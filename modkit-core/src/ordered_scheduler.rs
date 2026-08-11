use std::any::Any;
use std::collections::BTreeMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use anyhow::anyhow;
use crossbeam_channel::{bounded, unbounded, Receiver, Sender};

#[derive(Clone)]
struct FailureState {
    first: Arc<Mutex<Option<anyhow::Error>>>,
    signal: Sender<()>,
}

impl FailureState {
    fn report(&self, error: anyhow::Error) {
        let mut first =
            self.first.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
        if first.is_none() {
            *first = Some(error);
            let _ = self.signal.try_send(());
        }
    }

    fn has_failed(&self) -> bool {
        self.first
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .is_some()
    }

    fn take(&self) -> Option<anyhow::Error> {
        self.first
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .take()
    }
}

fn panic_message(panic: Box<dyn Any + Send>) -> String {
    if let Some(message) = panic.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = panic.downcast_ref::<String>() {
        message.clone()
    } else {
        "unknown panic payload".to_string()
    }
}

fn run_stage<R, F>(
    label: &str,
    failure: FailureState,
    mut resources: R,
    task: F,
) where
    F: FnOnce(&mut R) -> anyhow::Result<()>,
{
    // Keep every channel in `resources` alive until the originating failure is
    // recorded, so a derived disconnection cannot replace the root cause.
    let outcome = catch_unwind(AssertUnwindSafe(|| task(&mut resources)));
    match outcome {
        Ok(Ok(())) => {}
        Ok(Err(error)) => failure.report(anyhow!("{label} failed: {error:#}")),
        Err(panic) => failure
            .report(anyhow!("{label} panicked: {}", panic_message(panic))),
    }
}

fn join_stage(label: &str, handle: JoinHandle<()>, failure: &FailureState) {
    if let Err(panic) = handle.join() {
        failure.report(anyhow!(
            "{label} panicked outside its scheduler guard: {}",
            panic_message(panic)
        ));
    }
}

pub(crate) trait OrderedWorker<J, B>: Send {
    fn process(&mut self, job: J, buffer: B) -> anyhow::Result<B>;
}

pub(crate) fn run_ordered_scheduler<I, J, B, W, C>(
    stage_name: &'static str,
    feeder: I,
    workers: Vec<W>,
    empty_buffers: Receiver<B>,
    output_queue_size: usize,
    mut consume: C,
) -> anyhow::Result<()>
where
    I: Iterator<Item = anyhow::Result<J>> + Send + 'static,
    J: Send + 'static,
    B: Send + 'static,
    W: OrderedWorker<J, B> + 'static,
    C: FnMut(B) -> anyhow::Result<()>,
{
    if workers.is_empty() {
        return Err(anyhow!(
            "{stage_name} scheduler requires at least one worker"
        ));
    }

    let (jobs_tx, jobs_rx) = bounded(workers.len() * 2);
    let (results_tx, results_rx) = unbounded();
    let (records_tx, records_rx) = bounded(output_queue_size);
    // Dropping the only sender broadcasts cancellation to every blocked
    // receiver without consuming a one-shot cancellation message.
    let (cancel_tx, cancel_rx) = bounded::<()>(0);
    let (failure_signal_tx, failure_signal_rx) = bounded(1);
    let failure = FailureState {
        first: Arc::new(Mutex::new(None)),
        signal: failure_signal_tx,
    };

    let source_failure = failure.clone();
    let source_cancel = cancel_rx.clone();
    let source = std::thread::spawn(move || {
        run_stage(
            &format!("{stage_name} source"),
            source_failure,
            (feeder, jobs_tx, empty_buffers, source_cancel),
            |resources| {
                let (feeder, jobs_tx, empty_buffers, cancel) = resources;
                let mut seq = 0usize;
                loop {
                    let job = match feeder.next() {
                        Some(Ok(job)) => job,
                        Some(Err(error)) => {
                            return Err(anyhow!("feeder error: {error:#}"))
                        }
                        None => return Ok(()),
                    };
                    let buffer = crossbeam_channel::select_biased! {
                        recv(cancel) -> _ => return Ok(()),
                        recv(empty_buffers) -> result => result.map_err(|_| {
                            anyhow!("buffer-credit channel disconnected")
                        })?,
                    };
                    let sent = crossbeam_channel::select_biased! {
                        recv(cancel) -> _ => return Ok(()),
                        send(jobs_tx, (seq, job, buffer)) -> result => result,
                    };
                    sent.map_err(|_| anyhow!("job channel disconnected"))?;
                    seq = seq.wrapping_add(1);
                }
            },
        );
    });

    let mut worker_handles = Vec::with_capacity(workers.len());
    for (index, worker) in workers.into_iter().enumerate() {
        let worker_failure = failure.clone();
        let worker_cancel = cancel_rx.clone();
        let worker_jobs = jobs_rx.clone();
        let worker_results = results_tx.clone();
        worker_handles.push(std::thread::spawn(move || {
            let label = format!("{stage_name} worker {index}");
            run_stage(
                &label,
                worker_failure,
                (worker, worker_jobs, worker_results, worker_cancel),
                |resources| {
                    let (worker, jobs, results, cancel) = resources;
                    loop {
                        let received = crossbeam_channel::select_biased! {
                            recv(cancel) -> _ => return Ok(()),
                            recv(jobs) -> result => result,
                        };
                        let (seq, job, buffer) = match received {
                            Ok(job) => job,
                            Err(_) => return Ok(()),
                        };
                        let output = worker.process(job, buffer)?;
                        let sent = crossbeam_channel::select_biased! {
                            recv(cancel) -> _ => return Ok(()),
                            send(results, (seq, output)) -> result => result,
                        };
                        sent.map_err(|_| {
                            anyhow!("results channel disconnected")
                        })?;
                    }
                },
            );
        }));
    }
    drop(results_tx);
    drop(jobs_rx);

    let aggregator_failure = failure.clone();
    let aggregator_cancel = cancel_rx;
    let aggregator = std::thread::spawn(move || {
        run_stage(
            &format!("{stage_name} aggregator"),
            aggregator_failure,
            (results_rx, records_tx, aggregator_cancel),
            |resources| {
                let (results, records, cancel) = resources;
                let mut next_seq = 0usize;
                let mut buffer = BTreeMap::new();
                loop {
                    let received = crossbeam_channel::select_biased! {
                        recv(cancel) -> _ => return Ok(()),
                        recv(results) -> result => result,
                    };
                    let (seq, output) = match received {
                        Ok(result) => result,
                        Err(_) if buffer.is_empty() => return Ok(()),
                        Err(_) => {
                            return Err(anyhow!(
                                "results channel closed before sequence \
                                 {next_seq}"
                            ))
                        }
                    };
                    if buffer.insert(seq, output).is_some() {
                        return Err(anyhow!("duplicate result sequence {seq}"));
                    }
                    while let Some(output) = buffer.remove(&next_seq) {
                        let sent = crossbeam_channel::select_biased! {
                            recv(cancel) -> _ => return Ok(()),
                            send(records, output) -> result => result,
                        };
                        sent.map_err(|_| {
                            anyhow!("ordered-output channel disconnected")
                        })?;
                        next_seq = next_seq.wrapping_add(1);
                    }
                }
            },
        );
    });

    let mut cancel_tx = Some(cancel_tx);
    loop {
        crossbeam_channel::select_biased! {
            recv(failure_signal_rx) -> _ => {
                drop(cancel_tx.take());
                break;
            },
            recv(records_rx) -> result => match result {
                Ok(output) => {
                    match catch_unwind(AssertUnwindSafe(|| consume(output))) {
                        Ok(Ok(())) => {}
                        Ok(Err(error)) => {
                            failure.report(anyhow!(
                                "{stage_name} output consumer failed: {error:#}"
                            ));
                            drop(cancel_tx.take());
                            break;
                        }
                        Err(panic) => {
                            failure.report(anyhow!(
                                "{stage_name} output consumer panicked: {}",
                                panic_message(panic)
                            ));
                            drop(cancel_tx.take());
                            break;
                        }
                    }
                }
                Err(_) => {
                    if failure.has_failed() {
                        drop(cancel_tx.take());
                    }
                    break;
                }
            },
        }
    }
    drop(records_rx);
    // Cancellation makes every channel-blocked stage wake before joining.
    drop(cancel_tx.take());

    join_stage(&format!("{stage_name} source"), source, &failure);
    for (index, worker) in worker_handles.into_iter().enumerate() {
        join_stage(&format!("{stage_name} worker {index}"), worker, &failure);
    }
    join_stage(&format!("{stage_name} aggregator"), aggregator, &failure);

    match failure.take() {
        Some(error) => Err(error),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;
    use std::panic::{catch_unwind, AssertUnwindSafe};
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Barrier, Mutex};
    use std::thread;
    use std::time::Duration;

    use anyhow::{anyhow, bail};
    use crossbeam_channel::{bounded, unbounded, Receiver, Sender};

    use super::{run_ordered_scheduler, OrderedWorker};

    struct TestBuffer {
        value: usize,
        drops: Option<Arc<AtomicUsize>>,
    }

    impl Drop for TestBuffer {
        fn drop(&mut self) {
            if let Some(drops) = self.drops.as_ref() {
                drops.fetch_add(1, Ordering::SeqCst);
            }
        }
    }

    #[derive(Clone, Copy)]
    enum WorkerBehavior {
        Succeed,
        ErrorAt(usize),
        PanicAt(usize),
    }

    struct TestWorker {
        behavior: WorkerBehavior,
        drops: Option<Arc<AtomicUsize>>,
    }

    impl Drop for TestWorker {
        fn drop(&mut self) {
            if let Some(drops) = self.drops.as_ref() {
                drops.fetch_add(1, Ordering::SeqCst);
            }
        }
    }

    impl OrderedWorker<usize, TestBuffer> for TestWorker {
        fn process(
            &mut self,
            job: usize,
            mut buffer: TestBuffer,
        ) -> anyhow::Result<TestBuffer> {
            match self.behavior {
                WorkerBehavior::Succeed => {}
                WorkerBehavior::ErrorAt(target) if job == target => {
                    bail!("single worker failure at job {job}")
                }
                WorkerBehavior::PanicAt(target) if job == target => {
                    panic!("scripted worker panic at job {job}")
                }
                WorkerBehavior::ErrorAt(_) | WorkerBehavior::PanicAt(_) => {}
            }
            buffer.value = job;
            Ok(buffer)
        }
    }

    struct BlockingJoinWorker {
        job_two_started: Sender<()>,
        wait_for_job_two: Receiver<()>,
        release_job_two: Receiver<()>,
        failure_ready: Sender<()>,
        drops: Arc<AtomicUsize>,
    }

    impl Drop for BlockingJoinWorker {
        fn drop(&mut self) {
            self.drops.fetch_add(1, Ordering::SeqCst);
        }
    }

    impl OrderedWorker<usize, TestBuffer> for BlockingJoinWorker {
        fn process(
            &mut self,
            job: usize,
            mut buffer: TestBuffer,
        ) -> anyhow::Result<TestBuffer> {
            match job {
                0 => {
                    self.wait_for_job_two
                        .recv_timeout(Duration::from_secs(1))
                        .map_err(|_| anyhow!("job 2 never started"))?;
                    let _ = self.failure_ready.send(());
                    bail!("cancellation trigger at job 0")
                }
                2 => {
                    let _ = self.job_two_started.send(());
                    self.release_job_two
                        .recv_timeout(Duration::from_secs(1))
                        .map_err(|_| anyhow!("job 2 was not released"))?;
                }
                _ => {}
            }
            buffer.value = job;
            Ok(buffer)
        }
    }

    struct CoordinatedErrorWorker {
        ready: Arc<Barrier>,
    }

    impl OrderedWorker<usize, TestBuffer> for CoordinatedErrorWorker {
        fn process(
            &mut self,
            job: usize,
            _buffer: TestBuffer,
        ) -> anyhow::Result<TestBuffer> {
            self.ready.wait();
            bail!("scripted root failure at job {job}")
        }
    }

    struct ConsumerJoinWorker {
        job_one_started: Sender<()>,
        wait_for_job_one: Receiver<()>,
        release_job_one: Receiver<()>,
        drops: Arc<AtomicUsize>,
    }

    impl Drop for ConsumerJoinWorker {
        fn drop(&mut self) {
            self.drops.fetch_add(1, Ordering::SeqCst);
        }
    }

    impl OrderedWorker<usize, TestBuffer> for ConsumerJoinWorker {
        fn process(
            &mut self,
            job: usize,
            mut buffer: TestBuffer,
        ) -> anyhow::Result<TestBuffer> {
            match job {
                0 => {
                    self.wait_for_job_one
                        .recv_timeout(Duration::from_secs(1))
                        .map_err(|_| anyhow!("job 1 never started"))?;
                }
                1 => {
                    let _ = self.job_one_started.send(());
                    self.release_job_one
                        .recv_timeout(Duration::from_secs(1))
                        .map_err(|_| anyhow!("job 1 was not released"))?;
                }
                _ => {}
            }
            buffer.value = job;
            Ok(buffer)
        }
    }

    struct OrderedCompletionWorker {
        job_one_completed: Sender<()>,
        release_job_zero: Receiver<()>,
        completions: Arc<Mutex<Vec<usize>>>,
    }

    impl OrderedWorker<usize, TestBuffer> for OrderedCompletionWorker {
        fn process(
            &mut self,
            job: usize,
            mut buffer: TestBuffer,
        ) -> anyhow::Result<TestBuffer> {
            if job == 0 {
                self.release_job_zero
                    .recv_timeout(Duration::from_secs(1))
                    .map_err(|_| anyhow!("job 0 was not released"))?;
            }
            self.completions.lock().unwrap().push(job);
            if job == 1 {
                let _ = self.job_one_completed.send(());
            }
            buffer.value = job;
            Ok(buffer)
        }
    }

    #[derive(Clone, Copy)]
    enum ConsumerFailure {
        Error,
        Panic,
    }

    struct TestFeeder {
        items: VecDeque<anyhow::Result<usize>>,
        drops: Option<Arc<AtomicUsize>>,
    }

    impl TestFeeder {
        fn jobs(jobs: impl IntoIterator<Item = usize>) -> Self {
            Self { items: jobs.into_iter().map(Ok).collect(), drops: None }
        }
    }

    impl Iterator for TestFeeder {
        type Item = anyhow::Result<usize>;

        fn next(&mut self) -> Option<Self::Item> {
            self.items.pop_front()
        }
    }

    impl Drop for TestFeeder {
        fn drop(&mut self) {
            if let Some(drops) = self.drops.as_ref() {
                drops.fetch_add(1, Ordering::SeqCst);
            }
        }
    }

    fn make_buffer_pool(
        count: usize,
        drops: Option<Arc<AtomicUsize>>,
    ) -> (Sender<TestBuffer>, crossbeam_channel::Receiver<TestBuffer>) {
        let (sender, receiver) = unbounded();
        for _ in 0..count {
            sender
                .send(TestBuffer { value: usize::MAX, drops: drops.clone() })
                .unwrap();
        }
        (sender, receiver)
    }

    fn recycle_buffer(
        sender: &Sender<TestBuffer>,
        buffer: TestBuffer,
    ) -> anyhow::Result<()> {
        let _ = sender.send(buffer);
        Ok(())
    }

    fn run_scenario(
        feeder: TestFeeder,
        behavior: WorkerBehavior,
        worker_count: usize,
        buffer_count: usize,
    ) -> anyhow::Result<()> {
        let workers = (0..worker_count)
            .map(|_| TestWorker { behavior, drops: None })
            .collect();
        let (empty_sender, empty_buffers) =
            make_buffer_pool(buffer_count, None);
        let recycle = empty_sender.clone();
        run_ordered_scheduler(
            "pileup",
            feeder,
            workers,
            empty_buffers,
            worker_count * 2,
            move |buffer| recycle_buffer(&recycle, buffer),
        )
    }

    fn run_with_watchdog<F>(timeout: Duration, task: F) -> anyhow::Result<()>
    where
        F: FnOnce() -> anyhow::Result<()> + Send + 'static,
    {
        let (sender, receiver) = bounded(1);
        thread::spawn(move || {
            let result =
                catch_unwind(AssertUnwindSafe(task)).map_err(|panic| {
                    if let Some(message) = panic.downcast_ref::<&str>() {
                        anyhow!("scheduler panicked: {message}")
                    } else if let Some(message) = panic.downcast_ref::<String>()
                    {
                        anyhow!("scheduler panicked: {message}")
                    } else {
                        anyhow!("scheduler panicked")
                    }
                });
            let _ = sender.send(result.and_then(|result| result));
        });
        receiver.recv_timeout(timeout).map_err(|_| {
            anyhow!("scheduler watchdog expired after {timeout:?}")
        })?
    }

    fn assert_consumer_failure_joins_all_threads(mode: ConsumerFailure) {
        let worker_drops = Arc::new(AtomicUsize::new(0));
        let feeder_drops = Arc::new(AtomicUsize::new(0));
        let buffer_drops = Arc::new(AtomicUsize::new(0));
        let (job_one_started, wait_for_job_one) = bounded(1);
        let (release_job_one, job_one_release) = bounded(1);
        let (consumer_called, wait_for_consumer) = bounded(1);
        let (finished, wait_for_finish) = bounded(1);
        let worker_drop_probe = worker_drops.clone();
        let feeder_drop_probe = feeder_drops.clone();
        let buffer_drop_probe = buffer_drops.clone();

        thread::spawn(move || {
            let feeder = TestFeeder {
                items: (0..64).map(Ok).collect(),
                drops: Some(feeder_drop_probe),
            };
            let workers = (0..2)
                .map(|_| ConsumerJoinWorker {
                    job_one_started: job_one_started.clone(),
                    wait_for_job_one: wait_for_job_one.clone(),
                    release_job_one: job_one_release.clone(),
                    drops: worker_drop_probe.clone(),
                })
                .collect();
            let (empty_sender, empty_buffers) =
                make_buffer_pool(4, Some(buffer_drop_probe));
            let outcome = catch_unwind(AssertUnwindSafe(|| {
                run_ordered_scheduler(
                    "pileup",
                    feeder,
                    workers,
                    empty_buffers,
                    0,
                    move |_buffer| {
                        let _ = consumer_called.send(());
                        match mode {
                            ConsumerFailure::Error => {
                                bail!("scripted consumer failure")
                            }
                            ConsumerFailure::Panic => {
                                panic!("scripted consumer panic")
                            }
                        }
                    },
                )
            }));
            drop(empty_sender);
            let _ = finished.send(outcome);
        });

        wait_for_consumer
            .recv_timeout(Duration::from_secs(1))
            .expect("consumer was never called");
        if wait_for_finish.recv_timeout(Duration::from_millis(50)).is_ok() {
            let _ = release_job_one.send(());
            panic!("scheduler returned before joining the blocked worker");
        }
        release_job_one.send(()).expect("failed to release blocked worker");
        let outcome = wait_for_finish
            .recv_timeout(Duration::from_secs(1))
            .expect("scheduler did not finish after cancellation");
        let result = match outcome {
            Ok(result) => result,
            Err(_) => panic!("consumer panic escaped the scheduler"),
        };
        let error = result.expect_err("consumer failure should be propagated");
        match mode {
            ConsumerFailure::Error => {
                assert!(error.to_string().contains("scripted consumer failure"))
            }
            ConsumerFailure::Panic => {
                assert!(error.to_string().contains(
                    "pileup output consumer panicked: scripted consumer panic"
                ))
            }
        }
        assert_eq!(worker_drops.load(Ordering::SeqCst), 2);
        assert_eq!(feeder_drops.load(Ordering::SeqCst), 1);
        assert_eq!(buffer_drops.load(Ordering::SeqCst), 4);
    }

    #[test]
    fn zero_workers_are_rejected_without_feeding_or_consuming() {
        let feeder_calls = Arc::new(AtomicUsize::new(0));
        let consumer_calls = Arc::new(AtomicUsize::new(0));
        let feeder_probe = feeder_calls.clone();
        let consumer_probe = consumer_calls.clone();
        let error = run_with_watchdog(Duration::from_secs(1), move || {
            let feeder = std::iter::from_fn(move || {
                feeder_probe.fetch_add(1, Ordering::SeqCst);
                Some(Ok(0))
            });
            let (_empty_sender, empty_buffers) = unbounded();
            run_ordered_scheduler(
                "pileup",
                feeder,
                Vec::<TestWorker>::new(),
                empty_buffers,
                0,
                move |_buffer: TestBuffer| {
                    consumer_probe.fetch_add(1, Ordering::SeqCst);
                    Ok(())
                },
            )
        })
        .expect_err("an empty worker set should be rejected");

        assert!(error.to_string().contains("at least one worker"));
        assert_eq!(feeder_calls.load(Ordering::SeqCst), 0);
        assert_eq!(consumer_calls.load(Ordering::SeqCst), 0);
    }

    #[test]
    fn propagates_single_worker_error() {
        let error = run_with_watchdog(Duration::from_secs(1), || {
            run_scenario(
                TestFeeder::jobs([0]),
                WorkerBehavior::ErrorAt(0),
                1,
                2,
            )
        })
        .expect_err("worker error should be propagated");
        assert!(error.to_string().contains("single worker failure at job 0"));
    }

    #[test]
    fn concurrent_worker_errors_preserve_a_scripted_root_cause() {
        let error = run_with_watchdog(Duration::from_secs(1), || {
            let ready = Arc::new(Barrier::new(2));
            let workers = (0..2)
                .map(|_| CoordinatedErrorWorker { ready: ready.clone() })
                .collect();
            let (empty_sender, empty_buffers) = make_buffer_pool(4, None);
            let recycle = empty_sender.clone();
            run_ordered_scheduler(
                "pileup",
                TestFeeder::jobs(0..64),
                workers,
                empty_buffers,
                4,
                move |buffer| recycle_buffer(&recycle, buffer),
            )
        })
        .expect_err("a coordinated worker error should be propagated");
        let message = error.to_string();
        assert!(
            message.contains("scripted root failure at job 0")
                || message.contains("scripted root failure at job 1"),
            "unexpected error: {message}"
        );
        assert!(
            !message.contains("disconnected"),
            "derived error won: {message}"
        );
    }

    #[test]
    fn worker_panic_is_propagated_without_sequence_hole_deadlock() {
        let error = run_with_watchdog(Duration::from_secs(1), || {
            run_scenario(
                TestFeeder::jobs(0..64),
                WorkerBehavior::PanicAt(0),
                1,
                2,
            )
        })
        .expect_err("worker panic should be propagated");
        assert!(error.to_string().contains("scripted worker panic at job 0"));
    }

    #[test]
    fn feeder_error_is_propagated() {
        let error = run_with_watchdog(Duration::from_secs(1), || {
            run_scenario(
                TestFeeder {
                    items: vec![Err(anyhow!("scripted feeder failure"))].into(),
                    drops: None,
                },
                WorkerBehavior::Succeed,
                1,
                2,
            )
        })
        .expect_err("feeder error should be propagated");
        assert!(error.to_string().contains("scripted feeder failure"));
    }

    #[test]
    fn disconnected_buffer_credit_channel_is_propagated() {
        let error = run_with_watchdog(Duration::from_secs(1), || {
            let (empty_sender, empty_buffers) = unbounded::<TestBuffer>();
            drop(empty_sender);
            run_ordered_scheduler(
                "pileup",
                TestFeeder::jobs([0]),
                vec![TestWorker {
                    behavior: WorkerBehavior::Succeed,
                    drops: None,
                }],
                empty_buffers,
                2,
                |_| Ok(()),
            )
        })
        .expect_err("buffer-credit channel failure should be propagated");
        assert!(error.to_string().contains("buffer"));
    }

    #[test]
    fn failure_cancels_blocked_stages_and_joins_all_threads() {
        let worker_drops = Arc::new(AtomicUsize::new(0));
        let feeder_drops = Arc::new(AtomicUsize::new(0));
        let buffer_drops = Arc::new(AtomicUsize::new(0));
        let (job_two_started, wait_for_job_two) = bounded(1);
        let (release_job_two, job_two_release) = bounded(1);
        let (failure_ready, wait_for_failure) = bounded(1);
        let (finished, wait_for_finish) = bounded(1);
        let worker_drop_probe = worker_drops.clone();
        let feeder_drop_probe = feeder_drops.clone();
        let buffer_drop_probe = buffer_drops.clone();

        thread::spawn(move || {
            let feeder = TestFeeder {
                items: (0..64).map(Ok).collect(),
                drops: Some(feeder_drop_probe),
            };
            let workers = (0..2)
                .map(|_| BlockingJoinWorker {
                    job_two_started: job_two_started.clone(),
                    wait_for_job_two: wait_for_job_two.clone(),
                    release_job_two: job_two_release.clone(),
                    failure_ready: failure_ready.clone(),
                    drops: worker_drop_probe.clone(),
                })
                .collect();
            let (empty_sender, empty_buffers) =
                make_buffer_pool(4, Some(buffer_drop_probe));
            let recycle = empty_sender.clone();
            let result = run_ordered_scheduler(
                "pileup",
                feeder,
                workers,
                empty_buffers,
                1,
                move |buffer| recycle_buffer(&recycle, buffer),
            );
            drop(empty_sender);
            let _ = finished.send(result);
        });

        wait_for_failure
            .recv_timeout(Duration::from_secs(1))
            .expect("scripted worker failure never became ready");
        assert!(
            wait_for_finish.recv_timeout(Duration::from_millis(50)).is_err(),
            "scheduler returned before joining the blocked worker"
        );
        release_job_two.send(()).expect("failed to release blocked worker");
        let error = wait_for_finish
            .recv_timeout(Duration::from_secs(1))
            .expect("scheduler did not finish after cancellation")
            .expect_err("worker failure should cancel the pipeline");
        assert!(error.to_string().contains("cancellation trigger at job 0"));
        assert_eq!(worker_drops.load(Ordering::SeqCst), 2);
        assert_eq!(feeder_drops.load(Ordering::SeqCst), 1);
        assert_eq!(buffer_drops.load(Ordering::SeqCst), 4);
    }

    #[test]
    fn consumer_error_cancels_and_joins_all_threads() {
        assert_consumer_failure_joins_all_threads(ConsumerFailure::Error);
    }

    #[test]
    fn consumer_panic_is_converted_and_joins_all_threads() {
        assert_consumer_failure_joins_all_threads(ConsumerFailure::Panic);
    }

    #[test]
    fn out_of_order_completions_are_consumed_in_feeder_order() {
        let completions = Arc::new(Mutex::new(Vec::new()));
        let consumed = Arc::new(Mutex::new(Vec::new()));
        let (job_one_completed, wait_for_job_one) = bounded(1);
        let (release_job_zero, job_zero_release) = bounded(1);
        let (finished, wait_for_finish) = bounded(1);
        let completion_probe = completions.clone();
        let consumed_probe = consumed.clone();

        thread::spawn(move || {
            let workers = (0..2)
                .map(|_| OrderedCompletionWorker {
                    job_one_completed: job_one_completed.clone(),
                    release_job_zero: job_zero_release.clone(),
                    completions: completion_probe.clone(),
                })
                .collect();
            let (empty_sender, empty_buffers) = make_buffer_pool(2, None);
            let recycle = empty_sender.clone();
            let result = run_ordered_scheduler(
                "pileup",
                TestFeeder::jobs(0..2),
                workers,
                empty_buffers,
                0,
                move |buffer| {
                    consumed_probe.lock().unwrap().push(buffer.value);
                    recycle_buffer(&recycle, buffer)
                },
            );
            drop(empty_sender);
            let _ = finished.send(result);
        });

        wait_for_job_one
            .recv_timeout(Duration::from_secs(1))
            .expect("job 1 did not complete");
        assert_eq!(*completions.lock().unwrap(), vec![1]);
        assert!(consumed.lock().unwrap().is_empty());
        assert!(wait_for_finish
            .recv_timeout(Duration::from_millis(50))
            .is_err());
        release_job_zero.send(()).expect("failed to release job 0");
        wait_for_finish
            .recv_timeout(Duration::from_secs(1))
            .expect("scheduler did not finish")
            .expect("successful scheduler run should complete");

        assert_eq!(*completions.lock().unwrap(), vec![1, 0]);
        assert_eq!(*consumed.lock().unwrap(), vec![0, 1]);
    }
}
