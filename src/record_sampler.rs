use crate::util::{get_master_progress_bar, get_spinner};
use indicatif::ProgressBar;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// A utility data structure that when used in an interator allows
/// to randomly sample either a preset number of reads or a fraction
/// of reads. If sampling a preset number of reads, say N reads, then
/// the first N reads are taken.
pub struct RecordSampler {
    pub(crate) num_reads: Option<usize>,
    pub(crate) sample_frac: Option<f64>,
    rng: StdRng,
    reads_sampled: usize,
}

impl RecordSampler {
    pub(crate) fn new_num_reads(num_reads: usize) -> Self {
        Self {
            num_reads: Some(num_reads),
            sample_frac: None,
            rng: StdRng::from_entropy(),
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_sample_frac(sample_frac: f64, seed: Option<u64>) -> Self {
        let rng = seed
            .map(|s| StdRng::seed_from_u64(s))
            .unwrap_or(StdRng::from_entropy());
        Self {
            num_reads: None,
            sample_frac: Some(sample_frac),
            rng,
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_passthrough() -> Self {
        Self {
            num_reads: None,
            sample_frac: None,
            rng: StdRng::from_entropy(),
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_from_options(
        sample_frac: Option<f64>,
        num_reads: Option<usize>,
        seed: Option<u64>,
    ) -> Self {
        match (sample_frac, num_reads) {
            (_, Some(num_reads)) => RecordSampler::new_num_reads(num_reads),
            (Some(f), _) => RecordSampler::new_sample_frac(f, seed),
            (None, None) => RecordSampler::new_passthrough(),
        }
    }

    pub(crate) fn get_progress_bar(&self) -> ProgressBar {
        let spinner = if let Some(num) = self.num_reads {
            get_master_progress_bar(num)
        } else {
            get_spinner()
        };
        spinner.set_message("records sampled");
        spinner
    }

    fn check_num_reads(&self) -> Indicator {
        let indicator = if self.reads_sampled >= self.num_reads.unwrap() {
            Indicator::Done
        } else {
            Indicator::Use
        };
        indicator
    }

    fn check_sample_frac(&mut self) -> Indicator {
        if self.rng.gen_bool(self.sample_frac.unwrap()) {
            Indicator::Use
        } else {
            Indicator::Skip
        }
    }

    pub(crate) fn ask(&mut self) -> Indicator {
        match (self.num_reads, self.sample_frac) {
            (Some(_nr), _) => self.check_num_reads(),
            (_, Some(_sample_frac)) => self.check_sample_frac(),
            (None, None) => Indicator::Use,
        }
    }
}

pub(crate) enum Indicator {
    Use,
    Skip,
    Done,
}
