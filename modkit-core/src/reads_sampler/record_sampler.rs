use crate::util::{get_master_progress_bar, get_ticker};

use indicatif::ProgressBar;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rust_htslib::bam;

use super::deterministic_sampler::DeterministicFractionSampler;

/// A utility data structure that when used in an interator allows
/// to randomly sample either a preset number of reads or a fraction
/// of reads. If sampling a preset number of reads, say N reads, then
/// the first N reads are taken.
pub struct RecordSampler {
    pub(crate) num_reads: Option<usize>,
    pub(crate) sample_frac: Option<f64>,
    deterministic_fraction_sampler: Option<DeterministicFractionSampler>,
    rng: StdRng,
    reads_sampled: usize,
}

impl RecordSampler {
    pub(crate) fn new_num_reads(num_reads: usize) -> Self {
        Self {
            num_reads: Some(num_reads),
            sample_frac: None,
            deterministic_fraction_sampler: None,
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
            deterministic_fraction_sampler: None,
            rng,
            reads_sampled: 0,
        }
    }

    /// Construct a stateless fractional sampler for indexed record fetching.
    /// Every worker in the indexed job must receive the same resolved seed.
    pub(crate) fn new_deterministic_sample_frac(
        sample_frac: f64,
        master_seed: u64,
    ) -> Self {
        let deterministic_fraction_sampler = Some(
            DeterministicFractionSampler::new(master_seed, sample_frac)
                .expect("sampling fraction should already be validated"),
        );
        Self {
            num_reads: None,
            sample_frac: Some(sample_frac),
            deterministic_fraction_sampler,
            rng: StdRng::seed_from_u64(master_seed),
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_passthrough() -> Self {
        Self {
            num_reads: None,
            sample_frac: None,
            deterministic_fraction_sampler: None,
            rng: StdRng::from_entropy(),
            reads_sampled: 0,
        }
    }

    pub fn new_from_options(
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
            get_ticker()
        };
        spinner.set_message("records sampled");
        spinner
    }

    fn check_num_reads(&mut self) -> Indicator {
        if self.reads_sampled >= self.num_reads.unwrap() {
            Indicator::Done
        } else {
            Indicator::Use(Token)
        }
    }

    fn check_sample_frac(&mut self) -> Indicator {
        if self.rng.gen_bool(self.sample_frac.unwrap()) {
            Indicator::Use(Token)
        } else {
            Indicator::Skip
        }
    }

    pub(crate) fn ask(&mut self) -> Indicator {
        debug_assert!(
            self.deterministic_fraction_sampler.is_none(),
            "deterministic sampling decisions require a BAM record"
        );
        match (self.num_reads, self.sample_frac) {
            (Some(_nr), _) => self.check_num_reads(),
            (_, Some(_sample_frac)) => self.check_sample_frac(),
            (None, None) => Indicator::Use(Token),
        }
    }

    pub(crate) fn ask_record(&mut self, record: &bam::Record) -> Indicator {
        if let Some(sampler) = self.deterministic_fraction_sampler {
            if sampler.include(record) {
                Indicator::Use(Token)
            } else {
                Indicator::Skip
            }
        } else {
            self.ask()
        }
    }

    pub(crate) fn used(&mut self, _token: Token) {
        self.reads_sampled += 1;
    }
}

pub(crate) struct Token;

pub(crate) enum Indicator {
    Use(Token),
    Skip,
    Done,
}
