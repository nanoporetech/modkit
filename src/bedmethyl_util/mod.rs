use std::io::BufRead;

use anyhow::bail;
use bigtools::{
    bed::bedparser::{BedValueError, StreamingBedValues},
    Value,
};
use indicatif::ProgressBar;
use log::debug;
use rustc_hash::FxHashSet;

use crate::{
    dmr::bedmethyl::BedMethylLine, mod_base_code::ModCodeRepr, util::StrandRule,
};

pub mod subcommands;
struct BedMethylStream<R: BufRead> {
    in_stream: R,
    buf: String,
    mod_codes: FxHashSet<ModCodeRepr>,
    negative_strand_values: bool,
    curr_record: Option<BedMethylLine>,
    next_record: Option<BedMethylLine>,
    counter: ProgressBar,
}

impl<R: BufRead> BedMethylStream<R> {
    fn new(
        mut in_stream: R,
        mod_codes: FxHashSet<ModCodeRepr>,
        negative_strand_values: bool,
        counter: ProgressBar,
    ) -> anyhow::Result<Self> {
        let mut first_record = None;
        let mut buf = String::new();
        loop {
            buf.clear();
            let b = in_stream.read_line(&mut buf)?;
            if b == 0 {
                break;
            }
            let bm_record = BedMethylLine::parse(&buf)?;
            let keep = mod_codes.contains(&bm_record.raw_mod_code);
            if keep {
                first_record = Some(bm_record);
                break;
            } else {
                continue;
            }
        }
        if first_record.is_none() {
            bail!("no bedmethyl lines")
        } else {
            Ok(Self {
                in_stream,
                buf,
                mod_codes,
                next_record: first_record,
                curr_record: None,
                negative_strand_values,
                counter,
            })
        }
    }
}

impl<R: BufRead> BedMethylStream<R> {
    fn get_next(&mut self) -> Result<Option<(&str, Value)>, BedValueError> {
        if let Some(mut record) = self.next_record.take() {
            loop {
                // slurp up the next bedmethyl line
                self.buf.clear();
                let b = self.in_stream.read_line(&mut self.buf)?;
                if b == 0 {
                    // finished
                    assert!(self.next_record.is_none());
                    break;
                }
                let next = BedMethylLine::parse(&self.buf)
                    .map_err(|e| BedValueError::InvalidInput(e.to_string()))?;

                let keep = self.mod_codes.contains(&next.raw_mod_code);
                if !keep {
                    continue;
                }

                if record.same_position_and_strand_as(&next) {
                    if next.raw_mod_code == record.raw_mod_code {
                        return Err(BedValueError::InvalidInput(format!(
                            "duplicated record at {next}"
                        )));
                    }
                    if next.valid_coverage != record.valid_coverage {
                        return Err(BedValueError::InvalidInput(format!(
                            "invalid overlapping records at {next}, valid \
                             coverage should be the same if they apply to the \
                             same primary base"
                        )));
                    }
                    if (next.valid_coverage == record.valid_coverage)
                        && (next.count_canonical == record.count_canonical)
                    {
                        debug!(
                            "combining modification counts for {record} and \
                             {next}"
                        );
                        record.count_methylated += next.count_methylated;
                    }
                } else if record.same_position_as(&next) {
                    // same genomic position, different strand
                    if next.valid_coverage == record.valid_coverage {
                        debug!(
                            "ambiguous data at {record}, {next}, cannot \
                             determine correct strand to use"
                        );
                    }
                    if next.valid_coverage > record.valid_coverage {
                        debug!("replacing {record} with {next}");
                        record = next;
                    }
                } else {
                    self.next_record = Some(next);
                    break;
                }
            }
            let fact = if self.negative_strand_values
                && record.strand.overlaps(&StrandRule::Negative)
            {
                -100f32
            } else {
                100f32
            };

            let v = Value {
                start: record.start() as u32,
                end: record.stop() as u32,
                value: record.frac_modified() * fact,
            };
            self.curr_record = Some(record);
            self.counter.inc(1);
            Ok(Some((self.curr_record.as_ref().unwrap().chrom.as_str(), v)))
        } else {
            Ok(None)
        }
    }
}

impl<R: BufRead> StreamingBedValues for BedMethylStream<R> {
    type Value = Value;

    fn next(&mut self) -> Option<Result<(&str, Self::Value), BedValueError>> {
        self.get_next().transpose()
    }
}

#[cfg(test)]
mod bedmethylutil_tests {
    use std::io::BufReader;

    use rustc_hash::FxHashSet;

    use crate::bedmethyl_util::BedMethylStream;
    use crate::util::get_ticker;
    use crate::{
        dmr::bedmethyl::BedMethylLine, mod_base_code::ModCodeRepr,
        position_filter::Iv, tabix::ParseBedLine, util::StrandRule,
    };

    #[test]
    fn test_bedmethyl_stream_semantics() {
        let lines = vec![
            BedMethylLine::new(
                "chrom".to_string(),
                Iv { start: 0, stop: 1, val: () },
                ModCodeRepr::Code('m'),
                StrandRule::Positive,
                5,
                10,
                5,
                0,
                0,
                0,
                0,
                0,
            ),
            BedMethylLine::new(
                "chrom".to_string(),
                Iv { start: 0, stop: 1, val: () },
                ModCodeRepr::Code('a'),
                StrandRule::Positive,
                0,
                5,
                10,
                0,
                0,
                0,
                0,
                0,
            ),
        ];
        let lines = lines.into_iter().fold(String::new(), |mut acc, next| {
            acc.push_str(&next.to_line());
            acc
        });
        let reader = BufReader::new(lines.as_bytes());
        let mod_codes = FxHashSet::<ModCodeRepr>::from_iter(
            ['a'.into(), 'm'.into()].into_iter(),
        );
        let mut stream =
            BedMethylStream::new(reader, mod_codes, false, get_ticker())
                .unwrap();
        let _conflict_err = stream.get_next().expect_err("shouldn't work");
    }
}
