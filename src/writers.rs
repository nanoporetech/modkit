use crate::mod_pileup::ModBasePileup;
use std::fs::File;
use std::io::{BufWriter, Write};

pub trait OutWriter<T> {
    fn write(&mut self, item: T) -> Result<u64, String>;
}

pub struct BEDWriter {
    buf_writer: BufWriter<File>,
}

impl BEDWriter {
    pub fn new(buf_writer: BufWriter<File>) -> Self {
        Self { buf_writer }
    }
}

impl OutWriter<ModBasePileup> for BEDWriter {
    fn write(&mut self, item: ModBasePileup) -> Result<u64, String> {
        let mut rows_written = 0;
        let sep = '\t';
        for (pos, feature_counts) in item.iter_counts() {
            // chrom
            // pos
            // pos+1
            // raw_code
            // filtered cov (score)
            // strand
            // pos
            // pos+1
            // 0,0,0
            // filtered_cov
            // percent mod
            // n canonical
            // n modified
            // n other mod
            // n delete
            // n filtered
            // n diff
            for feature_count in feature_counts {
                let row = format!(
                    "{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}\
                    {}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}\n",
                    item.chrom_name,
                    pos,
                    pos + 1,
                    feature_count.raw_mod_code,
                    feature_count.filtered_coverage,
                    feature_count.strand.to_char(),
                    pos,
                    pos + 1,
                    "0,0,0",
                    feature_count.filtered_coverage,
                    format!("{:.0}", feature_count.fraction_modified * 100f32),
                    feature_count.n_modified,
                    feature_count.n_canonical,
                    feature_count.n_other_modified,
                    feature_count.n_delete,
                    feature_count.n_filtered,
                    feature_count.n_diff
                );
                self.buf_writer
                    .write(row.as_bytes())
                    .map_err(|e| e.to_string())?;
                rows_written += 1;
            }
        }
        Ok(rows_written)
    }
}
