use rust_htslib::bam::{self, HeaderView};
use rust_htslib::htslib::bam1_t;
use std::rc::Rc;

pub struct SafeRecord {
    inner: bam1_t,
}

unsafe impl Send for SafeRecord {}
unsafe impl Sync for SafeRecord {}

impl SafeRecord {
    pub fn ingest(read: &mut bam::Record) -> Self {
        let dummy = bam::Record::new();
        let inner = std::mem::replace(&mut read.inner, dummy.inner);
        Self { inner }
    }

    pub fn to_hts_record(self, header: Option<Rc<HeaderView>>) -> bam::Record {
        let mut record: bam::Record = self.into();
        if let Some(header) = header {
            record.set_header(header);
        }
        record
    }
}

impl From<bam::Record> for SafeRecord {
    fn from(mut value: bam::Record) -> Self {
        // necessary to manually move the data out otherwise `bam::Record` will
        // try to free the data when it's dropped. If the `own`
        // attribute was exposed, this wouldn't be necessary
        let dummy = bam::Record::new();
        let inner = std::mem::replace(&mut value.inner, dummy.inner);
        Self { inner }
    }
}

impl Into<bam::Record> for SafeRecord {
    fn into(self) -> bam::Record {
        let mut record = bam::Record::new();
        let data = record.inner_mut();
        let _ = std::mem::replace(data, self.inner);
        record
    }
}

#[cfg(test)]
pub mod saferecord_tests {
    use crate::SafeRecord;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;

    pub fn get_n_tests() -> usize {
        std::env::var("SAFE_RECORD_N_TESTS")
            .map_err(|e| e.to_string())
            .and_then(|x| x.parse::<usize>().map_err(|e| e.to_string()))
            .unwrap_or(20)
    }

    #[test]
    fn test_channel_reading() {
        let n_tests = get_n_tests();
        for _i in 0..n_tests {
            let mut reader =
                bam::Reader::from_path("tests/resources/10kreads.bam").unwrap();
            reader.set_threads(4).unwrap();
            let (snd, rcv) = std::sync::mpsc::channel();

            std::thread::spawn(move || {
                for record in
                    reader.records().filter_map(|r| r.ok()).take(10000)
                {
                    let safe_record = SafeRecord::from(record);
                    snd.send(safe_record).unwrap();
                }
            });

            let mut agg = 0;
            let mut n_recs = 0;
            for record in rcv {
                let record: bam::Record = record.into();
                agg += record.seq_len();
                n_recs += 1;
            }

            assert_eq!(agg, 88098576);
            assert_eq!(n_recs, 10000);
        }
    }

    #[test]
    fn test_record_buffer() {
        for i in 0..get_n_tests() {
            let mut reader =
                bam::Reader::from_path("tests/resources/10kreads.bam").unwrap();
            reader.set_threads(4).unwrap();
            let (snd, rcv) = std::sync::mpsc::channel();

            std::thread::spawn(move || {
                let mut record = bam::Record::new();
                for i in 0..10000 {
                    if let Some(_) = reader.read(&mut record) {
                        let safe_record = SafeRecord::ingest(&mut record);
                        snd.send(safe_record).unwrap();
                    }
                }
            });

            let mut agg = 0;
            let mut n_recs = 0;
            for record in rcv {
                let record: bam::Record = record.into();
                agg += record.seq_len();
                n_recs += 1;
            }

            assert_eq!(agg, 88098576);
            assert_eq!(n_recs, 10000);
        }
    }
}
