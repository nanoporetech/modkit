use log::LevelFilter;
use log4rs::append::console::{ConsoleAppender, Target};
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use log4rs::filter::threshold::ThresholdFilter;
use log4rs::{Config, Handle};
use std::path::PathBuf;

pub fn init_logging(log_fp: Option<&PathBuf>) -> Handle {
    let level = LevelFilter::Info;

    let file_endcoder = Box::new(PatternEncoder::new("[{f}::{L}][{l}] {m}{n}"));
    let console_encoder = Box::new(PatternEncoder::new("> {m}{n}"));
    let stderr = ConsoleAppender::builder()
        .encoder(console_encoder)
        .target(Target::Stderr)
        .build();

    let config = if let Some(fp) = log_fp {
        let logfile = FileAppender::builder()
            .encoder(file_endcoder)
            .build(fp)
            .unwrap();
        Config::builder()
            .appender(
                Appender::builder()
                    .filter(Box::new(ThresholdFilter::new(level)))
                    .build("stderr", Box::new(stderr)),
            )
            .appender(Appender::builder().build("logfile", Box::new(logfile)))
            .build(
                Root::builder()
                    .appender("logfile")
                    .appender("stderr")
                    .build(LevelFilter::Trace),
            )
            .unwrap()
    } else {
        Config::builder()
            .appender(
                Appender::builder()
                    .filter(Box::new(ThresholdFilter::new(level)))
                    .build("stderr", Box::new(stderr)),
            )
            .build(Root::builder().appender("stderr").build(LevelFilter::Trace))
            .unwrap()
    };

    let handle = log4rs::init_config(config).expect("failed to init logging");

    handle
}
