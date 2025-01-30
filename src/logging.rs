use log::{debug, LevelFilter};
use log4rs::append::console::{ConsoleAppender, Target};
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use log4rs::filter::threshold::ThresholdFilter;
use log4rs::{Config, Handle};
use std::path::PathBuf;

pub fn init_logging_smart(
    log_fp: Option<&PathBuf>,
    quiet_stdout: bool,
) -> Handle {
    let level = LevelFilter::Info;

    let file_endcoder = Box::new(PatternEncoder::new(
        "[{f}::{L}][{d(%Y-%m-%d %H:%M:%S)}][{l}] {m}{n}",
    ));
    let console_encoder = Box::new(PatternEncoder::new("{h(>)} {m}{n}"));
    let stderr = ConsoleAppender::builder()
        .encoder(console_encoder)
        .target(Target::Stderr)
        .build();

    let config = if let Some(fp) = log_fp {
        let logfile =
            FileAppender::builder().encoder(file_endcoder).build(fp).unwrap();
        let mut config = Config::builder();
        let logfile_appender =
            Appender::builder().build("logfile", Box::new(logfile));
        config = config.appender(logfile_appender);
        let mut root_logger = Root::builder().appender("logfile");
        if !quiet_stdout {
            config = config.appender(
                Appender::builder()
                    .filter(Box::new(ThresholdFilter::new(level)))
                    .build("stderr", Box::new(stderr)),
            );
            root_logger = root_logger.appender("stderr");
        }

        config.build(root_logger.build(LevelFilter::Trace)).unwrap()
    } else if !quiet_stdout {
        Config::builder()
            .appender(
                Appender::builder()
                    .filter(Box::new(ThresholdFilter::new(level)))
                    .build("stderr", Box::new(stderr)),
            )
            .build(Root::builder().appender("stderr").build(LevelFilter::Trace))
            .unwrap()
    } else {
        Config::builder()
            .build(Root::builder().build(LevelFilter::Trace))
            .unwrap()
    };

    let handle = log4rs::init_config(config).expect("failed to init logging");
    let command_line = std::env::args().collect::<Vec<String>>().join(" ");
    debug!("command line: {command_line}");
    handle
}

pub fn init_logging(log_fp: Option<&PathBuf>) -> Handle {
    init_logging_smart(log_fp, false)
}
