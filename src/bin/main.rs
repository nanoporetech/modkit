use clap::Parser;
use log::error;
use mod_kit::commands::Commands;

#[derive(Parser)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    unsafe {
        rust_htslib::htslib::hts_set_log_level(
            rust_htslib::htslib::htsLogLevel_HTS_LOG_OFF,
        );
    }
    if let Err(err) = cli.command.run() {
        error!("Error! {err}");
        std::process::exit(1);
    }
    Ok(())
}
