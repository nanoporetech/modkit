use clap::Parser;
use log::error;
use mod_kit::commands::Commands;

#[derive(Parser)]
#[command(version)]
/// Modkit is a bioinformatics tool for working with modified bases from Oxford
/// Nanopore
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
        for cause in err.chain().skip(1) {
            error!(" caused by {cause}")
        }
        std::process::exit(1);
    }
    Ok(())
}
