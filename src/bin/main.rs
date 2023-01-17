use clap::Parser;
use log::error;
use mod_kit::commands::Commands;

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    if let Err(err) = cli.command.run() {
        error!("Error! {err}");
        std::process::exit(1);
    }
    Ok(())
}
