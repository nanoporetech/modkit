use clap::Parser;
use mod_kit::commands::Commands;

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    if let Err(err) = cli.command.run() {
        eprintln!("Error! {err}");
        std::process::exit(1);
    }
    Ok(())
}
