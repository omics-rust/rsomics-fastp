mod cli;
use clap::Parser;
use rsomics_common::Tool;

fn main() -> std::process::ExitCode {
    let args = cli::Cli::parse();
    args.run()
}
