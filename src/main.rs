// Copyright 2021-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

mod cli;
mod utils;

use clap::Parser;
use log::info;

use anyhow::Context;
use std::io;

use std::time::Instant;

fn main() -> anyhow::Result<()> {
    let start_time = Instant::now();
    let stderr = io::stderr();
    let mut ehandle = stderr.lock();

    let cli = cli::Args::parse();
    utils::setup_logging(cli.quiet)?;

    // Handle input file/STDIN
    let infile = utils::handle_input(&cli.file, &mut ehandle)?;

    // Check and handle output files
    utils::handle_output_files(&cli.prefix, cli.force, &mut ehandle)?;

    // Process primers
    let primers = utils::process_primers(&cli, &mut ehandle)?;
    utils::validate_mismatch(&primers, cli.mismatch)?;

    // Log program information
    utils::log_program_info(cli.mismatch, cli.force, &cli.prefix);

    // Core processing
    utils::get_hypervar_regions(&infile, primers, &cli.prefix, cli.mismatch)
        .context("Failed to extract hypervariable regions")?;
    info!("Done getting hypervariable regions");

    // Cleanup and final logging
    utils::cleanup_and_log(start_time)?;

    Ok(())
}
