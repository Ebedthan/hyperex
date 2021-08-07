// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate anyhow;
extern crate clap;
extern crate log;
extern crate whoami;

mod app;
mod utils;

use anyhow::{Context, Result};
use clap::crate_version;
use log::{info, warn};

use std::env;
use std::fs;
use std::path::Path;

fn main() -> Result<()> {
    // Define command-line arguments ----------------------------------------
    let matches = app::build_app().get_matches_from(env::args_os());

    // Read command-line arguments ------------------------------------------
    let infile = matches.value_of("FILE").unwrap();
    let outfile = matches.value_of("output").unwrap();
    let mut primers: Vec<Vec<&str>> = Vec::new();
    if matches.is_present("forward_primer") && primers.is_empty() {
        primers = vec![vec![
            matches.value_of("forward_primer").unwrap(),
            matches.value_of("reverse_primer").unwrap(),
        ]];
    } else if matches.is_present("region") {
        primers = utils::region_to_primer(matches.value_of("region").unwrap())?;
    } else {
        primers = utils::region_to_primer("all")?;
    }

    let force = matches.is_present("force");
    let verbosity: u64 = matches.occurrences_of("verbose");

    // Set up program logging -----------------------------------------------
    utils::setup_logging(verbosity)?;

    info!("This is hyvrex v{}", crate_version!());
    info!("Written by Anicet Ebou");
    info!("Available at https://github.com/Ebedthan/hyvrex.git");
    info!("Localtime is {}", chrono::Local::now().format("[%H:%M:%S]"));
    info!("You are {}", whoami::username());
    info!("Operating system is {}", whoami::platform());

    // Handling creation of file
    if force {
        warn!("Reusing file {}", outfile);
        fs::remove_file(Path::new(outfile)).with_context(|| format!("Could not remove file '{}'. Do you have the appropriate permission?", outfile))?;
    }

    utils::process_fa(infile, primers, outfile)?;

    Ok(())
}
