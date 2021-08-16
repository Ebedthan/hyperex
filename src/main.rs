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

use anyhow::Result;
use bio::io::fasta;
use clap::crate_version;
use log::{error, info, warn};

use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::process;
use std::time::Instant;

// TODO: make --region to accept primer file like forward_primer\treverse_primer\n
fn main() -> Result<()> {
    // SET/GET REQUIREMENTS -------------------------------------------------
    // Starting up the Walltime chrono
    let startime = Instant::now();

    // Get command-line arguments (see app.rs)
    let matches = app::build_app().get_matches_from(env::args_os());

    // is --quiet option specified by the user?
    let quiet = matches.is_present("quiet");
    utils::setup_logging(quiet)?; // Settting up logging

    // Reading input data
    // This can be a piped data or a filename
    // So we match the value to '-' or some other value and read it
    let infile = match matches.value_of("FILE") {
        Some(val) => {
            if val != "-" {
                val
            } else {
                let mut writer = fasta::Writer::to_file("infile.fa")?;
                let mut records = fasta::Reader::new(io::stdin()).records();
                while let Some(Ok(record)) = records.next() {
                    writer.write_record(&record)?;
                }

                "infile.fa"
            }
        }
        // empty value that should never been reached as clap will always check
        // that infile is supplied
        _ => "",
    };

    // Check that the supplied file exists
    if infile != "infile.fa" {
        match Path::new(infile).exists() {
            true => (),
            false => {
                writeln!(std::io::stderr(), "error: No such file or directory. Is the path correct? Do you have permission to read the file?")?;
                process::exit(1);
            }
        }
    }

    // Read prefix for output files
    let prefix = matches.value_of("prefix").unwrap();

    // Get primers from command-line as a list of primer can be specified
    let mut primers: Vec<Vec<&str>> = Vec::new();
    if matches.is_present("forward_primer") && primers.is_empty() {
        // Read supplied forward and reverse primers
        let first = matches
            .values_of("forward_primer")
            .unwrap()
            .collect::<Vec<_>>();
        let second = matches
            .values_of("reverse_primer")
            .unwrap()
            .collect::<Vec<_>>();

        if (first.len() % 2) == 0 && first.len() != second.len() {
            error!(
                "Supplied forward and reverse primers are not multiple of 2. Please check specified primers"
            );
            process::exit(1);
        }
        // Combine both vec into one big vec
        primers = utils::combine_vec(first, second);
    } else if matches.is_present("region") {
        // Get supplied region names
        let regions = matches.values_of("region").unwrap().collect::<Vec<_>>();
        primers = regions
            .iter()
            .map(|x| utils::region_to_primer(x).unwrap())
            .collect::<Vec<_>>();
    } else {
        // Case when no region or primer is supplied, all the built-in regions are
        // extracted
        let all = vec![
            "v1v2", "v1v3", "v1v9", "v3v4", "v3v5", "v4", "v4v5", "v5v7",
            "v6v9", "v7v9",
        ];
        primers = all
            .iter()
            .map(|x| utils::region_to_primer(x).unwrap())
            .collect::<Vec<_>>();
    }
    let mis = matches.value_of("mismatch").unwrap().to_string();
    let mismatch = mis.parse::<u8>()?;

    // STARTING CORE PROGRAM ------------------------------------------------
    info!("This is hyperex v{}", crate_version!());
    info!("Written by Anicet Ebou");
    info!("Available at https://github.com/Ebedthan/hyperex.git");
    info!("Localtime is {}", chrono::Local::now().format("%H:%M:%S"));
    info!("You are {}", whoami::username());
    info!("Operating system is {}", whoami::platform());

    if mismatch != 0 {
        warn!(
            "You have allowed {} mismatch in the primer sequence",
            mismatch
        );
    }

    utils::get_hypervar_regions(infile, primers, prefix, mismatch)?;
    info!("Done getting hypervariable regions");

    // FINISHING ------------------------------------------------------------
    // Cleaning around
    if Path::new("infile.fa").exists() {
        fs::remove_file("infile.fa")?;
    }
    let duration = startime.elapsed();
    let y = 60 * 60 * 1000;
    let hours = duration.as_millis() / y;
    let minutes = (duration.as_millis() - (hours * y)) / (y / 60);
    let seconds =
        (duration.as_millis() - (hours * y) - (minutes * (y / 60))) / 1000;
    let milliseconds = duration.as_millis()
        - (hours * y)
        - (minutes * (y / 60))
        - (seconds * 1000);

    info!(
        "{}",
        format!(
            "Walltime: {}h:{}m:{}s.{}ms",
            hours, minutes, seconds, milliseconds
        )
    );
    info!("Enjoy. Share. Come back again!");

    Ok(())
}
