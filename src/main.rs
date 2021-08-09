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
use std::io;
use std::path::Path;
use std::time::Instant;

// TODO: make --region to accept primer file like forward_primer\treverse_primer\n
fn main() -> Result<()> {
    // Starting up the chrono
    let startime = Instant::now();

    // Define command-line arguments ----------------------------------------
    let matches = app::build_app().get_matches_from(env::args_os());

    // Read command-line arguments ------------------------------------------
    // Either read the input through stdin or use the supplied value as filename
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
        // dump value that should never been reached as clap will always check
        // that infile is supplied
        _ => "no",
    };

    // Check that the supplied file is valid
    if infile != "infile.fa" {
        utils::is_fasta(infile).unwrap();
    }

    let prefix = matches.value_of("prefix").unwrap();

    let quiet = matches.is_present("quiet");

    // Set up program logging -----------------------------------------------
    utils::setup_logging()?;

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

        if first.len() != second.len() {
            error!(
                "Supplied forward and reverse primers are not multiple of 2"
            );
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
        // Case no region or primer is supplied, all the built-in regions are
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

    if !quiet {
        info!("This is hyvrex v{}", crate_version!());
        info!("Written by Anicet Ebou");
        info!("Available at https://github.com/Ebedthan/hyvrex.git");
        info!("Localtime is {}", chrono::Local::now().format("%H:%M:%S"));
        info!("You are {}", whoami::username());
        info!("Operating system is {}", whoami::platform());
    }
    if mismatch != 0 {
        warn!(
            "You have allowed {} mismatch in the primer sequence",
            mismatch
        );
    }
    utils::process_fa(infile, primers, prefix, mismatch)?;

    // Finishing
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

    if !quiet {
        info!(
            "{}",
            format!(
                "Walltime: {}h:{}m:{}s.{}ms",
                hours, minutes, seconds, milliseconds
            )
        );
        info!("Done getting hypervariable regions");
        info!("Enjoy. Share. Come back again!");
    }

    Ok(())
}
