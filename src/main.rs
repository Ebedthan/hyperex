// Copyright 2021-2022 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate anyhow;
extern crate clap;
extern crate log;
extern crate sysinfo;

mod app;
mod utils;

use anyhow::Result;
use bio::io::fasta;
use clap::crate_version;
use log::{error, info, warn};
use sysinfo::{System, SystemExt};

use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::process;
use std::time::Instant;

fn main() -> Result<()> {
    // SET/GET REQUIREMENTS -------------------------------------------------
    // Setup sysinfo
    let mut sys = System::new_all();
    sys.refresh_all();

    // Starting up the Walltime chrono
    let startime = Instant::now();
    let stderr = std::io::stderr();
    let mut ehandle = stderr.lock();

    // Get command-line arguments (see app.rs)
    let matches = app::build_app().get_matches_from(env::args_os());

    // is --quiet option specified by the user?
    let quiet = matches.is_present("quiet");
    utils::setup_logging(quiet)?; // Settting up logging

    // Reading input data
    // This can be a piped data or a filename
    // So we match the value to '-' or some other value and read it
    let infile = match matches.value_of("FILE") {
        // Read from file if passed arg is not '-', otherwise read from stdin
        Some(value) => {
            if value == "-" {
                let mut writer = fasta::Writer::to_file("infile.fa")?;
                let mut records = fasta::Reader::new(io::stdin()).records();
                while let Some(Ok(record)) = records.next() {
                    writer.write_record(&record)?;
                }
                "infile.fa"
            } else {
                value
            }
        }
        // Read from STDIN
        None => {
            let mut writer = fasta::Writer::to_file("infile.fa")?;
            let mut records = fasta::Reader::new(io::stdin()).records();
            while let Some(Ok(record)) = records.next() {
                writer.write_record(&record)?;
            }
            "infile.fa"
        }
    };

    // Check that the supplied file exists
    if infile != "infile.fa" {
        match Path::new(infile).exists() {
            true => (),
            false => {
                writeln!(ehandle, "error: No such file or directory. Is the path correct? Do you have permission to read the file?")?;
                process::exit(1);
            }
        }
    }

    // Read prefix for output files
    let prefix = matches.value_of("prefix").unwrap();
    let force = matches.is_present("force");
    if !force {
        if Path::new(format!("{}.fa", prefix).as_str()).exists()
            || Path::new(format!("{}.gff", prefix).as_str()).exists()
        {
            writeln!(std::io::stderr(), "error: file already exists. Please change it using --prefix option or use --force to overwrite it")?;
            process::exit(1);
        }
    } else if force {
        fs::remove_file(format!("{}.fa", prefix).as_str())?;
        fs::remove_file(format!("{}.gff", prefix).as_str())?;
    }

    // Get primers from command-line as a list of primer can be specified
    let mut primers: Vec<Vec<String>> = Vec::new();
    let all = vec![
        "v1v2", "v1v3", "v1v9", "v3v4", "v3v5", "v4", "v4v5", "v5v7", "v6v9",
        "v7v9",
    ];

    // Case the user go for -f and -r options
    if matches.is_present("forward_primer") && primers.is_empty() {
        // Read supplied forward and reverse primers
        let first: Vec<String> = matches.values_of_t("forward_primer")?;
        let second: Vec<String> = matches.values_of_t("reverse_primer")?;

        // Primers should be in pairs!
        if (first.len() % 2) == 0 && first.len() != second.len() {
            writeln!(ehandle,
                "Supplied forward and reverse primers are not multiple of 2. Please check specified primers"
            )?;
            process::exit(1);
        }

        // Combine both Vec<String> into Vec<Vec<String>>
        primers = utils::combine_vec(first, second);

    // Case user goes for --region option
    } else if matches.is_present("region") {
        // Get supplied region names which can be multiple
        let regions: Vec<String> = matches.values_of_t("region")?;

        // Check if its a file that have been supplied or region name
        if Path::new(&regions[0]).is_file() {
            // We will consider in this case that the region name is a file
            primers = utils::file_to_vec(&regions[0]).unwrap();
        // Check that the region name is supported
        } else if regions.iter().all(|x| all.contains(&x.as_str())) {
            primers = regions
                .iter()
                .map(|x| utils::region_to_primer(x).unwrap())
                .collect::<Vec<_>>();
        } else {
            writeln!(ehandle, "Supplied region is not a correct file name nor a supported region name")?;
            process::exit(1);
        }
    } else {
        // Case when no region or primer is supplied, all the built-in regions are
        // extracted
        primers = all
            .iter()
            .map(|x| utils::region_to_primer(x).unwrap())
            .collect::<Vec<_>>();
    }

    let mismatch: u8 = matches.value_of_t("mismatch")?;

    // STARTING CORE PROGRAM ------------------------------------------------
    info!("This is hyperex v{}", crate_version!());
    info!("Written by Anicet Ebou");
    info!("Available at https://github.com/Ebedthan/hyperex.git");
    info!("Localtime is {}", chrono::Local::now().format("%H:%M:%S"));
    info!(
        "You are {}",
        sys.host_name()
            .unwrap_or_else(|| String::from("not teeling me who you are"))
    );
    info!(
        "Operating system is {}",
        sys.name().unwrap_or_else(|| String::from("not known"))
    );

    if mismatch != 0 {
        warn!(
            "You have allowed {} mismatch in the primer sequence",
            mismatch
        );
    }

    if force {
        warn!("Overwriting {}.fa and {}.gff files", prefix, prefix);
    }

    // Check that required number of mismatch is not greater than
    // the length of the longest primer
    let cp_primers = primers.clone();
    let longest_primer_length =
        cp_primers.into_iter().flatten().map(|x| x.len()).max();

    match longest_primer_length {
        Some(l) => {
            if mismatch as usize > l {
                error!("Supplied mismatch is greater that length of primer");
                error!("Aborting...");
                process::exit(1);
            }
        }
        None => {
            error!("No primer sequence detected");
            error!("Aborting...");
            process::exit(1);
        }
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
        "Walltime: {}h:{}m:{}s.{}ms",
        hours, minutes, seconds, milliseconds
    );
    info!("Memory used: {} KB", sys.used_memory());
    info!("Enjoy. Share. Come back again!");

    Ok(())
}
