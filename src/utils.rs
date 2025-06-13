// Copyright 2021-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use super::cli;

use anyhow::{anyhow, Context, Result};
use bio::io::fasta;
use bio::pattern_matching::myers::MyersBuilder;
use fern::colors::{Color, ColoredLevelConfig};
use log::{error, info, warn};
use phf::phf_map;
use std::{
    fs::{self, File, OpenOptions},
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

// Constants for configuration
const TEMP_FILE: &str = "infile.fa";
const SUPPORTED_REGIONS: &[&str] = &[
    "v1v2", "v1v3", "v1v9", "v3v4", "v3v5", "v4", "v4v5", "v5v7", "v6v9", "v7v9",
];
const MIN_SEQ_LENGTH: usize = 1500;

// Primer data using phf maps for fast lookup
static PRIMER_TO_REGION: phf::Map<&'static str, &'static str> = phf_map! {
    "AGAGTTTGATCMTGGCTCAG" => "v1",
    "CCTACGGGNGGCWGCAG" => "v3",
    "GTGCCAGCMGCCGCGGTAA" => "v4",
    "GTGYCAGCMGCCGCGGTAA" => "v4",
    "AACMGGATTAGATACCCKG" => "v5",
    "TAAAACTYAAAKGAATTGACGGGG" => "v6",
    "YAACGAGCGCAACCC" => "v7",
    "ACTGCTGCSYCCCGTAGGAGTCT" => "v2",
    "ATTACCGCGGCTGCTGG" => "v3",
    "GACTACHVGGGTATCTAATCC" => "v4",
    "CCGTCAATTYMTTTRAGT" => "v5",
    "GGACTACHVGGGTWTCTAAT" => "v4",
    "CCCCGYCAATTCMTTTRAGT" => "v5",
    "ACGTCATCCCCACCTTCC" => "v7",
    "TACGGYTACCTTGTTAYGACTT" => "v9"
};

static FORWARD_PRIMERS: phf::Map<&'static str, &'static str> = phf_map! {
    "27F" => "AGAGTTTGATCMTGGCTCAG",
    "341F" => "CCTACGGGNGGCWGCAG",
    "515F" => "GTGCCAGCMGCCGCGGTAA",
    "515F-Y" => "GTGYCAGCMGCCGCGGTAA",
    "799F" => "AACMGGATTAGATACCCKG",
    "928F" => "TAAAACTYAAAKGAATTGACGGGG",
    "1100F" => "YAACGAGCGCAACCC",
};

static REVERSE_PRIMERS: phf::Map<&'static str, &'static str> = phf_map! {
    "336R" => "ACTGCTGCSYCCCGTAGGAGTCT",
    "534R" => "ATTACCGCGGCTGCTGG",
    "805R" => "GACTACHVGGGTATCTAATCC",
    "926Rb" => "CCGTCAATTYMTTTRAGT",
    "806R" => "GGACTACHVGGGTWTCTAAT",
    "909-928R" => "CCCCGYCAATTCMTTTRAGT",
    "1193R" => "ACGTCATCCCCACCTTCC",
    "1492Rmod" => "TACGGYTACCTTGTTAYGACTT",
};

// Improved logging setup with better color configuration
pub fn setup_logging(quiet: bool) -> Result<(), fern::InitError> {
    let colors = ColoredLevelConfig::new()
        .debug(Color::Blue)
        .info(Color::Green)
        .warn(Color::Yellow)
        .error(Color::Red);

    let level_filter = if quiet {
        log::LevelFilter::Warn
    } else {
        log::LevelFilter::Debug
    };

    fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "[{}][{}] {}",
                chrono::Local::now().format("%H:%M:%S"),
                colors.color(record.level()),
                message
            ))
        })
        .level(level_filter)
        .chain(io::stdout())
        .chain(
            fern::Dispatch::new()
                .format(|out, message, record| {
                    out.finish(format_args!(
                        "{}[{}][{}] {}",
                        chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                        record.target(),
                        record.level(),
                        message
                    ))
                })
                .chain(fern::log_file("hyperex.log")?),
        )
        .apply()?;

    Ok(())
}

// Region to primer with compile-time checks
pub fn region_to_primer(region: &str) -> Result<Vec<String>> {
    let (f_key, r_key) = match region {
        "v1v2" => ("27F", "336R"),
        "v1v3" => ("27F", "534R"),
        "v1v9" => ("27F", "1492Rmod"),
        "v3v4" => ("341F", "805R"),
        "v3v5" => ("341F", "926Rb"),
        "v4" => ("515F", "806R"),
        "v4v5" => ("515F-Y", "909-928R"),
        "v5v7" => ("799F", "1193R"),
        "v6v9" => ("928F", "1492Rmod"),
        "v7v9" => ("1100F", "1492Rmod"),
        _ => return Ok(Vec::new()),
    };

    Ok(vec![
        FORWARD_PRIMERS[f_key].to_string(),
        REVERSE_PRIMERS[r_key].to_string(),
    ])
}

// File parsing with better error handling
pub fn file_to_vec(filename: &str) -> Result<Vec<Vec<String>>> {
    fs::read_to_string(filename)?
        .lines()
        .map(|line| {
            if line.contains(',') {
                Ok(line.split(',').map(String::from).collect())
            } else {
                Err(anyhow!("Primer file must be comma-separated"))
            }
        })
        .collect()
}

// Vector combination
pub fn combine_vec(first: Vec<String>, second: Vec<String>) -> Vec<Vec<String>> {
    first
        .into_iter()
        .zip(second)
        .map(|(f, r)| vec![f, r])
        .collect()
}

// File reading
fn read_file(filename: &str) -> Result<(Box<dyn io::Read>, niffler::compression::Format)> {
    let file =
        File::open(filename).with_context(|| format!("Failed to open file: {}", filename))?;
    let reader = Box::new(io::BufReader::new(file));
    niffler::get_reader(reader)
        .with_context(|| format!("Failed to determine compression for: {}", filename))
}

// Primer to region mapping
fn primers_to_region(primers: &[String]) -> String {
    let first = PRIMER_TO_REGION.get(&primers[0]).unwrap_or(&"");
    let second = PRIMER_TO_REGION.get(&primers[1]).unwrap_or(&"");

    if *first == "v4" && *second == "v4" {
        first.to_string()
    } else {
        format!("{}{}", first, second)
    }
}

// Efficient complement generation
fn to_complement(primer: &str, alphabet: &str) -> String {
    primer
        .chars()
        .map(|c| match (c, alphabet) {
            ('A', "dna") => 'T',
            ('T', "dna") => 'A',
            ('C', "dna") => 'G',
            ('G', "dna") => 'C',
            ('A', "rna") => 'U',
            ('U', "rna") => 'A',
            ('C', "rna") => 'G',
            ('G', "rna") => 'C',
            // IUPAC ambiguity codes
            ('R', _) => 'Y',
            ('Y', _) => 'R',
            ('K', _) => 'M',
            ('M', _) => 'K',
            ('B', _) => 'V',
            ('V', _) => 'B',
            ('D', _) => 'H',
            ('H', _) => 'D',
            // Default cases
            _ => c,
        })
        .collect()
}

// Reverse complement
fn to_reverse_complement(primer: &str, alphabet: &str) -> String {
    to_complement(primer, alphabet).chars().rev().collect()
}

#[derive(Debug, PartialEq)]
pub enum Alphabet {
    Dna,
    Rna,
}

pub fn sequence_type(sequence: &str) -> Option<Alphabet> {
    let has_u = sequence.contains('U');
    let has_t = sequence.contains('T');

    if has_u && !has_t {
        Some(Alphabet::Rna)
    } else if has_t && !has_u {
        Some(Alphabet::Dna)
    } else if !has_u && !has_t {
        // Could be either, check other IUPAC codes
        if sequence.chars().all(|c| "ACGTRYSWKMBDHVN".contains(c)) {
            Some(Alphabet::Dna)
        } else {
            None
        }
    } else {
        None
    }
}

// Core processing function
pub fn get_hypervar_regions(
    file: &str,
    primers: Vec<Vec<String>>,
    prefix: &str,
    mismatch: u8,
) -> Result<()> {
    let (reader, _) = read_file(file)?;
    let records = fasta::Reader::new(reader).records();

    // Initialize output writers
    let mut fasta_writer = fasta::Writer::to_file(format!("{}.fa", prefix))?;
    let gff_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(format!("{}.gff", prefix))?;
    let mut gff_writer = BufWriter::new(gff_file);
    gff_writer.write_all(b"##gff-version 3\n")?;

    // Configure Myers builder with IUPAC ambiguities
    // Using Vec<u8> to handle variable length ambigs
    let ambigs = vec![
        (b'M', b"AC".to_vec()),
        (b'R', b"AG".to_vec()),
        (b'W', b"AT".to_vec()),
        (b'S', b"CG".to_vec()),
        (b'Y', b"CT".to_vec()),
        (b'K', b"GT".to_vec()),
        (b'V', b"ACGMRS".to_vec()),
        (b'H', b"ACTMWY".to_vec()),
        (b'D', b"AGTRWK".to_vec()),
        (b'B', b"CGTSYK".to_vec()),
        (b'N', b"ACGTMRWSYKVHDB".to_vec()),
    ];

    let mut builder = MyersBuilder::new();
    for (base, equivalents) in ambigs {
        builder.ambig(base, equivalents.as_slice());
    }

    for record in records.flatten() {
        let seq = record.seq();
        let seq_str = std::str::from_utf8(seq)?;
        let alphabet = match sequence_type(seq_str) {
            Some(Alphabet::Dna) => "dna",
            Some(Alphabet::Rna) => "rna",
            None => {
                error!("Unrecognized sequence type in record: {}", record.id());
                continue;
            }
        };

        if seq.len() <= MIN_SEQ_LENGTH {
            warn!(
                "Sequence {} is short ({} bp). Some regions may not be found.",
                record.id(),
                seq.len()
            );
        }

        for primer_pair in &primers {
            let region = primers_to_region(primer_pair);
            let rev_comp = to_reverse_complement(&primer_pair[1], alphabet);

            let mut forward_myers = builder.build_64(primer_pair[0].as_bytes());
            let mut reverse_myers = builder.build_64(rev_comp.as_bytes());

            // Find all matches for both primers
            let forward_matches: Vec<_> = forward_myers.find_all_lazy(seq, mismatch).collect();

            let reverse_matches: Vec<_> = reverse_myers.find_all_lazy(seq, mismatch).collect();

            // Find best matches (lowest distance)
            let forward_best = forward_matches.iter().min_by_key(|&&(_, dist)| dist);
            let reverse_best = reverse_matches.iter().min_by_key(|&&(_, dist)| dist);

            match (forward_best, reverse_best) {
                (Some(&(f_end, _)), Some(&(r_end, _))) => {
                    // Calculate positions - Myers returns end positions
                    let f_start = f_end - primer_pair[0].len() + 1;
                    //let r_start = r_end - primer_pair[1].len() + 1;

                    // Write FASTA record
                    let desc = if region.is_empty() {
                        format!("forward={} reverse={}", primer_pair[0], primer_pair[1])
                    } else {
                        format!(
                            "region={} forward={} reverse={}",
                            region, primer_pair[0], primer_pair[1]
                        )
                    };

                    fasta_writer.write_record(&fasta::Record::with_attrs(
                        record.id(),
                        Some(&desc),
                        &seq[f_start..r_end + 1], // +1 to include the last base
                    ))?;

                    // Write GFF record (1-based coordinates)
                    writeln!(
                        gff_writer,
                        "{}\thyperex\tregion\t{}\t{}\t.\t.\t.\tNote=Hypervariable region {}",
                        record.id(),
                        f_start + 1,
                        r_end + 1,
                        region
                    )?;
                }
                (None, Some(_)) => warn!(
                    "Forward primer {} not found for region {} in sequence {}",
                    primer_pair[0],
                    region,
                    record.id()
                ),
                (Some(_), None) => warn!(
                    "Reverse primer {} not found for region {} in sequence {}",
                    primer_pair[1],
                    region,
                    record.id()
                ),
                (None, None) => warn!(
                    "Both primers not found for region {} in sequence {}",
                    region,
                    record.id()
                ),
            }
        }
    }

    Ok(())
}

pub fn handle_input(file_arg: &Option<String>, ehandle: &mut io::StderrLock) -> Result<String> {
    match file_arg {
        Some(value) if value == "-" => read_stdin_to_temp_file(),
        Some(value) => {
            if !Path::new(value).exists() {
                writeln!(ehandle, "error: No such file or directory. Is the path correct? Do you have permission to read the file?")?;
                std::process::exit(1);
            }
            Ok(value.clone())
        }
        None => read_stdin_to_temp_file(),
    }
}

pub fn read_stdin_to_temp_file() -> Result<String> {
    let mut writer = fasta::Writer::to_file(TEMP_FILE)?;
    let reader = fasta::Reader::new(io::stdin());

    for record in reader.records().flatten() {
        writer.write_record(&record)?;
    }

    Ok(TEMP_FILE.to_string())
}

pub fn handle_output_files(prefix: &str, force: bool, ehandle: &mut io::StderrLock) -> Result<()> {
    let output_files = [
        PathBuf::from(format!("{}.fa", prefix)),
        PathBuf::from(format!("{}.gff", prefix)),
    ];

    if !force {
        for file in &output_files {
            if file.exists() {
                writeln!(
                    ehandle,
                    "error: file {} already exists. Please change it using --prefix option or use --force to overwrite it",
                    file.display()
                )?;
                std::process::exit(1);
            }
        }
    } else {
        for file in &output_files {
            if file.exists() {
                fs::remove_file(file)?;
            }
        }
    }

    Ok(())
}

pub fn process_primers(cli: &cli::Args, ehandle: &mut io::StderrLock) -> Result<Vec<Vec<String>>> {
    match (&cli.forward, &cli.reverse, &cli.region) {
        (Some(forward), Some(reverse), None) => {
            if forward.len() != reverse.len() {
                writeln!(
                    ehandle,
                    "Supplied forward and reverse primers must have the same number of elements"
                )?;
                std::process::exit(1);
            }
            Ok(combine_vec(forward.clone(), reverse.clone()))
        }
        (None, None, Some(regions)) => {
            if regions.is_empty() {
                return Ok(SUPPORTED_REGIONS
                    .iter()
                    .map(|x| region_to_primer(x).unwrap())
                    .collect());
            }

            if Path::new(&regions[0].to_string()).is_file() {
                file_to_vec(&regions[0].to_string())
            } else if regions
                .iter()
                .all(|x| SUPPORTED_REGIONS.contains(&x.to_string().as_str()))
            {
                Ok(regions
                    .iter()
                    .map(|x| region_to_primer(&x.to_string()).unwrap())
                    .collect())
            } else {
                writeln!(
                    ehandle,
                    "Supplied region is not a correct file name nor a supported region name"
                )?;
                std::process::exit(1);
            }
        }
        _ => {
            // Default case: use all built-in regions
            Ok(SUPPORTED_REGIONS
                .iter()
                .map(|x| region_to_primer(x).unwrap())
                .collect())
        }
    }
}

pub fn validate_mismatch(primers: &[Vec<String>], mismatch: u8) -> Result<()> {
    let max_primer_len = primers
        .iter()
        .flatten()
        .map(|s| s.len())
        .max()
        .ok_or_else(|| anyhow::anyhow!("No primer sequences detected"))?;

    if mismatch as usize > max_primer_len {
        error!(
            "Supplied mismatch ({}) is greater than length of longest primer ({})",
            mismatch, max_primer_len
        );
        error!("Aborting...");
        std::process::exit(1);
    }

    Ok(())
}

pub fn log_program_info(mismatch: u8, force: bool, prefix: &str) {
    info!("This is hyperex v0.2");
    info!("Written by Anicet Ebou");
    info!("Available at https://github.com/Ebedthan/hyperex.git");
    info!("Localtime is {}", chrono::Local::now().format("%H:%M:%S"));

    if mismatch != 0 {
        warn!(
            "You have allowed {} mismatch in the primer sequence",
            mismatch
        );
    }

    if force {
        warn!("Overwriting {}.fa and {}.gff files", prefix, prefix);
    }
}

pub fn cleanup_and_log(start_time: Instant) -> Result<()> {
    if Path::new(TEMP_FILE).exists() {
        fs::remove_file(TEMP_FILE)?;
    }

    let duration = start_time.elapsed();
    log_duration(duration);
    info!("Enjoy. Share. Come back again!");

    Ok(())
}

pub fn log_duration(duration: Duration) {
    let total_ms = duration.as_millis();
    let (hours, remaining) = (total_ms / 3_600_000, total_ms % 3_600_000);
    let (minutes, remaining) = (remaining / 60_000, remaining % 60_000);
    let (seconds, milliseconds) = (remaining / 1_000, remaining % 1_000);

    info!(
        "Walltime: {}h:{}m:{}s.{}ms",
        hours, minutes, seconds, milliseconds
    );
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_complement_dna() {
        assert_eq!(
            to_complement("ATCGATCGATCGATCGRYKBVDHX", "dna"),
            String::from("TAGCTAGCTAGCTAGCYRMVBHDX")
        );
    }

    #[test]
    fn test_complement_rna() {
        assert_eq!(
            to_complement("AUCGAUCGAUCGAUCGRYKBVDHMXN", "rna"),
            String::from("UAGCUAGCUAGCUAGCYRMVBHDKXN")
        );
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(
            to_reverse_complement("GTGCCAGCMGCCGCGGTAAN", "dna"),
            "NTTACCGCGGCKGCTGGCAC"
        );
    }

    #[test]
    fn test_sequence_type_dna_ok() {
        assert_eq!(sequence_type("ATCGATCGATCG"), Some(Alphabet::Dna));
    }

    #[test]
    fn test_sequence_type_dna_iupac_ok() {
        assert_eq!(sequence_type("ATCGMTGCAATCG"), Some(Alphabet::Dna));
    }

    #[test]
    fn test_sequence_type_rna_ok() {
        assert_eq!(sequence_type("AGCUUUGCA"), Some(Alphabet::Rna));
    }

    #[test]
    fn test_sequence_type_rna_iupac_ok() {
        assert_eq!(sequence_type("GUUUUAACCCAAM"), Some(Alphabet::Rna));
    }

    #[test]
    fn test_sequence_type_err() {
        assert_eq!(sequence_type("ATCXXXRMGU"), None);
    }

    #[test]
    fn test_region_to_primer_ok() {
        assert_eq!(
            region_to_primer("v1v2").unwrap(),
            vec!["AGAGTTTGATCMTGGCTCAG", "ACTGCTGCSYCCCGTAGGAGTCT"]
        );
        assert_eq!(
            region_to_primer("v1v3").unwrap(),
            vec!["AGAGTTTGATCMTGGCTCAG", "ATTACCGCGGCTGCTGG"]
        );
        assert_eq!(
            region_to_primer("v1v9").unwrap(),
            vec!["AGAGTTTGATCMTGGCTCAG", "TACGGYTACCTTGTTAYGACTT"]
        );
        assert_eq!(
            region_to_primer("v3v4").unwrap(),
            vec!["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"]
        );
        assert_eq!(
            region_to_primer("v3v5").unwrap(),
            vec!["CCTACGGGNGGCWGCAG", "CCGTCAATTYMTTTRAGT"]
        );
        assert_eq!(
            region_to_primer("v4").unwrap(),
            vec!["GTGCCAGCMGCCGCGGTAA", "GGACTACHVGGGTWTCTAAT"]
        );
        assert_eq!(
            region_to_primer("v4v5").unwrap(),
            vec!["GTGYCAGCMGCCGCGGTAA", "CCCCGYCAATTCMTTTRAGT"]
        );
        assert_eq!(
            region_to_primer("v5v7").unwrap(),
            vec!["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"]
        );
        assert_eq!(
            region_to_primer("v6v9").unwrap(),
            vec!["TAAAACTYAAAKGAATTGACGGGG", "TACGGYTACCTTGTTAYGACTT"]
        );
        assert_eq!(
            region_to_primer("v7v9").unwrap(),
            vec!["YAACGAGCGCAACCC", "TACGGYTACCTTGTTAYGACTT"]
        );
    }

    #[test]
    fn test_write_fa_ok2() {
        let mut tmpfile = NamedTempFile::new().expect("Cannot create temp file");
        writeln!(tmpfile, ">id_str desc\nATCGCCG").expect("Cannot write to tmp file");

        let mut fa_records = fasta::Reader::from_file(tmpfile)
            .expect("Cannot read file.")
            .records();

        while let Some(Ok(rec)) = fa_records.next() {
            assert_eq!(rec.id(), "id_str");
            assert_eq!(rec.desc(), Some("desc"));
            assert_eq!(rec.seq(), b"ATCGCCG");
        }
    }

    #[test]
    fn test_combine_vec() {
        let first = vec!["ab".to_string(), "cd".to_string(), "ef".to_string()];
        let second = vec!["cd".to_string(), "ef".to_string(), "gh".to_string()];
        assert_eq!(
            combine_vec(first, second),
            vec![
                vec!["ab".to_string(), "cd".to_string()],
                vec!["cd".to_string(), "ef".to_string()],
                vec!["ef".to_string(), "gh".to_string()]
            ]
        );
    }

    #[test]
    fn test_combine_vec_not_ok() {
        let first = vec!["ab".to_string(), "cd".to_string(), "ef".to_string()];
        let second = vec!["ab".to_string()];
        assert_ne!(
            combine_vec(first, second),
            vec![
                vec!["ab".to_string(), "cd".to_string()],
                vec!["cd".to_string(), "ef".to_string()],
                vec!["ef".to_string(), "gh".to_string()]
            ]
        );
    }

    #[test]
    fn test_get_hypervar_regions() {
        assert!(get_hypervar_regions(
            "tests/test.fa.gz",
            vec![vec![
                "AGAGTTTGATCMTGGCTCAG".to_string(),
                "TACGGYTACCTTGTTAYGACTT".to_string()
            ]],
            "hyperex",
            0
        )
        .is_ok());
        fs::remove_file("hyperex.fa").expect("cannot delete file");
        fs::remove_file("hyperex.gff").expect("cannot delete file");
    }

    #[test]
    fn test_setup_logging() {
        assert!(setup_logging(false).is_ok());
    }

    #[test]
    fn test_read_file() {
        let myfile = "tests/test.fa.gz";
        assert!(read_file(myfile).is_ok());
    }

    #[test]
    fn test_file_to_vec() {
        assert_eq!(
            file_to_vec("tests/primers.txt").unwrap(),
            vec![
                vec![
                    "CCTACGGGNGGCWGCAG".to_string(),
                    "ATTACCGCGGCTGCTGG".to_string()
                ],
                vec![
                    "GTGCCAGCMGCCGCGGTAA".to_string(),
                    "GACTACHVGGGTATCTAATCC".to_string()
                ]
            ]
        );
    }

    #[test]
    fn test_file_to_vec_no_ok() {
        assert!(file_to_vec("test.fa").is_err());
    }
}
