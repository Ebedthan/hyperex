// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate chrono;
extern crate fern;
extern crate log;
extern crate niffler;
extern crate phf;
extern crate regex;

use anyhow::{anyhow, Context, Result};
use bio::io::fasta;
use fern::colors::ColoredLevelConfig;
use log::warn;
use phf::phf_map;
use regex::Regex;

use std::fs::File;
use std::fs::OpenOptions;
use std::io;

pub fn setup_logging(verbosity: u64) -> Result<(), fern::InitError> {
    let colors = ColoredLevelConfig::default();

    let mut base_config = fern::Dispatch::new();

    base_config = match verbosity {
        0 => base_config
            .level(log::LevelFilter::Info)
            .level_for("overly-verbose-target", log::LevelFilter::Warn),
        1 => base_config
            .level(log::LevelFilter::Debug)
            .level_for("overly-verbose-target", log::LevelFilter::Info),
        2 => base_config.level(log::LevelFilter::Debug),
        _3_or_more => base_config.level(log::LevelFilter::Trace),
    };

    // Separate file config so we can include year, month and day in file logs
    let file_config = fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                record.target(),
                record.level(),
                message
            ))
        })
        .chain(fern::log_file("hyvrex.log")?);

    let stdout_config = fern::Dispatch::new()
        .format(move |out, message, record| {
            // special format for debug messages coming from our own crate.
            if record.level() > log::LevelFilter::Info
                && record.target() == "hyvrex"
            {
                out.finish(format_args!(
                    "---\nDEBUG: {}: {}\n---",
                    chrono::Local::now().format("%H:%M:%S"),
                    message
                ))
            } else {
                out.finish(format_args!(
                    "[{}][{}] {}",
                    chrono::Local::now().format("%H:%M:%S"),
                    colors.color(record.level()),
                    message
                ))
            }
        })
        .chain(io::stdout());

    base_config
        .chain(file_config)
        .chain(stdout_config)
        .apply()?;

    Ok(())
}

// Primers data
static PRIMER_TO_REGION: phf::Map<&'static str, &'static str> = phf_map! {
    "AGAGTTTGATCMTGGCTCAG" => "v1",
    "CCTACGGGNGGCWGCAG" => "v3",
    "GTGCCAGCMGCCGCGGTAA" => "v4",
    "GTGYCAGCMGCCGCGGTAA" => "v4",
    "AACMGGATTAGATACCCKG" => "v5",
    "TAAAACTYAAAKGAATTGACGGGG" => "v6",
    "YAACGAGCGCAACCC" => "v7",
    "CYIACTGCTGCCTCCCGTAG" => "v2",
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
    "337R" => "CYIACTGCTGCCTCCCGTAG",
    "534R" => "ATTACCGCGGCTGCTGG",
    "805R" => "GACTACHVGGGTATCTAATCC",
    "926Rb" => "CCGTCAATTYMTTTRAGT",
    "806R" => "GGACTACHVGGGTWTCTAAT",
    "909-928R" => "CCCCGYCAATTCMTTTRAGT",
    "1193R" => "ACGTCATCCCCACCTTCC",
    "1492Rmod" => "TACGGYTACCTTGTTAYGACTT",
};

pub fn region_to_primer(region: &str) -> Result<Vec<Vec<&str>>> {
    match region {
        "v1v2" => {
            Ok(vec![vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["337R"]]])
        }
        "v1v3" => {
            Ok(vec![vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["534R"]]])
        }
        "v1v9" => Ok(vec![vec![
            FORWARD_PRIMERS["27F"],
            REVERSE_PRIMERS["1492Rmod"],
        ]]),
        "v3v4" => {
            Ok(vec![vec![FORWARD_PRIMERS["341F"], REVERSE_PRIMERS["805R"]]])
        }
        "v3v5" => Ok(vec![vec![
            FORWARD_PRIMERS["341F"],
            REVERSE_PRIMERS["926Rb"],
        ]]),
        "v4" => {
            Ok(vec![vec![FORWARD_PRIMERS["515F"], REVERSE_PRIMERS["806R"]]])
        }
        "v4v5" => Ok(vec![vec![
            FORWARD_PRIMERS["515F-Y"],
            REVERSE_PRIMERS["909-928R"],
        ]]),
        "v5v7" => Ok(vec![vec![
            FORWARD_PRIMERS["799F"],
            REVERSE_PRIMERS["1193R"],
        ]]),
        "v6v9" => Ok(vec![vec![
            FORWARD_PRIMERS["928F"],
            REVERSE_PRIMERS["1492Rmod"],
        ]]),
        "v7v9" => Ok(vec![vec![
            FORWARD_PRIMERS["1100F"],
            REVERSE_PRIMERS["1492Rmod"],
        ]]),
        _ => Ok(vec![
            vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["337R"]],
            vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["534R"]],
            vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["1492Rmod"]],
            vec![FORWARD_PRIMERS["341F"], REVERSE_PRIMERS["805R"]],
            vec![FORWARD_PRIMERS["341F"], REVERSE_PRIMERS["926Rb"]],
            vec![FORWARD_PRIMERS["515F"], REVERSE_PRIMERS["806R"]],
            vec![FORWARD_PRIMERS["515F-Y"], REVERSE_PRIMERS["909-928R"]],
            vec![FORWARD_PRIMERS["799F"], REVERSE_PRIMERS["1193R"]],
            vec![FORWARD_PRIMERS["928F"], REVERSE_PRIMERS["1492Rmod"]],
            vec![FORWARD_PRIMERS["1100F"], REVERSE_PRIMERS["1492Rmod"]],
        ]),
    }
}

// read_file function -------------------------------------------------------

/// Get reader and compression format of file
///
/// # Example
/// ```rust
/// # use std::path::Path;
///
/// let path = Path::new("path/to/file");
/// let (reader, compression) = read_file(&path);
/// ```
///
fn read_file(
    filename: &str,
) -> Result<(Box<dyn io::Read>, niffler::compression::Format)> {
    let raw_in = Box::new(io::BufReader::new(
        File::open(filename).with_context(|| "Cannot read a file")?,
    ));

    niffler::get_reader(raw_in).with_context(|| {
        anyhow!("Could not detect compression of file '{}'", filename)
    })
}

// write_to_fa function -----------------------------------------------------

/// Write to provided data to a fasta file in append mode
///
/// # Example
/// ```rust
///
/// # use bio::io::fasta;
/// let filename = "myfile.fa";
/// let record = fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
/// write_to_fa(filename, &record);
/// ```
///
fn write_fa(file: &std::fs::File, record: fasta::Record) -> Result<()> {
    let handle = io::BufWriter::new(file);
    let mut writer = fasta::Writer::new(handle);
    let write_result = writer.write_record(&record)?;

    Ok(write_result)
}

fn primers_to_region(primers: Vec<&str>) -> String {
    let mut first_part = "";
    let mut second_part = "";

    if PRIMER_TO_REGION.contains_key(primers[0]) {
        first_part = PRIMER_TO_REGION[primers[0]];
    }

    if PRIMER_TO_REGION.contains_key(primers[1]) {
        second_part = PRIMER_TO_REGION[primers[1]];
    }

    if first_part == "v4" && second_part == "v4" {
        first_part.to_string()
    } else {
        format!("{}{}", first_part, second_part)
    }
}

fn from_iupac_to_regex(primer: &str) -> String {
    let clean_regex: String = primer
        .chars()
        .map(|x| match x {
            'R' => "[AG]".to_string(),
            'Y' => "[CT]".to_string(),
            'S' => "[GC]".to_string(),
            'W' => "[AT]".to_string(),
            'K' => "[GT]".to_string(),
            'M' => "[AC]".to_string(),
            'B' => "[CGT]".to_string(),
            'D' => "[AGT]".to_string(),
            'H' => "[ACT]".to_string(),
            'V' => "[ACG]".to_string(),
            'N' => ".".to_string(),
            _ => String::from(x),
        })
        .collect();

    clean_regex
}

fn to_complement(primer: &str) -> String {
    let complement: String = primer
        .chars()
        .map(|x| match x {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'R' => 'Y',
            'Y' => 'R',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            _ => x,
        })
        .collect();

    complement
}

fn to_reverse_complement(primer: &str) -> String {
    let complement = to_complement(primer);
    let reverse_complement = complement.chars().rev().collect();

    reverse_complement
}

pub fn process_fa(
    file: &str,
    primers: Vec<Vec<&str>>,
    outfile: &str,
) -> Result<()> {
    let (reader, mut _compression) =
        read_file(file).with_context(|| "Cannot read file")?;

    let mut records = fasta::Reader::new(reader).records();

    let out = OpenOptions::new().create(true).append(true).open(outfile)?;

    while let Some(Ok(record)) = records.next() {
        let seq = record.seq();
        if seq.len() <= 1400 {
            warn!("Sequence length is less than 1400 bp. We may not be able to find some regions");
        }
        for primer_pair in primers.iter() {
            let forward_re = Regex::new(&from_iupac_to_regex(primer_pair[0]))?;
            let reverse_re = Regex::new(&from_iupac_to_regex(
                &to_reverse_complement(primer_pair[1]),
            ))?;

            let region = primers_to_region(primer_pair.to_vec());

            let forward_match = forward_re.find(std::str::from_utf8(seq)?);
            let reverse_match = reverse_re.find(std::str::from_utf8(seq)?);

            match forward_match {
                Some(fmat) => {
                    match reverse_match {
                        Some(rmat) => {
                            if !region.is_empty() {
                                write_fa(
                                    &out,
                                    fasta::Record::with_attrs(
                                        record.id(),
                                        Some(
                                            format!(
                                                "region={} forward={} reverse={}",
                                                region, primer_pair[0], primer_pair[1]
                                            )
                                            .as_str(),
                                        ),
                                        &seq[fmat.start()..rmat.end()],
                                    ),
                                )?;
                            } else {
                                write_fa(
                                    &out,
                                    fasta::Record::with_attrs(
                                        record.id(),
                                        Some(
                                            format!(
                                                "forward={} reverse={}",
                                                primer_pair[0], primer_pair[1]
                                            )
                                            .as_str(),
                                        ),
                                        &seq[fmat.start()..rmat.end()],
                                    ),
                                )?;
                            }
                        },
                        None => {
                            warn!("Region {} not found because primer {} was not found in the sequence", region, primer_pair[1]);
                        }
                    }
                },
                None => {
                    match reverse_match {
                        Some(_) => {
                            warn!("Region {} not found because primer {} was not found in the sequence", region, primer_pair[0]);
                        },
                        None => warn!("Region {} not found because primers {}, {} was not found in the sequence", region, primer_pair[0], primer_pair[1])
                    }
                }
            }
        }
    }

    Ok(())
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primers_to_region_ok() {
        assert_eq!(
            primers_to_region(vec!["CCTACGGGNGGCWGCAG", "GTGCCAGCMGCCGCGGTAA"]),
            "v3v4".to_string()
        );
    }

    #[test]
    fn test_primers_to_region_ok2() {
        assert_eq!(
            primers_to_region(vec![
                "GTGCCAGCMGCCGCGGTAA",
                "GTGCCAGCMGCCGCGGTAA"
            ]),
            "v4".to_string()
        );
    }

    #[test]
    fn test_primers_to_region_empty() {
        assert_eq!(primers_to_region(vec!["ZZZZZ", "AAAAAA"]), "".to_string());
    }

    #[test]
    fn test_complement_ok() {
        assert_eq!(
            to_complement("ATCGATCGATCGATCG"),
            String::from("TAGCTAGCTAGCTAGC")
        );
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(
            to_reverse_complement("GTGCCAGCMGCCGCGGTAA"),
            "TTACCGCGGCKGCTGGCAC"
        );
    }

    #[test]
    fn test_from_iupac_to_regex() {
        assert_eq!(from_iupac_to_regex("ATCGYATCH"), "ATCG[CT]ATC[ACT]");
    }

    #[test]
    fn test_from_iupac_to_regex_2() {
        assert_eq!(from_iupac_to_regex("GGG.AAA"), "GGG.AAA");
    }
}
