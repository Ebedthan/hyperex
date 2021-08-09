// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate chrono;
extern crate fern;
extern crate log;
extern crate niffler;
extern crate phf;

use anyhow::{anyhow, Context, Result};
use bio::io::fasta;
use bio::pattern_matching::myers::MyersBuilder;
use fern::colors::ColoredLevelConfig;
use log::{error, info, warn};
use phf::phf_map;

use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::path::Path;

pub fn setup_logging() -> Result<(), fern::InitError> {
    let colors = ColoredLevelConfig::default();

    fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "{}[{}] {}",
                chrono::Local::now().format("[%H:%M:%S]"),
                colors.color(record.level()),
                message
            ))
        })
        .level(log::LevelFilter::Debug)
        .chain(std::io::stdout())
        .chain(fern::log_file("hyvrex.log")?)
        .apply()?;

    Ok(())
}

pub fn is_fasta(filename: &str) -> Result<(), &str> {
    if Path::new(&filename).exists()
        && (filename.contains(".fasta")
            || filename.contains(".fa")
            || filename.contains(".fas")
            || filename.contains(".fna"))
    {
        Ok(())
    } else {
        Err("Is file path correct? with file extension (fa|fas|fasta|fna) clearly stated with appropriate permission?")
    }
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

pub fn region_to_primer(region: &str) -> Result<Vec<&str>> {
    match region {
        "v1v2" => Ok(vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["337R"]]),
        "v1v3" => Ok(vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["534R"]]),
        "v1v9" => Ok(vec![FORWARD_PRIMERS["27F"], REVERSE_PRIMERS["1492Rmod"]]),
        "v3v4" => Ok(vec![FORWARD_PRIMERS["341F"], REVERSE_PRIMERS["805R"]]),
        "v3v5" => Ok(vec![FORWARD_PRIMERS["341F"], REVERSE_PRIMERS["926Rb"]]),
        "v4" => Ok(vec![FORWARD_PRIMERS["515F"], REVERSE_PRIMERS["806R"]]),
        "v4v5" => {
            Ok(vec![FORWARD_PRIMERS["515F-Y"], REVERSE_PRIMERS["909-928R"]])
        }
        "v5v7" => Ok(vec![FORWARD_PRIMERS["799F"], REVERSE_PRIMERS["1193R"]]),
        "v6v9" => {
            Ok(vec![FORWARD_PRIMERS["928F"], REVERSE_PRIMERS["1492Rmod"]])
        }
        "v7v9" => {
            Ok(vec![FORWARD_PRIMERS["1100F"], REVERSE_PRIMERS["1492Rmod"]])
        }
        _ => Ok(vec![""]),
    }
}

pub fn combine_vec<'a>(
    first: Vec<&'a str>,
    second: Vec<&'a str>,
) -> Vec<Vec<&'a str>> {
    let mut newvec: Vec<Vec<&str>> = Vec::new();

    for i in 0..first.len() {
        newvec.push(vec![first[i], second[i]]);
    }

    newvec
}

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

fn to_complement(primer: &str, alphabet: &str) -> String {
    let mut complement = String::new();

    if alphabet == "dna" {
        complement = primer
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
    } else if alphabet == "rna" {
        complement = primer
            .chars()
            .map(|x| match x {
                'A' => 'U',
                'U' => 'A',
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
    }

    complement
}

fn to_reverse_complement(primer: &str, alphabet: &str) -> String {
    let complement = to_complement(primer, alphabet);
    let reverse_complement = complement.chars().rev().collect();

    reverse_complement
}

#[derive(Debug, PartialEq)]
pub enum Alphabet {
    Dna,
    Rna,
}

pub fn sequence_type(sequence: &str) -> Option<Alphabet> {
    let valid_dna_iupac = "ACGTRYSWKMBDHVN";
    let valid_rna_iupac = "ACGURYSWKMBDHVN";

    if sequence.chars().all(|x| valid_dna_iupac.contains(x)) {
        Some(Alphabet::Dna)
    } else if sequence.chars().all(|x| valid_rna_iupac.contains(x)) {
        Some(Alphabet::Rna)
    } else {
        None
    }
}

pub fn process_fa(
    file: &str,
    primers: Vec<Vec<&str>>,
    outfile: &str,
    mismatch: u8,
) -> Result<()> {
    let (reader, mut _compression) =
        read_file(file).with_context(|| "Cannot read file")?;

    let mut records = fasta::Reader::new(reader).records();

    let out = OpenOptions::new().create(true).append(true).open(outfile)?;

    // Build Myers with IUPAC ambiguities in patterns
    let ambigs = [
        (b'M', &b"AC"[..]),
        (b'R', &b"AG"[..]),
        (b'W', &b"AT"[..]),
        (b'S', &b"CG"[..]),
        (b'Y', &b"CT"[..]),
        (b'K', &b"GT"[..]),
        (b'V', &b"ACGMRS"[..]),
        (b'H', &b"ACTMWY"[..]),
        (b'D', &b"AGTRWK"[..]),
        (b'B', &b"CGTSYK"[..]),
        (b'N', &b"ACGTMRWSYKVHDB"[..]),
    ];

    let mut builder = MyersBuilder::new();

    for &(base, equivalents) in &ambigs {
        builder.ambig(base, equivalents);
    }

    while let Some(Ok(record)) = records.next() {
        let seq = record.seq();
        let mut alphabet = "";
        match sequence_type(std::str::from_utf8(seq)?) {
            Some(alp) => {
                if alp == Alphabet::Dna {
                    info!("Sequence type is DNA");
                    alphabet = "dna";
                } else if alp == Alphabet::Rna {
                    info!("Sequence type is RNA");
                    alphabet = "rna";
                }
            }
            None => error!("Sequence type is not recognized as DNA or RNA"),
        }
        if seq.len() <= 1400 {
            warn!("Sequence length is less than 1400 bp. We may not be able to find some regions");
        }

        for primer_pair in primers.iter() {
            let region = primers_to_region(primer_pair.to_vec());

            let mut forward_myers = builder.build_64(primer_pair[0].as_bytes());
            let mut reverse_myers = builder.build_64(
                to_reverse_complement(primer_pair[1], alphabet).as_bytes(),
            );

            let mut forward_matches =
                forward_myers.find_all_lazy(seq, mismatch);
            let mut reverse_matches =
                reverse_myers.find_all_lazy(seq, mismatch);

            // Get the best hit
            let forward_best_hit =
                forward_matches.by_ref().min_by_key(|&(_, dist)| dist);
            let reverse_best_hit =
                reverse_matches.by_ref().min_by_key(|&(_, dist)| dist);

            match forward_best_hit {
                Some((forward_best_hit_end, _)) => {
                    match reverse_best_hit {
                        Some((reverse_best_hit_end, _)) => {
                            // Get match start position of forward primer
                            let (forward_start, _) = forward_matches
                                .hit_at(forward_best_hit_end)
                                .unwrap();
                            // Get match start position of reverse primer
                            let (reverse_start, _) = reverse_matches
                                .hit_at(reverse_best_hit_end)
                                .unwrap();

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
                                        &seq[forward_start
                                            ..reverse_start + primer_pair[1].len()],
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
                                        &seq[forward_start
                                            ..reverse_start
                                                + primer_pair[1].len()],
                                    ),
                                )?;
                            }
                        }
                        None => {
                            warn!("Region {} not found because primer {} was not found in the sequence", region, primer_pair[1])
                        }
                    }
                }
                None => match reverse_best_hit {
                    Some((_, _)) => {
                        warn!("Region {} not found because primer {} was not found in the sequence", region, primer_pair[0]);
                    }
                    None => {
                        warn!("Region {} not found because primers {}, {} was not found in the sequence", region, primer_pair[0], primer_pair[1])
                    }
                },
            }
        }
    }

    Ok(())
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{tempfile, NamedTempFile};

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
    fn test_complement_dna() {
        assert_eq!(
            to_complement("ATCGATCGATCGATCGRYKBVDH", "dna"),
            String::from("TAGCTAGCTAGCTAGCYRMVBHD")
        );
    }

    #[test]
    fn test_complement_rna() {
        assert_eq!(
            to_complement("AUCGAUCGAUCGAUCGRYKBVDHM", "rna"),
            String::from("UAGCUAGCUAGCUAGCYRMVBHDK")
        );
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(
            to_reverse_complement("GTGCCAGCMGCCGCGGTAA", "dna"),
            "TTACCGCGGCKGCTGGCAC"
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
            vec![vec!["AGAGTTTGATCMTGGCTCAG", "CYIACTGCTGCCTCCCGTAG"]]
        );
        assert_eq!(
            region_to_primer("v1v3").unwrap(),
            vec![vec!["AGAGTTTGATCMTGGCTCAG", "ATTACCGCGGCTGCTGG"]]
        );
        assert_eq!(
            region_to_primer("v1v9").unwrap(),
            vec![vec!["AGAGTTTGATCMTGGCTCAG", "TACGGYTACCTTGTTAYGACTT"]]
        );
        assert_eq!(
            region_to_primer("v3v4").unwrap(),
            vec![vec!["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"]]
        );
        assert_eq!(
            region_to_primer("v3v5").unwrap(),
            vec![vec!["CCTACGGGNGGCWGCAG", "CCGTCAATTYMTTTRAGT"]]
        );
        assert_eq!(
            region_to_primer("v4").unwrap(),
            vec![vec!["GTGCCAGCMGCCGCGGTAA", "GGACTACHVGGGTWTCTAAT"]]
        );
        assert_eq!(
            region_to_primer("v4v5").unwrap(),
            vec![vec!["GTGYCAGCMGCCGCGGTAA", "CCCCGYCAATTCMTTTRAGT"]]
        );
        assert_eq!(
            region_to_primer("v5v7").unwrap(),
            vec![vec!["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"]]
        );
        assert_eq!(
            region_to_primer("v6v9").unwrap(),
            vec![vec!["TAAAACTYAAAKGAATTGACGGGG", "TACGGYTACCTTGTTAYGACTT"]]
        );
        assert_eq!(
            region_to_primer("v7v9").unwrap(),
            vec![vec!["YAACGAGCGCAACCC", "TACGGYTACCTTGTTAYGACTT"]]
        );
        assert_eq!(
            region_to_primer("").unwrap(),
            vec![
                vec!["AGAGTTTGATCMTGGCTCAG", "CYIACTGCTGCCTCCCGTAG"],
                vec!["AGAGTTTGATCMTGGCTCAG", "ATTACCGCGGCTGCTGG"],
                vec!["AGAGTTTGATCMTGGCTCAG", "TACGGYTACCTTGTTAYGACTT"],
                vec!["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"],
                vec!["CCTACGGGNGGCWGCAG", "CCGTCAATTYMTTTRAGT"],
                vec!["GTGCCAGCMGCCGCGGTAA", "GGACTACHVGGGTWTCTAAT"],
                vec!["GTGYCAGCMGCCGCGGTAA", "CCCCGYCAATTCMTTTRAGT"],
                vec!["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"],
                vec!["TAAAACTYAAAKGAATTGACGGGG", "TACGGYTACCTTGTTAYGACTT"],
                vec!["YAACGAGCGCAACCC", "TACGGYTACCTTGTTAYGACTT"]
            ]
        );
    }

    #[test]
    fn test_write_fa_ok() {
        let record =
            fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
        let file = tempfile().expect("Cannot create temp file");

        assert!((write_fa(&file, record)).is_ok());
    }

    #[test]
    fn test_write_fa_ok2() {
        let mut tmpfile =
            NamedTempFile::new().expect("Cannot create temp file");
        writeln!(tmpfile, ">id_str desc\nATCGCCG")
            .expect("Cannot write to tmp file");

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
        let first = vec!["ab", "cd", "ef"];
        let second = vec!["cd", "ef", "gh"];
        assert_eq!(
            combine_vec(first, second),
            [["ab", "cd"], ["cd", "ef"], ["ef", "gh"]]
        );
    }

    fn test_combine_vec_not_ok() {
        let first = vec!["ab", "cd", "ef"];
        let second = vec!["ab"];
    }
}
