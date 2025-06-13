// Copyright 2021-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{Parser, ValueEnum};

#[derive(Parser, Debug)]
#[command(
    name = "hyperex",
    version,
    author = "Anicet Ebou <anicet.ebou@gmail.com>",
    about = "Hypervariable region primer-based extractor",
    override_usage = "hyperex [options] [<FILE>]"
)]
pub struct Args {
    /// Input fasta file or stdin
    #[arg(
        long_help = "input fasta file. With no FILE, or when FILE is -, read standard input. Input data can be gzip'd, xz'd or bzip'd"
    )]
    pub file: Option<String>,

    /// Forward primer sequence
    #[arg(
        short = 'f',
        long = "forward",
        long_help = "Specifies forward primer sequence which can contains IUPAC ambiguities",
        conflicts_with = "region",
        requires = "reverse",
        value_name = "STR",
        num_args = 1..
    )]
    pub forward: Option<Vec<String>>,

    /// Reverse primer sequence
    #[arg(
        short = 'r',
        long = "reverse",
        long_help = "Specifies reverse primer sequence which can contains IUPAC ambiguities",
        conflicts_with = "region",
        value_name = "STR",
        num_args = 1..
    )]
    pub reverse: Option<Vec<String>>,

    /// Hypervariable region name
    #[arg(
        long = "region",
        long_help = "Specifies 16S rRNA region name wanted. Supported values are\nv1v1, v1v3, v1v9, v3v4, v3v5, v4, v4v5, v5v7, v6v9, v7v9",
        value_name = "STR",
        hide_possible_values = true,
        num_args = 1..
    )]
    pub region: Option<Vec<Region>>,

    /// Number of allowed mismatch
    #[arg(
        short = 'm',
        long = "mismatch",
        long_help = "Specifies the number of allowed mismatch. This cannot\nbe greate than the length of the lengthest primer",
        value_name = "N",
        hide_possible_values = true,
        default_value = "0"
    )]
    pub mismatch: u8,

    /// Prefix of output files
    #[arg(
        short = 'p',
        long = "prefix",
        long_help = "Specifies the prefix for output files. Paths are supported",
        value_name = "PATH",
        default_value = "hyperex_out"
    )]
    pub prefix: String,

    /// Overwrite output
    #[arg(long = "force")]
    pub force: bool,

    /// Decreases program verbosity
    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,
}

#[derive(ValueEnum, Clone, Debug)]
pub enum Region {
    V1V2,
    V1V3,
    V1V9,
    V3V4,
    V3V5,
    V4,
    V4V5,
    V5V7,
    V6V9,
    V7V9,
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Region::V1V2 => write!(f, "v1v2"),
            Region::V1V3 => write!(f, "v1v3"),
            Region::V1V9 => write!(f, "v1v9"),
            Region::V3V4 => write!(f, "v3v4"),
            Region::V3V5 => write!(f, "v3v5"),
            Region::V4 => write!(f, "v4"),
            Region::V4V5 => write!(f, "v4v5"),
            Region::V5V7 => write!(f, "v5v7"),
            Region::V6V9 => write!(f, "v6v9"),
            Region::V7V9 => write!(f, "v7v9"),
        }
    }
}
