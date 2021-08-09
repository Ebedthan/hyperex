// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{crate_version, App, AppSettings, Arg};

use std::path::Path;

pub fn build_app() -> App<'static, 'static> {
    let clap_color_setting = if std::env::var_os("NO_COLOR").is_none() {
        AppSettings::ColoredHelp
    } else {
        AppSettings::ColorNever
    };

    let app = App::new("hyvrex")
        .version(crate_version!())
        .usage(
            "hyvrex [FLAGS/OPTIONS] <FILE>\n\n\
            EXAMPLES:\n\
            \t# With built-in 16S primer names\n\
            \thyvrex -f 27F -r 337R file.fa.gz\n\n\
            \t# With built-in 16S region names\n\
            \thyvrex --region v3v4 file.fa.xz\n\n\
            \t# With custom primer sequences\n\
            \thyvrex -o outfile --forward-primer ATCG --reverse-primer TYAATG file.fa.bz2"
        )
        .setting(clap_color_setting)
        .setting(AppSettings::DeriveDisplayOrder)
        .after_help(
            "Note: `hyvrex -h` prints a short and concise overview while `hyvrex --help` gives all \
                 details.",
        )
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("Hypervariable region primer-based extractor")
        .arg(
            Arg::with_name("FILE")
                .help("Input fasta file. Can be gzip'd, xz'd or bzip'd")
                .long_help("Input fasta file. Can be gzip'd, xz'd or bzip'd")
                .required(true)
                .index(1)
                .validator(is_fasta),
        )
        .arg(
            Arg::with_name("forward_primer")
                .short("f")
                .long("forward-primer")
                .help("Specifies forward primer sequence. Can be a\nsequence, a regexp or a prosite pattern\nwith degenerate bases")
                .conflicts_with("region")
                .requires("reverse_primer")
                .value_name("PRIMER")
        )
        .arg(
            Arg::with_name("reverse_primer")
                .short("r")
                .long("reverse-primer")
                .help("Specifies reverse primer sequence. Can be a\nsequence, a regexp or a prosite pattern\nwith degenerate bases")
                .conflicts_with("region")
                .value_name("PRIMER")
        )
        .arg(
            Arg::with_name("region")
                .long("region")
                .help("Specifies a hypervariable region to extract")
                .possible_values(&["v1v2", "v1v3", "v1v9", "v3v4", "v3v5", "v4", "v4v5", "v5v7", "v6v9", "v7v9"])
                .hide_possible_values(true)
                .value_name("REGION")
        )
        .arg(
            Arg::with_name("mismatch")
                .help("Specifies number of allowed mismatch")
                .long("mismatch")
                .short("m")
                .value_name("N")
                .possible_values(&["0", "1", "2", "3", "4"])
                .hide_possible_values(true)
                .default_value("0")
                .takes_value(true)
        )
        .arg(
            Arg::with_name("output")
                .help("Specifies the ouput file")
                .short("o")
                .long("out")
                .value_name("FILE")
                .default_value("hyvrex_out.fa")
                .validator(already_exists),
        )
        .arg(
            Arg::with_name("verbose")
                .long_help("Increases program verbosity each use for up to 3 times")
                .short("v")
                .long("verbose")
                .multiple(true)
                .takes_value(false),
        );

    app
}

fn is_fasta(filename: String) -> Result<(), String> {
    if Path::new(&filename).exists()
        && (filename.contains(".fasta")
            || filename.contains(".fa")
            || filename.contains(".fas")
            || filename.contains(".fna"))
    {
        Ok(())
    } else {
        Err(String::from("Is file path correct? with file extension (fa|fas|fasta|fna) clearly stated with appropriate permission?"))
    }
}

fn already_exists(filename: String) -> Result<(), String> {
    if !Path::new(&filename).exists() {
        Ok(())
    } else {
        Err(String::from(
            "Selected file already exists. Please change it using --out option",
        ))
    }
}
