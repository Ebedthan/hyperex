// Copyright 2021-2022 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{crate_version, Command, ColorChoice, AppSettings, Arg};

pub fn build_app() -> Command<'static> {
    let clap_color_setting = if std::env::var_os("NO_COLOR").is_none() {
        ColorChoice::Always
    } else {
        ColorChoice::Never
    };

    let app = Command::new("hyperex")
        .version(crate_version!())
        .override_usage(
            "hyperex [options] [<FILE>]"
        )
        .color(clap_color_setting)
        .setting(AppSettings::DeriveDisplayOrder)
        .after_help(
            "Note: `hyperex -h` prints a short and concise overview while `hyperex --help` gives all \
                 details.",
        )
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("Hypervariable region primer-based extractor")
        .arg(
            Arg::new("FILE")
                .help("input fasta file or stdin")
                .long_help("input fasta file. With no FILE, or when FILE is -, read standard input. Input data can be gzip'd, xz'd or bzip'd")
                .index(1),
        )
        .arg(
            Arg::new("forward_primer")
                .short('f')
                .long("forward-primer")
                .help("forward primer sequence")
                .long_help("Specifies forward primer sequence which can contains IUPAC ambiguities")
                .conflicts_with("region")
                .requires("reverse_primer")
                .multiple_occurrences(true)
                .number_of_values(1)
                .value_name("STR")
        )
        .arg(
            Arg::new("reverse_primer")
                .short('r')
                .long("reverse-primer")
                .help("reverse primer sequence")
                .long_help("Specifies reverse primer sequence which can contains IUPAC ambiguities")
                .conflicts_with("region")
                .multiple_occurrences(true)
                .number_of_values(1)
                .value_name("STR")
        )
        .arg(
            Arg::new("region")
                .long("region")
                .help("hypervariable region name")
                .long_help(
                    "Specifies 16S rRNA region name wanted. Supported values are\n\
                    v1v2, v1v3, v1v9, v3v4, v3v5, v4, v4v5, v5v7, v6v9, v7v9"
                )
                .possible_values(&["v1v2", "v1v3", "v1v9", "v3v4", "v3v5", "v4", "v4v5", "v5v7", "v6v9", "v7v9"])
                .hide_possible_values(true)
                .multiple_occurrences(true)
                .number_of_values(1)
                .value_name("STR")
        )
        .arg(
            Arg::new("mismatch")
                .help("number of allowed mismatch")
                .long_help(
                    "Specifies the number of allowed mismatch. This cannot\n\
                    be greater than the length of the lengthest primer"
                )
                .long("mismatch")
                .short('m')
                .value_name("N")
                .hide_possible_values(true)
                .default_value("0")
                .takes_value(true)
        )
        .arg(
            Arg::new("prefix")
                .help("prefix of output files")
                .long_help("Specifies the prefix for output files. Paths are supported")
                .short('p')
                .long("prefix")
                .value_name("PATH")
                .default_value("hyperex_out"),
        )
        .arg(
            Arg::new("force")
                .help("overwrite output")
                .long("force")
                .takes_value(false),
        )
        .arg(
            Arg::new("quiet")
                .long_help("decreases program verbosity")
                .short('q')
                .long("quiet")
                .takes_value(false),
        );

    app
}
