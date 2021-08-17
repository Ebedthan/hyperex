// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{crate_version, App, AppSettings, Arg};

pub fn build_app() -> App<'static, 'static> {
    let clap_color_setting = if std::env::var_os("NO_COLOR").is_none() {
        AppSettings::ColoredHelp
    } else {
        AppSettings::ColorNever
    };

    let app = App::new("hyperex")
        .version(crate_version!())
        .usage(
            "hyperex [FLAGS/OPTIONS] <FILE>\n\n\
            EXAMPLES:\n\
            \t# With built-in 16S primer names\n\
            \thyperex -f 27F -r 337R file.fa.gz\n\n\
            \t# With built-in 16S region names\n\
            \thyperex --region v3v4 file.fa.xz\n\n\
            \t# With custom primer sequences\n\
            \thyperex -o outfile --forward-primer ATCG --reverse-primer TYAATG file.fa.bz2"
        )
        .setting(clap_color_setting)
        .setting(AppSettings::DeriveDisplayOrder)
        .after_help(
            "Note: `hyperex -h` prints a short and concise overview while `hyperex --help` gives all \
                 details.",
        )
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("Hypervariable region primer-based extractor")
        .arg(
            Arg::with_name("FILE")
                .help("Input fasta file. Can be gzip'd, xz'd or bzip'd")
                .long_help("Input fasta file. Can be gzip'd, xz'd or bzip'd")
                .index(1),
        )
        .arg(
            Arg::with_name("forward_primer")
                .short("f")
                .long("forward-primer")
                .help("Specifies forward primer sequence. Can be a\nsequence with degenerate bases")
                .conflicts_with("region")
                .requires("reverse_primer")
                .multiple(true)
                .number_of_values(1)
                .value_name("PRIMER")
        )
        .arg(
            Arg::with_name("reverse_primer")
                .short("r")
                .long("reverse-primer")
                .help("Specifies reverse primer sequence. Can be a\nsequence with degenerate bases")
                .conflicts_with("region")
                .multiple(true)
                .number_of_values(1)
                .value_name("PRIMER")
        )
        .arg(
            Arg::with_name("region")
                .long("region")
                .help("Specifies a hypervariable region to extract")
                .possible_values(&["v1v2", "v1v3", "v1v9", "v3v4", "v3v5", "v4", "v4v5", "v5v7", "v6v9", "v7v9"])
                .hide_possible_values(true)
                .multiple(true)
                .number_of_values(1)
                .value_name("REGION")
        )
        .arg(
            Arg::with_name("mismatch")
                .help("Specifies number of allowed mismatch")
                .long("mismatch")
                .short("m")
                .value_name("N")
                .hide_possible_values(true)
                .default_value("0")
                .takes_value(true)
        )
        .arg(
            Arg::with_name("prefix")
                .help("Specifies the prefix for the output files")
                .short("p")
                .long("prefix")
                .value_name("PATH")
                .default_value("hyperex_out"),
        )
        .arg(
            Arg::with_name("force")
                .help("Force output overwritting")
                .long("force")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("quiet")
                .long_help("Decreases program verbosity")
                .short("q")
                .long("quiet")
                .takes_value(false),
        );

    app
}
