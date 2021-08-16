```
                                 _    _                                
                                | |  | |                               
                                | |__| |_   _ _ __   ___ _ __ _____  __
                                |  __  | | | | '_ \ / _ \ '__/ _ \ \/ /
                                | |  | | |_| | |_) |  __/ | |  __/>  < 
                                |_|  |_|\__, | .__/ \___|_|  \___/_/\_\
                                         __/ | |                       
                                        |___/|_|   
```

<p align="center">
    <br/>
    <a href="https://github.com/Ebedthan/hyperex/blob/master/LICENSE">
        <img src="https://img.shields.io/badge/license-MIT-blue?style=flat">
    </a>
    <a href="https://github.com/Ebedthan/hyperex/workflows/CI">
        <img src="https://github.com/Ebedthan/hyperex/workflows/CI/badge.svg">
    </a>
    <a href="https://codecov.io/gh/Ebedthan/hyperex">
        <img src="https://codecov.io/gh/Ebedthan/hyperex/branch/main/graph/badge.svg">
    </a>
</p>   


## About

Hyperex (hypervariable region extractor) is a tool that extract hypervariable region from 16S/18S/23S rRNA based on supplied primer sequences. 
Hyperex have built-in primer sequences for 16S rRNA region that can be specified using their names (with `--forward-primer`, `--reverse-primer`) or region names (with `--region`). The user can also use custom primer with `--forward-primer` and `--reverse-primer` options.

## Some examples

```
# With built-in 16S primer names
hyperex -f 27F -r 337R file.fa.gz

# With built-in 16S region names
hyperex --region v3v4 file.fa.xz

# With custom primer sequences
hyperex -o outfile --forward-primer ATCG --reverse-primer TYAATG file.fa.bz2
```

## Command-line arguments

```
USAGE:
    hyperex [FLAGS/OPTIONS] <FILE>

FLAGS:
    -v, --verbose    Increases program verbosity each use for up to 3 times
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --forward-primer <PRIMER>    Specifies forward primer sequence. Can be a
                                     sequence, a regexp or a prosite pattern
                                     with degenerate bases
    -r, --reverse-primer <PRIMER>    Specifies reverse primer sequence. Can be a
                                     sequence, a regexp or a prosite pattern
                                     with degenerate bases
        --region <REGION>            Specifies a hypervariable region to extract
    -o, --out <FILE>                 Specifies the ouput file [default: hyperex_out.fa]

ARGS:
    <FILE>    Input fasta file. Can be gzip'd, xz'd or bzip'd
```

## Requirements

- [Rust](https://rust-lang.org) in stable channel
- libgz for gz file support
- liblzma for xz file support
- libbzip2 for bzip2 file support


## Installation

## From crates.io
If you already have a functional rust installation do:

```
cargo install hyperex
```

## From source
```
git clone https://github.com/Ebedthan/hyperex.git
cd hyperex

cargo build --release
cargo test

# To install hyperex in current directory
cargo install --path .
```

## Note
hyperex use colored output in help, nevertheless hyperex honors [NO_COLORS](https://no-color.org/) environment variable.

## Bugs
Submit problems or requests to the [Issue Tracker](https://github.com/Ebedthan/hyperex/issues).

## License
Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
