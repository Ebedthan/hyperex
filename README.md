[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![CI](https://github.com/Ebedthan/hyvrex/workflows/CI/badge.svg)
[![CodeCov](https://codecov.io/gh/Ebedthan/hyvrex/branch/main/graph/badge.svg)](https://codecov.io/gh/Ebedthan/hyvrex)


```text
                     _
                    | |                                              
                    | |__  _   ___   ___ __ _____  __
                    | '_ \| | | \ \ / / '__/ _ \ \/ /
                    | | | | |_| |\ V /| | |  __/>  < 
                    |_| |_|\__, | \_/ |_|  \___/_/\_\  v0.2.0
                            __/ |                    
                           |___/ 
                
                   HyperVariable Primer-based Region Extractor
```

## Introduction

Hyvrex (HyperVariable Region EXtractor) is a tool that extract hypervariable region from 16S/18S/23S rRNA based on supplied primer sequences. 
Hyvrex have built-in primer sequences for 16S rRNA region that can be specified using their names (with `--forward-primer`, `--reverse-primer`) or region names (with `--region`). The user can also use custom primer with `--forward-primer` and `--reverse-primer` options.

## Command-line arguments

```
USAGE:
    hyvrex [FLAGS/OPTIONS] <FILE>

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
    -o, --out <FILE>                 Specifies the ouput file [default: hyvrex_out.fa]

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
cargo install hyvrex
```

## From source
```
git clone https://github.com/Ebedthan/hyvrex.git
cd hyvrex

cargo build --release
cargo test

# To install hyvrex in current directory
cargo install --path .
```

## Note
Hyvrex use colored output in help, nevertheless hyvrex honors [NO_COLORS](https://no-color.org/) environment variable.

## Bugs
Submit problems or requests to the [Issue Tracker](https://github.com/Ebedthan/hyvrex/issues).

## License
Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.