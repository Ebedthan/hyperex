[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![CI](https://github.com/Ebedthan/hyvrex/workflows/CI/badge.svg)
[![CodeCov](https://codecov.io/gh/Ebedthan/hyvrex/branch/main/graph/badge.svg)](https://codecov.io/gh/Ebedthan/hyvrex)

# hyvrex: HyperVariable Primer-based Region Extractor

## Introduction

Hyvrex is a tool that extract hypervariable region from 16S rRNA. 
Provided with a fasta file, hyvrex search all hypervariable regions 
in the gene. The gene name can be changed using the --gene option.


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

## Authors
* [Anicet Ebou](https://orcid.org/0000-0003-4005-177X)