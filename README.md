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

Hyperex (pronounced "hype rex" for hypervariable region extractor) is a tool that extracts 16S ribosomal RNA (rRNA) hypervariable region based on a set of primers. By default when no option is specified, hyperex extracts all hypervariable region from the supplied sequences assuming 16S rRNA sequences. To do this it has a set of built-in primer sequences which are universal 16S primers sequences.
Nevertheless, the user can choose to specify the wanted region by specifying the `--region` option or by providing the primer sequences using `--forward-primer` and `--reverse-primer`. The `--region` option takes only the region names like "v1v2" or "v4v5" while the `--forward-primer` and `--reverse-primer` takes only the sequences which can contains IUPAC ambiguities.  
For more than one needed region, one can use multiple time the `--region`, `--forward-primer`, `reverse-primer` options to specify the wanted region. Theses option takes only one argument, but can be repeat multiple time (see Examples below).

For more praticability, the user can also provide a supplied file containing primer sequences to extract the wanted region using the `--region` option. The primer sequences file should be a no header tab separated value file like:
```
FORWARD_PRIMER_1\tREVERSE_PRIMER_1
FORWARD_PRIMER_2\tREVERSE_PRIMER_2
...
```

Moreover, one can allow a number of mismatch in the primer sequence using the `--mismatch` option.

The outputs are a fasta file containing the extracted regions and a GFF3 file indicating the extracted regions positions.

## How to run hyperex ?

### By default with no options

```
hyperex file.fa

cat file.fa | hyperex
```

### Using built-in 16S primer names

```
hyperex -f 27F -r 337R file.fa.gz

zcat file.fa.gz | hyperex -f 27F -r 337R
```

### Using built-in 16S region names

```
hyperex --region v3v4 file.fa.xz

xzcat file.fa.xz | hyperex --region v3v4
```

### Using custom primer sequences

```
hyperex -p prefix --forward-primer ATCG --reverse-primer TYAATG file.fa.bz2

bzcat file.fa.bz2 | hyperex -p prefix --forward-primer ATCG --reverse-primer TYAATG
```

### Using custom list of primers: primers.txt

```
hyperex --region primers.txt file.fa
```

### Using multiple primers

```
hyperex --region v1v2 --region v3v4 file.fa

hyperex -f ATC -f YGA -r GGCC -r TTRC file.fa
```

## Hyperex command-line arguments

```
FLAGS:
        --force      Force output overwritting
    -q, --quiet      Decreases program verbosity
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --forward-primer <PRIMER>...    Specifies forward primer sequence. Can be a
                                        sequence, a regexp or a prosite pattern
                                        with degenerate bases
    -r, --reverse-primer <PRIMER>...    Specifies reverse primer sequence. Can be a
                                        sequence, a regexp or a prosite pattern
                                        with degenerate bases
        --region <REGION>...            Specifies a hypervariable region to extract
    -m, --mismatch <N>                  Specifies number of allowed mismatch [default: 0]
    -p, --prefix <PATH>                 Specifies the prefix for the output files [default: hyperex_out]

ARGS:
    <FILE>    Input fasta file. Can be gzip'd, xz'd or bzip'd

Note: `hyperex -h` prints a short and concise overview while `hyperex --help` gives all details.
```

## Requirements

### Mandatory
- [Rust](https://rust-lang.org) in stable channel

### Optional
- libgz for gz file support
- liblzma for xz file support
- libbzip2 for bzip2 file support


## Installation

### Using prebuilt binaries

Please see the [release page](https://github.com/Ebedthan/hyperex/releases) for prebuilt binaries for major operating system

### From crates.io
If you already have a functional rust installation you can easily do:

```
cargo install hyperex
```

And that's all!

### From source
```
git clone https://github.com/Ebedthan/hyperex.git
cd hyperex

cargo build --release
cargo test

# To install hyperex in current directory
cargo install --path .
```

And you are good to go!

## Note
hyperex use colored output in help, nevertheless hyperex honors [NO_COLORS](https://no-color.org/) environment variable.

## Bugs
Submit problems or requests to the [Issue Tracker](https://github.com/Ebedthan/hyperex/issues).
