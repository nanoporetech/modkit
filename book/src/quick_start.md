# Basic Usage

### `modkit` is a bioinformatics tool for working with modified bases from Oxford Nanopore.

![ONT_logo](./images/ONT_logo_590x106.png)

## Installation

Pre-compiled binaries are provided for Linux from the [release
page](https://github.com/nanoporetech/modkit/releases). We recommend the use of these in
most circumstances. As a rust-based project, `modkit` can also be installed with 
[cargo](https://www.rust-lang.org/learn/get-started).
```bash
git clone https://github.com/nanoporetech/modkit.git
cd modkit
cargo install --path .
# or
cargo install --git https://github.com/nanoporetech/modkit.git
```

## Common Use Cases
1. [Creating a bedMethyl table with `pileup`](./intro_bedmethyl.md)
1. [Updating and Adjusting MM tags with `adjust-mods` and `update-tags`](./intro_adjust.md)
1. [Summarizing a modBAM with `summary`](./intro_summary.md)
1. [Making a motif BED file with `motif-bed`](./intro_motif_bed.md)

## Notes and troubleshooting
1. [General troubleshooting](./troubleshooting.md)

