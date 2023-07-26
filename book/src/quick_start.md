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
2. [Updating and Adjusting MM tags with `adjust-mods` and `update-tags`](./intro_adjust.md)
3. [Summarizing a modBAM with `summary`](./intro_summary.md)
4. [Making a motif BED file with `motif-bed`](./intro_motif_bed.md)
5. [Extracting per-read base modification data into a table](./intro_extract.md)
6. [Convert modification probabilities into hard calls](./intro_call_mods.md)
7. [Removing base modification calls at the ends of reads](./intro_edge_filter.md)
8. [Narrow analysis to only specific positions with a BED file](./intro_include_bed.md)
9. [Repairing/adding MM/ML tags to reads with clipped sequences](./intro_repair.md)

## Notes and troubleshooting
1. [General troubleshooting](./troubleshooting.md)
2. [Threshold evaluation examples](./filtering_details.md) (for advanced users)

