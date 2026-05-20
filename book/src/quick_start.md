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
cargo install --path modkit
# or
cargo install --git https://github.com/nanoporetech/modkit.git
```

## Common Use Cases
1. [Creating a bedMethyl table with `pileup`](./intro_pileup.md)
1. [Summarizing a modBAM with `modbam summary`](./intro_summary.md)
1. [Extracting per-read base modification data into a table](./intro_extract.md)
1. [Checking modified base tags in a modBAM](./intro_modbam_check_tags.md)
1. [Making a motif BED file with `motif-bed`](./intro_motif_bed.md)
1. [Performing differential methylation scoring with `dmr`](./intro_dmr.md)
1. [Investigating differential methylation in direct RNA](./intro_rna.md)
1. [Convert bedMethyl files to bigWig for visualization](./intro_bedmethyl_merge.md#convert-bedmethyl-to-bigwig)
1. [Predict regions of open chromatin on MTase-treated DNA](./intro_open_chromatin.md)
1. [Updating and Adjusting MM tags with `adjust-mods` and `update-tags`](./intro_adjust.md)
1. [Convert modification probabilities into hard calls](./intro_call_mods.md)
1. [Removing base modification calls at the ends of reads](./intro_edge_filter.md)
1. [Narrow analysis to only specific positions with a BED file](./intro_include_bed.md)
1. [Repairing/adding MM/ML tags to reads with clipped sequences](./intro_repair.md)
1. [Creating hemi-methylation pattern bedMethyl tables with `pileup-hemi`](./intro_pileup_hemi.md)

## Notes and troubleshooting
1. [General troubleshooting](./troubleshooting.md)
1. [Threshold evaluation examples](./filtering_details.md) (for advanced users)
1. [Querying the logs in `motif search`](./motif_search_structured_logging.md)

