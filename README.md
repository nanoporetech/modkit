![Oxford Nanopore Technologies logo](ONT_logo_590x106.png)

# Modkit

A bioinformatics tool for working with modified bases from Oxford Nanopore.
Specifically for quickly converting modBAM to bedMethyl files using best practices, but also manipulating modBAM files and generating summary statistics.
Detailed documentation and quick-start can be found in the [online documentation](https://nanoporetech.github.io/modkit/).

## Installation

Pre-compiled binaries are provided for Linux from the [release page](https://github.com/nanoporetech/modkit/releases). We recommend the use of these in most circumstances.

### Building from source

The provided packages should be used where possible. We understand that some users may wish to compile the software from its source code. To build `modkit` from source [cargo](https://www.rust-lang.org/learn/get-started) should be used.

```bash
git clone https://github.com/nanoporetech/modkit.git
cd modkit
cargo install --path modkit
# or
cargo install --git https://github.com/nanoporetech/modkit.git
```

### macOS (Apple Silicon)

A script is provided to compile modkit on Apple Silicon Macs with Metal GPU (MPS) acceleration.
Just download the script and run it with the desired installation directory, modkit version, and Python provider (system, conda, pyenv, or uv).

```bash
bash mac_compile_modkit.sh ~/tools
```

This installs all dependencies (Homebrew, Rust, PyTorch) and compiles modkit automatically. See the [macOS installation guide](./book/src/mac_compile_modkit.md) for full details, Python version control options, and troubleshooting.

## Usage

Modkit comprises a suite of tools for manipulating modified-base data stored in [BAM](http://www.htslib.org/) files. Modified base information is stored in the `MM` and `ML` tags (see section 1.7 of the [SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf) specification). These tags are produced by contemporary basecallers of data from Oxford Nanopore Technologies sequencing platforms.

### Constructing bedMethyl tables

A primary use of `modkit` is to create summary counts of modified and unmodified bases in an extended [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) format. bedMethyl files tabulate the counts of base modifications from every sequencing read over each reference genomic position.

For best performance use the `--modified-bases` option with the base modifications you intend to analyze. 

For example:
```bash
modkit pileup \
  path/to/reads.bam \
  path/to/output.bed.gz \
  --modified-bases 5mC 5hmC \
  --reference path/to/reference.fasta \
  --log path/to/log.txt \ # optional, recommended
  --bgzf \ # optional
```

A single file (described [below](#description-of-bedmethyl-output)) with base count summaries will be created. The final argument here specifies an optional log file output.

The program performs best-practices filtering and manipulation of the raw data stored in the input file. For further details see [filtering modified-base calls](./book/src/filtering.md).

For user convenience the counting process can be modulated using several additional transforms and filters. The most basic of these is to report only counts from reference CpG dinucleotides. This option requires a reference sequence in order to locate the CpGs in the reference:

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed \
  --cpg \
  --modified-bases 5mC 5hmC \
  --ref path/to/reference.fasta
```

For more details on `pileup` and other commands, please see the [online documentation](https://nanoporetech.github.io/modkit/).

For more information on the individual options see the [Advanced Usage](./book/src/advanced_usage.md) help document.

## Advanced usage examples

For complete usage instructions please see the command-line help of the program or the [Advanced usage](./book/src/advanced_usage.md) help documentation. Some more commonly required examples are provided below.
To combine multiple base modification calls into one, for example to combine basecalls for both 5hmC and 5mC into a count for "all cytosine modifications" (with code `C`) the `--combine-mods` option can be used:

```bash
modkit pileup \
  path/to/reads.bam \
  output/path/pileup.bed \
  --modified-bases 5mC 5hmC \
  --combine-mods \
  --ref path/to/reference.fasta \
  [--cpg] \  # optional
  [--combine-strands] \ # optional
```

If you have a modBAM with phased reads containing a `HP` tag. These can be partitioned into separate bedMethyl files on output by passing the `--phased` flag.

```bash
modkit pileup \
  path/to/reads.bam \
  output/directory/ \
  --cpg \
  --modified-bases 5mC 5hmC \
  --phased \
  --ref <reference.fasta>
```

The output will be 3 files: `hp1.bedmethyl`, `hp2.bedmethyl`, and `combined.bedmethyl`.
hp1.bedmethyl and hp2.bedmethyl contain counts for records with `HP=1` and `HP=2` tags, respectively. combined.bedmethyl contains counts for all modBAM records.


## Description of bedMethyl output

Below is a description of the bedMethyl columns generated by `modkit pileup`. A brief description of the
bedMethyl specification can be found on [Encode](https://www.encodeproject.org/data-standards/wgbs/).

### Definitions:

* N<sub>mod</sub> - Number of calls passing filters that were classified as a residue with a specified base modification.
* N<sub>canonical</sub> - Number of calls passing filters were classified as the canonical base rather than modified. The
exact base must be inferred by the modification code. For example, if the modification code is `m` (5mC) then
the canonical base is cytosine. If the modification code is `a`, the canonical base is adenosine.
* N<sub>other mod</sub> - Number of calls passing filters that were classified as modified, but where the modification is different from the listed base (and the corresponding canonical base is equal). For example, for a given cytosine there may be 3 reads with
`h` calls, 1 with a canonical call, and 2 with `m` calls. In the bedMethyl row for `h` N<sub>other_mod</sub> would be 2. In the
`m` row N<sub>other_mod</sub> would be 3.
* N<sub>valid_cov</sub> - the valid coverage. N<sub>valid_cov</sub> = N<sub>mod</sub> + N<sub>other_mod</sub> + N<sub>canonical</sub>, also used as the `score` in the bedMethyl
* N<sub>diff</sub> - Number of reads with a base other than the canonical base for this modification. For example, in a row
for `h` the canonical base is cytosine, if there are 2 reads with C->A substitutions, N<sub>diff</sub> will be 2.
* N<sub>delete</sub> - Number of reads with a deletion at this reference position
* N<sub>fail</sub> - Number of calls where the probability of the call was below the threshold. The threshold can be
set on the command line or computed from the data (usually failing the lowest 10th percentile of calls).
* N<sub>nocall</sub> - Number of reads aligned to this reference position, with the correct canonical base, but without a base
modification call. This can happen, for example, if the model requires a CpG dinucleotide and the read has a
CG->CH substitution such that no modification call was produced by the basecaller.

### bedMethyl column descriptions

| column | name                  | description                                                                    | type  |
|--------|-----------------------|--------------------------------------------------------------------------------|-------|
| 1      | chrom                 | name of reference sequence from BAM header                                     | str   |
| 2      | start position        | 0-based start position                                                         | int   |
| 3      | end position          | 0-based exclusive end position                                                 | int   |
| 4      | modified base code    | single letter code for modified base                                           | str   |
| 5      | score                 | Equal to N<sub>valid_cov</sub>.                                                | int   |
| 6      | strand                | '+' for positive strand '-' for negative strand, '.' when strands are combined | str   |
| 7      | start position        | included for compatibility                                                     | int   |
| 8      | end position          | included for compatibility                                                     | int   |
| 9      | color                 | included for compatibility, always 255,0,0                                     | str   |
| 10     | N<sub>valid_cov</sub> | See definitions above.                                                         | int   |
| 11     | fraction modified     | N<sub>mod</sub> / N<sub>valid_cov</sub>                                        | float |
| 12     | N<sub>mod</sub>       | See definitions above.                                                         | int   |
| 13     | N<sub>canonical</sub> | See definitions above.                                                         | int   |
| 14     | N<sub>other_mod</sub> | See definitions above.                                                         | int   |
| 15     | N<sub>delete</sub>    | See definitions above.                                                         | int   |
| 16     | N<sub>fail</sub>      | See definitions above.                                                         | int   |
| 17     | N<sub>diff</sub>      | See definitions above.                                                         | int   |
| 18     | N<sub>nocall</sub>    | See definitions above.                                                         | int   |


## Description of columns in `modkit summary`:
### Totals table
The lines of the totals table are prefixed with a `#` character.

| row | name                    | description                                                             | type   |
|-----|-------------------------|-------------------------------------------------------------------------|--------|
| 1   | bases                   | comma-separated list of canonical bases with modification calls.        | str    |
| 2   | total_reads_used        | total number of reads from which base modification calls were extracted | int    |
| 3+  | count_reads_{base}      | total number of reads that contained base modifications for {base}      | int    |
| 4+  | filter_threshold_{base} | filter threshold used for {base}                                        | float  |

### Modification calls table
The modification calls table follows immediately after the totals table.

| column | name       | description                                                                              | type  |
|--------|------------|------------------------------------------------------------------------------------------|-------|
| 1      | base       | canonical base with modification call                                                    | char  |
| 2      | code       | base modification code, or `-` for canonical                                             | char  |
| 3      | pass_count | total number of passing (confidence >= threshold) calls for the modification in column 2 | int   |
| 4      | pass_frac  | fraction of passing (>= threshold) calls for the modification in column 2                | float |
| 5      | all_count  | total number of calls for the modification code in column 2                              | int   |
| 6      | all_frac   | fraction of all calls for the modification in column 2                                   | float |




**Licence and Copyright**

(c) 2023 Oxford Nanopore Technologies Plc.

Modkit is distributed under the terms of the Oxford Nanopore Technologies, Ltd. Public License, v. 1.0.
If a copy of the License was not distributed with this file, You can obtain one at http://nanoporetech.com
