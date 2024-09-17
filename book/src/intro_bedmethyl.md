# Constructing bedMethyl tables.

A primary use of `modkit` is to create summary counts of modified and unmodified bases in
an extended [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) format.
bedMethyl files tabulate the counts of base modifications from every sequencing read over
each aligned reference genomic position. In order to create a bedMethyl table, your modBAM
must be aligned to a reference genome. The genome sequence is only required if you are using
the `--cpg` flag or `traditional` preset. Only **primary alignments** are used in generating 
the table, it is recommended to mark duplicate alignments before running as multiple primary
alignments can be double counted (but the behavior is logged). See [limitations](./limitations.md)
for details.

## Basic usage

In its simplest form `modkit pileup` creates a bedMethyl file using the following:

```text
modkit pileup path/to/reads.bam output/path/pileup.bed --log-filepath pileup.log
```

No reference sequence is required. A single file (described
[below](#description-of-bedmethyl-output)) with base count summaries will be created. The
final argument here specifies an optional log file output.

The program performs best-practices filtering and manipulation of the raw data stored in
the input file. For further details see [filtering modified-base calls](./filtering.md).

### Narrowing output to CpG dinucleotides

For user convenience, the counting process can be modulated using several additional
transforms and filters. The most basic of these is to report only counts from reference
CpG dinucleotides. This option requires a reference sequence in order to locate the CpGs
in the reference:

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed --cpg --ref path/to/reference.fasta
```

**Note** that when passing a reference with `--ref` a FASTA index `.fai` file is required to be at `path/to/reference.fasta.fai`.

To restrict output to only certain CpGs, pass the `--include-bed` option with the CpGs to be used, 
see [this page](./intro_include_bed.md) for more details.

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed \
  --cpg \
  --ref path/to/reference.fasta \
  --include-bed path/to/my_cpgs.bed
```

The program also contains preset which combine several options for ease of use. The
`traditional` preset,

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed \
  --ref path/to/reference.fasta \
  --preset traditional
```

performs three transforms:
* restricts output to locations where there is a CG dinucleotide in the reference,
* reports only a C and 5mC counts, using procedures to take into account counts of other
  forms of cytosine modification (notably 5hmC), and
* aggregates data across strands. The strand field of the output will be marked as '.'
  indicating that the strand information has been lost.

Using this option is equivalent to running with the options:

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed --cpg --ref <reference.fasta> --ignore h --combine-strands
```

### Narrowing output to specific motifs

By default, `modkit` will output a BED row for all genomic positions where
there is at least one base modification in the input modBAM. We define a motif
as a short DNA sequence potentially containing [degenerate
codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation). To ease downstream
analysis, the `--motif <Motif> <offset, 0-based>` option can be used to
pre-filter and annotate the bedMethyl rows. The `--cpg` flag is a alias for
`--motif CG 0` where the sequence motif is `CG` and the offset is `0`, meaning
pileup base modification counts for the first `C` in the motif on the top
strand the second `C` (complement to `G`) on the bottom strand. Another example
may be `--motif GATC 1`, signaling to pileup counts for the `A` in the second
position on the top strand and the `A` in the third position on the bottom
strand.

When multiple motifs are specified the `name` column ([column
4](#bedmethyl-column-descriptions)), will indicate which motif the counts are
tabulated for. For example, if `--motif CGCG 2 --motif CG 0` are passed you may
see lines such as:

```text
oligo_741_adapters  39 40 m,CG,0   4	-	39	40	255,0,0	4 100.00 4 0 0 0 0 0 0
oligo_741_adapters  39 40 m,CGCG,2 4	-	39	40	255,0,0	4 100.00 4 0 0 0 0 0 0

```

The `--combine-strands` flag can be combined with `--motif` however all motifs
must be reverse-complement palindromic (`CG` _is_ a palindrome but `CHH` is
not).


### Partitioning reads based on SAM tag values

If have a modBAM with reads from different conditions are other SAM tag annotations (for example `RG` or `HP`) you 
can pass the `--partition-tag` option and `modkit` will output a separate bedMethyl with counts for only the reads 
with that tag value. For example, if you have haplotype-annotated reads with the `HP` tag, you could use a command
like the following:

```bash
modkit pileup path/to/reads.bam output/directory/ --cpg --ref <reference.fasta> --partition-tag HP --prefix haplotyped
```
The output will be multiple files in placed in `output/directory/haplotyped_<1|2|etc>.bed`, multiple `--partition-tag`
options can be passed and the output files will correspond to the observed combinations of tags found in the modBAM. 
For example if `--partition-tag RG` and `--partition-tag HP` are passed:

```bash
outdir/
  <prefix>_<RG_value_1>_<HP_value_1>.bed
  <prefix>_<RG_value_2>_<HP_value_1>.bed
  <prefix>_<RG_value_1>_<HP_value_2>.bed
  <prefix>_<RG_value_2>_<HP_value_2>.bed
  # ... etc
```

Note that only tag values that can be easily turned into strings will be considered valid (e.g. numbers, characters,
strings, etc.), array values will not be used, and will result in `missing` being used. Reads missing all of the 
SAM tags will be put in `ungrouped.bed`.


For more information on the individual options see the [Advanced Usage](./advanced_usage.md) help document.



## Description of bedMethyl output.

Below is a description of the bedMethyl columns generated by `modkit pileup`. A brief description of the
bedMethyl specification can be found on [Encode](https://www.encodeproject.org/data-standards/wgbs/).

### Definitions:

* N<sub>mod</sub> - Number of calls passing filters that were classified as a residue with a specified base modification.
* N<sub>canonical</sub> - Number of calls passing filters were classified as the canonical base rather than modified. The
exact base must be inferred by the modification code. For example, if the modification code is `m` (5mC) then
the canonical base is cytosine. If the modification code is `a`, the canonical base is adenine.
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

### bedMethyl column descriptions.

| column | name                         | description                                                                     | type  |
|--------|------------------------------|---------------------------------------------------------------------------------|-------|
| 1      | chrom                        | name of reference sequence from BAM header                                      | str   |
| 2      | start position               | 0-based start position                                                          | int   |
| 3      | end position                 | 0-based exclusive end position                                                  | int   |
| 4      | modified base code and motif | single letter code for modified base and motif when more than one motif is used | str   |
| 5      | score                        | equal to N<sub>valid_cov</sub>                                                  | int   |
| 6      | strand                       | '+' for positive strand '-' for negative strand, '.' when strands are combined  | str   |
| 7      | start position               | included for compatibility                                                      | int   |
| 8      | end position                 | included for compatibility                                                      | int   |
| 9      | color                        | included for compatibility, always 255,0,0                                      | str   |
| 10     | N<sub>valid_cov</sub>        | see definitions above.                                                          | int   |
| 11     | percent modified             | (N<sub>mod</sub> / N<sub>valid_cov</sub>) * 100                                 | float |
| 12     | N<sub>mod</sub>              | see definitions above                                                           | int   |
| 13     | N<sub>canonical</sub>        | see definitions above                                                           | int   |
| 14     | N<sub>other_mod</sub>        | see definitions above                                                           | int   |
| 15     | N<sub>delete</sub>           | see definitions above                                                           | int   |
| 16     | N<sub>fail</sub>             | see definitions above                                                           | int   |
| 17     | N<sub>diff</sub>             | see definitions above                                                           | int   |
| 18     | N<sub>nocall</sub>           | see definitions above                                                           | int   |

## Performance considerations

The `--interval-size`, `--threads`, `--chunk-size`, and `--max-depth` parameters can be used to tweak the parallelism and memory consumption of `modkit pileup`.
The defaults should be suitable for most use cases, for more details see [performance considerations](./perf_considerations.md) sections.
