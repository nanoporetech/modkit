# Constructing bedMethyl tables.

A primary use of `modkit` is to create summary counts of modified and unmodified bases in an extended [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) format.
bedMethyl files tabulate the counts of base modifications from every sequencing read over each aligned reference genomic position.
In order to create a bedMethyl table, your modBAM must be aligned to a reference genome or transcriptome.
Although not required, providing a reference and using the `--modified-bases` option will provide the clearest results with the best computational performance.
Only **primary alignments** are used in generating the table, it is recommended to mark duplicate alignments before running as multiple primary alignments can be double counted.

## Recommended usage

> [!IMPORTANT]
> Changes for v0.6.0+

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

> [!TIP]
> Note when using a **transcriptome-aligned** modBAM (for direct RNA), pass the `--preload-references` flag to increase performance.

Passing `--modified-bases` not required, but directs Modkit to use optimized routines which will lead to better efficiency.
This option _does_ require a FASTA reference, and will only emit bedmethyl records for base modifications to the primary sequence base in the reference.
For example, a command with the option `--modified-bases 5mC 5hmC 6mA` will _not_ have 6mA records on reference Cytosine locations where a read has a C>A mismatch nor 5mC/5hmC records at reference Adenine locations.
For more details on how this option changes results from previous versions see [migrating to v0.6.0+](./migrating_060.md).

A subset of the base modifications present in the modBAM can specified
For example, passing the option `--modified-bases 5mC` when the modBAM contains 5mC and 5hmC calls will only produce 5mC records.

> [!NOTE]
> If the reference FASTA does not have an index at `path/to/reference.fasta.fai` one will be created.

### Syntax of `--modified-bases`
You may pass the "long name" such as "5mC" or a primary base and the "short name" separated by a `:`.
For example, `--modified-bases 5mC 5hmC 6mA` and `--modified-bases C:m C:h A:a` are equivalent.
ChEBI codes can also be used.
For example, a to make a pileup for a modBAM with direct RNA reads you may use: `--modified-bases A:17596 A:69426 A:a C:m C:19228 G:19229`.

> [!TIP]
> If you don't know which modified bases are present in your modBAM run: `modkit modbam check-tags $bam --head 100`.

> [!IMPORTANT]
> If you have a modBAM with very high sequencing depth (>60,000X) use the `--high-depth` argument. See [limitations](./limitations.md) for details.

### Running without a reference
A reference and `--modified-bases` is not required, for example:

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed --log-filepath pileup.log
```

This command will produce a bedMethyl record for every position in the reference for which there is at least one base modification call (either modified or unmodified).
This invocation may be slower than using `--modified-bases` with a reference.

In both cases, a single file (described [below](#description-of-bedmethyl-output)) with base count summaries will be created.

The program performs best-practices filtering and manipulation of the raw data stored in the input file.
For further details see [filtering modified-base calls](./filtering.md).


## Common options that change how base modifications are tabulated

### Narrowing output to CpG dinucleotides

For user convenience, the counting process can be modulated using several additional transforms and filters.
The most basic of these is to report only counts from reference CpG dinucleotides, this behavior is enabled by passing the `--cpg` flag.
Any motif can be used see [below](#narrowing-output-to-specific-motifs).

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed \
  --cpg \
  --modified-bases 5mC 5hmC \
  --ref path/to/reference.fasta
```

### Specifying `--combine-mods` and/or `--combine-strands` to simplify output

With palindromic motifs, base modification counts can be summed across strands with `--combine-strands`:

```bash
modkit pileup \
  path/to/reads.bam \
  output/path/pileup.bed \
  --modified-bases 5mC 5hmC \
  --cpg \
  --combine-strands \
  --ref path/to/reference.fasta \
```
Note that the strand field of the output will be marked as '.' indicating that the strand information has been lost.

Finally, `--combine-mods` can be used to produce a bedMethyl with binary "modified/not-modified" counts and percentages:

```bash
modkit pileup \
  path/to/reads.bam \
  output/path/pileup.bed \
  --modified-bases 5mC 5hmC \
  --combine-mods \
  --ref path/to/reference.fasta \
  [--cpg] \  # optional
  [--combine-strands] \ #optional
```

As can be seen in the above example, these flags can be combined to generate the following output types (flags in `[braces]` indicate that it is optional):
1. Nothing: All modifications at each, stranded, position.
1. `--cpg`: Narrow output to only CpGs.
1. `--combine-mods [--cpg]`: sum counts of all modification types together.
1. `--combine-strands --cpg [--combine-mods]`: Sum (+)-strand and (-)-strand counts onto the (+)-strand position.

To restrict output to only certain CpGs, pass the `--include-bed` option with the CpGs or regions to be used, see [this page](./intro_include_bed.md) for more details.

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed \
  --modified-bases 5mC 5hmC \
  --cpg \
  --ref path/to/reference.fasta \
  --include-bed path/to/my_cpgs.bed
```

### Narrowing output to specific motifs

By default, `modkit` will output a BED row for all genomic positions where there is at least one base modification for the reference base in the input modBAM.
We define a motif as a short DNA sequence potentially containing [degenerate codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation).
To ease downstream analysis, the `--motif <Motif> <offset, 0-based>` option can be used to pre-filter and annotate the bedMethyl rows.
The `--cpg` flag is a alias for `--motif CG 0` where the sequence motif is `CG` and the offset is `0`, meaning pileup base modification counts for the first `C` in the motif on the top strand the second `C` (complement to `G`) on the bottom strand.
Another example may be `--motif GATC 1`, signaling to pileup counts for the `A` in the second position on the top strand and the `A` in the third position on the bottom strand.

When multiple motifs are specified the `name` column ([column 4](#bedmethyl-column-descriptions)), will indicate which motif the counts are tabulated for.
For example, if `--motif CGCG 2 --motif CG 0` are passed you may see lines such as:

```text
oligo_741_adapters  39 40 m,CG,0   4	-	39	40	255,0,0	4 100.00 4 0 0 0 0 0 0
oligo_741_adapters  39 40 m,CGCG,2 4	-	39	40	255,0,0	4 100.00 4 0 0 0 0 0 0
```

Pileup supports at most eight motifs, including motifs added by `--cpg` or `--modified-bases`.

The `--combine-strands` flag requires exactly one motif, supplied with `--motif` or `--cpg`, and that motif must be reverse-complement palindromic (`CG` _is_ a palindrome but `CHH` is not).
The reverse-strand anchor must not precede the forward-strand anchor; for example, `--motif CGCG 0 --combine-strands` is supported, while `--motif CGCG 2 --combine-strands` is not currently supported.
See [limitations](./limitations.md) for more details.

## Partitioning reads based on phasing information with `--phased`

If you have a modBAM with phased reads containing a `HP` tag. These can be partitioned into separate bedMethyl files on output by passing the `--phased` flag.

```bash
modkit pileup path/to/reads.bam output/directory/ \
  --cpg --modified-bases 5mC 5hmC --ref <reference.fasta> --phased
```
The output will be 3 files: `hp1.bedmethyl`, `hp2.bedmethyl`, and `combined.bedmethyl`.
hp1.bedmethyl and hp2.bedmethyl contain counts for records with `HP=1` and `HP=2` tags, respectively. combined.bedmethyl contains counts for all modBAM records.

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

## bedRMod output

> [!IMPORTANT]
> Changes for v0.6.1

The bedRMod specification includes multiple header lines which all start with the `#` character.
Enabling this output can be done by passing `--bedrmod` on the command line, keep in mind that this flag requires `--modified-bases` and `--ref`.
Although this output specification is intended for direct RNA base modifications, some of the header fields may be generally useful.
The output table below the header is the same with the exception that column 10 (N<sub>valid_cov</sub>), is substituted with total coverage.
Total coverage is calculated as N<sub>valid_cov</sub> + N<sub>fail</sub>.
The `score` column remains the same: N<sub>valid_cov</sub>.

| column | name           | description                              | type |
|--------|----------------|------------------------------------------|------|
| 10     | Total coverage | N<sub>valid_cov</sub> + N<sub>fail</sub> | int  |

More details about bedRMod can be found [here](https://dieterich-lab.github.io/euf-specs/bedRModv2.pdf).


## Performance considerations 
> [!IMPORTANT]
> Changes for v0.6.0

As of version 0.6.0 the efficiency of `pileup` has been greatly improved.
As a result, fewer threads are often required to get high throughput.
For example, increasing threads beyond 8 a 40X human genome BAM will often not yield much benefit.
