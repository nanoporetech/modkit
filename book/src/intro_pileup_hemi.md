# Make hemi-methylation bedMethyl tables with `pileup-hemi`

Base modifications in DNA are inherently single-stranded, they (usually [^1]) don't change the base
pairing of the modified base. However, it may be of interest to know the correspondance
between the methylation state of a single base and another nearby base on the opposite strand - 
on the same molecule. In CpG dinucleotides, this is called "hemi-methylation", when one cytosine
is methylated and the neighbor on the opposite strand is not:

```text
     m
5'GATCGTACA
  CTAGCATGT
      -
```

In the above diagram, the cytosine in the fourth position on the positive strand is methylated (5mC) and the 
cytosine in the fifth position is canonical (-), indicating a "hemi-methylation".

In the case of 5mC and canonical, there are 4 "patterns" of methylation:

```text
m,m (5mC, 5mC)
-,m (canonical, 5mC)
m,- (5mC, canonical)
-,- (canonical, canonical)
```

These are all measured at the _single molecule_ level, meaning each molecule must report on both strands (as
is the case with [duplex](https://www.youtube.com/watch?v=8DVMG7FEBys) reads). For CpGs in the example above the
`MM` tags would be `C+m?` and `G-m?` for the top-strand and bottom-strand cytosines, respectively.

The `modkit pileup-hemi` will perform an aggregation of the methylation "patterns" at genomic positions. An example
command to perform hemi-methylation analysis at CpGs would be

```bash
modkit pileup-hemi \
  /path/to/duplex_reads.bam \
  --cpg \
  -r /path/to/reference.fasta \
  -o hemi_pileup.bed \
  --log modkit.log
```

Many of the `pileup` options are available in `pileup-hemi` with a couple differences: :

1. A motif must be provided. The `--cpg` flag is a preset to aggregate CpG hemi-methylation patterns as shown above.
2. A reference must be provided.
3. Both the positive strand base modification probability and the negative strand base modification probability must be above the pass threshold.

See [Advanced Usage](./advanced_usage.md) for details on all the options.


## Description of hemi-methylation patterns
The `modkit pileup-hemi` command aggregates a pair of base modification calls at each reference motif position
for each double-stranded DNA molecule. The base modification "pattern" indicates the methylation state on each base 
in 5-prime to 3-prime order, using the base modification code to indicate the identity of the base modification and 
`-` to indicate canonical (unmodified). For example `m,-,C` would mean the first base (from the 5' direction) is 5mC 
and the second base is unmodified and the primary base is cytosone. Similarly, `h,m,C` indicates the first base is 
5hmC and the second base is 5mC. The primary base called by the read is included to help disambiguate the unmodified
patterns (`-,-`). All patterns recognized at a location will be reported in the bedMethyl output.

### Definitions:

* N<sub>pattern</sub> - Number of call-pairs passing filters that had the pattern and primary base in column 4. E.g. `m,-,C` 
  indicates the first base in the 5' to 3' direction is 5mC, the second base is unmodified and the primary base in the reads was C.
* N<sub>canonical</sub> - Number of call-pairs passing filters that were classified as unmodified (i.e. the pattern is `-,-`).
* N<sub>other_pattern</sub> - Number of call-pairs passing filters where the pattern is different from the pattern in 
  column 4, but where the primary read base is the same. This count includes the unmodified pattern (`-,-`). **Note** this 
  differs from `pileup` where N<sub>other</sub> does not contain the canonical counts.
* N<sub>valid_cov</sub> - the valid coverage, total number of valid call-pairs.
* N<sub>diff</sub> - Number of reads with a primary base other than the primary base in column 4.
* N<sub>delete</sub> - Number of reads with a deletion at this reference position.
* N<sub>fail</sub> - Number of call-pairs where the probability of the at least one of the calls in the pair was below 
  the pass threshold. The threshold can be set on the command line or computed from the data (usually failing the 
  lowest 10th percentile of calls).
* N<sub>nocall</sub> - Number of reads where either one or both of the base modification calls was not present in the read.

### bedMethyl column descriptions.

| column | name                      | description                                                                                       | type  |
|--------|---------------------------|---------------------------------------------------------------------------------------------------|-------|
| 1      | chrom                     | name of reference sequence from BAM header                                                        | str   |
| 2      | start position            | 0-based start position                                                                            | int   |
| 3      | end position              | 0-based exclusive end position                                                                    | int   |
| 4      | methylation pattern       | comma-separated pair of modification codes `-` means canonical, followed by the primary read base | str   |
| 5      | score                     | Equal to N<sub>valid_cov</sub>.                                                                   | int   |
| 6      | strand                    | always '.' because strand information is combined                                                 | str   |
| 7      | start position            | included for compatibility                                                                        | int   |
| 8      | end position              | included for compatibility                                                                        | int   |
| 9      | color                     | included for compatibility, always 255,0,0                                                        | str   |
| 10     | N<sub>valid_cov</sub>     | See definitions above.                                                                            | int   |
| 11     | fraction modified         | N<sub>pattern</sub> / N<sub>valid_cov</sub>                                                       | float |
| 12     | N<sub>pattern</sub>       | See definitions above.                                                                            | int   |
| 13     | N<sub>canonical</sub>     | See definitions above.                                                                            | int   |
| 14     | N<sub>other_pattern</sub> | See definitions above.                                                                            | int   |
| 15     | N<sub>delete</sub>        | See definitions above.                                                                            | int   |
| 16     | N<sub>fail</sub>          | See definitions above.                                                                            | int   |
| 17     | N<sub>diff</sub>          | See definitions above.                                                                            | int   |
| 18     | N<sub>nocall</sub>        | See definitions above.                                                                            | int   |


## Limitations
1. Only one motif can be used at a time, this limitation may be removed in a later version.
2. Partitioning on tag key:value pairs is not currently supported.

[^1] In biology, there are almost always exceptions to every rule!
