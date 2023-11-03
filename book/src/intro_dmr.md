# Perform differential methylation scoring

The `modkit dmr` command contains two subcommands, `pair` and `multi`, that will compare a pair
of samples and multiple samples, respectively. The details of `multi` are the same as `pair` (
it simply does all the pairwise comparisons), so most of the description below will focus on how
to run `pair` and how to interpret the outputs.

## Preparing the input data
The inputs to `modkit dmr` are two or more bedMethyl files (created by `modkit pileup`) that have
been compressed with [bgzip](https://www.htslib.org/doc/bgzip.html) and indexed with 
[tabix](https://www.htslib.org/doc/tabix.html). An example workflow to generate the input data is shown below:

```bash
ref=grch38.fasta
threads=32

norm=normal_sample.bam
norm_pileup=normal_pileup.bed

modkit pileup ${norm} ${norm_pileup} \
  --cpg \
  --ref ${ref} \
  --threads ${threads} \
  --log-filepath log.txt

bgzip ${norm_pileup}
tabix ${norm_pileup}.gz

tumor=tumor_sample.bam
tumor_pileup=tumor_pileup.bed

modkit pileup ${tumor} ${tumor_pileup} \
  --cpg \
  --ref ${ref} \
  --threads ${threads} \
  --log-filepath log.txt 

bgzip ${tumor_pileup}
tabix ${tumor_pileup}.gz
```

## Running differential methylation scoring
Once you have the two (or more) samples to be compared in the appropriate format, the final piece necessary 
is a BED file of the regions to be compared. The `modkit dmr` functionality does not "segment" or otherwise
discover regions, it scores the differences between user-provided regions. To continue with the above example
we can get CpG Islands from the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables). The data may not 
always be appropriate input for `modkit`. For example, the CpG Islands track has extra columns and a header line:

```text
#bin    chrom   chromStart      chromEnd        name       length  cpgNum  gcNum   perCpg  perGc   obsExp
1065    chr20   63004819        63007703        CpG: 272   2884    272     1869    18.9    64.8    0.9
1065    chr20   63009128        63009816        CpG: 48    688     48      432     14      62.8    0.71
```

Therefore we'd need to transform the data with `awk` or similar, such as:
```bash 
awk 'BEGIN{FS="\t"; OFS="\t"} NR>1 {print $2, $3, $4, $5}' cpg_islands_ucsc.bed \
  | bedtools sort -i - >  cpg_islands_ucsc_cleaned.bed
```

Keeping the `name` column is optional. Sorting the regions isn't _strictly_ necessary, the output will
be in the same order as the regions file. Below is an example command to produce the scored output
(continuing from the top example). The `--base` option tells `modkit dmr` which bases to use for scoring
the differences, the argument should be a canonical nucleotide (`A`, `C`, `G`, or `T`) whichever primary 
sequence base has the modifications you're interested in capturing. For example, for CpG islands the base
we're interested in is `C`.

```bash
regions=cpg_islands_ucsc_cleaned.bed
dmr_result=cpg_islands_tumor_normal.bed

modkit dmr pair \
  -a ${norm_pileup}.gz \
  --index-a ${norm_pileup}.gz.tbi \ # optional
  -b ${tumor_pileup}.gz \
  --index-b ${tumor_pileup}.gz.tbi \ # optional
  -o ${dmr_result} \ # output to stdout if not present
  -r ${regions} \
  --ref ${ref} \
  --base C \  # may be repeated if multiple modifications are being used
  --threads ${threads} \
  --log-filepath dmr.log
```

### Running multiple samples
The `modkit dmr multi` command runs all pairwise comparisons for more than two samples.
The preparation of the data is identical to that for `dmr pair` (for each sample, of course). 
An example command could be:
```bash
modkit dmr multi \
  -s ${norm_pileup_1}.gz norm1 \
  -s ${tumor_pileup_1}.gz tumor1 \
  -s ${norm_pileup_2}.gz norm2 \
  -s ${tumor_pileup_2}.gz tumor2 \
  -o ${dmr_dir} \ # required for multi
  -r ${cpg_islands} \
  --ref ${ref} \
  --base C \
  -t 10 \
  -f \
  --log-filepath dmr_multi.log
```

For example the samples could be haplotype-partitioned bedMethyl tables or biological replicates.

## Differential methylation output format
The output from `modkit dmr pair` (and for each pairwise comparison with `modkit dmr multi`) is (roughly)
a BED file with the following schema:

| column | name                         | description                                                                               | type  |
|--------|------------------------------|-------------------------------------------------------------------------------------------|-------|
| 1      | chrom                        | name of reference sequence from bedMethyl input samples                                   | str   |
| 2      | start position               | 0-based start position, from `--regions` argument                                         | int   |
| 3      | end position                 | 0-based exclusive end position, from `--regions` argument                                 | int   |
| 4      | name                         | `name` column from `--regions` BED, or `chr:start-stop` if absent                         | str   |
| 5      | score                        | Difference score, more positive values have increased difference                          | float |
| 6      | sample<sub>a</sub> counts    | Counts of each base modification in the region, comma-separated, for sample A             | str   |
| 7      | sample<sub>a</sub> total     | Total number of base modification calls in the region, including unmodified, for sample A | str   |
| 8      | sample<sub>b</sub> counts    | Counts of each base modification in the region, comma-separated, for sample B             | str   |
| 9      | sample<sub>b</sub> total     | Total number of base modification calls in the region, including unmodified, for sample B | str   |
| 10     | sample<sub>a</sub> fractions | Fraction of calls for each base modification in the region, comma-separated, for sample A | str   |
| 11     | sample<sub>b</sub> fractions | Fraction of calls for each base modification in the region, comma-separated, for sample B | str   |

an example of the output is given below:
```text
chr10   73861   74083   chr10:73861-74083       -0.5007740865394226     h:7,m:18        950     h:8,m:16        802     h:0.74,m:1.89   h:1.00,m:2.00
chr10   74090   74289   chr10:74090-74289       0.5533780473006118      h:8,m:5         936     h:3,m:7         853     h:0.85,m:0.53   h:0.35,m:0.82
chr10   76139   76313   chr10:76139-76313       1.334274110255592       h:6,m:46        507     h:13,m:35       446     h:1.18,m:9.07   h:2.91,m:7.85
```

## Scoring details
The aim of `modkit dmr` is to enable exploratory data analysis of methylation patterns. To that aim, the approach to 
scoring methylation differences is intended to be simple and interpretable. For every region provided, within a sample, 
we model each potentially methylated base as arising from the same distribution. In other words, we discard the relative 
ordering of the base modification calls within a region. We then define a model for the frequency of observing each base 
modification state. In the case of methylated versus unmodified (5mC vs C, or 6mA vs A), we use the binomial distribution: 

\\[
    \mathbf{X}|p \sim \text{Bin}(n, p)
\\]
\\[
    p \sim \text{Beta}(\alpha, \beta)
\\]

where \\(n\\) is the number of potentially methylated bases reported on in the 
region and \\(\mathbf{X}\\) is the vector of counts (canonical and methylated). In the case where there are more than two
states (for example, 5hmC, 5mC, and unmodified C) we use a multinomial distribution: 
\\[
    \mathbf{X}|\pi \sim \text{Mult}(n, \pi)
\\]

\\[
    \pi \sim \text{Dir}(\alpha)
\\]

Let \\(\theta\\) be the maximum a posteriori (MAP) parameters of the model ( \\( \alpha, \beta \\) for the binary case, 
and \\(\alpha \\) in the general case). The `score` reported is the result of a likelihood ratio test:

\\[
\text{score} = \text{log}(\frac{l( \mathbf{X_a} | \theta_{a}) l(\mathbf{X_b} | \theta_{b})}{l(\mathbf{X_{a+b}} | \theta_{a+b})})
\\]

Where \\(\theta_a\\) and \\(\theta_b\\) are the MAP parameters of the model with the two
conditions modeled separately, and \\(\theta_{a+b}\\) are the MLE parameters when the two
conditions are modeled together. For all cases, we use [Jeffrey's prior](https://en.wikipedia.org/wiki/Jeffreys_prior) 
as the prior distribution.
