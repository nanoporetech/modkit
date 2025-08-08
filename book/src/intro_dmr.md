# Perform differential methylation scoring

The `modkit dmr` command contains two subcommands, `pair` and `multi`, that will compare pairwise conditions and multiple conditions.
The `pair` command can be used to perform differential methylation detection on single genome positions (for example CpGs) or regions provided as a BED file.
On the other hand, `multi` can only be used to compare regions (such as CpG islands), provided as a BED file.
There are essentially three differential methylation workflows:
1. Perform differential methylation scoring with a pair of samples on regions of the genome.
1. Perform differential methylation scoring across all pairs of samples on regions of the genome.
1. Perform base-level differential modification detection for a pair of conditions.

Each application is explained below. For details on the scoping of these applications see the [limitations](./limitations.md).

## Preparing the input data
The inputs to all `modkit dmr` commands are two or more bedMethyl files (created by `modkit pileup`) that have been compressed with [bgzip](https://www.htslib.org/doc/bgzip.html) and indexed with [tabix](https://www.htslib.org/doc/tabix.html).
An example of how to generate the input data is shown below:

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

bgzip -k ${norm_pileup}
tabix -p bed ${norm_pileup}.gz

# pileup and compression can also be done in one step
tumor=tumor_sample.bam
tumor_pileup=tumor_pileup.bed.gz

modkit pileup ${tumor} - \
  --cpg \
  --ref ${ref} \
  --threads ${threads} \
  --log-filepath log.txt | ${bgzip} -c > ${tumor_pileup}

tabix -p bed ${tumor_pileup}
```

## 1. Perform differential methylation scoring of genomic regions for a pair of samples.
Once you have the two samples to be compared in the appropriate format, the final piece necessary is a BED file of the regions to be compared.
To continue with our example we can get CpG Islands from the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables).
The data may not always be appropriate input for `modkit`.
For example, the CpG Islands track has extra columns and a header line:

```text
#bin  chrom  chromStart  chromEnd  name          length  cpgNum  gcNum  perCpg  perGc  obsExp
660   chr20  9838623     9839213   CpG:  47      590     47     383     15.9   64.9    0.76
661   chr20  10034962    10035266  CpG:  35      304     35     228     23     75      0.85
```

Therefore, we need to transform the data with `awk` or similar, such as:
```bash 
awk 'BEGIN{FS="\t"; OFS="\t"} NR>1 {print $2, $3, $4, $5}' cpg_islands_ucsc.bed \
  | bedtools sort -i - >  cpg_islands_ucsc_cleaned.bed
```

Keeping the `name` column is optional.
Sorting the regions isn't _strictly_ necessary, the output will be in the same order as the regions file.
Below is an example command to produce the scored output.
The `--base` option tells `modkit dmr` which bases to use for scoring the differences, the argument should be a canonical nucleotide (`A`, `C`, `G`, or `T`) whichever primary sequence base has the modifications you're interested in capturing.
For example, for CpG islands the base we're interested in is `C`.

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

The full schema is described [below](#differential-methylation-output-format) with an example output.

## 2. Perform differential methylation detection on all pairs of samples over regions from the genome.
The `modkit dmr multi` command runs all pairwise comparisons for more than two samples for all regions provided in the regions BED file.
The preparation of the data is identical to that for the [previous section](#preparing-the-input-data) (for each sample, of course).

**Note** that if multiple samples are given the same name, they will be combined.

An example command could be:

```bash
modkit dmr multi \
  -s ${norm_pileup_1}.gz norm1 \
  -s ${tumor_pileup_1}.gz tumor1 \
  -s ${norm_pileup_2}.gz norm2 \
  -s ${tumor_pileup_2}.gz tumor2 \
  -o ${dmr_dir} \ # required for multi
  -r ${cpg_islands} \ # skip this option to perform base-level DMR
  --ref ${ref} \
  --base C \
  -t 10 \
  -f \
  --log-filepath dmr_multi.log
```

For example the samples could be haplotype-partitioned bedMethyl tables or biological replicates.
Unlike for `modkit dmr pair` a sample name (e.g. `norm1` and `tumor1` above) must be provided for each input
sample.

## 3. Detecting differential modification at single base positions
The `modkit dmr pair` command has the ability to score individual bases (e.g. differentially methylated CpGs).
To run single-base analysis on one or more paired samples, simply omit the `--regions` (`-r`) option when running `modkit dmr pair`.
When performing single-base analysis the likelihood ratio score and a MAP-based p-value are available.
For details on the likelihood ratio score and the MAP-based p-value, see the [scoring details section](./dmr_scoring_details.md).
For example the above command becomes:

```bash
dmr_result=single_base_haplotype_dmr.bed

modkit dmr pair \
  -a ${hp1_pileup}.gz \
  -b ${hp2_pileup}.gz \
  -o ${dmr_result} \
  --ref ${ref} \
  --base C \
  --threads ${threads} \
  --log-filepath dmr.log
```

Multiple replicates can be provided as well by repeating the `-a` and `-b` options, such as:

```bash
dmr_result=tumor_normal_single_base_replicates.bed

modkit dmr pair \
  -a ${norm_pileup_1}.gz \
  -a ${norm_pileup_2}.gz \
  -b ${tumor_pileup_1}.gz \
  -b ${tumor_pileup_2}.gz \
  -o ${dmr_result_replicates} \
  --ref ${ref} \
  --base C \
  --threads ${threads} \
  --log-filepath dmr.log
```

Keep in mind that the MAP-based p-value provided in single-site analysis is based on a "modified" vs "unmodified" model, see the [scoring section](./dmr_scoring_details.md) and [limitations](./limitations.md) for additional details.

### Note about modification codes
The `modkit dmr` commands require the `--base` option to determine which genome positions to compare, i.e. `--base C` tells `modkit` to compare methylation at cytosine bases.
You may use this option multiple times to compare methylation at multiple primary sequence bases.
It is possible that, during `pileup` a read will have a mismatch and a modification call, such as a C->A mismatch and a 6mA call on that A, and you may not want to use that 6mA call when calculating the differential methylation metrics.
To filter out bedMethyl records like this, `modkit` uses the [SAM specification](https://samtools.github.io/hts-specs/SAMtags.pdf) (page 9) of modification codes to determine which modification codes apply to which primary sequence bases.
For example, `h` is 5hmC and applies to cytosine bases, `a` is 6mA and applies to adenine bases.
However, `modkit pileup` does not require that you use modification codes only in the specification.
If your bedMethyl has records with custom modification codes or codes that aren't in the specification yet, use `--assign-code <mod_code>:<primary_base>` to indicate the code applies to a given primary sequence base.


## Differential methylation output format
The output from `modkit dmr pair` (and for each pairwise comparison with `modkit dmr multi`) is (roughly) a BED file with the following schema:

| column | name                                 | description                                                                                                | type  |
|--------|--------------------------------------|------------------------------------------------------------------------------------------------------------|-------|
| 1      | chrom                                | name of reference sequence from bedMethyl input samples                                                    | str   |
| 2      | start position                       | 0-based start position, from `--regions` argument                                                          | int   |
| 3      | end position                         | 0-based exclusive end position, from `--regions` argument                                                  | int   |
| 4      | name                                 | `name` column from `--regions` BED, or `chr:start-stop` if absent, "." for single sites                    | str   |
| 5      | score                                | difference score, more positive values have increased difference                                           | float |
| 6      | strand                               | strand for the region or single-base position                                                              | str   |
| 7      | sample<sub>a</sub> counts            | counts of each base modification in the region, comma-separated, for sample A                              | str   |
| 8      | sample<sub>a</sub> total             | total number of base modification calls in the region, including unmodified, for sample A                  | int   |
| 9      | sample<sub>b</sub> counts            | counts of each base modification in the region, comma-separated, for sample B                              | str   |
| 10     | sample<sub>b</sub> total             | total number of base modification calls in the region, including unmodified, for sample B                  | int   |
| 11     | sample<sub>a</sub> percents          | percent of calls for each base modification in the region, comma-separated, for sample A                   | str   |
| 12     | sample<sub>b</sub> percents          | percent of calls for each base modification in the region, comma-separated, for sample B                   | str   |
| 13     | sample<sub>a</sub> fraction modified | fraction modification (of any kind) in sample A                                                            | float |
| 14     | sample<sub>b</sub> fraction modified | fraction modification (of any kind) in sample B                                                            | float |
| (15)   | cohen_h                              | Cohen's h [statistic](https://en.wikipedia.org/wiki/Cohen%27s_h) (useful with regions and high-depth runs) | float |
| (16)   | cohen_h_low                          | 95% confidence interval lower bound                                                                        | float |
| (17)   | cohen_h_high                         | 95% confidence interval upper bound                                                                        | float |

an example of the output is given below:

```tsv
#chrom  start     end       name         score                 strand    a_counts    a_total   b_counts    b_total  a_mod_percentages  b_mod_percentages  a_pct_modified  b_pct_modified  effect_size   cohen_h               cohen_h_low           cohen_h_high
chr20   9838623   9839213   CpG:  47     266.4940337502594     .         h:12,m:45   1777      h:40,m:569  2101     h:0.68,m:2.53      h:1.90,m:27.08     0.032076534     0.28986198      -0.25778544   -0.7769052105982954   0.7137374830784644    0.8400729381181264
chr20   10034962  10035266  CpG:  35     5.451198465947527     .         h:7,m:0     1513      h:8,m:8     1352     h:0.46,m:0.00      h:0.59,m:0.59      0.00462657      0.01183432      -0.00720775   -0.08185985199348272  0.008509433779133421  0.155210270207832
chr20   10172120  10172545  CpG:  35     7.8346692400409665    .         h:15,m:28   1229      h:10,m:60   1091     h:1.22,m:2.28      h:0.92,m:5.50      0.034987796     0.06416132      -0.029173527  -0.1358646714824175   0.05433723875835786   0.21739210420647712
chr20   10217487  10218336  CpG:  59     185.80326974114905    .         h:29,m:108  2339      h:38,m:453  1852     h:1.24,m:4.62      h:2.05,m:24.46     0.05857204      0.26511878      -0.20654674   -0.5928872447820701   0.5319235444885967    0.6538509450755435
chr20   10433628  10434345  CpG:  71     -0.06490108458092436  .         h:29,m:2    2750      h:32,m:5    3737     h:1.05,m:0.07      h:0.86,m:0.13      0.0112727275    0.00990099      0.0013717376  0.013409890828678855  -0.03583284167510216  0.06265262333245987
chr20   10671925  10674963  CpG:  255    6.372129285562551     .         h:43,m:24   9461      h:101,m:54  12880    h:0.45,m:0.25      h:0.78,m:0.42      0.0070817037    0.012034162     -0.004952458  -0.05133810350234183  0.024799813076449084  0.07787639392823457
```

**n.b.** Columns 15, 16, and 17 are present when the `--regions` option is passed, but these columns are on the right side of the table when performing single-site analysis (below).
It is generally recommended to use the `--header` flag and standard CSV parsing to make sure the schema's between experiments are maintained.

When performing single-site analysis, the following additional columns are added:

| column | name                       | description                                                                                           | type  |
|--------|----------------------------|-------------------------------------------------------------------------------------------------------|-------|
| 15     | MAP-based p-value          | ratio of the posterior probability of observing the effect size over zero effect size                 | float |
| 16     | effect size                | percent modified in sample A (col 12) minus percent modified in sample B (col 13)                     | float |
| 17     | balanced MAP-based p-value | MAP-based p-value when all replicates are balanced                                                    | float |
| 18     | balanced effect size       | effect size when all replicates are balanced                                                          | float |
| 19     | pct_a_samples              | percent of 'a' samples used in statistical test                                                       | float |
| 20     | pct_b_samples              | percent of 'b' samples used in statistical test                                                       | float |
| 21     | per-replicate p-values     | MAP-based p-values for matched replicate pairs                                                        | float |
| 22     | per-replicate effect sizes | effect sizes matched replicate pairs                                                                  | float |
| 23     | cohen_h                    | Cohen's h [statistic](https://en.wikipedia.org/wiki/Cohen%27s_h) (useful with regions and high-depth runs) | float |
| 24     | cohen_h_low                | 95% confidence interval lower bound                                                                   | float |
| 25     | cohen_h_high               | 95% confidence interval upper bound                                                                   | float |


Columns 16-19 are only produced when multiple samples are provided, columns 20 and 21 are only produced when there is an equal number of 'a' and 'b' samples.
When using multiple samples, it is possible that not every sample will have a modification fraction at a position. 
When this happens, the statistical test is still performed and the values of `pct_a_samples` and `pct_b_samples` reflect the percent of samples from each condition used in the test.

Columns 20 and 21 have the replicate pairwise MAP-based p-values and effect sizes which are calculated based on their order provided on the command line.
For example in the abbreviated command below:

```bash
modkit dmr pair \
  -a ${norm_pileup_1}.gz \
  -a ${norm_pileup_2}.gz \
  -b ${tumor_pileup_1}.gz \
  -b ${tumor_pileup_2}.gz \
  ...
```

Column 20 will contain the MAP-based p-value comparing `norm_pileup_1` versus `tumor_pileup_1` and `norm_pileup_2` versus `norm_pileup_2`.
Column 21 will contain the effect sizes, values are comma-separated.
If you have a different number of samples for each condition, such as:

```bash
modkit dmr pair \
  -a ${norm_pileup_1}.gz \
  -a ${norm_pileup_2}.gz \
  -a ${norm_pileup_3}.gz \
  -b ${tumor_pileup_1}.gz \
  -b ${tumor_pileup_2}.gz \
```
these columns will not be present.

The full table when performing single-site analysis with equal numbers of samples in groups, when running `modkit dmr pair`, will have the following schema:

| column | name                                 | description                                                                                                | type  |
|--------|--------------------------------------|------------------------------------------------------------------------------------------------------------|-------|
| 1      | chrom                                | name of reference sequence from bedMethyl input samples                                                    | str   |
| 2      | start position                       | 0-based start position, from `--regions` argument                                                          | int   |
| 3      | end position                         | 0-based exclusive end position, from `--regions` argument                                                  | int   |
| 4      | name                                 | `name` column from `--regions` BED, or `chr:start-stop` if absent, "." for single sites                    | str   |
| 5      | score                                | difference score, more positive values have increased difference                                           | float |
| 6      | strand                               | strand for the region or single-base position                                                              | str   |
| 7      | sample<sub>a</sub> counts            | counts of each base modification in the region, comma-separated, for sample A                              | str   |
| 8      | sample<sub>a</sub> total             | total number of base modification calls in the region, including unmodified, for sample A                  | int   |
| 9      | sample<sub>b</sub> counts            | counts of each base modification in the region, comma-separated, for sample B                              | str   |
| 10     | sample<sub>b</sub> total             | total number of base modification calls in the region, including unmodified, for sample B                  | int   |
| 11     | sample<sub>a</sub> percents          | percent of calls for each base modification in the region, comma-separated, for sample A                   | str   |
| 12     | sample<sub>b</sub> percents          | percent of calls for each base modification in the region, comma-separated, for sample B                   | str   |
| 13     | sample<sub>a</sub> fraction modified | fraction modification (of any kind) in sample A                                                            | float |
| 14     | sample<sub>b</sub> fraction modified | fraction modification (of any kind) in sample B                                                            | float |
| 15     | MAP-based p-value          | ratio of the posterior probability of observing the effect size over zero effect size                 | float |
| 16     | effect size                | percent modified in sample A (col 12) minus percent modified in sample B (col 13)                     | float |
| 17     | balanced MAP-based p-value | MAP-based p-value when all replicates are balanced                                                    | float |
| 18     | balanced effect size       | effect size when all replicates are balanced                                                          | float |

## Segmenting on differential methylation

When running `modkit dmr` without `--regions` (i.e. [single-site analysis](#3-detecting-differential-modification-at-single-base-positions)) you can generate regions of differential methylation on-the-fly using the segmenting [hidden Markov model](./dmr_scoring_details.html#dmr-segmentation-hidden-markov-model) (HMM).
To run segmenting on the fly, add the `--segment $segments_bed_fp` option to the command such as:


```bash
dmr_result=single_base_haplotype_dmr.bed
dmr_segments=single_base_segements.bed

modkit dmr pair \
  -a ${hp1_pileup}.gz \
  -b ${hp2_pileup}.gz \
  -o ${dmr_result} \
  --segment ${dmr_segments} \ # indicates to run segmentation
  --ref ${ref} \
  --base C \
  --threads ${threads} \
  --log-filepath dmr.log
```

The default settings for the HMM are to run in "coarse-grained" mode which will more eagerly join neighboring sites, potentially at the cost of including sites that are not differentially modified within "Different" blocks.
To activate "fine-grained" mode, pass the `--fine-grained` flag. 

The output schema for the segments is:

| column | name                                 | description                                                                                           | type  |
|--------|--------------------------------------|-------------------------------------------------------------------------------------------------------|-------|
| 1      | chrom                                | name of reference sequence from bedMethyl input samples                                               | str   |
| 2      | start position                       | 0-based start position, from `--regions` argument                                                     | int   |
| 3      | end position                         | 0-based exclusive end position, from `--regions` argument                                             | int   |
| 4      | state-name                           | "different" when sites are differentially modified, "same" otherwise                                  | str   |
| 5      | score                                | difference score, more positive values have increased difference                                      | float |
| 6      | N-sites                              | number of sites (bedmethyl records) in the segment                                                    | float |
| 7      | sample<sub>a</sub> counts            | counts of each base modification in the region, comma-separated, for sample A                         | str   |
| 8      | sample<sub>b</sub> counts            | counts of each base modification in the region, comma-separated, for sample B                         | str   |
| 9      | sample<sub>a</sub> percents          | percent of calls for each base modification in the region, comma-separated, for sample A              | str   |
| 10     | sample<sub>b</sub> percents          | percent of calls for each base modification in the region, comma-separated, for sample B              | str   |
| 11     | sample<sub>a</sub> fraction modified | percent modification (of any kind) in sample A                                                        | float |
| 12     | sample<sub>b</sub> fraction modified | percent modification (of any kind) in sample B                                                        | float |
| 13     | effect size                          | percent modified in sample A (col 11) minus percent modified in sample B (col 12)                     | float |
| 14     | cohen_h                              | Cohen's h [statistic](https://en.wikipedia.org/wiki/Cohen%27s_h) (useful with regions and high-depth runs) | float |
| 15     | cohen_h_low                          | 95% confidence interval lower bound                                                                   | float |
| 16     | cohen_h_high                         | 95% confidence interval upper bound                                                                   | float |

