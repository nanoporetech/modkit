# Summarizing a modBAM.

The `modkit modbam summary` sub-command is intended for collecting read-level or genome-level statistics from a modBAM.
It can optionally sample a fraction of the reads, sample a target number of reads, and subset to a region.

> [!IMPORTANT]
> Changes for v0.6.0+
> The `--no-sampling` option has been removed, using the entire modBAM is now the default.
The default behavior is to gather a summary of the entire modBAM.
The fastest run time will be to use `--num-reads`, followed by `--sample-frac`.

## Summarize the base modification calls in a modBAM.

```
modkit modbam summary input.bam 
```

will output a table similar to this

```text
> not subsampling, using all reads
> calculating threshold that removes lowest 10 percentile of calls
> collecting base modification calls from unmapped records
# bases                     C,A
# total_reads_used          8340404
# count_reads_C             8340404
# count_reads_A             8340404
# pass_threshold_A          0.7128906
# pass_threshold_C          0.7089844
# modification_codes_for_A  a
# modification_codes_for_C  h,m
 base  code  pass_count   pass_frac     all_count    all_frac
 A     -     38789589231  0.9522225     41493100163  0.91872257
 A     a     1946256088   0.04777748    3670806258   0.08127743
 C     -     26729831806  0.95418394    28992497768  0.9339372
 C     h     113554922    0.0040536085  374670368    0.01206928
 C     m     1169904115   0.041762467   1676137832   0.053993538`
```

> [!IMPORTANT]
> Changes for v0.6.0+
## Summarize only the matched base modification calls

When a reference is passed to `modkit modbam summary` the program can be set to only record base modifications on reads where the base in the read matches the reference base. 
This may sound intuitive, but the alternative is to record base modifications on mismatches, the default (e.g. a 5mC call on a A>C mismatch).

To enable this behavior, pass the `--matched-only` flag.
When this flag is passed, only mapped records are used.
Another option is to use `--motif`, for example `--motif DRACH 2`.
This has the effect of only calculating statistics at certain reference motifs _and_ the base in the read must match.

## Description of columns in `modkit summary`:
### Totals table
The lines of the totals table are prefixed with a `#` character.

| row | name                    | description                                                             | type   |
|-----|-------------------------|-------------------------------------------------------------------------|--------|
| 1   | bases                   | comma-separated list of canonical bases with modification calls.        | str    |
| 2   | total_reads_used        | total number of reads from which base modification calls were extracted | int    |
| 3+  | count_reads_{base}      | total number of reads that contained base modifications for {base}      | int    |
| 4+  | modification_codes_for_{base} | comma-separated list of modification codes found for this base | str| 
| 5+  | filter_threshold_{base} | filter threshold used for {base}                                        | float  |

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


For more details on thresholds see [filtering base modification calls](./filtering.md).


There are  `--filter-percentile`, and `--filter-threshold` options that
can be used with or without sampling.

### Passing a threshold directly.

To estimate the pass thresholds on a subset of reads, but then summarize _all_ of the
reads, there is a two-step process. First, determine the thresholds with `modkit modbam
sample-probs` (see [usage](./advanced_usage.html#sample-probs) for more details). Then run
`modkit modbam summary` with the threshold value specified.

```
modkit modbam sample-probs input.bam [--sampling-frac <frac> | --num-reads <num>]
```

This command will output a table like this:

```
> sampling 10042 reads from BAM
 base  percentile  threshold
 C     10          0.6972656
 C     50          0.96484375
 C     90          0.9941406
```

You can then use pass this threshold directly to `modkit modbam summary`:

```
modkit modbam summary input.bam --filter-threshold 0.6972656
```

## Sampling the modBAM
> [!IMPORTANT]
> Changes for v0.6.0+

Often summarizing the entire modBAM is unnecessary to get a quick answer.
There are 3 ways to sub-sample the modBAM:
1. Use `--region` (with aligned BAM only) to look at only a specific region
2. Use `--num-reads` to sample approximately this many reads from the BAM.
 See the details on [sampling](./sampling.md) for how this algorithm works.
3. Use `--sampling-frac` (optionally with `-seed`).

Using `--num-reads` is often the fastest.

## Performance considerations
When using `--ignore-index` or an unaligned modBAM, memory consumption should be stable and low.
In this case the execution time is bounded by the speed with which the program can read the BAM (set `--io-threads`) and parse the MM tag.
If you have an sorted, indexed modBAM the work can be divided up into intervals and processed in parallel.
Adding more `--threads` will generally speed up execution time proportionally, however it will require more memory.
Each worker ("thread") holds a histogram of probabilities that it uses compiles as it scans over it's section of the modBAM.
There is always a stable number of these histograms initialized at the start of execution, but there will be more of them initialized with more threads.
