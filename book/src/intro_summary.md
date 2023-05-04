# Summarizing a modBAM.

The `modkit summary` sub-command is intended for collecting read-level statistics on
either a sample of reads, a region, or an entire modBam.

## Summarize the base modification calls in a modBAM.

```
modkit summary input.bam 
```

will output a table similar to this

```bash
> parsing region chr20  # only present if --region option is provided
> sampling 10042 reads from BAM # modulated with --num-reads
> calculating threshold at 10% percentile # modulated with --filter-percentile
> calculated thresholds: C: 0.7167969 # calculated per-canonical base, on the fly
# bases             C
# total_reads_used  9989
# count_reads_C     9989
# pass_threshold_C  0.7167969
# region            chr20:0-64444167
 base  code  pass_count  pass_frac   all_count  all_frac
 C     m     1192533     0.58716166  1305956    0.5790408
 C     h     119937      0.0590528   195335     0.086608544
 C     -     718543      0.3537855   754087     0.33435062
```

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


For more details on thresholds see [filtering base modification calls](./filtering.md).


By default `modkit summary` will only use ten thousand reads when generating the summary
(or fewer if the modBAM has fewer than that). To use all of the reads in the modBAM set
the `--no-sampling` flag.

```
modkit summary input.bam --no-sampling
```

There are `--no-filtering`, `--filter-percentile`, and `--filter-threshold` options that
can be used with or without sampling.

### Passing a threshold directly.

To estimate the pass thresholds on a subset of reads, but then summarize _all_ of the
reads, there is a two-step process. First, determine the thresholds with `modkit
sample-probs` (see [usage](./advanced_usage.html#sample-probs) for more details). Then run
`modkit summary` with the threshold value specified.

```
modkit sample-probs input.bam [--sampling-frac <frac> | --num-reads <num>]
```

This command will output a table like this:

```
> sampling 10042 reads from BAM
 base  percentile  threshold
 C     10          0.6972656
 C     50          0.96484375
 C     90          0.9941406
```

You can then use pass this threshold directly to `modkit summary`:

```
modkit summary input.bam \
    --filter-threshold 0.6972656 \ # filter 10% lowest confidence calls
    --no-sampling
```

