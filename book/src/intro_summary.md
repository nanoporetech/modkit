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

The `pass_count` and `pass_frac` columns are the statistics for calls with confidence
greater than or equal to the `pass_threshold` for that canonical base's calls. For more
details on thresholds see [filtering base modification calls](./filtering.md).

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

