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
 base  code  all_count  all_frac     pass_count  pass_frac
 C     h     195335     0.086608544  119937      0.0590528
 C     m     1305956    0.5790408    1192533     0.58716166
 C     -     754087     0.33435062   718543      0.3537855
```

The `pass_count` and `pass_frac` columns are the statistics for calls with confidence
greater than or equal to the `pass_threshold` for that canonical base's calls. For more
details on thresholds see [filtering base modification calls](./filtering.md).
