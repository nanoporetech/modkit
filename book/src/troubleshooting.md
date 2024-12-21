# Troubleshooting

It's recommended to run all `modkit` commands with the `--log-filepath <path-to-file>`
option set. When unexpected outputs are produced inspecting this file will often indicate
the reason.


## Missing secondary and supplementary alignments in output

As of v0.2.4 secondary and supplementary alignments are supported in `adjust-mods`, `update-tags`, `call-mods`, and (optionally) in `extract`.
However, in order to use these alignment records correctly, the `MN` tag must be present and correct in the record. 
The `MN` tag indicates the length of the sequence corresponding to the `MM` and `ML` tags.
As of dorado v0.5.0 the `MN` tag is output when modified base calls are produced.
If the aligner has hard-clipped the sequence, this number will not match the sequence length and the record cannot be used. 
Similarly, if the SEQ field is empty (sequence length zero), the record cannot be used. 
One way to use supplementary alignments is to specify the `-Y` flag when using [dorado](https://github.com/nanoporetech/dorado/) or [minimap2](https://lh3.github.io/minimap2/minimap2.html). 
For these programs, when `-Y` is specified, the sequence will not be hardclipped in supplementary alignments and will be present in secondary alignments. 
Other mapping algorithms that are "MM tag-aware" may allow hard-clipping and update the `MM` and `ML` tags, `modkit` will accept these records as long as the `MN` tag indicates the correct sequence length.

## No rows in `modkit pileup` output.

First, check the logfile, there may be many lines with a variant of

```text
record 905b4cb4-15e2-4060-9c98-866d0aff78bb has un-allowed mode
(ImplicitProbModified), use '--force-allow-implicit' or 'modkit update-tags --mode
ambiguous
```

As suggested, using [update](./intro_adjust.md##updating-the-flag--and-) or setting the
`--force-allow-implicit` flag should produce output from these records.

If the MM tag contains an un-supported mod-code (a [limitation](./limitations.md) that
will be removed in a future release), these errors will be logged. To process the modBAM,
first [convert](./intro_adjust.md#changing-the-base-modification-code) the tags.

In general, if there are MM/ML tags in the modBAM and the bedMethyl is still empty the
log file will contain a line explaining why each record is being skipped.


## Not sampling enough reads to estimate threshold.

If you have previously downsampled the modBAM to a specific region, for example with a
command like
```
samtools view -bh input.sorted.bam chr1
```
modkit will still try and sample reads evenly across the entire genome but fail to find
any aligned to other contigs. To remedy this, either pass the same `--region chr1` option
or set the `--sampling-frac 1.0`. The former will set `modkit` to only use the region
present in the BAM. For example,

```
modkit pileup input.bam output.bed --region chr1
# or
modkit pileup input.bam output.bed --sample-region chr1
```

The latter will sample all of the reads in the BAM, which may slow
down the process, so in general the best advice is to either use the same `--region
<contig>` when using `modkit` or don't downsample the BAM ahead of time (and use `--region
<contig>` in order to get the desired bedMethyl).

The `summary`, `sample-probs`, and `pileup` subcommands all use the same `--region`
option.

## CG positions are missing from output bedMethyl.

When running with `--preset traditional` or `--cpg` the resultant bedMethyl file will only
contains CG positions. However, it will not include positions for which the pass coverage
is zero (see [the column
descriptions](./intro_pileup.md#description-of-bedmethyl-output)). This is to be
expected.
