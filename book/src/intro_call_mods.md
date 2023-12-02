# Calling mods in a modBAM

The `call-mods` subcommand in `modkit` transforms one modBAM into another
modBAM where the base modification probabilities have been clamped to 100% and
0%. If the `--filter-threshold` and/or `--mod-threshold`
[options](./advanced_usage.md#call-mods) are provided, base modification calls
failing the threshold will be removed prior to changing the probabilities. The
output modBAM can be used for visualization, `pileup`, or other applications.
If alignment information is present, only the **primary alignment** is used,
and supplementary alignments will not be in the output (see 
[limitations](./limitations.md)).

A modBAM that has been transformed with `call-mods` using `--filter-threshold`
and/or `--mod-threshold` cannot be re-transformed with different thresholds.

Note on `pileup` with clamped probabilities: `modkit pileup` will attempt to
estimate the threshold probability by default, but it is unnecessary if the
modBAM is the result of `call-mods`. The threshold probabilities will be
artificially high (i.e. not representative of the model's output
probabilities). Similarly, specifying `--filter-threshold` and
`--mod-threshold` is not useful because all the ML probabilities have been set
to 0 and 100%.

## Example usages

### Estimate the threshold on the fly, apply to modBAM and clamp the modification calls to certainty.
```
modkit call-mods <in.bam> <out.bam>
```

### Specify a filter threshold for your use-case
```
modkit call-mods <in.bam> <out.bam> --filter-threshold A:0.9 --mod-threshold a:0.95 --filter-threshold C:0.97
```

### Call mods with the estimated threshold and ignore modification calls within 100 base pairs of the ends of the reds
```
modkit call-mods <in.bam> <out.bam> --edge-filter 100
```
