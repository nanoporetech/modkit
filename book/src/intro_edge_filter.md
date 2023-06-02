# Removing modification calls at the ends of reads

If you have reads where you know base modifications near the ends should not be used (for example, if they are in adapters), you can use the `--edge-filter <n_basepairs>` option. This option is available in the following subcommands:

1. `pileup`, will ignore base modification calls that are `<n_basepairs>` from the ends.
2. `adjust-mods`, will **remove** base mofification calls that are `<n_basepairs>` from the ends from the resultant output modBAM.
3. `summary`, will ignore base modification calls that are `<n_basepairs>` from the ends.
4. `sample-probs`, will ignore base modification calls that are `<n_basepairs>` from the ends.
5. `call-mods`, will **remove** base mofification calls that are `<n_basepairs>` from the ends from the resultant output modBAM.
6. `extract`, will ignore base modification calls that are `<n_basepairs>` from the ends.

In `pileup` and `call-mods` the edge-filter is also respected when estimating the pass-thresholds.

## Example usages

### call mods with the estimated threshold and ignore modification calls within 100 base pairs of the ends of the reads
```
modkit call-mods <in.bam> <out.bam> --edge-filter 100
```

### Perform pileup, ignoring base modification calls within 100 base pairs of the ends of the reads
```
modkit puleup <in.bam> <out.bed> --edge-filter 100
```
