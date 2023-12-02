# Removing modification calls at the ends of reads

If you have reads where you know base modifications near the ends should not be used 
(for example, if they are in adapters), you can use the `--edge-filter <n_basepairs>` option.
Two comma-separated values may be provided to asymmetrically filter out
base modification calls from the start and end of reads. For example, 4,8 will
filter out base modification calls in the first 4 and last 8 bases of the read. One value
will filter symmetrically. 

1. `pileup`, will ignore base modification calls that are `<n_basepairs>` from the ends.
2. `adjust-mods`, will **remove** base modification calls that are `<n_basepairs>` from the ends 
    from the resultant output modBAM.
3. `summary`, will ignore base modification calls that are `<n_basepairs>` from the ends.
4. `sample-probs`, will ignore base modification calls that are `<n_basepairs>` from the ends.
5. `call-mods`, will **remove** base modification calls that are `<n_basepairs>` from the ends from
    the resultant output modBAM.
6. `extract`, will ignore base modification calls that are `<n_basepairs>` from the ends, this also applies
    when making the read-calls table (see [intro to extract](./intro_extract.md)).

In `pileup`, `call-mods`, and `extract` the edge-filter is also respected when estimating the pass-thresholds.
All commands have the flag `--invert-edge-filter` that will _keep_ only base modification probabilities within
`<n_basepairs>` of the ends of the reads.

## Example usages

### call mods with the estimated threshold and ignore modification calls within 100 base pairs of the ends of the reads
```
modkit call-mods <in.bam> <out.bam> --edge-filter 100
```

### Perform pileup, ignoring base modification calls within 100 base pairs of the ends of the reads
```
modkit pileup <in.bam> <out.bed> --edge-filter 100
```
Filter out base modification calls within the first 25 bases or the last 10 bases.
```
modkit pileup <in.bam> <out.bed> --edge-filter 25,10
```

