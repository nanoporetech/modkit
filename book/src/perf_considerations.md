# Performance considerations

## Sharding a large modBAM by region.

The `--region` option in `pileup`, `summary`, and `sample-probs` can be used to operate on
a subset of records in a large BAM. If you're working in a distributed environment, the
genome could be sharded into large sections which are specified to `modkit` in concurrent
processes and merged afterward in a "map-reduce" pattern.

## Setting the `--interval-size` and `--sampling-interval-size`.

Whenever operating on a sorted, indexed BAM, `modkit` will operate in parallel on disjoint
spans of the genome. The length of these spans (i.e. intervals) can be determined by the
`--interval-size` (for general `pileup` and other use cases) or the
`--sampling-interval-size` (for the sampling algorithm only). The defaults for these
parameters works well for genomes such as the human genome. For smaller genomes with high
coverage, you may decide to _decrease_ the interval size in order to take advantage of
parallelism. If you're operating in a memory-constrained environment, _decreasing_ the
interval size will reduce the memory usage per-thread.
