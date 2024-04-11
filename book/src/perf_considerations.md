# Performance considerations

## Sharding a large modBAM by region.

The `--region` option in `pileup`, `summary`, and `sample-probs` can be used to
operate on a subset of records in a large BAM. If you're working in a
distributed environment, the genome could be sharded into large sections which
are specified to `modkit` in concurrent processes and merged afterward in a
"map-reduce" pattern.

## Setting the `--interval-size` and `--chunk-size` (`pileup`).

Whenever operating on a sorted, indexed BAM, `modkit` will operate in parallel
on disjoint spans of the genome. The length of these spans (i.e. intervals) can
be determined by the `--interval-size`  or the `--sampling-interval-size` (for
the sampling algorithm only). The defaults for these parameters works well for
genomes such as the human genome. For smaller genomes with high coverage, you
may decide to _decrease_ the interval size in order to take advantage of
parallelism. The `pileup` subcommand also has a `--chunk-size` option that will
limit the total number of _intervals_ computed on in parallel. By default,
`modkit` will set this parameter to be 50% larger than the number of threads.
In general, this is a good setting for balancing parallelism and memory usage.
Increasing the `--chunk-size` can increase parallelism (and decrease run time)
but will consume more memory.

## Memory usage in `modkit extract`.

Transforming reads into a table with `modkit extract` can produce large files (especially with long reads).
Before the data can be written to disk, however, it is enqueued in memory and can potentially create a large memory burden.
There are a few ways to decrease the amount of memory `modkit extract` will use in these cases:
1. Lower the `--queue-size`, this decreased the number of batches that will be held in flight.
2. Use `--ignore-index` this will force `modkit extract` to run a serial scan of the mod-BAM.
3. Decrease the `--interval-size`, this will decrease the size of the batches.
