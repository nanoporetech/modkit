# Performance considerations

## Sharding a large modBAM by region.

The `--region` option in `pileup`, `summary`, and `sample-probs` can be used to
operate on a subset of records in a large BAM. If you're working in a
distributed environment, the genome could be sharded into large sections which
are specified to `modkit` in concurrent processes and merged afterward in a
"map-reduce" pattern.

## Setting the `--interval-size` and `--chunk-size`.

Whenever operating on a sorted, indexed BAM, `modkit` will operate in parallel on disjoint spans of the genome.
The length of these spans (i.e. intervals) can be determined by the `--interval-size`  or the `--sampling-interval-size` (for the sampling algorithm only).
The defaults for these parameters works well for genomes such as the human genome usually at CpG contexts.
For smaller genomes with high coverage, you may decide to _decrease_ the interval size in order to take advantage of parallelism.
When using transcriptomes multiple reference transcripts will be grouped - extra parameters are not required.
The `pileup` subcommand also has a `--chunk-size` option that will limit the total number of _intervals_ computed on in parallel.
By default, `modkit` will set this parameter to be 50% larger than the number of threads.
In general, this is a good setting for balancing parallelism and memory usage.
Increasing the `--chunk-size` can increase parallelism (and decrease run time) but will consume more memory.

One note is that using "all-context" base modification models where a base modification call may be at every cytosine, adenine, or both, will require more memory than CpG models.
In these cases, `modkit pileup` must simply hold more base modification records in memory for a given `interval-size`.
It is possible to reduce memory usage (usually at no run-time cost) by decreasing the `--interval-size` parameter.

## Memory usage in `modkit extract`.

Transforming reads into a table with `modkit extract` can produce large files (especially with long reads).
Passing the `--bgzf` flag will produce [compressed](https://www.htslib.org/doc/bgzip.html) (via [gzp](https://docs.rs/gzp/latest/gzp/deflate/struct.Bgzf.html#)) outputs and can dramatically reduce disk usage.

Before the data can be written to disk, however, it is enqueued in memory and can potentially create a large memory burden.
There are a few ways to decrease the amount of memory `modkit extract` will use in these cases:
1. Lower the `--queue-size`, this decreased the number of batches that will be held in flight.
2. Use `--ignore-index` this will force `modkit extract` to run a serial scan of the mod-BAM.
3. Decrease the `--interval-size`, this will decrease the size of the batches.
4. Use `--include-bed` on regions of interest to cull reads overlapping only regions.
5. Use `--motif` and/or `--ignore-implicit` to reduce the number of records written per read.


## Parallelism in `motif search`, `motif refine`, and when to `--skip-search`
The [search algorithm](./intro_find_motifs.md#simple-description-of-the-search-algorithm) takes advantage of parallelism at nearly every step and therefore hugely benefits from running with as many threads as possible (specified with `--threads`).
This horizontal scalability is most easily seen in the secondary search step where (by default) `129536` individual "seed sequences" are evaluated for potential refinement.
If you find that this search is taking a very long time (indicated by the progress bar message "<mod_code> seeds searched") you may consider one of the following:

- Increase the `--exhaustive-seed-min-log-odds` parameter, this will decrease the number of seeds passed on to the refinement step (which is more computationally expensive).
- Decrease the `--exhaustive-seed-len` to 2 or decrease the `--context-size`, this will exponentially decrease the number of seeds to be searched.

You may also decide to run `--skip-search` first and inspect the results.

