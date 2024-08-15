# Inspecting base modification probabilities

> For details on how base modification probabilities are calculated, see the [FAQ page](./faq.html#how-are-base-modification-probabilities-calculated)

For most use cases the automatic filtering enabled in `modkit` will produce nearly ideal results.
However, in some cases such as exotic organisms or specialized assays, you may want to interrogate the base modification probabilities directly and tune the pass thresholds.
The `modkit sample-probs` command is designed for this task.
There are two ways to use this command, first by simply running `modkit sample-probs $mod_bam` to get a tab-separated file of threshold values for each modified base. 
This can save time in downstream steps where you wish to re-use the threshold value by passing `--filter-threshold` and skip re-estimating the value.
To generate more advanced output, add `--hist --out-dir $output_dir` to the command and generate per-modification histograms of the output probabilities.
Using the command this way produces 3 files in the `$output_dir`:
1. An HTML document containing a histogram of the total counts of each probability emitted for each modification code (including canonical) in the sampled reads.
1. Another HTML document containing the proportions of each probability emitted.
1. A tab-separated table with the same information as the histograms and the percentile rank of each probability value.

The schema of the table is as follows:

| column | name            | description                                                                                  | type   |
|--------|-----------------|----------------------------------------------------------------------------------------------|--------|
| 1      | code            | modification code or '-' for canonical                                                       | string |
| 2      | primary base    | the primary DNA base for which the code applies                                              | string |
| 3      | range_start     | the inclusive start probability of the bin                                                   | float  |
| 4      | range_end       | the exclusive end probability of the bin                                                     | float  |
| 5      | count           | the total count of probabilities falling in this bin                                         | int    |
| 6      | frac            | the fraction of the total calls for this code/primary base in this bin                       | float  |
| 7      | percentile_rank | the [percentile rank](https://en.wikipedia.org/wiki/Percentile_rank) of this probability bin | float  |

From these plots and tables you can decide on a pass threshold per-modification code and use `--mod-threshold`/`--filter-threshold` [accordingly](./filtering.md).
