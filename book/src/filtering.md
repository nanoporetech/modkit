# Partitioning pass and fail base modification calls.

Base modification calls can be removed if they are low confidence based on the predicted modification probabilities. 
In general, `modkit` will estimate a pass confidence threshold value based on the input data. 
Threshold values for modifications on a primary sequence base can be specified on the command line with the `--filter-threshold` option. 
For example to set a threshold for cytosine modifications at 0.8 and adenine modifications at 0.9 provide `--filter-threshold C:0.8 --filter-threshold A:0.9`.
Pass threshold values per base modification can also be specified.
For example, to specify a threshold for canonical adenine at 0.8 and 6mA at 0.9 use `--filter-threshold A:0.8 --mod-thresholds a:0.9`.
Or to specify a threshold of 0.8 for 5mC, 0.9 for 5hmC, and 0.85 for canonical cytosine: `--filter-threshold C:0.85 --mod-thresholds m:0.8 --mod-thresholds h:0.9`

Keep in mind that the `--mod-threshold` option will treat `A`, `C`, `G`, and `T` and "any-mod" as per the [specification](https://samtools.github.io/hts-specs/SAMtags.pdf).

## Further details
1. [Examples of how thresholds affect base modification calls.](./filtering_details.md)
1. [Numerical details of how thresholds are calculated on the fly.](./filtering_numeric_details.md)
