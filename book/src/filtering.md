# Partitioning pass and fail base modification calls.

Base modification calls can be removed if they are low confidence based on the modification probabilities. In
general, `modkit` will estimate a pass confidence threshold value based on the input data. Threshold values
for modifications on a primary sequence base can be specified on the command line with the
`--filter-threshold` option. For example to set a threshold for cytosine modifications at 0.8 and adenosine
modifications at 0.9 provide `--filter-percentile C:0.8 --filter-percentile A:0.9`. Pass threshold values per
base modification can also be specified. For example, to specify a threshold for canonical adenosine at 0.8
and 6mA at 0.9 use `--mod-thresholds A:0.8 --mod-thresholds a:0.9`.

## Further details
1. [Examples of how thresholds affect base modification calls.](./filtering_details.md)
1. [Numerical details of how thresholds are calculated on the fly.](./filtering_numeric_details.md)
