# Narrow output to specific positions

The `pileup`, `sample-probs`, `summary`, and `extract` sub commands have a `--include-bed` (or `--include-positions`) option that will restrict analysis to only positions that overlap with the intervals contained within the BED file. 

In the case of `pileup`, `summary`, and `sample-probs`, the pass-threshold will be estimated with only base modification probabilities that are aligned to positions overlapping intervals in the BED. In the case of `pileup` and `extract` only positions will be reported if they overlap intervals in the BED.
