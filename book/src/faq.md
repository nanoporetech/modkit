# Frequently asked questions

## How are base modification probabilities calculated?

Base modifications are assigned a probability reflecting the confidence the base modification detection algorithm has in making a decision about the modification state of the molecule at a particular position.
The probabilities are parsed from the `ML` tag in the BAM record. These values reflect the probability of the base having a specific modification, `modkit` uses these values and calculates the probability for each modification as well as the probability the base is canonical:

\\[ \
P_{\text{canonical}} = 1 - \sum_{m \in \textbf{M}} P_{m} \
\\]

where \\(\textbf{M}\\) is the set of all of the potential modifications for the base.

For example, consider using a m6A model that predicts m6A or canonical bases at adenine residues, if the \\( P_{\text{m6A}} = 0.9 \\) then the probability of canonical \\( \text{A} \\) is  \\( P_{\text{canonical}} = 1 - P_{\text{m6A}} = 0.1 \\).
Or considering a typical case for cytosine modifications where the model predicts 5hmC, 5mC, and canonical cytosine:

\\[
P_{\text{5mC}} = 0.7, \\\\
P_{\text{5hmC}} = 0.2, \\\\
P_{\text{canonical}} = 1 - P_{\text{5mC}} + P_{\text{5hmC}} = 0.1, \\\\
\\]

A potential confusion is that `modkit` does not assume a base is canonical if the probability of modification is close to \\( \frac{1}{N_{\text{classes}}} \\), the lowest probability the algorithm may assign.

## What value for `--filter-threshold` should I use?

The same way that you may remove low quality data as a first step to any processing, `modkit` will filter out the lowest confidence base modification probabilities.
The filter threshold (or pass threshold) defines the minimum probability required for a read's base modification information at a particular position to be used in a downstream step.
This does not remove the whole read from consideration, just the base modification information attributed to a particular position in the read will be removed.
The most common place to encounter filtering is in `pileup`, where base modification probabilities falling below the pass threshold will be tabulated in the \\( \text{N}\_{\text{Fail}} \\)  column instead of the \\( \text{N}\_{\text{valid}} \\) column.
For highest accuracy, the general recommendation is to let `modkit` estimate this value for you based on the input data.
The value is calculated by first taking a sample of the base modification probabilities from the input dataset and determining the \\(10^{\text{th}}\\) percentile probability value.
This percentile can be changed with the `--filter-percentile` option.
Passing a value to `--filter-threshold` and/or `--mod-threshold` that is higher or lower than the estimated value will have the effect of excluding or including more probabilities, respectively.
It may be a good idea to inspect the distribution of probability values in your data, the `modkit sample-probs` [command](./intro_sample_probs.md) is designed for this task. 
Use the `--hist` and `--out-dir` options to collect a histogram of the prediction probabilities for each canonical base and modification.



<!-- ## How can I perform differential methylation analysis? -->
