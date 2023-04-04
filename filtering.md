# Filtering base modification calls

Modified base calls are qualified with a probability that is contained in the ML tag (see the
[specification](https://samtools.github.io/hts-specs/SAMtags.pdf)). We calculate the confidence that the model
has in the base modification prediction as $`\mathcal{q} = argmax(\textbf{P})`$ where $`\textbf{P}`$ is the
vector of probabilities for each modification. For example, given a model that can classify canonical
cytosine, 5mC, and 5hmC, $`\textbf{P}`$ is $`[P_{C}, P_m, P_h]`$, and $`\mathcal{q}`$ will be $`\mathcal{q} =
argmax(P_{C}, P_m, P_h)`$, the maximum of the three probabilities.  Filtering in `modkit` is performed by
first determining the value of $`\mathcal{q}`$ for the lowest n-th percentile of calls (10th percentile by
default).  The threshold value is typically an estimate because the base modification probabilities are
sampled from a subset of the reads. In practice, 10-50 thousand reads are sufficient to make this estimate.
All calls can be used to calculate the exact value, but the approximation usually gives the same value. When
using `pileup`, a region to sample the reads from can be specified with the `--sample-region` option. The
`sample-probs` sub-command is specifically taylored to investigate model confidence values at different
percentiles.

Once a threshold value has been determined base modification calls with a confidence value less than this
number will not be counted.  Determination of the threshold value can be performed on the fly (by sampling,
described above) or the threshold value can be specified on the command line with the `--filter-threshold`
flag.
