# Filtering base modification calls.

Modified base calls are qualitied with a probability that is contained in the ML tag (see the
[specification](https://samtools.github.io/hts-specs/SAMtags.pdf)). We calculate the confidence that the model
has in the base modification prediction as \\(\mathcal{q} = argmax(\textbf{P})\\) where \\(\textbf{P}\\) is the
vector of probabilities for each modification. For example, given a model that can classify canonical
cytosine, 5mC, and 5hmC, \\(\textbf{P}\\) is \\([P_{C}, P_m, P_h]\\), and \\(\mathcal{q}\\) will be \\(\mathcal{q} =
argmax(P_{C}, P_m, P_h)\\), the maximum of the three probabilities.  Filtering in `modkit` is performed by
first determining the value of \\(\mathcal{q}\\) for the lowest n-th percentile of calls (10th percentile by
default).  The threshold value is typically an estimate because the base modification probabilities are
sampled from a subset of the reads. All calls can be used to calculate the exact value, but in practice the
approximation gives the same value. Base modification calls with a confidence value less than this number will
not be counted.  Determination of the threshold value can be performed on the fly (by sampling) or the
threshold value can be specified on the command line with the `--filter_threshold` flag. The `sample-probs`
command can be used to quickly estimate the value of \\(\mathcal{q}\\) at various percentiles.
