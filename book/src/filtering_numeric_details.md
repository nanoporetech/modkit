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

# A note on the "probability of modification" (`.`) MM flag.
The `.` flag (e.g. `C+m.`) indicates that primary sequence bases without base modification calls can be inferred to be
canonical ([SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf)). Some base modification callers,
for example [`dorado`](https://github.com/nanoporetech/dorado/) have a default threshold, below which a base modification
probability will be omitted (meaning it will be inferred to be a canonical/unmodified base). In general, omitting
the base modification probabilities when they are very confidently canonical will not change the results from `modkit`.
However, since the base modification probabilities are effectively changed to 100% canonical, can affect the
estimation of the pass threshold.
