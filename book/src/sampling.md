# Details on how reads are sampled with `--num-reads` and `--sample-frac`

## Sampling approximately `--num-reads` from an aligned modBAM

This algorithm is meant be work quickly when the modBAM has generally uniform coverage over the aligned reference.
It works by dividing the reference into "intervals" of length `--interval-size` (or sometimes `--sampling-interval-size` is available as an override) and taking the first \\( \textbf{N} \\) reads from that interval.
The number of reads to take from that interval, \\( \textbf{N} \\), is determined by first proportionally distributing the number of requested reads over the contigs in the reference.
Then for each contig in the reference, the number of reads to be sampled from that contig are divided evenly over the intervals that tile that reference.

Here is a worked example:

contig_1, length = 100
contig_2, length = 50
contig_3, length = 10

Total reference length = 100 + 50 + 10 = 160

Requested reads: 1000 (`--num-reads 1000`)

contig_1 is \\( \frac{100}{160}\\) of the total, so it will get \\( \lceil 1000 * \frac{100}{160} \rceil \\).
A small "fudge factor" of 25% is added.
A table showing how many reads per interval and per contig is printed to the debug log (set with `--log <file_name>`).
The same arithmetic is done for contig_2 and contig_3.

After the number of reads per interval is calculated for each contig, Modkit will sample this many reads from the start of each interval:

```text

contig_1 >------|-----|-----|-----|---<
         ^      ^     ^     ^     ^
N reads are taken from the start of each "^"
```

This algorithm is the "roughest" of the options, however, with a typical whole genome or whole transcriptome alignment provides a robust estimate.

## Sampling `--sampling-frac` fraction of reads.

This is a more traditional option for sampling reads from the modBAM.
In this scenario each reader (more on that below) only keeps a fraction of the records for processing, and discards the others.
A "reader" can be a worker assigned to one interval of the reference (like above) or a single reader assigned to scan the entire modBAM file.
The latter is the case when no index is available (unaligned modBAM) or `--ignore-index` is passed.
The random state in the readers are seeded, so the same command will always sample the same records.
The following will change this determinism:
1. No options or flags.
  Deterministic unless the `--threads` or `--interval-size` changes.
1. Set `--seed`, deterministic regardless of the number of workers (`--threads`), not deterministic when `--interval-size` changes!



