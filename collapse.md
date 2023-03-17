# Removing modification calls from BAMs.

There may be multiple possible base modifications to a specific canonical base, for example cytosine has four
known modified variants. As technologies are increasingly able to detect these chemical moieties, not all
downstream tools may be capable of accepting them, or you may want to convert your data to make it comparable
to other data with less specific resolution. To address this issue, `modkit` implements methods for removing
one or more DNA base modification call from a BAM as well as a command line flag to simply combine all
modification calls when performing a pileup to generate a bedMethyl file.


## Removing DNA base modification probabilities.

BAM records containing probabilities for multiple base modifications at each residue can be transformed to
ignore one (or more) of the predicted base modifications. For example, if a BAM file contains 5hmC and 5mC
probabilities, the total probability is necessarily distributed over 5mC $`p_m`$, 5hmC $`p_h`$, and canonical
C, $`p_{C}`$. Described below are the two implemented methods to reduce the dimensionality of the probability
distribution by one or more predicted base modification class.

### Distribute method: `dist`

This is the preferred method.  Take for example a 3-way modification call for C, 5mC, and 5hmC, the
probabilities of each being $`p_{C}`$, $`p_{m}`$, $`p_{h}`$, respectively.  To remove the 5hmC probability,
$`p_{h}`$ is distributed evenly across the other two options. In this example, the updates are $`p_C
\leftarrow p_C + (\frac{p_h}{2})`$ and $`p_m \leftarrow p_m + (\frac{p_h}{2})`$.

### Normalize method: `norm`

Continuing with the above example and notation.  To remove 5hmC calls the updates to the canonical and 5mC
probabilities are $`p_m \leftarrow \frac{p_m}{p_m + p_C}`$ and $`p_C \leftarrow \frac{p_C}{p_m + p_C}`$,
respectively.


## Combining multiple base modifications into a single count.

### Combine: `--combine`

In `modkit pileup` the combine method can be used to produce binary modified/unmodified counts. Continuing
with the above example, the counts for number of modified reads, `N_mod`, becomes the number of reads with
5hmC calls plus the number of reads with 5mC calls. `N_other_mod` will always be `0`. The `raw_mod_code` will
be reported as the ambiguous mod code for the canonical base. See
https://samtools.github.io/hts-specs/SAMtags.pdf for details on base modification codes and `schema.yaml` in
this project for details on the notation for counts.
