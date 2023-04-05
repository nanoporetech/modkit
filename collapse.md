# Augmenting modification calls from BAMs.

There may be multiple possible base modifications to a specific canonical base, for example cytosine has four
known modified variants. As technologies are increasingly able to detect these chemical moieties, not all
downstream tools may be capable of accepting them, or you may want to convert your data to make it comparable
to other data with less specific resolution. To address this issue, `modkit` implements a method for removing
one or more DNA base modifications from a BAM. A command-line flag to simply combine all
modification calls when performing a pileup to generate a bedMethyl file is provided also.


## Manipulating DNA base modification probabilities.

BAM records containing probabilities for multiple base modifications at each residue can be transformed to
ignore one (or more) of the predicted base modifications. For example, if a BAM file contains 5hmC and 5mC
probabilities, the total probability is necessarily distributed over 5mC $p_m$, 5hmC $p_h$, and canonical
C, $p_{C}$.

Take for example a 3-way modification call for C, 5mC, and 5hmC, the probabilities of each being $`p_{C}`$,
$`p_{m}`$, $`p_{h}`$, respectively. To remove the 5hmC probability one method is to distribute evenly the probability mass of $p_{h}$ across the other two bases. In this example, the updates are $p_C \leftarrow p_C + (\frac{p_h}{2})$ and $p_m\leftarrow p_m + (\frac{p_h}{2})$.


### Combining multiple base modifications into a single count.

In `modkit pileup` the `--combine` option can be used to produce binary modified vs. unmodified counts. Continuing
with the above example, the counts for number of modified reads, N<sub>mod</sub>, becomes the number of reads with
5hmC calls plus the number of reads with 5mC calls. N<sub>other_mod</sub> will always be zero. The modified base code (column 4 of the bedMethyl output) will
be reported as the ambiguous modified base code for the canonical base. See
the [SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf) specification for details on base modification codes and `schema.yaml` in
this project for details on the notation for counts.
