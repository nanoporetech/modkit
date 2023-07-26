# Repair MM/ML tags on trimmed reads

The `modkit repair` command is useful when you have a BAM with reads where the
canonical sequences have been altered in some way that either renders the `MM`
and `ML` tags invalid (for example, trimmed or hard-clipped) or the data has
been lost completely. This command requires that you have the original base
modification calls for each read you want to repair, and it will project these
base modification calls onto the sequences in the altered BAM.

The command uses two arguments called the "donor" and the "acceptor". The
donor, contains the original, correct, `MM` and `ML` tags and the acceptor is
either missing `MM` and `ML` tags or they are invalid (they will be discarded
either way). The reads in the donor must be a superset of the reads in the
acceptor, meaning you can have extra reads in the donor BAM if some reads have
been removed or filtered earlier in the workflow. Both the donor and the
acceptor must be sorted by read name prior to running `modkit repair`.
Duplicate reads in the acceptor are allowed so long as they have valid SEQ
fields. Lastly, `modkit repair` only works on reads that have been trimmed,
other kinds of alteration such as run-length-encoding are not currently
supported. Split reads, or other derived transformations, are not currently
repairable with this command.

For example a typical workflow may look like this:
```text
# original base modification calls
basecalls_5mC_5hmC.bam

# basecalls that have been trimmed
trimmed.bam # could also be fastq, but would require conversion to BAM

# the two BAM files need to be sorted
samtools -n trimmed.bam -O BAM > trimed_read_sort.bam 
samtools -n basecalls_5mC_5hmC.bam -O BAM > basecalls_5mC_5hmC_read_sort.bam

modkit repair \
    --donor-bam basecalls_5mC_5hmC_read_sort.bam \
    --acceptor-bam trimed_read_sort.bam \
    --log-filepath modkit_repair.log \
    --output-bam trimmed_repaired.bam
```
