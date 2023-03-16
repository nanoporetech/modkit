![Oxford Nanopore Technologies logo](https://github.com/epi2me-labs/modbam2bed/raw/master/images/ONT_logo_590x106.png)

# Modkit

A bioinformatics tool for working with modified bases from Oxford Nanopore. Specifically for converting modBAM
to bedMethyl files using best practices, but also manipulating modBAM files and generating summary statistics.

## Creating a bedMethyl pileup from a modBam

The most typical use case, take a BAM with modified bases (as MM/ML or Mm/Ml tags) and sum the calls from
every read over each genomic position (a pileup). 

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed 
```

No reference sequence is required. A single file (description below) with pileup calls will be created.
Modification filtering will be performed for you.

Some typical options:

1. Only emit counts from reference CpG dinucleotides. This option requires a reference sequence in order to
   locate the CpGs.

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed --cpg --ref path/to/reference.fasta
```
2. Prepare bedMethyl for direct comparison to whole genome bisulfite data.

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed --cpg --ref path/to/reference.fasta --bisulfite
```

## bedMethyl output description

### Some definitions:

**N_mod**: Number of filtered calls that classified a residue as with a specific base modification.  For
example, if the base modification is `h` (5hmC) then this number is the number of filtered reads with a 5hmC
call aligned to this position.

**N_canonical**: Number of filtered calls that classified a residue as canonical as opposed to modified. The
exact base must be inferred by the modification code. For example, if the modification code is `m` (5mC) then
the canonical base is cytosine. If the modification code is `a`, the canonical base is adenosine.

**N_other_mod**: Number of filtered calls that classified a residue as modified where the canonical base is the
same, but the actual modification is different. For example, for a given cytosine there may be 3 reads with
`h` calls, 1 with a canonical call, and 2 with `m` calls. In the row for `h` N_other_mod would be 2 and in the
`m` row N_other_mod would be 3.

**filtered_coverage**: N_mod + N_other_mod + N_canonical, also used as the `score` in the bedMethyl

**N_diff**: Number of reads with a base other than the canonical base for this modification. For example, in a row
for `h` the canonical base is cytosine, if there are 2 reads with C->A substitutions, N_diff will be 2.

**N_delete**: Number of reads with a delete at this position

**N_filtered**: Number of calls where the probability of the call was below the threshold. The threshold can be
set on the command line or computed from the data (usually filtering out the lowest 10th percentile of calls).

**N_nocall**: Number of reads aligned to this position, with the correct canonical base, but without a base
modification call. This can happen, for example, if the model requires a CpG dinucleotide and the read has a
CG->CH substitution.

### Columns in the bedMethy

```yaml
chrom:
  type: str
  description: name of reference sequence from BAM header
start_pos:
  type: int
  description: 0-based index of modified base
end_pos:
  type: int
  description: start_pos + 1
raw_mod_code:
  type: str
  description: single letter code of modified base
score:
  type: int
  description: filtered_coverage
strand:
  type: str
  description: + for positive strand - for negative strand
start_pos:
  type: int
  description: included for compatibility 
end_pos:
  type: int
  description: included for compatibility 
color:
  type: str
  description: included for compatibility, always 255,0,0
filtered_coverage:
  type: int
  description: see definitions
percent_modified:
  type: float
  description: N_mod / filtered_coverage
N_mod:
  type: int
  description: Number of filtered calls for raw_mod_code.
N_canonical:
  type: int
  description: Number of filtered calls for a canonical residue.
N_other_mod:
  type: int
  description: Number of filtered calls for a modification other than raw_mod_code.
N_delete:
  type: int
  description: Number of reads with a delete at this position.
N_filtered:
  type: int
  description: Number of calls that were filtered out.
N_diff:
  type: int
  description: Number of reads with a base other than the canonical base corresponding to raw_mod_code.
N_nocall:
  type: int
  description: Number of reads with no base modification information at this position.
```

## Advanced usage examples

1. Combine multiple base modification calls into one, for example if your data has 5hmC and 5mC
   this will combine the counts into a `C` (any mod) count.

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed --combine-mods
```

2. Combine CpG calls on opposite strands. CpG motifs are reverse complement equivalent, so you may want to
   combine the calls from the positive stand C with the negative strand C (reference G). This operation
   _requires_ that you use the `--cpg` flag and specify a reference sequence.

```bash
modkit pileup path/to/reads.bam output/path/pileup.bed --cpg --ref path/to/reference.fasta \
    --combine-strands <--bisulfite> 
```

3. Produce a bedGraph for each modification in the BAM file file. Counts for the positive and negative strands
   will be put in separate files. Can also be combined with `--cpg` and `--combine-strands` options.

```bash
modkit pileup path/to/reads.bam output/directory/path --bedgraph
```


## Terms and Licence

This is a research release provided under the terms of the Oxford Nanopore Technologies' Public Licence.
Research releases are provided as technology demonstrators to provide early access to features or stimulate
Community development of tools.  Support for this software will be minimal and is only provided directly by
the developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking
and pull requests.  Much as we would like to rectify every issue, the developers may have limited resource for
support of this software.  Research releases may be unstable and subject to rapid change by Oxford Nanopore
Technologies.

© 2023 Oxford Nanopore Technologies Ltd.  Remora is distributed under the terms of the Oxford Nanopore
Technologies' Public Licence.

## Research Release

Research releases are provided as technology demonstrators to provide early access to features or stimulate
Community development of tools. Support for this software will be minimal and is only provided directly by the
developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking and
pull requests. However much as we would like to rectify every issue and piece of feedback users may have, the
developers may have limited resource for support of this software. Research releases may be unstable and
subject to rapid iteration by Oxford Nanopore Technologies.

