# Making a motif BED file.

Downstream analysis may require a [BED](https://en.wikipedia.org/wiki/BED_(file_format))
file to select motifs of interest. For example, selecting GATC motifs in _E. coli_. This
command requires a reference sequence in
[FASTA](https://en.wikipedia.org/wiki/FASTA_format) a motif to find, which can include
[IUPAC ambiguous bases](https://en.wikipedia.org/wiki/Nucleic_acid_notation) and a
position within the motif.

The following command would make a BED file for CG motifs.

```
modkit motif-bed reference.fasta CG 0 1> cg_modifs.bed
```

The output is directed to standard out.

