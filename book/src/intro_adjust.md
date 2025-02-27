# Updating and Adjusting MM tags.

The `adjust-mods` subcommand can be used to manipulate MM (and corresponding ML) tags in a
modBam. In general, these simple commands are run prior to `pileup`, visualization, or
other analysis. For `adjust-mods` and `update-tags`, if a correct `MN` tag is found, secondary and supplementary
alignments will be output. See [troubleshooting](./troubleshooting.md) for details.


## Ignoring a modification class.

To remove a base modification class from a modBAM and produce a new modBAM, use the
`--ignore` option for `adjust-mods`.

```
modkit adjust-mods input.bam output.adjust.bam --ignore <mod_code_to_ignore>
```
For example the command below will remove 5hmC calls, leaving just 5mC calls.

```
modkit adjust-mods input.bam output.adjust.bam --ignore h
```
For technical details on the transformation see [Removing modification calls from
BAMs](./collapse.md#removing-dna-base-modification-probabilities).

## Combining base modification probabilities.

Combining base modification probabilities may be desirable for downstream analysis or
visualization. Unlike `--ignore` which removes the probability of a class, `--convert`
will sum the probability of one class with another if the second class already exists. For
example, the command below will convert probabilities associated with `h` probability into
`m` probability. If `m` already exists, the probabilities will be summed.  As described in
[changing the modification code](./intro_adjust.md#changing-the-base-modification-code),
if the second base modification code doesn't exist, the probabilities are left unchanged.

```
modkit adjust-mods input.bam output.convert.bam --convert h m
```


## Updating the flag (`?` and `.`).
The [specification](https://samtools.github.io/hts-specs/SAMtags.pdf) (Section 1.7) allows
for omission of the MM flag, however this may not be the intent of missing base
modification probabilities for some models. The command below will add or change the `?` flag to a modBAM.

```
modkit adjust-mods input.bam output.bam --mode ambiguous
```

Another option is to set the flag to `.`, the "implicitly canonical" mode:

```
modkit adjust-mods input.bam output.bam --mode implicit
```

## Changing the base modification code.
Some functions in `modkit` or other tools may require the mod-codes in the MM tag be in
the [specification](https://samtools.github.io/hts-specs/SAMtags.pdf). 

For example, the following command will change `C+Z,` tags to `C+m,` tags.

```
modkit adjust-mods input.bam output.bam --convert Z m
```

## Filtering to specific motifs
If you want to remove base modification calls that don't match a specific basecall sequence motif, you can use the `--motif` in `adjust-mods` retain only base modification calls that match the motif. 
The format for specifying the motif is `<sequence> <offset>` [IUPAC](https://www.bioinformatics.org/sms/iupac.html) codes are allowed in the motif sequence and `<offset>` specifies the 0-based offset into the sequence for the primary base carrying the modification.
For example for CpG dinucleotides `--motif CG 0`.
