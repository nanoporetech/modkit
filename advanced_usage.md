# modkit, subcommand documentation

The goal of `modkit` is to enable easy and correct manipulation of BAM files containing modified base
information (modBAMs). The various sub-commands and tools available in `modkit` are described below.
This information can be obtained by invoking the long help (`--help`) for each command.

> Advanced usage information.


```bash
Usage: modkit <COMMAND>

Commands:
  pileup        Tabulates base modification calls across genomic positions. This command produces a bedMethyl formatted file. Schema and description of fields can be found in the README.
  adjust-mods   Performs various operations on BAM files containing base modification information, such as converting base modification codes and ignoring modification calls. Produces a BAM output file.
  update-tags   Renames Mm/Ml to tags to MM/ML. Also allows changing the the mode flag from silent '.' to '?' or '.'.
  sample-probs  Calculate an estimate of the distribution of mod-base prediction probabilities distribution.
  summary       Summarize the mod tags present in a BAM and get basic statistics.
  motif-bed     Create BED file with all locations of a motif.
  help          Print this message or the help of the given subcommand(s).

Options:
  -h, --help  Print help information
```

## pileup
```bash
Tabulates base modification calls across genomic positions. This command produces a bedMethyl formatted file. Schema and description of fields can be found in the README.

Usage: modkit pileup [OPTIONS] <IN_BAM> <OUT_BED>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index available.

  <OUT_BED>
          Output file to write results into.

Options:
  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently.
          
          [default: 4]

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is recommended.

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size to process concurrently. Smaller interval chunk sizes will use less memory but incur more overhead.
          
          [default: 100000]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the `filter-percentile`. In practice, 50-100 thousand reads is sufficient to estimate the model output distribution and determine the filtering threshold.
          
          [default: 0.1]

      --seed <SEED>
          Set a random seed for deterministic running, the default is non-deterministic.

      --no-filtering
          Do not perform any filtering, include all mod base calls in output. See filtering.md for details on filtering.

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out mod-calls where the probability of the predicted variant is below this percentile. For example, 0.1 will filter out the 10% lowest confidence modification calls.
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Filter threshold, drop calls below this probability.

      --only-tabs
          For bedMethyl output, separate columns with only tabs. The default is to use tabs for the first 10 fields and spaces thereafter. The default behavior is more likely to be compatible with genome viewers. Enabling this option may make it easier to parse the output with tabular data handlers that expect a single kind of separator.

      --bedgraph
          Output bedGraph format, see https://genome.ucsc.edu/goldenPath/help/bedgraph.html. For this setting, specify a directory for output files to be make in. Two files for each modification will be produced, one for the positive strand and one for the negative strand. So for 5mC (m) and 5hmC (h) there will be 4 files produced.

      --force-allow-implicit
          Force allow implicit-canonical mode. By default modkit does not allow pileup with the implicit mode ('.', or silent). The `update-tags` subcommand is provided to update tags to the new mode. This option allows the interpretation of implicit mode tags: residues without modified base probability will be interpreted as being the non-modified base. We do not recommend using this option.

      --region <REGION>
          Process only the specified region of the BAM when performing pileup. Format should be <chrom_name>:<start>-<end>.

      --cpg
          Only output counts at CpG motifs. Requires a reference sequence to be provided.

      --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format. Required for CpG motif filtering.

      --bisulfite
          Prepare data for comparison to whole genome bisulfite sequencing runs. This setting will ignore modification calls other than 5mC, by removing them using the redistribute method described in collapse.md.

      --combine-mods
          Combine mod calls, all counts of modified bases are summed together. See collapse.md for details.

      --combine-strands
          When performing CpG analysis, sum the counts from the positive and negative strands into the counts for the positive strand.

  -h, --help
          Print help information (use `-h` for a summary).
```

## adjust-mods
```bash
Performs various operations on BAM files containing base modification information, such as converting base modification codes and ignoring modification calls. Produces a BAM output file.

Usage: modkit adjust-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM file to collapse mod call from.
  <OUT_BAM>  File path to new BAM file to be created.

Options:
      --ignore <IGNORE>              Modified base code to ignore/remove, see https://samtools.github.io/hts-specs/SAMtags.pdf for details on the modified base codes. [default: h]
  -t, --threads <THREADS>            Number of threads to use. [default: 4]
  -f, --ff                           Fast fail, stop processing at the first invalid sequence record. Default behavior is to continue and report failed/skipped records at the end.
      --method <METHOD>              Method to use to collapse mod calls, 'norm', 'dist'. A full description of the methods can be found in collapse.md. [default: dist]
      --convert <CONVERT> <CONVERT>  Convert one mod-tag to another, summing the probabilities together if the retained mod tag is already present.
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path.
  -h, --help                         Print help information.
```

## update-tags
```bash
Renames Mm/Ml to tags to MM/ML. Also allows changing the the mode flag from silent '.' to '?' or '.'.

Usage: modkit update-tags [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM file to update modified base tags in.
  <OUT_BAM>  File path to new BAM file to be created.

Options:
  -m, --mode <MODE>                  Mode, change mode to this value, options {'ambiguous', 'implicit'}. See spec at: https://samtools.github.io/hts-specs/SAMtags.pdf. 'ambiguous' ('?') means residues without explicit modification probabilities will not be assumed canonical or modified. 'implicit' means residues without explicit modification probabilities are assumed to be canonical [possible values: ambiguous, implicit].
  -t, --threads <THREADS>            Number of threads to use. [default: 4]
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path.
  -h, --help                         Print help information.
```

## sample-probs
```bash
Calculate an estimate of the distribution of mod-base prediction probabilities distribution.

Usage: modkit sample-probs [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>  Input BAM with modified base tags.

Options:
  -f, --sampling-frac <SAMPLING_FRAC>  Fraction of reads to sample, for example 0.1 will sample 1/10th of the reads. [default: 0.1]
  -t, --threads <THREADS>              Number of threads to use reading BAM. [default: 4]
  -s, --seed <SEED>                    Random seed for deterministic running, the default is non-deterministic.
  -p, --percentiles <PERCENTILES>      Percentiles to calculate, a space separated list of floats. [default: 0.1,0.5,0.9]
  -h, --help                           Print help information.
```

## summary
```bash
Summarize the mod tags present in a BAM and get basic statistics.

Usage: modkit summary [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>  Input ModBam file.

Options:
  -t, --threads <THREADS>  Number of threads to use reading BAM. [default: 4]
  -h, --help               Print help information.
```

## motif-bed
```bash
Create BED file with all locations of a motif.

Usage: modkit motif-bed <FASTA> <MOTIF> <OFFSET>

Arguments:
  <FASTA>   Input FASTA file.
  <MOTIF>   Motif to search for within FASTA.
  <OFFSET>  Offset within motif.

Options:
  -h, --help  Print help information.
```
