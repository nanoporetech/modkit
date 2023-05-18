# modkit, subcommand documentation

The goal of `modkit` is to enable best-practices manipulation of BAM files containing
modified base information (modBAMs). The various sub-commands and tools available in
`modkit` are described below.  This information can be obtained by invoking the long help
(`--help`) for each command.

> Advanced usage information.


```text
Usage: modkit <COMMAND>

Commands:
  pileup        Tabulates base modification calls across genomic positions. This command
                    produces a bedMethyl formatted file. Schema and description of fields can be
                    found in the README.
  adjust-mods   Performs various operations on BAM files containing base modification
                    information, such as converting base modification codes and ignoring
                    modification calls. Produces a BAM output file.
  update-tags   Renames Mm/Ml to tags to MM/ML. Also allows changing the the mode flag from
                    silent '.' to explicitly '?' or '.'.
  sample-probs  Calculate an estimate of the base modification probability distribution.
  summary       Summarize the mod tags present in a BAM and get basic statistics. The default
                    output is a totals table (designated by '#' lines) and a modification calls
                    table. Descriptions of the columns can be found in the README
  call-mods     Call mods from a modbam, discarding the probabilities.
  motif-bed     Create BED file with all locations of a sequence motif.
                    Example: modkit motif-bed CG 0
  help          Print this message or the help of the given subcommand(s).

Options:
  -h, --help     Print help information.
  -V, --version  Print version information.
```

## pileup
```text
Tabulates base modification calls across genomic positions. This command produces a bedMethyl
formatted file. Schema and description of fields can be found in the README.

Usage: modkit pileup [OPTIONS] <IN_BAM> <OUT_BED>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index available.

  <OUT_BED>
          Output file (or directory with --bedgraph option) to write results into. Specify "-" or
          "stdout" to direct output to stdout.

Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended.

      --region <REGION>
          Process only the specified region of the BAM when performing pileup. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>.

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently.
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size to process concurrently. Smaller interval chunk sizes will use less
          memory but incur more overhead.
          
          [default: 100000]

  -n, --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. Reads will be sampled
          evenly across aligned genome. If a region is specified, either with the --region option or
          the --sample-region option, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand reads is sufficient
          to estimate the model output distribution and determine the filtering threshold.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the filter-percentile. In practice,
          50-100 thousand reads is sufficient to estimate the model output distribution and
          determine the filtering threshold. See filtering.md for details on filtering.

      --seed <SEED>
          Set a random seed for deterministic running, the default is non-deterministic.

      --no-filtering
          Do not perform any filtering, include all mod base calls in output. See filtering.md for
          details on filtering.

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted variant is below
          this confidence percentile. For example, 0.1 will filter out the 10% lowest confidence
          modification calls.
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter threshold can be
          specified with by a decimal number (e.g. 0.75). Per-base thresholds can be specified by
          colon-separated values, for example C:0.75 specifies a threshold value of 0.75 for
          cytosine modification calls. Additional per-base thresholds can be specified by repeating
          the option: for example --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a
          single base option and a default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenosine and 0.9 for
          all other base modification calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and 5mC unless the `--filter-threshold` option is
          also passed. See the online documentation for more details.

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold probability. If
          this option is not provided, but --region is provided, the genomic interval passed to
          --region will be used. Format should be <chrom_name>:<start>-<end> or <chrom_name>.

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size to process concurrently when estimating the threshold probability, can
          be larger than the pileup processing interval.
          
          [default: 1000000]

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --force-allow-implicit
          Force allow implicit-canonical mode. By default modkit does not allow pileup with the
          implicit mode ('.', or silent). The `update-tags` subcommand is provided to update tags to
          the new mode. This option allows the interpretation of implicit mode tags: residues
          without modified base probability will be interpreted as being the non-modified base. We
          do not recommend using this option.

      --cpg
          Only output counts at CpG motifs. Requires a reference sequence to be provided.

      --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format. Required for CpG motif filtering.

  -k, --mask
          Respect soft masking in the reference FASTA.

      --preset <PRESET>
          Optional preset options for specific applications. traditional: Prepares bedMethyl
          analogous to that generated from other technologies for the analysis of 5mC modified
          bases. Shorthand for --cpg --combine-strands --ignore h.
          
          [possible values: traditional]

      --combine-mods
          Combine base modification calls, all counts of modified bases are summed together. See
          collapse.md for details.

      --combine-strands
          When performing CpG analysis, sum the counts from the positive and negative strands into
          the counts for the positive strand.

      --only-tabs
          For bedMethyl output, separate columns with only tabs. The default is to use tabs for the
          first 10 fields and spaces thereafter. The default behavior is more likely to be
          compatible with genome viewers. Enabling this option may make it easier to parse the
          output with tabular data handlers that expect a single kind of separator.

      --bedgraph
          Output bedGraph format, see https://genome.ucsc.edu/goldenPath/help/bedgraph.html. For
          this setting, specify a directory for output files to be make in. Two files for each
          modification will be produced, one for the positive strand and one for the negative
          strand. So for 5mC (m) and 5hmC (h) there will be 4 files produced.

      --prefix <PREFIX>
          Prefix to prepend on bedgraph output file names. Without this option the files will be
          <mod_code>_<strand>.bedgraph.

  -h, --help
          Print help information (use `-h` for a summary).
```

## adjust-mods
```text
Performs various operations on BAM files containing base modification information, such as
converting base modification codes and ignoring modification calls. Produces a BAM output file.

Usage: modkit adjust-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM file to collapse mod call from.
  <OUT_BAM>  File path to new BAM file to be created.

Options:
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path.
      --ignore <IGNORE>              Modified base code to ignore/remove, see
                                     https://samtools.github.io/hts-specs/SAMtags.pdf for details on
                                     the modified base codes. [default: h]
  -t, --threads <THREADS>            Number of threads to use. [default: 4]
  -f, --ff                           Fast fail, stop processing at the first invalid sequence
                                     record. Default behavior is to continue and report
                                     failed/skipped records at the end.
      --convert <CONVERT> <CONVERT>  Convert one mod-tag to another, summing the probabilities
                                     together if the retained mod tag is already present.
  -h, --help                         Print help information.
```

## update-tags
```text
Renames Mm/Ml to tags to MM/ML. Also allows changing the the mode flag from silent '.' to explicitly
'?' or '.'.

Usage: modkit update-tags [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM file to update modified base tags in.
  <OUT_BAM>  File path to new BAM file to be created.

Options:
  -m, --mode <MODE>                  Mode, change mode to this value, options {'ambiguous',
                                     'implicit'}. See spec at:
                                     https://samtools.github.io/hts-specs/SAMtags.pdf. 'ambiguous'
                                     ('?') means residues without explicit modification
                                     probabilities will not be assumed canonical or modified.
                                     'implicit' means residues without explicit modification
                                     probabilities are assumed to be canonical. [possible values:
                                     ambiguous, implicit]
  -t, --threads <THREADS>            Number of threads to use. [default: 4]
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path.
  -h, --help                         Print help information.
```

## sample-probs
```text
Calculate an estimate of the base modification probability distribution.

Usage: modkit sample-probs [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input BAM with modified base tags. If a index is found reads will be sampled evenly across
          the length of the reference sequence.

Options:
  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended.

  -p, --percentiles <PERCENTILES>
          Percentiles to calculate, a space separated list of floats.
          
          [default: 0.1,0.5,0.9]

  -o, --out-dir <OUT_DIR>
          Directory to deposit result tables into. Required for model probability histogram output.
          Creates two files probabilities.tsv and probabilities.txt The .txt contains
          ASCII-histograms and the .tsv contains tab-separated variable data represented by the
          histograms.

      --prefix <PREFIX>
          Label to prefix output files with. E.g. 'foo' will output foo_thresholds.tsv,
          foo_probabilities.tsv, and foo_probabilities.txt.

      --force
          Overwrite results if present.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --hist
          Output histogram of base modification prediction probabilities.

      --buckets <BUCKETS>
          Number of buckets for the histogram, if used.
          
          [default: 128]

  -n, --num-reads <NUM_READS>
          Max number of reads to use, especially recommended when using a large BAM without an
          index. If an indexed BAM is provided, the reads will be sampled evenly over the length of
          the aligned reference. If a region is passed with the --region option, they will be
          sampled over the genomic region.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of reads to sample, for
          example 0.1 will sample 1/10th of the reads.

      --no-sampling
          No sampling, use all of the reads to calculate the filter thresholds.

  -s, --seed <SEED>
          Random seed for deterministic running, the default is non-deterministic.

      --region <REGION>
          Process only the specified region of the BAM when collecting probabilities. Format should
          be <chrom_name>:<start>-<end> or <chrom_name>.

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size to process concurrently. Smaller interval chunk sizes will use less
          memory but incur more overhead. Only used when sampling probs from an indexed bam.
          
          [default: 1000000]

  -h, --help
          Print help information (use `-h` for a summary).
```

## summary
```text
Summarize the mod tags present in a BAM and get basic statistics. The default output is a totals
table (designated by '#' lines) and a modification calls table. Descriptions of the columns can be
found in the README.

Usage: modkit summary [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>  Input ModBam file.

Options:
  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended.

      --tsv
          Output summary as a tab-separated variables stdout instead of a table.

  -n, --num-reads <NUM_READS>
          Max number of reads to use for estimating the filter threshold and generating the summary,
          especially recommended when using a large BAM without an index. If an indexed BAM is
          provided, the reads will be sampled evenly over the length of the aligned reference. If a
          region is passed with the --region option, they will be sampled over the genomic region.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of reads to sample when
          estimating the filter threshold. For example 0.1 will sample 1/10th of the reads.

      --no-sampling
          No sampling, use all of the reads to calculate the filter thresholds and generating the
          summary.

  -s, --seed <SEED>
          Sets a random seed for deterministic running (when using --sample-frac), the default is
          non-deterministic.

      --no-filtering
          Do not perform any filtering, include all base modification calls in the summary. See
          filtering.md for details on filtering.

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted variant is below
          this confidence percentile. For example, 0.1 will filter out the 10% lowest confidence
          base modification calls.
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter threshold can be
          specified with by a decimal number (e.g. 0.75). Per-base thresholds can be specified by
          colon-separated values, for example C:0.75 specifies a threshold value of 0.75 for
          cytosine modification calls. Additional per-base thresholds can be specified by repeating
          the option: for example --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a
          single base option and a default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenosine and 0.9 for
          all other base modification calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and 5mC unless the `--filter-threshold` option is
          also passed. See the online documentation for more details.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --region <REGION>
          Process only the specified region of the BAM when collecting probabilities. Format should
          be <chrom_name>:<start>-<end> or <chrom_name>.

  -i, --interval-size <INTERVAL_SIZE>
          When using regions, interval chunk size to process concurrently. Smaller interval chunk
          sizes will use less memory but incur more overhead.
          
          [default: 1000000]

  -h, --help
          Print help information (use `-h` for a summary).
```

## motif-bed
```text
Create BED file with all locations of a sequence motif. Example: modkit motif-bed CG 0

Usage: modkit motif-bed [OPTIONS] <FASTA> <MOTIF> <OFFSET>

Arguments:
  <FASTA>   Input FASTA file.
  <MOTIF>   Motif to search for within FASTA, e.g. CG.
  <OFFSET>  Offset within motif, e.g. 0.

Options:
  -k, --mask  Respect soft masking in the reference FASTA.
  -h, --help  Print help information.
```

## call-mods
```text
Call mods from a modbam, discarding the probabilities.

Usage: modkit call-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index available.

  <OUT_BAM>
          Output BAM filepath.

Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended.

      --ff
          Fast fail, stop processing at the first invalid sequence record. Default behavior is to
          continue and report failed/skipped records at the end.

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently.
          
          [default: 4]

  -n, --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. If alignments are present
          reads will be sampled evenly across aligned genome. If a region is specified, either with
          the --region option or the --sample-region option, then reads will be sampled evenly
          across the region given. This option is useful for large BAM files. In practice, 10-50
          thousand reads is sufficient to estimate the model output distribution and determine the
          filtering threshold.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the filter-percentile. In practice,
          50-100 thousand reads is sufficient to estimate the model output distribution and
          determine the filtering threshold. See filtering.md for details on filtering.

      --seed <SEED>
          Set a random seed for deterministic running, the default is non-deterministic.

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold probability. If
          this option is not provided, but --region is provided, the genomic interval passed to
          --region will be used. Format should be <chrom_name>:<start>-<end> or <chrom_name>.

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size to process concurrently when estimating the threshold probability, can
          be larger than the pileup processing interval.
          
          [default: 1000000]

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted variant is below
          this confidence percentile. For example, 0.1 will filter out the 10% lowest confidence
          modification calls.
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per primary base. A global filter threshold can
          be specified with by a decimal number (e.g. 0.75). Per-base thresholds can be specified by
          colon-separated values, for example C:0.75 specifies a threshold value of 0.75 for
          cytosine modification calls. Additional per-base thresholds can be specified by repeating
          the option: for example --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a
          single base option and a default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenosine and 0.9 for
          all other base modification calls.

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and 5mC unless the `--filter-threshold` option is
          also passed. See the online documentation for more details.

      --no-filtering
          Don't filter base modification calls, assign each base modification to the highest
          probability prediction.

  -h, --help
          Print help information (use `-h` for a summary).
```
