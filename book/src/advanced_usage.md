# modkit, subcommand documentation

The goal of `modkit` is to enable best-practices manipulation of BAM files containing
modified base information (modBAMs). The various sub-commands and tools available in
`modkit` are described below.  This information can be obtained by invoking the long help
(`--help`) for each command.

> Advanced usage information.


```text
Modkit is a bioinformatics tool for working with modified bases from Oxford Nanopore

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
                    table. Descriptions of the columns can be found in the README.
  call-mods     Call mods from a modbam, creates a new modbam with probabilities set to 100% if
                    a base modification is called or 0% if called canonical.
  motif-bed     Create BED file with all locations of a sequence motif. 
                    Example: modkit motif-bed CG 0
  extract       Extract read-level base modification information from a modBAM into a
                    tab-separated values table.
  repair        Repair MM and ML tags in one bam with the correct tags from another. To use this
                    command, both modBAMs _must_ be sorted by read name. The "donor" modBAM's reads
                    must be a superset of the acceptor's reads. Extra reads in the donor are
                    allowed, and multiple reads with the same name (secondary, etc.) are allowed in
                    the acceptor. Reads with an empty SEQ field cannot be repaired and will be
                    rejected. Reads where there is an ambiguous alignment of the acceptor to the
                    donor will be rejected (and logged). See the full documentation for details.
  dmr           Perform DMR test on a set of regions. Output a BED file of regions with the
                    score column indicating the magnitude of the difference. Find the schema and
                    description of fields can in the README as well as a description of the model
                    and method. See subcommand help for additional details.
  pileup-hemi   Tabulates double-stranded base modification patters (such as hemi-methylation)
                    across genomic motif positions. This command produces a bedMethyl file, the
                    schema can be found in the online documentation.
  validate      Validate results from a set of mod-BAM files and associated BED files containing
                    the ground truth modified base status at reference positions.
  find-motifs   Find sequence motifs in a bedMethyl pileup that are enriched for base
                    modification.
  entropy       Use a mod-BAM to calculate methylation entropy over genomic windows.
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
          recommended. (alias: log)

      --region <REGION>
          Process only the specified region of the BAM when performing pileup. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>. Commas are allowed.

      --max-depth <MAX_DEPTH>
          Maximum number of records to use when calculating pileup. This argument is passed to the
          pileup engine. If you have high depth data, consider increasing this value substantially.
          Must be less than 2147483647 or an error will be raised.
          
          [default: 8000]

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently.
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller interval chunk sizes
          will use less memory but incur more overhead.
          
          [default: 100000]

      --chunk-size <CHUNK_SIZE>
          Break contigs into chunks containing this many intervals (see `interval_size`). This
          option can be used to help prevent excessive memory usage, usually with no performance
          penalty. By default, modkit will set this value to 1.5x the number of threads specified,
          so if 4 threads are specified the chunk_size will be 6. A warning will be shown if this
          option is less than the number of threads specified.

      --suppress-progress
          Hide the progress bar.

  -n, --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. Reads will be sampled
          evenly across aligned genome. If a region is specified, either with the --region option or
          the --sample-region option, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand reads is sufficient
          to estimate the model output distribution and determine the filtering threshold.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the pass-threshold. In practice, 10-100
          thousand reads is sufficient to estimate the model output distribution and determine the
          filtering threshold. See filtering.md for details on filtering.

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
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenine and 0.9 for all
          other base modification calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless the
          `--filter-threshold` option is also passed. See the online documentation for more details.

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold probability. If
          this option is not provided, but --region is provided, the genomic interval passed to
          --region will be used. Format should be <chrom_name>:<start>-<end> or <chrom_name>.

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently when estimating the threshold
          probability, can be larger than the pileup processing interval.
          
          [default: 1000000]

      --include-bed <INCLUDE_BED>
          BED file that will restrict threshold estimation and pileup results to positions
          overlapping intervals in the file. (alias: include-positions)

      --include-unmapped
          Include unmapped base modifications when estimating the pass threshold.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --force-allow-implicit
          Force allow implicit-canonical mode. By default modkit does not allow pileup with the
          implicit mode (e.g. C+m, no '.' or '?'). The `update-tags` subcommand is provided to
          update tags to the new mode. This option allows the interpretation of implicit mode tags:
          residues without modified base probability will be interpreted as being the non-modified
          base.

      --motif <MOTIF> <MOTIF>
          Output pileup counts for only sequence motifs provided. The first argument should be the
          sequence motif and the second argument is the 0-based offset to the base to pileup base
          modification counts for. For example: --motif CGCG 0 indicates to pileup counts for the
          first C on the top strand and the last C (complement to G) on the bottom strand. The --cpg
          argument is short hand for --motif CG 0.
          
          This argument can be passed multiple times. When more than one motif is used, the
          resulting output BED file will indicate the motif in the "name" field as
          <mod_code>,<motif>,<offset>. For example, given `--motif CGCG 2 --motif CG 0` there will
          be output lines with name fields such as "m,CG,0" and "m,CGCG,2". To use this option with
          `--combine-strands`, all motifs must be reverse-complement palindromic or an error will be
          raised.

      --cpg
          Only output counts at CpG motifs. Requires a reference sequence to be provided.

  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format. Required for motif (e.g. CpG) filtering.

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
          When performing motif analysis (such as CpG), sum the counts from the positive and
          negative strands into the counts for the positive strand position.

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

      --only-tabs
          **Deprecated** The default output has all tab-delimiters. For bedMethyl output, separate
          columns with only tabs. The default is to use tabs for the first 10 fields and spaces
          thereafter. The default behavior is more likely to be compatible with genome viewers.
          Enabling this option may make it easier to parse the output with tabular data handlers
          that expect a single kind of separator.

      --mixed-delim
          Output bedMethyl where the delimiter of columns past column 10 are space-delimited instead
          of tab-delimited. This option can be useful for some browsers and parsers that don't
          expect the extra columns of the bedMethyl format.

      --bedgraph
          Output bedGraph format, see https://genome.ucsc.edu/goldenPath/help/bedgraph.html. For
          this setting, specify a directory for output files to be make in. Two files for each
          modification will be produced, one for the positive strand and one for the negative
          strand. So for 5mC (m) and 5hmC (h) there will be 4 files produced.

      --with-header
          Output a header with the bedMethyl.

      --prefix <PREFIX>
          Prefix to prepend on bedgraph output file names. Without this option the files will be
          <mod_code>_<strand>.bedgraph.

      --partition-tag <PARTITION_TAG>
          Partition output into multiple bedMethyl files based on tag-value pairs. The output will
          be multiple bedMethyl files with the format
          `<prefix>_<tag_value_1>_<tag_value_2>_<tag_value_n>.bed` prefix is optional and set with
          the `--prefix` flag.

  -h, --help
          Print help information (use `-h` for a summary)
```

## adjust-mods
```text
Performs various operations on BAM files containing base modification information, such as
converting base modification codes and ignoring modification calls. Produces a BAM output file

Usage: modkit adjust-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
          BAM file to collapse mod call from. Can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input.

  <OUT_BAM>
          File path to new BAM file to be created. Can be a path to a file or one of `-` or `stdin`
          to specify a stream from standard output.

Options:
      --log-filepath <LOG_FILEPATH>
          Output debug logs to file at this path.

      --ignore <IGNORE>
          Modified base code to ignore/remove, see https://samtools.github.io/hts-specs/SAMtags.pdf
          for details on the modified base codes.

  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

  -f, --ff
          Fast fail, stop processing at the first invalid sequence record. Default behavior is to
          continue and report failed/skipped records at the end.

      --convert <CONVERT> <CONVERT>
          Convert one mod-tag to another, summing the probabilities together if the retained mod tag
          is already present.

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

      --output-sam
          Output SAM format instead of BAM.

      --suppress-progress
          Hide the progress bar

  -h, --help
          Print help information (use `-h` for a summary).
```

## update-tags
```text
Renames Mm/Ml to tags to MM/ML. Also allows changing the the mode flag from silent '.' to explicitly
'?' or '.'.

Usage: modkit update-tags [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM to update modified base tags in. Can be a path to a file or one of `-` or `stdin`
             to specify a stream from standard input.
  <OUT_BAM>  File to new BAM file to be created or one of `-` or `stdin` to specify a stream from
             standard output.

Options:
  -m, --mode <MODE>                  Mode, change mode to this value, options {'explicit',
                                     'implicit'}. See spec at:
                                     https://samtools.github.io/hts-specs/SAMtags.pdf. 'explicit'
                                     ('?') means residues without modification probabilities will
                                     not be assumed canonical or modified. 'implicit' means residues
                                     without explicit modification probabilities are assumed to be
                                     canonical. [possible values: explicit, implicit]
  -t, --threads <THREADS>            Number of threads to use [default: 4]
      --no-implicit-probs            Don't add implicit canonical calls. This flag is important when
                                     converting from one of the implicit modes ( `.` or `""`) to
                                     explicit mode (`?`). By passing this flag, the bases without
                                     associated base modification probabilities will not be assumed
                                     to be canonical. No base modification probability will be
                                     written for these bases, meaning there is no information. The
                                     mode will automatically be set to the explicit mode `?`.
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path.
      --output-sam                   Output SAM format instead of BAM.
  -h, --help                         Print help information.
```

## sample-probs
```text
Calculate an estimate of the base modification probability distribution.

Usage: modkit sample-probs [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input BAM with modified base tags. If a index is found reads will be sampled evenly across
          the length of the reference sequence. Can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input.

Options:
  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended.

      --suppress-progress
          Hide the progress bar.

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

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

      --hist
          Output histogram of base modification prediction probabilities.

      --buckets <BUCKETS>
          Number of buckets for the histogram, if used.
          
          [default: 128]

  -n, --num-reads <NUM_READS>
          Approximate maximum number of reads to use, especially recommended when using a large BAM
          without an index. If an indexed BAM is provided, the reads will be sampled evenly over the
          length of the aligned reference. If a region is passed with the --region option, they will
          be sampled over the genomic region. Actual number of reads used may deviate slightly from
          this number.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of reads to sample, for
          example 0.1 will sample 1/10th of the reads.

      --no-sampling
          No sampling, use all of the reads to calculate the filter thresholds.

  -s, --seed <SEED>
          Random seed for deterministic running, the default is non-deterministic, only used when no
          BAM index is provided.

      --region <REGION>
          Process only the specified region of the BAM when collecting probabilities. Format should
          be <chrom_name>:<start>-<end> or <chrom_name>.

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller interval chunk sizes
          will use less memory but incur more overhead. Only used when sampling probs from an
          indexed bam.
          
          [default: 1000000]

      --include-bed <INCLUDE_BED>
          Only sample base modification probabilities that are aligned to the positions in this BED
          file. (alias: include-positions)

      --only-mapped
          Only use base modification probabilities that are aligned (i.e. ignore soft-clipped, and
          inserted bases).

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
  <IN_BAM>
          Input modBam, can be a path to a file or one of `-` or `stdin` to specify a stream from
          standard input.

Options:
  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended.

      --tsv
          Output summary as a tab-separated variables stdout instead of a table.

      --suppress-progress
          Hide the progress bar.

  -n, --num-reads <NUM_READS>
          Approximate maximum number of reads to use, especially recommended when using a large BAM
          without an index. If an indexed BAM is provided, the reads will be sampled evenly over the
          length of the aligned reference. If a region is passed with the --region option, they will
          be sampled over the genomic region. Actual number of reads used may deviate slightly from
          this number.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of reads to sample when
          estimating the filter threshold. For example 0.1 will sample 1/10th of the reads.

      --no-sampling
          No sampling, use all of the reads to calculate the filter thresholds and generating the
          summary.

  -s, --seed <SEED>
          Sets a random seed for deterministic running (when using --sample-frac), the default is
          non-deterministic, only used when no BAM index is provided.

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
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenine and 0.9 for all
          other base modification calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless the
          `--filter-threshold` option is also passed. See the online documentation for more details.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

      --include-bed <INCLUDE_BED>
          Only summarize base modification probabilities that are aligned to the positions in this
          BED file. (alias: include-positions)

      --only-mapped
          Only use base modification probabilities that are aligned (i.e. ignore soft-clipped, and
          inserted bases).

      --region <REGION>
          Process only the specified region of the BAM when collecting probabilities. Format should
          be <chrom_name>:<start>-<end> or <chrom_name>.

  -i, --interval-size <INTERVAL_SIZE>
          When using regions, interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead.
          
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
Call mods from a modBam, creates a new modBam with probabilities set to 100% if a base modification
is called or 0% if called canonical.

Usage: modkit call-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
          Input BAM, may be sorted and have associated index available. Can be a path to a file or
          one of `-` or `stdin` to specify a stream from standard input.

  <OUT_BAM>
          Output BAM, can be a path to a file or one of `-` or `stdin` to specify a stream from
          standard input.

Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended.

      --ff
          Fast fail, stop processing at the first invalid sequence record. Default behavior is to
          continue and report failed/skipped records at the end.

      --suppress-progress
          Hide the progress bar.

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently.
          
          [default: 4]

  -n, --num-reads <NUM_READS>
          Sample approximately this many reads when estimating the filtering threshold. If
          alignments are present reads will be sampled evenly across aligned genome. If a region is
          specified, either with the --region option or the --sample-region option, then reads will
          be sampled evenly across the region given. This option is useful for large BAM files. In
          practice, 10-50 thousand reads is sufficient to estimate the model output distribution and
          determine the filtering threshold.
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the filter-percentile. In practice,
          50-100 thousand reads is sufficient to estimate the model output distribution and
          determine the filtering threshold. See filtering.md for details on filtering.

      --seed <SEED>
          Set a random seed for deterministic running, the default is non-deterministic, only used
          when no BAM index is provided.

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
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenine and 0.9 for all
          other base modification calls.

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless the
          `--filter-threshold` option is also passed. See the online documentation for more details.

      --no-filtering
          Don't filter base modification calls, assign each base modification to the highest
          probability prediction.

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

      --output-sam
          Output SAM format instead of BAM.

  -h, --help
          Print help information (use `-h` for a summary).
```

## extract
```text
Extract read-level base modification information from a modBAM into a tab-separated values table.

Usage: modkit extract [OPTIONS] <IN_BAM> <OUT_PATH>

Arguments:
  <IN_BAM>
          Path to modBAM file to extract read-level information from, or one of `-` or `stdin` to
          specify a stream from standard input. If a file is used it may be sorted and have
          associated index.

  <OUT_PATH>
          Path to output file, "stdout" or "-" will direct output to standard out. Specifying "null"
          will not output the extract table (useful if all you need is the `--read-calls` output
          table).

Options:
  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

  -q, --queue-size <QUEUE_SIZE>
          Number of reads that can be in memory at a time. Increasing this value will increase
          thread usage, at the cost of memory usage.
          
          [default: 10000]

      --log-filepath <LOG_FILEPATH>
          Path to file to write run log, setting this file is recommended.

      --mapped-only
          Include only mapped bases in output. (alias: mapped)

      --allow-non-primary
          Output aligned secondary and supplementary base modification probabilities as additional
          rows. The primary alignment will have all of the base modification probabilities
          (including soft-clipped ones, unless --mapped-only is used). The non-primary alignments
          will only have mapped bases in the output.

      --num-reads <NUM_READS>
          Number of reads to use. Note that when using a sorted, indexed modBAM that the sampling
          algorithm will attempt to sample records evenly over the length of the reference sequence.
          The result is the final number of records used may be slightly more or less than the
          requested number. When piping from stdin or using a modBAM without an index, the requested
          number of reads will be exact.

      --region <REGION>
          Process only reads that are aligned to a specified region of the BAM. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>.

      --force
          Force overwrite of output file.

      --suppress-progress
          Hide the progress bar.

      --kmer-size <KMER_SIZE>
          Set the query and reference k-mer size (if a reference is provided). Maximum number for
          this value is 50.
          
          [default: 5]

      --ignore-index
          Ignore the BAM index (if it exists) and default to a serial scan of the BAM.

      --read-calls-path <READ_CALLS_PATH>
          Produce a table of read-level base modification calls. This table has, for each read, one
          row for each base modification call in that read using the same thresholding algorithm as
          in pileup, or summary (see online documentation for details on thresholds). Passing this
          option will cause `modkit` to estimate the pass thresholds from the data unless a
          `--filter-threshold` value is passed to the command. Use 'stdout' to stream this table to
          stdout, but note that you cannot stream this table and the raw extract table to stdout.

      --pass-only
          Only output base modification calls that pass the minimum confidence threshold. (alias:
          pass)

      --no-headers
          Don't print the header lines in the output tables.

      --reference <REFERENCE>
          Path to reference FASTA to extract reference context information from. If no reference is
          provided, `ref_kmer` column will be "." in the output. (alias: ref)

      --include-bed <INCLUDE_BED>
          BED file with regions to include (alias: include-positions). Implicitly only includes
          mapped sites

  -v, --exclude-bed <EXCLUDE_BED>
          BED file with regions to _exclude_ (alias: exclude).

      --motif <MOTIF> <MOTIF>
          Output read-level base modification probabilities restricted to the reference sequence
          motifs provided. The first argument should be the sequence motif and the second argument
          is the 0-based offset to the base to pileup base modification counts for. For example:
          --motif CGCG 0 indicates include base modifications for which the read is aligned to the
          first C on the top strand and the last C (complement to G) on the bottom strand. The --cpg
          argument is short hand for --motif CG 0. This argument can be passed multiple times.

      --cpg
          Only output counts at CpG motifs. Requires a reference sequence to be provided.

  -k, --mask
          When using motifs, respect soft masking in the reference sequence.

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter threshold can be
          specified with by a decimal number (e.g. 0.75). Per-base thresholds can be specified by
          colon-separated values, for example C:0.75 specifies a threshold value of 0.75 for
          cytosine modification calls. Additional per-base thresholds can be specified by repeating
          the option: for example --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a
          single base option and a default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenine and 0.9 for all
          other base modification calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless the
          `--filter-threshold` option is also passed. See the online documentation for more details.

      --no-filtering
          Do not perform any filtering, include all mod base calls in output. See filtering.md for
          details on filtering.

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently when estimating the threshold
          probability.
          
          [default: 1000000]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the pass-threshold. In practice, 10-100
          thousand reads is sufficient to estimate the model output distribution and determine the
          filtering threshold. See filtering.md for details on filtering.

  -n, --sample-num-reads <SAMPLE_NUM_READS>
          Sample this many reads when estimating the filtering threshold. If a sorted, indexed
          modBAM is provided reads will be sampled evenly across aligned genome. If a region is
          specified, with the --region, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand reads is sufficient
          to estimate the model output distribution and determine the filtering threshold.
          
          [default: 10042]

      --seed <SEED>
          Set a random seed for deterministic running, the default is non-deterministic.

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted variant is below
          this confidence percentile. For example, 0.1 will filter out the 10% lowest confidence
          modification calls.
          
          [default: 0.1]

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller interval chunk sizes
          will use less memory but incur more overhead. Only used when an indexed modBam is provided.
          
          [default: 100000]

      --ignore-implicit
          Ignore implicitly canonical base modification calls. When the `.` flag is used in the MM
          tag, this implies that bases missing a base modification probability are to be assumed
          canonical. Set this flag to omit those base modifications from the output. For additional
          details see the SAM spec: https://samtools.github.io/hts-specs/SAMtags.pdf.

  -h, --help
          Print help information (use `-h` for a summary).
```

## repair
```text
Repair MM and ML tags in one bam with the correct tags from another. To use this command, both
modBAMs _must_ be sorted by read name. The "donor" modBAM's reads must be a superset of the
acceptor's reads. Extra reads in the donor are allowed, and multiple reads with the same name
(secondary, etc.) are allowed in the acceptor. Reads with an empty SEQ field cannot be repaired and
will be rejected. Reads where there is an ambiguous alignment of the acceptor to the donor will be
rejected (and logged). See the full documentation for details.

Usage: modkit repair [OPTIONS] --donor-bam <DONOR_BAM> --acceptor-bam <ACCEPTOR_BAM> --output-bam <OUTPUT_BAM>

Options:
  -d, --donor-bam <DONOR_BAM>        Donor modBAM with original MM/ML tags. Must be sorted by read
                                     name.
  -a, --acceptor-bam <ACCEPTOR_BAM>  Acceptor modBAM with reads to have MM/ML base modification data
                                     projected on to. Must be sorted by read name.
  -o, --output-bam <OUTPUT_BAM>      output modBAM location.
      --log-filepath <LOG_FILEPATH>  File to write logs to, it is recommended to use this option as
                                     some reads may be rejected and logged here
  -t, --threads <THREADS>            The number of threads to use [default: 4]
  -h, --help                         Print help information.
```

## validate
```text
Validate results from a set of mod-BAM files and associated BED files containing the ground truth
modified base status at reference positions.

Usage: modkit validate [OPTIONS]

Options:
      --bam-and-bed <BAM> <BED>
          Argument accepts 2 values. The first value is the BAM file path with modified base tags.
          The second is a bed file with ground truth reference positions. The name field in the
          ground truth bed file should be the short name (single letter code or ChEBI ID) for a
          modified base or `-` to specify a canonical base ground truth position. This argument can
          be provided more than once for multiple samples.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

  -c, --canonical-base <CANONICAL_BASE>
          Canonical base to evaluate. By default, this will be derived from mod codes in ground
          truth BED files. For ground truth with only canonical sites and/or ChEBI codes this values
          must be set.
          
          [possible values: A, C, G, T]

      --min-identity <MIN_ALIGNMENT_IDENTITY>
          Only use reads with alignment identity >= this number, in Q-space (phred score).

      --min-length <MIN_ALIGNMENT_LENGTH>
          Remove reads with fewer aligned reference bases than this threshold.

  -q, --filter-quantile <FILTER_QUANTILE>
          Filter out modified base calls where the probability of the predicted variant is below
          this confidence percentile. For example, 0.1 will filter out the 10% lowest confidence
          modification calls.
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify modified base probability filter threshold value. If specified, --filter-threshold
          will override --filter-quantile.

  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

      --suppress-progress
          Hide the progress bar.

  -o, --out-filepath <OUT_FILEPATH>
          Specify a file for machine parseable output.

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended. (alias: log)

  -h, --help
          Print help information (use `-h` for a summary).
```

## pileup-hemi
```text
Tabulates double-stranded base modification patters (such as hemi-methylation) across genomic motif
positions. This command produces a bedMethyl file, the schema can be found in the online
documentation.

Usage: modkit pileup-hemi [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index available.

Options:
  -o, --out-bed <OUT_BED>
          Output file to write results into. Will write to stdout if not provided.

      --cpg
          Aggregate double-stranded base modifications for CpG dinucleotides. This flag is
          short-hand for --motif CG 0.

      --motif <MOTIF> <MOTIF>
          Specify the sequence motif to pileup double-stranded base modification pattern counts for.
          The first argument should be the sequence motif and the second argument is the 0-based
          offset to the base to pileup base modification counts for. For example: --motif CG 0
          indicates to generate pattern counts for the C on the top strand and the following C
          (opposite to G) on the negative strand. The motif must be reverse-complement palindromic
          or an error will be raised. See the documentation for more examples and details.

  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format.

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended. (alias: log)

      --region <REGION>
          Process only the specified region of the BAM when performing pileup. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>. Commas are allowed.

      --max-depth <MAX_DEPTH>
          Maximum number of records to use when calculating pileup. This argument is passed to the
          pileup engine. If you have high depth data, consider increasing this value substantially.
          Must be less than 2147483647 or an error will be raised.
          
          [default: 8000]

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently.
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller interval chunk sizes
          will use less memory but incur more overhead.
          
          [default: 100000]

      --chunk-size <CHUNK_SIZE>
          Break contigs into chunks containing this many intervals (see `interval_size`). This
          option can be used to help prevent excessive memory usage, usually with no performance
          penalty. By default, modkit will set this value to 1.5x the number of threads specified,
          so if 4 threads are specified the chunk_size will be 6. A warning will be shown if this
          option is less than the number of threads specified.

      --suppress-progress
          Hide the progress bar.

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
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenine and 0.9 for all
          other base modification calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless the
          `--filter-threshold` option is also passed. See the online documentation for more details.

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold probability. If
          this option is not provided, but --region is provided, the genomic interval passed to
          --region will be used. Format should be <chrom_name>:<start>-<end> or <chrom_name>.

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently when estimating the threshold
          probability, can be larger than the pileup processing interval.
          
          [default: 1000000]

      --include-bed <INCLUDE_BED>
          BED file that will restrict threshold estimation and pileup results to positions
          overlapping intervals in the file. (alias: include-positions)

      --include-unmapped
          Include unmapped base modifications when estimating the pass threshold.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --force-allow-implicit
          Force allow implicit-canonical mode. By default modkit does not allow pileup with the
          implicit mode (e.g. C+m, no '.' or '?'). The `update-tags` subcommand is provided to
          update tags to the new mode. This option allows the interpretation of implicit mode tags:
          residues without modified base probability will be interpreted as being the non-modified
          base.

  -k, --mask
          Respect soft masking in the reference FASTA.

      --combine-mods
          Combine base modification calls, all counts of modified bases are summed together. See
          collapse.md for details.

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. Two comma-separated values may be provided to asymmetrically filter out base
          modification calls from the start and end of the reads. For example, 4,8 will filter out
          base modification calls in the first 4 and last 8 bases of the read.

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification calls at the ends of
          reads, only _keep_ base modification calls at the ends of reads. E.g. if usually, "4,8"
          would remove (i.e. filter out) base modification calls in the first 4 and last 8 bases of
          the read, using this flag will keep only base modification calls in the first 4 and last 8
          bases.

      --only-tabs
          **Deprecated** The default output has all tab-delimiters. For bedMethyl output, separate
          columns with only tabs. The default is to use tabs for the first 10 fields and spaces
          thereafter. The default behavior is more likely to be compatible with genome viewers.
          Enabling this option may make it easier to parse the output with tabular data handlers
          that expect a single kind of separator.

      --mixed-delim
          Output bedMethyl where the delimiter of columns past column 10 are space-delimited instead
          of tab-delimited. This option can be useful for some browsers and parsers that don't
          expect the extra columns of the bedMethyl format.

  -h, --help
          Print help information (use `-h` for a summary).
```

## find-motifs
```text
Find sequence motifs in a bedMethyl pileup that are enriched for base modification

Usage: modkit find-motifs [OPTIONS] --in-bedmethyl <IN_BEDMETHYL> --ref <REFERENCE_FASTA>

Options:
  -i, --in-bedmethyl <IN_BEDMETHYL>
          Input bedmethyl table, can be used directly from modkit pileup.
  -o, --out-table <OUT_TABLE>
          Optionally output a machine-parsable TSV (human-readable table will always be output to
          the log).
      --known-motifs-table <OUT_KNOWN_TABLE>
          Optionally output machine parsable table with known motif modification frequencies that
          were not found during search.
  -t, --threads <THREADS>
          Number of threads to use. [default: 4]
  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format used for the pileup. Reference sequence in FASTA format
          used for the pileup.
      --low-thresh <LOW_THRESHOLD>
          Fraction modified threshold below which consider a genome location to be "low
          modification". [default: 0.2]
      --high-thresh <HIGH_THRESHOLD>
          Fraction modified threshold above which consider a genome location to be "high
          modification" or enriched for modification. [default: 0.6]
      --min-log-odds <MIN_LOG_ODDS>
          Minimum log-odds to consider a motif sequence to be enriched. [default: 1.5]
      --exhaustive-seed-min-log-odds <EXHAUSTIVE_SEED_MIN_LOG_ODDS>
          Minimum log-odds to consider a motif sequence to be enriched when performing exhaustive
          search. [default: 2.5]
      --exhaustive-seed-len <EXHAUSTIVE_SEED_LEN>
          Exhaustive search seed length, increasing this value increases computational time. [default:
          3]
      --skip-search
          Skip the exhaustive search phase, saves time but the results may be less sensitive.
      --min-coverage <MIN_COVERAGE>
          Minimum coverage in the bedMethyl to consider a record valid. [default: 5]
      --context-size <CONTEXT_SIZE> <CONTEXT_SIZE>
          Upstream and downstream number of bases to search for a motif sequence around a modified
          base. Example: --context-size 12 12 [default: 12 12]
      --init-context-size <INIT_CONTEXT_SIZE> <INIT_CONTEXT_SIZE>
          Initial "fixed" seed window size in base pairs around the modified base. Example:
          --init-context-size 2 2 [default: 2 2]
      --min-sites <MIN_SITES>
          Minimum number of total sites in the genome required for a motif to be considered. [default:
          300]
      --min-frac-mod <FRAC_SITES_THRESH>
          Minimum fraction of sites in the genome to be "high-modification" for a motif to be
          considered. [default: 0.85]
      --known-motif <KNOWN_MOTIFS> <KNOWN_MOTIFS> <KNOWN_MOTIFS>
          Gather enrichment information for a known motif as well as compare to discovered motifs.
          Format should be <sequence> <offset> <mod_code>.
      --mod-code <MOD_CODES>
          Specify which modification codes to process, default will process all modification codes
          found in the input bedMethyl file.
      --force-override-spec
          Force override SAM specification of association of modification codes to primary sequence
          bases.
      --log-filepath <LOG_FILEPATH>
          Output log to this file.
      --suppress-progress
          Disable the progress bars.
  -h, --help
          Print help information.
```

## entropy
```text
Use a mod-BAM to calculate methylation entropy over genomic windows.

Usage: modkit entropy [OPTIONS] --in-bam <IN_BAMS> --ref <REFERENCE_FASTA>

Options:
  -s, --in-bam <IN_BAMS>
          Input mod-BAM, may be repeated multiple times to calculate entropy across all input
          mod-BAMs.

  -o, --out-bed <OUT_BED>
          Output BED file, if using `--region` this must be a directory.

      --prefix <PREFIX>
          Only used with `--regions`, prefix files in output directory with this string.

  -n, --num-positions <NUM_POSITIONS>
          Number of modified positions to consider at a time.
          
          [default: 4]

  -w, --window-size <WINDOW_SIZE>
          Maximum length interval that "num_positions" modified bases can occur in. The maximum
          window size decides how dense the positions are packed. For example, consider that the
          num_positions is equal to 4, the motif is CpG, and the window_size is equal to 8, this
          configuration would require that the modified positions are immediately adjacent to each
          other, "CGCGCGCG". On the other hand, if the window_size was set to 12, then multiple
          sequences with various patterns of other bases can be used CGACGATCGGCG.
          
          [default: 50]

      --no-filtering
          Do not perform any filtering, include all mod base calls in output.

      --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. Reads will be sampled
          evenly across aligned genome. If a region is specified, either with the --region option or
          the --sample-region option, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand reads is sufficient
          to estimate the model output distribution and determine the filtering threshold.
          
          [default: 10042]

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted variant is below
          this confidence percentile. For example, 0.1 will filter out the 10% lowest confidence
          modification calls.
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or for the canonical calls. When specified, base
          modification call probabilities will be required to be greater than or equal to this
          number. If `--mod-thresholds` is also specified, _this_ value will be used for canonical
          calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless the
          `--filter-threshold` option is also passed. See the online documentation for more details.

  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of BAM-reading threads to use.

      --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format.

      --mask
          Respect soft masking in the reference FASTA.

      --motif <MOTIF> <MOTIF>
          Motif to use for entropy calculation, default will be CpG.

      --cpg
          Use CpG motifs. Short hand for --motif CG 0 --combine-strands

      --base <BASE>
          Primary sequence base to calculate modification entropy on.
          
          [possible values: A, C, G, T]

      --regions <REGIONS_FP>
          Regions over which to calculate descriptive statistics.

      --combine-strands
          Combine modification counts on the positive and negative strands and report entropy on
          just the positive strand.

      --min-coverage <MIN_VALID_COVERAGE>
          Minimum coverage required at each position in the window. Windows without at least this
          many valid reads will be skipped, but positions within the window with enough coverage can
          be used by neighboring windows.
          
          [default: 3]

      --log-filepath <LOG_FILEPATH>
          Send debug logs to this file, setting this file is recommended.

      --suppress-progress
          Hide progress bars.

      --force
          Force overwrite output.

      --header
          Write a header line.

      --drop-zeros
          Omit windows with zero entropy.

      --max-filtered-positions <MAX_FILTERED_POSITIONS>
          Maximum number of filtered positions a read is allowed to have in a window, more than this
          number and the read will be discarded. Default will be 50% of `num_positions`.

  -h, --help
          Print help information (use `-h` for a summary).
```

## dmr pair
```text
Compare regions in a pair of samples (for example, tumor and normal or control and experiment). A
sample is input as a bgzip pileup bedMethyl (produced by pileup, for example) that has an associated
tabix index. Output is a BED file with the score column indicating the magnitude of the difference
in methylation between the two samples. See the online documentation for additional details.

Usage: modkit dmr pair [OPTIONS] -a <CONTROL_BED_METHYL> -b <EXP_BED_METHYL> --ref <REFERENCE_FASTA>

Options:
  -a <CONTROL_BED_METHYL>
          Bgzipped bedMethyl file for the first (usually control) sample. There should be a tabix
          index with the same name and .tbi next to this file or the --index-a option must be
          provided.

  -b <EXP_BED_METHYL>
          Bgzipped bedMethyl file for the second (usually experimental) sample. There should be a
          tabix index with the same name and .tbi next to this file or the --index-b option must be
          provided.

  -o, --out-path <OUT_PATH>
          Path to file to direct output, optional, no argument will direct output to stdout

      --header
          Include header in output.

  -r, --regions-bed <REGIONS_BED>
          BED file of regions over which to compare methylation levels. Should be tab-separated
          (spaces allowed in the "name" column). Requires chrom, chromStart and chromEnd. The Name
          column is optional. Strand is currently ignored. When omitted, methylation levels are
          compared at each site.

      --ref <REFERENCE_FASTA>
          Path to reference fasta for used in the pileup/alignment

      --segment <SEGMENTATION_FP>
          Run segmentation, output segmented differentially methylated regions to this file.

      --max-gap-size <MAX_GAP_SIZE>
          Maximum number of base pairs between modified bases for them to be segmented together.
          
          [default: 5000]

      --dmr-prior <DMR_PRIOR>
          Prior probability of a differentially methylated position.
          
          [default: 0.1]

      --diff-stay <DIFF_STAY>
          Maximum probability of continuing a differentially methylated block, decay will be dynamic
          based on proximity to the next position.
          
          [default: 0.9]

      --significance-factor <SIGNIFICANCE_FACTOR>
          Significance factor, effective p-value necessary to favor the "Different" state.
          
          [default: 0.01]

      --log-transition-decay
          Use logarithmic decay for "Different" stay probability.

      --decay-distance <DECAY_DISTANCE>
          After this many base pairs, the transition probability will become the prior probability
          of encountering a differentially modified position.
          
          [default: 500]

      --fine-grained
          Preset HMM segmentation parameters for higher propensity to switch from "Same" to
          "Different" state. Results will be shorter segments, but potentially higher sensitivity.

  -m, --base <MODIFIED_BASES>
          Bases to use to calculate DMR, may be multiple. For example, to calculate differentially
          methylated regions using only cytosine modifications use --base C.

      --assign-code <MOD_CODE_ASSIGNMENTS>
          Extra assignments of modification codes to their respective primary bases. In general,
          modkit dmr will use the SAM specification to know which modification codes are appropriate
          to use for a given primary base. For example "h" is the code for 5hmC, so is appropriate
          for cytosine bases, but not adenine bases. However, if your bedMethyl file contains custom
          codes or codes that are not part of the specification, you can specify which primary base
          they belong to here with --assign-code x:C meaning associate modification code "x" with
          cytosine (C) primary sequence bases. If a code is encountered that is not part of the
          specification, the bedMethyl record will not be used, this will be logged.

      --log-filepath <LOG_FILEPATH>
          File to write logs to, it's recommended to use this option.

  -t, --threads <THREADS>
          Number of threads to use.
          
          [default: 4]

      --batch-size <BATCH_SIZE>
          Control the  batch size. The batch size is the number of regions to load at a time. Each
          region will be processed concurrently. Loading more regions at a time will decrease IO to
          load data, but will use more memory. Default will be 50% more than the number of threads
          assigned.

  -k, --mask
          Respect soft masking in the reference FASTA.

      --suppress-progress
          Don't show progress bars.

  -f, --force
          Force overwrite of output file, if it already exists.

      --index-a <INDEX_A>
          Path to tabix index associated with -a (--control-bed-methyl) bedMethyl file.

      --index-b <INDEX_B>
          Path to tabix index associated with -b (--exp-bed-methyl) bedMethyl file.

      --missing <HANDLE_MISSING>
          How to handle regions found in the `--regions` BED file. quiet => ignore regions that are
          not found in the tabix header warn => log (debug) regions that are missing fatal => log
          (error) and exit the program when a region is missing.
          
          [default: warn]
          [possible values: quiet, warn, fail]

      --min-valid-coverage <MIN_VALID_COVERAGE>
          Minimum valid coverage required to use an entry from a bedMethyl. See the help for pileup
          for the specification and description of valid coverage.
          
          [default: 0]

      --prior <PRIOR> <PRIOR>
          Prior distribution for estimating MAP-based p-value. Should be two arguments for alpha and
          beta (e.g. 1.0 1.0). See `dmr_scoring_details.md` for additional details on how the metric
          is calculated.

      --delta <DELTA>
          Consider only effect sizes greater than this when calculating the MAP-based p-value.
          
          [default: 0.05]

  -N, --n-sample-records <N_SAMPLE_RECORDS>
          Sample this many reads when estimating the max coverage thresholds.
          
          [default: 10042]

      --max-coverages <MAX_COVERAGES> <MAX_COVERAGES>
          Max coverages to enforce when calculating estimated MAP-based p-value.

      --cap-coverages
          When using replicates, cap coverage to be equal to the maximum coverage for a single
          sample. For example, if there are 3 replicates with max_coverage of 30, the total coverage
          would normally be 90. Using --cap-coverages will down sample the data to 30X.

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller interval chunk sizes
          will use less memory but incur more overhead.
          
          [default: 100000]

  -h, --help
          Print help information (use `-h` for a summary).
```

## dmr multi
```text
Compare regions between all pairs of samples (for example a trio sample set or haplotyped trio
sample set). As with `pair` all inputs must be bgzip compressed bedMethyl files with associated
tabix indices. Each sample must be assigned a name. Output is a directory of BED files with the
score column indicating the magnitude of the difference in methylation between the two samples
indicated in the file name. See the online documentation for additional details

Usage: modkit dmr multi [OPTIONS] --regions-bed <REGIONS_BED> --out-dir <OUT_DIR> --ref <REFERENCE_FASTA>

Options:
  -s, --sample <SAMPLES> <SAMPLES>
          Two or more named samples to compare. Two arguments are required <path> <name>. This
          option should be repeated at least two times.
  -i, --index <INDICES> <INDICES>
          Optional, paths to tabix indices associated with named samples. Two arguments are required
          <path> <name> where <name> corresponds to the name of the sample given to the -s/--sample
          argument.
  -r, --regions-bed <REGIONS_BED>
          BED file of regions over which to compare methylation levels. Should be tab-separated
          (spaces allowed in the "name" column). Requires chrom, chromStart and chromEnd. The Name
          column is optional. Strand is currently ignored.
      --header
          Include header in output.
  -o, --out-dir <OUT_DIR>
          Directory to place output DMR results in BED format.
  -p, --prefix <PREFIX>
          Prefix files in directory with this label.
      --ref <REFERENCE_FASTA>
          Path to reference fasta for the pileup.
  -m, --base <MODIFIED_BASES>
          Bases to use to calculate DMR, may be multiple. For example, to calculate differentially
          methylated regions using only cytosine modifications use --base C.
      --assign-code <MOD_CODE_ASSIGNMENTS>
          Extra assignments of modification codes to their respective primary bases. In general,
          modkit dmr will use the SAM specification to know which modification codes are appropriate
          to use for a given primary base. For example "h" is the code for 5hmC, so is appropriate
          for cytosine bases, but not adenine bases. However, if your bedMethyl file contains custom
          codes or codes that are not part of the specification, you can specify which primary base
          they belong to here with --assign-code x:C meaning associate modification code "x" with
          cytosine (C) primary sequence bases. If a code is encountered that is not part of the
          specification, the bedMethyl record will not be used, this will be logged.
      --log-filepath <LOG_FILEPATH>
          File to write logs to, it's recommended to use this option.
  -t, --threads <THREADS>
          Number of threads to use. [default: 4]
  -k, --mask
          Respect soft masking in the reference FASTA.
      --suppress-progress
          Don't show progress bars.
  -f, --force
          Force overwrite of output file, if it already exists.
      --missing <HANDLE_MISSING>
          How to handle regions found in the `--regions` BED file. quiet => ignore regions that are
          not found in the tabix header warn => log (debug) regions that are missing fatal => log
          (error) and exit the program when a region is missing. [default: warn] [possible values:
          quiet, warn, fail]
      --min-valid-coverage <MIN_VALID_COVERAGE>
          Minimum valid coverage required to use an entry from a bedMethyl. See the help for pileup
          for the specification and description of valid coverage. [default: 0]
  -h, --help
          Print help information.
```
