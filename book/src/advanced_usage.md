# modkit, subcommand documentation

The goal of `modkit` is to enable best-practices manipulation of BAM files containing
modified base information (modBAMs). The various sub-commands and tools available in
`modkit` are described below.  This information can be obtained by invoking the long help
(`--help`) for each command.

> Advanced usage information.


```text
Modkit is a bioinformatics tool for working with modified bases from Oxford
Nanopore

Usage: modkit <COMMAND>

Commands:
  pileup          Tabulates base modification calls across genomic positions.
                  This command produces a bedMethyl formatted file. Schema and
                  description of fields can be found in the README
  extract         Extract read-level base modification information from a modBAM
                  into a tab-separated values table
  dmr             Perform DMR test on a set of regions. Output a BED file of
                  regions with the score column indicating the magnitude of the
                  difference. Find the schema and description of fields can in
                  the README as well as a description of the model and method.
                  See subcommand help for additional details
  pileup-hemi     Tabulates double-stranded base modification patters (such as
                  hemi-methylation) across genomic motif positions. This command
                  produces a bedMethyl file, the schema can be found in the
                  online documentation
  validate        Validate results from a set of mod-BAM files and associated
                  BED files containing the ground truth modified base status at
                  reference positions
  motif           Various commands to search for, evaluate, or further regine
                  sequence motifs enriched for base modification. Also can
                  generate BED files of motif locations
  entropy         Use a mod-BAM to calculate methylation entropy over genomic
                  windows
  localize        Investigate patterns of base modifications, by aggregating
                  pileup counts "localized" around genomic features of interest
  stats           Calculate base modification levels over regions
  bedmethyl       Utilities to work with bedMethyl files
  modbam          Utilities to work with modBAM files
  open-chromatin  Identify regions of open chromatin based on exogenous 6mA
                  signal
  help            Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

## pileup
```text
Tabulates base modification calls across genomic positions. This command
produces a bedMethyl formatted file. Schema and description of fields can be
found in the README

Usage: modkit pileup [OPTIONS] <IN_BAM> <OUT_BED>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index available

  <OUT_BED>
          Output file to write results into. Specify "-" or "stdout" to direct
          output to stdout

Options:
      --duplex
          Specify that modBAM contains duplex modification calls and these
          should be used

      --modified-bases <MODIFIED_BASES>...
          Specify which modified bases to tablulate counts for. These can be the
          "long name" such as '5mC', '6mA' (or 'm6A' for RNA), or 'inosine'. You
          can also pass <primary_base>:<mod_code>, such as 'C:m'. Finally, when
          running with --combine-mods the reference base can be passed alone,
          such as 'C' or 'A'

  -h, --help
          Print help (see a summary with '-h')

Output Options:
      --bgzf
          Output bgzf-compressed bedmethyl files suitable for Tabix-indexing

      --header
          Output bedMethyl where the delimiter of columns past column 10 are
          space-delimited instead of tab-delimited. This option can be useful
          for some browsers and parsers that don't expect the extra columns of
          the bedMethyl format. Output a header with the bedMethyl

      --prefix <PREFIX>
          Prefix to prepend on phased output file names. Without this option the
          files will be <mod_code>_<strand>.bedmethyl

      --phased
          Produce separate bedMethyl tables for haplotype-labeled modBAM
          records. Currently only supports diploid genomes. Produces
          hp1.bedmethyl, hp2.bedmethyl, and combined.bedmethyl files (optionally
          compressed). hp1.bedmethyl and hp2.bedmethyl contain counts for
          records with HP=1 and HP=2 tags, respectively. combined.bedmethyl
          contains counts for all modBAM records

      --bedrmod
          Output BedRModV2 header counts based on V2 Specification. Details can
          be found at https://dieterich-lab.github.io/euf-specs/bedRModv2.pdf

Compute Options:
      --bgzf-threads <BGZF_THREADS>
          Specify the number of threads to dedicate to BGZF compression
          
          [default: 4]

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of IO/decompression threads to use
          
          [default: 4]

  -s, --sampling-threads <SAMPLING_THREADS>
          Use this many threads when performing threshold estimation, setting
          this to a higher value can speed up execution

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead
          
          [default: 1000000]

      --queue-size <QUEUE_SIZE>
          Size of queue for writing records, default will be the number of
          threads

      --max-depth <MAX_DEPTH>
          Maximum number of reads to use whilest tabulating counts at a given
          position
          
          [default: 8000]

      --high-depth
          Flag to indicate that modBAM has exceptionally high depth (>65,000X)
          and should be downsampled during pileup

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended. (alias: log)

      --suppress-progress
          Hide the progress bar

Selection Options:
      --region <REGION>
          Process only the specified region of the BAM when performing pileup.
          Format should be <chrom_name>:<start>-<end> or <chrom_name>. Commas
          are allowed

      --include-bed <INCLUDE_BED>
          BED file that will restrict threshold estimation and pileup results to
          positions overlapping intervals in the file. (alias:
          include-positions)

      --include-unmapped
          Include unmapped base modifications when estimating the pass threshold

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

Sampling Options:
  -n, --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. Reads
          will be sampled evenly across aligned genome. If a region is
          specified, either with the --region option or the --sample-region
          option, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand
          reads is sufficient to estimate the model output distribution and
          determine the filtering threshold
          
          [default: 10042]

      --strict-sampling
          Fail if we cannot sample at least the requested number of reads when
          estimating the pass threshold

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the pass-threshold.
          In practice, 10-100 thousand reads is sufficient to estimate the model
          output distribution and determine the filtering threshold. See
          filtering.md for details on filtering

      --seed <SEED>
          Set a random seed for deterministic running, the default is
          non-deterministic

Filtering Options:
      --no-filtering
          Do not perform any filtering, include all mod base calls in output.
          See filtering.md for details on filtering

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter
          threshold can be specified with by a decimal number (e.g. 0.75).
          Per-base thresholds can be specified by colon-separated values, for
          example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

      --max-filter-threshold <MAX_FILTER_THRESHOLD>
          Maximum threhsold value to use, if the estimated threshold value for
          any primary sequence base is greater than this value, use this value
          instead. To set a default value for all primary bases, use
          --max-threshold 0.8, for example. Or use per-base threshold values
          like --max-threshold C:0.8

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold
          probability. If this option is not provided, but --region is provided,
          the genomic interval passed to --region will be used. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently when
          estimating the threshold probability, can be larger than the pileup
          processing interval
          
          [default: 1000000]

Modified Base Options:
      --motif <MOTIF> <MOTIF>
          Output pileup counts for only sequence motifs provided. The first
          argument should be the sequence motif and the second argument is the
          0-based offset to the base to pileup base modification counts for. For
          example: --motif CGCG 0 indicates to pileup counts for the first C on
          the top strand and the last C (complement to G) on the bottom strand.
          The --cpg argument is short hand for --motif CG 0.
          
          This argument can be passed multiple times. When more than one motif
          is used, the resulting output BED file will indicate the motif in the
          "name" field as <mod_code>,<motif>,<offset>. For example, given
          `--motif CGCG 2 --motif CG 0` there will be output lines with name
          fields such as "m,CG,0" and "m,CGCG,2". To use this option with
          `--combine-strands`, all motifs must be reverse-complement palindromic
          or an error will be raised.

      --cpg
          Only output counts at CpG motifs

  -r, --reference <REFERENCE_FASTA>
          Reference sequence in FASTA format. Required for motif (e.g. CpG)
          filtering, requires FAI fasta index to be pre-generated. (alias:
          'ref')

      --preload-references
          Preload the reference sequences, useful when working with many, short
          reference sequences such as a transcriptome

  -k, --mask
          Respect soft masking in the reference FASTA

      --combine-mods
          Combine base modification calls, all counts of modified bases are
          summed together. See collapse.md for details

      --allow-non-primary
          Use non-primary alignments in pileup (e.g. secondary and
          supplementaty). Keep in mind these mappings will be used _in addition_
          to the primary mappings

      --combine-strands
          When performing motif analysis (such as CpG), sum the counts from the
          positive and negative strands into the counts for the positive strand
          position. Requires that the BAM records have correct MN tags

BedRMod Options:
      --organism <ORGANISM>
          NCBI Taxonomic identifier, details at:
          https://doi.org/10.1093/database/baaa062
          
          [default: unknown]

      --modification-type <MODIFICATION_TYPE>
          A valid RNA type
          
          [default: RNA]

      --assembly <ASSEMBLY>
          Genome assembly
          
          [default: unknown]

      --annotation-source <ANNOTATION_SOURCE>
          Annotation source
          
          [default: unknown]

      --annotation-version <ANNOTATION_VERSION>
          Annotation version
          
          [default: unknown]

      --sequencing-platform <SEQUENCING_PLATFORM>
          Sequencing platform
          
          [default: ont]

      --basecalling <BASECALLING>
          Basecalling model, override model set in BAM header

      --bioinformatics-workflow <BIOINFORMATICS_WORKFLOW>
          Link to bioinformatics workflow; program name, version, and/or call;
          information relevant to score, coverage, or frequency calculation; etc

      --experiment <EXPERIMENT>
          Information about or link to experimental protocol and design

      --external-source <EXTERNAL_SOURCE>
          Databank:ID of data

      --modomics-code <MODOMICS_CODE> <MODOMICS_CODE> <MODOMICS_CODE>
          User-specified code, format: <Code> <MODOMICS Short Name> <Primary
          Base>. <Code> is the base modification code described in the SAMtags
          or a numeric ChEBI code present in the BAM file. For example 'a m6A A'
          would be for <Code> 'a', <MODOMICS Short Name> m6A, and <Primary Base>
          A
```

## adjust-mods
```text
Performs various operations on BAM files containing base modification
information, such as converting base modification codes and ignoring
modification calls. Produces a BAM output file

Usage: modkit adjust-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
          Input BAM file, can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input

  <OUT_BAM>
          File path to new BAM file to be created. Can be a path to a file or
          one of `-` or `stdin` to specify a stream from standard output

Options:
  -f, --ff
          Fast fail, stop processing at the first invalid sequence record.
          Default behavior is to continue and report failed/skipped records at
          the end

  -h, --help
          Print help (see a summary with '-h')

Output Options:
      --log-filepath <LOG_FILEPATH>
          Output debug logs to file at this path

      --output-sam
          Output SAM format instead of BAM

Modified Base Options:
      --ignore <IGNORE>
          Modified base code to ignore/remove, see
          https://samtools.github.io/hts-specs/SAMtags.pdf for details on the
          modified base codes

      --convert <CONVERT> <CONVERT>
          Convert one mod-tag to another, summing the probabilities together if
          the retained mod tag is already present

      --motif <MOTIF> <MOTIF>
          Filter out any base modification call that isn't part of a basecall
          sequence motif. This argument can be passed multiple times. Format is
          <motif_sequence> <offset>. For example the argument to match CpG
          dinucleotides is `--motif CG 0`, or to match CG[5mC]G the argument
          would be `--motif CGCG 2`. Single bases can be used as motifs to keep
          only base modification calls for a specific primary base, for example
          `--motif C 0`

      --cpg
          Shorthand for --motif CG 0

      --discard-motifs
          Discard base modification calls that match the provided motifs
          (instead of keeping them)

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

Selection Options:
      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

      --filter-probs
          Filter out the lowest confidence base modification probabilities

      --mapped-only
          Only use base modification probabilities from bases that are aligned
          when estimating the filter threshold (i.e. ignore soft-clipped, and
          inserted bases)

Sampling Options:
  -n, --num-reads <NUM_READS>
          Sample approximately this many reads when estimating the filtering
          threshold. If alignments are present reads will be sampled evenly
          across aligned genome. If a region is specified, either with the
          --region option or the --sample-region option, then reads will be
          sampled evenly across the region given. This option is useful for
          large BAM files. In practice, 10-50 thousand reads is sufficient to
          estimate the model output distribution and determine the filtering
          threshold
          
          [default: 10042]

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold
          probability. If this option is not provided, but --region is provided,
          the genomic interval passed to --region will be used. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size to process concurrently when estimating the
          threshold probability, can be larger than the pileup processing
          interval
          
          [default: 1000000]

Filtering Options:
  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per primary base. A global
          filter threshold can be specified with by a decimal number (e.g.
          0.75). Per-base thresholds can be specified by colon-separated values,
          for example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

Logging Options:
      --suppress-progress
          Hide the progress bar
```

## update-tags
```text
Renames Mm/Ml to tags to MM/ML. Also allows changing the mode flag from silent
'.' to explicitly '?' or '.'

Usage: modkit update-tags [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM to update modified base tags in. Can be a path to a file or one
             of `-` or `stdin` to specify a stream from standard input
  <OUT_BAM>  File to new BAM file to be created or one of `-` or `stdin` to
             specify a stream from standard output

Options:
  -m, --mode <MODE>        Mode, change mode to this value, options {'explicit',
                           'implicit'}. See spec at:
                           https://samtools.github.io/hts-specs/SAMtags.pdf.
                           'explicit' ('?') means residues without modification
                           probabilities will not be assumed canonical or
                           modified. 'implicit' means residues without explicit
                           modification probabilities are assumed to be
                           canonical [possible values: explicit, implicit]
      --no-implicit-probs  Don't add implicit canonical calls. This flag is
                           important when converting from one of the implicit
                           modes ( `.` or `""`) to explicit mode (`?`). By
                           passing this flag, the bases without associated base
                           modification probabilities will not be assumed to be
                           canonical. No base modification probability will be
                           written for these bases, meaning there is no
                           information. The mode will automatically be set to
                           the explicit mode `?`
  -h, --help               Print help

Compute Options:
  -t, --threads <THREADS>  Number of threads to use [default: 4]

Logging Options:
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path

Output Options:
      --output-sam  Output SAM format instead of BAM
```

## sample-probs
```text
Calculate an estimate of the base modification probability distribution

Usage: modkit sample-probs [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input BAM with modified base tags. If a index is found reads will be
          sampled evenly across the length of the reference sequence. Can be a
          path to a file or one of `-` or `stdin` to specify a stream from
          standard input

Options:
      --ignore-index
          Ignore the BAM index (if it exists) perform a serial scan of the BAM,
          otherwise if an index is found, multiple workers will be used to read
          the BAM in parallel. Conflicts with --threads since there will only be
          one reader

  -r, --reference <REFERENCE_FASTA>
          Reference sequence in FASTA format. (alias: 'ref')

      --preload-references
          Preload the reference sequences, useful when working with many, short
          reference sequences such as a transcriptome

  -k, --mask
          Respect soft masking in the reference FASTA

  -h, --help
          Print help (see a summary with '-h')

Compute Options:
  -t, --threads <THREADS>
          Number of workers to run concurrently. When an aligned, indexed BAM is
          found, multiple workers (threads) will be spawned to read disjoint
          sections of the BAM in parallel
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of IO/decompression threads to use, these are shared across
          workers (if an indexed BAM is found) or used for a single reader if
          --ignore-index is passed, stdin is used, or an index is not found
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead.
          Only used when sampling probs from an indexed bam
          
          [default: 1000000]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended

      --suppress-progress
          Hide the progress bar

Output Options:
  -p, --percentiles <PERCENTILES>
          Percentiles to calculate, a space separated list of floats (for
          example: "10,20,90")
          
          [default: 10,50,90]

  -q, --quantiles <QUANTILES>
          Quantiles to calculate, a space separated list of floats (for example:
          "0.1,0.2,0.9")

  -o, --out-dir <OUT_DIR>
          Directory to deposit result tables into. Required for model
          probability histogram output

      --prefix <PREFIX>
          Label to prefix output files with

      --force
          Overwrite results if present

      --hist
          Output histogram of base modification prediction probabilities

      --dna-color <PRIMARY_BASE_COLORS> <PRIMARY_BASE_COLORS>
          Set colors of primary bases in histogram, should be RGB format, e.g.
          "#0000FF" is defailt for canonical cytosine

      --mod-color <MOD_BASE_COLORS> <MOD_BASE_COLORS>
          Set colors of modified bases in histogram, should be RGB format, e.g.
          "#FF00FF" is default for 5hmC

Selection Options:
      --motif <MOTIF> <MOTIF>
          Tabulate base modification probabilities at motifs (or single bases).
          The first argument should be the sequence motif and the second
          argument is the 0-based offset to the base to pileup base modification
          counts for. For example: --motif CGCG 0 indicates to collect base
          mofification probabilities for the first C on the top strand and the
          last C (complement to G) on the bottom strand. For single bases,
          --motif C 0 or --motif A 0 would restrict to probabilities on C or A
          reference bases, respectively. This argument can be passed multiple
          times

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --region <REGION>
          Process only the specified region of the BAM when collecting
          probabilities. Format should be <chrom_name>:<start>-<end> or
          <chrom_name>

      --include-bed <INCLUDE_BED>
          Only sample base modification probabilities that are aligned to the
          positions in this BED file. (alias: include-positions)

      --use-non-primary
          Use non-primary mappings, may cause double-counting of reads with
          primary and one or more non-primary mappings. Requires MN tag to be
          correct. See documentation for help on proper alignment

Sampling Options:
  -n, --num-reads <NUM_READS>
          

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of
          reads to sample, for example 0.1 will sample 1/10th of the reads

      --no-sampling
          No sampling, use all of the reads to calculate the filter thresholds

  -s, --seed <SEED>
          Random seed for deterministic running, the default is
          non-deterministic, only used when no BAM index is provided
```

## summary
```text
Summarize the mod tags present in a BAM and get basic statistics. The default
output is a totals table (designated by '#' lines) and a modification calls
table. Descriptions of the columns can be found in the README

Usage: modkit summary [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input modBam, can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input

Options:
      --ignore-index
          Ignore the BAM index (if it exists) perform a serial scan of the BAM,
          otherwise if an index is found, multiple workers will be used to read
          the BAM in parallel. Conflicts with --threads since there will only be
          one reader

  -r, --reference <REFERENCE_FASTA>
          Reference sequence in FASTA format. (alias: 'ref')

      --preload-references
          Preload the reference sequences, useful when working with many, short
          reference sequences such as a transcriptome

  -k, --mask
          Respect soft masking in the reference FASTA

      --cpg
          Shorthand for --motif CG 0. Requires a sorted, indexed modBAM

  -h, --help
          Print help (see a summary with '-h')

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of IO/decompression threads to use, these are shared across
          workers (if an indexed BAM is found) or used for a single reader if
          --ignore-index is passed, stdin is used, or an index is not found
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          When using regions, interval chunk size in base pairs to process
          concurrently. Smaller interval chunk sizes will use less memory but
          incur more overhead
          
          [default: 1000000]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended

      --suppress-progress
          Hide the progress bar

Output Options:
      --tsv
          Output summary as a tab-separated variables stdout instead of a table

      --combine-mods
          Combine modification counts per primary sequence base when calculating
          summary statistics

Sampling Options:
  -n, --num-reads <NUM_READS>
          Approximate maximum number of reads to use. If an indexed BAM is
          provided, the reads will be sampled evenly over the length of the
          aligned reference. If a region is passed with the --region option,
          they will be sampled over the region. Actual number of reads used may
          deviate slightly from this number

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of
          reads to sample when estimating the filter threshold. For example 0.1
          will sample 1/10th of the reads

  -s, --seed <SEED>
          Sets a random seed for deterministic running (when using
          --sample-frac)

Filtering Options:
  -p, --filter-quantile <FILTER_QUANTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence base modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter
          threshold can be specified with by a decimal number (e.g. 0.75).
          Per-base thresholds can be specified by colon-separated values, for
          example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

      --max-filter-threshold <MAX_FILTER_THRESHOLD>
          Maximum threhsold value to use, if the estimated threshold value for
          any primary sequence base is greater than this value, use this value
          instead. To set a default value for all primary bases, use
          --max-threshold 0.8, for example. Or use per-base threshold values
          like --max-threshold C:0.8

Selection Options:
      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --include-bed <INCLUDE_BED>
          Only summarize base modification probabilities that are aligned to the
          positions in this BED file. (alias: include-positions)

      --mapped-only
          Only use modBAM records that are mapped, don't include unmapped
          records

      --matched-only
          Only summarize base modifications that are mapped and are matches to
          the appropriate primary sequence base. For example, require that 5mC
          base modification calls be aligned to a reference C. Requires a
          reference FASTA and an sorted, indexed modBAM

      --motif <MOTIF> <MOTIF>
          Calculate base modification statistics at reference motifs. This
          option has the same behavior as "--matched-only" except it only
          includes specified motifs. The first argument should be the sequence
          motif and the second argument is the 0-based offset to the base to
          pileup base modification counts for. For example: --motif CGCG 0
          indicates to collect base mofification probabilities for the first C
          on the top strand and the last C (complement to G) on the bottom
          strand. For single bases, --motif C 0 or --motif A 0 would restrict to
          probabilities on C or A reference bases, respectively. This argument
          can be passed multiple times. Requires a sorted, indexed modBAM

      --region <REGION>
          Process only the specified region of the BAM when collecting
          probabilities. Format should be <chrom_name>:<start>-<end> or
          <chrom_name>
```

## call-mods
```text
Call mods from a modbam, creates a new modbam with probabilities set to 100% if
a base modification is called or 0% if called canonical

Usage: modkit call-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
          Input BAM, may be sorted and have associated index available. Can be a
          path to a file or one of `-` or `stdin` to specify a stream from
          standard input

  <OUT_BAM>
          Output BAM, can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input

Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended

      --ff
          Fast fail, stop processing at the first invalid sequence record.
          Default behavior is to continue and report failed/skipped records at
          the end

      --suppress-progress
          Hide the progress bar

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently
          
          [default: 4]

  -n, --num-reads <NUM_READS>
          Sample approximately this many reads when estimating the filtering
          threshold. If alignments are present reads will be sampled evenly
          across aligned genome. If a region is specified, either with the
          --region option or the --sample-region option, then reads will be
          sampled evenly across the region given. This option is useful for
          large BAM files. In practice, 10-50 thousand reads is sufficient to
          estimate the model output distribution and determine the filtering
          threshold
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the
          filter-percentile. In practice, 50-100 thousand reads is sufficient to
          estimate the model output distribution and determine the filtering
          threshold. See filtering.md for details on filtering

      --seed <SEED>
          Set a random seed for deterministic running, the default is
          non-deterministic, only used when no BAM index is provided

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold
          probability. If this option is not provided, but --region is provided,
          the genomic interval passed to --region will be used. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size to process concurrently when estimating the
          threshold probability, can be larger than the pileup processing
          interval
          
          [default: 1000000]

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per primary base. A global
          filter threshold can be specified with by a decimal number (e.g.
          0.75). Per-base thresholds can be specified by colon-separated values,
          for example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

      --no-filtering
          Don't filter base modification calls, assign each base modification to
          the highest probability prediction

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

      --motif <MOTIF> <MOTIF>
          Filter out any base modification call that isn't part of a basecall
          sequence motif This argument can be passed multiple times. Format is
          <motif_sequence> <offset>. For example the argument to match CpG
          dinucleotides is `--motif CG 0`, or to match CG[5mC]G the argument
          would be `--motif CGCG 2`

      --cpg
          Shorthand for --motif CG 0

      --discard-motifs
          Discard base modification calls that match the provided motifs
          (instead of keeping them)

      --output-sam
          Output SAM format instead of BAM

  -h, --help
          Print help (see a summary with '-h')
```

## repair
```text
Repair MM and ML tags in one bam with the correct tags from another. To use this
command, both modBAMs _must_ be sorted by read name. The "donor" modBAM's reads
must be a superset of the acceptor's reads. Extra reads in the donor are
allowed, and multiple reads with the same name (secondary, etc.) are allowed in
the acceptor. Reads with an empty SEQ field cannot be repaired and will be
rejected. Reads where there is an ambiguous alignment of the acceptor to the
donor will be rejected (and logged). See the full documentation for details

Usage: modkit repair [OPTIONS] --donor-bam <DONOR_BAM> --acceptor-bam <ACCEPTOR_BAM> --output-bam <OUTPUT_BAM>

Options:
  -d, --donor-bam <DONOR_BAM>
          Donor modBAM with original MM/ML tags. Must be sorted by read name
  -a, --acceptor-bam <ACCEPTOR_BAM>
          Acceptor modBAM with reads to have MM/ML base modification data
          projected on to. Must be sorted by read name
  -o, --output-bam <OUTPUT_BAM>
          output modBAM location
      --log-filepath <LOG_FILEPATH>
          File to write logs to, it is recommended to use this option as some
          reads may be rejected and logged here
  -t, --threads <THREADS>
          The number of threads to use [default: 4]
  -h, --help
          Print help
```

## validate
```text
Validate results from a set of mod-BAM files and associated BED files containing
the ground truth modified base status at reference positions

Usage: modkit validate [OPTIONS]

Options:
  -h, --help
          Print help (see a summary with '-h')

Sample Options:
      --bam-and-bed <BAM> <BED>
          Argument accepts 2 values. The first value is the BAM file path with
          modified base tags. The second is a bed file with ground truth
          reference positions. The name field in the ground truth bed file
          should be the short name (single letter code or ChEBI ID) for a
          modified base or `-` to specify a canonical base ground truth
          position. This argument can be provided more than once for multiple
          samples

  -c, --canonical-base <CANONICAL_BASE>
          Canonical base to evaluate. By default, this will be derived from mod
          codes in ground truth BED files. For ground truth with only canonical
          sites and/or ChEBI codes this values must be set
          
          [possible values: A, C, G, T]

Modified Base Options:
      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base
          modification probability equally across other options. For example, if
          collapsing 'h', with 'm' and canonical options, half of the
          probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md

Selection Options:
      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

      --min-identity <MIN_ALIGNMENT_IDENTITY>
          Only use reads with alignment identity >= this number, in Q-space
          (phred score)

      --min-length <MIN_ALIGNMENT_LENGTH>
          Remove reads with fewer aligned reference bases than this threshold

Filtering Options:
  -p, --filter-quantile <FILTER_QUANTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify modified base probability filter threshold value. If
          specified, --filter-threshold will override --filter-quantile

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

Logging Options:
      --suppress-progress
          Hide the progress bar

      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended. (alias: log)

Output Options:
  -o, --out-filepath <OUT_FILEPATH>
          Specify a file for machine parseable output
```

## pileup-hemi
```text
Tabulates double-stranded base modification patters (such as hemi-methylation)
across genomic motif positions. This command produces a bedMethyl file, the
schema can be found in the online documentation

Usage: modkit pileup-hemi [OPTIONS] --ref <REFERENCE_FASTA> <IN_BAM>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index available

Options:
  -o, --out-bed <OUT_BED>
          Output file to write results into. Will write to stdout if not
          provided

  -h, --help
          Print help (see a summary with '-h')

Modified Base Options:
      --cpg
          Aggregate double-stranded base modifications for CpG dinucleotides.
          This flag is short-hand for --motif CG 0

      --motif <MOTIF> <MOTIF>
          Specify the sequence motif to pileup double-stranded base modification
          pattern counts for. The first argument should be the sequence motif
          and the second argument is the 0-based offset to the base to pileup
          base modification counts for. For example: --motif CG 0 indicates to
          generate pattern counts for the C on the top strand and the following
          C (opposite to G) on the negative strand. The motif must be
          reverse-complement palindromic or an error will be raised. See the
          documentation for more examples and details

  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base
          modification probability equally across other options. For example, if
          collapsing 'h', with 'm' and canonical options, half of the
          probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md

      --force-allow-implicit
          Force allow implicit-canonical mode. By default modkit does not allow
          pileup with the implicit mode (e.g. C+m, no '.' or '?'). The
          `update-tags` subcommand is provided to update tags to the new mode.
          This option allows the interpretation of implicit mode tags: residues
          without modified base probability will be interpreted as being the
          non-modified base

  -k, --mask
          Respect soft masking in the reference FASTA

      --combine-mods
          Combine base modification calls, all counts of modified bases are
          summed together. See collapse.md for details

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended. (alias: log)

      --suppress-progress
          Hide the progress bar

Selection Options:
      --region <REGION>
          Process only the specified region of the BAM when performing pileup.
          Format should be <chrom_name>:<start>-<end> or <chrom_name>. Commas
          are allowed

      --max-depth <MAX_DEPTH>
          Maximum number of records to use when calculating pileup. This
          argument is passed to the pileup engine. If you have high depth data,
          consider increasing this value substantially. Must be less than
          2147483647 or an error will be raised
          
          [default: 8000]

      --include-bed <INCLUDE_BED>
          BED file that will restrict threshold estimation and pileup results to
          positions overlapping intervals in the file. (alias:
          include-positions)

      --include-unmapped
          Include unmapped base modifications when estimating the pass threshold

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead
          
          [default: 100000]

      --queue-size <QUEUE_SIZE>
          Size of queue for writing records
          
          [default: 1000]

      --chunk-size <CHUNK_SIZE>
          Break contigs into chunks containing this many intervals (see
          `interval_size`). This option can be used to help prevent excessive
          memory usage, usually with no performance penalty. By default, modkit
          will set this value to 1.5x the number of threads specified, so if 4
          threads are specified the chunk_size will be 6. A warning will be
          shown if this option is less than the number of threads specified

Sampling Options:
  -n, --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. Reads
          will be sampled evenly across aligned genome. If a region is
          specified, either with the --region option or the --sample-region
          option, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand
          reads is sufficient to estimate the model output distribution and
          determine the filtering threshold
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the
          filter-percentile. In practice, 50-100 thousand reads is sufficient to
          estimate the model output distribution and determine the filtering
          threshold. See filtering.md for details on filtering

      --seed <SEED>
          Set a random seed for deterministic running, the default is
          non-deterministic

Filtering Options:
      --no-filtering
          Do not perform any filtering, include all mod base calls in output.
          See filtering.md for details on filtering

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter
          threshold can be specified with by a decimal number (e.g. 0.75).
          Per-base thresholds can be specified by colon-separated values, for
          example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold
          probability. If this option is not provided, but --region is provided,
          the genomic interval passed to --region will be used. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently when
          estimating the threshold probability, can be larger than the pileup
          processing interval
          
          [default: 1000000]

Output Options:
      --only-tabs
          **Deprecated** The default output has all tab-delimiters. For
          bedMethyl output, separate columns with only tabs. The default is to
          use tabs for the first 10 fields and spaces thereafter. The default
          behavior is more likely to be compatible with genome viewers. Enabling
          this option may make it easier to parse the output with tabular data
          handlers that expect a single kind of separator

      --mixed-delim
          Output bedMethyl where the delimiter of columns past column 10 are
          space-delimited instead of tab-delimited. This option can be useful
          for some browsers and parsers that don't expect the extra columns of
          the bedMethyl format
```

## entropy
```text
Use a mod-BAM to calculate methylation entropy over genomic windows

Usage: modkit entropy [OPTIONS] --in-bam <IN_BAMS> --ref <REFERENCE_FASTA>

Options:
  -s, --in-bam <IN_BAMS>
          Input mod-BAM, may be repeated multiple times to calculate entropy
          across all input mod-BAMs

  -n, --num-positions <NUM_POSITIONS>
          Number of modified positions to consider at a time
          
          [default: 4]

  -w, --window-size <WINDOW_SIZE>
          Maximum length interval that "num_positions" modified bases can occur
          in. The maximum window size decides how dense the positions are
          packed. For example, consider that the num_positions is equal to 4,
          the motif is CpG, and the window_size is equal to 8, this
          configuration would require that the modified positions are
          immediately adjacent to each other, "CGCGCGCG". On the other hand, if
          the window_size was set to 12, then multiple sequences with various
          patterns of other bases can be used CGACGATCGGCG
          
          [default: 50]

      --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format

      --mask
          Respect soft masking in the reference FASTA

      --motif <MOTIF> <MOTIF>
          Motif to use for entropy calculation, multiple motifs can be used by
          repeating this option. When multiple motifs are used that specify
          different modified primary bases, all modification possibilities will
          be used in the calculation

      --cpg
          Use CpG motifs. Short hand for --motif CG 0 --combine-strands

      --base <BASE>
          Primary sequence base to calculate modification entropy on
          
          [possible values: A, C, G, T]

      --regions <REGIONS_FP>
          Regions over which to calculate descriptive statistics

      --combine-strands
          Combine modification counts on the positive and negative strands and
          report entropy on just the positive strand

      --min-coverage <MIN_VALID_COVERAGE>
          Minimum coverage required at each position in the window. Windows
          without at least this many valid reads will be skipped, but positions
          within the window with enough coverage can be used by neighboring
          windows
          
          [default: 3]

      --max-filtered-positions <MAX_FILTERED_POSITIONS>
          Maximum number of filtered positions a read is allowed to have in a
          window, more than this number and the read will be discarded. Default
          will be 50% of `num_positions`

  -h, --help
          Print help (see a summary with '-h')

Output Options:
  -o, --out-bed <OUT_BED>
          Output BED file, if using `--region` this must be a directory

      --prefix <PREFIX>
          Only used with `--regions`, prefix files in output directory with this
          string

      --force
          Force overwrite output

      --header
          Write a header line

      --drop-zeros
          Omit windows with zero entropy

Filtering Options:
      --no-filtering
          Do not perform any filtering, include all mod base calls in output

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or for the canonical calls. When
          specified, base modification call probabilities will be required to be
          greater than or equal to this number. If `--mod-threshold` is also
          specified, _this_ value will be used for canonical calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

Sampling Options:
      --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. Reads
          will be sampled evenly across aligned genome. If a region is
          specified, either with the --region option or the --sample-region
          option, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand
          reads is sufficient to estimate the model output distribution and
          determine the filtering threshold
          
          [default: 10042]

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of BAM-reading threads to use

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Send debug logs to this file, setting this file is recommended

      --verbose-logging
          Log regions that have zero or insufficient coverage. Requires log file

      --suppress-progress
          Hide progress bars
```

## localize
```text
Investigate patterns of base modifications, by aggregating pileup counts
"localized" around genomic features of interest

Usage: modkit localize [OPTIONS] --regions <REGIONS> --genome-sizes <GENOME_SIZES> <IN_BEDMETHYL>

Arguments:
  <IN_BEDMETHYL>
          Input bedMethyl table. Should be bgzip-compressed and have an
          associated Tabix index. The tabix index will be assumed to be
          $this_file.tbi

Options:
      --regions <REGIONS>
          BED file of regions to calculate enrichment around. These BED records
          serve as the points from which the `--window` number of bases is
          centered

  -w, --window <EXPAND_WINDOW>
          Number of base pairs to search around, for example if your BED region
          records are single positions, a window of 500 will look 500 base pairs
          upstream and downstream of that position. If your region BED records
          are larger regions, this will expand from the midpoint of that region
          
          [default: 2000]

  -s, --stranded <STRANDED>
          Whether to only keep bedMethyl records on the "same" strand or
          "opposite" strand
          
          [possible values: same, opposite]

      --stranded-features <STRANDED_FEATURES>
          Force use bedMethyl records from a particular strand, default is to
          use the strand as given in the BED file (will use BOTH for BED3)
          
          [possible values: positive, negative, both]

      --min-coverage <MIN_COVERAGE>
          Minimum valid coverage to use a bedMethyl record
          
          [default: 3]

  -r, --genome-sizes <GENOME_SIZES>
          TSV of genome sizes, should be <chrom>\t<size_in_bp>

  -o, --out-file <OUT_FILE>
          Optionally specify a file to write output to, default is stdout

  -h, --help
          Print help (see a summary with '-h')

Output Options:
      --chart <CHART_FILEPATH>
          Create plots showing %-modification vs. offset. Argument should be a
          path to a file

      --name <CHART_NAME>
          Give the HTML document and chart a name

  -f, --force
          Force overwrite of existing output file

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file to write debug logs to

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of tabix/bgzf IO threads to use
          
          [default: 2]

      --batch-size <BATCH_SIZE_BP>
          [default: 500000]
```

## stats
```text
Calculate base modification levels over regions

Usage: modkit stats [OPTIONS] --regions <REGIONS> --out-table <OUT_TABLE> <IN_BEDMETHYL>

Arguments:
  <IN_BEDMETHYL>  Input bedMethyl table. Should be bgzip-compressed and have an
                  associated Tabix index. The tabix index will be assumed to be
                  $this_file.tbi

Options:
      --regions <REGIONS>
          BED file of regions to aggregate base modification over
  -c, --mod-codes <MOD_CODES>
          Specify which base modification codes to use. Default will report
          information on all base modification codes encountered
  -m, --min-coverage <MIN_COVERAGE>
          Only use records with at least this much valid coverage [default: 1]
  -h, --help
          Print help

Output Options:
  -o, --out-table <OUT_TABLE>  Specify the output file to write the results
                               table, or "-"/"stdout" for stdout
      --force                  Force overwrite the output file
      --no-header              Don't add the header describing the columns to
                               the output

Logging Options:
      --log-filepath <LOG_FILEPATH>  Specify a file to write debug logs to

Compute Options:
  -t, --threads <THREADS>        Number of threads to use [default: 4]
      --io-threads <IO_THREADS>  Number of tabix/bgzf threads to use [default:
                                 2]
```

## open-chromatin predict
```text
Usage: modkit open-chromatin predict [OPTIONS] --model <MODEL_PATH> --output <OUTPUT> <IN_BAM>

Arguments:
  <IN_BAM>
          Input modBAM with 6mA base modification calls

Options:
      --model <MODEL_PATH>
          Path to directory with open-chromatin model

      --batch-size <BATCH_SIZE>
          Collect this many windows of data for each run through the model
          
          [default: 1024]

      --super-batch-size <SUPER_BATCH_SIZE>
          Number of "batches" to collect at once, see documentation for exact
          calculation or run with --dryrun to see output
          
          [default: 100]

      --step-size <STEP_SIZE>
          Number of base pairs to step along the genome or genomic intervals to
          make predictions. Smaller numbers will result in better resolution but
          take more computation
          
          [default: 25]

  -t, --min-coverage <MIN_COVERAGE>
          Require this many reads covering each prediction window
          
          [default: 5]

      --threshold <FILTER_THRESHOLD>
          Only emit records/windows where the open chromatin probability is
          greater than or equal to this value

  -o, --output <OUTPUT>
          Output bedGraph file, or "-"/"stdout" to write to stdout

      --force
          Force overwrite output file

      --log-filepath <LOG_FILEPATH>
          Path to file to write logs to, setting this parameter is recommended

      --include-bed <INCLUDE_BED>
          BED file of regions over which to make predictions

      --region <REGION>
          Region in the format <chrom>:start-end over which to make predictions

      --suppress-progress
          Hide the progress bar

      --no-header
          Don't print header in output

      --dryrun
          Perform setup and validations, but stop before performing inference

  -h, --help
          Print help (see a summary with '-h')
```

## extract full
```text
Transform the probabilities from the MM/ML tags in a modBAM into a table

Usage: modkit extract full [OPTIONS] <IN_BAM> <OUT_PATH>

Arguments:
  <IN_BAM>
          Path to modBAM file to extract read-level information from, or one of
          `-` or `stdin` to specify a stream from standard input. If a file is
          used it may be sorted and have associated index

  <OUT_PATH>
          Path to output file, "stdout" or "-" will direct output to standard
          out

Options:
      --reference <REFERENCE>
          Path to reference FASTA to extract reference context information from.
          Required for motif selection

  -h, --help
          Print help (see a summary with '-h')

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of IO/decompression threads to use
          
          [default: 4]

      --out-threads <OUT_THREADS>
          Number of threads to use for parallel bgzf writing
          
          [default: 4]

  -q, --queue-size <QUEUE_SIZE>
          Number of reads that can be in memory at a time. Increasing this value
          will increase thread usage, at the cost of memory usage
          
          [default: 10000]

      --ignore-index
          Ignore the BAM index (if it exists) and default to a serial scan of
          the BAM

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead.
          Only used when an indexed modBAM is provided
          
          [default: 100000]

Output Options:
      --bgzf
          Write output as BGZF compressed file

      --force
          Force overwrite of output file

      --kmer-size <KMER_SIZE>
          Set the query and reference k-mer size (if a reference is provided).
          Maximum number for this value is 50
          
          [default: 5]

      --no-headers
          Don't print the header lines in the output tables

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Path to file to write run log

      --suppress-progress
          Hide the progress bar

Selection Options:
      --mapped-only
          Include only mapped bases in output (alias: mapped)

      --allow-non-primary
          Output aligned secondary and supplementary base modification
          probabilities as additional rows. The primary alignment will have all
          of the base modification probabilities (including soft-clipped ones,
          unless --mapped-only is used). The non-primary alignments will only
          have mapped bases in the output

      --num-reads <NUM_READS>
          Number of reads to use. Note that when using a sorted, indexed modBAM
          that the sampling algorithm will attempt to sample records evenly over
          the length of the reference sequence. The result is the final number
          of records used may be slightly more or less than the requested
          number. When piping from stdin or using a modBAM without an index, the
          requested number of reads will be the first `num_reads` records

      --region <REGION>
          Process only reads that are aligned to a specified region of the BAM.
          Format should be <chrom_name>:<start>-<end> or <chrom_name>

      --include-bed <INCLUDE_BED>
          BED file with regions to include (alias: include-positions).
          Implicitly only includes mapped sites

  -v, --exclude-bed <EXCLUDE_BED>
          BED file with regions to _exclude_ (alias: exclude)

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

      --ignore-implicit
          Ignore implicitly canonical base modification calls. When the `.` flag
          is used in the MM tag, this implies that bases missing a base
          modification probability are to be assumed canonical. Set this flag to
          omit those base modifications from the output. For additional details
          see the SAM spec: https://samtools.github.io/hts-specs/SAMtags.pdf

Modified Base Options:
      --motif <MOTIF> <MOTIF>
          Output read-level base modification probabilities restricted to the
          reference sequence motifs provided. The first argument should be the
          sequence motif and the second argument is the 0-based offset to the
          base to pileup base modification counts for. For example: --motif CGCG
          0 indicates include base modifications for which the read is aligned
          to the first C on the top strand and the last C (complement to G) on
          the bottom strand. The --cpg argument is short hand for --motif CG 0.
          This argument can be passed multiple times

      --annotate-motifs
          When used with `--motif` or `--cpg` emit all modified base alignment
          information even if it does not align to a reference motif, but
          annotate which aligned positions match which motifs in the "motifs"
          column. "." will be used when an aligned position does not match a
          motif

      --cpg
          Only output counts at CpG motifs. Requires a reference sequence to be
          provided

  -k, --mask
          When using motifs, respect soft masking in the reference sequence

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base
          modification probability equally across other options. For example, if
          collapsing 'h', with 'm' and canonical options, half of the
          probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md
```

## extract calls
```text
Produce a table of read-level base modification calls. This table has, for each
read, one row for each base modification call in that read using the same
thresholding algorithm as in pileup, or summary (see online documentation for
details on thresholds)

Usage: modkit extract calls [OPTIONS] <IN_BAM> <OUT_PATH>

Arguments:
  <IN_BAM>
          Path to modBAM file to extract read-level information from, or one of
          `-` or `stdin` to specify a stream from standard input. If a file is
          used it may be sorted and have associated index

  <OUT_PATH>
          Path to output file, "stdout" or "-" will direct output to standard
          out

Options:
      --reference <REFERENCE>
          Path to reference FASTA to extract reference context information from.
          If no reference is provided, `ref_kmer` column will be "." in the
          output. (alias: ref)

      --preload-references
          Preload the reference sequences, useful when working with many, short
          reference sequences such as a transcriptome

  -h, --help
          Print help (see a summary with '-h')

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of IO/decompression threads to use
          
          [default: 4]

      --out-threads <OUT_THREADS>
          Number of threads to use for parallel bgzf writing
          
          [default: 4]

  -q, --queue-size <QUEUE_SIZE>
          Number of reads that can be in memory at a time. Increasing this value
          will increase thread usage, at the cost of memory usage
          
          [default: 10000]

      --ignore-index
          Ignore the BAM index (if it exists) and default to a serial scan of
          the BAM

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead.
          Only used when an indexed modBAM is provided
          
          [default: 100000]

Output Options:
      --bgzf
          Write output as BGZF compressed file

      --force
          Force overwrite of output file

      --kmer-size <KMER_SIZE>
          Set the query and reference k-mer size (if a reference is provided).
          Maximum number for this value is 50
          
          [default: 5]

      --no-headers
          Don't print the header lines in the output tables

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Path to file to write run log

      --suppress-progress
          Hide the progress bar

Selection Options:
      --mapped-only
          Include only mapped bases in output (alias: mapped)

      --allow-non-primary
          Output aligned secondary and supplementary base modification
          probabilities as additional rows. The primary alignment will have all
          of the base modification probabilities (including soft-clipped ones,
          unless --mapped-only is used). The non-primary alignments will only
          have mapped bases in the output

      --num-reads <NUM_READS>
          Number of reads to use. Note that when using a sorted, indexed modBAM
          that the sampling algorithm will attempt to sample records evenly over
          the length of the reference sequence. The result is the final number
          of records used may be slightly more or less than the requested
          number. When piping from stdin or using a modBAM without an index, the
          requested number of reads will be the first `num_reads` records

      --region <REGION>
          Process only reads that are aligned to a specified region of the BAM.
          Format should be <chrom_name>:<start>-<end> or <chrom_name>

      --include-bed <INCLUDE_BED>
          BED file with regions to include (alias: include-positions).
          Implicitly only includes mapped sites

  -v, --exclude-bed <EXCLUDE_BED>
          BED file with regions to _exclude_ (alias: exclude)

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

      --ignore-implicit
          Ignore implicitly canonical base modification calls. When the `.` flag
          is used in the MM tag, this implies that bases missing a base
          modification probability are to be assumed canonical. Set this flag to
          omit those base modifications from the output. For additional details
          see the SAM spec: https://samtools.github.io/hts-specs/SAMtags.pdf

      --pass-only
          Only output base modification calls that pass the minimum confidence
          threshold. (alias: pass)

Modified Base Options:
      --motif <MOTIF> <MOTIF>
          Output read-level base modification probabilities restricted to the
          reference sequence motifs provided. The first argument should be the
          sequence motif and the second argument is the 0-based offset to the
          base to pileup base modification counts for. For example: --motif CGCG
          0 indicates include base modifications for which the read is aligned
          to the first C on the top strand and the last C (complement to G) on
          the bottom strand. The --cpg argument is short hand for --motif CG 0.
          This argument can be passed multiple times

      --annotate-motifs
          When used with `--motif` or `--cpg` emit all modified base alignment
          information even if it does not align to a reference motif, but
          annotate which aligned positions match which motifs in the "motifs"
          column. "." will be used when an aligned position does not match a
          motif

      --cpg
          Only output counts at CpG motifs. Requires a reference sequence to be
          provided

  -k, --mask
          When using motifs, respect soft masking in the reference sequence

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base
          modification probability equally across other options. For example, if
          collapsing 'h', with 'm' and canonical options, half of the
          probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md

Filtering Options:
      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter
          threshold can be specified with by a decimal number (e.g. 0.75).
          Per-base thresholds can be specified by colon-separated values, for
          example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

      --max-filter-threshold <MAX_FILTER_THRESHOLD>
          Maximum threhsold value to use, if the estimated threshold value for
          any primary sequence base is greater than this value, use this value
          instead. To set a default value for all primary bases, use
          --max-threshold 0.8, for example. Or use per-base threshold values
          like --max-threshold C:0.8

      --no-filtering
          Don't estimate the pass threshold, all calls will "pass"

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

Sampling Options:
      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently when
          estimating the threshold probability
          
          [default: 1000000]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the pass-threshold.
          In practice, 10-100 thousand reads is sufficient to estimate the model
          output distribution and determine the filtering threshold. See
          filtering.md for details on filtering

  -n, --sample-num-reads <SAMPLE_NUM_READS>
          Sample this many reads when estimating the filtering threshold. If a
          sorted, indexed modBAM is provided reads will be sampled evenly across
          aligned genome. If a region is specified, with the --region, then
          reads will be sampled evenly across the region given. This option is
          useful for large BAM files. In practice, 10-50 thousand reads is
          sufficient to estimate the model output distribution and determine the
          filtering threshold
          
          [default: 10042]

      --seed <SEED>
          Set a random seed for deterministic running, the default is
          non-deterministic when using `sampling_frac`. When using `num_reads`
          the output is still deterministic
```

## extract read-stats
```text
Produce a table where modification counts are summarized on the read level. This
table will have one record per valid read and count the number of modified and
unmodified bases for each base modification requested

Usage: modkit extract read-stats [OPTIONS] <IN_BAM> <OUT_PATH>

Arguments:
  <IN_BAM>
          Input modBAM. Can be a path to a file or "-"/"stdin" for standard
          input

  <OUT_PATH>
          Path to file for output, default will be to stdout

Options:
      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter
          threshold can be specified with by a decimal number (e.g. 0.75).
          Per-base thresholds can be specified by colon-separated values, for
          example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed

      --mod-codes <MOD_CODES>...
          Specify which modcodes to tabulate, should be
          <primary_base>:<mod_code>

      --queue-size <QUEUE_SIZE>
          Queue size, maximum number of reads to be held in memory waiting for
          records to be written
          
          [default: 10000]

      --io-threads <IO_THREADS>
          [default: 4]

  -h, --help
          Print help (see a summary with '-h')

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Path to file to write run log

      --suppress-progress
          Hide the progress bar
```

## motif bed
```text
Create BED file with all locations of a sequence motif. Example: modkit motif
bed CG 0

Usage: modkit motif bed [OPTIONS] <FASTA> <MOTIF> <OFFSET>

Arguments:
  <FASTA>   Input FASTA file
  <MOTIF>   Motif to search for within FASTA, e.g. CG
  <OFFSET>  Offset within motif, e.g. 0

Options:
  -k, --mask  Respect soft masking in the reference FASTA
  -h, --help  Print help
```

## motif search
```text
Search for modification-enriched subsequences in a reference genome

Usage: modkit motif search [OPTIONS] --in-bedmethyl <IN_BEDMETHYL> --ref <REFERENCE_FASTA>

Options:
      --force-override-spec
          Force override SAM specification of association of modification codes
          to primary sequence bases

  -h, --help
          Print help (see a summary with '-h')

Input Options:
  -i, --in-bedmethyl <IN_BEDMETHYL>
          Input bedmethyl table, can be used directly from modkit pileup

  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format used for the pileup

      --contig <CONTIG>
          Use only bedMethyl records from this contig, requires that the
          bedMethyl be BGZIP-compressed and tabix-indexed

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of tabix/bgzf IO threads to use
          
          [default: 2]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Output log to this file

      --suppress-progress
          Disable the progress bars

Search Options:
      --low-thresh <LOW_THRESHOLD>
          Fraction modified threshold below which consider a genome location to
          be "low modification"
          
          [default: 0.2]

      --high-thresh <HIGH_THRESHOLD>
          Fraction modified threshold above which consider a genome location to
          be "high modification" or enriched for modification
          
          [default: 0.6]

      --min-frac-mod <FRAC_SITES_THRESH>
          Minimum fraction of sites in the genome to be "high-modification" for
          a motif to be considered
          
          [default: 0.85]

      --context-size <CONTEXT_SIZE> <CONTEXT_SIZE>
          Upstream and downstream number of bases to search for a motif sequence
          around a modified base. Example: --context-size 12 12
          
          [default: 12 12]

      --min-coverage <MIN_COVERAGE>
          Minimum valid coverage in the bedMethyl to consider a record valid
          
          [default: 5]

      --min-sites <MIN_SITES>
          Minimum number of total sites in the genome required for a motif to be
          considered
          
          [default: 300]

      --min-log-odds <MIN_LOG_ODDS>
          Minimum log-odds to consider a motif sequence to be enriched
          
          [default: 1.5]

      --init-context-size <INIT_CONTEXT_SIZE> <INIT_CONTEXT_SIZE>
          Initial "fixed" seed window size in base pairs around the modified
          base. Example: --init-context-size 2 2
          
          [default: 2 2]

      --mod-code <MOD_CODES>
          Specify which modification codes to process, default will process all
          modification codes found in the input bedMethyl file

Output Options:
  -o, --out-table <OUT_TABLE>
          Optionally output a machine-parsable TSV (human-readable table will
          always be output to the log)

      --known-motif <KNOWN_MOTIFS> <KNOWN_MOTIFS> <KNOWN_MOTIFS>
          Include statistics on a suspected or known motif. Format should be
          <sequence> <offset> <mod_code>

      --known-motifs-table <KNOWN_MOTIFS_TABLE>
          Path to known motifs in tabular format. Tab-separated values:
          <mod_code>\t<motif_seq>\t<offset>. May have the same header as the
          output table from this command

      --eval-motifs-table <OUT_KNOWN_TABLE>
          Optionally output machine parsable table with known motif modification
          frequencies that were not found during search

Exhaustive Search Options:
      --exhaustive-seed-min-log-odds <EXHAUSTIVE_SEED_MIN_LOG_ODDS>
          Minimum log-odds to consider a motif seed sequence to be enriched when
          performing exhaustive search, decreasing this number will increase the
          number of seeds searched and thus computational time
          
          [default: 2.5]

      --exhaustive-seed-len <EXHAUSTIVE_SEED_LEN>
          Exhaustive search seed length, increasing this value increases
          computational time
          
          [default: 3]

      --skip-search
          Skip the exhaustive search phase, saves time but the results may be
          less sensitive

      --search-top-pct <SEARCH_TOP_PCT>
          During exhaustive search, instead of searching all seeds with log-odds
          above `exhaustive_seed_min_log_odds`, only search the top X-percent of
          seeds. Can be used with `min_exhaustive_seeds` and
          `max_exhaustive_seeds`

      --narrow-search
          When used in conjunction with `search_top_pct`, search the top
          X-percent of seeds, and then narrow the search space by removing
          contexts matching any motifs found. Then iterate until zero additional
          motifs are found or another stopping condition is reached

      --search-timeout <SEARCH_TIMEOUT>
          A stopping condition when using `--narrow-search`, stop once exaustive
          search for a modification code has been worked on for this long

      --search-batch-size <SEARCH_BATCH_SIZE>
          Set the batch size when performing a simple timeout on search. At
          least this many seeds will be evaluated
          
          [default: 100]

      --max-exhaustive-seeds <MAX_EXHAUSTIVE_SEEDS>
          Set the maximum number of exhaustive seeds to be searched in a batch.
          Overrides the X-percent of seeds to be searched when that number
          exceeds this setting

      --min-exhaustive-seeds <MIN_EXHAUSTIVE_SEEDS>
          Search at least this many seeds. Overrides the X-percent of seeds to
          be searched when that number is less than this setting
          
          [default: 20]

      --max-narrow-iters <MAX_NARROW_ITERS>
          Stopping condition when using `--narrow-search` and
          `--search-top-pct`, stop after this many iterations regardless if the
          timeout is provided and has been reached. Exaustive search will still
          stop when once no more motifs are found
```

## motif evaluate
```text
Calculate enrichment statistics on a set of motifs from a bedMethyl table

Usage: modkit motif evaluate [OPTIONS] --in-bedmethyl <IN_BEDMETHYL> --ref <REFERENCE_FASTA>

Options:
      --force-override-spec
          Force override SAM specification of association of modification codes
          to primary sequence bases

      --suppress-table
          Don't print final table to stderr (will still go to log file)

  -h, --help
          Print help (see a summary with '-h')

Input Options:
  -i, --in-bedmethyl <IN_BEDMETHYL>
          Input bedmethyl table, can be used directly from modkit pileup

  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format used for the pileup

      --contig <CONTIG>
          Use only bedMethyl records from this contig, requires that the
          bedMethyl be BGZIP-compressed and tabix-indexed

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of tabix/bgzf IO threads to use
          
          [default: 2]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Output log to this file

      --suppress-progress
          Disable the progress bars

Output Options:
      --known-motif <KNOWN_MOTIFS> <KNOWN_MOTIFS> <KNOWN_MOTIFS>
          Format should be <sequence> <offset> <mod_code>

      --known-motifs-table <KNOWN_MOTIFS_TABLE>
          Path to known motifs in tabular format. Tab-separated values:
          <mod_code>\t<motif_seq>\t<offset>. May have the same header as the
          output table from this command

      --out <OUT_TABLE>
          Machine-parsable table of refined motifs. Human-readable table always
          printed to stderr and log

Search Options:
      --min-coverage <MIN_COVERAGE>
          Minimum valid coverage in the bedMethyl to consider a record valid
          
          [default: 5]

      --context-size <CONTEXT_SIZE> <CONTEXT_SIZE>
          Upstream and downstream number of bases to search for a motif sequence
          around a modified base. Example: --context-size 12 12
          
          [default: 12 12]

      --low-thresh <LOW_THRESHOLD>
          Fraction modified threshold below which consider a genome location to
          be "low modification"
          
          [default: 0.2]

      --high-thresh <HIGH_THRESHOLD>
          Fraction modified threshold above which consider a genome location to
          be "high modification" or enriched for modification
          
          [default: 0.6]
```

## motif refine
```text
Use a previously defined list of motif sequences and further refine them with a
bedMethyl table

Usage: modkit motif refine [OPTIONS] --in-bedmethyl <IN_BEDMETHYL> --ref <REFERENCE_FASTA>

Options:
      --min-refine-frac-mod <MIN_REFINE_FRAC_MODIFIED>
          Minimum fraction of sites in the genome to be "high-modification" for
          a motif to be further refined, otherwise it will be discarded
          
          [default: 0.6]

      --min-refine-sites <MIN_REFINE_SITES>
          Minimum number of total sites in the genome required for a motif to be
          further refined, otherwise it will be discarded
          
          [default: 300]

      --force-override-spec
          Force override SAM specification of association of modification codes
          to primary sequence bases

  -h, --help
          Print help (see a summary with '-h')

Input Options:
  -i, --in-bedmethyl <IN_BEDMETHYL>
          Input bedmethyl table, can be used directly from modkit pileup

  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format used for the pileup

      --contig <CONTIG>
          Use only bedMethyl records from this contig, requires that the
          bedMethyl be BGZIP-compressed and tabix-indexed

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of tabix/bgzf IO threads to use
          
          [default: 2]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Output log to this file

      --suppress-progress
          Disable the progress bars

Output Options:
      --known-motif <KNOWN_MOTIFS> <KNOWN_MOTIFS> <KNOWN_MOTIFS>
          Format should be <sequence> <offset> <mod_code>

      --known-motifs-table <KNOWN_MOTIFS_TABLE>
          Path to known motifs in tabular format. Tab-separated values:
          <mod_code>\t<motif_seq>\t<offset>. May have the same header as the
          output table from this command

      --out <OUT_TABLE>
          Machine-parsable table of refined motifs. Human-readable table always
          printed to stderr and log

Search Options:
      --low-thresh <LOW_THRESHOLD>
          Fraction modified threshold below which consider a genome location to
          be "low modification"
          
          [default: 0.2]

      --high-thresh <HIGH_THRESHOLD>
          Fraction modified threshold above which consider a genome location to
          be "high modification" or enriched for modification
          
          [default: 0.6]

      --min-frac-mod <FRAC_SITES_THRESH>
          Minimum fraction of sites in the genome to be "high-modification" for
          a motif to be considered
          
          [default: 0.85]

      --context-size <CONTEXT_SIZE> <CONTEXT_SIZE>
          Upstream and downstream number of bases to search for a motif sequence
          around a modified base. Example: --context-size 12 12
          
          [default: 12 12]

      --min-coverage <MIN_COVERAGE>
          Minimum valid coverage in the bedMethyl to consider a record valid
          
          [default: 5]

      --min-sites <MIN_SITES>
          Minimum number of total sites in the genome required for a motif to be
          considered
          
          [default: 300]

      --min-log-odds <MIN_LOG_ODDS>
          Minimum log-odds to consider a motif sequence to be enriched
          
          [default: 1.5]
```

## dmr pair
```text
Compare regions in a pair of samples (for example, tumor and normal or control
and experiment). A sample is input as a bgzip pileup bedMethyl (produced by
pileup, for example) that has an associated tabix index. Output is a BED file
with the score column indicating the magnitude of the difference in methylation
between the two samples. See the online documentation for additional details

Usage: modkit dmr pair [OPTIONS] --ref <REFERENCE_FASTA>

Options:
  -r, --regions-bed <REGIONS_BED>
          BED file of regions over which to compare methylation levels. Should
          be tab-separated (spaces allowed in the "name" column). Requires
          chrom, chromStart and chromEnd. The Name column is optional. Strand is
          currently ignored. When omitted, methylation levels are compared at
          each site

      --ref <REFERENCE_FASTA>
          Path to reference fasta for used in the pileup/alignment

  -h, --help
          Print help (see a summary with '-h')

Sample Options:
  -a <CONTROL_BED_METHYL>
          Bgzipped bedMethyl file for the first (usually control) sample. There
          should be a tabix index with the same name and .tbi next to this file
          or the --index-a option must be provided

  -b <EXP_BED_METHYL>
          Bgzipped bedMethyl file for the second (usually experimental) sample.
          There should be a tabix index with the same name and .tbi next to this
          file or the --index-b option must be provided

  -m, --base <MODIFIED_BASES>
          Bases to use to calculate DMR, may be multiple. For example, to
          calculate differentially methylated regions using only cytosine
          modifications use --base C

      --assign-code <MOD_CODE_ASSIGNMENTS>
          Extra assignments of modification codes to their respective primary
          bases. In general, modkit dmr will use the SAM specification to know
          which modification codes are appropriate to use for a given primary
          base. For example "h" is the code for 5hmC, so is appropriate for
          cytosine bases, but not adenine bases. However, if your bedMethyl file
          contains custom codes or codes that are not part of the specification,
          you can specify which primary base they belong to here with
          --assign-code x:C meaning associate modification code "x" with
          cytosine (C) primary sequence bases. If a code is encountered that is
          not part of the specification, the bedMethyl record will not be used,
          this will be logged

      --single-code <SINGLE_MODIFICATION_CODE>
          Calculate differential methylation for this modification code in
          isolation

  -k, --mask
          Respect soft masking in the reference FASTA

      --min-valid-coverage <MIN_VALID_COVERAGE>
          Minimum valid coverage required to use an entry from a bedMethyl. See
          the help for pileup for the specification and description of valid
          coverage
          
          [default: 0]

Output Options:
  -o, --out-path <OUT_PATH>
          Path to file to direct output, optional, no argument will direct
          output to stdout

      --header
          Include header in output

Segmentation Options:
      --segment <SEGMENTATION_FP>
          Run segmentation, output segmented differentially methylated regions
          to this file

      --max-gap-size <MAX_GAP_SIZE>
          Maximum number of base pairs between modified bases for them to be
          segmented together
          
          [default: 5000]

      --dmr-prior <DMR_PRIOR>
          Prior probability of a differentially methylated position
          
          [default: 0.1]

      --diff-stay <DIFF_STAY>
          Maximum probability of continuing a differentially methylated block,
          decay will be dynamic based on proximity to the next position
          
          [default: 0.9]

      --significance-factor <SIGNIFICANCE_FACTOR>
          Significance factor, effective p-value necessary to favor the
          "Different" state
          
          [default: 0.01]

      --log-transition-decay
          Use logarithmic decay for "Different" stay probability

      --decay-distance <DECAY_DISTANCE>
          After this many base pairs, the transition probability will become the
          prior probability of encountering a differentially modified position
          
          [default: 500]

      --fine-grained
          Preset HMM segmentation parameters for higher propensity to switch
          from "Same" to "Different" state. Results will be shorter segments,
          but potentially higher sensitivity

Logging Options:
      --careful
          Log out which sequences are in common between the samples and the
          reference FASTA, useful for debugging

      --log-filepath <LOG_FILEPATH>
          File to write logs to, it's recommended to use this option

      --suppress-progress
          Don't show progress bars

      --missing <HANDLE_MISSING>
          How to handle regions found in the `--regions` BED file. quiet =>
          ignore regions that are not found in the tabix header warn => log
          (debug) regions that are missing fatal => log (error) and exit the
          program when a region is missing
          
          [default: quiet]
          [possible values: quiet, warn, fail]

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of threads to use when for decompression
          
          [default: 4]

      --batch-size <BATCH_SIZE>
          Control the  batch size. The batch size is the number of regions to
          load at a time. Each region will be processed concurrently. Loading
          more regions at a time will decrease IO to load data, but will use
          more memory. Default will be 50% more than the number of threads
          assigned

  -f, --force
          Force overwrite of output file, if it already exists

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead
          
          [default: 100000]

Single-site Options:
      --prior <PRIOR> <PRIOR>
          Prior distribution for estimating MAP-based p-value. Should be two
          arguments for alpha and beta (e.g. 1.0 1.0). See
          `dmr_scoring_details.md` for additional details on how the metric is
          calculated. Alpha and beta must each be positive, and their sum must
          be >= 1.0

      --delta <DELTA>
          Consider only effect sizes greater than this when calculating the
          MAP-based p-value
          
          [default: 0.05]

  -N, --n-sample-records <N_SAMPLE_RECORDS>
          Sample this many reads when estimating the max coverage thresholds
          
          [default: 10042]

      --max-coverages <MAX_COVERAGES> <MAX_COVERAGES>
          Max coverages to enforce when calculating estimated MAP-based p-value

      --cap-coverages
          When using replicates, cap coverage to be equal to the maximum
          coverage for a single sample. For example, if there are 3 replicates
          with max_coverage of 30, the total coverage would normally be 90.
          Using --cap-coverages will down sample the data to 30X
```

## dmr multi
```text
Compare regions between all pairs of samples (for example a trio sample set or
haplotyped trio sample set). As with `pair` all inputs must be bgzip compressed
bedMethyl files with associated tabix indices. Each sample must be assigned a
name. Output is a directory of BED files with the score column indicating the
magnitude of the difference in methylation between the two samples indicated in
the file name. See the online documentation for additional details

Usage: modkit dmr multi [OPTIONS] --regions-bed <REGIONS_BED> --out-dir <OUT_DIR> --ref <REFERENCE_FASTA>

Options:
  -h, --help  Print help

Sample Options:
  -s, --sample <SAMPLES> <SAMPLES>
          Two or more named samples to compare. Two arguments are required
          <path> <name>. This option should be repeated at least two times. When
          two samples have the same name, they will be combined
  -r, --regions-bed <REGIONS_BED>
          BED file of regions over which to compare methylation levels. Should
          be tab-separated (spaces allowed in the "name" column). Requires
          chrom, chromStart and chromEnd. The Name column is optional. Strand is
          currently ignored
      --ref <REFERENCE_FASTA>
          Path to reference fasta for the pileup
  -m, --base <MODIFIED_BASES>
          Bases to use to calculate DMR, may be multiple. For example, to
          calculate differentially methylated regions using only cytosine
          modifications use --base C
      --assign-code <MOD_CODE_ASSIGNMENTS>
          Extra assignments of modification codes to their respective primary
          bases. In general, modkit dmr will use the SAM specification to know
          which modification codes are appropriate to use for a given primary
          base. For example "h" is the code for 5hmC, so is appropriate for
          cytosine bases, but not adenine bases. However, if your bedMethyl file
          contains custom codes or codes that are not part of the specification,
          you can specify which primary base they belong to here with
          --assign-code x:C meaning associate modification code "x" with
          cytosine (C) primary sequence bases. If a code is encountered that is
          not part of the specification, the bedMethyl record will not be used,
          this will be logged
      --single-code <SINGLE_MODIFICATION_CODE>
          Calculate differential methylation for this modification code in
          isolation
  -k, --mask
          Respect soft masking in the reference FASTA
      --min-valid-coverage <MIN_VALID_COVERAGE>
          Minimum valid coverage required to use an entry from a bedMethyl. See
          the help for pileup for the specification and description of valid
          coverage [default: 0]
      --ref-sample <REF_SAMPLE>
          Only performs comparisons between this sample and all other samples

Output Options:
      --header             Include header in output
  -o, --out-dir <OUT_DIR>  Directory to place output DMR results in BED format
  -p, --prefix <PREFIX>    Prefix files in directory with this label
  -f, --force              Force overwrite of output file, if it already exists

Logging Options:
      --log-filepath <LOG_FILEPATH>
          File to write logs to, it's recommended to use this option
      --suppress-progress
          Don't show progress bars
      --missing <HANDLE_MISSING>
          How to handle regions found in the `--regions` BED file. quiet =>
          ignore regions that are not found in the tabix header warn => log
          (debug) regions that are missing fatal => log (error) and exit the
          program when a region is missing [default: quiet] [possible values:
          quiet, warn, fail]

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use [default: 4]
      --io-threads <IO_THREADS>
          Number of threads to use when for decompression [default: 4]
```

## dmr isoform
```text
Use a transcriptome-aligned bedMethyl table and GTF file to compare methylation
across gene transcript isoforms. Run analysis on the entire transcriptome or on
a single gene. Optionally also make plots for single genes

Usage: modkit dmr isoform [OPTIONS] --gtf <GTF> <INPUT_BEDMETHYL> <OUTPUT_TABLE>

Arguments:
  <INPUT_BEDMETHYL>
          Input bedmethyl table. Expected to be bgzip-compressed and have a .tbi
          tabix index. This bedMethyl table must be generated from a
          transcriptome-aligned modBAM, NOT a genome-aligned modBAM. Expected
          that this table was generated by pileup with `--modified-bases`
          argument

  <OUTPUT_TABLE>
          Path to output BED file or '-' for standard output

Options:
      --gtf <GTF>
          Path to GTF file to use for transcript models. May be compressed

      --min-valid-coverage <MIN_VALID_COVERAGE>
          Minimum valid coverage required to use an entry from a bedMethyl. See
          the help for pileup for the specification and description of valid
          coverage
          
          [default: 1]

      --single-code <SINGLE_MODIFICATION_CODE>
          Calculate differential methylation for this modification code in
          isolation. Unmodified and other modification counts will be combined

      --full
          Include per-modification combined proportions per position and
          per-isoform proportions in output. Adding this flag will make the
          output file substantially larger

      --ignore-version
          Ignore gene ID version and transcript ID version in GTF, if there are
          multiple identical transcript IDs or gene IDs with different versions
          the last one will be kept

      --log-filepath <LOG_FILEPATH>
          Path to optional debug log (recommended)

  -g, --gene-id <GENE_ID>
          Only process this Gene Id

  -n, --gene-name <GENE_NAME>
          Only process this Gene Name

      --gene-strand <GENE_STRAND>
          Specify a strand for the Gene Name (or GeneID) if there are multiple
          chromosomes with this gene

      --gene-chrom <GENE_CHROM>
          Specify a chromosome for the Gene Name (or GeneID) if there are
          multiple chromosomes with this gene

      --suppress-progress
          Don't show the progress bar

  -t, --threads <THREADS>
          Number of concurrent processors to use. Only used when running on all
          genes
          
          [default: 4]

      --plot <PLOT>
          Generate an SVG plot displaying the methylation positions on each
          isoform (filtered by --min-score or --max-pval) and the common exons
          of the gene. Also displays negative log of the p-value of the DMR
          metric. See the online documentation for an example

      --marker-positions <MARKER_POSITIONS>...
          Comma-separated list of genomic positions to draw a vertical dotted
          line and a label

      --width <WIDTH>
          SVG width in px
          
          [default: 1600]

      --track-height <TRACK_HEIGHT>
          Height per transcript track in px
          
          [default: 26]

      --track-gap <TRACK_GAP>
          Gap between transcript tracks in px
          
          [default: 18]

      --left-margin <LEFT_MARGIN>
          Left margin in px
          
          [default: 220]

      --right-margin <RIGHT_MARGIN>
          Right margin in px
          
          [default: 40]

      --top-margin <TOP_MARGIN>
          Top margin in px
          
          [default: 210]

      --bottom-margin <BOTTOM_MARGIN>
          Bottom margin in px
          
          [default: 50]

      --intron-px <INTRON_PX>
          Length, in px, of the compressed intron segments
          
          [default: 24]

      --top-bar-height <TOP_BAR_HEIGHT>
          Height in px of the top bar
          
          [default: 90]

      --min-score <MIN_SCORE>
          Only plot methylation markers with negative log p-value greater than
          this value

      --max-pval <MAX_PVAL>
          Only plot methylation markers with p-value less than this value

  -h, --help
          Print help (see a summary with '-h')
```

## dmr compare-tx-sites
```text
Uses two transcriptome-aligned bedmMethyl tables to compare methylation for two
conditions at transcriptome sites aggregated at the gene level. This command
will combine the methylation counts across isoforms into a single summary
statistic per-site for the gene, then compare those counts between conditions at
the gene level and report results on genomic coordinates. The transcript models
are provided by a GTF file. Optionally make a volcano plot of the top-k
differentially methylated sites per gene

Usage: modkit dmr compare-tx-sites [OPTIONS] -a <COND_A> -b <COND_B> --out <OUTPUT_TABLE> --gtf <GTF>

Options:
  -a <COND_A>
          Input bedmethyl table for first condition (condition 'A'). Expected to
          be bgzip-compressed and have a .tbi tabix index. This bedMethyl table
          must be generated from a transcriptome-aligned modBAM, NOT a
          genome-aligned modBAM. Expected that this table was generated by
          pileup with `--modified-bases` argument
  -b <COND_B>
          Input bedmethyl table for first condition (condition 'B')
  -o, --out <OUTPUT_TABLE>
          Path to output BED file or '-' for stdout
      --full
          Include combined proportions and counts for each condition. This will
          increase the size of the output file
      --gtf <GTF>
          Path to GTF file to use for transcript models. Can be compressed
      --ignore-version
          Ignore gene ID version and transcript ID version in GTF, if there are
          multiple identical transcript IDs or gene IDs with different versions
          the last one will be kept
      --min-valid-coverage <MIN_VALID_COVERAGE>
          Minimum valid coverage required to use an entry from a bedMethyl. See
          the help for pileup for the specification and description of valid
          coverage [default: 1]
      --single-code <SINGLE_MODIFICATION_CODE>
          Calculate differential methylation for this modification code in
          isolation. Unmodified and other modification counts will be combined
      --log-filepath <LOG_FILEPATH>
          Path to optional debug log (recommended)
      --suppress-progress
          Don't show the progress bar
  -t, --threads <THREADS>
          Number of concurrent processors to use [default: 4]
      --plot <PLOT>
          Path to file to create an SVG volcano plot where the x-axis is log2
          fold change and the y-axis is negative log p-value. See the online
          documentation for and example
      --top-k <TOP_K>
          For each gene, sort the positions in increasing p-value and take this
          many for plotting. Setting a large number will increase the size of
          the SVG plot [default: 1]
      --label-top-k-genes <LABEL_TOP_K_GENES>
          Sort the points to be plotted by increasing p-value, then label the
          points from the top-k genes
      --min-effect-size <MIN_EFFECT_SIZE>
          Discard positions with an effect size less than this before
          considering them for plotting [default: 0.1]
      --plot-title <PLOT_TITLE>
          Add a title to the plot
      --gene-labels <GENE_LABELS>
          Path to plain text file of a list of gene names, one per line, to mark
          in the volcano plot
      --sort-by-effect-size
          Determine which points to plot by sorting to decreasing effect size
          instead of increasing p-value
  -h, --help
          Print help
```

## bedmethyl merge
```text
Perform an outer join on two or more bedMethyl files, summing their counts for
records that overlap

Usage: modkit bedmethyl merge [OPTIONS] --out-bed <OUT_BED> --genome-sizes <GENOME_SIZES> [IN_BEDMETHYL] [IN_BEDMETHYL]...

Arguments:
  [IN_BEDMETHYL] [IN_BEDMETHYL]...
          Input bedMethyl table(s). Should be bgzip-compressed and have an
          associated Tabix index. The tabix index will be assumed to be
          $this_file.tbi

Options:
  -g, --genome-sizes <GENOME_SIZES>
          TSV of genome sizes, should be <chrom>\t<size_in_bp>

  -h, --help
          Print help (see a summary with '-h')

Output Options:
  -o, --out-bed <OUT_BED>
          Specify the output file to write the results table

      --force
          Force overwrite the output file

      --header
          Output a header with the bedMethyl

      --mixed-delim
          Output bedMethyl where the delimiter of columns past column 10 are
          space-delimited instead of tab-delimited. This option can be useful
          for some browsers and parsers that don't expect the extra columns of
          the bedMethyl format

Compute Options:
      --chunk-size <CHUNK_SIZE>
          Chunk size for how many start..end regions for each chromosome to
          read. Larger values will lead to faster merging at the expense of
          memory usage, while smaller values will be slower with lower memory
          usage. This option will only impact large bedmethyl files

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead
          
          [default: 100000]

  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --queue-size <QUEUE_SIZE>
          Number of batches (of size chunk size) allowed to be in a pre-written
          state at once. Increasing this number will increase memory usage
          
          [default: 30]

      --io-threads <IO_THREADS>
          Number of tabix/bgzf threads to use
          
          [default: 2]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file to write debug logs to

Filtering Options:
      --min-samples <MIN_SAMPLES>
          Only output a position present in at least this many input bedMethyl
          files. Accepts an integer greater than 1, or "all" to require the
          position in every input (an inner join across replicates)

      --min-sample-coverage <MIN_SAMPLE_COVERAGE>
          Minimum valid coverage for an input's record to count towards a
          position, for both the --min-samples tally and the summed counts
```

## bedmethyl tobigwig
```text
Make a BigWig track from a bedMethyl file or stream. For details on the BigWig
format see https://doi.org/10.1093/bioinformatics/btq351

Usage: modkit bedmethyl tobigwig [OPTIONS] --mod-codes <MOD_CODES> <--sizes <CHROMSIZES>|--header <INPUT_BAM>> <IN_BEDMETHYL> <OUT_FP>

Arguments:
  <IN_BEDMETHYL>
          Input bedmethyl, uncompressed, "-" or "stdin" indicates an input
          stream

  <OUT_FP>
          Output bigWig filename

Options:
  -g, --sizes <CHROMSIZES>
          A chromosome sizes file. Each line should be a chromosome and its size
          in bases, separated by whitespace. A fasta index (.fai) works as well.
          Use instead of the bam header

  -b, --header <INPUT_BAM>
          modBAM from which the pileup was generated. Chromosome sizes are
          gathered from the header. Use instead of the chromosome sizes file

  -m, --mod-codes <MOD_CODES>
          Make a bigWig track where the values are the percent of bases with
          this modification, use multiple comma-separated codes to combine
          counts. For example --mod-code m makes a track of the 5mC percentages
          and --mod-codes h,m will make a track of the combined counts from 5hmC
          and 5mC. Combining counts for different primary bases will cause an
          error (e.g. --mod-codes a,h will error)

  -h, --help
          Print help (see a summary with '-h')

Output Options:
      --negative-strand-values
          Report the percentages on the negative strand as negative values. The
          data range will be [-100, 100]

  -z, --nzooms <NZOOMS>
          Set the maximum of zooms to create
          
          [default: 10]

      --zooms <ZOOMS>...
          Set the zoom resolutions to use (overrides the --nzooms argument)

  -u, --uncompressed
          Don't use compression

      --block-size <BLOCK_SIZE>
          Number of items to bundle in r-tree
          
          [default: 256]

      --items-per-slot <ITEMS_PER_SLOT>
          Number of data points bundled at lowest level
          
          [default: 1024]

Compute Options:
  -t, --nthreads <NTHREADS>
          Set the number of threads to use. This tool will typically use ~225%
          CPU on a HDD. SDDs may be higher. (IO bound)
          
          [default: 6]

      --inmemory
          Do not create temporary files for intermediate data

      --force-chromosome-ordering
          If input bedMethyl has sorting of the same scheme as `sort`, this
          option may speed up conversion

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended. (alias: log)

      --suppress-progress
          Hide the progress bar
```

## bedmethyl map-to-genome
```text
Map a transcriptome-aligned bedMethyl file to genome-coordinates based on a GTF

Usage: modkit bedmethyl map-to-genome [OPTIONS] --transcript-id <TRANSCRIPT_ID> --gtf <GTF> <INPUT_BEDMETHYL> <OUTPUT_BEDMETHYL>

Arguments:
  <INPUT_BEDMETHYL>   
  <OUTPUT_BEDMETHYL>  

Options:
  -t, --transcript-id <TRANSCRIPT_ID>
          
      --transcript-version <TRANSCRIPT_VERSION>
          
      --ignore-version
          Ignore gene ID version and transcript ID version in GTF, if there are
          multiple identical transcript IDs or gene IDs with different versions
          the last one will be kept
      --gtf <GTF>
          
      --header
          
  -h, --help
          Print help
```

## modbam adjust-mods
```text
Performs various operations on BAM files containing base modification
information, such as converting base modification codes and ignoring
modification calls. Produces a BAM output file

Usage: modkit modbam adjust-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
          Input BAM file, can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input

  <OUT_BAM>
          File path to new BAM file to be created. Can be a path to a file or
          one of `-` or `stdin` to specify a stream from standard output

Options:
  -f, --ff
          Fast fail, stop processing at the first invalid sequence record.
          Default behavior is to continue and report failed/skipped records at
          the end

  -h, --help
          Print help (see a summary with '-h')

Output Options:
      --log-filepath <LOG_FILEPATH>
          Output debug logs to file at this path

      --output-sam
          Output SAM format instead of BAM

Modified Base Options:
      --ignore <IGNORE>
          Modified base code to ignore/remove, see
          https://samtools.github.io/hts-specs/SAMtags.pdf for details on the
          modified base codes

      --convert <CONVERT> <CONVERT>
          Convert one mod-tag to another, summing the probabilities together if
          the retained mod tag is already present

      --motif <MOTIF> <MOTIF>
          Filter out any base modification call that isn't part of a basecall
          sequence motif. This argument can be passed multiple times. Format is
          <motif_sequence> <offset>. For example the argument to match CpG
          dinucleotides is `--motif CG 0`, or to match CG[5mC]G the argument
          would be `--motif CGCG 2`. Single bases can be used as motifs to keep
          only base modification calls for a specific primary base, for example
          `--motif C 0`

      --cpg
          Shorthand for --motif CG 0

      --discard-motifs
          Discard base modification calls that match the provided motifs
          (instead of keeping them)

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

Selection Options:
      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

      --filter-probs
          Filter out the lowest confidence base modification probabilities

      --mapped-only
          Only use base modification probabilities from bases that are aligned
          when estimating the filter threshold (i.e. ignore soft-clipped, and
          inserted bases)

Sampling Options:
  -n, --num-reads <NUM_READS>
          Sample approximately this many reads when estimating the filtering
          threshold. If alignments are present reads will be sampled evenly
          across aligned genome. If a region is specified, either with the
          --region option or the --sample-region option, then reads will be
          sampled evenly across the region given. This option is useful for
          large BAM files. In practice, 10-50 thousand reads is sufficient to
          estimate the model output distribution and determine the filtering
          threshold
          
          [default: 10042]

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold
          probability. If this option is not provided, but --region is provided,
          the genomic interval passed to --region will be used. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size to process concurrently when estimating the
          threshold probability, can be larger than the pileup processing
          interval
          
          [default: 1000000]

Filtering Options:
  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per primary base. A global
          filter threshold can be specified with by a decimal number (e.g.
          0.75). Per-base thresholds can be specified by colon-separated values,
          for example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

Logging Options:
      --suppress-progress
          Hide the progress bar
```

## modbam update-tags
```text
Renames Mm/Ml to tags to MM/ML. Also allows changing the mode flag from silent
'.' to explicitly '?' or '.'

Usage: modkit modbam update-tags [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM to update modified base tags in. Can be a path to a file or one
             of `-` or `stdin` to specify a stream from standard input
  <OUT_BAM>  File to new BAM file to be created or one of `-` or `stdin` to
             specify a stream from standard output

Options:
  -m, --mode <MODE>        Mode, change mode to this value, options {'explicit',
                           'implicit'}. See spec at:
                           https://samtools.github.io/hts-specs/SAMtags.pdf.
                           'explicit' ('?') means residues without modification
                           probabilities will not be assumed canonical or
                           modified. 'implicit' means residues without explicit
                           modification probabilities are assumed to be
                           canonical [possible values: explicit, implicit]
      --no-implicit-probs  Don't add implicit canonical calls. This flag is
                           important when converting from one of the implicit
                           modes ( `.` or `""`) to explicit mode (`?`). By
                           passing this flag, the bases without associated base
                           modification probabilities will not be assumed to be
                           canonical. No base modification probability will be
                           written for these bases, meaning there is no
                           information. The mode will automatically be set to
                           the explicit mode `?`
  -h, --help               Print help

Compute Options:
  -t, --threads <THREADS>  Number of threads to use [default: 4]

Logging Options:
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path

Output Options:
      --output-sam  Output SAM format instead of BAM
```

## modbam sample-probs
```text
Calculate an estimate of the base modification probability distribution

Usage: modkit modbam sample-probs [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input BAM with modified base tags. If a index is found reads will be
          sampled evenly across the length of the reference sequence. Can be a
          path to a file or one of `-` or `stdin` to specify a stream from
          standard input

Options:
      --ignore-index
          Ignore the BAM index (if it exists) perform a serial scan of the BAM,
          otherwise if an index is found, multiple workers will be used to read
          the BAM in parallel. Conflicts with --threads since there will only be
          one reader

  -r, --reference <REFERENCE_FASTA>
          Reference sequence in FASTA format. (alias: 'ref')

      --preload-references
          Preload the reference sequences, useful when working with many, short
          reference sequences such as a transcriptome

  -k, --mask
          Respect soft masking in the reference FASTA

  -h, --help
          Print help (see a summary with '-h')

Compute Options:
  -t, --threads <THREADS>
          Number of workers to run concurrently. When an aligned, indexed BAM is
          found, multiple workers (threads) will be spawned to read disjoint
          sections of the BAM in parallel
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of IO/decompression threads to use, these are shared across
          workers (if an indexed BAM is found) or used for a single reader if
          --ignore-index is passed, stdin is used, or an index is not found
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller
          interval chunk sizes will use less memory but incur more overhead.
          Only used when sampling probs from an indexed bam
          
          [default: 1000000]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended

      --suppress-progress
          Hide the progress bar

Output Options:
  -p, --percentiles <PERCENTILES>
          Percentiles to calculate, a space separated list of floats (for
          example: "10,20,90")
          
          [default: 10,50,90]

  -q, --quantiles <QUANTILES>
          Quantiles to calculate, a space separated list of floats (for example:
          "0.1,0.2,0.9")

  -o, --out-dir <OUT_DIR>
          Directory to deposit result tables into. Required for model
          probability histogram output

      --prefix <PREFIX>
          Label to prefix output files with

      --force
          Overwrite results if present

      --hist
          Output histogram of base modification prediction probabilities

      --dna-color <PRIMARY_BASE_COLORS> <PRIMARY_BASE_COLORS>
          Set colors of primary bases in histogram, should be RGB format, e.g.
          "#0000FF" is defailt for canonical cytosine

      --mod-color <MOD_BASE_COLORS> <MOD_BASE_COLORS>
          Set colors of modified bases in histogram, should be RGB format, e.g.
          "#FF00FF" is default for 5hmC

Selection Options:
      --motif <MOTIF> <MOTIF>
          Tabulate base modification probabilities at motifs (or single bases).
          The first argument should be the sequence motif and the second
          argument is the 0-based offset to the base to pileup base modification
          counts for. For example: --motif CGCG 0 indicates to collect base
          mofification probabilities for the first C on the top strand and the
          last C (complement to G) on the bottom strand. For single bases,
          --motif C 0 or --motif A 0 would restrict to probabilities on C or A
          reference bases, respectively. This argument can be passed multiple
          times

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --region <REGION>
          Process only the specified region of the BAM when collecting
          probabilities. Format should be <chrom_name>:<start>-<end> or
          <chrom_name>

      --include-bed <INCLUDE_BED>
          Only sample base modification probabilities that are aligned to the
          positions in this BED file. (alias: include-positions)

      --use-non-primary
          Use non-primary mappings, may cause double-counting of reads with
          primary and one or more non-primary mappings. Requires MN tag to be
          correct. See documentation for help on proper alignment

Sampling Options:
  -n, --num-reads <NUM_READS>
          

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of
          reads to sample, for example 0.1 will sample 1/10th of the reads

      --no-sampling
          No sampling, use all of the reads to calculate the filter thresholds

  -s, --seed <SEED>
          Random seed for deterministic running, the default is
          non-deterministic, only used when no BAM index is provided
```

## modbam summary
```text
Summarize the mod tags present in a BAM and get basic statistics. The default
output is a totals table (designated by '#' lines) and a modification calls
table. Descriptions of the columns can be found in the README

Usage: modkit modbam summary [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input modBam, can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input

Options:
      --ignore-index
          Ignore the BAM index (if it exists) perform a serial scan of the BAM,
          otherwise if an index is found, multiple workers will be used to read
          the BAM in parallel. Conflicts with --threads since there will only be
          one reader

  -r, --reference <REFERENCE_FASTA>
          Reference sequence in FASTA format. (alias: 'ref')

      --preload-references
          Preload the reference sequences, useful when working with many, short
          reference sequences such as a transcriptome

  -k, --mask
          Respect soft masking in the reference FASTA

      --cpg
          Shorthand for --motif CG 0. Requires a sorted, indexed modBAM

  -h, --help
          Print help (see a summary with '-h')

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --io-threads <IO_THREADS>
          Number of IO/decompression threads to use, these are shared across
          workers (if an indexed BAM is found) or used for a single reader if
          --ignore-index is passed, stdin is used, or an index is not found
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          When using regions, interval chunk size in base pairs to process
          concurrently. Smaller interval chunk sizes will use less memory but
          incur more overhead
          
          [default: 1000000]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended

      --suppress-progress
          Hide the progress bar

Output Options:
      --tsv
          Output summary as a tab-separated variables stdout instead of a table

      --combine-mods
          Combine modification counts per primary sequence base when calculating
          summary statistics

Sampling Options:
  -n, --num-reads <NUM_READS>
          Approximate maximum number of reads to use. If an indexed BAM is
          provided, the reads will be sampled evenly over the length of the
          aligned reference. If a region is passed with the --region option,
          they will be sampled over the region. Actual number of reads used may
          deviate slightly from this number

  -f, --sampling-frac <SAMPLING_FRAC>
          Instead of using a defined number of reads, specify a fraction of
          reads to sample when estimating the filter threshold. For example 0.1
          will sample 1/10th of the reads

  -s, --seed <SEED>
          Sets a random seed for deterministic running (when using
          --sample-frac)

Filtering Options:
  -p, --filter-quantile <FILTER_QUANTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence base modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter
          threshold can be specified with by a decimal number (e.g. 0.75).
          Per-base thresholds can be specified by colon-separated values, for
          example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

      --max-filter-threshold <MAX_FILTER_THRESHOLD>
          Maximum threhsold value to use, if the estimated threshold value for
          any primary sequence base is greater than this value, use this value
          instead. To set a default value for all primary bases, use
          --max-threshold 0.8, for example. Or use per-base threshold values
          like --max-threshold C:0.8

Selection Options:
      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --include-bed <INCLUDE_BED>
          Only summarize base modification probabilities that are aligned to the
          positions in this BED file. (alias: include-positions)

      --mapped-only
          Only use modBAM records that are mapped, don't include unmapped
          records

      --matched-only
          Only summarize base modifications that are mapped and are matches to
          the appropriate primary sequence base. For example, require that 5mC
          base modification calls be aligned to a reference C. Requires a
          reference FASTA and an sorted, indexed modBAM

      --motif <MOTIF> <MOTIF>
          Calculate base modification statistics at reference motifs. This
          option has the same behavior as "--matched-only" except it only
          includes specified motifs. The first argument should be the sequence
          motif and the second argument is the 0-based offset to the base to
          pileup base modification counts for. For example: --motif CGCG 0
          indicates to collect base mofification probabilities for the first C
          on the top strand and the last C (complement to G) on the bottom
          strand. For single bases, --motif C 0 or --motif A 0 would restrict to
          probabilities on C or A reference bases, respectively. This argument
          can be passed multiple times. Requires a sorted, indexed modBAM

      --region <REGION>
          Process only the specified region of the BAM when collecting
          probabilities. Format should be <chrom_name>:<start>-<end> or
          <chrom_name>
```

## modbam call-mods
```text
Call mods from a modbam, creates a new modbam with probabilities set to 100% if
a base modification is called or 0% if called canonical

Usage: modkit modbam call-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
          Input BAM, may be sorted and have associated index available. Can be a
          path to a file or one of `-` or `stdin` to specify a stream from
          standard input

  <OUT_BAM>
          Output BAM, can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input

Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended

      --ff
          Fast fail, stop processing at the first invalid sequence record.
          Default behavior is to continue and report failed/skipped records at
          the end

      --suppress-progress
          Hide the progress bar

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently
          
          [default: 4]

  -n, --num-reads <NUM_READS>
          Sample approximately this many reads when estimating the filtering
          threshold. If alignments are present reads will be sampled evenly
          across aligned genome. If a region is specified, either with the
          --region option or the --sample-region option, then reads will be
          sampled evenly across the region given. This option is useful for
          large BAM files. In practice, 10-50 thousand reads is sufficient to
          estimate the model output distribution and determine the filtering
          threshold
          
          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the
          filter-percentile. In practice, 50-100 thousand reads is sufficient to
          estimate the model output distribution and determine the filtering
          threshold. See filtering.md for details on filtering

      --seed <SEED>
          Set a random seed for deterministic running, the default is
          non-deterministic, only used when no BAM index is provided

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold
          probability. If this option is not provided, but --region is provided,
          the genomic interval passed to --region will be used. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size to process concurrently when estimating the
          threshold probability, can be larger than the pileup processing
          interval
          
          [default: 1000000]

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted
          variant is below this confidence percentile. For example, 0.1 will
          filter out the 10% lowest confidence modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per primary base. A global
          filter threshold can be specified with by a decimal number (e.g.
          0.75). Per-base thresholds can be specified by colon-separated values,
          for example C:0.75 specifies a threshold value of 0.75 for cytosine
          modification calls. Additional per-base thresholds can be specified by
          repeating the option: for example --filter-threshold C:0.75
          --filter-threshold A:0.70 or specify a single base option and a
          default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for
          adenine and 0.9 for all other base modification calls

      --mod-threshold <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification,
          independent of the threshold for the primary sequence base or the
          default. For example, to set the pass threshold for 5hmC to 0.8 use
          `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and other modifications unless
          the `--filter-threshold` option is also passed. See the online
          documentation for more details

      --no-filtering
          Don't filter base modification calls, assign each base modification to
          the highest probability prediction

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the
          start or the end of the read. Two comma-separated values may be
          provided to asymmetrically filter out base modification calls from the
          start and end of the reads. For example, 4,8 will filter out base
          modification calls in the first 4 and last 8 bases of the read

      --invert-edge-filter
          Invert the edge filter, instead of filtering out base modification
          calls at the ends of reads, only _keep_ base modification calls at the
          ends of reads. E.g. if usually, "4,8" would remove (i.e. filter out)
          base modification calls in the first 4 and last 8 bases of the read,
          using this flag will keep only base modification calls in the first 4
          and last 8 bases

      --motif <MOTIF> <MOTIF>
          Filter out any base modification call that isn't part of a basecall
          sequence motif This argument can be passed multiple times. Format is
          <motif_sequence> <offset>. For example the argument to match CpG
          dinucleotides is `--motif CG 0`, or to match CG[5mC]G the argument
          would be `--motif CGCG 2`

      --cpg
          Shorthand for --motif CG 0

      --discard-motifs
          Discard base modification calls that match the provided motifs
          (instead of keeping them)

      --output-sam
          Output SAM format instead of BAM

  -h, --help
          Print help (see a summary with '-h')
```

## modbam repair
```text
Repair MM and ML tags in one bam with the correct tags from another. To use this
command, both modBAMs _must_ be sorted by read name. The "donor" modBAM's reads
must be a superset of the acceptor's reads. Extra reads in the donor are
allowed, and multiple reads with the same name (secondary, etc.) are allowed in
the acceptor. Reads with an empty SEQ field cannot be repaired and will be
rejected. Reads where there is an ambiguous alignment of the acceptor to the
donor will be rejected (and logged). See the full documentation for details

Usage: modkit modbam repair [OPTIONS] --donor-bam <DONOR_BAM> --acceptor-bam <ACCEPTOR_BAM> --output-bam <OUTPUT_BAM>

Options:
  -d, --donor-bam <DONOR_BAM>
          Donor modBAM with original MM/ML tags. Must be sorted by read name
  -a, --acceptor-bam <ACCEPTOR_BAM>
          Acceptor modBAM with reads to have MM/ML base modification data
          projected on to. Must be sorted by read name
  -o, --output-bam <OUTPUT_BAM>
          output modBAM location
      --log-filepath <LOG_FILEPATH>
          File to write logs to, it is recommended to use this option as some
          reads may be rejected and logged here
  -t, --threads <THREADS>
          The number of threads to use [default: 4]
  -h, --help
          Print help
```

## modbam check-tags
```text
Usage: modkit modbam check-tags [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>
          Input modBam, can be a path to a file or one of `-` or `stdin` to
          specify a stream from standard input

Options:
      --permissive
          Don't exit 1 when invalid records are found in the input

  -h, --help
          Print help (see a summary with '-h')

IO Options:
  -o, --out-dir <OUT_DIR>
          Write output tables into this directory. The directory will be created
          if it doesn't exist

  -f, --force
          Force overwrite of previous output

      --prefix <PREFIX>
          Prefix output files with this string

Compute Options:
  -t, --threads <THREADS>
          Number of threads to use
          
          [default: 4]

      --ignore-index
          Perform a linear scan of the modBAM even if the index is found

  -i, --interval-size <INTERVAL_SIZE>
          When using regions, interval chunk size in base pairs to process
          concurrently. Smaller interval chunk sizes will use less memory but
          incur more overhead
          
          [default: 5000000]

Logging Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them.
          Setting a file is recommended

      --suppress-progress
          Hide the progress bar

Selection Options:
  -n, --num-reads <NUM_READS>
          Approximate maximum number of reads to use, especially recommended
          when using a large BAM without an index. If an indexed BAM is
          provided, the reads will be sampled evenly over the length of the
          aligned reference. If a region is passed with the --region option,
          they will be sampled over the region. Actual number of reads used may
          deviate slightly from this number

      --allow-non-primary
          Check tags on non-primary alignments as well. Keep in mind this may
          incur a double-counting of the read with its primary mapping

      --mapped-only
          Only check alignments that are mapped

      --region <REGION>
          Process only the specified region of the BAM when collecting
          probabilities. Format should be <chrom_name>:<start>-<end> or
          <chrom_name>

      --head <HEAD>
          Use the first N reads without sampling. Shorthand for --ignore-index
          --num-reads N
```
