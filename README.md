# modkit

> Tools for handling modBAM files and pileups

```bash
Usage: modkit <COMMAND>

Commands:
  adjust-mods   Collapse N-way base modification calls to (N-1)-way
  update-tags   Update mod tags, changes Mm/Ml-style tags to MM/ML-style. Also allows to change the mode to '?' or '.' instead of implicitly '.'
  pileup        Pileup (combine) mod calls across genomic positions
  sample-probs  Get an estimate of the distribution of mod-base prediction probabilities
  summary       Summarize the mod tags present in a BAM and get basic statistics
  motif-bed     Create BED file with all locations of a motif
  help          Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help information

```

## adjust-mods
```bash
Collapse N-way base modification calls to (N-1)-way

Usage: modkit adjust-mods [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM file to collapse mod call from
  <OUT_BAM>  File path to new BAM file

Options:
      --ignore <IGNORE>              mod base code to ignore/remove [default: h]
  -t, --threads <THREADS>            number of threads to use [default: 4]
  -f, --ff                           Fast fail, stop processing at the first invalid sequence record. Default behavior is to continue and report failed/skipped records at the end
      --method <METHOD>              Method to use to collapse mod calls, 'norm', 'dist'. A full description of the methods can be found in collapse.md [default: norm]
      --convert <CONVERT> <CONVERT>  Convert one mod-tag to another, summing the probabilities together if the retained mod tag is already present
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path
  -h, --help                         Print help information
```

## update-tags
```bash
Update mod tags, changes Mm/Ml-style tags to MM/ML-style. Also allows to change the mode to '?' or '.' instead of implicitly '.'

Usage: modkit update-tags [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>   BAM file to collapse mod call from
  <OUT_BAM>  File path to new BAM file

Options:
  -m, --mode <MODE>                  Mode, change mode to this value, options {'ambiguous', 'implicit'}. See spec at: https://samtools.github.io/hts-specs/SAMtags.pdf. 'ambiguous' ('?') means residues without explicit modification probabilities will not be assumed canonical or modified. 'implicit' means residues without explicit modification probabilities are assumed to be canonical [possible values: ambiguous, implicit]
  -t, --threads <THREADS>            number of threads to use [default: 4]
      --log-filepath <LOG_FILEPATH>  Output debug logs to file at this path
  -h, --help                         Print help information
```

## pileup
```bash
Pileup (combine) mod calls across genomic positions. Produces bedMethyl formatted file. Schema and description of fields can be found in schema.yaml

Usage: modkit pileup [OPTIONS] <IN_BAM> <OUT_BED>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index

  <OUT_BED>
          Output file

Options:
  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently
          
          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size to process concurrently. Smaller interval chunk sizes will use less memory but incur more overhead
          
          [default: 100000]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the `filter-percentile`. In practice, 50-100 thousand reads is sufficient to estimate the model output distribution and determine the filtering threshold
          
          [default: 0.1]

      --seed <SEED>
          random seed for deterministic running, default is non-deterministic

      --no-filtering
          Do not perform any filtering, include all mod base calls in output

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter (remove) mod-calls where the probability of the predicted variant is below this percentile. For example, 0.1 will filter out the lowest 10% of modification calls
          
          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Filter threshold, drop calls below this probability

      --log-filepath <LOG_FILEPATH>
          Output debug logs to file at this path

      --combine
          Combine mod calls, all counts of modified bases are summed together

      --collapse <COLLAPSE>
          Collapse _in_situ_. Arg is the method to use {'norm', 'dist'}

      --method <METHOD>
          Method to use to collapse mod calls, 'norm', 'dist'. A full description of the methods can be found in collapse.md
          
          [default: norm]

      --force-allow-implicit
          Force allow implicit-canonical mode. By default modkit does not allow pileup with the implicit mode ('.', or omitted). The `update-tags` subcommand is provided to update tags to the new mode, however if the user would like to assume that residues with no probability associated canonical, this option will allow that behavior

  -h, --help
          Print help information (use `-h` for a summary)
```

## sample-probs
```bash
Get an estimate of the distribution of mod-base prediction probabilities

Usage: modkit sample-probs [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>  Input BAM, should be sorted and have associated index

Options:
  -f, --sampling-frac <SAMPLING_FRAC>  Sample fraction [default: 0.1]
  -t, --threads <THREADS>              number of threads to use reading BAM [default: 4]
  -s, --seed <SEED>                    random seed for deterministic running, default is non-deterministic
  -p, --percentiles <PERCENTILES>      Percentiles to calculate, space separated list [default: 0.1,0.5,0.9]
  -h, --help                           Print help information
```

## summary
```bash
Summarize the mod tags present in a BAM and get basic statistics

Usage: modkit summary [OPTIONS] <IN_BAM>

Arguments:
  <IN_BAM>  Input ModBam file

Options:
  -t, --threads <THREADS>  number of threads to use reading BAM [default: 4]
  -h, --help               Print help information
```

## motif-bed
```bash
Create BED file with all locations of a motif

Usage: modkit motif-bed <FASTA> <MOTIF> <OFFSET>

Arguments:
  <FASTA>   Input FASTA file
  <MOTIF>   Motif to search for within FASTA
  <OFFSET>  Offset within motif

Options:
  -h, --help  Print help information
```
