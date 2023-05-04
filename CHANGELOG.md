# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.1.5]
### Changes
- Table output is default in `modkit summary`, columns are all_counts/frac and pass_counts/frac instead of counts/frac and filt_counts/frac. See help for details.

### Fixes
- The `modkit summary` command will not use all of the reads in the input modBAM when the `--no-filtering` flag is used. To disable sampling (and use all of the reads) use `--no-sampling`. See docs for details.

### Adds
- The `sample-probs` command will now output histograms of the base modification probabilities with the `--hist` option.
- "Book-style" docs.


## [v0.1.4]
### Changes
- Implicitly canonical mode `.` See [SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf) does not require `--force-allow-implicit`. However no mode at all still does.
- `modkit adjust-mods` will allow conversion of base modification codes that aren't in the [specification](https://samtools.github.io/hts-specs/SAMtags.pdf) (e.g. `Z`). In order to process these with `modkit pileup` you must first run `modkit adjust-mods in.bam out.bam --convert Z m` (for example). 
- Per-canonical base thresholds will be automatically inferred, e.g. different thresholds will be chosen for A mods and C mods. This is also true for `modkit summary` and `modkit sample-probs`.
- `modkit summary` will now allow specification of a `--region` when using an indexed BAM and will sample evenly over the aligned contigs when the region is omitted.
### Adds
- `modkit summary` has `--table` option that produces a human-readable table.
### Deprecated
- `modkit summary` tab-separated format will not be default in the next release. The table output (currently specified with the `--table` option will become the default).
### Fixes
- Restructuring of the code allows for more parallelism when estimating threshold values, speeding up this process at the start of `modkit pileup` or `modkit summary`.

## [v0.1.4-rc1]

### Changes
- Filter base modification calls that are < the threshold value instead of <= the threshold value.


## [v0.1.3]
### Adds
- `--mask` flag in `pileup` and `motif-bed`

### Fixes
- Default behavior will find CGs in masked regions of the reference.


## [v0.1.2]
### Adds
- First release.
    
### Fixes

