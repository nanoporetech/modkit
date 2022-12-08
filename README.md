# mod_flatten

```bash
flattens N-way base-modification calls into N-1

Usage: mod_flatten [OPTIONS] <IN_BAM> <OUT_BAM>

Arguments:
  <IN_BAM>
  <OUT_BAM>

Options:
  -b, --base <BASE>          canonical base to flatten calls for [default: C]
  -m, --mod-base <MOD_BASE>  mod base code to flatten/remove [default: h]
  -t, --threads <THREADS>    number of threads to use [default: 1]
  -f, --ff                   exit on bad reads, otherwise continue
  -h, --help                 Print help information
```

## Example
```bash
cargo run --release -- tests/resources/ref_out_5mC_auto.bam flattened_5mC.bam

# or if installed

mod_flatten tests/resources/ref_out_5mC_auto.bam flattened_5mC.bam
```

## TODO
- [ ] C-API
- [ ] python lib
