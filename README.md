# mod_kit

> Tools for handling BAMs with mod bases


# `modkit `

```bash
Usage: modkit <COMMAND>

Commands:
  collapse  Collapse N-way base modification calls to (N-1)-way
  pileup    Pileup (combine) mod calls across genomic positions
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help information
```
## TODO (all)
- [ ] CI/CD that runs tests and builds artifacts
- [ ] python lib

## TODO (`collapse`)
- [ ] Don't remove sup/dup/secondary alignments
- [ ] Handle removing mod-prob in ambiguous mode (`.`).
- [ ] duplex

## TODO (`pileup`)
- [x] Proper command line
- [ ] logging?
  - Should log to a file
- [ ] collapse
- [x] estimation of threshold
- [ ] duplex
- [ ] GC in read cache
- [x] test on larger dataset
- [ ] bigwig output
