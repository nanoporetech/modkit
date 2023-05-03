# Current limitations

Known limitations and forecasts for when they will be removed.

1. ChEBI codes are not supported at all, only mod-codes in the
[specification](https://samtools.github.io/hts-specs/SAMtags.pdf) (top of page 9) are
supported in `pileup`.
    - This limitation will be removed by version 0.2.0

1. Only one MM-flag (`.`, `?`) per-canonical base is supported within a read.
    - This limitation may be removed in the future.
