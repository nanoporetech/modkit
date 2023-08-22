# Current limitations

Known limitations and forecasts for when they will be removed.

1. ChEBI codes are not supported at all, only mod-codes in the
   [specification](https://samtools.github.io/hts-specs/SAMtags.pdf) (top of page 9) are supported in
   `pileup`.
    - This limitation will be removed by version 0.2.0
2. Ambiguous DNA bases in ML tags are not supported (for example `N+m?`).
   - This limitation will be removed in version 0.2.z
3. Only one MM-flag (`.`, `?`) per-canonical base is supported within a read.
    - This limitation may be removed in the future.
    - If you have previously run a transformation that removes base modification calls, for example `modkit
      adjust-mods` with the `--edge-filter` option. Following `modkit extract` commands should be run with 
      the same `--edge-filter` option or else the trimmed base modifications will be output as canonical. 
4. During `modkit pileup`, it is assumed that each read should only have one primary alignment. If a read name
is detected more than once, the occurance is logged but both alignments will be used. This limitation may be
removed in the future with a form of dynamic de-duplication.
