# Current limitations

Known limitations and forecasts for when they will be removed.

1. Ambiguous DNA bases in ML tags are not supported (for example `N+m?`).
   - This limitation will be removed in version 0.2.z
2. During `modkit pileup`, it is assumed that each read should only have one primary alignment. If a read name
   is detected more than once, the occurance is logged but both alignments will be used. This limitation may be
   removed in the future with a form of dynamic de-duplication.
3. Only one MM-flag (`.`, `?`) per-canonical base is supported within a read.
    - This limitation may be removed in the future.
