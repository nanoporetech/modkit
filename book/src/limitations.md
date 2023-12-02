# Current limitations

Known limitations and forecasts for when they will be removed.

1. Ambiguous DNA bases in ML tags are not supported (for example `N+m?`).
   - This limitation will be removed in version 0.2.z
2. During `modkit pileup`, it is assumed that each read should only have one primary alignment. If a read name
   is detected more than once, the occurrence is logged but both alignments will be used. This limitation may be
   removed in the future with a form of dynamic de-duplication.
3. Only one MM-flag (`.`, `?`) per-canonical base is supported within a read.
    - This limitation may be removed in the future.
4. Functions that transform a modBAM into another modBAM (and manipulate the MM and ML tags) can only do so 
   with the primary alignments. Supplementary and secondary alignments will not be present in the output.
   There are plans to remove this limitation in the near future.
