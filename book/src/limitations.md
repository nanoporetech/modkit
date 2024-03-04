# Current limitations

Known limitations and forecasts for when they will be removed.

1. Ambiguous DNA bases in ML tags are not supported (for example `N+m?`).
   - This limitation will be removed in version 0.2.z
2. During `modkit pileup`, it is assumed that each read should only have one primary alignment. If a read name
   is detected more than once, the occurrence is logged but both alignments will be used. This limitation may be
   removed in the future with a form of dynamic de-duplication.
3. Only one MM-flag (`.`, `?`) per-canonical base is supported within a read.
    - This limitation may be removed in the future.
4. The MAP-based p-value metric ([details](./dmr_scoring_details.md#map-based-p-value)) performs a test that there is a difference in modification (of any kind) between two conditions.
If a position has multiple base modification calls (such as 5hmC and 5mC) the calls are summed together into a single "modified" count. 
If a position differs only in the modification _type_ (such as one condition has more 5hmC and the other has more 5mC) this effect will not be captured in the MAP-based p-value. The [likelihood ratio](./dmr_scoring_details.md#likelihood-ratio-scoring-details) test _does_ capture changes in modification type.
5. The MAP-based p-value is not available when performing DMR on regions.
This is because there is potentially large variability in the number of modified bases and coverage over regions. This variability translates into varying degrees of statistical power and makes comparisons difficult to interpret.
