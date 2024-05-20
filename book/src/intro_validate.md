# Validating ground truth results.

The `modkit validate` sub-command is intended for validating results in a uniform manner from samples with known modified base content. Specifically the modified base status at any annotated reference location should be known.

## Validating from modBAM reads and BED reference annotation.

The input to the `modkit validate` command will be pairs of modBAM and BED files.
modBAM files should contain modified base calls in the MM/ML tag as input to most modkit commands.
BED files paired to each input modBAM file describe the ground truth modified base status at reference positions.
This ground truth status should be known by the researcher due to previous experimental conditions and cannot be derived by modkit.

```
modkit validate \
    --bam-and-bed sample1.bam sample1_annotation.bed \
    --bam-and-bed sample2.bam sample2_annotation.bed
```

This will produce output such as the following:

```text
> Parsing BED at /path/to/sample1_annotation.bed
> Processed 10 BED lines
> Parsing BED at /path/to/sample2_annotation.bed
> Processed 10 BED lines
> Canonical base: C
> Parsing mapping at /path/to/sample1.bam
> Processed 10 mapping recrods
> Parsing mapping at /path/to/sample2.bam
> Processed 10 mapping recrods
> Raw counts summary
              Called Base
         ┌───┬───────┬───────┬───┬───┬───┬──────────┐
         │   │ -     │ a     │ C │ G │ T │ Deletion │
         ├───┼───────┼───────┼───┼───┼───┼──────────┤
 Ground  │ - │ 9,900 │   100 │ 1 │ 1 │ 1 │        2 │
 Truth   │ a │   100 │ 9,900 │ 1 │ 1 │ 1 │        2 │
         └───┴───────┴───────┴───┴───┴───┴──────────┘
> Balancing ground truth call totals
> Raw accuracy: 99.00%
> Raw modified base calls contingency table
              Called Base
         ┌───┬────────┬────────┐
         │   │ -      │ a      │
         ├───┼────────┼────────┤
 Ground  │ - │ 99.00% │  1.00% │
 Truth   │ a │  1.00% │ 99.00% │
         └───┴────────┴────────┘
> Call probability threshold: 0.95
> Filtered accuracy: 99.90%
> Filtered modified base calls contingency table
              Called Base
         ┌───┬────────┬────────┐
         │   │ -      │ a      │
         ├───┼────────┼────────┤
 Ground  │ - │ 99.90% │  0.10% │
 Truth   │ a │  0.10% │ 99.90% │
         └───┴────────┴────────┘
```

The filtering threshold is computed in the same manner as in all other modkit commands.
Currently only a defined percentage of input data filtering threshold estimation is implemented.
The default value is 10% (as in other modkit commands) and can be adjusted with the `--filter-quantile` argument.
Other methods (including user-defined thresholds) will be implemented in a future version.

The `Call probability threshold` is intended a value to be used for user-defined thresholds for other modkit commands.

## BED ground truth annotation file:

A BED file is a tab-delimited file.
For this command the first 6 fields are processed.
These fields are as follows:

| column | name           | description                    |     |
|--------|----------------|--------------------------------|-----|
| 1      | chrom          | name of reference sequence     | str |
| 2      | start position | 0-based start position         | int |
| 3      | end position   | 0-based exclusive end position | int |
| 4      | mod code       | modified base code             | str |
| 6      | strand         | strand (e.g. +,-,.)            | str |

The 5th column is ignored in the validate command.

The 4th column represents the modified base code annotating the status at this reference position (or range of reference positions).
This value can be `-` representing a canonical base (note that this differs from the `remora validate` annotation), a single letter code as defined in the modBAM tag specification, or any ChEBI code.
The validate command will assume that any base from the associated modBAM file overlapping these positions should match this annotation.

## Output file

The `--out-filepath` option is provided to allow persistent storage of results in a machine-parseable format without other logging lines.
This format outputs all contingency tables in a machine-parseable format.
For example this contingency table `[["ground_truth_label","-","a"],["-",9900,100],["a",100,9900]]` would be produced from the above example results.
