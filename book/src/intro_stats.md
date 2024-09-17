# Calculating modification statistics in regions

There are many analysis operations available in `modkit` once you've generated a bedMethyl table.
One such operation is to calculate aggregation statistics on specific regions, for example in CpG islands or gene promoters.
The `modkit stats` command is designed for this purpose. 

```bash
# these files can be found in the modkit repository
cpgs=tests/resources/cpg_chr20_with_orig_names_selection.bed
sample=tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz
modkit stats ${sample} --regions ${cpgs} -o ./stats.tsv [--mod-codes "h,m"]
```

> Note that the argument `--mod-codes` can alternatively be passed multiple times, e.g. this is equivalent: <br />
> `--mod-codes c --mod-codes h`

The output TSV has the following schema:

| column | Name           | Description                                                                   | type  |
|--------|----------------|-------------------------------------------------------------------------------|-------|
| 1      | chrom          | name of reference sequence from BAM header                                    | str   |
| 2      | start position | 0-based start position                                                        | int   |
| 3      | end position   | 0-based exclusive end position                                                | int   |
| 4      | name           | name of the region from input BED (`.` if not provided)                       | str   |
| 5      | strand         | Strand (`+`, `-`, `.`) from the input BED (`.` assumed for when not provided) | str   |
| 6+     | count_x        | total number of `x` base modification codes in the region                     | int   |
| 7+     | count_valid_x  | total valid calls for the primary base modified by code `x`                   | int   |
| 8+     | percent_x      | `count_x` / `count_vali_x` * 100                                              | float |

Columns 6, 7, and 8 are repeated for each modification code found in the bedMethyl file or provided with `--mod-codes` argument.

An example output:

```text
chrom  start     end       name     strand   count_h        count_valid_h  percent_h   count_m        count_valid_m  percent_m
chr20  9838623   9839213   CpG: 47  .        12             1777           0.6752954   45             1777           2.532358
chr20  10034962  10035266  CpG: 35  .        7              1513           0.46265697  0              1513           0
chr20  10172120  10172545  CpG: 35  .        15             1229           1.2205045   28             1229           2.278275
chr20  10217487  10218336  CpG: 59  .        29             2339           1.2398461   108            2339           4.617358
chr20  10433628  10434345  CpG: 71  .        29             2750           1.0545455   2              2750           0.07272727
chr20  10671925  10674963  CpG: 255 .        43             9461           0.45449743  24             9461           0.25367296
```

