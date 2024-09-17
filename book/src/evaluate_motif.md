# Evaluate a table of known motifs

The `modkit search` command has an option to provide any number of known motifs with `--know-motif`.
If you already have a list of candidate motifs (e.f. from a previous run of `modkit motif search`) you can check these motifs quickly against a bedMethyl table with `modkit motif evaluate`.

```bash
modkit motif evaluate -i ${bedmethyl} --known-motifs-table motifs.tsv -r ${ref}
```

Similarly, the search [algorithm](./intro_find_motifs.md#simple-description-of-the-search-algorithm) can be run using known motifs as seeds:

```bash
modkit motif refine -i ${bedmethyl} --known-motifs-table motifs.tsv -r ${ref}
```

The output tables to both of these commands have the same schema:

| column | name       | description                                                                                     | type  |
|--------|------------|-------------------------------------------------------------------------------------------------|-------|
| 1      | mod_code   | code specifying the modification found in the motif                                             | str   |
| 2      | motif      | sequence of identified motif using [IUPAC](https://www.bioinformatics.org/sms/iupac.html) codes | str   |
| 3      | offset     | 0-based offset into the motif sequence of the modified base                                     | int   |
| 4      | frac_mod   | fraction of time this sequence is found in the _high modified_ set col-5 / (col-5 + col-6)      | float |
| 5      | high_count | number of occurances of this sequence in the _high-modified_ set                                | int   |
| 6      | low_count  | number of occurances of this sequence in the _low-modified_ set                                 | int   |
| 7      | mid_count  | number of occurances of this sequence in the _mid-modified_ set                                 | int   |
| 8      | log_odds   | log2 odds of the motif being in the high-modified set                                           | int   |

In the human-readable table columns (1) and (2) are merged to show the modification code in the motif sequence context, the rest of the columns are the same as the machine-readable table.

