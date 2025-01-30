# ModBam Utilities

## Checking MM/ML tags

In general, when Modkit encounters an record that it cannot use due to the modified base information being invalid or not conforming to the [specification](https://samtools.github.io/hts-specs/SAMtags.pdf) it will log the error and ignore the record.
However, it may be important to know precisely how many records are unusable and for what reason.
The `modkit modbam check-tags` command is for this purpose.
It will tabulate the error conditions in the modBAM and report out how many of each problem it encounters.
When there are any records that cause an error, the command will exit 1 (unless the `--permissive` flag is given).
You can use this command as part of a workflow to exit early if a modBAM is invalid.

### Example usage
```bash
modkit modbam check-tags ${modbam} -t 20 --log ${outdir}/check_tags.log
# alternatively use the abbreviated subcommand
modkit mb check-tags ${modbam} -t 20 --log ${outdir}/check_tags.log
```

The command outputs the following tables

1. Error counts, each error encountered and the number of records with the error. The `pct` column shows the percent of errors. 
```text
+-----------------------------------------+-------+-------+
| error                                   | count | pct   |
+-----------------------------------------+-------+-------+
| conflict-explicit-and-inferred          | 3792  | 92.71 |
| conflict-explicit-prob-greater-than-one | 298   | 7.29  |
| total                                   | 4090  | 100   |
+-----------------------------------------+-------+-------+
```

2. Tag-header tables, counts of how many of each MM-tag "header" are found, e.g. "C+m." means cytosine is the fundamental base, on the positive strand, and the "mode" is '.' meaning that cytosines without explicit probabilities are inferred to be un-modified.
See the specification for more details on the "mode".
There is one table for each the valid and invalid records.

```text
+------------+-------+
| tag_header | count |
+------------+-------+
| C+m.       | 4090  |
| A+a.       | 4090  |
| N+e?       | 4090  |
| C+h.       | 4090  |
+------------+-------+
```

3. Modified bases, a table of which primary sequence bases have which modification codes assigned to them as well as the "mode" in the tag corresponding to that assignment.

```text
+--------+--------------+----------+------+
| strand | primary_base | mod_code | mode |
+--------+--------------+----------+------+
| +      | T            | b        | ?    |
| +      | T            | e        | ?    |
| +      | G            | b        | ?    |
| +      | G            | e        | ?    |
| +      | C            | b        | .    |
| +      | C            | e        | .    |
| +      | C            | h        | .    |
| +      | C            | m        | .    |
| +      | A            | a        | .    |
| +      | A            | b        | .    |
| +      | A            | e        | .    |
+--------+--------------+----------+------+
```
