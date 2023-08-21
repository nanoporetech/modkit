"""Script was used to update the test data in test_extract.rs
"""
import os
import shutil
import pandas as pd


if __name__ == '__main__':
    files = [
        "tests/resources/bc_anchored_10_reads.sorted.methylprofile_ignoreh.tsv",
        "tests/resources/bc_anchored_10_reads.sorted.methylprofile.tsv",
        "tests/resources/bc_anchored_10_reads.sorted.methylprofile_ref.tsv",
        "tests/resources/duplex_sorted.tsv",
        "tests/resources/test_extract_include_sites_duplex_regression_expected.bed",
    ]
    for fp in files:
        print(f"> updating {fp}")
        assert os.path.exists(fp)
        backup = f"{fp}.bak"
        shutil.copy(fp, backup)
        assert os.path.exists(backup)

        X = pd.read_csv(backup, delim_whitespace=True)
        X["inferred"] = "false"
        X.to_csv(fp, sep="\t", index=None)

    print("> done")

