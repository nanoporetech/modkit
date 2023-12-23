#!/bin/bash

abort() {
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "Test failure" >&2
    exit 1
}


function realpath () {
    python -c "import os; print(os.path.realpath('${1}'))"
}

trap 'abort' 0

set -eu

cargo check

#    test_extract_correct_output
cargo run extract \
  tests/resources/bc_anchored_10_reads.sorted.bam \
  tests/resources/bc_anchored_10_reads.sorted.methylprofile.tsv \
  -i 25 --force

#    test_extract_collapse_correct_output
cargo run extract \
  tests/resources/bc_anchored_10_reads.sorted.bam \
  tests/resources/bc_anchored_10_reads.sorted.methylprofile_ignoreh.tsv \
  -i 25 --force --ignore h

#    test_extract_calls_regression
cargo run extract tests/resources/2_reads_all_context.bam null \
  --read-calls tests/resources/test_read_calls_estimate_thresh.tsv \
  --ref tests/resources/CGI_ladder_3.6kb_ref.fa \
  --force

#    test_extract_correct_output_with_ref
cargo run extract \
  tests/resources/bc_anchored_10_reads.sorted.bam \
  tests/resources/bc_anchored_10_reads.sorted.methylprofile_ref.tsv \
  -i 25 --force --ref tests/resources/CGI_ladder_3.6kb_ref.fa

#    test_extract_duplex_correct_output
cargo run extract \
  tests/resources/duplex_modbam.sorted.bam \
  tests/resources/duplex_sorted.tsv \
  --force --region chr17

#    test_extract_implicit_mod_calls
cargo run extract \
  tests/resources/implicit_mod_tags.bam \
  tests/resources/extract_with_implicit.tsv \
  --force

#    test_extract_include_sites_duplex_regression
cargo run extract \
  tests/resources/duplex_modbam.sorted.bam \
  tests/resources/test_extract_include_sites_duplex_regression_expected.tsv \
  --ignore-index --force \
  --include-bed tests/resources/hg38_chr17_CG0_snip.bed

trap : 0
echo "> all done"

