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

test_dir=/tmp/scratch
mkdir -p ${test_dir}
test_file=/tmp/scratch/test_extract.tsv
control_file=/tmp/scratch/control_extract.tsv

#    test_extract_correct_output
cargo run extract tests/resources/bc_anchored_10_reads.sorted.bam - -i 25 --force \
  | awk -v OFS="\t" 'NR>1' | cut -f 1-18 | sort -k 1 > ${test_file}
awk -v OFS="\t" 'NR>1' tests/resources/bc_anchored_10_reads.sorted.methylprofile.tsv | sort -k 1 > ${control_file}
diff ${test_file} ${control_file}

#    test_extract_collapse_correct_output
cargo run extract tests/resources/bc_anchored_10_reads.sorted.bam - \
  -i 25 --force --ignore h \
  | awk -v OFS="\t" 'NR>1' | cut -f 1-18 | sort -k 1 > ${test_file}
awk -v OFS="\t" 'NR>1' tests/resources/bc_anchored_10_reads.sorted.methylprofile_ignoreh.tsv | sort -k 1 > ${control_file}
diff ${test_file} ${control_file}


#    test_extract_calls_regression
cargo run extract tests/resources/2_reads_all_context.bam null \
  --read-calls /tmp/scratch/test_calls_tmp.tsv \
  --ref tests/resources/CGI_ladder_3.6kb_ref.fa \
  --force
awk -v OFS="\t" 'NR>1' /tmp/scratch/test_calls_tmp.tsv | cut -f 1-20 | sort -k 1 > ${test_file}
awk -v OFS="\t" 'NR>1' tests/resources/test_read_calls_estimate_thresh.tsv | sort -k 1 > ${control_file}
diff ${test_file} ${control_file}


#    test_extract_correct_output_with_ref
cargo run extract tests/resources/bc_anchored_10_reads.sorted.bam - \
  -i 25 --force --ref tests/resources/CGI_ladder_3.6kb_ref.fa \
  | awk -v OFS="\t" 'NR>1' | cut -f 1-18 | sort -k 1 > ${test_file}

awk -v OFS="\t" 'NR>1' tests/resources/bc_anchored_10_reads.sorted.methylprofile_ref.tsv | sort -k 1 > ${control_file}
diff ${test_file} ${control_file}

#    test_extract_duplex_correct_output
cargo run extract tests/resources/duplex_modbam.sorted.bam - \
  --region chr17 \
  | awk -v OFS="\t" 'NR>1' | cut -f 1-18 | sort -k 1 > ${test_file}
awk -v OFS="\t" 'NR>1' tests/resources/duplex_sorted.tsv | sort -k 1 > ${control_file}
diff ${test_file} ${control_file}

#    test_extract_implicit_mod_calls
cargo run extract tests/resources/implicit_mod_tags.bam - --force \
  | awk -v OFS="\t" 'NR>1' | cut -f 1-18 | sort -k 1 > ${test_file}
awk -v OFS="\t" 'NR>1' tests/resources/extract_with_implicit.tsv | sort -k 1 > ${control_file}
diff ${test_file} ${control_file}

#    test_extract_include_sites_duplex_regression
cargo run extract tests/resources/duplex_modbam.sorted.bam - \
  --ignore-index --force --include-bed tests/resources/hg38_chr17_CG0_snip.bed \
  | awk -v OFS="\t" 'NR>1' | cut -f 1-18 | sort -k 1 > ${test_file}
awk -v OFS="\t" 'NR>1' tests/resources/test_extract_include_sites_duplex_regression_expected.tsv | sort -k 1 > ${control_file}
diff ${test_file} ${control_file}


trap : 0
echo "> all good!"

