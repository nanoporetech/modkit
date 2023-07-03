#!/bin/bash


if [ -z $1 ]; then
  echo "USAGE: $0 old_tsv"
  exit 1
fi

set -eu

function realpath () {
    python -c "import os; print(os.path.realpath('${1}'))"
}

old=$(realpath $1)

pushd ../
cargo b --release
popd

modkit=$(realpath ../target/release/modkit)

${modkit} \
  extract \
  /Users/art.rand/projects/mod_flatten/tests/resources/duplex_modbam.sorted.bam \
  a.tsv \
  --include-bed ../tests/resources/hg38_chr17_CG0_snip.bed \
  --force


awk 'NR>1 {print}' ${old} | sort -k 1,2 -s > old_sort.tsv

awk 'NR>1 {print}' a.tsv | cut -f 1-6,8-16 | sort -k 1,2 -s > b.tsv

echo "checking $old vs b.tsv"

diff b.tsv old_sort.tsv && echo "all good" || echo "they're different"


