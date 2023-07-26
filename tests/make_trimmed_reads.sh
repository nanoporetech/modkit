#!/bin/bash

# this script makes some reads that have their ends trimmed

set -eu

python ./trim_reads.py \
  resources/bc_anchored_10_reads.sorted.bam \
  ./trimmed.bam \
  --head-trim 25 \
  --tail-trim 30 

samtools fastq trimmed.bam > trimmed.fastq

minimap2 -ax map-ont \
  ./resources/CGI_ladder_3.6kb_ref.fa \
  trimmed.fastq \
  > trimmed.mapped.sam

samtools sort -n trimmed.mapped.sam -O BAM > trimmed_read_sort.mapped.bam 
samtools sort -O BAM -n resources/bc_anchored_10_reads.sorted.bam > donor_read_sort.bam