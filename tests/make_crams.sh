#!/usr/bin/env bash
set -e

pushd resources
samtools view -C -T CGI_ladder_3.6kb_ref.fa bc_anchored_10_reads.sorted.bam > bc_anchored_10_reads.sorted.cram
samtools index bc_anchored_10_reads.sorted.cram

samtools view -C -T CGI_ladder_3.6kb_ref.fa bc_anchored_10_reads.unmapped.bam > bc_anchored_10_reads_unmapped.cram
samtools index bc_anchored_10_reads_unmapped.cram
popd