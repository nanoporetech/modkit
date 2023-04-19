#!/bin/bash

set -e

samtools view -h resources/bc_anchored_10_reads.sorted.bam \
    | sed s/C+m?,/C+Z,/g \
    | sed s/C+h?,/C+Y,/g \
    | samtools view -b -o resources/bc_anchored_10_reads_old_tags.bam

