#!/usr/bin/env bash
# Create test paired-end FASTQs containing only reads that map to chr1 (or a chr1 region).
# Use these outputs in the pipeline's fastq/ dir (with naming that matches config) to test the workflow quickly.
#
# Usage:
#   ./scripts/make_chr1_test_fastq.sh <input.bam> <output_prefix> [region]
#
# Arguments:
#   input.bam     BAM file (e.g. from a previous run: alignment/bam/sample.bowtie2.bam)
#   output_prefix Base name for outputs; will create output_prefix.1.fq.gz and output_prefix.2.fq.gz
#   region        Optional. SAM region to keep (default: chr1). Examples: chr1  or  chr1:1-5000000
#
# Requires: samtools
#
# Example - chr1 only, write into pipeline fastq dir:
#   ./scripts/make_chr1_test_fastq.sh alignment/bam/WT_CTCF_Rep1.bowtie2.bam fastq/test_chr1
#
# Example - first 5 Mb of chr1 (smaller test):
#   ./scripts/make_chr1_test_fastq.sh alignment/bam/WT_CTCF_Rep1.bowtie2.bam fastq/test_chr1_small chr1:1-5000000
#
# Then add the sample to config (e.g. samples: ["test_chr1"]) or ensure fastq/test_chr1.1.fq.gz and .2.fq.gz exist for discovery.

set -euo pipefail
BAM="${1:?Usage: $0 input.bam output_prefix [region]}"
PREFIX="${2:?Usage: $0 input.bam output_prefix [region]}"
REGION="${3:-chr1}"

echo "Extracting region=$REGION from $BAM -> ${PREFIX}.1.fq.gz / ${PREFIX}.2.fq.gz"
samtools view -b -F 0x904 "$BAM" "$REGION" | \
  samtools sort -n -o - | \
  samtools fastq -1 "${PREFIX}.1.fq.gz" -2 "${PREFIX}.2.fq.gz" -
echo "Done."
