#!/bin/bash

# Bash script authors: Priit Paluoja, Priit Palta, Hindrek Teder, Kaarel Krjutshkov. (CCHT)
#
# Usage:
# 1. Provide system with Bowtie 2, Samtools and R.
# 2. Download and modify seqFF to accept piped sam file.
# 3. Download Bowtie 2 hg19 pre-built index file and update bash script with corresponding location.
# 4. All sample X files must be in one folder.

INPUT="$1"  # (un)compressed FASTQ
OUTPUT="$2"  # output directory
CPUS="$3"  # number of CPUs
MEM="$4"M  # memory per CPU in megabytes
REF=software/SeqFF/reference/hg19  # pre-built Bowtie 2 index (must be GRCh37/hg19)
SAMPLE="$(basename "$INPUT" .fastq)"  # sample name


bowtie2 -x "$REF" -U "$INPUT" --very-sensitive -p "$CPUS" |  # aligne
samtools view -u -q 30 - |  # convert and filter
samtools sort -m "$MEM" -o "$OUTPUT"/aligned.filtered.sorted.bam -  # sort
samtools index "$OUTPUT"/aligned.filtered.sorted.bam  # index

# Get statistics; keep chromosome names and read counts; chromosomes 13, 18, 21, X; last line / (sum of all lines - last line) -> ratio of Y and remaining chromosomes
GENDER="$(samtools idxstats "$OUTPUT"/aligned.filtered.sorted.bam | cut -f 1,3 | head -n 24 | sed '/chr13\|chr18\|chr21\|chrX/d' | awk '{ sum += $2 } END { print $2 / (sum - $2) }')"

# Calculate fetal fraction + append gender ratio
samtools view "$OUTPUT"/aligned.filtered.sorted.bam | cut -f 3,4 | Rscript software/SeqFF/seqff.sam.R | (printf "$SAMPLE\t$GENDER\t" && cat) > "$OUTPUT"/fetal_fraction.tsv
