#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=priitpaluoja@gmail.com
#SBATCH --mem=5000
#SBATCH -J SeqSAM1
# The job requires 1 compute node
#SBATCH -N 1
# The job requires 1 task per node
#SBATCH --ntasks-per-node=1
# Number of CPU cores per task
#SBATCH --cpus-per-task=10

# Bash script authors: Priit Paluoja, Priit Palta, Hindrek Teder, Kaarel Krjutshkov. (CCHT)
#
# Usage:
# 1. Provide system with R and Anaconda with source "ff" which contains bowtie2 and samtools. (conda install bowtie2 and conda install samtools)
# 2. Download and modify seqFF to accept piped sam file.
# 3. Download bowtie2 hg19 pre-built index file and update bash script with corresponding location
# 4. Fill file "locations" with folder paths to samples. All sample X files must be in one folder. Locations can contain any number of samples.
# 5. Run the script. Last tested on Rocket cluster (rocket.hpc.ut.ee) on 26.05.2018.


module load R-3.2.0
module load python-3.6.3

source activate ff # Source most provide bowtie2 and samtools

# Pre-built Bowtie 2 index (must be GRCh37/hg19)
REF=/gpfs/hpchome/ppaluoja/software/databases/hg19
# File which provides path(s) to fastq/fastq.gz files. For each sample, all the sample files must be in the same directory. Locations file can contain multiple samples.
locations="locations1.txt"

# Read locations from file and calculate FF.
# 1. Concatenate fastq
# 2. Map
# 3. Filter
# 4. Calculate FF
while IFS='' read -r line || [[ -n "$line" ]]; do
    location=$line
    sample=${location##*/}
    echo $sample
    # Align and filter
    zcat --force $location/*.fastq* | bowtie2 --very-sensitive -X 500 -q - --norc -x $REF --no-unal -p 10 | samtools view -q 30 -S - > $sample.aligned.sam
    # Sort and index for gender determination
    samtools view -Sb $sample.aligned.sam | samtools sort > $sample.aligned.sorted.bam
    samtools index $sample.aligned.sorted.bam
    # Get statistics; keep chromosome names and read counts; chromosomes 13,18,21,X; last line / (sum of all lines - last line) -> ratio of Y and remaining chromosomes
    gender=$(cat $sample.aligned.sorted.bam | samtools idxstats | cut -f1,3 | head -n 24 | sed '/chr13\|chr18\|chr21\|chrX/d' |  awk '{ sum += $2 } END { print $2/(sum-$2) }')
    # Calculate fetal fraction + append gender ratio
    cat $sample.aligned.sam | cut -f3,4 | Rscript seqff.sam.R | (printf "$sample\t$gender\t" && cat) >> results.sam.tsv
    rm $sample.aligned.sam
    rm $sample.aligned.sorted.bam
    rm $sample.aligned.sorted.bai
done < "$locations"


echo "DONE!"
