#!/bin/bash
#SBATCH --time=08:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=priitpaluoja@gmail.com
#SBATCH --mem=14000
#SBATCH -J SeqFFMM
# The job requires 1 compute node
#SBATCH -N 1
# The job requires 1 task per node
#SBATCH --ntasks-per-node=1
# Number of CPU cores per task
#SBATCH --cpus-per-task=10

module load samtools-0.1.19
module load python-3.6.3
module load R-3.2.0

# Build index
# bowtie2-build GCA_000001405.15_GRCh38_genomic.fna bwtie38
SOFTBWT=/gpfs/hpchome/ppaluoja/software/bowtie2-2.3.3.1/bowtie2
REF=/gpfs/rocket/samba/CCHT/BelgiaNIPT/fastq/bwtie38
locations="locations"

# Read locations from file and calculate FF.
while IFS='' read -r line || [[ -n "$line" ]]; do
    location=$line
    sample=${location##*/}
    echo $sample
    zcat --force $location/*.fastq* | $SOFTBWT --very-sensitive -X 500 -q - --norc -x $REF --no-unal -p 10 --mm | samtools view -q 10 -S - | python3 separator.py $sample "50000" | Rscript seqff.R >> results.all.10.ff.tsv
done < "$locations"

echo "DONE!"
