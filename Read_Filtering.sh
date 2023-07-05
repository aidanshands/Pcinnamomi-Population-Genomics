#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 8
#SBATCH --time=5-00:00:00
##SBATCH --output=my.stdout
#SBATCH --mail-type=ALL
#SBATCH --job-name="Fastq-Filter"
#SBATCH -p batch


source activate ea-utils

cd $SLURM_SUBMIT_DIR

SAMPFILE=Isolates.csv
RAWREADDIR=Raw_Reads
STATDIR=Filtered_Stats
FILTERREADDIR=Filtered_Reads

IFS=,
# Filering the reads & Getting the stats
tail -n +2 $SAMPFILE | while read ISOLATE
do
  fastq-mcf \
  -q 20 -o $FILTERREADDIR/$ISOLATE.F.1.fq.gz -o $FILTERREADDIR/$ISOLATE.F.2.fq.gz \
  Adapters.fasta \
  $RAWREADDIR/$ISOLATE.1.fq.gz \
  $RAWREADDIR/$ISOLATE.2.fq.gz
  # filtered stats
  fastq-stats $FILTERREADDIR/$ISOLATE.F.1.fq.gz > $STATDIR/$ISOLATE.stats.1.txt
  fastq-stats $FILTERREADDIR/$ISOLATE.F.2.fq.gz > $STATDIR/$ISOLATE.stats.2.txt
done
