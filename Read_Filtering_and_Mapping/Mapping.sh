#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=200G
#SBATCH -c 64
#SBATCH --time=12-00:00:00
##SBATCH --output=my.stdout
#SBATCH --mail-type=ALL
#SBATCH --job-name="BWA"
#SBATCH -p batch

# Load modules
module load samtools/1.14
module load bwa/0.7.17

cd $SLURM_SUBMIT_DIR
# Setting variables
SAMPFILE=Isolates.csv
FILTERREADDIR=Filtered_Reads
SORTEDCRAM=Sorted_CRAMs
# Indexing genome
bwa index -p Pc2113T1 Pc2113T1_genome.fasta
samtools faidx Pc2113T1_genome.fasta

IFS=,
# running BWA mem then piping the output to samtools to sort bam file then convert that to a cram file
tail -n +2 $SAMPFILE | while read ISOLATE
do
  # BWA
  bwa mem -t 64 Pc2113T1 \
  -M \
  -R $(echo "@RG\tID:"$ISOLATE"\tSM:"$ISOLATE"\tLB:"$ISOLATE"\tPL:ILLUMINA") \
  $FILTERREADDIR/$ISOLATE.F.1.fq.gz \
  $FILTERREADDIR/$ISOLATE.F.2.fq.gz \
  | samtools sort -@ 32 -O bam -l 0 -T /tmp - | \
  samtools view -@ 32 -T Pc2113T1_genome.fasta -C -o $SORTEDCRAM/$ISOLATE.sort.cram -
done
