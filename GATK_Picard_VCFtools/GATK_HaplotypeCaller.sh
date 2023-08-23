#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH --time=20-00:00:00
##SBATCH --output=my.stdout
#SBATCH --job-name="GATK Diploid"
#SBATCH --array=1-136
#SBATCH -p batch
# Load modules
module load picard
module load java/17.0.2
module load gatk/4.2.5.0
module load bcftools
module load samtools/1.14

MEM=32g
SAMPFILE=Isolate_Names.csv

if [ -f config.txt ]; then
    source config.txt
fi

DICT=$(echo $REFGENOME | sed 's/fasta$/dict/')

if [ ! -f $DICT ]; then
        picard CreateSequenceDictionary R=$GENOMEIDX O=$DICT
fi
mkdir -p $GVCFFOLDER
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "need to provide a number by --array slurm or on the cmdline"
 exit
fi

hostname
date
IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SAMPID
do
  # BEGIN THIS PART IS PROJECT SPECIFIC LIKELY
  # END THIS PART IS PROJECT SPECIFIC LIKELY
  echo "STRAIN is $STRAIN"
  GVCF=$GVCFFOLDER/$STRAIN.g.vcf.gz
  ALNFILE=$ALNFOLDER/$STRAIN.RG.mkdup.bam
  if [ -s $GVCF ]; then
    echo "Skipping $STRAIN - Already called $STRAIN.g.vcf.gz"
    exit
  fi
  if [[ ! -f $GVCF || $ALNFILE -nt $GVCF ]]; then
      time gatk --java-options -Xmx${MEM} HaplotypeCaller \
          -ERC GVCF --sample-ploidy 2 \
          -I $ALNFILE -R $REFGENOME \
          -O $GVCF --native-pair-hmm-threads $CPU
 fi
done
