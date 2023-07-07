#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH --time=20-00:00:00
##SBATCH --output=my.stdout
#SBATCH --job-name="GATK Diploid"
#SBATCH -p batch
#SBATCH --array 1-735
# Load modules

module load picard
module load java/17.0.2
module load gatk/4.2.5.0
module load bcftools
module load samtools/1.14

MEM=32g

if [ -f config.txt ]; then
    source config.txt
fi


out_name="final_vcfs/genotyped_contig${SLURM_ARRAY_TASK_ID}.vcf.gz"

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
 N=$1
fi

i=$(expr $SLURM_ARRAY_TASK_ID - 1)

FAI=( `cut -f 1 "${REFGENOME}.fai" `)
read contig <<< "${FAI[$i]}"

if [ ! $N ]; then
 echo "need to provide a number by --array slurm or on the cmdline"
 exit
fi

hostname
date
gatk --java-options -Xmx${MEM} GenotypeGVCFs \
-R $REFGENOME -L $contig -V Combined136.g.vcf.gz -O $out_name
