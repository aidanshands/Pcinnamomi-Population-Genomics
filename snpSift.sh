#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH --time=20-00:00:00
##SBATCH --output=my.stdout
#SBATCH --job-name="SNPeff"
#SBATCH -p batch
#SBATCH --array 1-136


# Load modules
module load bcftools
module load samtools/1.14

cd $SLURM_SUBMIT_DIR
VCFDIR=snpEff_VCFs
SAMPFILE=Isolate_Names_Pc136.csv
VCFOUTDIR=snpSift_VCFs
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

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN
do
  java -Xmx32g -jar /rhome/ashan007/snpEff/SnpSift.jar \
  filter 'isVariant( GEN[0] )' \
  -f $VCFDIR/$STRAIN.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.snpEff.vcf > \
  $VCFOUTDIR/$STRAIN.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.snpSift.vcf
done
