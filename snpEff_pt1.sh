#!/bin/bash -l

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
# working dir: /bigdata/manosalvalab/pmanosal/Pc_PopGen_2022/variant_calling_diploid/Effector_Selection/SNPeff
VCFDIR=Sample_VCFs
VCF=Pc_only136.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.recode.vcf.gz
SAMPFILE=Isolate_Names.csv
VCFOUTDIR=snpEff_VCFs
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
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN
do
  java -Xmx32g -jar /rhome/ashan007/snpEff/snpEff.jar -v -stats \
  $STRAIN.html Pc2113 \
  $VCFDIR/$STRAIN.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.vcf > \
  $VCFOUTDIR/$STRAIN.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.snpEff.vcf
done
