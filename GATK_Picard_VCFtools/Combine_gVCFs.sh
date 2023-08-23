#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=12-00:00:00
##SBATCH --output=my.stdout
#SBATCH --mail-type=ALL
#SBATCH --job-name="CombineVCFs"
#SBATCH -p batch

module load picard
module load java/17.0.2
module load gatk/4.2.5.0
module load bcftools
module load samtools/1.14

picard MergeVcfs \
I=Diploid_vcf.list \
O=Combined142.vcf.gz
