#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=12-00:00:00
##SBATCH --output=my.stdout
#SBATCH --job-name="SelectVariants"
#SBATCH -p batch

module load picard
module load java/17.0.2
module load gatk/4.2.5.0
module load bcftools
module load samtools/1.14

MEM=32g

if [ -f config.txt ]; then
    source config.txt
fi

gatk --java-options -Xmx${MEM} SelectVariants \
-R $REFGENOME \
-V Combined142.vcf.gz \
--restrict-alleles-to BIALLELIC \
--remove-unused-alternates true \
--select-type-to-include SNP \
--exclude-non-variants true \
--sample-name Pc_Isolate_Names.list \
-O Pc_only136.p2.f1SV.vcf.gz
