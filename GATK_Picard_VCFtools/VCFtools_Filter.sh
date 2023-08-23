#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=12-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="filtering"
#SBATCH -p batch

module load picard
module load java/17.0.2
module load gatk/4.2.5.0
module load bcftools
module load samtools/1.14
module load miniconda3
module load vcftools/0.1.16-18

#-------------------------------------------------------------------------------
# Global settings
# setting dirs
wd=variant_calling
# setting variables and names
vcf=Pc_only136.p2.f1SV.vcf.gz
fil_wd=Light_Filtering
#-------------------------------------------------------------------------------
# Identify Singletons
singletons=$fil_wd/singletons
singletons_out=$fil_wd/singletons.singletons
singletons_Pos=$fil_wd/Singletons_Pos.singletons
# create singletons file and extract position
CMD="vcftools --gzvcf $vcf --singletons --out $singletons"
echo $CMD
#eval $CMD
awk 'NF{NF-=3};3' < $singletons_out > $singletons_Pos
#-------------------------------------------------------------------------------
# Light filtering part 1
max_missing_vars_1=1 # No missing data
min_dp=4
min_Q=20
remove_sms_out="$fil_wd/""$(basename --suffix ".vcf.gz" $vcf)"".vcft_dp${min_dp}_Q${min_Q}_NoSingleton_noMD"
CMD="vcftools --gzvcf $vcf --out $remove_sms_out --recode --recode-INFO-all \
--max-missing $max_missing_vars_1 --minDP $min_dp --minQ $min_Q --exclude-positions $singletons_Pos"
echo $CMD
#eval $CMD
#-------------------------------------------------------------------------------
# Light filtering part 2
vcf2=Light_Filtering/Pc_only136.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD.recode.vcf
bed=Pc2113T1_genome.repeats.v2.bed
F2_out="$fil_wd/Pc_only136.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats"
CMD="vcftools --vcf $vcf2 --out $F2_out --recode --recode-INFO-all --exclude-bed $bed"
echo $CMD
#eval $CMD
