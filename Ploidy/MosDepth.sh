#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=12-00:00:00
##SBATCH --output=my.stdout
#SBATCH --job-name="Mosdepth"
#SBATCH -p batch

module load mosdepth/0.3.3
module load samtools/1.14


bam_list=Final_Bams.list
mean_dps_out=stats/mean-dp_PcPopAll.mosdepth.txt
mapQ=20
touch $mean_dps_out
while read bam; do
  sm_name=$(basename --suffix ".RG.mkdup.bam" $bam)
  mosdepth --no-per-base -t 32 --mapq 20 ${sm_name}_mapQ$mapQ $bam
  mean_cov_mosdp=$(tail -n +2 stats/${sm_name}_mapQ${mapQ}.mosdepth.summary.txt | awk '{w = w + $2; e = e + $4 * $2;} END {print e/w}')
  printf "${sm_name}\t${mean_cov_mosdp}\n" >> $mean_dps_out
done < $bam_list
