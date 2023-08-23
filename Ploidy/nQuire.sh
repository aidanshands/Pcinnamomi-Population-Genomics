#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=5-00:00:00
##SBATCH --output=my.stdout
#SBATCH --mail-type=ALL
#SBATCH --job-name="nQuire"
#SBATCH -p batch

cd $SLURM_SUBMIT_DIR

SAMPFILE=PcOnly.mosdepth_MeanCov_Rounded.csv
BINDIR=Bins
DENOISE_BINDIR=Denoised_Bins
LRD_Out_FILE="lrdmodel_Output.txt"
LRD_DN_Out_FILE="lrdmodel_Denoised_Output.txt"

# nQuire Create
IFS=,
tail -n +2 $SAMPFILE | while read ISOLATE COV
do
  ~/nQuire/nQuire create \
  -b /bigdata/manosalvalab/pmanosal/Pc_PopGen_2022/mapping/Final_BAMs_RGs/$ISOLATE.RG.mkdup.bam \
  -o $BINDIR/$ISOLATE \
  -q 20 -c $COV
done;

# nQuire lrdmodel
for file in $BINDIR/*.bin; do
    ~/nQuire/nQuire lrdmodel -t 32 "$file" >> "$LRD_Out_FILE"
done

# Denoise Bins
tail -n +2 $SAMPFILE | while read ISOLATE COV
do
  ~/nQuire/nQuire denoise $BINDIR/$ISOLATE.bin -o $DENOISE_BINDIR/$ISOLATE
done;

IFS=,
# nQuire lrdmodel on Denooised Bins
for file in $DENOISE_BINDIR/*.bin; do
    ~/nQuire/nQuire lrdmodel -t 32 "$file" >> "$LRD_DN_Out_FILE"
done
