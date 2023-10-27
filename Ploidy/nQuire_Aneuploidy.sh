#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=5-00:00:00
##SBATCH --output=my.stdout
#SBATCH --mail-user=ashan007@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="nQuire"
#SBATCH -p batch

# wd = '/bigdata/manosalvalab/pmanosal/Pc_PopGen_2022/nQuire_Ploidy/Cov_Adj_nQuire'
cd $SLURM_SUBMIT_DIR

PC2113COV=Pc2113_mapQ20.mosdepth.ContigCov.csv
ZITRCOV=ZITR3_3_mapQ20.mosdepth.ContigCov.csv
ZITRBINDIR=ZITR3-3_Bins
ZITRDENOISE_BINDIR=ZITR3-3_Denoised_Bins
PC2113BINDIR=Pc2113_Bins
PC2113DENOISE_BINDIR=Pc2113_Denoised_Bins
PC2113LRD_Out_FILE="Pc2113_lrdmodel_Output.txt"
PC2113LRD_DN_Out_FILE="Pc2113_lrdmodel_Denoised_Output.txt"
ZITRLRD_Out_FILE="ZITR3-3_lrdmodel_Output.txt"
ZITRLRD_DN_Out_FILE="ZITR3-3_lrdmodel_Denoised_Output.txt"
BEDDIR=/bigdata/manosalvalab/pmanosal/Pc_PopGen_2022/nQuire_Ploidy/nQuire_Aneuploidy/Pc2113_Contig_BEDs
#-------------------------------------------------------------------------------
# Pc2113
#-------------------------------------------------------------------------------
# nQuire Create
IFS=,
tail -n +2 $PC2113COV | while read CONTIG COV
do
  ~/nQuire/nQuire create \
  -b /bigdata/manosalvalab/pmanosal/Pc_PopGen_2022/mapping/Final_BAMs_RGs/Pc2113.RG.mkdup.bam \
  -o $PC2113BINDIR/Pc2113 \
  -q 20 -c $COV -r $BEDDIR/$CONTIG.bed
done;

# nQuire lrdmodel
for file in $PC2113BINDIR/*.bin; do
    ~/nQuire/nQuire lrdmodel -t 32 "$file" >> "$PC2113LRD_Out_FILE"
done;

# Denoise Bins
tail -n +2 $PC2113COV | while read CONTIG COV
do
  ~/nQuire/nQuire denoise $PC2113BINDIR/Pc2113-$CONTIG.bin -o $PC2113DENOISE_BINDIR/$CONTIG
done;

# nQuire lrdmodel on Denooised Bins
for file in $PC2113DENOISE_BINDIR/*.bin; do
    ~/nQuire/nQuire lrdmodel -t 32 "$file" >> "$PC2113LRD_DN_Out_FILE"
done;
#-------------------------------------------------------------------------------
# ZITR3-3
#-------------------------------------------------------------------------------
tail -n +2 $ZITRCOV | while read CONTIG COV
do
  ~/nQuire/nQuire create \
  -b /bigdata/manosalvalab/pmanosal/Pc_PopGen_2022/mapping/Final_BAMs_RGs/ZITR3_3.RG.mkdup.bam \
  -o $ZITRBINDIR/ZITR3_3 \
  -q 20 -c $COV -r $BEDDIR/$CONTIG.bed
done;

# nQuire lrdmodel
for file in $ZITRBINDIR/*.bin; do
    ~/nQuire/nQuire lrdmodel -t 32 "$file" >> "$ZITRLRD_Out_FILE"
done;

# Denoise Bins
tail -n +2 $ZITRCOV | while read CONTIG COV
do
  ~/nQuire/nQuire denoise $ZITRBINDIR/ZITR3_3-$CONTIG.bin -o $ZITRDENOISE_BINDIR/$CONTIG
done;

IFS=,
# nQuire lrdmodel on Denooised Bins
for file in $ZITRDENOISE_BINDIR/*.bin; do
    ~/nQuire/nQuire lrdmodel -t 32 "$file" >> "$ZITRLRD_DN_Out_FILE"
done;
