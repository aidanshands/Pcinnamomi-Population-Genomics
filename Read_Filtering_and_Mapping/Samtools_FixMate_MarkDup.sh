#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH--mem=64G
#SBATCH -c 32
#SBATCH --time=12-00:00:00
##SBATCH --output=my.stdout
#SBATCH --mail-type=ALL
#SBATCH --job-name="Samtools"
#SBATCH -p batch

# Load modules
module load samtools/1.14

cd $SLURM_SUBMIT_DIR

# Setting variables
SAMPFILE=Isolates.csv
IFS=,


# Setting genome-specific variables
INCRAMDIR=Sorted_CRAMs
OUTCRAMDIR=Name_Sorted_CRAMs
FIXMATEDIR=Fixmate_CRAMs
MKDUPSTATSDIR=Sorted_mkdup_CRAM_Stats
FINALBAMDIR=Final_Bams
FLAGSTATSDIR=Bam_Flagstats
# running samtools sort to sort by name (required by samtools fixmate) and then running samtools fixmate
tail -n +2 $SAMPFILE | while read ISOLATE
do
  samtools sort --reference Pc2113T1_genome.fasta -T $ISOLATE -@ 32 -n $INCRAMDIR/$ISOLATE.sort.cram -o $OUTCRAMDIR/$ISOLATE.namesort.cram
  samtools fixmate -m -@ 32 -O cram $OUTCRAMDIR/$ISOLATE.namesort.cram $FIXMATEDIR/$ISOLATE.namesort.fixmate.cram
done

# Sorting the fixmate CRAM by coordinate (required by markdup), then running samtools markdup, then generating flagstats
tail -n +2 $SAMPFILE | while read ISOLATE
do
  samtools sort --reference Pc2113T1_genome.fasta -O bam -@ 32 -T $ISOLATE $FIXMATEDIR/$ISOLATE.namesort.fixmate.cram | \
  samtools markdup -T $ISOLATE -f $MKDUPSTATSDIR/$ISOLATE.mkdup.stats.txt - $FINALBAMDIR/$ISOLATE.mkdup.bam
  samtools flagstat $FINALBAMDIR/$ISOLATE.mkdup.bam > $FLAGSTATSDIR/$ISOLATE.mkdup.flagstat.txt
done

# Here we are looking at our flagstats outputs and creating a nice summary file with all samples and the % Mapped
# Kindly shared by Nick Carelson (https://github.com/Neato-Nick)
tail -n +2 $SAMPFILE | while read ISOLATE
do
  mapped_pct=$(grep "mapped" $FLAGSTATSDIR/$ISOLATE.mkdup.flagstat.txt | awk -F "[(|%]" '{print $2}' | head -n 1)
  printf "$ISOLATE\t${mapped_pct}\n" >> $FLAGSTATSDIR/SUM.mkdup.flagstat.txt
done
