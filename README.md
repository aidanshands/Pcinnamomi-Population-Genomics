# _Phytophthora cinnamomi_ Population Genomics

## Read Trimming & Filtering

Read filtering was performed using fastq-mcf v1.05 (Aronesty, 2011). Script: Read_Filtering.sh. 

``` bash
SAMPFILE=Isolates.csv
RAWREADDIR=Raw_Reads
STATDIR=Filtered_Stats
FILTERREADDIR=Filtered_Reads

IFS=,
# Filering the reads & Getting the stats
tail -n +2 $SAMPFILE | while read ISOLATE
do
  fastq-mcf \
  -q 20 -o $FILTERREADDIR/$ISOLATE.F.1.fq.gz -o $FILTERREADDIR/$ISOLATE.F.2.fq.gz \
  Adapters.fasta \
  $RAWREADDIR/$ISOLATE.1.fq.gz \
  $RAWREADDIR/$ISOLATE.2.fq.gz
  # filtered stats
  fastq-stats $FILTERREADDIR/$ISOLATE.F.1.fq.gz > $STATDIR/$ISOLATE.stats.1.txt
  fastq-stats $FILTERREADDIR/$ISOLATE.F.2.fq.gz > $STATDIR/$ISOLATE.stats.2.txt
done
```

## Read Mapping & BAM Processing

Reads were mapped to _P. cinnamomi_ isolate 2113 (Shands et al., 2023) (GenBank:XXXXX) using BWA-mem v. 0.7.17 (Li 2013) with the following settings: -M, -R, -t 64. The respective BWA-mem output was piped into samtools (Li et al., 2009) and the respective BAM files were sorted by position for samtools flagstats. Next, the BAM files were sorted by name for samtools fixmate and the respective output sorted by position and marked for duplicates using samtools mkdup. These respective BAM files served as the input for variant calling. 

**Mapping**

Script: Mapping.sh

``` bash
SAMPFILE=Isolates.csv
FILTERREADDIR=Filtered_Reads
SORTEDCRAM=Sorted_CRAMs
# Indexing genome
bwa index -p Pc2113 Pc2113T1_genome.fasta
samtools faidx Pc2113T1_genome.fasta

IFS=,
# running BWA mem then piping the output to samtools to sort bam file then convert that to a cram file
tail -n +2 $SAMPFILE | while read ISOLATE
do
  # BWA
  bwa mem -t 64 Pc2113T1_genome.fasta \
  -M \
  -R $(echo "@RG\tID:"$ISOLATE"\tSM:"$ISOLATE"\tLB:"$ISOLATE"\tPL:ILLUMINA") \
  $FILTERREADDIR/$ISOLATE.Filtered.1.fq.gz \
  $FILTERREADDIR/$ISOLATE.Filtered.2.fq.gz \
  | samtools sort -@ 32 -O bam -l 0 -T /tmp - | \
  samtools view -@ 32 -T Pc2113T1_genome.fasta -C -o $SORTEDCRAM/$ISOLATE.sort.cram -
done
```

**BAM Processing & Mapping Statistics**

Script: Samtools_FixMate_MarkDup.sh

``` bash
# Setting genome-specific variables
SAMPFILE=Isolates.csv
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
```

