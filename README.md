# _Phytophthora cinnamomi_ Population Genomics

## Read Trimming & Filtering

Read filtering was performed using fastq-mcf v1.05 (Aronesty, 2011). 

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
