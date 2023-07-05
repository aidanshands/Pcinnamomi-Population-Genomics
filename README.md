# _Phytophthora cinnamomi_ Population Genomics

The workflow and materials used for the _P. cinnamomi_ population genomic study in Shands _et al._ (2023). Code utilizing array jobs were adopted from my dissertation committee member Dr. Jason Stajich ([https://github.com/hyphaltip/PopGenomics_Aureobasidium_pullulans/tree/ma](https://github.com/hyphaltip/PopGenomics_Aureobasidium_pullulans/tree/main)), thank you Jason. 

## Read Trimming & Filtering

Read filtering was performed using fastq-mcf v1.05 (Aronesty, 2011). Script: Read_Filtering.sh. 


## Read Mapping & BAM Processing

Reads were mapped to _P. cinnamomi_ isolate 2113 (Shands et al., 2023) (GenBank:XXXXX) using BWA-mem v. 0.7.17 (Li 2013). The respective BWA-mem output was piped into samtools (Li et al., 2009) and the respective BAM files were sorted by position for samtools flagstats. Next, the BAM files were sorted by name for samtools fixmate and the respective output sorted by position and marked for duplicates using samtools mkdup. These respective BAM files served as the input for variant calling. 

**Mapping**

Script: Mapping.sh

**BAM Processing & Mapping Statistics**

Script: Samtools_FixMate_MarkDup.sh

## Haplotype Calling via GATK 

Script: GATK_HaplotypeCaller.sh

