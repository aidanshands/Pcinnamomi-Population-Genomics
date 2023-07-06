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

## Variant Calling, Filtering & Annotation
Single-nucleotide polymorphisms (SNPs) were called from each isolate using Genome Analysis Toolkit (GATK; version 4.2.5.0.) (DePristo et al., 2011) in two steps. First, variants were called HaplotypeCaller using the -ERC GVCF option with the ploidy level set to 2. Following, the variant files were merged using GATK’s CombineGVCFs function and the merged gvcf file was subjected to variant calling using GATK’s GenotypeGVCFs function. The resulting VCF files were then merged using PICARD v.2.26.11 MergeVcfs function. The merged vcf was subjected to light filtering using GATK’s SelectVariants where only the biallelic SNPs were retained, unused alternatives removed, and non-variants were excluded.

**GATK HaplotypeCaller**

Script: GATK_HaplotypeCaller.sh

**Picard MergeVcfs**

Script: Combine_gVCF.sh

**GATK GenotypeGVCF**

Script: GATK_GenotypeGVCF.sh

**GATK SelectVariants**

Script: GATK_SelectVariants.sh

**VCFtools Filtering**

Script: VCFtools_Filter.sh

