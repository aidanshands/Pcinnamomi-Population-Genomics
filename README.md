# _Phytophthora cinnamomi_ Population Genomics

The workflow and materials used for the _P. cinnamomi_ population genomic study in Shands _et al._ (2023). Code utilizing array jobs were adopted from my dissertation committee member Dr. Jason Stajich ([https://github.com/hyphaltip/PopGenomics_Aureobasidium_pullulans/tree/ma](https://github.com/hyphaltip/PopGenomics_Aureobasidium_pullulans/tree/main)), thank you Jason. 


## Read Filtering, Mapping & BAM Processing

See Read Filtering and Mapping in Methods section in Shands _et al._ (2023). 

**Read Filtering**

Script: Read_Filtering.sh. 

**Mapping**

Script: Mapping.sh

**BAM Processing & Mapping Statistics**

Script: Samtools_FixMate_MarkDup.sh

## Variant Calling, Filtering & Annotation

See Variant Calling, Filtering & Annotation in Methods section in Shands _et al._ (2023). 

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

**SNPeff Annotation**

Scripts: snpEff_pt1.sh, snpSift.sh & snpEff_pt2.sh

## Population Genomic Analyses

See Population Genomic Analyses in Methods section in Shands _et al._ (2023). 

**Population Genomics in R**

Script: PopGen_Analysis.R

**FastStructure**

Scripts: FastStructure_SP_k1-10.sh & FastStructure_LP_k2-5.sh

**Mosdepth**

Script: MosDepth.sh

**nQuire**

Script: nQuire.sh

**vcfR Allele Balance**

Script: Allele_Balance.R

## Caluclating EC50

See _In Vitro_ Sensitivity to Potassium Phosphite in Methods section in Shands _et al._ (2023). 

Script: Calculate_EC50.py



