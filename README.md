# _Phytophthora cinnamomi_ Population Genomics

Please note that this repository is under development, thank you.

The workflow and materials used for the _P. cinnamomi_ population genomic study in Shands _et al._ (2024). Code utilizing array jobs were adopted from my dissertation committee member Dr. Jason Stajich ([https://github.com/hyphaltip/PopGenomics_Aureobasidium_pullulans/tree/ma](https://github.com/hyphaltip/PopGenomics_Aureobasidium_pullulans/tree/main)), thank you Jason. 


## Read Filtering, Mapping & BAM Processing

See Read Filtering and Mapping in Methods section in Shands _et al._ (2024). 

**Read Filtering**

Script: Read_Filtering.sh. 

**Mapping**

Script: Mapping.sh

**BAM Processing & Mapping Statistics**

Script: Samtools_FixMate_MarkDup.sh

## Variant Calling, Filtering & Annotation

See Variant Calling, Filtering & Annotation in Methods section in Shands _et al._ (2024). 

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


## Population Genomic Analyses

See Population Genomic Analyses in Methods section in Shands _et al._ (2024). 

**Population Genomics in R**

Scripts: PopGen_Analysis.R & Calculate_Ia_Simulation.v2.R

**FastStructure**

Scripts: FastStructure_SP_k1-10.sh & FastStructure_LP_k2-5.sh

**Mosdepth**

Script: MosDepth.sh

**nQuire**

Script: nQuire.sh

**vcfR Allele Balance**

Script: Allele_Balance.R

## Caluclating EC50

See _In Vitro_ Sensitivity to Potassium Phosphite in Methods section in Shands _et al._ (2024). 

Script: Calculate_EC50.py



