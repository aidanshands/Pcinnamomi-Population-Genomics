# 7/27/23
# Aidan Shands

library(adegenet)
library(poppr)
library(vcfR)
library(ape)
library(RColorBrewer)
library(ape)
library(dartR)
library(mmod)
library(ggplot2)
library(agricolae)
library(ggthemes)
#setwd('/bigdata/manosalvalab/pmanosal/Pc_PopGen_2022/variant_calling_diploid/PopGen_Stats_R/R_Analysis_on_HPCC_PAPER')
setwd('/Users/manosalvalab/Desktop/Aidan/Pc_PopGen_2022/R_Analysis_on_HPCC_PAPER')
# helpful resource : https://github.com/grunwaldlab/GBS-Pcinnamomi/blob/master/By_group.Rmd

#-------------------------------------------------------------------------------
# Importing & Setting global variables
#-------------------------------------------------------------------------------
# Annotations 
annot = read.csv('annotations/Pc_Annot.PAPER.csv', header = T)
# VCF 
vcf_file_OG = 'Pc_only136.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.recode.vcf.gz'
vcf_OG<- read.vcfR(vcf_file_OG, verbose = FALSE)
vcf_OG

# creating Genlight and GenInd objects 
vcf_GL_OG<- vcfR2genlight(vcf_OG)
vcf_GL_OG

vcf_GL_OG$ind.names
annot$ID
ploidy(vcf_GL_OG) <- 2
vcf_GL_OG$pop <- as.factor(annot$location_Simple)
vcf_GL_OG$other <- as.list(annot$ID)

set.seed(100)

DatasetCols <- c("Asia/Oceania" = "#1B9E77", "Mexico" = "#E7298A", 
                  "CA-South" = "#7570B3", "CA-North" = "#D95F02",
                 "Clonal" = "azure3", "Mostly-clonal" = "azure3" , 
                 "Semi-clonal" = "azure3", "Sexual" = "azure3" )

#-------------------------------------------------------------------------------
# Defining Functions
#-------------------------------------------------------------------------------
# Simulation and Ia function
SimulateGL_CaluclateIA = function(Genlight, n_reps, n_SNPs, Genlight_Name){
  # Simulate Genlights
  Clone = glSim(nInd(Genlight),
                n.snp.nonstruc = floor(0.1*nLoc(Genlight)),
                n.snp.struc=ceiling(0.9*nLoc(Genlight)), ploidy=2, LD = T)
  Most_clone <- glSim(nInd(Genlight),
                      n.snp.nonstruc = ceiling(nLoc(Genlight)/3),
                      n.snp.struc= 2*nLoc(Genlight)/3, ploidy=2, LD=T)
  Semi_clone <- glSim(nInd(Genlight),
                      n.snp.nonstruc = 0.5*nLoc(Genlight),
                      n.snp.struc= 0.5*nLoc(Genlight), ploidy=2, LD=T)
  Sexual <- glSim(n.ind = nInd(Genlight), 
                        n.snp.nonstruc = ceiling(0.9*nLoc(Genlight)), 
                        n.snp.struc = floor(0.1*nLoc(Genlight)), ploidy=2, LD=T)
  # Calculate Ia using samp.ia()
  ia.clone <- samp.ia(Clone, quiet = T, reps = n_reps, n.snp = n_SNPs)
  ia.most <- samp.ia(Most_clone, quiet = T, reps = n_reps, n.snp = n_SNPs)
  ia.semi <- samp.ia(CA_South_semi_clone, quiet = T,reps = n_reps, n.snp = n_SNPs)
  ia.sex <- samp.ia(Sexual, quiet = T, reps = n_reps, n.snp = n_SNPs)
  ia.Pc <- samp.ia(Genlight,  reps = n_reps, quiet = T, n.snp = n_SNPs)
  # cast results to df
  df1 <- data.frame(ia.Pc, rep(Genlight_Name, length(ia.Pc)))
  df2 <- data.frame(ia.sex, rep("Sexual", length(ia.sex)))
  df3 <- data.frame(ia.clone, rep("Clonal", length(ia.clone)))
  df4 <- data.frame(ia.semi, rep("Semi-clonal", length(ia.semi)))
  df5 <- data.frame(ia.most, rep("Mostly-clonal", length(ia.most)))
  # rename df columns
  colnames(df1) <- c("ia","dataset")
  colnames(df2) <- c("ia","dataset")
  colnames(df3) <- c("ia","dataset")
  colnames(df4) <- c("ia","dataset")
  colnames(df5) <- c("ia","dataset")
  ia.total <- rbind(df3, df5, df4, df2, df1)
  return(ia.total)
}

# Shapiro Wilks test 
QC_Stats = function(data, Pc_Data){
  Sexual = subset(data, dataset == "Sexual")
  Clonal = subset(data, dataset == "Clonal")
  Semiclonal = subset(data, dataset == "Semi-clonal")
  Mostlyclonal = subset(data, dataset == "Mostly-clonal")
  Pc = subset(data, dataset == Pc_Data)
  frames <- list(as.data.frame(Sexual), as.data.frame(Clonal), 
                 as.data.frame(Semiclonal), as.data.frame(Mostlyclonal), 
                 as.data.frame(Pc))
  normality <- list()
  for (i in 1:length(frames)){
    normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
  }
  return(normality)
}

# Performs One-way ANOVA, Tukey's HSD or Fisher's LSD
Stats = function(data, variable, value, test){
  dataname = deparse(substitute(data)) # get name of data we are working on 
  print(paste("Stats for ", dataname, sep=""))
  print(paste("One-Way ANOVA for: ", dataname, sep=""))
  # generate a formula (value~variable)
  formula <- reformulate(variable, response = value)
  AOV = aov(lm(formula, data = data)) # run one-way anova
  outnameAOV = summary(AOV)# generate summary
  print(outnameAOV) # print Summary
  if (test == "Fisher"){
    print(paste("Fisher's LSD Groups for: ", dataname, sep=""))
    TestResults = (LSD.test(AOV, factor(variable), group=TRUE))# run Fisher's LSD
  } else if (test == "Tukey"){
    print(paste("Tukey HSD Groups for: ", dataname, sep=""))
    TestResults = (HSD.test(AOV, factor(variable), group=TRUE))# run Tukey HSD
  }
  # generate a df with Tukey HSD or Fisher's results
  results <- data.frame(
    Isolate = rownames(TestResults$groups),
    Mean = TestResults$groups[, value],
    Groups = TestResults$groups[, "groups"],
    stringsAsFactors = FALSE
  )
  return(results)
}

#-------------------------------------------------------------------------------
#  Raw Subsets 
#-------------------------------------------------------------------------------
Pops <- seppop(vcf_GL_OG)

CA_South = as.snpclone(Pops$'CA-South')
CA_North = as.snpclone(Pops$'CA-North')
Asia = as.snpclone(Pops$'Asia/Oceania')
Mexico = as.snpclone(Pops$'Mexico')

#-------------------------------------------------------------------------------
#  All Data Uncorrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
All_IA = SimulateGL_CaluclateIA(vcf_GL_OG, n_reps = 100, n_SNPs = 10000, "All Isolates")
write_csv(All_IA, "All_IA.csv")
All_Normality = QC_Stats(All_IA, "All Isolates")
All_Normality
# Analysis of variance
All_Results = Stats(All_IA, "dataset", "ia", "Tukey")
colnames(All_Results) <- c("dataset","Mean", "Groups")
All_Results
write_csv(All_Results, "All_IA_Results.csv")
# Visualization
All_BP = ggplot(All_IA, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

All_BP = All_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
All_BP


#-------------------------------------------------------------------------------
#  CA-South Uncorrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
CA_South_IA = SimulateGL_CaluclateIA(CA_South, n_reps = 100, n_SNPs = 10000, "CA-South")
write_csv(CA_South_IA, "CA_South_IA.csv")
CAS_Normality = QC_Stats(CA_South_IA, "CA-South")
CAS_Normality
# Analysis of variance
CAS_Results = Stats(CA_South_IA, "dataset", "ia", "Tukey")
colnames(CAS_Results) <- c("dataset","Mean", "Groups")
CAS_Results
write_csv(CAS_Results, "CA_South_IA_Results.csv")

# Visualization
CAS_BP = ggplot(CA_South_IA, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

CAS_BP = CAS_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
CAS_BP

#-------------------------------------------------------------------------------
#  CA-North Uncorrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
CA_North_IA = SimulateGL_CaluclateIA(CA_North, n_reps = 100, n_SNPs = 10000, "CA-North")
write_csv(CA_North_IA, "CA_North_IA.csv")
CAN_Normality = QC_Stats(CA_North_IA, "CA-North")
CAN_Normality
# Analysis of variance
CAN_Results = Stats(CA_North_IA, "dataset", "ia", "Tukey")
colnames(CAN_Results) <- c("dataset","Mean", "Groups")
CAN_Results
write_csv(CAN_Results, "CA_North_IA_Results.csv")
# Visualization

CAN_BP = ggplot(CA_North_IA, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

CAN_BP = CAN_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
CAN_BP

#-------------------------------------------------------------------------------
#  Mexico Uncorrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
Mexico_IA = SimulateGL_CaluclateIA(MX, n_reps = 100, n_SNPs = 10000, "Mexico")
write_csv(Mexico_IA, "Mexico_IA.csv")
MX_Normality = QC_Stats(Mexico_IA, "Mexico")
MX_Normality
# Analysis of variance
Mexico_Results = Stats(Mexico_IA, "dataset", "ia", "Tukey")
colnames(Mexico_Results) <- c("dataset","Mean", "Groups")
Mexico_Results
write_csv(Mexico_Results, "Mexico_IA_Results.csv")
# Visualization
MX_BP = ggplot(Mexico_IA, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

MX_BP = MX_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
MX_BP

#-------------------------------------------------------------------------------
#  Asia/Oceania Uncorrected 
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
Asia_IA = SimulateGL_CaluclateIA(Asia, n_reps = 100, n_SNPs = 10000, "Asia/Oceania")
write_csv(Asia_IA, "Asia_Oceania_IA.csv")
Asia_Normality = QC_Stats(Asia_IA, "Asia/Oceania")
Asia_Normality
# Analysis of variance
Asia_Results = Stats(Asia_IA, "dataset", "ia", "Tukey")
colnames(Asia_Results) <- c("dataset","Mean", "Groups")
Asia_Results
write_csv(Asia_Results, "Asia_Oceania_IA_Results.csv")
# Visualization
AO_BP = ggplot(Asia_IA, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

AO_BP = AO_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
AO_BP


#-------------------------------------------------------------------------------
#  Visualizing All Uncorrected 
#-------------------------------------------------------------------------------
# Since the simulated values (Clonal, Semi-clonal, & Sexual) for each population 
# do not change much, I'm going to plot the simulated values of the overall 
# genlight (All_IA_Results). I'm going to combine the results from each pop and
# plot them together. 
Combined_Uncorrected = data.frame()
Combined_Uncorrected = rbind(
  Asia_IA[Asia_IA$dataset == "Asia/Oceania", ],
  Mexico_IA[Mexico_IA$dataset == "Mexico", ],
  CA_South_IA[CA_South_IA$dataset == "CA-South", ],
  CA_North_IA[CA_North_IA$dataset == "CA-North", ],
  All_IA
)

# Custom Order
custom_order = c("Clonal", "Mostly-clonal", "Semi-clonal", "Sexual",
                  "All Isolates", "Asia/Oceania", "CA-North", "CA-South", "Mexico")
plot_order = c("Clonal", "Semi-clonal", "Sexual",
              "All Isolates", "Asia/Oceania", "CA-North", "CA-South", "Mexico")
Combined_Uncorrected$dataset <- factor(Combined_Uncorrected$dataset, levels = custom_order)
write_csv(Combined_Uncorrected, "Combined_Uncorrected_IA.csv")

Combined_Uncorrected = read.csv("IA_Simulation_Outputs/Combined_Uncorrected_IA.csv")
# remove Mostly-clonal
Combined_Uncorrected = subset(Combined_Uncorrected, dataset != "Mostly-clonal")
# Visualization
Combined_BP = ggplot(Combined_Uncorrected, aes(x = factor(dataset, levels = plot_order), y=ia, fill=as.factor(dataset)))+ 
  geom_boxplot(lwd = 0.3) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

Combined_BP = Combined_BP + theme_few() + 
  labs(x="", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "none")
Combined_BP

#ggsave(Combined_BP,filename = "IA_BP.FINAL.pdf", width = 8, height = 6, units = "in", dpi = 300)
ggsave("IA_BP.FINAL.v2.pdf", Combined_BP, width=4.28, height=3.78, units="in")

# Remove the All Isolates sample 
Combined_Uncorrected_2 <- subset(Combined_Uncorrected, dataset != "All Isolates")

# Visualization Only Pops 
Combined_BP2 = ggplot(Combined_Uncorrected_2, aes(x = factor(dataset, levels = plot_order), y=ia, fill=as.factor(dataset))) + 
  geom_boxplot(lwd = 0.1, outlier.size = 0.3) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

Combined_BP2 = Combined_BP2 + theme_few() + 
  labs(x="", y = "Index of association")+
  theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 0.55, .1), limits = c(0, 0.55))
Combined_BP2

# ggsave(Combined_BP2,filename = "IA_BP_OnlyPops.FINAL.pdf", width = 8, height = 6, units = "in", dpi = 300)
ggsave("IA_BP_OnlyPops.FINAL.v2.pdf", Combined_BP2, width=4.28, height=3.78, units="in")

#-------------------------------------------------------------------------------
#  Clone Correction & Subsets 
#-------------------------------------------------------------------------------
# convert to snpclone
SC = as.snpclone(vcf_GL_OG)
pop(SC) <- as.factor(annot$location_Simple) # set pop by location
strata(SC) = annot # set strata 

# for MLG filter of 0.01
mlg.filter(SC, algorithm = "average_neighbor", distance = "bitwise.dist") <- 0.01
SC

# Clone correct 
SC_CC = clonecorrect(SC, strata= ~location_Simple)
SC_CC

#-------------------------------------------------------------------------------
#  Clone-Corrected Subsets 
#-------------------------------------------------------------------------------
Pops_CC <- seppop(SC_CC)

CA_South_CC = as.snpclone(Pops_CC$'CA-South')
CA_North_CC = as.snpclone(Pops_CC$'CA-North')
Asia_CC = as.snpclone(Pops_CC$'Asia/Oceania')
Mexico_CC = as.snpclone(Pops_CC$'Mexico')

#-------------------------------------------------------------------------------
#  All Data Clone-corrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
All_IA_CC = SimulateGL_CaluclateIA(SC_CC, n_reps = 100, n_SNPs = 10000, "All Isolates")
write_csv(All_IA_CC, "All_IA_CC.csv")
All_Normality_CC = QC_Stats(All_IA_CC, "All Isolates")
All_Normality_CC
# Analysis of variance
All_Results_CC = Stats(All_IA_CC, "dataset", "ia", "Tukey")
colnames(All_Results_CC) <- c("dataset","Mean", "Groups")
All_Results_CC
write_csv(All_Results_CC, "All_IA_CC_Results.csv")
# Visualization
All_BP_CC = ggplot(All_IA_CC, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

All_BP_CC = All_BP_CC + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
All_BP_CC


#-------------------------------------------------------------------------------
#  CA-South Clone-corrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
CA_South_IA_CC= SimulateGL_CaluclateIA(CA_South_CC, n_reps = 100, n_SNPs = 10000, "CA-South")
write_csv(CA_South_IA_CC, "CA_South_IA_CC.csv")
CAS_Normality_CC = QC_Stats(CA_South_IA_CC, "CA-South")
CAS_Normality_CC
# Analysis of variance
CAS_Results_CC = Stats(CA_South_IA_CC, "dataset", "ia", "Tukey")
colnames(CAS_Results_CC) <- c("dataset","Mean", "Groups")
CAS_Results_CC
write_csv(CAS_Results_CC, "CA_South_IA_CC_Results.csv")

# Visualization
CASCC_BP = ggplot(CA_South_IA_CC, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

CASCC_BP = CASCC_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
CASCC_BP

#-------------------------------------------------------------------------------
#  CA-North Clone-corrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
CA_North_IA_CC = SimulateGL_CaluclateIA(CA_North_CC, n_reps = 100, n_SNPs = 10000, "CA-North")
write_csv(CA_North_IA_CC, "CA_North_IA_CC.csv")
CAN_Normality_CC = QC_Stats(CA_North_IA_CC, "CA-North")
CAN_Normality_CC
# Analysis of variance
CAN_Results_CC = Stats(CA_North_IA_CC, "dataset", "ia", "Tukey")
colnames(CAN_Results_CC) <- c("dataset","Mean", "Groups")
CAN_Results_CC
write_csv(CAN_Results_CC, "CA_North_IA_CC_Results.csv")
# Visualization

CANCC_BP = ggplot(CA_North_IA_CC, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

CANCC_BP = CANCC_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
CANCC_BP


#-------------------------------------------------------------------------------
#  Mexico Clone-corrected
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
Mexico_IA_CC = SimulateGL_CaluclateIA(Mexico_CC, n_reps = 100, n_SNPs = 10000, "Mexico")
write_csv(Mexico_IA_CC, "Mexico_IA_CC.csv")
MX_Normality_CC = QC_Stats(Mexico_IA_CC, "Mexico")
MX_Normality_CC
# Analysis of variance
Mexico_Results_CC = Stats(Mexico_IA_CC, "dataset", "ia", "Tukey")
colnames(Mexico_Results_CC) <- c("dataset","Mean", "Groups")
Mexico_Results_CC
write_csv(Mexico_Results_CC, "Mexico_IA_CC_Results.csv")
# Visualization
MXCC_BP = ggplot(Mexico_IA_CC, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

MXCC_BP = MXCC_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
MXCC_BP

#-------------------------------------------------------------------------------
#  Asia/Oceania Clone-corrected 
#-------------------------------------------------------------------------------
# Genlight Simulation & Ia Calculation
Asia_IA_CC = SimulateGL_CaluclateIA(Asia_CC, n_reps = 100, n_SNPs = 10000, "Asia/Oceania")
write_csv(Asia_IA_CC, "Asia_Oceania_IA_CC.csv")
Asia_Normality_CC = QC_Stats(Asia_IA_CC, "Asia/Oceania")
Asia_Normality_CC
# Analysis of variance
Asia_Results_CC = Stats(Asia_IA_CC, "dataset", "ia", "Tukey")
colnames(Asia_Results_CC) <- c("dataset","Mean", "Groups")
Asia_Results_CC
write_csv(Asia_Results_CC, "Asia_Oceania_IA_CC_Results.csv")
# Visualization
AOCC_BP = ggplot(Asia_IA_CC, aes(x = dataset, y=ia, fill=as.factor(dataset))) + 
  geom_boxplot() +
  #geom_text(data = CAS_Results, 
  #aes(label = Groups, y = max(CAS_ia.total$ia) + 0.02),
  #vjust = -0.5) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

AOCC_BP = AOCC_BP + theme_classic() + 
  labs(x="Dataset", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))
AOCC_BP






#-------------------------------------------------------------------------------
#  Visualizing All Clone-corrected  
#-------------------------------------------------------------------------------
# Since the simulated values (Clonal, Semi-clonal, & Sexual) for each population 
# do not change much, I'm going to plot the simulated values of the overall 
# genlight (All_IA_Results). I'm going to combine the results from each pop and
# plot them together. 
Combined_Corrected = data.frame()
Combined_Corrected = rbind(
  Asia_IA_CC[Asia_IA_CC$dataset == "Asia/Oceania", ],
  Mexico_IA_CC[Mexico_IA_CC$dataset == "Mexico", ],
  CA_South_IA_CC[CA_South_IA_CC$dataset == "CA-South", ],
  CA_North_IA_CC[CA_North_IA_CC$dataset == "CA-North", ],
  All_IA_CC
)

# Custom Order
Combined_Corrected$dataset <- factor(Combined_Corrected$dataset, levels = custom_order)
write_csv(Combined_Corrected, "Combined_Corrected_IA.csv")

Combined_Corrected = read.csv("IA_Simulation_Outputs/Combined_Corrected_IA.csv")
# remove Mostly-clonal
Combined_Corrected = subset(Combined_Corrected, dataset != "Mostly-clonal")
# Visualization All samples 
Combined_BP_CC1 = ggplot(Combined_Corrected, aes(x = factor(dataset, levels = plot_order), y=ia, fill=as.factor(dataset))) + 
  geom_boxplot(lwd = 0.3) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

Combined_BP_CC1 = Combined_BP_CC1 + theme_bw() + 
  labs(x="", y = "Index of Association")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "none")
Combined_BP_CC1

ggsave(Combined_BP_CC1,filename = "IA_CC_BP.FINAL.pdf", 
       width = 8, height = 6, units = "in", dpi = 300)

# Remove the All Isolates sample 
Combined_Corrected_2 <- subset(Combined_Corrected, dataset != "All Isolates")

# Visualization Only Pops 
Combined_BP_CC2 = ggplot(Combined_Corrected_2, aes(x = factor(dataset, levels = plot_order), y=ia, fill=as.factor(dataset))) + 
  geom_boxplot(lwd = 0.1, outlier.size = 0.3) +
  scale_fill_manual(values=DatasetCols, name="Dataset") 

Combined_BP_CC2 = Combined_BP_CC2 + theme_few() + 
  labs(x="", y = "Index of association")+
  theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10),
        legend.position = "none")+
  scale_y_continuous(breaks = seq(0, 0.55, .1), limits = c(0, 0.55))
Combined_BP_CC2

#ggsave(Combined_BP_CC2,filename = "IA_CC_BP_OnlyPops.FINAL.pdf", width = 8, height = 6, units = "in", dpi = 300)
ggsave("IA_CC_BP_OnlyPops.FINAL.v2.pdf", Combined_BP_CC2, width=4.28, height=3.78, units="in")

#-------------------------------------------------------------------------------
#  Stats Excluding Mostly-Clonal 
#-------------------------------------------------------------------------------

# Shapiro Wilks test 
QC_Stats2 = function(data, Pc_Data){
  Sexual = subset(data, dataset == "Sexual")
  Clonal = subset(data, dataset == "Clonal")
  Semiclonal = subset(data, dataset == "Semi-clonal")
  Pc = subset(data, dataset == Pc_Data)
  frames <- list(as.data.frame(Sexual), as.data.frame(Clonal), 
                 as.data.frame(Semiclonal), as.data.frame(Pc))
  normality <- list()
  for (i in 1:length(frames)){
    normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
  }
  return(normality)
}

#-------------
#  Uncorrected 
#-------------
# All pops 
All_IA = read.csv("IA_Simulation_Outputs/All_IA.csv")
All_IA = subset(All_IA, dataset != "Mostly-clonal")
All_Normality = QC_Stats2(All_IA, "All Isolates")
All_Normality
# Analysis of variance
All_Results = Stats(All_IA, "dataset", "ia", "Tukey")
colnames(All_Results) <- c("dataset","Mean", "Groups")
All_Results
write.csv(All_Results, "IA_Simulation_Outputs/All_IA_Results_NoMC.csv")


# CA-South
CA_South_IA = read.csv("IA_Simulation_Outputs/CA_South_IA.csv")
CA_South_IA = subset(CA_South_IA, dataset != "Mostly-clonal")
CA_South_Normality = QC_Stats2(CA_South_IA, "CA-South")
CA_South_Normality
# Analysis of variance
CA_South_Results = Stats(CA_South_IA, "dataset", "ia", "Tukey")
colnames(CA_South_Results) <- c("dataset","Mean", "Groups")
CA_South_Results
write.csv(CA_South_Results, "IA_Simulation_Outputs/CA_South_IA_Results_NoMC.csv")

# CA-North
CA_North_IA = read.csv("IA_Simulation_Outputs/CA_North_IA.csv")
CA_North_IA = subset(CA_North_IA, dataset != "Mostly-clonal")
CA_North_Normality = QC_Stats2(CA_North_IA, "CA-North")
CA_North_Normality
# Analysis of variance
CA_North_Results = Stats(CA_North_IA, "dataset", "ia", "Tukey")
colnames(CA_North_Results) <- c("dataset","Mean", "Groups")
CA_North_Results
write.csv(CA_North_Results, "IA_Simulation_Outputs/CA_North_IA_Results_NoMC.csv")

# Mexico
Mexico_IA = read.csv("IA_Simulation_Outputs/Mexico_IA.csv")
Mexico_IA = subset(Mexico_IA, dataset != "Mostly-clonal")
Mexico_Normality = QC_Stats2(Mexico_IA, "Mexico")
Mexico_Normality
# Analysis of variance
Mexico_Results = Stats(Mexico_IA, "dataset", "ia", "Tukey")
colnames(Mexico_Results) <- c("dataset","Mean", "Groups")
Mexico_Results
write.csv(Mexico_Results, "IA_Simulation_Outputs/Mexico_IA_Results_NoMC.csv")

# Asia/Oceania
Asia_IA = read.csv("IA_Simulation_Outputs/Asia_Oceania_IA.csv")
Asia_IA = subset(Asia_IA, dataset != "Mostly-clonal")
Asia_Normality = QC_Stats2(Asia_IA, "Asia/Oceania")
Asia_Normality
# Analysis of variance
Asia_Results = Stats(Asia_IA, "dataset", "ia", "Tukey")
colnames(Asia_Results) <- c("dataset","Mean", "Groups")
Asia_Results
write.csv(Asia_Results, "IA_Simulation_Outputs/Asia_Oceania_IA_Results_NoMC.csv")

#-----------------
#  Clone-corrected 
#-----------------

# All pops 
All_IA_CC = read.csv("IA_Simulation_Outputs/All_IA_CC.csv")
All_IA_CC = subset(All_IA_CC, dataset != "Mostly-clonal")
All_Normality_CC = QC_Stats2(All_IA_CC, "All Isolates")
All_Normality_CC
# Analysis of variance
All_Results_CC = Stats(All_IA_CC, "dataset", "ia", "Tukey")
colnames(All_Results_CC) <- c("dataset","Mean", "Groups")
All_Results
write.csv(All_Results_CC, "IA_Simulation_Outputs/All_IA_Results_CC_NoMC.csv")


# CA-South
CA_South_IA_CC = read.csv("IA_Simulation_Outputs/CA_South_IA_CC.csv")
CA_South_IA_CC = subset(CA_South_IA_CC, dataset != "Mostly-clonal")
CA_South_Normality_CC = QC_Stats2(CA_South_IA_CC, "CA-South")
CA_South_Normality_CC
# Analysis of variance
CA_South_Results_CC = Stats(CA_South_IA_CC, "dataset", "ia", "Tukey")
colnames(CA_South_Results_CC) <- c("dataset","Mean", "Groups")
CA_South_Results_CC
write.csv(CA_South_Results_CC, "IA_Simulation_Outputs/CA_South_IA_Results_CC_NoMC.csv")

# CA-North
CA_North_IA_CC = read.csv("IA_Simulation_Outputs/CA_North_IA_CC.csv")
CA_North_IA_CC = subset(CA_North_IA_CC, dataset != "Mostly-clonal")
CA_North_Normality_CC = QC_Stats2(CA_North_IA_CC, "CA-North")
CA_North_Normality_CC
# Analysis of variance
CA_North_Results_CC = Stats(CA_North_IA_CC, "dataset", "ia", "Tukey")
colnames(CA_North_Results_CC) <- c("dataset","Mean", "Groups")
CA_North_Results_CC
write.csv(CA_North_Results_CC, "IA_Simulation_Outputs/CA_North_IA_CC_Results_NoMC.csv")

# Mexico
Mexico_IA_CC = read.csv("IA_Simulation_Outputs/Mexico_IA_CC.csv")
Mexico_IA_CC = subset(Mexico_IA_CC, dataset != "Mostly-clonal")
Mexico_Normality_CC = QC_Stats2(Mexico_IA_CC, "Mexico")
Mexico_Normality_CC
# Analysis of variance
Mexico_Results_CC = Stats(Mexico_IA_CC, "dataset", "ia", "Tukey")
colnames(Mexico_Results_CC) <- c("dataset","Mean", "Groups")
Mexico_Results_CC
write.csv(Mexico_Results_CC, "IA_Simulation_Outputs/Mexico_IA_Results_CC_NoMC.csv")

# Asia/Oceania
Asia_IA_CC = read.csv("IA_Simulation_Outputs/Asia_Oceania_IA_CC.csv")
Asia_IA_CC = subset(Asia_IA_CC, dataset != "Mostly-clonal")
Asia_Normality_CC = QC_Stats2(Asia_IA_CC, "Asia/Oceania")
Asia_Normality_CC
# Analysis of variance
Asia_Results_CC = Stats(Asia_IA_CC, "dataset", "ia", "Tukey")
colnames(Asia_Results_CC) <- c("dataset","Mean", "Groups")
Asia_Results_CC
write.csv(Asia_Results_CC, "IA_Simulation_Outputs/Asia_Oceania_IA_Results_CC_NoMC.csv")


