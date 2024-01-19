# Aidan Shands

### ----------------------------------------------------------------------------
# Load packages
### ----------------------------------------------------------------------------
library(ggplot2)
library(agricolae)
library(RColorBrewer)
library(car)
library(cowplot)
library(ggpubr)
library(multcompView)
library(fdth)
library(MASS)
library(janitor)
library(colorspace)
library(unikn)
library(grid)
library(gridExtra)
library(gtable)
library(FSA)
library(dplyr)
library(tidyverse)
library(rcompanion)
### ----------------------------------------------------------------------------
# Importing data, & Setting factors 
### ----------------------------------------------------------------------------
setwd("/Users/manosalvalab/Desktop/Aidan/Pc_PopGen_2022/PopGen_Scripts_Github/GrowthRate_Phenotyping")

OGT_Data <- read.csv("GR_Data.csv", header = T)
# Setting factors 
OGT_Data$Isolate <- factor(OGT_Data$Isolate)
OGT_Data$RGPD_22 <- as.numeric(OGT_Data$RGPD_22)
OGT_Data$RGPD_25 <- as.numeric(OGT_Data$RGPD_25)
OGT_Data$RGPD_28 <- as.numeric(OGT_Data$RGPD_28)
OGT_Data$County <- factor(OGT_Data$County)
OGT_Data$State <- factor(OGT_Data$State)
OGT_Data$Country <- factor(OGT_Data$Country)
OGT_Data$Replicate <- factor(OGT_Data$Replicate)
OGT_Data$Growing_Region <- factor(OGT_Data$Growing_Region)

# Setting colors 
Cols <- c("South" = "#7570B3", "North" = "#D95F02")

### ----------------------------------------------------------------------------
# Defining Functions
### ----------------------------------------------------------------------------
# Perform the Levene's test, Bartlett's test, Shapiro-wilks test and QQplot
QC_Stats = function(data, variable){
  dataname = deparse(substitute(data))
  col = data[[variable]]
  print(paste("Levene's Test by Replicate for: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Replicate), data = data))
  print(paste("Levene's Test by Isolate for: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Isolate), data = data))
  print(paste("Bartlett's Test by Replicate for: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Replicate), data = data)$p.value)
  print(paste("Bartlett's Test by Isolate for: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Isolate), data = data)$p.value)
  print(paste("Shapiro Wilks test for: ", dataname, sep=""))
  print(shapiro.test(col))
  print(paste("QQplot: ", dataname, sep=""))
  qqPlot(col)
}
# Performs One-way ANOVA, Tukey's HSD or Fisher's LSD (Parametric)
P_Stats = function(data, variable, value, test){
  dataname = deparse(substitute(data)) # get name of data we are working on 
  print(paste("Stats for ", dataname, sep=""))
  print(paste("One-Way ANOVA for: ", dataname, sep=""))
  # generate a formula (value~variable)
  formula = reformulate(variable, response = value)
  AOV = aov(formula, data = data) # run one-way anova
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

# Performs Kruskal-Wallis test and the Dunn's test (with letter groups assigned) (Non-parametric)
NP_Stats = function(data, variable, value){
  dataname = deparse(substitute(data)) # get name of data we are working on 
  print(paste("Stats for ", dataname, sep=""))
  print(paste("Kruskal-Wallis for: ", dataname, sep=""))
  # generate a formula (value~variable)
  formula =reformulate(variable, response = value)
  KW = kruskal.test(formula, data = data) # run one-way anova
  print(KW) # print Summary
  print(paste("Dunn's Test Groups for: ", dataname, sep=""))
  DTresult = dunnTest(formula, data = data, method="holm")$res
  DunnResults = cldList(P.adj ~ Comparison, data = DTresult, threshold = 0.05, remove.zero = FALSE)
  DunnResults = DunnResults %>%
    select(-MonoLetter) %>%
    rename(Isolate = Group, Dunn_Group = Letter)
  return(DunnResults)
}

### ----------------------------------------------------------------------------
# QC Stats 
### ----------------------------------------------------------------------------
# 22˚C  
QC_Stats(OGT_Data, "RGPD_22")
# 25˚C
QC_Stats(OGT_Data, "RGPD_25")
# 28˚C
QC_Stats(OGT_Data, "RGPD_28")

### ----------------------------------------------------------------------------
# One-way ANOVA & Tukey's HSD Test 
### ----------------------------------------------------------------------------
# ANOVA by mean radial growth  by isolate
# 22˚C 
ONA_TT_22 = P_Stats(OGT_Data, "Isolate", "RGPD_22", "Tukey")
# 25˚C 
ONA_TT_25 = P_Stats(OGT_Data, "Isolate", "RGPD_25", "Tukey")
# 28˚C 
ONA_TT_28 = P_Stats(OGT_Data, "Isolate", "RGPD_28", "Tukey")

### ----------------------------------------------------------------------------
# One-way ANOVA & Fisher's LSD Test 
### ----------------------------------------------------------------------------
# ANOVA by mean radial growth by isolate
# 22˚C 
ONA_FLSD_22 = P_Stats(OGT_Data, "Isolate", "RGPD_22", "Fisher")
# 25˚C 
ONA_FLSD_25 = P_Stats(OGT_Data, "Isolate", "RGPD_25", "Fisher")
# 28˚C 
ONA_FLSD_28 = P_Stats(OGT_Data, "Isolate", "RGPD_28", "Fisher")

### ----------------------------------------------------------------------------
# Kruskal-Wallis & Dunn's test 
### ----------------------------------------------------------------------------
# 22˚C 
KW_DUNN_22 = NP_Stats(OGT_Data, "Isolate", "RGPD_22")
# 25˚C 
KW_DUNN_25 = NP_Stats(OGT_Data, "Isolate", "RGPD_25")
# 28˚C 
KW_DUNN_28 = NP_Stats(OGT_Data, "Isolate", "RGPD_28")

### ----------------------------------------------------------------------------
# Summarizing data (Avg + SD)
### ----------------------------------------------------------------------------
SUM22 <- as.data.frame(OGT_Data %>% group_by(Isolate, County, State, Country, Growing_Region) %>% 
  dplyr::summarize(Mean = mean(RGPD_22, na.rm=TRUE), SD = sd(RGPD_22, na.rm=TRUE)))

SUM25 <- as.data.frame(OGT_Data %>% group_by(Isolate, County, State, Country, Growing_Region) %>% 
  dplyr::summarize(Mean = mean(RGPD_25, na.rm=TRUE), SD = sd(RGPD_25, na.rm=TRUE)))

SUM28 <- as.data.frame(OGT_Data %>% group_by(Isolate, County, State, Country, Growing_Region) %>% 
  dplyr::summarize(Mean = mean(RGPD_28, na.rm=TRUE), SD = sd(RGPD_28, na.rm=TRUE)))

# rounding the decimals to 2 decimals
SUM22 = SUM22 %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUM25 = SUM25 %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUM28 = SUM28 %>% mutate(across(where(is.numeric), ~ round(., 2)))
# renaming columns
colnames(SUM22) <- c('Isolate', 'County', 'State', 'Country', "Growing_Region", 'RGPD_22', 'RGPD_22_SD')
colnames(SUM25) <- c('Isolate', 'County', 'State', 'Country', "Growing_Region", 'RGPD_25', 'RGPD_25_SD')
colnames(SUM28) <- c('Isolate', 'County', 'State', 'Country', "Growing_Region", 'RGPD_28', 'RGPD_28_SD')
# exporting summaries
#dir.create("Summary_Outputs_Final")
#write.csv(SUM22, "Summary_Outputs_Final/All_Data_22SUM.csv")
#write.csv(SUM25, "Summary_Outputs_Final/All_Data_25SUM.csv")
#write.csv(SUM28, "Summary_Outputs_Final/All_Data_28SUM.csv")

# Creating summary df for plotting
SUM22_j = subset(SUM22, select = c(Isolate, County, State, Growing_Region, Country, RGPD_22, RGPD_22_SD))
SUM25_j = subset(SUM25, select = c(Isolate, RGPD_25, RGPD_25_SD))
SUM28_j = subset(SUM28, select = c(Isolate, RGPD_28, RGPD_28_SD))
# Tukey HSD Results
TT22_j = subset(ONA_TT_22, select = c(Isolate, Groups))
colnames(TT22_j) <- c('Isolate', 'TT_22_Groups')
TT25_j = subset(ONA_TT_25, select = c(Isolate, Groups))
colnames(TT25_j) <- c('Isolate', 'TT_25_Groups')
TT28_j = subset(ONA_TT_28, select = c(Isolate, Groups))
colnames(TT28_j) <- c('Isolate', 'TT_28_Groups')
# Fisher's LSD Results
FLSD22_j = subset(ONA_FLSD_22, select = c(Isolate, Groups))
colnames(FLSD22_j) <- c('Isolate', 'FLSD_22_Groups')
FLSD25_j = subset(ONA_FLSD_25, select = c(Isolate, Groups))
colnames(FLSD25_j) <- c('Isolate', 'FLSD_25_Groups')
FLSD28_j = subset(ONA_FLSD_28, select = c(Isolate, Groups))
colnames(FLSD28_j) <- c('Isolate', 'FLSD_28_Groups')
# Dunn's Test Results 
DUNN22_j = subset(KW_DUNN_22, select = c(Isolate, Dunn_Group))
colnames(DUNN22_j) <- c('Isolate', 'Dunn_22_Groups')
DUNN25_j = subset(KW_DUNN_25, select = c(Isolate, Dunn_Group))
colnames(DUNN25_j) <- c('Isolate', 'Dunn_25_Groups')
DUNN28_j = subset(KW_DUNN_28, select = c(Isolate, Dunn_Group))
colnames(DUNN28_j) <- c('Isolate', 'Dunn_28_Groups')

Sum_List = list(SUM22_j, SUM25_j, SUM28_j,TT22_j, TT25_j, TT28_j, 
                FLSD22_j, FLSD25_j, FLSD28_j, DUNN22_j, DUNN25_j, DUNN28_j)
# This data is used in Table S4
Comb_Summary = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List) 
head(Comb_Summary)
#write.csv(Comb_Summary,"Summary_Outputs_Final/All_Data_Summarized.csv")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### ----------------------------------------------------------------------------
# Subsetting CA only and looking at North vs South 
### ----------------------------------------------------------------------------
CA_Isolates = subset(OGT_Data, State == "CA")
table(CA_Isolates$Growing_Region)

### ----------------------------------------------------------------------------
# QC Stats 
### ----------------------------------------------------------------------------
# 22˚C  
QC_Stats(CA_Isolates, "RGPD_22")
# 25˚C
QC_Stats(CA_Isolates, "RGPD_25")
# 28˚C
QC_Stats(CA_Isolates, "RGPD_28")

#### ----------------------------------------------------------------------------
# One-way ANOVA & Tukey's HSD Test by growing region
### ----------------------------------------------------------------------------
# 22˚C 
CA_ONA_TT_22 = P_Stats(CA_Isolates, "Growing_Region", "RGPD_22", "Tukey")
print(CA_ONA_TT_22)
# 25˚C 
CA_ONA_TT_25 = P_Stats(CA_Isolates, "Growing_Region", "RGPD_25", "Tukey")
print(CA_ONA_TT_25)
# 28˚C 
CA_ONA_TT_28 = P_Stats(CA_Isolates, "Growing_Region", "RGPD_28", "Tukey")
print(CA_ONA_TT_28)

### ----------------------------------------------------------------------------
# One-way ANOVA & Fisher's LSD Test by growing region
### ----------------------------------------------------------------------------
# 22˚C 
CA_ONA_FLSD_22 = P_Stats(CA_Isolates, "Growing_Region", "RGPD_22", "Fisher")
print(CA_ONA_FLSD_22)
# 25˚C 
CA_ONA_FLSD_25 = P_Stats(CA_Isolates, "Growing_Region", "RGPD_25", "Fisher")
print(CA_ONA_FLSD_25)
# 28˚C 
CA_ONA_FLSD_28 = P_Stats(CA_Isolates, "Growing_Region", "RGPD_28", "Fisher")
print(CA_ONA_FLSD_28)

### ----------------------------------------------------------------------------
# Kruskal-Wallis test by growing region
### ----------------------------------------------------------------------------
# 22˚C 
print(kruskal.test(RGPD_22 ~Growing_Region, data = CA_Isolates)) # Significant
# 25˚C 
print(kruskal.test(RGPD_25 ~Growing_Region, data = CA_Isolates)) # Significant 
# 28˚C 
print(kruskal.test(RGPD_28 ~Growing_Region, data = CA_Isolates)) # Not Significant

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### ----------------------------------------------------------------------------
# Performing Temperature Comparisons
### ----------------------------------------------------------------------------

OGT_Data2 <- read.csv("GR_Temperature_Data.csv", header = T)
head(OGT_Data2)

# Setting factors 
OGT_Data2$Isolate <- factor(OGT_Data2$Isolate)
OGT_Data2$RGPD <- as.numeric(OGT_Data2$RGPD)
OGT_Data2$County <- factor(OGT_Data2$County)
OGT_Data2$State <- factor(OGT_Data2$State)
OGT_Data2$Country <- factor(OGT_Data2$Country)
OGT_Data2$Replicate <- factor(OGT_Data2$Replicate)
OGT_Data2$Growing_Region <- factor(OGT_Data2$Growing_Region)
OGT_Data2$Temperature <- factor(OGT_Data2$Temperature)

### ----------------------------------------------------------------------------
# QC Stats for Temperature
### ----------------------------------------------------------------------------
QC_Stats(OGT_Data2, "RGPD")

### ----------------------------------------------------------------------------
# One-way ANOVA & Tukey's HSD Test 
### ----------------------------------------------------------------------------
ONA_TT_Temp = P_Stats(OGT_Data2, "Temperature", "RGPD", "Tukey")
print(ONA_TT_Temp)
### ----------------------------------------------------------------------------
# One-way ANOVA & Fisher's LSD Test 
### ----------------------------------------------------------------------------
ONA_FLSD_Temp = P_Stats(OGT_Data2, "Temperature", "RGPD", "Fisher")
print(ONA_FLSD_Temp)
### ----------------------------------------------------------------------------
# Kruskal-Wallis & Dunn's test 
### ----------------------------------------------------------------------------
KW_DUNN_Temp = NP_Stats(OGT_Data2, "Temperature", "RGPD")
print(KW_DUNN_Temp)
### ----------------------------------------------------------------------------
# Subsetting for North & South 
### ----------------------------------------------------------------------------
CA_Isolates_So = subset(OGT_Data2, Growing_Region == "South")
CA_Isolates_No = subset(OGT_Data2, Growing_Region == "North")
CA_Isolates_All = subset(OGT_Data2, State == "CA")

### ----------------------------------------------------------------------------
# One-Way ANOVA/Kruskal-Wallis  looking at RGR by Temp in CA 
### ----------------------------------------------------------------------------
ONA_TT_Temp_CA = P_Stats(CA_Isolates_All, "Temperature", "RGPD", "Tukey")
print(ONA_TT_Temp_CA)

ONA_FLSD_Temp_CA = P_Stats(CA_Isolates_All, "Temperature", "RGPD", "Fisher")
print(ONA_FLSD_Temp_CA)

KW_DUNN_Temp_CA = NP_Stats(CA_Isolates_All, "Temperature", "RGPD")
print(KW_DUNN_Temp_CA)
### ----------------------------------------------------------------------------
# One-Way ANOVA/Kruskal-Wallis looking at RGR by Temp in CA-South
### ----------------------------------------------------------------------------
ONA_TT_Temp_CAS = P_Stats(CA_Isolates_So, "Temperature", "RGPD", "Tukey")
print(ONA_TT_Temp_CAS)

ONA_FLSD_Temp_CAS = P_Stats(CA_Isolates_So, "Temperature", "RGPD", "Fisher")
print(ONA_FLSD_Temp_CAS)

KW_DUNN_Temp_CAS = NP_Stats(CA_Isolates_So, "Temperature", "RGPD")
print(KW_DUNN_Temp_CAS)

### ----------------------------------------------------------------------------
# One-Way ANOVA/Kruskal-Wallis  looking at RGR by Temp in CA-North
### ----------------------------------------------------------------------------
ONA_TT_Temp_CAN = P_Stats(CA_Isolates_No, "Temperature", "RGPD", "Tukey")
print(ONA_TT_Temp_CAN)

ONA_FLSD_Temp_CAN = P_Stats(CA_Isolates_No, "Temperature", "RGPD", "Fisher")
print(ONA_FLSD_Temp_CAN)

KW_DUNN_Temp_CAN = NP_Stats(CA_Isolates_No, "Temperature", "RGPD")
print(KW_DUNN_Temp_CAN)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### ----------------------------------------------------------------------------
# Ratio of Growth at 22 vs 28 
### ----------------------------------------------------------------------------
OGT_Data3 <- read.csv("22_vs_28.csv", header = T)

# Setting factors 
OGT_Data3$Isolate <- factor(OGT_Data3$Isolate)
OGT_Data3$County <- factor(OGT_Data3$County)
OGT_Data3$State <- factor(OGT_Data3$State)
OGT_Data3$Country <- factor(OGT_Data3$Country)
OGT_Data3$Growing_Region <- factor(OGT_Data3$Growing_Region)
OGT_Data3$Percent_Growth <- as.numeric(OGT_Data3$Percent_Growth)
OGT_Data3$Percent_Growth_Affected_by_28 <- as.numeric(OGT_Data3$Percent_Growth_Affected_by_28)

# qqPlots 
qqPlot(OGT_Data3$Percent_Growth_Affected_by_28)
# Shapiro-Wilks test 
shapiro.test(OGT_Data3$Percent_Growth_Affected_by_28) 
# Bartlett's Test
print(bartlett.test(Percent_Growth_Affected_by_28 ~ Growing_Region, data = OGT_Data3)$p.value)

### ----------------------------------------------------------------------------
# One-Way ANOVA looking at Percent_Growth_Affected_by_28 by Growing_Region CA only
### ----------------------------------------------------------------------------
# Subset CA 
CA_Isolates2 = subset(OGT_Data3, State == "CA")
# qqPlots 
qqPlot(CA_Isolates2$Percent_Growth_Affected_by_28)
# Shapiro-Wilks test 
shapiro.test(CA_Isolates2$Percent_Growth_Affected_by_28) 
# Bartlett's Test
print(bartlett.test(Percent_Growth_Affected_by_28 ~ Growing_Region, data = CA_Isolates2)$p.value)
# (sig. diff. p<0.05)
ONA_TT_Temp_Ratio_CA = P_Stats(CA_Isolates2, "Growing_Region", "Percent_Growth_Affected_by_28", "Tukey")
print(ONA_TT_Temp_Ratio_CA)

ONA_FLSD_Temp_Ratio_CA = P_Stats(CA_Isolates2, "Growing_Region", "Percent_Growth_Affected_by_28", "Fisher")
print(ONA_FLSD_Temp_Ratio_CA)
### ----------------------------------------------------------------------------
# Generating Figures for Paper 
### ----------------------------------------------------------------------------

#----------
# Figure 2A
#----------
BP4 = ggplot(OGT_Data2, aes(x=Temperature, y=RGPD)) + geom_boxplot(outlier.size = 0.75)
BP4= BP4 + theme_classic() + 
  labs(x="Temperature (\u00B0C)", y = "Mean Radial Growth Rate (mm/day)")+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title=element_text(size=8,face="bold"))+
  guides(fill = guide_legend(title = "Growing Region"))
BP4
ggsave("OGTBoxplot_All_Paper.v2.pdf",BP4, width=4.28, height=3.78, units="in")

# Letter groupings are 22˚C-a, 25˚C-b, 28˚C-c and despite non-normality the 
# groupings are all the same between ANOVA + Tukey's HSD, ANOVA + Fisher's LSD,
# and Kruskal-Wallis + Dunn's test

#----------
# Figure 3A
#----------
BP8 = ggplot(CA_Isolates_All, aes(x=Temperature, y=RGPD, fill= Growing_Region)) + geom_boxplot(outlier.size = 0.75)
BP8= BP8 + theme_classic() + scale_fill_manual(values=Cols) + 
  labs(x="Temperature (\u00B0C)", y = "Mean Radial Growth Rate (mm/day)") + 
  theme(axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title=element_text(size=10,face="bold"),
      legend.position = "none")
BP8
#ggsave("OGTBoxplot_All_CA_Paper.v2.pdf",BP8, width=4.28, height=3.78, units="in")

# Dunn groupings of CA-North at each temp (letters above the box)
# 22˚C - a, 25˚C-b, 28˚C-c
# KW_DUNN_Temp_CA (correct method)

# Dunn groupings of CA-South at each temp (letters above the box)
# 22˚C - a, 25˚C-b, 28˚C-c
# see KW_DUNN_Temp_CAS (correct method)

# Groupings of North vs South @ each Temp (within the two boxes)
# We performed Kruskal-Wallis test for this and it does not yeild groups
# but there were significant differences at 22˚C, 25˚C but not at 28˚C. 
# I placed a or b within the two boxes for each region at each temp for the
# figure. 


#----------
# Figure 3B
#----------
PTemp<- ggplot(CA_Isolates2, aes(x=Isolate, y=Percent_Growth_Affected_by_28,  fill= Growing_Region)) + 
  geom_bar(stat="identity", 
           position=position_dodge(width = 0.9), 
           aes(reorder(Isolate,Percent_Growth_Affected_by_28),Percent_Growth_Affected_by_28), width = 0.5)

PTemp = PTemp+labs(x="Isolate", y = "Percent") +
  theme_classic() + scale_fill_manual(values=Cols) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10,face="bold"),
        legend.position = "none") + coord_flip()
PTemp
#ggsave("Percent_Temp_CA_Paper.v2.pdf",PTemp, width=6, height=8, units="in")

### ----------------------------------------------------------------------------
# Session info
### ----------------------------------------------------------------------------
sessionInfo()
### ----------------------------------------------------------------------------
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] rcompanion_2.4.30  lubridate_1.9.2    forcats_1.0.0      stringr_1.5.0      purrr_1.0.1        readr_2.1.4       
# [7] tidyr_1.3.0        tibble_3.2.1       tidyverse_2.0.0    dplyr_1.1.1        dunn.test_1.3.5    FSA_0.9.4         
# [13] gtable_0.3.3       gridExtra_2.3      unikn_0.8.0        colorspace_2.1-0   janitor_2.2.0      MASS_7.3-54       
# [19] fdth_1.2-6         multcompView_0.1-8 ggpubr_0.6.0       cowplot_1.1.1      car_3.1-2          carData_3.0-5     
# [25] plyr_1.8.8         RColorBrewer_1.1-3 agricolae_1.3-5    ggplot2_3.4.2     

# loaded via a namespace (and not attached):
#   [1] TH.data_1.1-2       ggsignif_0.6.4      ellipsis_0.3.2      class_7.3-19        modeltools_0.2-23  
# [6] snakecase_0.11.0    gld_2.6.6           rstudioapi_0.14     proxy_0.4-27        farver_2.1.1       
# [11] fansi_1.0.4         mvtnorm_1.1-3       coin_1.4-2          codetools_0.2-18    splines_4.1.2      
# [16] rootSolve_1.8.2.3   libcoin_1.0-9       broom_1.0.4         cluster_2.1.2       png_0.1-8          
# [21] shiny_1.7.4         compiler_4.1.2      httr_1.4.5          backports_1.4.1     Matrix_1.3-4       
# [26] fastmap_1.1.1       cli_3.6.1           later_1.3.0         htmltools_0.5.5     tools_4.1.2        
# [31] ggmap_3.0.2         glue_1.6.2          lmom_2.9            Rcpp_1.0.10         cellranger_1.1.0   
# [36] vctrs_0.6.1         nlme_3.1-153        lmtest_0.9-40       timechange_0.2.0    mime_0.12          
# [41] miniUI_0.1.1.1      lifecycle_1.0.3     rstatix_0.7.2       zoo_1.8-11          scales_1.2.1       
# [46] hms_1.1.3           promises_1.2.0.1    parallel_4.1.2      sandwich_3.0-2      expm_0.999-7       
# [51] Exact_3.2           labelled_2.10.0     stringi_1.7.12      highr_0.10          klaR_1.7-2         
# [56] AlgDesign_1.2.1     nortest_1.0-4       e1071_1.7-13        boot_1.3-28         RgoogleMaps_1.4.5.3
# [61] rlang_1.1.0         pkgconfig_2.0.3     bitops_1.0-7        matrixStats_0.63.0  lattice_0.20-45    
# [66] labeling_0.4.2      tidyselect_1.2.0    magrittr_2.0.3      R6_2.5.1            DescTools_0.99.48  
# [71] generics_0.1.3      multcomp_1.4-24     combinat_0.0-8      pillar_1.9.0        haven_2.5.2        
# [76] withr_2.5.0         survival_3.2-13     abind_1.4-5         sp_1.6-0            questionr_0.7.8    
# [81] utf8_1.2.3          tzdb_0.3.0          jpeg_0.1-10         readxl_1.4.2        data.table_1.14.8  
# [86] digest_0.6.31       xtable_1.8-4        httpuv_1.6.9        stats4_4.1.2        munsell_0.5.0  
