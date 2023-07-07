# Aidan Shands
### ----------------------------------------------------------------------------
# Load packages
### ----------------------------------------------------------------------------
library(ggplot2)
library(agricolae)
library(RColorBrewer)
library(plyr)
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

### ----------------------------------------------------------------------------
# Importing data, & Setting factors 
### ----------------------------------------------------------------------------
OGT_Data <- read.csv("GR_Data.csv", header = T)
# Setting factors 
OGT_Data$Isolate <- factor(OGT_Data$Isolate)
OGT_Data$User <- factor(OGT_Data$User)
OGT_Data$Date_D4 <- factor(OGT_Data$Date_D4)
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
  col <- data[[variable]]
  print(paste("Levene's Test for: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Replicate), data = data))
  print(paste("Bartlett's Test for: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Replicate), data = data)$p.value)
  print(paste("Shapiro Wilks test for: ", dataname, sep=""))
  print(shapiro.test(col))
  print(paste("QQplot: ", dataname, sep=""))
  qqPlot(col)
}
# Performs One-way ANOVA, Tukey's HSD or Fisher's LSD
Stats = function(data, variable, value, test){
  dataname = deparse(substitute(data)) # get name of data we are working on 
  print(paste("Stats for ", dataname, sep=""))
  print(paste("One-Way ANOVA for: ", dataname, sep=""))
  # generate a formula (value~variable)
  formula <- reformulate(variable, response = value)
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

# Summarize data (mean+sd)
detach(package:dplyr, unload=TRUE)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
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
ONA_TT_22 = Stats(OGT_Data, "Isolate", "RGPD_22", "Tukey")
# 25˚C 
ONA_TT_25 = Stats(OGT_Data, "Isolate", "RGPD_25", "Tukey")
# 28˚C 
ONA_TT_28 = Stats(OGT_Data, "Isolate", "RGPD_28", "Tukey")

### ----------------------------------------------------------------------------
# One-way ANOVA & Fisher's LSD Test 
### ----------------------------------------------------------------------------
# ANOVA by mean radial growth by isolate
# 22˚C 
ONA_FLSD_22 = Stats(OGT_Data, "Isolate", "RGPD_22", "Fisher")
# 25˚C 
ONA_FLSD_25 = Stats(OGT_Data, "Isolate", "RGPD_25", "Fisher")
# 28˚C 
ONA_FLSD_28 = Stats(OGT_Data, "Isolate", "RGPD_28", "Fisher")

### ----------------------------------------------------------------------------
# Summarizing data (Avg + SD)
### ----------------------------------------------------------------------------
# Summarizing 22 
SUMMARIZED_22 <- data_summary(OGT_Data, varname="RGPD_22", 
                              groupnames=c("Isolate", "County", "State", "Country", "Growing_Region"))
# Summarizing 25
SUMMARIZED_25 <- data_summary(OGT_Data, varname="RGPD_25", 
                              groupnames=c("Isolate", "County", "State", "Country", "Growing_Region"))
# Summarizing 28
SUMMARIZED_28 <- data_summary(OGT_Data, varname="RGPD_28", 
                              groupnames=c("Isolate", "County", "State", "Country", "Growing_Region"))

# rounding the decimals to 2 decimals
SUMMARIZED_22 = SUMMARIZED_22 %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUMMARIZED_25 = SUMMARIZED_25 %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUMMARIZED_28 = SUMMARIZED_28 %>% mutate(across(where(is.numeric), ~ round(., 2)))
# renaming columns
colnames(SUMMARIZED_22) <- c('Isolate', 'County', 'State', 'Country', "Growing_Region", 'RGPD_22', 'RGPD_22_SD')
colnames(SUMMARIZED_25) <- c('Isolate', 'County', 'State', 'Country', "Growing_Region", 'RGPD_25', 'RGPD_25_SD')
colnames(SUMMARIZED_28) <- c('Isolate', 'County', 'State', 'Country', "Growing_Region", 'RGPD_28', 'RGPD_28_SD')
# exporting summaries
#dir.create("Summary_Outputs_Final")
#write.csv(SUMMARIZED_22, "All_Data_22SUM.csv")
#write.csv(SUMMARIZED_25, "All_Data_25SUM.csv")
#write.csv(SUMMARIZED_28, "All_Data_28SUM.csv")

# Creating summary df for plotting
SUMMARIZED_22_j = subset(SUMMARIZED_22, select = c(Isolate, County, State, Growing_Region, Country, RGPD_22, RGPD_22_SD))
SUMMARIZED_25_j = subset(SUMMARIZED_25, select = c(Isolate, RGPD_25, RGPD_25_SD))
SUMMARIZED_28_j = subset(SUMMARIZED_28, select = c(Isolate, RGPD_28, RGPD_28_SD))
TT22_j = subset(ONA_TT_22, select = c(Isolate, Groups))
colnames(TT22_j) <- c('Isolate', 'TT_22_Groups')
TT25_j = subset(ONA_TT_25, select = c(Isolate, Groups))
colnames(TT25_j) <- c('Isolate', 'TT_25_Groups')
TT28_j = subset(ONA_TT_28, select = c(Isolate, Groups))
colnames(TT28_j) <- c('Isolate', 'TT_28_Groups')
FLSD22_j = subset(ONA_FLSD_22, select = c(Isolate, Groups))
colnames(FLSD22_j) <- c('Isolate', 'FLSD_22_Groups')
FLSD25_j = subset(ONA_FLSD_25, select = c(Isolate, Groups))
colnames(FLSD25_j) <- c('Isolate', 'FLSD_25_Groups')
FLSD28_j = subset(ONA_FLSD_28, select = c(Isolate, Groups))
colnames(FLSD28_j) <- c('Isolate', 'FLSD_28_Groups')
Sum_List = list(SUMMARIZED_22_j, SUMMARIZED_25_j, SUMMARIZED_28_j, TT22_j, TT25_j, TT28_j, FLSD22_j, FLSD25_j, FLSD28_j)
# This data is used in Table S4
Comb_Summary = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List) 
#write.csv(Comb_Summary,"All_Data_Summarized.csv")
### ----------------------------------------------------------------------------
# Subsetting CA only and looking at North vs South 
### ----------------------------------------------------------------------------
CA_Isolates = subset(OGT_Data, State == "CA")
table(CA_Isolates$Growing_Region)

# QC Stats 
# 22˚C  
QC_Stats(CA_Isolates, "RGPD_22")
# 25˚C
QC_Stats(CA_Isolates, "RGPD_25")
# 28˚C
QC_Stats(CA_Isolates, "RGPD_28")

# ANOVA & Tukey's HSD with RGR by Growing Region
# 22˚C 
CA_ONA_TT_22 = Stats(CA_Isolates, "Growing_Region", "RGPD_22", "Tukey")
# 25˚C 
CA_ONA_TT_25 = Stats(CA_Isolates, "Growing_Region", "RGPD_25", "Tukey")
# 28˚C 
CA_ONA_TT_28 = Stats(CA_Isolates, "Growing_Region", "RGPD_28", "Tukey")

# ANOVA & Fisher's HSD with RGR by Growing Region
# 22˚C 
CA_ONA_FLSD_22 = Stats(CA_Isolates, "Growing_Region", "RGPD_22", "Fisher")
# 25˚C 
CA_ONA_FLSD_25 = Stats(CA_Isolates, "Growing_Region", "RGPD_25", "Fisher")
# 28˚C 
CA_ONA_FLSD_28 = Stats(CA_Isolates, "Growing_Region", "RGPD_28", "Fisher")

# summarizing CA data =
SUMMARIZED_22_CA <- data_summary(CA_Isolates, varname="RGPD_22", 
                                 groupnames=c("Isolate", "County", "State", "Country", "Growing_Region"))
SUMMARIZED_25_CA <- data_summary(CA_Isolates, varname="RGPD_25", 
                                 groupnames=c("Isolate", "County", "State", "Country", "Growing_Region"))
SUMMARIZED_28_CA <- data_summary(CA_Isolates, varname="RGPD_28", 
                                 groupnames=c("Isolate", "County", "State", "Country", "Growing_Region"))
# rounding the decimals to 2 decimals
SUMMARIZED_22_CA = SUMMARIZED_22_CA %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUMMARIZED_25_CA = SUMMARIZED_25_CA %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUMMARIZED_28_CA = SUMMARIZED_28_CA %>% mutate(across(where(is.numeric), ~ round(., 2)))

### ----------------------------------------------------------------------------
# Performing Temperature Comparisons
### ----------------------------------------------------------------------------

OGT_Data2 <- read.csv("GR_Temperature_Data.csv", header = T)
head(OGT_Data2)

# Setting factors 
OGT_Data2$Isolate <- factor(OGT_Data2$Isolate)
OGT_Data2$User <- factor(OGT_Data2$User)
OGT_Data2$Date_D4 <- factor(OGT_Data2$Date_D4)
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
ONA_TT_Temp = Stats(OGT_Data2, "Temperature", "RGPD", "Tukey")

### ----------------------------------------------------------------------------
# One-way ANOVA & Fisher's LSD Test 
### ----------------------------------------------------------------------------
ONA_FLSD_Temp = Stats(OGT_Data2, "Temperature", "RGPD", "Fisher")

### ----------------------------------------------------------------------------
# Subsetting for North & South 
### ----------------------------------------------------------------------------
CA_Isolates_So = subset(OGT_Data2, Growing_Region == "South")
CA_Isolates_No = subset(OGT_Data2, Growing_Region == "North")
CA_Isolates_All = subset(OGT_Data2, State == "CA")

### ----------------------------------------------------------------------------
# One-Way ANOVA looking at RGR by Temp in CA 
### ----------------------------------------------------------------------------
ONA_TT_Temp_CA = Stats(CA_Isolates_All, "Temperature", "RGPD", "Tukey")

ONA_FLSD_Temp_CA = Stats(CA_Isolates_All, "Temperature", "RGPD", "Fisher")

### ----------------------------------------------------------------------------
# One-Way ANOVA looking at RGR by Temp in CA-South
### ----------------------------------------------------------------------------
ONA_TT_Temp_CAS = Stats(CA_Isolates_So, "Temperature", "RGPD", "Tukey")

ONA_FLSD_Temp_CAS = Stats(CA_Isolates_So, "Temperature", "RGPD", "Fisher")

### ----------------------------------------------------------------------------
# One-Way ANOVA looking at RGR by Temp in CA-North
### ----------------------------------------------------------------------------
ONA_TT_Temp_CAN = Stats(CA_Isolates_No, "Temperature", "RGPD", "Tukey")

ONA_FLSD_Temp_CAN = Stats(CA_Isolates_No, "Temperature", "RGPD", "Fisher")

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

### ----------------------------------------------------------------------------
# One-Way ANOVA looking at Percent_Growth_Affected_by_28 by Growing_Region
### ----------------------------------------------------------------------------
ONA_TT_Temp_Ratio = Stats(OGT_Data3, "Growing_Region", "Percent_Growth_Affected_by_28", "Tukey")

ONA_FLSD_Temp_Ratio = Stats(OGT_Data3, "Growing_Region", "Percent_Growth_Affected_by_28", "Fisher")

### ----------------------------------------------------------------------------
# One-Way ANOVA looking at Percent_Growth_Affected_by_28 by Growing_Region CA only
### ----------------------------------------------------------------------------
# Subset CA 
CA_Isolates2 = subset(OGT_Data3, State == "CA")

ONA_TT_Temp_Ratio_CA = Stats(CA_Isolates2, "Growing_Region", "Percent_Growth_Affected_by_28", "Tukey")

ONA_FLSD_Temp_Ratio_CA = Stats(CA_Isolates2, "Growing_Region", "Percent_Growth_Affected_by_28", "Fisher")

# Percent growth vs Growing region (sig. diff. p<0.05)
AOV_Model_PercentGrowth2 <- aov(Percent_Growth_Affected_by_28 ~ Growing_Region, data = CA_Isolates2)
summary(AOV_Model_PercentGrowth2)

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

#----------
# Figure 3B
#----------
PTemp<- ggplot(CA_Isolates2, aes(x=Isolate, y=Percent_Growth_Affected_by_28,  fill= Growing_Region)) + 
  geom_bar(stat="identity", 
           position=position_dodge(width = 0.9), 
           aes(reorder(Isolate,Percent_Growth_Affected_by_28),Percent_Growth_Affected_by_28), width = 0.5)

PTemp = PTemp+labs(x="Isolate", y = "% growth affected from 22 \u00B0C to 28 \u00B0C") +
  theme_classic() + scale_fill_manual(values=Cols) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10,face="bold"),
        legend.position = "none") + coord_flip()
PTemp

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
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] multcompView_0.1-8 cowplot_1.1.1      plyr_1.8.8         ggtree_3.2.1       colorspace_2.1-0   lubridate_1.9.2   
# [7] forcats_1.0.0      stringr_1.5.0      purrr_1.0.1        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1      
# [13] tidyverse_2.0.0    reshape2_1.4.4     MASS_7.3-54        fdth_1.2-6         ggpubr_0.6.0       car_3.1-2         
# [19] carData_3.0-5      RColorBrewer_1.1-3 agricolae_1.3-5    ggplot2_3.4.2     

# loaded via a namespace (and not attached):
#   [1] seqinr_4.2-23       ggsignif_0.6.4      ellipsis_0.3.2      aplot_0.1.10        rstudioapi_0.14    
# [6] farver_2.1.1        fansi_1.0.4         splines_4.1.2       knitr_1.42          pegas_1.2          
# [11] ade4_1.7-22         jsonlite_1.8.4      poppr_2.9.4         broom_1.0.4         cluster_2.1.2      
# [16] png_0.1-8           shiny_1.7.4         compiler_4.1.2      httr_1.4.5          backports_1.4.1    
# [21] lazyeval_0.2.2      Matrix_1.3-4        fastmap_1.1.1       cli_3.6.1           later_1.3.0        
# [26] htmltools_0.5.5     tools_4.1.2         ggmap_3.0.2         igraph_1.4.1        gtable_0.3.3       
# [31] glue_1.6.2          dplyr_1.1.1         Rcpp_1.0.10         vctrs_0.6.1         ape_5.7-1          
# [36] nlme_3.1-153        xfun_0.38           timechange_0.2.0    mime_0.12           miniUI_0.1.1.1     
# [41] lifecycle_1.0.3     rstatix_0.7.2       scales_1.2.1        hms_1.1.3           promises_1.2.0.1   
# [46] parallel_4.1.2      yaml_2.3.7          ggfun_0.0.9         yulab.utils_0.0.6   labelled_2.10.0    
# [51] stringi_1.7.12      highr_0.10          klaR_1.7-2          AlgDesign_1.2.1     tidytree_0.4.2     
# [56] permute_0.9-7       boot_1.3-28         RgoogleMaps_1.4.5.3 rlang_1.1.0         pkgconfig_2.0.3    
# [61] bitops_1.0-7        polysat_1.7-7       evaluate_0.20       lattice_0.20-45     labeling_0.4.2     
# [66] treeio_1.18.1       patchwork_1.1.2     tidyselect_1.2.0    magrittr_2.0.3      R6_2.5.1           
# [71] generics_0.1.3      combinat_0.0-8      DBI_1.1.3           pillar_1.9.0        haven_2.5.2        
# [76] withr_2.5.0         mgcv_1.8-38         abind_1.4-5         sp_1.6-0            questionr_0.7.8    
# [81] utf8_1.2.3          tzdb_0.3.0          rmarkdown_2.21      jpeg_0.1-10         adegenet_2.1.10    
# [86] grid_4.1.2          vegan_2.6-4         digest_0.6.31       xtable_1.8-4        httpuv_1.6.9       
# [91] gridGraphics_0.5-1  munsell_0.5.0       ggplotify_0.1.0  

