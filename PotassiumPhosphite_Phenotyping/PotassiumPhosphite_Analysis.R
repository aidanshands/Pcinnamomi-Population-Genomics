#-------------------------------------------------------------------------------
# Aidan Shands 
# Manosalva Lab @ University of California, Riverside
#-------------------------------------------------------------------------------
### ----------------------------------------------------------------------------
# Load packages
### ----------------------------------------------------------------------------
library(ggplot2)
library(agricolae)
library(car)
library(dplyr)
library(FSA)
library(tidyverse)
library(rcompanion)
library(ggthemes)
#-------------------------------------------------------------------------------
setwd("/Users/manosalvalab/Desktop/Aidan/Chapter2_Pc_Population_Genomics/Pc_PopGen_Paper/2024_Revisions/GitHub/PotassiumPhosphite_Phenotyping")
PPP_Data <- read.csv("PPEC50_Data.csv", header = T)

Cols <- c("South" = "#7570B3", "North" = "#D95F02")

PPP_Data$Experiment <- factor(PPP_Data$Experiment)
PPP_Data$Growing_Region <- factor(PPP_Data$Growing_Region)
PPP_Data$Isolate <- factor(PPP_Data$Isolate)
table(PPP_Data$Growing_Region)

#-------------------------------------------------------------------------------
# Defining Functions 
#-------------------------------------------------------------------------------
# QQ plot, Levene's test & Shapiro-Wilks test
QC_Stats = function(data, variable){
  dataname = deparse(substitute(data))
  col <- data[[variable]]
  print(paste("Levene's Test by replicate for: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Replicate), data = data))
  print(paste("Levene's Test by Isolate for: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Isolate), data = data))
  print(paste("Bartlett's Test by replicate for: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Replicate), data = data)$p.value)
  print(paste("Bartlett's Test by Isolate for: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Isolate), data = data)$p.value)
  print(paste("Shapiro Wilks test for: ", dataname, sep=""))
  print(shapiro.test(col))
  print(paste("QQplot: ", dataname, sep=""))
  qqPlot(col)
}

# Performs One-way ANOVA, Tukey's HSD or Fisher's LSD
P_Stats = function(data, variable, value, test){
  dataname = deparse(substitute(data)) # get name of data we are working on 
  print(paste("Stats for ", dataname, sep=""))
  print(paste("One-Way ANOVA for: ", dataname, sep=""))
  # generate a formula (value~variable)
  formula <- reformulate(variable, response = value)
  AOV = aov(formula, data = data) # run one-way anova
  outnameAOV = summary(AOV)# generate summary
  print(outnameAOV) # print Summary
  if (test == "Fisher"){
    TestResults = (LSD.test(AOV, factor(variable), group=TRUE))# run Fisher's LSD
    print(TestResults)
  } else if (test == "Tukey"){
    TestResults = (HSD.test(AOV, factor(variable), group=TRUE))# run Tukey HSD
    print(TestResults)
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

#-------------------------------------------------------------------------------
# QC Stats
#-------------------------------------------------------------------------------
QC_Stats(PPP_Data, "EC50")
QC_Stats(PPP_Data, "Log10_EC50")

#-------------------------------------------------------------------------------
# One-way ANOVA looking at EC50 by Isolate & Tukey's HSD/Fisher's LSD
#-------------------------------------------------------------------------------
IS_TT = P_Stats(PPP_Data, "Isolate", "EC50", "Tukey")
IS_FLSD = P_Stats(PPP_Data, "Isolate", "EC50", "Fisher")
IS_DUNN = NP_Stats(PPP_Data, "Isolate", "EC50")

#-------------------------------------------------------------------------------
# Subset CA isolates 
#-------------------------------------------------------------------------------
CA = subset(PPP_Data, State == "CA")

#-------------------------------------------------------------------------------
# One-way ANOVA looking at EC50 by Growing Region & Tukey's HSD/Fisher's LSD
#-------------------------------------------------------------------------------
IS_TT_CA = P_Stats(CA, "Growing_Region", "EC50", "Tukey")
print(IS_TT_CA)

IS_FLSD_CA = P_Stats(CA, "Growing_Region", "EC50", "Fisher")
print(IS_FLSD_CA)
print(kruskal.test(EC50 ~Growing_Region, data = CA)) # Significant p<0.05

#-------------------------------------------------------------------------------
# Summarize Data
#-------------------------------------------------------------------------------
PPP_SUM <- as.data.frame(PPP_Data %>% group_by(Isolate, County, State, Country, Growing_Region) %>% 
                         dplyr::summarize(Mean = mean(EC50, na.rm=TRUE), SD = sd(EC50, na.rm=TRUE)))

PPP_SUM = PPP_SUM %>% mutate(across(where(is.numeric), ~ round(., 2)))

# joining TT, FLSD & Dunn groups for each isolate to PPP_SUM
TT_j = subset(IS_TT, select = c(Isolate, Groups))
colnames(TT_j) <- c('Isolate', 'TT_Groups')
FLSD_j = subset(IS_FLSD, select = c(Isolate, Groups))
colnames(FLSD_j) <- c('Isolate', 'FLSD_Groups')
DUNN_j = subset(IS_DUNN, select = c(Isolate, Dunn_Group))
colnames(DUNN_j) <- c('Isolate', 'Dunn_Groups')

Sum_List = list(PPP_SUM, TT_j, FLSD_j, DUNN_j)
# Summary used in Table S5
Comb_Summary = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List) 
colnames(Comb_Summary) <- c('Isolate', 'County', 'State',  'Country', 'Growing_Region',
                            "Mean_EC50", 'SD', 'TT_Groups', 'FLSD_Groups', 'Dunn_Groups')
#write.csv(Comb_Summary, file = "EC50_Sum.csv")

#-------------------------------------------------------------------------------
# Subset CA from Summarized Data
#-------------------------------------------------------------------------------
CA_SUM = subset(Comb_Summary, State == "CA")
table(CA_SUM$Growing_Region)
#-------------------------------------------------------------------------------
# Visualizations
#-------------------------------------------------------------------------------
#-----------
# Figure 4B
#-----------
BP4 = ggplot(CA_SUM, aes(x=Growing_Region, y=Mean_EC50, fill= Growing_Region)) + geom_boxplot()
BP4= BP4 + theme_few() + scale_fill_manual(values=Cols) + 
  labs(x="Growing Region", y = expression(EC[50] ~ paste("(",mu, "g/ml)")))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10),
        legend.position = "none")
BP4
#ggsave("PPP_EC50_Boxplot_Region_Paper.v4.pdf",BP4, width=4.28, height=3.78, units="in")

# Since the Kruskal-wallis test was significant, I'm placing a, and b above the 
# boxes to indicate significance. 

#-----------
# Figure 4A
#-----------
# Scott's bin
SB = hist(Comb_Summary$Mean_EC50, breaks = "Scott", xlab = "EC50", main = "Scott's bin of Mean EC50 values")
SB_Breaks = SB$breaks
SB_df <- data.frame(x = SB$mids, freq = SB$counts)
SB_Final = ggplot(data = SB_df, aes(x = x, y = freq)) +
  geom_bar(stat = "identity") + theme_few()+
  scale_x_continuous(breaks = SB_Breaks) +
  labs(y = "Frequency", x= expression(EC[50] ~ paste("(",mu, "g/ml)")))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10))
SB_Final
#ggsave("PPP_EC50_Scotts_Bin_Paper.v4.pdf",SB_Final, width=4.28, height=3.78, units="in")
