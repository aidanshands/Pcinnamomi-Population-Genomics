# Aidan Shands

library(ggplot2)
library(agricolae)
library(car)
library(dplyr)
#-------------------------------------------------------------------------------
setwd("/Users/manosalvalab/Desktop/Aidan/Pc_PopGen_2022/PopGen_Scripts_Github/PotassiumPhosphite_Phenotyping")
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
    TestResults = (LSD.test(AOV, factor(variable), group=TRUE))# run Fisher's LSD
  } else if (test == "Tukey"){
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

#-------------------------------------------------------------------------------
# QC Stats
#-------------------------------------------------------------------------------
QC_Stats(PPP_Data, "EC50")
QC_Stats(PPP_Data, "Log10_EC50")

#-------------------------------------------------------------------------------
# One-way ANOVA looking at EC50 by Isolate & Tukey's HSD/Fisher's LSD
#-------------------------------------------------------------------------------
IS_TT = Stats(PPP_Data, "Isolate", "EC50", "Tukey")
IS_FLSD = Stats(PPP_Data, "Isolate", "EC50", "Fisher")

#-------------------------------------------------------------------------------
# Subset CA isolates 
#-------------------------------------------------------------------------------
CA = subset(PPP_Data, State == "CA")

#-------------------------------------------------------------------------------
# One-way ANOVA looking at EC50 by Growing Region & Tukey's HSD/Fisher's LSD
#-------------------------------------------------------------------------------
IS_TT_CA = Stats(CA, "Growing_Region", "EC50", "Tukey")
IS_FLSD_CA = Stats(CA, "Growing_Region", "EC50", "Fisher")

#-------------------------------------------------------------------------------
# Summarize Data
#-------------------------------------------------------------------------------
PPP_SUM <- data_summary(PPP_Data, varname="EC50", 
                        groupnames=c("Isolate", "County", "State", "Country", "Growing_Region"))

# joining TT groups for each isolate to PPP_SUM
TT_j = subset(IS_TT, select = c(Isolate, Groups))
colnames(TT_j) <- c('Isolate', 'TT_Groups')
Sum_List = list(PPP_SUM, TT_j)
# Summary used in Table S4
Comb_Summary = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List) 
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
# Figure 3C
#-----------
BP4 = ggplot(CA_SUM, aes(x=Growing_Region, y=EC50, fill= Growing_Region)) + geom_boxplot()
BP4= BP4 + theme_classic() + scale_fill_manual(values=Cols) + 
  labs(x="Growing Region", y = "EC50")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10,face="bold"),
        legend.position = "none")
BP4
#ggsave("PPP_EC50_Boxplot_Region.pdf",BP4, width=6, height=4, units="in")

#-----------
# Figure 2B
#-----------
# Scott's bin
SB = hist(PPP_SUM$EC50, breaks = "Scott", xlab = "EC50", main = "Scott's bin of Mean EC50 values")
SB_Breaks = SB$breaks
SB_df <- data.frame(x = SB$mids, freq = SB$counts)
SB_Final = ggplot(data = SB_df, aes(x = x, y = freq)) +
  geom_bar(stat = "identity") + theme_classic()+
  scale_x_continuous(breaks = SB_Breaks) +
  labs(y = "Frequency", x= expression(EC[50]))+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10,face="bold"))
SB_Final
#ggsave("PPP_EC50_Scotts_Bin.pdf",SB_Final, width=6, height=4, units="in")


