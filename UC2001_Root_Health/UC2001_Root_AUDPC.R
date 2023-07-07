# Aidan Shands

# UC2001 Root Health AUDPC

library(ggplot2)
library(agricolae)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(car)
library(cowplot)
library(ggpubr)
library(multcompView)
library(fdth)
library(MASS)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(colorspace)
library(ggtree)

#-------------------------------------------------------------------------------
# Import data 
Root_data = read.csv("Final_UC2001_Root_Rating_Percent_Diseased.csv", header=T)
Root_Health_data = read.csv("UC2001_Root_Health.csv", header=T)
#-------------------------------------------------------------------------------
# Defining Functions 
#-------------------------------------------------------------------------------
# QQ plot, Levene's test & Shapiro-Wilks test
QC_Stats = function(data, variable){
  dataname = deparse(substitute(data))
  col <- data[[variable]]
  print(paste("Levene's Test for: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Isolate), data = data))
  print(paste("Shapiro Wilks test for: ", dataname, sep=""))
  print(shapiro.test(col))
  print(paste("QQplot: ", dataname, sep=""))
  qqPlot(col)
}
# Calculate AUDPC ( need to specify Root Time or Health Time )
Calculate_AUDPC = function(input_df, data){
  Root_Time = c(2,3.5,5)
  Health_Time = c(0,1,2,3,4,5)
  # conditional to determine which time to use 
  if (data == "Root") {
    time <- Root_Time
  } else if (data == "Health") {
    time <- Health_Time
  } else {
    stop("Invalid time argument. Please provide 'Root_Time' or 'Health_Time'.")
  }
  colnames2 = c("Isolate", "AUDPC_Values")# Setting new col names 
  output_df = NULL # defining empty df to populate
  # Creating a list of Isolates_w_Rep for loop
  isolates = unlist(input_df$Isolate_w_Rep, recursive = FALSE) 
  # a loop that calculates the AUDPC for each isolate rep, and outputs to a df 
  for (i in isolates){
    name = toString(i) # get name as a string
    # get position of the list that will match with the input data (AUDPC_Data)
    position = match(i,isolates)
    # convert the row at 'position' to a vector
    vector = as.vector(t(input_df[position,1-2])) 
    AUDPC_Res = audpc(vector, time, type = "absolute") # generate the AUDPC
    AUDPC_Value = AUDPC_Res[[1]]# extract the AUDPC value we want
    # append to the df
    output_df = rbind(output_df, data.frame(name,AUDPC_Value))
  }
  colnames(output_df) = colnames2 # rename columns of new df & Inspect
  # Adding back in the Fruit_Number, Experiment and Isolate Columns 
  output_df$Isolate <- gsub("-\\d+$", "", output_df$Isolate)
  print("Success!")
  return(output_df)
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

#-------------------------------------------------------------------------------
# Root Health Analysis 
#-------------------------------------------------------------------------------
# Here I am calculating the AUDPC of the seedling root health ratings of the 
# treatments via visual assessment. In this experiment, I have 3 data points 
# from 3 separate trees, respectively, that were ultimately discarded after 
# visual assessment at 2 and 3.5 WPI. For the remaining 6 trees assessed at the 
# end of the experiment, I wrote a python script (Calculate_Random_Pairs.py) to 
# randomly pair the remaining 6 data points for each treatment and take the 
# average of the pairs. This process leaves 3 data points remaining for the 
# AUDPC calculation. 
#-------------------------------------------------------------------------------

# Calculating the AUDPC of the Overall Root Health  
Calculate_AUDPC(Root_data, "Root")
# Casting results to new data frame 
Root_AUDPC = Calculate_AUDPC(Root_data, "Root")

#-------------------------------------------------------------------------------
# Assessing Normality and Performing One-way ANOVA and Tukey HSD
#-------------------------------------------------------------------------------

# Normality Check 
QC_Stats(Root_AUDPC, "AUDPC_Values")
# 1-way ANOVA & Tukey HSD
Root_THSD_Res = Stats(Root_AUDPC, "Isolate", "AUDPC_Values", "Tukey")
Root_FLSD_Res = Stats(Root_AUDPC, "Isolate", "AUDPC_Values", "Fisher")

# Visualizing the Root AUDPC Results with Tukey's HSD Groups 
BP1 = ggplot(Root_AUDPC, aes(x=Isolate, y=AUDPC_Values)) + 
  geom_boxplot(aes(reorder(Isolate,AUDPC_Values),AUDPC_Values)) +
  geom_text(data = Root_THSD_Res, 
            aes(label = Groups, y = max(Root_AUDPC$AUDPC_Values) + 0.5),
            vjust = -0.5)

BP1 = BP1 + theme_classic() + 
  labs(x="Isolate", y = "AUDPC")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold")) 
BP1 = BP1+labs(title=expression(paste("AUDPC of % healthy roots of UC2001 seedlings infected with ", italic("P. cinnamomi"))))

BP1

# Exporting Fisher LSD and Tukey HSD results 
#write.csv(Root_THSD_Res, "UC2001_Root_AUDPC_THSD_Final.csv")
#write.csv(Root_FLSD_Res, "UC2001_Root_AUDPC_FLSD_Final.csv")
#write.csv(Root_AUDPC, "UC2001_Root_AUDPC_Final.csv")

print(BP1)
# ggsave(BP1,filename = "Root_AUDPC_THSD.pdf", width = 8, height = 6, units = "in")
#-------------------------------------------------------------------------------
# Visualizing the Root AUDPC Results with Fisher's LSD Groups 
#-------------------------------------------------------------------------------
#-----------
# Figure 5B
#-----------
BP2 = ggplot(Root_AUDPC, aes(x=Isolate, y=AUDPC_Values)) + 
  geom_boxplot(aes(reorder(Isolate,AUDPC_Values),AUDPC_Values)) +
  geom_text(data = Root_FLSD_Res, 
            aes(label = Groups, y = max(Root_AUDPC$AUDPC_Values) + 0.5),
            vjust = -0.5)

BP2 = BP2 + theme_classic() + 
  labs(x="Isolate", y = "AUDPC")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold")) 
BP2 = BP2+labs(title=expression(paste("AUDPC of % healthy roots of UC2001 seedlings infected with ", italic("P. cinnamomi"))))

print(BP2)
#ggsave(BP2,filename = "Root_AUDPC_FLSD.pdf", width = 8, height = 6, units = "in")

#-------------------------------------------------------------------------------
# Performing one-way ANOVA, Tukey's HSD, and Fisher's LSD of the overall 
# root health (as % healthy roots) at 5 weeks post inoculation 
#-------------------------------------------------------------------------------
# Normality Check 
QC_Stats(Root_Health_data, "Root_health")
Root_Health_THSD_Res = Stats(Root_Health_data, "Isolate", "Root_health", "Tukey")
Root_Health_FLSD_Res = Stats(Root_Health_data, "Isolate", "Root_health", "Fisher")
Root_Health_Sum = data_summary(Root_Health_data, varname="Root_health", 
                               groupnames=c("Isolate"))

#-------------------------------------------------------------------------------
# Visualizing the Root Health Results with Tukey's HSD Groups 
#-------------------------------------------------------------------------------
BP3 <- ggplot(Root_Health_Sum, aes(x=Isolate, y=Root_health)) + 
  geom_bar(stat="identity", position=position_dodge(width = 0.6), 
           aes(reorder(Isolate,Root_health),Root_health), width = 0.6) +
  geom_errorbar(aes(ymin=Root_health-sd, ymax=Root_health+sd), width=.2,
                position=position_dodge(0.6)) +
  geom_text(data = Root_Health_THSD_Res, 
            aes(label = Groups, y = max(Root_Health_Sum$Root_health) + 3),
            vjust = -1) + labs(x="Isolate", y = "Root Health (% healthy)") + 
  theme_classic()
# Finished bar plot
BP3= BP3+labs(title=expression(paste("Root health of UC2001 seedlings infected with ", italic("P. cinnamomi")))) +
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))

print(BP3)
#ggsave(BP3, filename = "Root_Health_THSD.pdf", width = 8, height = 6, units = "in")

#-------------------------------------------------------------------------------
# Visualizing the Root Health Results with Fisher's LSD Groups 
#-------------------------------------------------------------------------------
#-----------
# Figure 5A
#-----------
BP4 <- ggplot(Root_Health_Sum, aes(x=Isolate, y=Root_health)) + 
  geom_bar(stat="identity", position=position_dodge(width = 0.6), 
           aes(reorder(Isolate,Root_health),Root_health), width = 0.6) +
  geom_errorbar(aes(ymin=Root_health-sd, ymax=Root_health+sd), width=.2,
                position=position_dodge(0.6)) +
  geom_text(data = Root_Health_FLSD_Res, 
            aes(label = Groups, y = max(Root_Health_Sum$Root_health) + 3),
            vjust = -1) + labs(x="Isolate", y = "Root Health (% healthy)") + 
  theme_classic()
# Finished bar plot
BP4 = BP4+labs(title=expression(paste("Root health of UC2001 seedlings infected with ", italic("P. cinnamomi")))) +
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12,face="bold"))

print(BP4)
#ggsave(BP4, filename = "Root_Health_FLSD.pdf", width = 8, height = 6, units = "in")

