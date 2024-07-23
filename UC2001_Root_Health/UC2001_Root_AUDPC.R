# Aidan Shands

# UC2001 Root Health AUDPC

library(ggplot2)
library(agricolae)
library(RColorBrewer)
library(dplyr)
library(car)
library(cowplot)
library(ggpubr)
library(multcompView)
library(fdth)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(colorspace)
library(FSA)
library(rcompanion)

# Import data 
Root_data = read.csv("Final_UC2001_Root_Rating_Percent_Diseased_for_AUDPC.csv", header=T)
Root_Health_data = read.csv("UC2001_Root_Health_Percent_Diseased.csv", header=T)

# Southern Isolates
South = c("P336", "P467", "P482", "Pc2109", "Pc2120", "Pc5241")
North = c("P344", "P444", "Pc2110", "Pc2113", "Pc2114")
Mexican = c("SER5_6", "TANR3_5", "TANS5_22")

labsX = c("P344", "P444", "Pc2110", "Pc2113", "Pc2114",
          "P336", "P467", "P482", "Pc2109", "Pc2120", "Pc5241",
          "SER5_6", "TANR3_5", "TANS5_22", "H20")

# Colors
LocationCols <- c( "Mexico" = "#E7298A", "CA-South" = "#7570B3", 
                   "CA-North" = "#D95F02")
#-------------------------------------------------------------------------------
# Defining Functions 
#-------------------------------------------------------------------------------
# QQ plot, Levene's test & Shapiro-Wilks test
QC_Stats = function(data, variable){
  dataname = deparse(substitute(data))
  col <- data[[variable]]
  print(paste("Levene's Test for: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Isolate), data = data))
  print(paste("Bartlett's Test for: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Isolate), data = data)$p.value)
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
  Subset_df = subset(input_df, select = -c(Location, Region))
  # Creating a list of Isolates_w_Rep for loop
  isolates = unlist(Subset_df$Isolate_w_Rep, recursive = FALSE) 
  # a loop that calculates the AUDPC for each isolate rep, and outputs to a df 
  for (i in isolates){
    name = toString(i) # get name as a string
    # get position of the list that will match with the input data (AUDPC_Data)
    position = match(i,isolates)
    # convert the row at 'position' to a vector
    vector = as.vector(t(Subset_df[position,1-2])) 
    AUDPC_Res = audpc(vector, time, type = "absolute") # generate the AUDPC
    AUDPC_Value = AUDPC_Res[[1]]# extract the AUDPC value we want
    # append to the df
    output_df = rbind(output_df, data.frame(name,AUDPC_Value))
  }
  colnames(output_df) = colnames2 # rename columns of new df & Inspect
  # Adding back in the Fruit_Number, Experiment and Isolate Columns 
  output_df$Isolate <- gsub("-\\d+$", "", output_df$Isolate)
  output_df$Location <- input_df$Location
  output_df$Region <- input_df$Region
  print("Success!")
  invisible(output_df)
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
# Overall Root Health  
#-------------------------------------------------------------------------------
QC_Stats(Root_Health_data, "Root_health")
Root_Health_THSD = Stats(Root_Health_data, "Isolate", "Root_health", "Tukey")
Root_Health_FLS = Stats(Root_Health_data, "Isolate", "Root_health", "Fisher")
Root_Health_DUNN = NP_Stats(Root_Health_data, "Isolate", "Root_health")

Root_Health_Sum <- as.data.frame(Root_Health_data %>% group_by(Isolate, Location) %>% 
                                   dplyr::summarize(Mean = mean(Root_health, na.rm=F), 
                                                    SD = sd(Root_health, na.rm=F)))

Root_Health_Sum = Root_Health_Sum %>% mutate(across(where(is.numeric), ~ round(., 2)))
# Create Results summary 
RH_THSD_j = subset(Root_Health_THSD, select = c(Isolate, Groups))
colnames(RH_THSD_j) <- c('Isolate', 'THSD_Groups')
RH_FLSD_j = subset(Root_Health_FLS, select = c(Isolate, Groups))
colnames(RH_FLSD_j) <- c('Isolate', 'FLSD_Groups')
RH_DUNN_j = subset(Root_Health_DUNN, select = c(Isolate, Dunn_Group))
colnames(RH_DUNN_j) <- c('Isolate', 'Dunn_Groups')

Sum_List = list(Root_Health_Sum, RH_THSD_j, RH_FLSD_j, RH_DUNN_j)
Comb_Summary = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List) 
colnames(Comb_Summary) <- c('Isolate', 'Location',"Mean_Root_health_as_Percent_Diseased", 
                            'SD', 'THSD_Groups', 'FLSD_Groups','Dunn_Groups')
#write_csv(Comb_Summary, "Root_Health_Summary.csv")

dfs3 = list(Root_Health_THSD, Root_Health_FLS, Root_Health_DUNN)
modified_dfs3 <- list()

for (df in dfs3) {
  df$Location = ifelse(df$Isolate %in% South, "CA-South",
                       ifelse(df$Isolate %in% North, "CA-North",
                              ifelse(df$Isolate %in% Mexican, "Mexico", NA)))
  
  df$Region = ifelse(df$Isolate %in% South, "CA",
                     ifelse(df$Isolate %in% North, "CA",
                            ifelse(df$Isolate %in% Mexican, "MX", NA)))
  
  # Append the modified dataframe to the list
  modified_dfs3[[length(modified_dfs3) + 1]] = df
}

Root_Health_THSD = modified_dfs3[[1]]
Root_Health_FLS = modified_dfs3[[2]]
Root_Health_DUNN = modified_dfs3[[3]]

labsX = c("P344", "P444", "Pc2110", "Pc2113", "Pc2114",
          "P336", "P467", "P482", "Pc2109", "Pc2120", "Pc5241",
          "SER5_6", "TANR3_5", "TANS5_22", "H20")

# Bar plot 
Bar <- ggplot(Root_Health_Sum, aes(x = interaction(Isolate, Location), y=Mean, fill=as.factor(Location), group=as.factor(Location))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(0.9))+ 
  ylab('% Diseased Roots') + 
  xlab("Isolate") +
  geom_text(data = Root_Health_DUNN, aes(label = Dunn_Group, y = max(Root_Health_Sum$Mean) + 23), vjust = -.3)+
  scale_fill_manual(values=LocationCols, name="Location") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12),
        legend.position = "none") +
  scale_x_discrete(labels = labsX)

Bar = Bar + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 
Bar
#ggsave(Bar,filename = "UC2001_Barplot.Dunn.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

#-------------------------------------------------------------------------------
# Calculating the AUDPC of the Overall Root Health  
#-------------------------------------------------------------------------------
Root_AUDPC = Calculate_AUDPC(Root_data, "Root")
print(head(Root_AUDPC))
# Normality Check 
QC_Stats(Root_AUDPC, "AUDPC_Values")

# 1-way ANOVA, Tukey HSD & Fisher's LSD
AUDPC_TT = Stats(Root_AUDPC, "Isolate", "AUDPC_Values", "Tukey")
AUDPC_LSD = Stats(Root_AUDPC, "Isolate", "AUDPC_Values", "Fisher")

AUDPC_Sum <- as.data.frame(Root_AUDPC %>% group_by(Isolate, Location) %>% 
                         dplyr::summarize(Mean = mean(AUDPC_Values, na.rm=F), 
                                          SD = sd(AUDPC_Values, na.rm=F)))
AUDPC_Sum = AUDPC_Sum %>% mutate(across(where(is.numeric), ~ round(., 2)))
# Create Results summary 
AUDPC_THSD_j = subset(AUDPC_TT, select = c(Isolate, Groups))
colnames(AUDPC_THSD_j) <- c('Isolate', 'THSD_Groups')
AUDPC_FLSD_j = subset(AUDPC_LSD, select = c(Isolate, Groups))
colnames(AUDPC_FLSD_j) <- c('Isolate', 'FLSD_Groups')
Sum_List2 = list(AUDPC_Sum, AUDPC_THSD_j, AUDPC_FLSD_j)
Comb_Summary2 = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List2) 
colnames(Comb_Summary2) <- c('Isolate', 'Location',"Mean_AUDPC", 
                            'SD', 'THSD_Groups', 'FLSD_Groups')

#write_csv(Comb_Summary2, "UC2001_AUDPC_Summary.csv")

# Location Summary

AUDPC_Sum_Loc <- as.data.frame(Root_AUDPC %>% group_by(Location) %>% 
                             dplyr::summarize(Mean = mean(AUDPC_Values, na.rm=F), 
                                              SD = sd(AUDPC_Values, na.rm=F)))
AUDPC_Sum_Loc = AUDPC_Sum_Loc %>% mutate(across(where(is.numeric), ~ round(., 2)))
#write_csv(AUDPC_Sum_Loc, "UC2001_AUDPC_Location_Summary.csv")

AUDPC_dfs = list(AUDPC_TT, AUDPC_LSD)
modified_AUDPC_dfs <- list()

for (df in AUDPC_dfs) {
  df$Location = ifelse(df$Isolate %in% South, "CA-South",
                       ifelse(df$Isolate %in% North, "CA-North",
                              ifelse(df$Isolate %in% Mexican, "Mexico", NA)))
  
  df$Region = ifelse(df$Isolate %in% South, "CA",
                     ifelse(df$Isolate %in% North, "CA",
                            ifelse(df$Isolate %in% Mexican, "MX", NA)))
  
  # Append the modified dataframe to the list
  modified_AUDPC_dfs[[length(modified_AUDPC_dfs) + 1]] = df
}

AUDPC_TT = modified_AUDPC_dfs[[1]]
AUDPC_LSD = modified_AUDPC_dfs[[2]]

# Subsetting CA only For North v South comparisons
AUDPC_CA = subset(Root_AUDPC, Region=="CA")

QC_Stats(AUDPC_CA, "AUDPC_Values")
AUDPC_TT_CA = Stats(AUDPC_CA, "Location", "AUDPC_Values", "Tukey")
AUDPC_LSD_CA = Stats(AUDPC_CA, "Location", "AUDPC_Values", "Fisher")
colnames(AUDPC_TT_CA)[1] <- "Location"
colnames(AUDPC_LSD_CA)[1] <- "Location"


#-------------------------------------------------------------------------------
# Visualizing AUDPC 
#-------------------------------------------------------------------------------
# All boxplots
BP1 = ggplot(Root_AUDPC, aes(x = interaction(Isolate, Location), y=AUDPC_Values, fill=as.factor(Location))) + 
  geom_boxplot() +
  geom_text(data = AUDPC_LSD, 
            aes(label = Groups, y = max(Root_AUDPC$AUDPC_Values)+5),
            vjust = -0.5) +
  scale_fill_manual(values=LocationCols, name="Location") +
  scale_x_discrete(labels = labsX)

BP1 = BP1 + theme_classic() + 
  labs(x="Isolate", y = "AUDPC")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12),
        legend.position = "none") + coord_cartesian(clip = "off")
BP1
#ggsave(BP1,filename = "UC2001_AUDPC.FLSD.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

# CA Boxplots  
BP5 = ggplot(AUDPC_CA, aes(x = Location, y=AUDPC_Values, fill=as.factor(Location))) + 
  geom_boxplot() +
  geom_text(data = AUDPC_LSD_CA, 
            aes(label = Groups, y = max(AUDPC_CA$AUDPC_Values) + 0.5),
            vjust = -0.5) +
  scale_fill_manual(values=LocationCols, name="Location") 

BP5 = BP5 + theme_classic() + 
  labs(x="Location", y = "AUDPC")+
  theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12),
        legend.position = "none") + coord_cartesian(clip = "off")
BP5
#ggsave(BP5,filename = "UC2001_AUDPC_CA_North_v_South.FLSD.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)


