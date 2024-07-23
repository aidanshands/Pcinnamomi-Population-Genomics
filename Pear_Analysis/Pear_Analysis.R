# Aidan Shands 

library(ggplot2)
library(agricolae)
library(RColorBrewer)
library(car)
library(cowplot)
library(ggpubr)
library(multcompView)
library(fdth)
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
# Lesion area data 
Pear_Data = read.csv("AUDPC_14_Isolates_AllReps.csv", header = T)
print(head(Pear_Data))

# Log transformation of the lesion area values 

#Pear_Data$DPI_3LT = log(Pear_Data$DPI_3)
#Pear_Data$DPI_4LT = log(Pear_Data$DPI_4)
#Pear_Data$DPI_5LT = log(Pear_Data$DPI_5)

# Square-root transformation of the lesion area values
#Pear_Data$DPI_3SQRT = sqrt(Pear_Data$DPI_3)
#Pear_Data$DPI_4SQRT = sqrt(Pear_Data$DPI_4)
#Pear_Data$DPI_5SQRT = sqrt(Pear_Data$DPI_5)


# Setting factors 
Pear_Data$Isolate_w_Rep <- factor(Pear_Data$Isolate_w_Rep)
Pear_Data$Fruit_Number <- factor(Pear_Data$Fruit_Number)
Pear_Data$Experiment <- factor(Pear_Data$Experiment)
Pear_Data$Isolate<- factor(Pear_Data$Isolate)
Pear_Data$Location<- factor(Pear_Data$Location)
Pear_Data$Region<- factor(Pear_Data$Region)
Pear_Data$DPI_3 <- as.numeric(Pear_Data$DPI_3)
Pear_Data$DPI_4 <- as.numeric(Pear_Data$DPI_4)
Pear_Data$DPI_5 <- as.numeric(Pear_Data$DPI_5)

# Lesion area summary data 
Pear_Sum_Data = read.csv("Pear_EXP_Sum_14Isolates_PAPER.csv", header = T)
print(head(Pear_Sum_Data))

LocationCols <- c( "Mexico" = "#E7298A", "CA-South" = "#7570B3", 
                   "CA-North" = "#D95F02")

# Southern Isolates
South = c("P336", "P467", "P482", "Pc2109", "Pc2120", "Pc5241")
North = c("P344", "P444", "Pc2110", "Pc2113", "Pc2114")
Mexican = c("SER5_6", "TANR3_5", "TANS5_22")

### ----------------------------------------------------------------------------
# Defining Functions
### ----------------------------------------------------------------------------
# QC Stats 
QC_Stats = function(data, variable){
  dataname = deparse(substitute(data))
  col <- data[[variable]]
  print(paste("Stats for ", dataname, sep=""))
  print(shapiro.test(col))
  print(paste("QQplot: ", dataname, sep=""))
  qqPlot(col)
  print(paste("Welch Two-sample T-test for: ", dataname, sep=""))
  Ttest = t.test(col ~ factor(Experiment), data = data, var.equal = FALSE)
  print(Ttest$p.value)
  print(paste("Student's t-test for: ", dataname, sep=""))
  Ttest2 = t.test(col ~ factor(Experiment), data = data, var.equal = T)
  print(Ttest2$p.value)
  print(paste("Bartlett's Test for equal variance between experiments: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Experiment), data = data)$p.value)
  print(paste("Levene's Test for equal variance between experiments: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Experiment), data = data))
  print(paste("Bartlett's Test for equal variance between isolates: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Isolate), data = data)$p.value)
  print(paste("Levene's Test for equal variance between isolates: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Isolate), data = data))
}

# QC Stats 
QC_Stats2 = function(data, variable){
  dataname = deparse(substitute(data))
  col <- data[[variable]]
  print(paste("Stats for ", dataname, sep=""))
  print(shapiro.test(col))
  print(paste("QQplot: ", dataname, sep=""))
  qqPlot(col)
  print(paste("Bartlett's Test for equal variance between isolates: ", dataname, sep=""))
  print(bartlett.test(col ~ factor(Isolate), data = data)$p.value)
  print(paste("Levene's Test for equal variance between isolates: ", dataname, sep=""))
  print(leveneTest(col ~ factor(Isolate), data = data))
}

# Calculate AUDPC
Calculate_AUDPC = function(input_df){
  time = c(3,4,5) # Setting time 
  # subset columns of interest
  AllReps_for_AUDPC = subset(input_df, 
                             select = -c(Fruit_Number, Experiment, Isolate, Location, Region))
  colnames2 = c("Isolate_w_Rep", "AUDPC_Values")# Setting new col names 
  output_df = NULL # defining empty df to populate
  # Creating a list of Isolates_w_Rep for loop
  isolates = unlist(AllReps_for_AUDPC$Isolate_w_Rep, recursive = FALSE) 
  # a loop that calculates the AUDPC for each isolate and rep and 
  # generates a new df with the results
  for (i in isolates){
    name = toString(i) # get name as a string
    # get position of the list that will match with the input data (AUDPC_Data)
    position = match(i,isolates) 
    # convert the row at 'position' to a vector
    vector = as.vector(t(AllReps_for_AUDPC[position,1-2])) 
    AUDPC_Res = audpc(vector, time, type = "absolute") # generate the AUDPC
    AUDPC_Value = AUDPC_Res[[1]]# extract the AUDPC value we want
    # append to the df
    output_df = rbind(output_df, data.frame(name,AUDPC_Value)) 
  }
  colnames(output_df) = colnames2 # rename columns of new df & Inspect
  # Adding back in the Fruit_Number, Experiment and Isolate Columns 
  output_df$Fruit_Number <- as.numeric(sub(".*-(\\d+)$", "\\1", 
                                           output_df$Isolate_w_Rep))
  output_df$Experiment <- ifelse(output_df$Fruit_Number <= 6, "EXP1", "EXP2")
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

### ----------------------------------------------------------------------------
# Lesion Area Comparisons
### ----------------------------------------------------------------------------
CA_Only = subset(Pear_Data, Region =="CA" )
head(Pear_Data)
#------------------
# 3 DPI (Combined)
#------------------
# Raw data
QC_Stats(Pear_Data, "DPI_3")
DPI3_TT = Stats(Pear_Data, "Isolate", "DPI_3", "Tukey")
DPI3_LSD = Stats(Pear_Data, "Isolate", "DPI_3", "Fisher")
DPI3_Dunn = NP_Stats(Pear_Data, "Isolate", "DPI_3")
# North v South 
CA_3DPI_TT = Stats(CA_Only, "Location", "DPI_3", "Tukey")
CA_3DPI_LSD = Stats(CA_Only, "Location", "DPI_3", "Fisher")

# Log-transformed data
#QC_Stats(Pear_Data, "DPI_3LT")
#DPI3_TT_LT = Stats(Pear_Data, "Isolate", "DPI_3LT", "Tukey")
#DPI3_LSD_LT = Stats(Pear_Data, "Isolate", "DPI_3LT", "Fisher")
#DPI3_Dunn_LT = NP_Stats(Pear_Data, "Isolate", "DPI_3LT")
# North v South 
#CA_3DPI_TT_LT = Stats(CA_Only, "Location", "DPI_3LT", "Tukey")
#CA_3DPI_LSD_LT = Stats(CA_Only, "Location", "DPI_3LT", "Fisher")

# Square-root transformed data
#QC_Stats(Pear_Data, "DPI_3SQRT")
#DPI3_TT_SQRT = Stats(Pear_Data, "Isolate", "DPI_3SQRT", "Tukey")
#DPI3_LSD_SQRT = Stats(Pear_Data, "Isolate", "DPI_3SQRT", "Fisher")
#DPI3_Dunn_SQRT = NP_Stats(Pear_Data, "Isolate", "DPI_3SQRT")
# North v South 
#CA_3DPI_TT_SQRT = Stats(CA_Only, "Location", "DPI_3SQRT", "Tukey")
#CA_3DPI_LSD_SQRT = Stats(CA_Only, "Location", "DPI_3SQRT", "Fisher")


#------------------
# 4 DPI (Combined)
#------------------
QC_Stats(Pear_Data, "DPI_4")
DPI4_TT = Stats(Pear_Data, "Isolate", "DPI_4", "Tukey")
DPI4_LSD = Stats(Pear_Data, "Isolate", "DPI_4", "Fisher")
DPI4_Dunn = NP_Stats(Pear_Data, "Isolate", "DPI_4")
# North v South 
CA_4DPI_TT = Stats(CA_Only, "Location", "DPI_4", "Tukey")
CA_4DPI_LSD = Stats(CA_Only, "Location", "DPI_4", "Fisher")

# Log-transformed data
#QC_Stats(Pear_Data, "DPI_4LT")
#DPI4_TT_LT = Stats(Pear_Data, "Isolate", "DPI_4LT", "Tukey")
#DPI4_LSD_LT = Stats(Pear_Data, "Isolate", "DPI_4LT", "Fisher")
#DPI4_Dunn_LT = NP_Stats(Pear_Data, "Isolate", "DPI_4LT")
# North v South 
#CA_4DPI_TT_LT = Stats(CA_Only, "Location", "DPI_4LT", "Tukey")
#CA_4DPI_LSD_LT = Stats(CA_Only, "Location", "DPI_4LT", "Fisher")

# Square-root transformed data
#QC_Stats(Pear_Data, "DPI_4SQRT")
#DPI4_TT_SQRT = Stats(Pear_Data, "Isolate", "DPI_4SQRT", "Tukey")
#DPI4_LSD_SQRT = Stats(Pear_Data, "Isolate", "DPI_4SQRT", "Fisher")
#DPI4_Dunn_SQRT = NP_Stats(Pear_Data, "Isolate", "DPI_4SQRT")
# North v South 
#CA_4DPI_TT_SQRT = Stats(CA_Only, "Location", "DPI_4SQRT", "Tukey")
#CA_4DPI_LSD_SQRT = Stats(CA_Only, "Location", "DPI_4SQRT", "Fisher")

#------------------
# 5 DPI (Combined)
#------------------
QC_Stats(Pear_Data, "DPI_5")
DPI5_TT = Stats(Pear_Data, "Isolate", "DPI_5", "Tukey")
DPI5_LSD = Stats(Pear_Data, "Isolate", "DPI_5", "Fisher")
DPI5_Dunn = NP_Stats(Pear_Data, "Isolate", "DPI_5")

# North v South 
CA_5DPI_TT = Stats(CA_Only, "Location", "DPI_5", "Tukey")
CA_5DPI_LSD = Stats(CA_Only, "Location", "DPI_5", "Fisher")

# Log-transformed data
#QC_Stats(Pear_Data, "DPI_5LT")
#DPI5_TT_LT = Stats(Pear_Data, "Isolate", "DPI_5LT", "Tukey")
#DPI5_LSD_LT = Stats(Pear_Data, "Isolate", "DPI_5LT", "Fisher")
#DPI5_Dunn_LT = NP_Stats(Pear_Data, "Isolate", "DPI_5LT")
# North v South 
#CA_5DPI_TT_LT = Stats(CA_Only, "Location", "DPI_5LT", "Tukey")
#CA_5DPI_LSD_LT = Stats(CA_Only, "Location", "DPI_5LT", "Fisher")
#--------------------
# Lesion Area Summary
#--------------------
# 3DPI
SUM3DPI <- as.data.frame(Pear_Data %>% group_by(Isolate, Location) %>% 
                         dplyr::summarize(DPI3_Mean = mean(DPI_3, na.rm=TRUE), DPI3_SD = sd(DPI_3, na.rm=TRUE)))
# 4DPI
SUM4DPI <- as.data.frame(Pear_Data %>% group_by(Isolate, Location) %>% 
                           dplyr::summarize(DPI4_Mean = mean(DPI_4, na.rm=TRUE), DPI4_SD = sd(DPI_4, na.rm=TRUE)))
# 5DPI
SUM5DPI <- as.data.frame(Pear_Data %>% group_by(Isolate, Location) %>% 
                           dplyr::summarize(DPI5_Mean = mean(DPI_5, na.rm=TRUE), DPI5_SD = sd(DPI_5, na.rm=TRUE)))

# rounding the decimals to 2 decimals
SUM3DPI = SUM3DPI %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUM4DPI = SUM4DPI %>% mutate(across(where(is.numeric), ~ round(., 2)))
SUM5DPI = SUM5DPI %>% mutate(across(where(is.numeric), ~ round(., 2)))

# Subsetting data
SUM4DPI_j = subset(SUM4DPI, select = c(Isolate, DPI4_Mean, DPI4_SD))
SUM5DPI_j = subset(SUM5DPI, select = c(Isolate, DPI5_Mean, DPI5_SD))
# 3dpi
DPI3_TT_j = subset(DPI3_TT, select = c(Isolate, Groups))
colnames(DPI3_TT_j) <- c('Isolate', 'TT_DPI3_Groups')
DPI3_FLSD_j = subset(DPI3_LSD, select = c(Isolate, Groups))
colnames(DPI3_FLSD_j) <- c('Isolate', 'FLSD_DPI3_Groups')
DPI3_Dunn_j = subset(DPI3_Dunn, select = c(Isolate, Dunn_Group))
colnames(DPI3_Dunn_j) <- c('Isolate', 'DUNN_DPI3_Groups')
# 4DPI
DPI4_TT_j = subset(DPI4_TT, select = c(Isolate, Groups))
colnames(DPI4_TT_j) <- c('Isolate', 'TT_DPI4_Groups')
DPI4_FLSD_j = subset(DPI4_LSD, select = c(Isolate, Groups))
colnames(DPI4_FLSD_j) <- c('Isolate', 'FLSD_DPI4_Groups')
DPI4_Dunn_j = subset(DPI4_Dunn, select = c(Isolate, Dunn_Group))
colnames(DPI4_Dunn_j) <- c('Isolate', 'DUNN_DPI4_Groups')
# 5DPI
DPI5_TT_j = subset(DPI5_TT, select = c(Isolate, Groups))
colnames(DPI5_TT_j) <- c('Isolate', 'TT_DPI5_Groups')
DPI5_FLSD_j = subset(DPI5_LSD, select = c(Isolate, Groups))
colnames(DPI5_FLSD_j) <- c('Isolate', 'FLSD_DPI5_Groups')
DPI5_Dunn_j = subset(DPI5_Dunn, select = c(Isolate, Dunn_Group))
colnames(DPI5_Dunn_j) <- c('Isolate', 'DUNN_DPI5_Groups')

# Create list of dfs
Sum_List = list(SUM3DPI, DPI3_TT_j, DPI3_FLSD_j, DPI3_Dunn_j, SUM4DPI_j, DPI4_TT_j,
                DPI4_FLSD_j, DPI4_Dunn_j, SUM5DPI_j, DPI5_TT_j, DPI5_FLSD_j, DPI5_Dunn_j)
# This data is used in Table S4
Comb_Summary = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List) 
#write.csv(Comb_Summary,"Pear_Lesion_Area_Summarized.csv")

### ----------------------------------------------------------------------------
# AUDPC Analysis
### ----------------------------------------------------------------------------
Pear_AUDPC = Calculate_AUDPC(Pear_Data)
#write.csv(Pear_AUDPC, "Pear_AUDPC_Results.csv")
head(Pear_AUDPC)
QC_Stats(Pear_AUDPC, "AUDPC_Values")
AUDPC_TT = Stats(Pear_AUDPC, "Isolate", "AUDPC_Values", "Tukey")
AUDPC_LSD = Stats(Pear_AUDPC, "Isolate", "AUDPC_Values", "Fisher")
AUDPC_DUNN = NP_Stats(Pear_AUDPC, "Isolate", "AUDPC_Values")
# Subset CA 
Pear_AUDPC_CA = subset(Pear_AUDPC, Region=="CA")

# Log transformed AUDPC values 
Pear_AUDPC$AUDPC_Values_LT = log(Pear_AUDPC$AUDPC_Values)
QC_Stats(Pear_AUDPC, "AUDPC_Values_LT")
AUDPC_TT_LT = Stats(Pear_AUDPC, "Isolate", "AUDPC_Values_LT", "Tukey")
AUDPC_LSD_LT = Stats(Pear_AUDPC, "Isolate", "AUDPC_Values_LT", "Fisher")
AUDPC_DUNN_LT = NP_Stats(Pear_AUDPC, "Isolate", "AUDPC_Values_LT")

# Square-root transformed AUDPC values
Pear_AUDPC$AUDPC_Values_SQRT = sqrt(Pear_AUDPC$AUDPC_Values)
QC_Stats(Pear_AUDPC, "AUDPC_Values_SQRT")
AUDPC_TT_SQRT = Stats(Pear_AUDPC, "Isolate", "AUDPC_Values_SQRT", "Tukey")
AUDPC_LSD_SQRT = Stats(Pear_AUDPC, "Isolate", "AUDPC_Values_SQRT", "Fisher")
AUDPC_DUNN_SQRT = NP_Stats(Pear_AUDPC, "Isolate", "AUDPC_Values_SQRT")


# Looking at Difference between CA-North & CA-South
QC_Stats(Pear_AUDPC_CA, "AUDPC_Values")
AUDPC_TT_CA = Stats(Pear_AUDPC_CA, "Location", "AUDPC_Values", "Tukey")
AUDPC_LSD_CA = Stats(Pear_AUDPC_CA, "Location", "AUDPC_Values", "Fisher")


colnames(AUDPC_TT_CA)[1] <- "Location"
colnames(AUDPC_LSD_CA)[1] <- "Location"

AUDPC_dfs = list(AUDPC_TT, AUDPC_LSD, AUDPC_DUNN)

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
AUDPC_DUNN = modified_AUDPC_dfs[[3]]

#--------------
# AUDPC Summary
#--------------
AUDPC_SUM <- as.data.frame(Pear_AUDPC %>% group_by(Isolate, Location) %>% 
                           dplyr::summarize(AUDPC_Mean = mean(AUDPC_Values, na.rm=TRUE), AUDPC_SD = sd(AUDPC_Values, na.rm=TRUE)))
AUDPC_SUM = AUDPC_SUM %>% mutate(across(where(is.numeric), ~ round(., 2)))

AUDPC_TT_j = subset(AUDPC_TT, select = c(Isolate, Groups))
colnames(AUDPC_TT_j) <- c('Isolate', 'TT_Groups')
AUDPC_FLSD_j = subset(AUDPC_LSD, select = c(Isolate, Groups))
colnames(AUDPC_FLSD_j) <- c('Isolate', 'FLSD_Groups')
# Create list of dfs
Sum_List2 = list(AUDPC_SUM, AUDPC_TT_j, AUDPC_FLSD_j)
# This data is used in Table S4
Comb_Summary2 = Reduce(function(x, y) merge(x, y, all=TRUE), Sum_List2) 
#write.csv(Comb_Summary2,"Pear_AUDPC_Summary.csv")

# Location Summary 
AUDPC_Loc_SUM <- as.data.frame(Pear_AUDPC %>% group_by(Location) %>% 
                             dplyr::summarize(AUDPC_Mean = mean(AUDPC_Values, na.rm=TRUE), AUDPC_SD = sd(AUDPC_Values, na.rm=TRUE)))
AUDPC_Loc_SUM = AUDPC_Loc_SUM %>% mutate(across(where(is.numeric), ~ round(., 2)))
#write.csv(AUDPC_Loc_SUM,"Pear_Location_AUDPC_Summary.csv")
#-------------------------------------
# Visualizing AUDPC with Combined Data 
#-------------------------------------
# All isolates with Dunn's Groups 
BP4 = ggplot(Pear_AUDPC, aes(x = interaction(Isolate, Location), y=AUDPC_Values, fill=as.factor(Location))) + 
  geom_boxplot() +
  geom_text(data = AUDPC_DUNN, 
            aes(label = Dunn_Group, y = max(Pear_AUDPC$AUDPC_Values) + 0.5),
            vjust = -0.5) +
  scale_fill_manual(values=LocationCols, name="Location") +
  scale_x_discrete(labels = function(x) sub(".CA-North|.CA-South|.Mexico", "", x))

BP4 = BP4 + theme_classic() + 
  labs(x="Isolate", y = "AUDPC")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12),
        legend.position = "none") + coord_cartesian(clip = "off")
BP4
#ggsave(BP4,filename = "Pear_AUDPC_Combined.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(BP4,filename = "Pear_AUDPC_Combined.Dunn.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

# CA Boxplots with Fisher's LSD but the Kruskal-Wallis indicates significance too and the group labels would be the same
BP5 = ggplot(Pear_AUDPC_CA, aes(x = Location, y=AUDPC_Values, fill=as.factor(Location))) + 
  geom_boxplot() +
  geom_text(data = AUDPC_LSD_CA, 
            aes(label = Groups, y = max(Pear_AUDPC_CA$AUDPC_Values) + 0.5),
            vjust = -0.5) +
  scale_fill_manual(values=LocationCols, name="Location") 

BP5 = BP5 + theme_classic() + 
  labs(x="Isolate", y = "AUDPC")+
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title=element_text(size=12),
        legend.position = "none") + coord_cartesian(clip = "off")
#BP5 = BP5+labs(title=expression(paste("AUDPC of California ", italic("P. cinnamomi "), "infecting pear (Fisher's LSD groups)")))
BP5
#ggsave(BP5,filename = "Pear_AUDPC_Combined_CA_North_v_South.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(BP5,filename = "Pear_AUDPC_Combined_CA_North_v_South.FLSD.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)




### ----------------------------------------------------------------------------
# Lesion Area Visualization with Fisher's LSD 
### ----------------------------------------------------------------------------
#-------------------------
# Visualizing Lesion Area
#   and Formatting DFs
#-------------------------

# Visualizing Lesion area 
# adding HPI column used for geom_text
DPI3_LSD$HPI = "72"
DPI4_LSD$HPI = "96"
DPI5_LSD$HPI = "120"

# merging the std data from Pear_Sum_Data to the HPIXXLSD dataframe for geom_text
DPI3_LSD <- merge(DPI3_LSD, Pear_Sum_Data[Pear_Sum_Data$HPI == 72, c("Isolate", "std")], by = "Isolate", all.x = TRUE)
DPI4_LSD <- merge(DPI4_LSD, Pear_Sum_Data[Pear_Sum_Data$HPI == 96, c("Isolate", "std")], by = "Isolate", all.x = TRUE)
DPI5_LSD <- merge(DPI5_LSD, Pear_Sum_Data[Pear_Sum_Data$HPI == 120, c("Isolate", "std")], by = "Isolate", all.x = TRUE)

# Southern Isolates
South = c("P336", "P467", "P482", "Pc2109", "Pc2120", "Pc5241")
North = c("P344", "P444", "Pc2110", "Pc2113", "Pc2114")
Mexican = c("SER5_6", "TANR3_5", "TANS5_22")

dfs = list(DPI3_LSD, DPI4_LSD, DPI5_LSD)

modified_dfs <- list()

for (df in dfs) {
  df$Location = ifelse(df$Isolate %in% South, "CA-South",
                        ifelse(df$Isolate %in% North, "CA-North",
                               ifelse(df$Isolate %in% Mexican, "Mexico", NA)))
  
  df$Region = ifelse(df$Isolate %in% South, "CA",
                      ifelse(df$Isolate %in% North, "CA",
                             ifelse(df$Isolate %in% Mexican, "MX", NA)))
  
  # Append the modified dataframe to the list
  modified_dfs[[length(modified_dfs) + 1]] = df
}

DPI3_LSD = modified_dfs[[1]]
DPI4_LSD = modified_dfs[[2]]
DPI5_LSD = modified_dfs[[3]]


# 72HPI 
Bar72 <- ggplot(DPI3_LSD, aes(x = interaction(Isolate, Location), y=Mean, fill=as.factor(Location), group=as.factor(Location))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=Mean-std, ymax=Mean+std), width=.2, position=position_dodge(0.9))+ 
  ylab(bquote('Lesion area '(cm^2))) + 
  xlab("Isolate") +
  geom_text(data = DPI3_LSD, aes(label = Groups, y = max(Mean) + 1), vjust = -.3)+
  scale_fill_manual(values=LocationCols, name="Location") + 
  #labs(title=expression(paste("Average lesion size of ", italic("P. cinnamomi "), "infecting pear at 72HPI (Fisher's LSD Groups)")))+
  theme_classic() +
  scale_x_discrete(labels = function(x) sub(".CA-North|.CA-South|.Mexico", "", x))

Bar72 = Bar72 + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                      axis.text.y = element_text(size = 10),
                      axis.title=element_text(size=12,face="bold"),
                      legend.position = "none") + coord_cartesian(clip = "off")
Bar72
#ggsave(Bar72,filename = "Pear_72HPI_Bars.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(Bar72,filename = "Pear_72HPI_Bars.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

# 96HPI 
Bar96 <- ggplot(DPI4_LSD, aes(x = interaction(Isolate, Location), y=Mean, fill=as.factor(Location), group=as.factor(Location))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=Mean-std, ymax=Mean+std), width=.2, position=position_dodge(0.9))+ 
  ylab(bquote('Lesion area '(cm^2))) + 
  xlab("Isolate") +
  geom_text(data = DPI4_LSD, aes(label = Groups, y = max(Mean) + 4.2), vjust = -1)+
  scale_fill_manual(values=LocationCols, name="Location") + 
  #labs(title=expression(paste("Average lesion size of ", italic("P. cinnamomi "), "infecting pear at 96HPI (Fisher's LSD Groups)")))+
  theme_classic() +
  scale_x_discrete(labels = function(x) sub(".CA-North|.CA-South|.Mexico", "", x))

Bar96 = Bar96 + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                      axis.text.y = element_text(size = 10),
                      axis.title=element_text(size=12,face="bold"),
                      legend.position = "none") + coord_cartesian(clip = "off")
Bar96
#ggsave(Bar96,filename = "Pear_96HPI_Bars.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(Bar96,filename = "Pear_96HPI_Bars.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

# 120HPI 
Bar120 <- ggplot(DPI5_LSD, aes(x = interaction(Isolate, Location), y=Mean, fill=as.factor(Location), group=as.factor(Location))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=Mean-std, ymax=Mean+std), width=.2, position=position_dodge(0.9))+ 
  ylab(bquote('Lesion area '(cm^2))) + 
  xlab("Isolate") +
  geom_text(data = DPI5_LSD, aes(label = Groups, y = max(Mean) + 18), vjust = -1)+
  scale_fill_manual(values=LocationCols, name="Location") + 
  #labs(title=expression(paste("Average lesion size of ", italic("P. cinnamomi "), "infecting pear at 120HPI (Fisher's LSD Groups)")))+
  theme_classic() +
  scale_x_discrete(labels = function(x) sub(".CA-North|.CA-South|.Mexico", "", x))

Bar120 = Bar120 + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                        axis.text.y = element_text(size = 10),
                        axis.title=element_text(size=12,face="bold"),
                        legend.position = "none") + coord_cartesian(clip = "off")
Bar120

#ggsave(Bar120,filename = "Pear_120HPI_Bars.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(Bar120,filename = "Pear_120HPI_Bars.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

### ----------------------------------------------------------------------------
# Lesion Area Visualization with Dunn's Groups
### ----------------------------------------------------------------------------
#-------------------------
# Visualizing Lesion Area
#   and Formatting DFs
#-------------------------

# Visualizing Lesion area 
# adding HPI column used for geom_text
DPI3_Dunn$HPI = "72"
DPI4_Dunn$HPI = "96"
DPI5_Dunn$HPI = "120"

# merging the std data from Pear_Sum_Data to the HPIXXLSD dataframe for geom_text
DPI3_Dunn <- merge(DPI3_Dunn, Pear_Sum_Data[Pear_Sum_Data$HPI == 72, c("Isolate", "mean", "std")], by = "Isolate", all.x = TRUE)
DPI4_Dunn <- merge(DPI4_Dunn, Pear_Sum_Data[Pear_Sum_Data$HPI == 96, c("Isolate", "mean", "std")], by = "Isolate", all.x = TRUE)
DPI5_Dunn <- merge(DPI5_Dunn, Pear_Sum_Data[Pear_Sum_Data$HPI == 120, c("Isolate", "mean", "std")], by = "Isolate", all.x = TRUE)

# Southern Isolates
South = c("P336", "P467", "P482", "Pc2109", "Pc2120", "Pc5241")
North = c("P344", "P444", "Pc2110", "Pc2113", "Pc2114")
Mexican = c("SER5_6", "TANR3_5", "TANS5_22")

dfs = list(DPI3_Dunn, DPI4_Dunn, DPI5_Dunn)

modified_dfs <- list()

for (df in dfs) {
  df$Location = ifelse(df$Isolate %in% South, "CA-South",
                       ifelse(df$Isolate %in% North, "CA-North",
                              ifelse(df$Isolate %in% Mexican, "Mexico", NA)))
  
  df$Region = ifelse(df$Isolate %in% South, "CA",
                     ifelse(df$Isolate %in% North, "CA",
                            ifelse(df$Isolate %in% Mexican, "MX", NA)))
  
  # Append the modified dataframe to the list
  modified_dfs[[length(modified_dfs) + 1]] = df
}

DPI3_Dunn = modified_dfs[[1]]
DPI4_Dunn = modified_dfs[[2]]
DPI5_Dunn = modified_dfs[[3]]


# 72HPI 
Bar72D <- ggplot(DPI3_Dunn, aes(x = interaction(Isolate, Location), y=mean, fill=as.factor(Location), group=as.factor(Location))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2, position=position_dodge(0.9))+ 
  ylab(bquote('Lesion area '(cm^2))) + 
  xlab("Isolate") +
  geom_text(data = DPI3_Dunn, aes(label = Dunn_Group, y = max(mean) + 1), vjust = -.3)+
  scale_fill_manual(values=LocationCols, name="Location") + 
  #labs(title=expression(paste("Average lesion size of ", italic("P. cinnamomi "), "infecting pear at 72HPI (Fisher's LSD Groups)")))+
  theme_classic() +
  scale_x_discrete(labels = function(x) sub(".CA-North|.CA-South|.Mexico", "", x))

Bar72D = Bar72D + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                      axis.text.y = element_text(size = 10),
                      axis.title=element_text(size=12),
                      legend.position = "none") + coord_cartesian(clip = "off")
Bar72D
#ggsave(Bar72,filename = "Pear_72HPI_Bars.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(Bar72D,filename = "Pear_72HPI_Bars.Dunn.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

# 96HPI 
Bar96D <- ggplot(DPI4_Dunn, aes(x = interaction(Isolate, Location), y=mean, fill=as.factor(Location), group=as.factor(Location))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2, position=position_dodge(0.9))+ 
  ylab(bquote('Lesion area '(cm^2))) + 
  xlab("Isolate") +
  geom_text(data = DPI4_Dunn, aes(label = Dunn_Group, y = max(mean) + 4.2), vjust = -1)+
  scale_fill_manual(values=LocationCols, name="Location") + 
  #labs(title=expression(paste("Average lesion size of ", italic("P. cinnamomi "), "infecting pear at 96HPI (Fisher's LSD Groups)")))+
  theme_classic() +
  scale_x_discrete(labels = function(x) sub(".CA-North|.CA-South|.Mexico", "", x))

Bar96D = Bar96D + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                      axis.text.y = element_text(size = 10),
                      axis.title=element_text(size=12),
                      legend.position = "none") + coord_cartesian(clip = "off")
Bar96D
#ggsave(Bar96,filename = "Pear_96HPI_Bars.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(Bar96D,filename = "Pear_96HPI_Bars.Dunn.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

# 120HPI 
Bar120D <- ggplot(DPI5_Dunn, aes(x = interaction(Isolate, Location), y=mean, fill=as.factor(Location), group=as.factor(Location))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2, position=position_dodge(0.9))+ 
  ylab(bquote('Lesion area '(cm^2))) + 
  xlab("Isolate") +
  geom_text(data = DPI5_Dunn, aes(label = Dunn_Group, y = max(mean) + 18), vjust = -1)+
  scale_fill_manual(values=LocationCols, name="Location") + 
  #labs(title=expression(paste("Average lesion size of ", italic("P. cinnamomi "), "infecting pear at 120HPI (Fisher's LSD Groups)")))+
  theme_classic() +
  scale_x_discrete(labels = function(x) sub(".CA-North|.CA-South|.Mexico", "", x))

Bar120D = Bar120D + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                        axis.text.y = element_text(size = 10),
                        axis.title=element_text(size=12),
                        legend.position = "none") + coord_cartesian(clip = "off")
Bar120D

#ggsave(Bar120,filename = "Pear_120HPI_Bars.pdf", width = 8, height = 6, units = "in", limitsize = FALSE)
ggsave(Bar120D,filename = "Pear_120HPI_Bars.Dunn.v2.pdf", width = 6, height = 4, units = "in", limitsize = FALSE)

