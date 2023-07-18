#-------------------------------------------------------------------------------
# P. cinnamomi Population Genomic Analysis
# Aidan Shands 
# Manosalva Lab @ University of California, Riverside
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Load Packages
#-------------------------------------------------------------------------------
library(adegenet)
library(poppr)
library(vcfR)
library(ape)
library(ggtree)
library(phytools)
library(phyloseq)
library(ggplot2)
library(treeio)
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
library(RColorBrewer)
library(treeio)
library(ape)
library(tidytree)
library(TDbook)
library(dplyr)
library(aplot)
library(tidyr)
library(igraph)
library(dartR)
library(SNPRelate)
library(terra)
library(reshape2)
library(ggpubr)
library(ggstance)
library(mmod)
library(gplots)
library(HardyWeinberg)
library(ggtern)
library(colorspace)
library(ggrepel)
library(grid)
library(gridExtra)
library(gtable)
library(gridGraphics)
library(scatterpie)
library(jsonlite)
library(googleway)
library(ggmapstyles)


#-------------------------------------------------------------------------------
# Importing Annotation Data
#-------------------------------------------------------------------------------
# Annotations 
annot = read.csv('annotations/Annotations.csv', header = T)
# MLG 0.01 for 92 individuals 
annot_MLG0.01 = read.csv('annotations/Pc92_583k_MLG0.01_CC_Annot.csv', header = T)
# MLG annotations for all 136 isolates 
annot_MLG0.01_All = read.csv('annotations/MLG0.01_Isolates.csv', header = T)
# Importing FastStructure Log Prior where k=4 ()
FastStructurek4= read.csv("annotations/FastStructure_LP_K4.csv", header =T)
FastStructurek4 = as.data.frame(FastStructurek4)

#-------------------------------------------------------------------------------
#        Generating Genlight Object and Formatting
#-------------------------------------------------------------------------------
vcf_file_OG = 'Pc_only136.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.recode.vcf.gz'
vcf_OG<- read.vcfR(vcf_file_OG, verbose = FALSE)
vcf_OG

# Export path 
F583k_Path = '583K_outputs'

# creating Genlight and GenInd objects 
vcf_GL_OG<- vcfR2genlight(vcf_OG)
# Defining ploidy, populations (by geographical location), and other (name)
ploidy(vcf_GL_OG) <- 2
vcf_GL_OG$pop <- as.factor(annot$location_Simple)
vcf_GL_OG$other <- as.list(annot$ID)
# Setting colors
cols <- brewer.pal(n = nPop(vcf_GL_OG), name = "Dark2")
LocationCols <- c("Asia/Oceania" = "#1B9E77", "MX" = "#E7298A", 
                  "CA-south" = "#7570B3", "CA-north" = "#D95F02")
#-------------------------------------------------------------------------------
#             Principal Component Analysis (PCA)
#                                Figure S1A
#-------------------------------------------------------------------------------

# Determining the amount of PCs to use
# all pops
vcf_GL.pca <- glPca(vcf_GL_OG, nf = 4)
barplot(100*vcf_GL.pca$eig/sum(vcf_GL.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
#quartz.save("Pc136_583k_Eigen_Bars.pdf", type="pdf")

# Variance explained in the 4 PCAs = 90.36%
100 * round(sum(vcf_GL.pca$eig[1:4]/sum(vcf_GL.pca$eig)), digits = 5)

# with all pops
vcf_GL.pca.scores <- as.data.frame(vcf_GL.pca$scores)
vcf_GL.pca.scores$pop <- pop(vcf_GL_OG)

# get the variance of each axis
vcf_GL.pca.eigen = data.frame(vcf_GL.pca$eig)
vcf_GL.pca.eigen$Percent_Variance= 100*vcf_GL.pca.eigen/sum(vcf_GL.pca.eigen)
vcf_GL.pca.eigen # take the first 2 and they reflect my eigenbars

PC1_Variance_Percent_Label = '81.2%'
PC2_Variance_Percent_Label = '6.06%'

set.seed(9)
dev.off()
# plotting 
p <- ggplot(vcf_GL.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=0.5, position = position_jitter(width = 0.1, 
                                                         height = 0.1))
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw() + ggtitle("PCA on Pc SNP data (583,627 SNPs)") + 
  xlab("PC1 (81.2%)") + ylab("PC2 (6.06%)")
p

#ggsave("Pc136_583k_PCA.FINAL.pdf", plot = p, path=F583k_Path, width = 5, height = 5, units = "in")

vcf_GLpca_Scores = vcf_GL.pca$scores
#write.csv(vcf_GLpca_Scores, "583k_outputs/PCA_Scores_583k_136Pc.csv")

#-------------------------------------------------------------------------------
#               Discriminant Analysis of Principal Components
#                                Figure S2C
#-------------------------------------------------------------------------------
# Generating DAPC
DAPC <- dapc(vcf_GL_OG, pop=as.factor(unlist(vcf_GL_OG$pop)), n.pca=10, n.da=3,
             var.loadings=T, pca.info=T)
# DAPC Scatterplot 
scatter.dapc(DAPC, grp=as.factor(unlist(vcf_GL_OG$pop)), legend=T, 
             col=cols, cex.lab=0.25, clabel = 0, cleg = 1.2,
             cellipse = 1.5, ratio.da=.15, posi.leg=c(11,11),
             posi.da = "topright")
#quartz.save("Pc136_583k_DAPC_Scatter.FINAL.v2.pdf", type="pdf")

#-------------------------------------------------------------------------------
#     Unweighted Pair Group Method with Arithmetic Mean (UPGMA) Tree
#                         Tree Used in Figure 1
#-------------------------------------------------------------------------------

# All pop
tree <- aboot(vcf_GL_OG, tree = "upgma", distance = bitwise.dist, 
              sample = 1000, showtree = F, cutoff = 50, quiet = T)
plot.phylo(tree, cex = 0.28, font = 2, adj = 0, 
           tip.color =  cols[pop(vcf_GL_OG)])

nodelabels(tree$node.label, adj = c(1.65, -0.5), 
           frame = "n", cex = 0.3,font = 3, xpd = TRUE)

legend('topleft', legend = c("Asia","CA-north","CA-south", "MX"), 
       fill = cols, border = FALSE, bty = "n", cex = 0.5)

axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
#quartz.save("Pc136_583k_UPGMA_Tree_1000BSR.pdf", type="pdf")
#write.nexus(tree, file = "Pc136_583k_UPGMA_Tree_1000BSR.nex") 
#-------------------------------------------------------------------------------
#                         Minimum spanning Network (MSN)
#                                Figure S2B
#-------------------------------------------------------------------------------
# customize in Shiny app
# imsn()

vcf_GL_OG_sub <- popsub(vcf_GL_OG, exclude = character(0))
vcf_GL_OG_nomiss <- missingno(vcf_GL_OG, type = 'mean')
vcf_GL_OG_dist <- bitwise.dist(vcf_GL_OG_nomiss, percent = TRUE, 
                               mat = FALSE, missing_match = TRUE, 
                               scale_missing = FALSE, euclidean = FALSE, 
                               differences_only = FALSE, threads = 0)

min_span_net <- poppr.msn(vcf_GL_OG_sub, vcf_GL_OG_dist, 
                          showplot = FALSE, include.ties = TRUE)

set.seed(69)
plot_poppr_msn(vcf_GL_OG,
               min_span_net,
               inds = "NONE",
               mlg = FALSE,
               gadj = 6,
               nodescale = 6,
               palette = cols,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_nicely)

#-------------------------------------------------------------------------------
#         Export Genlight to STRUCTURE Format for FastStructure
#                   FastStructure Data Used in Figure 1
#-------------------------------------------------------------------------------
gl = gl.compliance.check(vcf_GL_OG)
gl2structure(vcf_GL_OG ,ploidy = 2, 
             outfile = "Pc_only136.p2.f1SV.vcft_dp4_Q20_NoSingleton_noMD_NoRepeats.583k.STRUCTURE.txt",
             outpath = '/Users/manosalvalab/Desktop/Aidan/Pc_PopGen_2022/R_Analysis_Final/583K_outputs')

#-------------------------------------------------------------------------------
#               Clone Correction with MLG Filter of 0.1    
#                     MLGs used in Tables S2 & S4
#-------------------------------------------------------------------------------
# convert to snpclone
vcf_GL_OG_SC2 = as.snpclone(vcf_GL_OG)
# set pop by location
pop(vcf_GL_OG_SC2) <- as.factor(annot$location_Simple)
# set strata 
strata(vcf_GL_OG_SC2) = annot 
# MLG filter of 0.01
mlg.filter(vcf_GL_OG_SC2, 
           algorithm = "average_neighbor", 
           distance = "bitwise.dist") <- 0.01
# Clone correct 
vcf_GL_OG_SC2_CC = clonecorrect(vcf_GL_OG_SC2, strata= ~location_Simple)

# Export MLGs and Isolate names 
vcf_GL_OG_SC2_Vec = mlg.vector(vcf_GL_OG_SC2)
names(vcf_GL_OG_SC2_Vec) <- vcf_GL_OG$ind.names
#write.csv(vcf_GL_OG_SC2_Vec, "vcf_GL_OG_SC2_MLG0.01_Isolates.csv")

#-------------------------------------------------------------------------------
#                 UPGMA Tree with Clone-Corrected Data
#-------------------------------------------------------------------------------

# Tree 
tree <- aboot(vcf_GL_OG_SC2_CC, tree = "upgma", distance = bitwise.dist, 
              sample = 1000, showtree = F, cutoff = 50, quiet = T)
plot.phylo(tree, cex = 0.28, font = 2, adj = 0, 
           tip.color =  cols[pop(vcf_GL_OG_SC2_CC)])
nodelabels(tree$node.label, adj = c(1.65, -0.5), frame = "n", 
           cex = 0.3,font = 3, xpd = TRUE)
legend('topleft', legend = c("MX","CA-north","Asia/Oceania", "CA-South"), 
       fill = cols, border = FALSE, bty = "n", cex = 0.5)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
#quartz.save("Pc136_583k_UPGMA_Tree_0.01MLG_CC_1000BS.pdf", type="pdf")
#write.nexus(tree, file = "Pc92_583k_UPGMA_Tree_0.01MLG_CC_1000BS.nex") 
#-------------------------------------------------------------------------------
#                      MSN with Clone-Corrected Data
#                               Figure S2C
#-------------------------------------------------------------------------------
vcf_GL_OG_SC2_CC_sub <- popsub(vcf_GL_OG_SC2_CC, exclude = character(0))
vcf_GL_OG_SC2_CC_nomiss <- missingno(vcf_GL_OG_SC2_CC, type = 'mean')
vcf_GL_OG_SC2_CC_dist <- bitwise.dist(vcf_GL_OG_SC2_CC_nomiss, percent = TRUE, 
                                      mat = FALSE, missing_match = TRUE, 
                                      scale_missing = FALSE, euclidean = FALSE, 
                                      differences_only = FALSE, threads = 0)

min_span_net <- poppr.msn(vcf_GL_OG_SC2_CC_sub, vcf_GL_OG_SC2_CC_dist, 
                          showplot = FALSE, include.ties = TRUE)

set.seed(69)
plot_poppr_msn(vcf_GL_OG_SC2_CC,
               min_span_net,
               inds = "NONE",
               mlg = FALSE,
               gadj = 6,
               nodescale = 25,
               palette = LocationCols,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_nicely)

#-------------------------------------------------------------------------------
#                       PCA with Clone-Corrected Data 
#                                Figure S2A
#-------------------------------------------------------------------------------

# Determining the amount of PCs to use
vcf_GL_CC.pca <- glPca(vcf_GL_OG_SC2_CC, nf = 4)
barplot(100*vcf_GL_CC.pca$eig/sum(vcf_GL_CC.pca$eig), 
        col = heat.colors(50), main="PCA Eigenvalues")

title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
#quartz.save("Pc136_583k_0.1MLGCC_Eigen_Bars.pdf", type="pdf")
# % variance explained by the 4 PCs
100 * round(sum(vcf_GL_CC.pca$eig[1:4]/sum(vcf_GL_CC.pca$eig)), digits = 5)

# with all pops
vcf_GL_CC.pca.scores <- as.data.frame(vcf_GL_CC.pca$scores)
vcf_GL_CC.pca.scores$pop <- pop(vcf_GL_OG_SC2_CC)

# get the variance of each axis
vcf_GL_CC.pca.eigen = data.frame(vcf_GL_CC.pca$eig)
vcf_GL_CC.pca.eigen$Percent_Variance= 100*vcf_GL_CC.pca.eigen/sum(vcf_GL_CC.pca.eigen)
vcf_GL_CC.pca.eigen # take the first 2 and they reflect my eigenbars

PC1_Variance_Percent_Label = '68.2%'
PC2_Variance_Percent_Label = '11.1%'

set.seed(9)
dev.off()
# plotting 
p <- ggplot(vcf_GL_CC.pca.scores, aes(x=PC1, y=PC2, color=pop)) 
p <- p + geom_point(size=0.5, position = position_jitter(width = 0.1, height = 0.1))
p <- p + scale_color_manual(values = LocationCols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw() + ggtitle("PCA on Pc SNP data 0.1CC (583,627 SNPs)") + 
  xlab("PC1 (68.2%)") + ylab("PC2 (11.1%)")
p

#ggsave("Pc136_583k_PCA_0.1MLGCC.FINAL.pdf", plot = p, path=F583k_Path, width = 5, height = 5, units = "in")
#write.csv(vcf_GL_CC.pca.scores, "583k_outputs/PCA_Scores_583k_136Pc_0.1MLGCC.csv")

#-------------------------------------------------------------------------------
#                       DAPC with Clone-Corrected Data 
#                                Figure S2B
#-------------------------------------------------------------------------------

# identify number of PCs to use for n.pca
#pramx_F3 <- xvalDapc(tab(vcf_GL3, NA.method = "mean"), pop(vcf_GL3))
# Generating DAPC
DAPC <- dapc(vcf_GL_OG_SC2_CC, pop=as.factor(unlist(vcf_GL_OG_SC2_CC$pop)), 
             n.pca=10, n.da=3, var.loadings=T, pca.info=T)
group_levels <- levels(DAPC$grp)
group_colors <- LocationCols[group_levels]
# DAPC Scatterplot 
scatter.dapc(DAPC, grp=as.factor(unlist(vcf_GL_OG_SC2_CC$pop)), legend=T, 
             col=group_colors, cex.lab=0.25, clabel = 0, cleg = 1.2,
             cellipse = 1.5, ratio.da=.15, posi.leg="topleft",
             posi.da = "topright")
#quartz.save("Pc136_583k_0.1CC_DAPC_Scatter.FINAL.v2.pdf", type="pdf")

#-------------------------------------------------------------------------------
#      Subset by Population for Uncorrected and Clone-Corrected data
#-------------------------------------------------------------------------------
# Uncorrected
# CA-North 
CAnorth = popsub(vcf_GL_OG, exclude=c('Asia','MX', 'CA-south'))
# Mexico 
MX = popsub(vcf_GL_OG, exclude=c('Asia','CA-north', 'CA-south'))
# CA South 
CAsouth = popsub(vcf_GL_OG, exclude=c('Asia','CA-north', 'MX'))
# Asia/Oceania
Asia = popsub(vcf_GL_OG, exclude=c('CA-south','CA-north', 'MX'))

# Clone-Corrected 
# CA-North 
CAnorthMLG0.01 = popsub(vcf_GL_OG_SC2_CC, exclude=c('Asia','MX', 'CA-south'))
# Mexico 
MX_MLG0.01 = popsub(vcf_GL_OG_SC2_CC, exclude=c('Asia','CA-north', 'CA-south'))
# CA South 
CAsouthMLG0.01 = popsub(vcf_GL_OG_SC2_CC, exclude=c('Asia','CA-north', 'MX'))
# Asua
AsiaMLG0.01 = popsub(vcf_GL_OG_SC2_CC, exclude=c('CA-south','CA-north', 'MX'))

#-------------------------------------------------------------------------------
# Functions to perform IA, Basic Stats, & Pairwise GST (Hedrick's & Nei's)
#-------------------------------------------------------------------------------
# These were all performed on the UCR HPCC because they take a long time. 

# Calculate Index of Association (IA) with window size of 500 (Poppr)
CalculateIA = function(GL){
  name = deparse(substitute(GL))
  print(paste0("Beginning analysis for ", name))
  IA = win.ia(GL, window = 500, threads=0)
  write.csv(IA, paste0("Pc136_583k.Ia_500win.",name,".csv"))
  png(paste0("Pc136_583k.Ia_500win.",name, ".png"), width = 1200, height = 750)
  hist(IA, breaks = "fd")
  dev.off()
}

# Basic Stats (dartR/hierfstat)
Basic_Stats = function(GL){
  Stats = utils.basic.stats(GL)
  name = deparse(substitute(GL))
  write.csv(Stats$Ho, paste0("Pc_583k.Ho.",name,".csv"))
  write.csv(Stats$Hs,  paste0("Pc_583k.Hs.",name,".csv"))
  write.csv(Stats$perloc,  paste0("Pc_583k.perloc.",name,".csv"))
  write.csv(Stats$Fis,  paste0("Pc_583k.Fis.",name,".csv"))
  write.csv(Stats$overall,  paste0("Pc_583k.overall.",name,".csv"))
}

# Hedrick's and Nei's Gst (mmod)
PW_Gst_CC = function(Gind){
  name = deparse(substitute(Gind))
  HenGst = pairwise_Gst_Hedrick(Gind) 
  HenGST_DF<-as.data.frame(as.matrix(HenGst))
  write.csv(HenGST_DF,file = paste0("PW_HenGst.",name,".csv"))
  NieGst = pairwise_Gst_Nei(Gind)
  NeiGST_DF<-as.data.frame(as.matrix(NieGst))
  write.csv(NeiGST_DF,file = paste0("PW_NeiGst.",name,".csv"))
}

#-------------------------------------------------------------------------------
#               Calculating the Index of Association (IA)
#                                Table 2
#-------------------------------------------------------------------------------
# Uncorrected 
CalculateIA(vcf_GL_OG)
CalculateIA(CAnorth)
CalculateIA(MX)
CalculateIA(CAsouth)
CalculateIA(Asia)
# Clone-Corrected
CalculateIA(vcf_GL_OG_SC2_CC)
CalculateIA(CAnorthMLG0.01)
CalculateIA(MX_MLG0.01)
CalculateIA(CAsouthMLG0.01)
CalculateIA(AsiaMLG0.01)

#-------------------------------------------------------------------------------
#                         Generating Basic Stats 
#                                Table 2
#-------------------------------------------------------------------------------
# Uncorrected 
Basic_Stats(vcf_GL_OG)
Basic_Stats(CAnorth)
Basic_Stats(MX)
Basic_Stats(CAsouth)
Basic_Stats(Asia)
# Clone-Corrected
Basic_Stats(vcf_GL_OG_SC2_CC)
Basic_Stats(CAnorthMLG0.01)
Basic_Stats(MX_MLG0.01)
Basic_Stats(CAsouthMLG0.01)
Basic_Stats(AsiaMLG0.01)

#-------------------------------------------------------------------------------
#                        Calculating Pairwise Gst 
#                                Table 3
#-------------------------------------------------------------------------------
# Uncorrected 
PW_Gst_CC(vcf_GL_OG)
# Clone-Corrected
PW_Gst_CC(vcf_GL_OG_SC2_CC)

#-------------------------------------------------------------------------------
#                         GGtree & FastStructure
#.                              Figure 1 
#-------------------------------------------------------------------------------
# Raw UPGMA tree with 1000 BS and adding the FastStructure LogPrior data (k=4)

# Trees
Tree = read.beast("Pc136_583K_UPGMA_Tree_1000BSR_FigtreeOut.nex") 

FastStructurek4.v2 = FastStructurek4 %>% transform(ID = as.character(ID),
                                   Cluster1=
                                     as.numeric(as.character(Cluster1)),
                                   Cluster2 =
                                     as.numeric(as.character(Cluster2)), 
                                   Cluster3 =
                                     as.numeric(as.character(Cluster3)), 
                                   Cluster4 =
                                     as.numeric(as.character(Cluster4)))

FastStructurek4.v4 <- melt(FastStructurek4.v2)
# Setting colors 
ClusterCols=c("Cluster1" = "#00A0B0", "Cluster2" = "#6A4A3C",
              "Cluster3" = "#CC333F", "Cluster4" = "#EB6841")

# Plotting Tree
p = ggtree(Tree, right = T) %<+% Annot +
  geom_tiplab(aes(label = name), linetype=3, 
              linesize=.1, size=1.5, hjust = -0.2)+ 
  geom_point(aes(color = location_Simple), size=1.5) +
  geom_point2(aes(subset=!isTip,
                  fill=cut(as.numeric(Bootstrap), c(0, 75, 90, 100))), shape=21, size=1)+ 
  scale_fill_manual(values=c("black", "grey", "white"), guide='legend', 
                    name='Bootstrap Percentage (BP)', 
                    breaks=c('(90,100]', '(75,90]', '(0,75]'), 
                    labels=expression(BP>=90, 75 <= BP * " < 90", BP < 75))+ 
  scale_color_manual(values=LocationCols, name="Location")+ 
  geom_treescale(width = 0.01, fontsize = 3) + 
  ggplot2::xlim(0, NA) + 
  theme(legend.position = c(0.2, 0.2)) + 
  coord_cartesian(clip="off") 
p=p + new_scale("fill") 

# Fast Structure  
p3 = facet_plot(p,
                panel = "Fast Structure (K=4)",
                data = FastStructurek4.v4,
                geom = geom_barh,
                mapping = aes(x = value, fill = variable),
                stat = "identity") + scale_fill_manual(values=ClusterCols, 
                                                       name="Genetic Cluster")

p3 = p3 + xlim_expand(c(0.15, NA), 'Tree')

p3 = p3 + theme_classic() + theme(legend.position = c(0.2, 0.3), 
                                  legend.box="vertical", 
                                  legend.margin=margin(),
                                  legend.key.size = unit(4, 'mm'),
                                  legend.text = element_text(size = 8),
                                  legend.title =element_text(size = 10)) 
p3

# Reducing the size of the fastStructure panel
gt = ggplot_gtable(ggplot_build(p3))
gt$layout$l[grep('Fast Structure (K=4)', gt$layout$name)] 
gt$widths[7] = 0.5*gt$widths[7] 
grid.draw(gt) 

#ggsave(p3,filename = "Pc583k_UPGMA1kBS_FastStructureK4.v2.pdf", 
       #width = 270, height = 360, units = "mm", limitsize = FALSE)

#-------------------------------------------------------------------------------
#                                GGMap
#.                              Figure 1 
#-------------------------------------------------------------------------------
# Google API key
ggmap::register_google(key = "acquire your own ;)")
apikey= "acquire your own ;)"

# Importing JSON with google map styles (applied to all regions) 
style_list <- fromJSON("annotations/MapStlye.json", simplifyVector = FALSE)
style_list <- as.list(style_list)
Style_String = ggmap_style_string(style_list)

#-------------
# Southern CA 
#-------------
# Importing annotations for pie charts 
CA_Annot = read.csv('annotations/GGMap_Pc_Cluster_Loc_Data_CA.csv', header = T)
CA_Annot = as.data.frame(CA_Annot)
CA_Annot2 = CA_Annot %>% transform(latitude =
                                     as.numeric(as.character(latitude)),
                                   longitude =
                                     as.numeric(as.character(longitude)),
                                   Cluster1=
                                     as.numeric(as.character(Cluster1)),
                                   Cluster2 =
                                     as.numeric(as.character(Cluster2)), 
                                   Cluster3 =
                                     as.numeric(as.character(Cluster3)), 
                                   Cluster4 =
                                     as.numeric(as.character(Cluster4)), 
                                   region = 
                                     as.numeric(as.character(region)))

# Define the bounding box coordinates
bbox_bound_CA <- c(left = -121.4, bottom = 32.4754703838965,
                   right = -115, top = 34.90665090116335)

# Calculate the center point of the bounding box
center_lon_CA <- mean(bbox_bound_CA[c("left", "right")])
center_lat_CA <- mean(bbox_bound_CA[c("bottom", "top")])

# Specify the zoom level
zoom <- 8

# Get the map using the bounding box coordinates and zoom level
CA3 <- get_googlemap(center = c(lon = center_lon_CA, lat = center_lat_CA), 
                     zoom = zoom, style = Style_String, bbox = bbox_bound_CA)

CA3M = ggmap(CA3)

P7 = CA3M + geom_scatterpie(aes(x=longitude, y=latitude, group=region), 
                            data=CA_Annot2, pie_scale = 0.4,
                            cols=colnames(CA_Annot2)[5:8], 
                            legend_name = "Genetic_Cluster") +
  scale_fill_manual(values=ClusterCols) +
  theme(legend.position = 'none')

P7
#ggsave(P7,filename = "CA.v5.BW.pdf", width = 12.59, height = 7.06, 
       #units = "in", limitsize = FALSE)

#-------------
# Michoacan MX
#-------------
# Importing annotations for pie charts 
MX_Annot = read.csv('annotations/GGMap_Pc_Cluster_Loc_Data_MX.csv', header = T)
MX_Annot = as.data.frame(MX_Annot)
MX_Annot2 = MX_Annot %>% transform(latitude =
                                     as.numeric(as.character(latitude)),
                                   longitude =
                                     as.numeric(as.character(longitude)),
                                   Cluster1=
                                     as.numeric(as.character(Cluster1)),
                                   Cluster2 =
                                     as.numeric(as.character(Cluster2)), 
                                   Cluster3 =
                                     as.numeric(as.character(Cluster3)), 
                                   Cluster4 =
                                     as.numeric(as.character(Cluster4)), 
                                   region = 
                                     as.numeric(as.character(region)))

# Define the bounding box coordinates
bbox_bound_MX <- c(left = -104, bottom = 18,
                   right = -99.8, top = 20.5)

# Calculate the center point of the bounding box
center_lon_MX <- mean(bbox_bound_MX[c("left", "right")])
center_lat_MX <- mean(bbox_bound_MX[c("bottom", "top")])

# Specify the zoom level
zoom <- 8

# Get the map using the bounding box coordinates and zoom level
MX3 <- get_googlemap(center = c(lon = center_lon_MX, lat = center_lat_MX), 
                     zoom = zoom, style = Style_String, bbox = bbox_bound_MX)

MX3M = ggmap(MX3)

P8 = MX3M + geom_scatterpie(aes(x=longitude, y=latitude, group=region), 
                            data=MX_Annot2, pie_scale = 0.5,
                            cols=colnames(MX_Annot2)[5:8], 
                            legend_name = "Genetic_Cluster") +
  scale_fill_manual(values=ClusterCols) +
  theme(legend.position = 'none')

P8
#ggsave(P8,filename = "MX.v5.BW.pdf", width = 12.59, height = 7.06, 
       #units = "in", limitsize = FALSE)

#-------------
# Asia/Oceania
#-------------
# Importing annotations for pie charts 
ASIA_Annot = read.csv('annotations/GGMap_Pc_Cluster_Loc_Data_Asia.csv', header = T)
ASIA_Annot = as.data.frame(ASIA_Annot)
ASIA_Annot2 = ASIA_Annot %>% transform(latitude =
                                         as.numeric(as.character(latitude)),
                                       longitude =
                                         as.numeric(as.character(longitude)),
                                       Cluster1=
                                         as.numeric(as.character(Cluster1)),
                                       Cluster2 =
                                         as.numeric(as.character(Cluster2)), 
                                       Cluster3 =
                                         as.numeric(as.character(Cluster3)), 
                                       Cluster4 =
                                         as.numeric(as.character(Cluster4)), 
                                       region = 
                                         as.numeric(as.character(region)))

# Define the bounding box coordinates
bbox_bound_AS <- c(left = 98, bottom = -10,
                   right = 150, top = 30)

# Calculate the center point of the bounding box
center_lon_AS <- mean(bbox_bound_AS[c("left", "right")])
center_lat_AS <- mean(bbox_bound_AS[c("bottom", "top")])

# Specify the zoom level
zoom <- 4

# Get the map using the bounding box coordinates and zoom level
AS3 <- get_googlemap(center = c(lon = center_lon_AS, lat = center_lat_AS), 
                     zoom = zoom, style = Style_String, bbox = bbox_bound_AS)

AS3M = ggmap(AS3)

P9 = AS3M + geom_scatterpie(aes(x=longitude, y=latitude, group=region), 
                            data=ASIA_Annot2, pie_scale = 0.4,
                            cols=colnames(ASIA_Annot2)[5:8], 
                            legend_name = "Genetic_Cluster") +
  scale_fill_manual(values=ClusterCols) +
  theme(legend.position = 'none')

P9
#ggsave(P9,filename = "Asia.v5.BW.pdf", width = 12.59, height = 7.06, 
       #units = "in", limitsize = FALSE)




