# HEADER -----------------------------------------------------------------------
#
# TITLE:   6-PCA_BY_HMO_CLUSTERS.R
#
# PURPOSE: do a PCA analysis of miRNAs by HMO cluster
#
# DATE:    March 14, 2022
#
# SET UP -----------------------------------------------------------------------

#clear the workspace
rm(list = ls())

#turn off scientific notation
options(scipen = 100)

#load libraries
pacman::p_load(knitr, tidyverse, magrittr, lme4, lmerTest, GGally, corrplot, 
               Hmisc, kableExtra, dplyr, plyr, janitor, lubridate, survminer, 
               ggplot2, here, readr, tableone, officer, flextable,finalfit,
               purrr, stringr, lme4, corrplot, pscl, stargazer, MASS, lmerTest,
               readxl, factoextra, ggbiplot) 

# REad in meta data
meta <- read.csv("input/meta_clean_EY.csv")
meta$dyad_id <- as.character(meta$dyad_id)

# ------------------------------------------------------------------------------
# Do miRNA PCA: ####
miRNA_cpm <- read.csv("input/miRNA_counts_EY.csv") # reading in for colnames
# Read in PCA data from script 2
miRNA_cpm_noQC <- read_csv("input/miRNA_cpm_noQC_EY.csv")

# calculate principle components
# miRNA.pca <- prcomp(miRNA_cpm_noQC, center = T, scale = T) # non numeric error
miRNA.pca <- prcomp(miRNA_cpm_noQC[,2:210], center = T, scale = T)

# check whether we can use eigenvalue > 1 as a metric
I((summary(miRNA.pca)$sdev)^2) #this might not be NB. Not great. Too many
# eigenvalue = standard deviation^2. Too many to use those > 1.
plot(miRNA.pca)
#save the PCA results in a data frame
PCA <- data.frame(miRNA.pca$x)

meta_miRNA <- left_join(meta, PCA, by = "dyad_id")
miRNA_cpm$dyad_id <- paste0("X", miRNA_cpm$dyad_id)

#meta_miRNA <- left_join(meta_miRNA, miRNA_cpm[211:217], by = "dyad_id")
#editting the line above:
miRNA_cpm_cols <- dplyr::select(miRNA_cpm, c("dyad_id","LibSize", 
                                             "LibSizeNormalized", "file_name","LL", 
                                             "Address","Barcode",
                                             "prop_rRNA","Vol_Supernatant","Date_Evs"))
meta_miRNA <- left_join(meta_miRNA, miRNA_cpm_cols, by = "dyad_id")

colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC1"] <- "PC1"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC2"] <- "PC2"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC3"] <- "PC3"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC4"] <- "PC4"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC5"] <- "PC5"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC6"] <- "PC6"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC7"] <- "PC7"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC8"] <- "PC8"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC9"] <- "PC9"

temp <- data.frame(t(summary(miRNA.pca)$importance))
temp$PC <- row.names(temp)
temp <- temp[1:10,]
temp$PC <- factor(temp$PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                      "PC7", "PC8", "PC9", "PC10"))

plot_PC <- ggplot(data = temp, aes(x = PC, y = Proportion.of.Variance)) + 
  geom_bar(stat = "identity") +
  ylab("Proportion of Variance") +
  xlab("")
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
  theme(legend.position = "none")

plot_PC

# ------------------------------------------------------------------------------
# Cluster HMOs: ####

# Read in HMO clustering data from script 3
meta_hmo_trim_ID <- read_csv("input/meta_hmo_trim_EY.csv")
# Only Secretos 
hmo_sec_trim_ID <- read_csv("input/hmo_sec_trim_EY.csv")

#make numeric
meta_hmo_trim <- sapply(meta_hmo_trim_ID, as.numeric)
hmo_sec_trim <- sapply(hmo_sec_trim_ID, as.numeric)

# scale 
meta_hmo_trim_scaled <- data.frame(scale(meta_hmo_trim))
hmo_sec_trim_scaled <- data.frame(scale(hmo_sec_trim))

# now calculate the distance matrix
dist_mat <- dist(meta_hmo_trim_scaled, method = 'euclidean')
dist_mat_sec <- dist(hmo_sec_trim_scaled, method = "euclidean")
#Why Euclidian?

# cluster
hclust_avg <- hclust(dist_mat, method = "complete")
hclust_sec <- hclust(dist_mat_sec, method = "complete")

hclust <- cutree(hclust_avg, 6)
names(hclust) <- miRNA_cpm$dyad_id
#hclust <- as.data.frame(hclust)

secretors <- meta$dyad_id[which(meta$Secretor == "Yes")]
miRNA_sec <- miRNA_cpm[which(miRNA_cpm$X %in% secretors),]
hclust_sec2 <- cutree(hclust_sec, 6)
names(hclust_sec2) <- miRNA_sec$dyad_id
#hclust_sec2 <- as.data.frame(hclust_sec2)

# get some info about the clustering
table(hclust)
#1   2   3   4   5   6 
#81 101  15   2   1   7 
table(hclust_sec2)
#1  2  3  4 
#99 76  1  7

#make a list of individuals in cluster 1 - cluster 3 for all: 
cluster_1 <- hclust[hclust == 1]
cluster_2 <- hclust[hclust == 2]
cluster_3 <- hclust[hclust == 3]
cluster_4 <- hclust[hclust == 4]
cluster_5 <- hclust[hclust == 5]
cluster_6 <- hclust[hclust == 6]

# add cluster as a variable in meta
meta$cluster <- NA
meta$cluster[row.names(meta) %in% names(cluster_1)] <- 1
meta$cluster[row.names(meta) %in% names(cluster_2)] <- 2
meta$cluster[row.names(meta) %in% names(cluster_3)] <- 3
meta$cluster[row.names(meta) %in% names(cluster_4)] <- 4
meta$cluster[row.names(meta) %in% names(cluster_5)] <- 5
meta$cluster[row.names(meta) %in% names(cluster_6)] <- 6

summary(meta$cluster)

hc <- meta[which(!is.na(meta$cluster)),]
hc$cluster <- as.character(hc$cluster)

#only keep individuals in meta who have miRNA data
PCA$dyad_id <- row.names(PCA)
hc_pca <- hc[hc$dyad_id %in% row.names(PCA),] #51 in cluster 1-3 with miRNA

# ------------------------------------------------------------------------------
# visualizations for HMO cluster
ggbiplot(miRNA.pca, group = hc$cluster, ellipse = T, 
         var.axes = F) +
  theme_minimal() +
  labs(color = "hmo_cluster") +
  ylim(-4,2.5) + xlim(-3,3) +
  ylab("PC2 (explained miRNA var.)") + xlab("PC1 (explained miRNA var.)")

#
