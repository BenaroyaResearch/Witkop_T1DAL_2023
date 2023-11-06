#### P362_1 T1DAL Correlation with Pathology Indicators ####

# This script assess correlation of the percent of cells in each cluster with T1D pathology indicators

#### Load Libraries ####

library(tidyverse)
library(monocle3)
library(Matrix)
library(ggpubr)
library(gridExtra)
library(ComplexHeatmap)
library(apird)

setwd("/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup")
set.seed(42) # same seed as QC script 

#### LOAD ALL SAVED DATA AND SET PATHS ####

# Set Paths 
baseDir <- "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup"
plotDir <- file.path(baseDir,'FIGURES')
resultDir <- file.path(baseDir,'SAVED_DATA')
annotationDir <- file.path(baseDir, "RAW_DATA")

# Load data with cluster 9 removed 
load( file = file.path(resultDir, "cds_no_MAIT_no_9.Rdata"))

# load anno data
P362_libs <- apird::getProjectLibs("P362-1")
P362_anno <- apird::getAnno(P362_libs) 

# Set colors
R_NR_colors <- c("#46c19a","#6d80d8")

### Load QR and c-peptide data from Matt D 
# Notes from Matt D on this file: This file has the C-peptide slopes and QR for all subjects, including the 12 that you listed.  
# I checked that the identifiers match, and all 12 subjects are present in the data I have.  For QR, I 
# included the expected and observed 2-hour AUC values (in log (mean nmol/L + 1)), as well as the QR itself. 
# The attached plot shows those values for the placebo- and active-treated groups - QR is just the vertical distance 
# from the 1:1 line.  One thing you’re probably already aware of, subject 10241 had an extremely rapid loss of C-peptide.
# So much so that it throws off the slope calculation when using the random-effects model.  
# So I have a variation of the model that excludes that subject.  The file includes those two models,
# as well as slopes from a model where each subject gets their own independent slope.

T1DAL_QR_cpeptide <- read.csv(file.path(annotationDir,"T1DAL_cpeptide_rates_QR_for_EWitkop.2023-09-29.csv"))
colnames(T1DAL_QR_cpeptide)[2] <-"Donor.ID" # change colnames to match the UMAP metadata
T1DAL_QR_cpeptide$Donor.ID <- as.character(T1DAL_QR_cpeptide$Donor.ID)

#### Calculate percent of cells in each cluster and join with metadata ####

# Calculate the percent of exhausted cells per cluster per donor
cds_no_MAIT_no_9_percent <- as.data.frame(colData(cds_no_MAIT_no_9)) %>% 
  # get total per donor
  group_by(Donor.ID) %>%
  mutate(total_donor = n()) %>%
  # get total per cluster per donor
  group_by(Cluster.Name, Donor.ID) %>% 
  mutate(total_per_cluster = n(),
         percent_per_cluster = total_per_cluster/total_donor*100) %>%
  distinct(Donor.ID, Response, Cluster.Name, total_donor, total_per_cluster, percent_per_cluster )

# Join with QR and cpeptide data - join by Donor ID
cds_no_MAIT_no_9_percent_path <- left_join(cds_no_MAIT_no_9_percent, T1DAL_QR_cpeptide)

#### ASSESS CORRELATION BETWEEN EXHAUSTED CLUSTERS AND PATHOLOGY INDICATORS ####

column_list <- c("cpep_slope_random", "cpep_slope_random_ex_outlier", "cpep_slope_fixed","auc2hr_expected_QR",
                 "auc2hr_observed_QR","QR")

# write function to plot all correlation plots
plot_cor <- function(x){
  variable <- rlang::sym(x)
  plot <-  cds_no_MAIT_no_9_percent_path %>%
    filter(Cluster.Name %in% c(5,6,7,8)) %>%
    ggplot(aes(y = !!variable, x = percent_per_cluster, color = Response))+
    geom_point() + theme_minimal() + 
    geom_smooth(method='lm', formula= y~x) + 
    stat_cor( label.y.npc = 0.15, label.x.npc  = 0.5, show.legend = FALSE, size = 6) +
    scale_color_manual(values= R_NR_colors) +
    labs(x= "% Cells per Tex Cluster",y=str_to_title(str_replace_all(x, "_"," ")),
         title = paste(str_to_title(str_replace_all(x, "_"," ")), "vs % Tex")) +
    theme(text = element_text(size = 16))
  
}

cpep_slope_random <- plot_cor("cpep_slope_random")
cpep_slope_random_ex_outlier <- plot_cor("cpep_slope_random_ex_outlier")
cpep_slope_fixed <- plot_cor("cpep_slope_fixed")
auc2hr_expected_QR <- plot_cor("auc2hr_expected_QR")
auc2hr_observed_QR <- plot_cor("auc2hr_observed_QR")
QR <- plot_cor("QR")

# export plots
all_Tex_cor_plot <- ggarrange(cpep_slope_random, cpep_slope_random_ex_outlier, cpep_slope_fixed, auc2hr_expected_QR, 
                              auc2hr_observed_QR, QR, nrow = 2, ncol = 3)

ggsave(all_Tex_cor_plot, file = file.path(plotDir, "all_Tex_cor_plot.pdf"), height = 12, width =20)

#### FIGURE 6A, B: PLOT CORRELATION WITH C PEPTIDE FIXED AS BAR PLOT PER CLUSTER ####

# Replot just c-peptide slppe fixed for each Tex cluster separately
cpeptide_fixed_facet <-  cds_no_MAIT_no_9_percent_path %>%
  filter(Cluster.Name %in% c(5,6,7,8)) %>%
  ggplot(aes(y = cpep_slope_fixed, x = percent_per_cluster))+
  geom_point() + theme_minimal() + 
  #scale_color_manual(values= R_NR_colors) +
  geom_smooth(method='lm', formula= y~x) + 
  stat_cor( label.y.npc = 0.15, label.x.npc  = 0.1, show.legend = FALSE, size = 4) +
  labs(x= "% Cells per Tex Cluster",y= "C-peptide Slope") +
  theme(text = element_text(size = 12)) + facet_grid(.~Cluster.Name)

## FIGURE 6A:
ggsave(cpeptide_fixed_facet, file = file.path(plotDir,"cpeptide_fixed_facet.pdf"), width = 8, height = 5)

## Plot top and bottom quartile of C-peptide fixed slope 

cpeptide_fixed_quartile <- cds_no_MAIT_no_9_percent_path %>%
  group_by(Cluster.Name) %>%
  mutate(cpeptide_slope_fixed_quartile = dplyr::ntile(cpep_slope_fixed, 4))

cpeptide_fixed_quartile$cpeptide_slope_fixed_quartile <- factor(cpeptide_fixed_quartile$cpeptide_slope_fixed_quartile,
                                                                levels = c("1","2","3","4"),
                                                                labels = c("Bottom 25%","Bottom higher","Top Lower","Top 25%"))

# filter for top and bottom quartile and plot
cpeptide_fixed_quartile_Tex_plot <- cpeptide_fixed_quartile %>% 
  filter(cpeptide_slope_fixed_quartile %in% c("Bottom 25%","Top 25%")) %>%
  filter(Cluster.Name %in% c(5,6,7,8)) %>%
  distinct(Donor.ID, Cluster.Name, percent_per_cluster, cpeptide_slope_fixed_quartile, Response) %>%
  # plot
  ggplot( aes(x = cpeptide_slope_fixed_quartile, y = percent_per_cluster, color=Response)) +
  geom_boxplot(outlier.fill = NA) + geom_point() +facet_grid(.~Cluster.Name)+
  scale_color_manual(values= R_NR_colors) +
  theme_minimal() + labs(y = "% Cells per Tex Cluster", x = "C-peptide Slope Quartile") +
  stat_compare_means(show.legend = FALSE) +
  theme(panel.border = element_rect(colour = "black",fill = NA),text = element_text(size = 12))

## FIGURE 6B:
ggsave(cpeptide_fixed_quartile_Tex_plot, file = file.path(plotDir,"cpeptide_fixed_quartile_Tex_plot.pdf"), width = 8, height = 5)
