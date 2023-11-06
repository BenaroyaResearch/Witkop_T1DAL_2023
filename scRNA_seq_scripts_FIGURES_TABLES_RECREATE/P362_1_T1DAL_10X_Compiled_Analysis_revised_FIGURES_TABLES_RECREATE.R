### This script produced scRNA-seq UMAPS and gene list heatmaps

#### LOAD LIBRARIES ####

library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(scater)
library(stringr)
library(ggpubr)
library(gridExtra)
library(car)
library(lme4)
library(garnett)
library(ggrepel)
library(ComplexHeatmap)
library(viridis)
library(tidyverse)
library(apird)
library(org.Hs.eg.db)
library(limma)
library(egg)
library(cowplot)


setwd("/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup")
set.seed(42) # same seed as QC script 

#### LOAD ALL SAVED DATA AND SET PATHS ####

# Set Paths 
baseDir <- "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup"
plotDir <- file.path(baseDir,'FIGURES')
resultDir <- file.path(baseDir,'SAVED_DATA')
annotationDir <- file.path(baseDir, "RAW_DATA")

# Load saved annotations 
load(file= file.path(resultDir, "P362-1_annotation.Rdata"))

# Load post QC no MAIT cell data
load(file= file.path(resultDir,"P362-1 T1DAL cds object - postQC no MAIT cells.Rdata"))

P362_anno_IDs <- P362_anno %>% dplyr::rename(Lib.ID = "libid", Donor.ID = "participantID") %>% dplyr::select(Lib.ID, Donor.ID)

# rename for downstream analysis to preserve the original
cds_data <- cds_no_MAIT
rm(cds_no_MAIT)

# set color pallette
pal = toupper(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffed6f','#b15928', "gray","black","blue","red"))    

## load ITN dictionary matching operational donor IDs to the correct public IDs
T1DAL_ITN_ID_dictionary <- readxl::read_xlsx(file.path(annotationDir, "t1dal_mask_id_EW_formatted.xlsx"))
colnames(T1DAL_ITN_ID_dictionary) <- c("operational_PID", "masked_public_PID","Donor.ID")
# remove T1DAL from public_PID column
T1DAL_ITN_ID_dictionary$masked_public_PID <- str_remove(T1DAL_ITN_ID_dictionary$masked_public_PID, "T1DAL_")

#### FIGURE 3B PLOT UMAP WITH PSEUDOTIME ####
  
# remove cluster 9 from plot
cds_data_cluster_no_MAIT_no_9 <- plot_cells(cds_data[,colnames(cds_data) %in%row.names(colData(cds_data) %>% as.data.frame() %>%  filter(Cluster.Name !="9")) ],
                                            color_cells_by="Cluster.Name", 
                                            show_trajectory_graph = T, cell_size = 2, alpha=0.5,label_cell_groups = F,
                                            label_branch_points = F, label_roots  = F, label_leaves = F) + 
  scale_color_manual(values=pal) +
  theme(text = element_text(size = 20))
# PLOT FIGURE 3B - note annotations were added afterward in Inkscape
ggsave(cds_data_cluster_no_MAIT_no_9, file = file.path(plotDir,"cds_data_cluster_no_MAIT_no_9.pdf"), device = "pdf", width= 8, height = 6)


#### FIGURE 3C: Z-score scaled heatmap of key marker genes ####

# final addition version
figure_ex_gene_list <- c("IL7R","HAVCR2","CTLA4","KLRG1","FCGR3A","TIGIT","LAG3","PDCD1","EOMES","TOX","KIR3DL1","KIR2DS4",
                         "TBX21","KLRD1","LILRB1","KIR2DL3","KIR2DL1","KIR3DL2", "CD160", "NCR1","NKG7","GZMB","GZMA","GZMH", "CD244", "IKZF2",
                         "CX3CR1","SLAMF6", "GZMK")


data_for_figure1 <- as.data.frame(t(exprs(cds_data[rowData(cds_data)$gene_short_name %in% figure_ex_gene_list,])) )
length(figure_ex_gene_list ) == ncol(data_for_figure1) # TRUE
data_for_figure1$Cell.ID <- row.names(data_for_figure1)
colData(cds_data)$Cluster <- colData(cds_data)$Cluster.Name

anno_data <- as.data.frame(colData(cds_data)) #%>% select(-CCR7, -CD28, -IL7R, -IFNG, -MKI67) # this was originally commented out by Kirsten
anno_data$Cell.ID <- row.names(colData(cds_data))

# Calculate the mean of genes of interest, removing cluster 9
data_for_figure1_anno <-  data_for_figure1 %>%
  merge(anno_data, by = "Cell.ID") %>%
  filter_at(vars(figure_ex_gene_list ), any_vars(. > 0)) %>%
  group_by(Cluster.Name, Timepoint) %>%
  summarise_at(vars(figure_ex_gene_list), ~mean(. > 0)) %>%
  # remove cluster 9 for similarity to Kirsten's plot
  filter(Cluster.Name != 9)

cluster_medians_figure1 <- data_for_figure1_anno

## Use scale to subtract the mean values and get zscore
# ?scale(): If center is a numeric-alike vector with length equal to the number of columns of x,
#then each column of x has the corresponding value from center subtracted from it.
col_scaled_median_figure1 <- scale(cluster_medians_figure1[,figure_ex_gene_list])
row.names(col_scaled_median_figure1) <- cluster_medians_figure1$Cluster.Name


IR_gene_list_zscore_heatmap <- Heatmap(t(col_scaled_median_figure1),row_names_side="right",name="Z-score\nExpression",
                                       show_parent_dend_line=F, cluster_rows = T, cluster_columns = T,rect_gp = gpar(col = "white", lwd = 1),
                                       column_names_gp = grid::gpar(fontsize = 16), row_names_gp = grid::gpar(fontsize = 16), 
                                       heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12), title_gp = grid::gpar(fontsize = 16)))

IR_gene_list_zscore_heatmap_grob <- grid.grabExpr(draw(IR_gene_list_zscore_heatmap))

# FIGURE 3C:
pdf(file = file.path(plotDir,"IR_gene_list_zscore_heatmap.pdf"), width = 5, height = 6)
IR_gene_list_zscore_heatmap
dev.off()

#### FIGURE S15A, FIGURE S15B: Plot heatmaps with published gene lists ####

# load gene lists of interest
Daniel_Wherry_2022_Figure2b_Texterm_TexKLR <- readxl::read_xlsx(file.path(annotationDir,"Daniel_Wherry_2022_Figure2b_Texterm_TexKLR.xlsx"))
Giles_Wherry_Figure_3b_gene_list <- readxl::read_xlsx(file.path(annotationDir,"Giles_Wherry_Figure_3b_gene_list.xlsx"))
Figure_S3A_heatmap_gene_set <- readxl::read_xlsx(file.path(annotationDir,"Figure_S3A_heatmap_gene_set.xlsx"))

## C: Daniel_Wherry_2022_Figure2b_Texterm_TexKLR

data_for_figure1 <- as.data.frame(t(exprs(cds_data[rowData(cds_data)$gene_short_name %in% Daniel_Wherry_2022_Figure2b_Texterm_TexKLR$gene_name,])) )
length(Daniel_Wherry_2022_Figure2b_Texterm_TexKLR$gene_name ) == ncol(data_for_figure1) # TRUE
data_for_figure1$Cell.ID <- row.names(data_for_figure1)
colData(cds_data)$Cluster <- colData(cds_data)$Cluster.Name

anno_data <- as.data.frame(colData(cds_data)) #%>% select(-CCR7, -CD28, -IL7R, -IFNG, -MKI67) # this was originally commented out by Kirsten
anno_data$Cell.ID <- row.names(colData(cds_data))

# Calculate the mean of genes of interest, removing cluster 9
data_for_figure1_anno <-  data_for_figure1 %>%
  merge(anno_data, by = "Cell.ID") %>%
  filter_at(vars(Daniel_Wherry_2022_Figure2b_Texterm_TexKLR$gene_name ), any_vars(. > 0)) %>%
  group_by(Cluster.Name, Timepoint) %>%
  summarise_at(vars(Daniel_Wherry_2022_Figure2b_Texterm_TexKLR$gene_name), ~mean(. > 0)) %>%
  # remove cluster 9 for similarity to Kirsten's plot
  filter(Cluster.Name != 9)

cluster_medians_figure1 <- data_for_figure1_anno

## Use scale to subtract the mean values and get zscore
# ?scale(): If center is a numeric-alike vector with length equal to the number of columns of x,
#then each column of x has the corresponding value from center subtracted from it.
col_scaled_median_figure1 <- scale(cluster_medians_figure1[,Daniel_Wherry_2022_Figure2b_Texterm_TexKLR$gene_name])
row.names(col_scaled_median_figure1) <- cluster_medians_figure1$Cluster.Name

# cols = c(brewer.pal(n=12,name="Paired"),brewer.pal(n=3,name="RdBu")[1],brewer.pal(n=3,name="PiYG")[1],"black","gray")
# cols = cols[sample(1:length(cols))]
cols = c("#FFFF99","#1F78B4","#FB9A99","#FDBF6F","#B2DF8A","#FF7F00","gray","#EF8A62","#CAB2D6","#E31A1C","#E9A3C9","black","#6A3D9A","#A6CEE3","#B15928", "#33A02C")
cluster_cols <- cols

Daniel_Wherry_2022_Figure2b_Texterm_TexKLR_gene_list_zscore_heatmap <- 
  Heatmap(t(col_scaled_median_figure1),row_names_side="right",name="Z-score\nExpression",
          show_parent_dend_line=F, cluster_rows = T, cluster_columns = T,rect_gp = gpar(col = "white", lwd = 1),
          column_names_gp = grid::gpar(fontsize = 16), row_names_gp = grid::gpar(fontsize = 16), 
          heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12), title_gp = grid::gpar(fontsize = 16)))


## PLOT FIGURE S15B
pdf(file = file.path(plotDir, "Daniel_Wherry_2022_Figure2b_Texterm_TexKLR_gene_list_zscore_heatmap.pdf"), width = 5, height = 9)
Daniel_Wherry_2022_Figure2b_Texterm_TexKLR_gene_list_zscore_heatmap
dev.off()

### PLOT Li et al gene signature
## Figure_S3A_heatmap_gene_set

data_for_figure1 <- as.data.frame(t(exprs(cds_data[rowData(cds_data)$gene_short_name %in% Figure_S3A_heatmap_gene_set$gene_name,])) )
length(Giles_Wherry_Figure_3b_gene_list$gene_name ) == ncol(data_for_figure1) # TRUE
data_for_figure1$Cell.ID <- row.names(data_for_figure1)
colData(cds_data)$Cluster <- colData(cds_data)$Cluster.Name

anno_data <- as.data.frame(colData(cds_data)) #%>% select(-CCR7, -CD28, -IL7R, -IFNG, -MKI67) # this was originally commented out by Kirsten
anno_data$Cell.ID <- row.names(colData(cds_data))

# Calculate the mean of genes of interest, removing cluster 9
data_for_figure1_anno <-  data_for_figure1 %>%
  merge(anno_data, by = "Cell.ID") %>%
  filter_at(vars(Figure_S3A_heatmap_gene_set$gene_name ), any_vars(. > 0)) %>%
  group_by(Cluster.Name, Timepoint) %>%
  summarise_at(vars(Figure_S3A_heatmap_gene_set$gene_name), ~mean(. > 0)) %>%
  # remove cluster 9 for similarity to Kirsten's plot
  filter(Cluster.Name != 9)

cluster_medians_figure1 <- data_for_figure1_anno

## Use scale to subtract the mean values and get zscore
# ?scale(): If center is a numeric-alike vector with length equal to the number of columns of x,
#then each column of x has the corresponding value from center subtracted from it.
col_scaled_median_figure1 <- scale(cluster_medians_figure1[,Figure_S3A_heatmap_gene_set$gene_name])
row.names(col_scaled_median_figure1) <- cluster_medians_figure1$Cluster.Name

# cols = c(brewer.pal(n=12,name="Paired"),brewer.pal(n=3,name="RdBu")[1],brewer.pal(n=3,name="PiYG")[1],"black","gray")
# cols = cols[sample(1:length(cols))]
cols = c("#FFFF99","#1F78B4","#FB9A99","#FDBF6F","#B2DF8A","#FF7F00","gray","#EF8A62","#CAB2D6","#E31A1C","#E9A3C9","black","#6A3D9A","#A6CEE3","#B15928", "#33A02C")
cluster_cols <- cols

Figure_S3A_heatmap_gene_set_zscore_heatmap <- 
  Heatmap(t(col_scaled_median_figure1),row_names_side="right",name="Z-score\nExpression",
          show_parent_dend_line=F, cluster_rows = T, cluster_columns = T,rect_gp = gpar(col = "white", lwd = 1),
          column_names_gp = grid::gpar(fontsize = 16), row_names_gp = grid::gpar(fontsize = 16), 
          heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12), title_gp = grid::gpar(fontsize = 16)))

### Plot FIGURE S15A
pdf(file = file.path(plotDir, "Figure_S3A_heatmap_gene_set_zscore_heatmap.pdf"), width = 5, height = 9)
Figure_S3A_heatmap_gene_set_zscore_heatmap
dev.off()



