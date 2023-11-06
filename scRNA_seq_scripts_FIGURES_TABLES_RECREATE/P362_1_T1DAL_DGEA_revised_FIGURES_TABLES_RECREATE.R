#### P362-1 T1DAL Monocle Cluster Differential Expression Analysis ####

## Goals: better characterize the transcriptional signature differentiating CD57, PD1 and DN clusters 

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
library(KEGGREST)
library(Matrix)

setwd("/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup")
set.seed(42) # same seed as QC script 

#### LOAD ALL SAVED DATA AND SET PATHS ####

# Set Paths 
baseDir <- "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup"
plotDir <- file.path(baseDir,'FIGURES')
resultDir <- file.path(baseDir,'SAVED_DATA')
annotationDir <- file.path(baseDir, "RAW_DATA")

#### LOAD ALL SAVED DATA ####

# Load saved annotations 
load(file= file.path(resultDir, "P362-1_annotation.Rdata"))

# Load post QC no MAIT cell data
load(file= file.path(resultDir,"P362-1 T1DAL cds object - postQC no MAIT cells.Rdata"))

# remove cells in cluster 9 from cds object
cds_no_MAIT_no_9 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name != 9))]
#save(cds_no_MAIT_no_9, file = file.path(resultDir, "cds_no_MAIT_no_9.Rdata"))

load( file = file.path(resultDir, "cds_no_MAIT_no_9.Rdata"))

## load ITN dictionary matching operational donor IDs to the correct public IDs
T1DAL_ITN_ID_dictionary <- readxl::read_xlsx(file.path(annotationDir, "t1dal_mask_id_EW_formatted.xlsx"))
colnames(T1DAL_ITN_ID_dictionary) <- c("operational_PID", "masked_public_PID","Donor.ID")
# remove T1DAL from public_PID column
T1DAL_ITN_ID_dictionary$masked_public_PID <- str_remove(T1DAL_ITN_ID_dictionary$masked_public_PID, "T1DAL_")

donor_colors <- c("#bb5595", "#8fb03d", "#6677db", "#c99f3b", "#583687", "#66a453", "#c26ac6", "#45c097", "#b64468", "#848bd3", "#ab7434", "#bb4c41")

#### FIGURE S8: Find top marker genes expressed by each cluster with cluster 9 removed ####

# Use top_markers to ask what genes make clusters different from one another
# only re-run top_markers if necessary, it takes a long time to run
#marker_test_res <- top_markers(cds_no_MAIT_no_9, 
#                               group_cells_by="Cluster.Name")
# marker_test_res contains a number of metrics for how specifically expressed each gene is in each partition. 
# save
#save(marker_test_res, file = file.path(resultDir, "marker_test_res.Rdata"))

# load previously run marker test results
load(file.path(resultDir, "marker_test_res.Rdata"))

# export all markers that meet the specificity threshold
top_specific_markers <- marker_test_res %>%
  filter(specificity >= 0.25) %>%
  group_by(cell_group) 
top_specific_markers_ids <- unique(top_specific_markers %>% pull(gene_short_name))
top_specific_markers %>% dplyr::count(cell_group)

# still there are a few markers that are specific for more than 1 cluster

# Now we can plot the expression and fraction of cells that express each 
# marker in each group with the plot_genes_by_group function
top_markers_expression <- plot_genes_by_group(cds_no_MAIT_no_9,
                                              top_specific_markers_ids,
                                              group_cells_by="Cluster.Name",
                                              ordering_type= "maximal_on_diagonal",
                                              max.size=5) +
  theme(text = element_text(size = 16))
# FIGURE S8A
ggsave(top_markers_expression, file = file.path(plotDir, "top_markers_expression_no9.pdf"), height = 10, width = 8)

#### TABLE S3 SHEET 4: export results of top markers ####
write.csv(top_specific_markers, file = file.path(resultDir, "top_specific_markers_updated_9_26_2023.csv"))
top_specific_markers <- read.csv(file = file.path(resultDir, "top_specific_markers_updated_9_26_2023.csv"))
View(top_specific_markers)

top_specific_markers %>% filter(cell_group %in% c(5,6,8)) %>% View()


#### FIGURE 4A, TABLE S4: Analyze PD1 vs CD57 differential expression by cluster using regression analysis ####

# use the data subset to only include those clusters after cluster 4
cds_no_MAIT_ex <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,6,7,8)))]

# create dummy variable grouping 5,6,8 clusters and cluster 7 separately
cluster_ex_group <-  as.data.frame(colData(cds_no_MAIT_ex)) %>% mutate(cluster_ex = case_when(
  Cluster.Name %in% c(5,6,8) ~ "CD57",
  Cluster.Name == 7 ~ "PD1"
))

# join onto coldata
colData(cds_no_MAIT_ex)$cluster_ex <- cluster_ex_group$cluster_ex
levels(as.factor(colData(cds_no_MAIT_ex)$cluster_ex)) # "CD57" "PD1" - meaning that it was PD1 vs CD57 since first term is usually basis for comparison

# negative = higher in CD57
# positive = higher in PD1

##Re-Running the following section is very time consuming, just load significant results below

# test for genes that differ between cluster ex
#gene_fits_ex <- fit_models(cds_no_MAIT_ex , model_formula_str = "~cluster_ex")
#save(gene_fits_ex ,file = file.path(resultDir, "gene_fits_cluster_ex.Rdata"))
#load(file = file.path(resultDir, "gene_fits_cluster_ex.Rdata"))
#
## extract table of coefficients
#fit_coefs <- coefficient_table(gene_fits_ex)
#
## extract the cluster term
#cluster_ex_terms <- fit_coefs %>% filter(term != "(Intercept)") 
#levels(as.factor(cluster_ex_terms$term)) #  "cluster_exPD1"
#
## filter for those genes significantly differing between cluster 7 and cluster 5,6,8
#cluster_ex_terms_sig <- cluster_ex_terms   %>% filter (q_value < 0.05) %>%
#  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_ex_terms_sig, file = file.path(resultDir, "cluster_ex_terms_sig.Rdata"))
load(file = file.path(resultDir, "cluster_ex_terms_sig.Rdata"))
rm(gene_fits_ex)

# join category for up or down
cluster_ex_terms_sig_plot <- cluster_ex_terms_sig %>% mutate(up_down = case_when(normalized_effect >=0 ~"PD-1-like",
                                                                                 normalized_effect <0 ~ "CD57-like"))

# filter for up and down
cluster_ex_terms_sig_up <- cluster_ex_terms_sig  %>% filter(normalized_effect >=0)
cluster_ex_terms_sig_down <- cluster_ex_terms_sig  %>% filter(normalized_effect <0)

# negative = higher in CD57
# positive = higher in PD1

# TABLE S4: write.csv
write.csv(cluster_ex_terms_sig_up, file = file.path(resultDir, "cluster_ex_terms_sig_up.csv"))
write.csv(cluster_ex_terms_sig_down, file = file.path(resultDir, "cluster_ex_terms_sig_down.csv"))

# save table for plotting
#save(cluster_ex_terms_sig_plot, file= file.path(resultDir, "cluster_ex_terms_sig_plot.Rdata"))

### Check overlap of PD1 vs CD57 with atac lists

# Load ATAC-seq contrast results for those all genes and not just those that are overlapping the start site

ATAC_file_list = c("ATACseqData_P452_3_norm_db_anno_CD57pos_DN_up.csv",
                   "ATACseqData_P452_3_norm_db_anno_CD57minus_DN_up.csv",
                   "ATACseqData_P452_3_norm_db_anno_up.csv",
                   "ATACseqData_P452_3_norm_db_anno_down.csv")

read_atac <- function(i) {
  
  df <- read.csv(paste0(resultDir,"/",i)) %>%
    mutate(atac_list = i) %>% mutate(atac_list = str_remove(atac_list, "ATACseqData_P452_3_norm_db_anno_")) %>%
    mutate(atac_list = str_remove(atac_list, ".csv"))
  df
}
atac_all_up <- lapply(ATAC_file_list, read_atac)  
names(atac_all_up) <- ATAC_file_list
atac_all_up_df <- bind_rows(atac_all_up)

# join significant genes with the ATAC list 
cluster_ex_terms_sig_atac <- cluster_ex_terms_sig %>% dplyr::rename(gene_name = gene_short_name) %>%
  left_join(atac_all_up_df) %>% filter(!is.na(atac_list)) %>% mutate(up_down = case_when(
    normalized_effect >=0 ~ "up_in_PD1",
    normalized_effect <0 ~"up_in_CD57"
  ))
cluster_ex_terms_sig_atac %>% dplyr::count(atac_list, up_down )
# atac_list       up_down        n
# <chr>           <chr>      <int>
# 1 CD57minus_DN_up up_in_CD57    40
# 2 CD57minus_DN_up up_in_PD1     39
# 3 CD57pos_DN_up   up_in_CD57    73
# 4 CD57pos_DN_up   up_in_PD1     64
# 5 down            up_in_PD1     12
# 6 up              up_in_CD57     5

cluster_ex_terms_sig_atac %>% filter(atac_list == "down" & up_down == "up_in_PD1") %>% View()
cluster_ex_terms_sig_atac %>% filter(atac_list == "down" & up_down == "up_in_PD1") %>% dplyr::select(gene_name)
#gene_name
#<chr>    
#  1 RCAN3    
#2 RGS1     
#3 DPP4     
#4 INPP4B   
#5 IL7R     
#6 IL7R     
#7 THEMIS   
#8 SGK1     
#9 BEX3     
#10 JAML     
#11 HAPLN3   
#12 R3HDM4  

cluster_ex_terms_sig_atac %>% filter(atac_list == "up" & up_down == "up_in_CD57") %>% View()
cluster_ex_terms_sig_atac %>% filter(atac_list == "up" & up_down == "up_in_CD57") %>%  dplyr::select(gene_name)
#gene_name
#<chr>    
#  1 RNF220   
#2 KLRF1    
#3 S1PR5    
#4 S1PR5    
#5 KIR3DL2 

#### FIGURE 4A: Plotting PD1 vs CD57 DEGs using curated gene list ####

load(file= file.path(resultDir, "cluster_ex_terms_sig_plot.Rdata"))

curated_label <- data.frame(gene_short_name = c("IFIT1",
                                                "TBX21", "KLRD1", "KIR3DL2",  
                                                "GNLY", "GZMB",
                                                "PFN1",
                                                "PDCD1","GZMA","LAG3","TIGIT",
                                                "CTLA4","CD2",
                                                "AIF1","NKFBIA", "KLRF1","IL7R","THEMIS","TCF7","LEF1"))
curated_label$label <- curated_label$gene_short_name                       



cluster_ex_terms_sig_plot_curated <- left_join(cluster_ex_terms_sig_plot,curated_label)

# format the volcano plot to highlight all genes
cluster_ex_terms_sig_plot_curated_highlight <- ggplot(cluster_ex_terms_sig_plot_curated, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                             label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 4.5, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#93a24eff", "#ba4d4cc7")) + 
  theme(legend.position="none", text = element_text(size = 20)) +
  labs(x = "log2 Fold Change",y = "-log10 q-value")

# add annotation under the columns

text_low <- textGrob("Up in CD57+", gp=gpar(fontsize=18, fontface="bold"))
text_high <- textGrob("Up in PD-1+", gp=gpar(fontsize=18, fontface="bold"))

cluster_ex_terms_sig_plot_curated_highlight_plot <- cluster_ex_terms_sig_plot_curated_highlight + 
  annotation_custom(text_high,xmin=3,xmax=3,ymin=-45,ymax=-45) + 
  annotation_custom(text_low,xmin=-2,xmax=-2,ymin=-45,ymax=-45) +
  coord_cartesian( clip="off")

# FIGURE 4A: export plot
ggsave(plot = cluster_ex_terms_sig_plot_curated_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_curated_highlight_plot.pdf"),device = "pdf", width = 7, height = 6)


#### TABLE S7: Analyze differential expression of CD57 clusters versus all other clusters ####

load( file = file.path(resultDir, "cds_no_MAIT_no_9.Rdata"))

# create dummy variable grouping 5,6,8 clusters as CD57 and all others separately
cluster_group <-  as.data.frame(colData(cds_no_MAIT_no_9)) %>% mutate(cluster_CD57 = case_when(
  Cluster.Name %in% c(5,6,8) ~ "CD57",
  !Cluster.Name %in% c(5,6,8) ~ "other"
))

# join onto coldata
colData(cds_no_MAIT_no_9)$cluster_CD57 <- cluster_group$cluster_CD57
levels(as.factor(colData(cds_no_MAIT_no_9)$cluster_CD57)) # "CD57"  "other" - meaning CD57 is the reference
colData(cds_no_MAIT_no_9)$cluster_CD57 <- factor(colData(cds_no_MAIT_no_9)$cluster_CD57, levels= c("other","CD57"))
levels(as.factor(colData(cds_no_MAIT_no_9)$cluster_CD57)) # "other" "CD57"  - other is the reference

## Re-runnning this code is very time consuming, load significant objects below code

# test for genes that differ between CD57 and all other
#gene_fits <- fit_models(cds_no_MAIT_no_9 , model_formula_str = "~cluster_CD57")
#
## extract table of coefficients
#fit_coefs <- coefficient_table(gene_fits)
#
## extract the response term
#cluster_terms <- fit_coefs %>% filter(term != "(Intercept)") 
#levels(as.factor(cluster_terms$term)) #  
##save(cluster_terms, file = file.path(resultDir, "cluster_terms_CD57_other.Rdata"))
#
## filter for those genes significantly differing between CD57 and all other
#cluster_terms_sig <- cluster_terms  %>% filter (q_value < 0.05) %>%
#  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
## normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_terms_sig, file = file.path(resultDir, "cluster_terms_sig_CD57_other.Rdata"))

load(file = file.path(resultDir, "cluster_terms_sig_CD57_other.Rdata"))
#rm(gene_fits)

## TABLE S7 export
write.csv(cluster_terms_sig, file = file.path(resultDir, "cluster_terms_sig_CD57_other.csv"))

#### TABLE S5: Analyze differential expression by response in PD1 or CD57 clusters ####

### Assess Response differential expression in PD1 cluster

# use the data subset to only include those clusters after cluster 4
cds_no_MAIT_ex_PD1 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(7)))]

## Re-runnning this code is very time consuming, load significant objects below code

# test for genes that differ between response
#gene_fits_response_PD1 <- fit_models(cds_no_MAIT_ex_PD1 , model_formula_str = "~Response")
#
## extract table of coefficients
#fit_coefs_response_PD1 <- coefficient_table(gene_fits_response_PD1)
#
## extract the response term
#cluster_response_terms_PD1 <- fit_coefs_response_PD1 %>% filter(term != "(Intercept)") 
#levels(as.factor(cluster_response_terms_PD1$term)) #  "ResponseR"
#
## filter for those genes significantly differing between R and NR in PD1 cluster
#cluster_response_terms_PD1_sig <- cluster_response_terms_PD1   %>% filter (q_value < 0.05) %>%
#  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
## normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_response_terms_PD1_sig, file = file.path(resultDir, "cluster_response_terms_PD1_sig.Rdata"))

load(file = file.path(resultDir, "cluster_response_terms_PD1_sig.Rdata"))
#rm(gene_fits_response_PD1)

## TABLE S5B
write.csv(cluster_response_terms_PD1_sig, file = file.path(resultDir, "cluster_response_terms_PD1_sig.csv"))
nrow(cluster_response_terms_PD1_sig)

# Repeat for CD57
# use the data subset to only include CD57 clusters
cds_no_MAIT_ex_CD57 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,6,8)))]
#save(cds_no_MAIT_ex_CD57, file = file.path(resultDir, "cds_no_MAIT_ex_CD57.Rdata"))

## Re-runnning this code is very time consuming, load significant objects below code

# test for genes that differ between response
#gene_fits_response_CD57 <- fit_models(cds_no_MAIT_ex_CD57 , model_formula_str = "~Response")
#
## extract table of coefficients
#fit_coefs_response_CD57 <- coefficient_table(gene_fits_response_CD57)
#
## extract the response term
#cluster_response_terms_CD57 <- fit_coefs_response_CD57 %>% filter(term != "(Intercept)") 
#levels(as.factor(cluster_response_terms_CD57$term)) #  "ResponseR"
#
## filter for those genes significantly differing between R and NR in CD57 cluster
#cluster_response_terms_CD57_sig <- cluster_response_terms_CD57   %>% filter (q_value < 0.05) %>%
#  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
## normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_response_terms_CD57_sig, file = file.path(resultDir, "cluster_response_terms_CD57_sig.Rdata"))

load(file = file.path(resultDir, "cluster_response_terms_CD57_sig.Rdata"))
#rm(gene_fits_response_CD57)

## TABLE S5A
write.csv(cluster_response_terms_CD57_sig, file = file.path(resultDir, "cluster_response_terms_CD57_sig.csv"))

#### FIGURE S9: Plot volcano plots of Differential expression by response for PD1 and CD57 clusters ####

# load DEG results
load(file = file.path(resultDir, "cluster_response_terms_PD1_sig.Rdata"))
load(file = file.path(resultDir, "cluster_response_terms_CD57_sig.Rdata"))

R_NR_colors <- c("#46c19a","#6d80d8")

PD1_curated_label <- data.frame(gene_short_name = c("PSME2", "IFNG", "PSME1", "OASL", 
                                                    "IFITM1", "IRF1", "ISG20", "IFI6", "IFNG", 
                                                    "TNFSF10", "ICOS", "SELL", "GZMK", "BHLHE40", "CXCR3", "LAG3"))

PD1_curated_label <- PD1_curated_label %>% mutate(label = gene_short_name)                        

cluster_response_terms_PD1_sig <- cluster_response_terms_PD1_sig %>% mutate(up_down = case_when(normalized_effect >= 0 ~ "up",
                                                                                                normalized_effect <0 ~ "down"))
PD1_cluster_ex_terms_sig_plot_curated <- left_join(cluster_response_terms_PD1_sig, PD1_curated_label)

# format the volcano plot to highlight all genes
cluster_ex_terms_PD1_sig_plot_curated_highlight <- ggplot(PD1_cluster_ex_terms_sig_plot_curated, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                                     label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 4.5, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#719944",
                                "#9750a1")) +
  theme(text = element_text(size = 20),
        legend.position="none") +
  labs(x = "log2 Fold Change",y = "-log10 q-value", title = "PD-1+ Cluster")

# add annotation under the columns

text_low <- textGrob("Up in NR", gp=gpar(fontsize=18, fontface="bold"))
text_high <- textGrob("Up in R", gp=gpar(fontsize=18, fontface="bold"))

cluster_ex_terms_PD1_sig_plot_curated_highlight_plot <- cluster_ex_terms_PD1_sig_plot_curated_highlight  + 
  annotation_custom(text_high,xmin=1,xmax=1,ymin=-8,ymax=-8) + 
  annotation_custom(text_low,xmin=-3.5,xmax=-3.5,ymin=-8,ymax=-8) +
  coord_cartesian( clip="off")

# FIGURE S9B: export plot
ggsave(plot = cluster_ex_terms_PD1_sig_plot_curated_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_PD1_sig_plot_curated_highlight_plot.pdf"),device = "pdf", width = 6, height = 6)

## Repeat for CD57 cluster plot
CD57_curated_label <- data.frame(gene_short_name = c( "KIR3DL2", "KIR2DL3",
                                                      "KLRC1", "KLRB1", "EIF1","TMSB4X","TGFB1","EIF5A","PTMA","XIST","NKG7",
                                                      "NFKBIA", "KRAS","TNFSF10","KLRB1","TRAF3","MAP2K1"))
CD57_curated_label$label <-CD57_curated_label$gene_short_name

cluster_response_terms_CD57_sig <- cluster_response_terms_CD57_sig %>% mutate(up_down = case_when(normalized_effect >= 0 ~ "up",
                                                                                                  normalized_effect <0 ~ "down"))
CD57_cluster_ex_terms_sig_plot_curated <- left_join(cluster_response_terms_CD57_sig, CD57_curated_label)

# format the volcano plot to highlight all genes
cluster_ex_terms_CD57_sig_plot_curated_highlight <- ggplot(CD57_cluster_ex_terms_sig_plot_curated, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                                       label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 4.5, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = R_NR_colors) +
  theme( text = element_text(size = 20),
         legend.position="none") +
  labs(x = "log2 Fold Change",y = "-log10 q-value", title = "CD57+ Clusters")

# add annotation under the columns

text_low <- textGrob("Up in NR", gp=gpar(fontsize=16, fontface="bold"))
text_high <- textGrob("Up in R", gp=gpar(fontsize=16, fontface="bold"))

cluster_ex_terms_CD57_sig_plot_curated_highlight_plot <- cluster_ex_terms_CD57_sig_plot_curated_highlight  + 
  annotation_custom(text_high,xmin=2,xmax=2,ymin=-25,ymax=-25) + 
  annotation_custom(text_low,xmin=-2.5,xmax=-2.5,ymin=-25,ymax=-25) +
  coord_cartesian( clip="off")

# FIGURE S9A: export plot
ggsave(plot = cluster_ex_terms_CD57_sig_plot_curated_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_CD57_sig_plot_curated_highlight_plot.pdf"),device = "pdf", width = 6, height = 6)

#### FIGURE S9C: Replot volcano plot highlighting those genes differing in PD1 vs CD57 that are also different between R and NR...####

### Compare whether DEGs between R and NR account for additional heterogeneity ###

# load previously generated objects
load(file = file.path(resultDir, "cluster_response_terms_CD57_sig.Rdata"))
load(file = file.path(resultDir, "cluster_ex_terms_sig.Rdata"))
load(file= file.path(resultDir, "cluster_ex_terms_sig_plot.Rdata"))

cluster_response_terms_CD57_sig_compare <- cluster_response_terms_CD57_sig %>% dplyr::rename(gene_name = gene_short_name)

length(unique(cluster_response_terms_CD57_sig_compare$gene_name)) # 635

## Compare R vs NR list to CD57 vs PD1 lists # 264 overlaps between these lists
CD57_R_vs_NR_PD1_vs_CD57_overlap <- cluster_response_terms_CD57_sig[cluster_response_terms_CD57_sig$gene_short_name %in% unique(cluster_ex_terms_sig $gene_short_name),]
CD57_R_vs_NR_PD1_vs_CD57_overlap$overlap <- "yes"
cluster_response_terms_CD57_sig[cluster_response_terms_CD57_sig$gene_short_name %in% unique(cluster_ex_terms_sig $gene_short_name),] %>% View()


cluster_ex_terms_sig_plot_R_NR_overlap <- cluster_ex_terms_sig_plot %>% left_join(.,unique(CD57_R_vs_NR_PD1_vs_CD57_overlap[,c("gene_short_name","overlap")])) %>%
  mutate(up_down_overlap = case_when(overlap == "yes" ~ "R vs. NR\nOverlap",
                                     up_down == "PD-1-like"    ~ "PD-1 Tex",
                                     up_down == "CD57-like"    ~ "CD57 Tex"))

curated_label_new <- data.frame(gene_short_name = c("IFIT1",
                                                    "TBX21", "KLRD1", "KIR3DL2",  
                                                    "GNLY", "GZMB",
                                                    "PFN1",
                                                    "PDCD1","GZMA","LAG3","TIGIT",
                                                    "CTLA4","CD2",
                                                    "AIF1","NKFBIA", "KLRF1","IL7R","THEMIS","TCF7","LEF1"))
curated_label_new$label <- curated_label_new$gene_short_name                       


cluster_ex_terms_sig_plot_curated_overlap <- left_join(cluster_ex_terms_sig_plot_R_NR_overlap,curated_label_new)

# format the volcano plot to highlight all genes
cluster_ex_terms_sig_plot_curated_highlight_overlap <- ggplot(cluster_ex_terms_sig_plot_curated_overlap, aes(x = normalized_effect, y = -log10(q_value), color = up_down_overlap,
                                                                                                             label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 4.5, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#ba4d4cc7", "#fdbf6fff", "gray")) +
  theme( text = element_text(size = 20)) +
  labs(x = "log2 Fold Change",y = "-log10 q-value", colour = "")

# add annotation under the columns

text_low <- textGrob("Up in CD57+", gp=gpar(fontsize=18, fontface="bold"))
text_high <- textGrob("Up in PD-1+", gp=gpar(fontsize=18, fontface="bold"))

cluster_ex_terms_sig_plot_curated_highlight_plot_overlap <- cluster_ex_terms_sig_plot_curated_highlight_overlap + 
  annotation_custom(text_high,xmin=3,xmax=3,ymin=-45,ymax=-45) + 
  annotation_custom(text_low,xmin=-2,xmax=-2,ymin=-45,ymax=-45) +
  coord_cartesian( clip="off")

# FIGURE S9C: export plot
ggsave(plot = cluster_ex_terms_sig_plot_curated_highlight_plot_overlap, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_curated_highlight_plot_overlap.pdf"),device = "pdf", width = 9, height = 6)

#### TABLE S3 SHEET 1,2,3: Explore differential expression between CD57 clusters ####

# assess differential expression between cluster 5 and 6
cds_no_MAIT_ex_5_6 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,6)))]

## Re-runnning this code is very time consuming, load significant objects below code

## test for genes that differ between cluster ex
##gene_fits_ex_5_6 <- fit_models(cds_no_MAIT_ex_5_6 , model_formula_str = "~Cluster.Name")
##save(gene_fits_ex_5_6,file = file.path(resultDir, "gene_fits_ex_5_6.Rdata"))
##load(file = file.path(resultDir, "gene_fits_ex_5_6.Rdata"))
#
## extract table of coefficients
#fit_coefs_5_6 <- coefficient_table(gene_fits_ex_5_6)
#
## extract the cluster term
#cluster_ex_terms_5_6 <- fit_coefs_5_6 %>% filter(term != "(Intercept)") 
#levels(as.factor(cluster_ex_terms_5_6$term)) # "Cluster.Name6"
#
## filter for those genes significantly differing
#cluster_ex_terms_5_6_select <- cluster_ex_terms_5_6 %>% dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
#cluster_ex_term_5_6_sig <- cluster_ex_terms_5_6_select  %>% filter (q_value < 0.05) 
## normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_ex_terms_5_6_select, file = file.path(resultDir, "cluster_ex_terms_5_6.Rdata"))
#save(cluster_ex_term_5_6_sig, file = file.path(resultDir, "cluster_ex_term_5_6_sig.Rdata"))

load(file = file.path(resultDir, "cluster_ex_term_5_6_sig.Rdata"))

### TABLE S3 SHEET 3
write.csv(cluster_ex_term_5_6_sig, file = file.path(resultDir, "cluster_ex_term_6_vs_5_sig.csv"), row.names = FALSE)
#rm(gene_fits_ex_5_6)
nrow(cluster_ex_term_5_6_sig) # 1124

## Repeat for 5,8
cds_no_MAIT_ex_5_8 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,8)))]

## Re-runnning this code is very time consuming, load significant objects below code

## test for genes that differ between cluster ex
#gene_fits_ex_5_8 <- fit_models(cds_no_MAIT_ex_5_8 , model_formula_str = "~Cluster.Name")
##save(gene_fits_ex_5_8,file = file.path(resultDir, "gene_fits_ex_5_8.Rdata"))
##load(file = file.path(resultDir, "gene_fits_ex_5_8.Rdata"))
#
## extract table of coefficients
#fit_coefs_5_8 <- coefficient_table(gene_fits_ex_5_8)
#
## extract the cluster term
#cluster_ex_terms_5_8 <- fit_coefs_5_8 %>% filter(term != "(Intercept)") 
#levels(as.factor(cluster_ex_terms_5_8$term)) #  "Cluster.Name8"
#
## filter for those genes significantly differing 
#cluster_ex_terms_5_8_select <- cluster_ex_terms_5_8 %>%
#  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
#cluster_ex_term_5_8_sig <- cluster_ex_terms_5_8_select %>% filter (q_value < 0.05) 
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_ex_terms_5_8_select, file = file.path(resultDir, "cluster_ex_terms_5_8.Rdata"))
#save(cluster_ex_term_5_8_sig, file = file.path(resultDir, "cluster_ex_term_5_8_sig.Rdata"))
load(file = file.path(resultDir, "cluster_ex_term_5_8_sig.Rdata"))
#load(file = file.path(resultDir, "cluster_ex_term_5_8.Rdata"))

## TABLE S3 SHEET 2
write.csv(cluster_ex_term_5_8_sig, file = file.path(resultDir, "cluster_ex_term_8_vs_5_sig.csv"), row.names = FALSE)
#rm(gene_fits_ex_5_8)
nrow(cluster_ex_term_5_8_sig) # 458

## Repeat for 6,8

cds_no_MAIT_ex_6_8 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(6,8)))]

## Re-runnning this code is very time consuming, load significant objects below code

## test for genes that differ between cluster ex
#gene_fits_ex_6_8 <- fit_models(cds_no_MAIT_ex_6_8 , model_formula_str = "~Cluster.Name")
##save(gene_fits_ex_6_8,file = file.path(resultDir, "gene_fits_ex_6_8.Rdata"))
##load(file = file.path(resultDir, "gene_fits_ex_6_8.Rdata"))
#
## extract table of coefficients
#fit_coefs_6_8 <- coefficient_table(gene_fits_ex_6_8)
#
## extract the cluster term
#cluster_ex_terms_6_8 <- fit_coefs_6_8 %>% filter(term != "(Intercept)") 
#levels(as.factor(cluster_ex_terms_6_8$term)) #  "Cluster.Name8"
#
## filter for those genes significantly differing 
#cluster_ex_terms_6_8_select <- cluster_ex_terms_6_8 %>% 
#  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
#cluster_ex_term_6_8_sig <- cluster_ex_terms_6_8_select  %>% filter (q_value < 0.05) 
## normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
##save(cluster_ex_term_6_8_select, file = file.path(resultDir, "cluster_ex_term_6_8.Rdata"))
##save(cluster_ex_term_6_8_sig, file = file.path(resultDir, "cluster_ex_term_6_8_sig.Rdata"))

## TABLE S3 SHEET 1
load(file = file.path(resultDir, "cluster_ex_term_6_8_sig.Rdata"))
write.csv(cluster_ex_term_6_8_sig, file = file.path(resultDir, "cluster_ex_term_8_vs_6_sig.csv"), row.names = FALSE)
#rm(gene_fits_ex_6_8)

nrow(cluster_ex_term_6_8_sig)  # 466

#### FIGURE S5A,B PLOTTING ####

## Add correct public masked ID and replot without cluster 9!
load( file = file.path(resultDir,"cds_no_MAIT_no_9.Rdata"))
cds_no_MAIT_no_9_meta <- as.data.frame(colData(cds_no_MAIT_no_9))
cds_no_MAIT_no_9_meta <- left_join(cds_no_MAIT_no_9_meta, T1DAL_ITN_ID_dictionary)
colData(cds_no_MAIT_no_9)$masked_public_PID <- cds_no_MAIT_no_9_meta$masked_public_PID

# replot
cds_no_MAIT_UMAP_donor_no_9_PID <- plot_cells(cds_no_MAIT_no_9, color_cells_by = "masked_public_PID", label_cell_groups = F) + #+ facet_grid(~ Response) + 
  scale_color_manual(values = donor_colors) + 
  #labs(color = "Donor ID") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size=4),title='Donor ID')) 

# FIGURE S5A
ggsave(cds_no_MAIT_UMAP_donor_no_9_PID, file = file.path(plotDir, "cds_no_MAIT_UMAP_Donor.ID_23_01_23.tiff"), device = "tiff",
       width = 6, height = 5)

# replot with R vs NR
cds_no_MAIT_UMAP_response_9_group <- plot_cells(cds_no_MAIT_no_9, color_cells_by = "Response", label_cell_groups = F) + #+ facet_grid(~ Response) + 
  scale_color_manual(values = c("#56ae6c","#b0457b")) + 
  theme(text = element_text(size = 20)) 
# FIGURE S5A 
ggsave(cds_no_MAIT_UMAP_response_9_group, file = file.path(plotDir, "cds_no_MAIT_UMAP_response_23_09_22.tiff"), device = "tiff",
       width = 6, height = 5)

# plotting the number of cells per donor per cluster as a stacked bar plot
Donor_ID_cell_number <- as.data.frame(colData(cds_no_MAIT)) %>% filter(Cluster.Name != 9)
Donor_ID_cell_number$Cluster.Name <- factor(Donor_ID_cell_number$Cluster.Name, levels = c("8","7","6","5","4","3","2","1"))
total_cells_per_cluster <- Donor_ID_cell_number %>%
  dplyr::count(Cluster.Name) %>%
  ggplot( aes(x = Cluster.Name, y = n)) + geom_col() + 
  theme_bw() + labs(x = "Cluster", y = "Cell Number") + 
  theme(text = element_text(size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

percent_cells_per_donor_per_cluster <- Donor_ID_cell_number %>%
  group_by(Cluster.Name) %>%
  mutate(total_cells_cluster = n()) %>%
  ungroup() %>%
  group_by(Cluster.Name, Donor.ID) %>%
  mutate(total_per_cluster = n()) %>%
  distinct(Donor.ID, Cluster.Name, total_cells_cluster, total_per_cluster, Response) %>%
  mutate(percent_per_cluster = total_per_cluster/total_cells_cluster*100) 

percent_cells_per_donor_per_cluster_plot <- percent_cells_per_donor_per_cluster %>%
  left_join(.,T1DAL_ITN_ID_dictionary ) %>%
  ggplot(aes(x = Cluster.Name, y = percent_per_cluster, fill =masked_public_PID)) + geom_col() + 
  scale_fill_manual(values = donor_colors) + labs(x = "Cluster",y = "Cells per Donor (%)", fill = "Donor ID") + 
  theme_bw() + theme(text= element_text(size = 16))

percent_cells_per_response_per_cluster_plot <- percent_cells_per_donor_per_cluster %>%
  left_join(.,T1DAL_ITN_ID_dictionary ) %>%
  ggplot(aes(x = Cluster.Name, y = percent_per_cluster, fill =Response)) + geom_col() + 
  scale_fill_manual(values = c("#56ae6c","#b0457b") ) + labs(x = "Cluster",y = "Cells per Donor (%)", fill = "Response") + 
  theme_bw() + theme(text= element_text(size = 16))

# Plot the percent cells of each donor in each cluster
cluster_pal <- rev(c("#a6cee3" ,"#1f78b4", "#b2df8a", "#33a02c", "#fb9a99" ,"#e31a1c", "#fdbf6f", "#ff7f00"))

percent_cells_per_cluster_per_donor <- Donor_ID_cell_number %>%
  dplyr::count(Cluster.Name, Donor.ID, Response) %>% group_by(Donor.ID) %>% 
  mutate(total_per_donor = sum(n), percent_donor_cluster = n/total_per_donor*100) 

percent_cells_per_cluster_per_donor_plot <- percent_cells_per_cluster_per_donor %>%
  # fix donor ID to be public ID
  left_join(.,T1DAL_ITN_ID_dictionary ) %>%
  ggplot(aes(x = masked_public_PID, y = percent_donor_cluster, fill = Cluster.Name)) + geom_col() + 
  facet_grid(.~Response, scales = "free") +
  scale_fill_manual(values = cluster_pal) + labs(x = "Donor",y = "Cells per Cluster (%)", fill = "Cluster") + 
  theme_bw() + theme(text= element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1))

# FIGURE S5B
pdf(file = file.path(plotDir,"combined_cluster_percent_total.pdf"), height = 12, width = 7)
egg::ggarrange(total_cells_per_cluster, percent_cells_per_donor_per_cluster_plot, 
               percent_cells_per_response_per_cluster_plot,
               percent_cells_per_cluster_per_donor_plot,
               ncol = 1, heights = c(0.3,0.5, 0.5, 0.5)) 
dev.off()

