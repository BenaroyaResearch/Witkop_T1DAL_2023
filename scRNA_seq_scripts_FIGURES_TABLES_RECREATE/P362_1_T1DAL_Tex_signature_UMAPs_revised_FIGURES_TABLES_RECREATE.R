#### P362_1_T1DAL_Tex_signature_UMAPs.R ####

#### Load Libraries ####

library(tidyverse)
library(monocle3)
library(Matrix)
library(ggpubr)
library(gridExtra)
library(ComplexHeatmap)
library(readxl)
library(ggrepel)

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

#### FIGURE 3D, FIGURE S7A: Plot log 10 gene expression of particular gene lists on UMAP ####

# Plot marker expression to show Tex populations 

## List potential gene sets to analyze
# Kirsten did not have annotation for what these gene sets are
gene_list <- c("TIGIT", "KLRG1","PDCD1","TBX21", "TOX","CTLA4","LAG3","CD160","CD244","HAVCR2")

# rename data to prevent altering original data
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,gene_list]))
cds_exprs[cds_exprs == -Inf] <- NA 
colData(cds_data) <- cbind(colData(cds_data),cds_exprs)

## PLOT FIGURE 3D
pdf(file.path(plotDir, "Tex_population_ID_marker_expression_UMAP.pdf"), height = 15, width = 14)
all_zscore_plots <- list()

for(i in 1:length(gene_list)){
  order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(!is.na(get(gene_list[i])),get(gene_list[i]))
  cds_all_ordered <- cds_data[,row.names(order_cds)]
  all_zscore_plots[[i]] <- plot_cells(cds_all_ordered, color_cells_by =gene_list[i],show_trajectory_graph = F,cell_size=1.5) + 
    theme(text = element_text(size = 20))
}

grid.arrange(grobs=all_zscore_plots,ncol=3)
dev.off()


# Repeat with KIR gene set
gene_list <- c("KIR3DL3","KIR2DL3", "KIR2DL1", "KIR2DL4", "KIR3DL1","KIR2DS4","KIR3DL2")

# rename data to prevent altering original data
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,gene_list]))
cds_exprs[cds_exprs == -Inf] <- NA 
colData(cds_data) <- cbind(colData(cds_data),cds_exprs)

## PLOT FIGURE S7A
pdf(file.path(plotDir, "Tex_population_KIR_expression_UMAP.pdf"), height = 10, width = 14)
all_zscore_plots <- list()

for(i in 1:length(gene_list)){
  order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(!is.na(get(gene_list[i])),get(gene_list[i]))
  cds_all_ordered <- cds_data[,row.names(order_cds)]
  all_zscore_plots[[i]] <- plot_cells(cds_all_ordered, color_cells_by =gene_list[i],show_trajectory_graph = F,cell_size=1.5) + 
    theme(text = element_text(size = 20))
}

grid.arrange(grobs=all_zscore_plots,ncol=3)
dev.off()


#### FIGURE S7B: Plot Tex gene set mean as heatmap ####

cds_data <- cds_no_MAIT_no_9
gene_list <- c("TIGIT", "KLRG1","PDCD1","TBX21", "TOX","CTLA4","LAG3","CD160","CD244","HAVCR2")
## Get subset of log transformed cds data with genes in set
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data))[,gene_list])))
cds_exprs[cds_exprs == -Inf] <- NA
# cds_exprs <- scale(cds_exprs)

## Format expression data to get gene set mean per cell
cds_subset_genes <- cbind(colData(cds_data),cds_exprs) %>%
  as.data.frame() %>%
  mutate(gene_mean = rowMeans(dplyr::select(.,!(colnames(colData(cds_data)))),na.rm=T)) # had to add explicitly call dplyr

## Add per-cell gene set means back to cds colData for plotting

colData(cds_data) <- cbind(colData(cds_data),  cds_subset_genes)

# calculate mean across cluster
gene_set_mean_vals <- colData(cds_data) %>% as.data.frame() %>% group_by(Cluster.Name) %>% 
  summarise_at(vars(gene_mean),mean, na.rm=T) 
row.names(gene_set_mean_vals) <- gene_set_mean_vals$Cluster.Name

# Use scale to again get the zscore - but using the median
plot_scaled_mean <- scale(gene_set_mean_vals[,-1]) %>% as.data.frame() 
row.names(plot_scaled_mean) <- gene_set_mean_vals$Cluster.Name
colnames(plot_scaled_mean) <- "Tex Gene Set"

## FIGURE S7B
# Plot the scaled heatmap - using median
pdf(file=file.path(plotDir, "Tex_zscore_gene_set_heatmap.pdf"), width = 2, height = 3)
Heatmap(as.matrix(plot_scaled_mean),
        name="Zscore\nGene Set\nexprs",
        show_parent_dend_line=F, 
        show_row_names = T, 
        cluster_rows = T, cluster_columns = T, 
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#### FIGURE S15D: PLOT GENE LISTS FROM OTHER PUBLISHED DATASETS ONTO UMAPs ####

# THE UMAPS FOR FIGURE S16D WERE SAVED SEPARATELY AND COMPILED FOR THE FIGURE, the first part of the figure was copied from Figure  3B in the main text

## load Li et al 2022 Science non-exhausted KIR gene set Table S3, sheet 1 that was a result of each cluster vs all other clusters
Li_KIR_non_ex_CD57 <- read_excel(file.path(annotationDir, "science.abi9591_table_s3.xlsx"),
                                 sheet=1) %>% 
  # filter for genes in the KIR+ effector cell set versus all other clusters
  filter(cluster == "KIR+ effector CD8") %>%
  # filter for those increased specifically in the KIR+ effector set in the MS patients
  filter(MS_avg_logFC > 0 )

colnames( Li_KIR_non_ex_CD57)[1] <- "gene_name"

# get just gene names to aggregate  
Li_KIR_non_ex_CD57_genes <-  Li_KIR_non_ex_CD57$gene_name

## Load Giles et al...Wherry Nat Immun 2022 exhausted dataset from Giles_Wherry_2022_Supplementary_Table_1_DEGs_by_cluster.xlsx
# Use the tab regarding Cl13 data, which is from the chronic infection. DEGs by cluster. DEGs were calculated with using Seurat FindAllMarkers two-sided Wilcoxon test using Bonferroni correction.
# this table lists both positive and negative DEGs, I want to only look at those that are up in their exhausted KLR gene set

Giles_CD57_Tex <- read_excel(file.path(annotationDir, "Giles_Wherry_2022_Supplementary_Table_1_DEGs_by_cluster.xlsx"),
                             # get sheet 3 only 
                             sheet = 3) %>%
  # filter for exhausted KLR gene set
  filter(cluster == "Exh-KLR") %>%
  # filter for only those genes that are increased in this comparison
  filter(avg_log2FC >0)
Giles_CD57_Tex_genes <-  Giles_CD57_Tex$gene
# make uppercase to convert to human
Giles_CD57_Tex_genes <- str_to_upper( Giles_CD57_Tex_genes)
length(Giles_CD57_Tex_genes) #180
# filter for genes that are in the rownames of our data
Giles_CD57_Tex_genes <-  Giles_CD57_Tex_genes[Giles_CD57_Tex_genes %in% rownames( cds_no_MAIT_no_9)] 
length(Giles_CD57_Tex_genes) # 153

## load Daniel et al, 2022 Nat Immun exhausted dataset https://www.nature.com/articles/s41590-022-01337-5 
Daniel_CD57_Tex <- read_excel(file.path(annotationDir, "Daniel_science_41590_2022_1337_MOESM2_ESM.xlsx"), 
                              # sheet one, supplementary table 1 lists the specific DEGs for each population
                              sheet =1) %>%
  # filter for TexKLR
  filter(cluster == "TexKLR") 
# log fold change is already positive here so I don't need to further filter
Daniel_CD57_Tex  <- str_to_upper( Daniel_CD57_Tex$gene)
length( Daniel_CD57_Tex) # 138
# filter for genes that are in the rownames of our data
Daniel_CD57_Tex_genes <-   Daniel_CD57_Tex[ Daniel_CD57_Tex %in% rownames( cds_no_MAIT_no_9)] 
length( Daniel_CD57_Tex_genes) # 126

## Load EOMES module from Alice Long paper - supplementary table 2 PMID: 28664195
EOMES_mod <- read_excel(file.path(annotationDir, "aai7793_Table S2.xlsx"))
# these are the top 800 genes associated with EOMES expression
EOMES_mod <- EOMES_mod$`Table S2. EOMES-associated genes, the top 800 genes correlated with EOMES expression, as described in Linsley et al PLoS One. 2014 Oct 14;9(10):e109760.`
EOMES_mod  <- EOMES_mod [ EOMES_mod  %in% rownames( cds_no_MAIT_no_9)] 

## Load canonical cancer Tex gene set from Zheng et al PMID: 34914499 
# load signature genes of meta clusters from supplementary Table 3, using sheet 2 that lists the CD8 cell meta-clusters
Zheng_S3_CD8 <- read_excel(file.path(annotationDir, "science.abe6474_table_s3.xlsx"), sheet = 2, skip = 1)
unique(Zheng_S3_CD8$cluster.name)
# extract unique signatures
Zheng_S3_CD8_KIR_EOMES_NK_like <- Zheng_S3_CD8 %>% filter(cluster.name == "CD8.c08(KIR+EOMES+ NK-like)")
Zheng_S3_CD8_KIR_EOMES_NK_like <- Zheng_S3_CD8_KIR_EOMES_NK_like$geneSymbol

Zheng_S3_CD8_terminal_Tex <-  Zheng_S3_CD8 %>% filter(cluster.name == "CD8.c12(terminal Tex)" )
Zheng_S3_CD8_terminal_Tex <- Zheng_S3_CD8_terminal_Tex$geneSymbol
Zheng_S3_CD8_terminal_Tex <- Zheng_S3_CD8_terminal_Tex[Zheng_S3_CD8_terminal_Tex %in% rownames( cds_no_MAIT_no_9)]

Zheng_S3_CD8_GZMK_Tex  <- Zheng_S3_CD8 %>% filter(cluster.name == "CD8.c11(GZMK+ Tex)")
Zheng_S3_CD8_GZMK_Tex <- Zheng_S3_CD8_GZMK_Tex$geneSymbol
Zheng_S3_CD8_GZMK_Tex <- Zheng_S3_CD8_GZMK_Tex[Zheng_S3_CD8_GZMK_Tex %in% rownames( cds_no_MAIT_no_9)]

Zheng_S3_CD8_TCF7_Tex  <- Zheng_S3_CD8 %>% filter(cluster.name == "CD8.c14(TCF7+ Tex)" )
Zheng_S3_CD8_TCF7_Tex <- Zheng_S3_CD8_TCF7_Tex$geneSymbol

Zheng_S3_CD8_KIR_TXK_NK_like <- Zheng_S3_CD8 %>% filter(cluster.name =="CD8.c09(KIR+TXK+ NK-like)" )
Zheng_S3_CD8_KIR_TXK_NK_like <- Zheng_S3_CD8_KIR_TXK_NK_like$geneSymbol
Zheng_S3_CD8_KIR_TXK_NK_like <- Zheng_S3_CD8_KIR_TXK_NK_like[Zheng_S3_CD8_KIR_TXK_NK_like%in% rownames( cds_no_MAIT_no_9)]

### Create UMAPs that plot the aggregated gene signature of each 

### Plot Li et al. gene set on UMAP
# rename data to prevent altering original data
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,   Li_KIR_non_ex_CD57_genes]))
cds_exprs[cds_exprs == -Inf] <- NA 
# sum across all genes and normalize by total set size
cds_exprs$gene_set_total <-  rowSums( cds_exprs, na.rm = TRUE)/length(Li_KIR_non_ex_CD57_genes)
cds_exprs_total <-  cds_exprs %>% dplyr::select(gene_set_total) 
cds_exprs_total_Li <- cds_exprs_total
colData(cds_data) <- cbind(colData(cds_data), cds_exprs_total)
colnames( colData(cds_data))

# order with highest values at the top
order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(gene_set_total)
cds_all_ordered <- cds_data[,row.names(order_cds)]

# plot as UMAP
colnames(colData(cds_all_ordered))[15] <- "Li et al., KIR+\nCD57 Gene Set"
Li_KIR_UMAP_gene_set <- plot_cells( cds_all_ordered , color_cells_by ="Li et al., KIR+\nCD57 Gene Set", show_trajectory_graph = F,cell_size=1.5) + 
  theme(text = element_text(size = 20)) +
  viridis::scale_color_viridis(option = "magma") +
  #limits = c(0,0.3)) +
  labs(colour = "Normlized Li KIR+\nGene Set")

ggsave(Li_KIR_UMAP_gene_set,file = file.path(plotDir, "Li_KIR_UMAP_gene_set.pdf"), height = 5, width = 7)

### Plot Giles CD57 Tex gene set on UMAP

# rename data to prevent altering original data
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,  Giles_CD57_Tex_genes]))
cds_exprs[cds_exprs == -Inf] <- NA 
# sum across all genes and normalize by set size
cds_exprs$gene_set_total <-  rowSums( cds_exprs, na.rm = TRUE)/length(Giles_CD57_Tex_genes)
cds_exprs_total <-  cds_exprs %>% dplyr::select(gene_set_total) 
cds_exprs_total_Giles <- cds_exprs_total
colData(cds_data) <- cbind(colData(cds_data), cds_exprs_total)
colnames( colData(cds_data))

# order with highest values at the top
order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(gene_set_total)
cds_all_ordered <- cds_data[,row.names(order_cds)]

# plot as UMAP
colnames(colData(cds_all_ordered))[15] <- "Giles Exh-KLR\nGene Set"
Giles_Exh_KLR_UMAP_gene_set <- plot_cells( cds_all_ordered , color_cells_by ="Giles Exh-KLR\nGene Set", show_trajectory_graph = F,cell_size=1.5) + 
  theme(text = element_text(size = 20)) +
  viridis::scale_color_viridis(option = "magma") + 
  #  limits = c(0,0.3)) +
  labs(colour = "Normalized Giles\nExh-KLR\nGene Set")

ggsave(Giles_Exh_KLR_UMAP_gene_set,file = file.path(plotDir, "Giles_Exh_KLR_Gene_Set_UMAP.pdf"), height = 5, width = 7)

### Plot with Daniel Tex-KLR dataset

# rename data to prevent altering original data
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,   Daniel_CD57_Tex_genes]))
cds_exprs[cds_exprs == -Inf] <- NA 
# sum across all genes normalize by set size 
cds_exprs$gene_set_total <-  rowSums( cds_exprs, na.rm = TRUE)/length(Daniel_CD57_Tex_genes)
cds_exprs_total <-  cds_exprs %>% dplyr::select(gene_set_total) 
cds_exprs_total_Daniel <- cds_exprs_total
colData(cds_data) <- cbind(colData(cds_data), cds_exprs_total)
colnames( colData(cds_data))

# order with highest values at the top
order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(gene_set_total)
cds_all_ordered <- cds_data[,row.names(order_cds)]

# plot as UMAP
colnames(colData(cds_all_ordered))[15] <- "Daniel TexKLR\nGene Set"
Daniel_TexKLR_UMAP_gene_set <- plot_cells( cds_all_ordered , color_cells_by ="Daniel TexKLR\nGene Set", show_trajectory_graph = F,cell_size=1.5) + 
  theme(text = element_text(size = 20)) +
  viridis::scale_color_viridis(option = "magma")+
  # limits = c(0,0.3)) +
  labs(colour = "Normalized Daniel\nTexKLR\nGene Set")

ggsave(Daniel_TexKLR_UMAP_gene_set,file = file.path(plotDir, "Daniel_TexKLR_UMAP_gene_set.pdf"), height = 5, width = 7)

### Plot EOMES mod 
# rename data to prevent altering original data
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,   EOMES_mod]))
cds_exprs[cds_exprs == -Inf] <- NA 
# sum across all genes normalize by set size 
cds_exprs$gene_set_total <-  rowSums( cds_exprs, na.rm = TRUE)/length(EOMES_mod)
cds_exprs_total <-  cds_exprs %>% dplyr::select(gene_set_total) 
colData(cds_data) <- cbind(colData(cds_data), cds_exprs_total)
colnames( colData(cds_data))

# order with highest values at the top
order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(gene_set_total)
cds_all_ordered <- cds_data[,row.names(order_cds)]

# plot as UMAP
colnames(colData(cds_all_ordered))[15] <- "EOMES Module\nGene Set"
EOMES_UMAP_gene_set <- plot_cells( cds_all_ordered , color_cells_by ="EOMES Module\nGene Set", show_trajectory_graph = F,cell_size=1.5) + 
  theme(text = element_text(size = 20)) +
  viridis::scale_color_viridis(option = "magma") + 
  # limits = c(0,0.3)) +
  labs(colour = "Normalized EOMES\nModule\nGene Set")

ggsave(EOMES_UMAP_gene_set,file = file.path(plotDir, "EOMES_UMAP_gene_set.pdf"), height = 5, width = 7)

### Plot Terminal Tex Zhang module
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,   Zheng_S3_CD8_terminal_Tex]))
cds_exprs[cds_exprs == -Inf] <- NA 
# sum across all genes normalize by set size 
cds_exprs$gene_set_total <-  rowSums( cds_exprs, na.rm = TRUE)/length(Zheng_S3_CD8_terminal_Tex)
cds_exprs_total <-  cds_exprs %>% dplyr::select(gene_set_total) 
colData(cds_data) <- cbind(colData(cds_data), cds_exprs_total)
colnames( colData(cds_data))

# order with highest values at the top
order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(gene_set_total)
cds_all_ordered <- cds_data[,row.names(order_cds)]

# plot as UMAP
colnames(colData(cds_all_ordered))[15] <- "Zheng Terminal Tex\nGene Set"
Zheng_terminal_tex_UMAP_gene_set <- plot_cells( cds_all_ordered , color_cells_by ="Zheng Terminal Tex\nGene Set", show_trajectory_graph = F,cell_size=1.5) + 
  theme(text = element_text(size = 20)) +
  viridis::scale_color_viridis(option = "magma") +
  #limits = c(0,0.3)) +
  labs(colour = "Normalized Zheng\nTerminal Tex\nGene Set")

ggsave(Zheng_terminal_tex_UMAP_gene_set ,file = file.path(plotDir, "Zheng_terminal_tex_UMAP_gene_set.pdf"), height = 5, width = 7)

### Repeat with Zhang GZMK Tex
cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,   Zheng_S3_CD8_GZMK_Tex]))
cds_exprs[cds_exprs == -Inf] <- NA 
# sum across all genes normalize by set size 
cds_exprs$gene_set_total <-  rowSums( cds_exprs, na.rm = TRUE)/length(Zheng_S3_CD8_GZMK_Tex)
cds_exprs_total <-  cds_exprs %>% dplyr::select(gene_set_total) 
colData(cds_data) <- cbind(colData(cds_data), cds_exprs_total)
colnames( colData(cds_data))

# order with highest values at the top
order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(gene_set_total)
cds_all_ordered <- cds_data[,row.names(order_cds)]

# plot as UMAP
colnames(colData(cds_all_ordered))[15] <- "Zheng GZMK+ Tex\nGene Set"
Zheng_GZMK_tex_UMAP_gene_set <- plot_cells( cds_all_ordered , color_cells_by ="Zheng GZMK+ Tex\nGene Set", show_trajectory_graph = F,cell_size=1.5) + 
  theme(text = element_text(size = 20)) +
  viridis::scale_color_viridis(option = "magma") +
  #limits = c(0,0.3)) +
  labs(colour = "Normalized Zheng\nGZMK+ Tex\nGene Set")

ggsave(Zheng_GZMK_tex_UMAP_gene_set ,file = file.path(plotDir, "Zheng_GZMK_tex_UMAP_gene_set.pdf"), height = 5, width = 7)

### Repeat with Zhang KIR TXK

cds_data <- cds_no_MAIT_no_9

# Log10 scale expression and subset for gene list
cds_exprs <- log10(as.data.frame(as.matrix(t(exprs(cds_data)))[,   Zheng_S3_CD8_KIR_TXK_NK_like]))
cds_exprs[cds_exprs == -Inf] <- NA 
# sum across all genes normalize by set size 
cds_exprs$gene_set_total <-  rowSums( cds_exprs, na.rm = TRUE)/length(Zheng_S3_CD8_KIR_TXK_NK_like)
cds_exprs_total <-  cds_exprs %>% dplyr::select(gene_set_total) 
colData(cds_data) <- cbind(colData(cds_data), cds_exprs_total)
colnames( colData(cds_data))

# order with highest values at the top
order_cds <- colData(cds_data) %>% as.data.frame() %>% arrange(gene_set_total)
cds_all_ordered <- cds_data[,row.names(order_cds)]

# plot as UMAP
colnames(colData(cds_all_ordered))[15] <- "Zheng_KIR_TXK_NK_like_Gene_Set"
Zheng_KIR_TXK_NK_like_UMAP_gene_set <- plot_cells( cds_all_ordered , color_cells_by ="Zheng_KIR_TXK_NK_like_Gene_Set", show_trajectory_graph = F,cell_size=1.5) + 
  theme(text = element_text(size = 20)) +
  viridis::scale_color_viridis(option = "magma") +
  #  limits = c(0,0.3)) +
  labs(colour = "Normalized Zheng\nKIR+TXK+ NK-like\nGene Set")

ggsave(Zheng_KIR_TXK_NK_like_UMAP_gene_set ,file = file.path(plotDir, "Zheng_KIR_TXK_NK_like_UMAP_gene_set.pdf"), height = 5, width = 7)

# AGAIN: ALL UMAPS ABOVE WERE COMPILED SEPARELY INTO A UMAP FOR FIGURE S16D

#### FIGURE S15C: Assess overlap between DGEA genes in CD57 vs all other clusters and loaded gene sets ####

# load signifcant genes from the CD57 analysis - this is all genes up and down
cluster_terms_sig_CD57 <- read.csv(file = file.path(resultDir, "cluster_terms_sig_CD57_other.csv")) %>% dplyr::select(-X)

# filter for top 1000 up genes 
cluster_terms_sig_CD57_up <- cluster_terms_sig_CD57 %>% filter(normalized_effect > 0) %>% top_n(normalized_effect,n=1000)
nrow(cluster_terms_sig_CD57_up) # 1000
cluster_terms_sig_CD57_up_genes <- cluster_terms_sig_CD57_up$gene_short_name

# make list of gene sets to check overlaps of 
gene_sets <- list("Zheng_S3_CD8_KIR_TXK_NK_like" = Zheng_S3_CD8_KIR_TXK_NK_like,
                  "Li_KIR_non_ex_CD57_genes" =Li_KIR_non_ex_CD57_genes,
                  "Giles_CD57_Tex_genes" = Giles_CD57_Tex_genes,
                  "Daniel_CD57_Tex_genes" = Daniel_CD57_Tex_genes,
                  "EOMES_mod" = EOMES_mod,
                  "Zheng_S3_CD8_terminal_Tex" = Zheng_S3_CD8_terminal_Tex,
                  "Zheng_S3_CD8_GZMK_Tex" =Zheng_S3_CD8_GZMK_Tex)
gene_sets_lengths <- lengths(gene_sets )
# variable lengths 
length(cluster_terms_sig_CD57_up_genes[cluster_terms_sig_CD57_up_genes %in% c(gene_sets[[1]])])

Tex_phyper_up <- data.frame()
for (row in 1:length(gene_sets)) {
  # group1 = our data
  # group2 = any set
  total = nrow(cds_no_MAIT_no_9) #which was the number of background genes in our data set
  intersect <- length(cluster_terms_sig_CD57_up_genes[cluster_terms_sig_CD57_up_genes %in% c(gene_sets[[row]])])
  pub_set_length <- gene_sets_lengths[[row]]
  print(pub_set_length)
  phyper_df <- phyper(intersect - 1, pub_set_length, (total - pub_set_length), length(cluster_terms_sig_CD57_up_genes), lower.tail = FALSE )
  phyper_df <- cbind(set  = names(gene_sets_lengths[row]), intersection = intersect, pub_length = pub_set_length, CD57_length =length(cluster_terms_sig_CD57_up_genes),
                     as.data.frame(phyper_df))
  print(phyper_df )
  Tex_phyper_up <- rbind(phyper_df, Tex_phyper_up)
  
}

View(Tex_phyper_up)

# add set group
Tex_phyper_up <- Tex_phyper_up %>% mutate(`Set Type` = case_when(
  set %in% c("Zheng_S3_CD8_terminal_Tex", "Zheng_S3_CD8_GZMK_Tex","EOMES_mod","Daniel_CD57_Tex_genes", "Giles_CD57_Tex_genes") ~ "Tex",
  set %in% c("Li_KIR_non_ex_CD57_genes","Zheng_S3_CD8_KIR_TXK_NK_like") ~"Non-Tex"
)) %>%
  mutate(labels = case_when(
    set == "Zheng_S3_CD8_GZMK_Tex"  ~ "Zheng CD8 GZMK Tex",
    set == "Zheng_S3_CD8_terminal_Tex" ~ "Zheng CD8 Terminal Tex",
    set == "EOMES_mod" ~ "EOMES Module",
    set == "Daniel_CD57_Tex_genes" ~ "Daniel Tex-KLR",
    set == "Giles_CD57_Tex_genes" ~ "Giles Exh-KLR",
    set == "Li_KIR_non_ex_CD57_genes" ~ "Li KIR+ Non-Ex",
    set == "Zheng_S3_CD8_KIR_TXK_NK_like" ~ "Zheng CD8 KIR TXK NK-like"
  ))

# Plot as a bubble plot
Tex_phyper_up_bubble <-
  ggplot(Tex_phyper_up, aes(y = intersection, x = -log10(phyper_df), label = labels, color = `Set Type`)) + geom_point(size = 2) +
  #scale_colour_gradientn(colours=c("gray","orange","red")) + 
  ggrepel::geom_text_repel(data = subset(Tex_phyper_up, set != "EOMES_mod"),
                           color = "black", size = 3) +
  ggrepel::geom_text_repel(data = subset(Tex_phyper_up, set == "EOMES_mod"),
                           color = "black", fontface = "bold", size = 4) +
  labs(x = "-log10 Hypergeometric P-value", y = "# Shared Genes with CD57+ Tex") + 
  theme_minimal() +
  scale_color_manual(values = c("#b84c7d",
                                "#9da140")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

## FIGURE S15C
ggsave(Tex_phyper_up_bubble, file = file.path(plotDir,"Tex_phyper_up_bubble.pdf"), height = 4, width = 4.5)

