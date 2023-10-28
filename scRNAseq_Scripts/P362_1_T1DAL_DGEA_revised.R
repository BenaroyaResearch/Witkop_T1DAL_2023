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

setwd("/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/")

set.seed(42) # same seed as QC script 

#### LOAD ALL SAVED DATA ####

# Load saved annotations (Kirsten erroneously called then P362-2)
load("./P362-1_annotation.Rdata")

# Load post QC no MAIT cell data
load("./EW_T1DAL_Results/P362-1 T1DAL cds object - postQC no MAIT cells.Rdata")

# set output directories
plotDir <- "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/FIGURES/DGEA"
resultDir <- "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/EW_T1DAL_Results"

# remove cells in cluster 9 from cds object
cds_no_MAIT_no_9 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name != 9))]
#save(cds_no_MAIT_no_9, file = file.path(resultDir, "cds_no_MAIT_no_9.Rdata"))

load( file = file.path(resultDir, "cds_no_MAIT_no_9.Rdata"))

#### Monocle DGEA Graph Autocorrelation Analysis ####

## This original code was edited from code in "P362_1_T1DAL_10X_Compiled_Analysis.R" 
## The graph-autocorrelation test shown below allows the user to identify genes that
# vary across clusters
# https://cole-trapnell-lab.github.io/monocle3/docs/differential/

# create new monocle object to ensure I don't mess up old one
subset_for_DGEA <- cds_no_MAIT_no_9
subset_for_DGEA <- reduce_dimension(subset_for_DGEA)

## Run monocle graph test for DEG analysis
# The function graph_test() uses the statistic Moran's I, which Cao & Spielmann et al showed to be effective
# in finding genes that vary in single-cell RNA-seq datasets.
# the Moran's I value indicates a gene's effect size, and genes with a  positive value are associated with a specific cluster

## graph test takes a very long time to run, load saved data in the future
graph_test_result_subset <- graph_test(subset_for_DGEA, neighbor_graph="knn", cores = 8)
# save result for later, preserving Kirsten's original name 
#save(graph_test_result_subset, file = "./EW_T1DAL_Results/P362-1_T1DAL_DGEA_rerun.Rdata")

## Load results of graph test previously computed since this code took a while
load(file = "./EW_T1DAL_Results/P362-1_T1DAL_DGEA_rerun.Rdata")

# Subset for significant genes 
subset_DEG_IDs <- row.names(subset(graph_test_result_subset, q_value < 0.05))
length(subset_DEG_IDs) # 8042

## Find modules of co-regulated genes
# find_gene_modules() essentially runs UMAP on the genes (as opposed to the cells) 
# and then groups them into modules using Louvain community analysis
# modules at this point are not related to previous clusters

gene_module_df <- find_gene_modules(subset_for_DGEA[subset_DEG_IDs,], resolution=0.001)

# Plot expression of the co-regulated gene modules with significant genes
coregulated_gene_modules <- plot_cells(subset_for_DGEA,
                                       genes=gene_module_df,
                                       label_cell_groups=TRUE,
                                       show_trajectory_graph=FALSE,
                                       cell_size=2)
ggsave(coregulated_gene_modules, file = file.path(plotDir, "coregulated_gene_module_expression.pdf"), height= 12, width = 15)


# plot cells without modules 3, 10, 1 removed since these have much higher expression scores and make other patterns difficult to see 
coregulated_gene_modules_3_10_1_rm <- plot_cells(subset_for_DGEA,
                                                 genes=gene_module_df[!(gene_module_df$module %in% c(3,10,1)),],
                                                 label_cell_groups=TRUE,
                                                 show_trajectory_graph=FALSE,
                                                 cell_size=2)
ggsave(coregulated_gene_modules_3_10_1_rm , file = file.path(plotDir, "coregulated_gene_module_expression_3_10_1_rm.pdf"), height= 12, width = 15)

### Plot module heatmap by cluster - to see which module genes are expressed in which clusters

cell_group_df <- tibble::tibble(cell=row.names(colData(subset_for_DGEA)), 
                                cell_group=colData(subset_for_DGEA)$Cluster.Name)

agg_mat <- aggregate_gene_expression(subset_for_DGEA, gene_module_df_0, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pdf(file =file.path(plotDir, "graph_auto_DEG_clusters.pdf"))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
dev.off()

#### Get module gene lists and run GO analysis ####
# Example: https://ucdavis-bioinformatics-training.github.io/2019_August_UCD_mRNAseq_Workshop/differential_expression/enrichment

## Get gene lists for each module
module_list <- unique(gene_module_df$module)
get_module_gene_list <- function(i, data) {
  data %>% dplyr::filter(module==i) %>% dplyr::select(id) %>% mutate(module = i )
  
}
module_gene_lists <- lapply(module_list, get_module_gene_list, gene_module_df)
names(module_gene_lists) <- paste("module",module_list, sep= "_")
list2env(module_gene_lists, envir = .GlobalEnv)

# get compiled module gene list
module_gene_lists_bind <- purrr::map_df(module_gene_lists , data.frame, .id = 'name')

# get background genes
bkgr_genes <- row.names(graph_test_result_subset)
bkgr_genes_entrez <- AnnotationDbi::select(org.Hs.eg.db,keys=bkgr_genes,columns=c("GENENAME","ENTREZID"),keytype="SYMBOL") %>%
  filter(!is.na(ENTREZID))

# annotate module gene lists
anno_gene_lists <- function(x, data ){
  
  mod_genes_entrez <- AnnotationDbi::select(org.Hs.eg.db,keys=x$id,columns=c("GENENAME","ENTREZID"),keytype="SYMBOL")
  # Run goana enrichment analysis on all modules
  goana_results <- goana(mod_genes_entrez$ENTREZID,FDR=0.05,universe=data$ENTREZID,trend=T)
  goana_sig <- goana_results[goana_results$P.DE<=0.05,]
  goana_sig[order(goana_sig$P.DE,decreasing=FALSE),]
  
}

sort(paste("module",module_list, sep= "_"))

module_names <-  list(
  "module_1"=  module_1,
  "module_2"=  module_2,
  "module_3"=  module_3,
  "module_4"=  module_4 ,
  "module_5"=  module_5,
  "module_6"=  module_6,
  "module_7"=  module_7,
  "module_8"=  module_8,
  "module_9" = module_9 ,
  "module_10" =module_10,
  "module_11" =module_11,
  "module_12" =module_12,
  "module_13" =module_13,
  "module_14"=module_14)

mod_genes_entrez <- lapply(module_names, anno_gene_lists, bkgr_genes_entrez)
names(mod_genes_entrez) <- paste("module_GO_enrichment",module_list, sep= "_")
list2env(mod_genes_entrez, envir = .GlobalEnv)

# output individual datasets and also output compiled data
mod_genes_entrez_bind <- purrr::map_df(mod_genes_entrez , data.frame, .id = 'name')

### Remove redundant GO terms 

## Matt D gave me his code to remove redundant GO terms
# this function removes redundant terms from the goana annotation enrichment results
# in current form it simply removes less strongly enriched ancestor terms
# it does this by ordering by significance of enrichment, then removing any later terms
# that are ancestor terms in GO.db
# it expects a column with GO_ID

filter_goana_output <- function(goana_result, id_col = "GO_ID") {
  
  # Change column name to be GO_ID
  goana_result <- goana_result %>% rownames_to_column(var = "GO_ID") 
  
  require(GO.db)
  go_ancestors_all <-
    c(as.list(GOBPANCESTOR), as.list(GOCCANCESTOR), as.list(GOMFANCESTOR))
  goana_result_filtered <-
    goana_result %>%
    dplyr::arrange(P.DE)
  
  rows_to_drop_list <- list()
  for (row.tmp in 1:nrow(goana_result_filtered)) {
    rows_to_drop.tmp <-
      which(
        goana_result_filtered[[id_col]] %in%
          go_ancestors_all[[
            goana_result_filtered[[id_col]][row.tmp]]]) %>%
      setdiff(1:row.tmp)
    if (length(rows_to_drop.tmp) > 0)
      rows_to_drop_list[[
        goana_result_filtered[[id_col]][row.tmp]]] <-
        rows_to_drop.tmp
  }
  
  rows_to_drop <- unique(unlist(rows_to_drop_list))
  if (length(rows_to_drop) > 0)
    goana_result_filtered <- goana_result_filtered[-rows_to_drop,]
  
  # put back in original order
  goana_result_filtered <-
    goana_result_filtered[
      order(
        match(goana_result_filtered[[id_col]], goana_result[[id_col]])),]
  
  return(goana_result_filtered)
}

# Run function on each dataset 
module_names_enrichment <- list(
  "module_GO_enrichment_1"= module_GO_enrichment_1,
  "module_GO_enrichment_2"= module_GO_enrichment_2,
  "module_GO_enrichment_3"= module_GO_enrichment_3,
  "module_GO_enrichment_4"= module_GO_enrichment_4,
  "module_GO_enrichment_5" = module_GO_enrichment_5 ,
  "module_GO_enrichment_6"= module_GO_enrichment_6,
  "module_GO_enrichment_7"= module_GO_enrichment_7,
  "module_GO_enrichment_8"= module_GO_enrichment_8,
  "module_GO_enrichment_9"= module_GO_enrichment_9,
  "module_GO_enrichment_10" = module_GO_enrichment_10,
  "module_GO_enrichment_11" = module_GO_enrichment_11,
  "module_GO_enrichment_12" = module_GO_enrichment_12,
  "module_GO_enrichment_13" = module_GO_enrichment_13,
  "module_GO_enrichment_14" = module_GO_enrichment_14
)

# lapply function to reduce the redundant GO terms
mod_genes_entrez_bind_reduced_GO <- lapply(module_names_enrichment,filter_goana_output )
names(mod_genes_entrez_bind_reduced_GO) <- paste("module_GO_enrichment_reduced",module_list, sep= "_")
# combine the results
mod_genes_entrez_bind_reduced_GO_bind <- purrr::map_df(mod_genes_entrez_bind_reduced_GO , data.frame, .id = 'name')

### Investigate results

# remove CC terms since I am less interested in these, and only keep the more significant terms
mod_genes_entrez_bind_reduced_GO_bind_no_CC <- mod_genes_entrez_bind_reduced_GO_bind %>% filter(Ont!="CC") %>% filter(P.DE <= 0.0001)
mod_genes_entrez_bind_reduced_GO_bind_no_CC$P.DE <- as.numeric(mod_genes_entrez_bind_reduced_GO_bind_no_CC$P.DE)

## Plot results of GO as heatmap to compare and contrast - with also filtering to keep only the top 50 terms
mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide <- mod_genes_entrez_bind_reduced_GO_bind_no_CC %>% ungroup() %>% 
  group_by(name) %>% 
  top_n(20) %>%
  dplyr::select(-N, -Term, -Ont, -DE) %>%
  pivot_wider(names_from = GO_ID, values_from = P.DE) %>% column_to_rownames(var = "name")
head(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide )
mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide[is.na(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide)] <- 0
head(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide )

pdf(file.path(plotDir, "GO_term_heatmap.pdf"), height = 12, width = 7)
ComplexHeatmap::Heatmap(t(as.matrix(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide)))
dev.off()

# View also with term as the column name
mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term <- mod_genes_entrez_bind_reduced_GO_bind_no_CC %>% ungroup() %>% 
  group_by(name) %>% 
  top_n(20) %>%
  dplyr::select(-N, -GO_ID, -Ont, -DE) %>%
  pivot_wider(names_from = Term, values_from = P.DE) %>% column_to_rownames(var = "name")
head(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term )
mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term[is.na(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term)] <- 0
head(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term )

pdf(file.path(plotDir, "GO_term_heatmap_term.pdf"), height = 20, width = 7)
ComplexHeatmap::Heatmap(t(as.matrix(mod_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term)),
                        row_names_gp = gpar(fontsize = 6))
dev.off()

### Export results
for (i in names(module_gene_lists)) {
  module = module_gene_lists_bind %>% filter(name == i)
  write.csv(module, file = file.path(resultDir, paste0(i, ".csv")))
}


#### Plot Heatmap of NK cell mediated cytotoxicity gene expression in each cluster ####


# Load KEGG pathway hsa04650 using KEGGREST
# download the NK cell mediated cytotoxicity pathway
NK_mediated_cyto <- keggGet(c("hsa04650"))

# the first line is the KEGG gene entry number and the second line is the gene name
NK_mediated_cyto_genes_entry <- NK_mediated_cyto[[1]]$GENE
NK_mediated_cyto_gene_names <- as.data.frame(NK_mediated_cyto_genes_entry) %>% filter(grepl("KO:", NK_mediated_cyto_genes_entry)) %>% 
  # separate by semi colon to get the gene symbol
  separate(NK_mediated_cyto_genes_entry, into = c("gene_name", "gene_full_name"), sep = ";") %>%
  # make copy column for joining
  mutate(hsa04650_NK_gene = gene_name)

module_6_NK_genes_list <- NK_mediated_cyto_gene_names$hsa04650_NK_gene

# Get mean expression across the gene set
cds_NK_gene_list <- as.data.frame(t(exprs(cds_no_MAIT_no_9[rowData(cds_no_MAIT_no_9)$gene_short_name %in% module_6_NK_genes_list ,])) )
cds_NK_gene_list$barcode <- row.names(cds_NK_gene_list)

# join with cluster anno
all(cds_NK_gene_list$barcode == rownames(colData(cds_no_MAIT_no_9))) # TRUE

cds_NK_gene_list$Cluster <- colData(cds_no_MAIT_no_9)$Cluster.Name

# Calculate the mean of genes of interest for each cluster
cds_NK_gene_list_mean <-  cds_NK_gene_list %>%
  # not all genes in the NK pathway were identified
  filter_at(vars(module_6_NK_genes_list[module_6_NK_genes_list %in% colnames(cds_NK_gene_list)]), any_vars(. > 0)) %>%
  group_by(Cluster) %>%
  summarise_at(vars(module_6_NK_genes_list[module_6_NK_genes_list %in% colnames(cds_NK_gene_list)]), ~mean(. > 0)) %>% column_to_rownames(var = "Cluster")

# plot heatmap
pdf(file=file.path(plotDir, "cds_no_9_NK_mean_expression.pdf"), width = 15, height = 8)
ComplexHeatmap::Heatmap(as.matrix(cds_NK_gene_list_mean[,-c(1:5)]))# remove HLA columns because they are the same across all clusters
dev.off()


#### Assess overlap between gene module lists and ATAC-seq contrasts ####

# Load ATAC-seq contrast results for those that overlap the start site 

# CD57pos vs CD57minus up start
CD57pos_vs_CD57minus_anno_up_start <- read.csv("./P452_3_SAVED_DATA/ATACseqData_P452_3_norm_db_anno_up_start.csv") %>%
  mutate(name = "CD57pos_vs_CD57minus_anno_up_start" ) %>%
  dplyr::select(name, gene_name)
nrow(CD57pos_vs_CD57minus_anno_up_start) # 29

# CD57pos vs DN up start
CD57pos_DN_up_start <- read.csv("./P452_3_SAVED_DATA/ATACseqData_P452_3_norm_db_anno_CD57pos_DN_up_start.csv") %>% 
  mutate(name = "CD57pos_DN_up_start" ) %>%
  dplyr::select(name, gene_name)
nrow(CD57pos_DN_up_start ) # 562

# CD57minus vs DN up start
CD57minus_DN_up_start <- read.csv("./P452_3_SAVED_DATA/ATACseqData_P452_3_norm_db_anno_CD57minus_DN_up_start.csv") %>% 
  mutate(name = "CD57minus_DN_up_start" ) %>%
  dplyr::select(name, gene_name)
nrow(CD57minus_DN_up_start ) # 267 

# sites commonly up between CD57pos vs DN and CD57minus vs DN - but not filtered for overlapping start
CD57pos_DN_common_vs_CD57minus_DN_UP <- read.csv("./P452_3_SAVED_DATA/ATACseqData_P452_3_norm_db_anno_CD57pos_DN_common_vs_DN_UP.csv") %>%
  mutate(name = "CD57pos_DN_common_vs_CD57minus_DN_UP" ) %>%
  dplyr::select(name, gene_name)
nrow(CD57pos_DN_common_vs_CD57minus_DN_UP) # 740

# change module column names for binding
module_gene_lists_bind_gene <- module_gene_lists_bind %>% dplyr::rename(gene_name = id) %>% dplyr::select(-module)
module_gene_lists_bind_gene %>% dplyr::count(name)
#name    n
#1   module_1 1072
#2  module_10  290
#3  module_11  199
#4  module_12  138
#5  module_13   41
#6  module_14   27
#7   module_2 1031
#8   module_3 1019
#9   module_4  982
#10  module_5  844
#11  module_6  748
#12  module_7  616
#13  module_8  559
#14  module_9  476

# rbind the ATAC lists and the module gene lists
DE_module_ATAC_comparison <- rbind(CD57pos_vs_CD57minus_anno_up_start, 
                                   CD57pos_DN_up_start,
                                   CD57minus_DN_up_start, 
                                   CD57pos_DN_common_vs_CD57minus_DN_UP,
                                   module_gene_lists_bind_gene )
# spread list into binary matrix for upset plot, make set the columns
DE_module_ATAC_comparison_spread <- DE_module_ATAC_comparison %>% distinct() %>%
  mutate(Test= 1) %>%
  tidyr::pivot_wider(values_from = Test, names_from = name) %>% 
  replace(is.na(.), 0) %>% column_to_rownames(var = "gene_name")
colnames(DE_module_ATAC_comparison_spread)

# Plot upset plot using distinct mode (default)
pdf(file.path(plotDir, "DE_module_ATAC_comparison_spread_upset.pdf"), height = 10, width = 10)
UpSetR::upset(DE_module_ATAC_comparison_spread[,-c(2:4)], nsets = 20) # keep only 1
UpSetR::upset(DE_module_ATAC_comparison_spread[,-c(1,3:4)], nsets = 20) # keep 2
UpSetR::upset(DE_module_ATAC_comparison_spread[,-c(1:2,4)], nsets = 20) # keep 3
UpSetR::upset(DE_module_ATAC_comparison_spread[,-c(1:3)], nsets = 20) # keep 4
dev.off()

pdf(file.path(plotDir, "DE_module_ATAC_comparison_spread_upset_modules_4_6.pdf"))
UpSetR::upset(DE_module_ATAC_comparison_spread[,c(1,10:11)], nsets = 20) # keep only 1
UpSetR::upset(DE_module_ATAC_comparison_spread[,c(2, 10:11)], nsets = 20) # keep 2
UpSetR::upset(DE_module_ATAC_comparison_spread[,c(3, 10:11)], nsets = 20) # keep 3
UpSetR::upset(DE_module_ATAC_comparison_spread[,c(4,10:11)], nsets = 20) # keep 4
dev.off()

# only ook at overlaps with modules 4 and 6 
ATAC_contrasts <- c("CD57pos_vs_CD57minus_anno_up_start","CD57pos_DN_up_start", "CD57minus_DN_up_start", "CD57pos_DN_common_vs_CD57minus_DN_UP")

filter_ATAC_module_4_6_overlap <- function(x, data) {
  data <- data %>% dplyr::select(x,module_4, module_6 )
  data4 <- data %>% filter(data[,1] == 1 & module_4 == 1) %>% mutate(id = paste(x,"module_4_overlap"))
  data6 <- data %>% filter(data[,1] == 1 & module_6 == 1) %>% mutate(id = paste(x,"module_6_overlap"))
  rbind(data4, data6)
  
}

ATAC_contrasts_4_6 <- lapply(ATAC_contrasts,filter_ATAC_module_4_6_overlap ,  DE_module_ATAC_comparison_spread)
names(ATAC_contrasts_4_6) <- paste0("module_4_6_", ATAC_contrasts)
list2env(ATAC_contrasts_4_6, envir = .GlobalEnv)

rownames(module_4_6_CD57pos_vs_CD57minus_anno_up_start)
module_4_6_CD57pos_DN_up_start            
module_4_6_CD57minus_DN_up_start
module_4_6_CD57pos_DN_common_vs_CD57minus_DN_UP


# get lists of overlaps

# plot all as heatmap
pdf(file.path(plotDir, "DE_module_ATAC_comparison_spread_modules_4_6_heatmap.pdf"))
ComplexHeatmap::Heatmap(as.matrix(DE_module_ATAC_comparison_spread_4_6_subset_CD57pos_vs_CD57minus_anno_up_start)) # keep only 1
ComplexHeatmap::Heatmap(as.matrix(DE_module_ATAC_comparison_spread_4_6_subset_CD57pos_DN_up_start)) # keep 2
ComplexHeatmap::Heatmap(as.matrix(DE_module_ATAC_comparison_spread_4_6_subset_CD57minus_DN_up_start)) # keep 3
ComplexHeatmap::Heatmap(as.matrix(DE_module_ATAC_comparison_spread_4_6_subset_CD57pos_DN_common_vs_CD57minus_DN_UP)) # keep 4
dev.off()

#### Fine Tune Parameters of Monocle DGEA Analysis ####

set.seed(42) # the seed makes a difference for how genes are clustered into modules!

## Load results of graph test previously computed since this code took a while
load(file = "./EW_T1DAL_Results/P362-1_T1DAL_DGEA_rerun.Rdata")

# plot distribution of p values and moran's I values
class(graph_test_result_subset$q_value) # numeric
graph_test_result_subset %>% arrange(q_value) %>% 
  ggplot( aes(x = morans_I, y = q_value)) + geom_point()

q_moran <- graph_test_result_subset %>% arrange(q_value) %>% top_n(-1000) %>%
  ggplot( aes(x = morans_I, y = q_value)) + geom_point()

ggsave(q_moran, file = file.path(plotDir, "q_moran_relationship.pdf"))

# this shows there are a lot of significant points where the moran's I is very low

# plot the distribution of q-values
q_value_dist <- graph_test_result_subset %>% arrange(q_value) %>% 
  filter(morans_I >=0) %>%
  ggplot( aes(x = q_value)) + geom_histogram() +
  geom_vline(xintercept=0.05, linetype="dashed", 
             color = "red", size=1)
ggsave(q_value_dist, file = file.path(plotDir, "q_value_dist.pdf"))

graph_test_result_subset %>% arrange(q_value) %>% 
  filter(morans_I <0) %>%
  ggplot( aes(x = q_value)) + geom_histogram() +
  geom_vline(xintercept=0.05, linetype="dashed", 
             color = "red", size=1)

graph_test_result_subset %>% arrange(q_value) %>% 
  filter(q_value <=0.05) %>%
  ggplot( aes(x = q_value)) + geom_histogram()

graph_test_result_subset %>% arrange(q_value) %>% 
  filter(q_value <=0.05) %>%
  ggplot( aes(x = morans_I)) + geom_histogram()

# subset for significant genes at q-value of 0 and 0.05
subset_DEG_IDs_pos_moran_0 <- row.names(subset(graph_test_result_subset, q_value ==0 & morans_I > 0))
length(subset_DEG_IDs_pos_moran_0) # 543
subset_DEG_IDs_pos_moran_0.05 <- row.names(subset(graph_test_result_subset, q_value <=0.05 & morans_I > 0))
length(subset_DEG_IDs_pos_moran_0.05) # 8007

## Find modules of co-regulated genes
# find_gene_modules() essentially runs UMAP on the genes (as opposed to the cells) 
# and then groups them into modules using Louvain community analysis
# modules at this point are not related to previous clusters
# I tested multiple different resolution values

gene_module_df_0 <- find_gene_modules(subset_for_DGEA[subset_DEG_IDs_pos_moran_0,], resolution = c( 0.05))
gene_module_df_0_res_0.01 <- find_gene_modules(subset_for_DGEA[subset_DEG_IDs_pos_moran_0,], resolution = c( 0.01))

#save(gene_module_df_0, file = file.path(resultDir,"gene_module_df_0.Rdata"))
#save(gene_module_df_0_res_0.01, file = file.path(resultDir,"gene_module_df_0_res_0.01.Rdata"))

# re-run with q-value filtering at 0.05 rather than 0
gene_module_df_0.05 <- find_gene_modules(subset_for_DGEA[subset_DEG_IDs_pos_moran_0.05,], resolution = c( 0.05))
gene_module_df_0.05_res_0.01 <- find_gene_modules(subset_for_DGEA[subset_DEG_IDs_pos_moran_0.05,], resolution = c( 0.01))

#save(gene_module_df_0.05, file = file.path(resultDir,"gene_module_df_0.05.Rdata"))
#save(gene_module_df_0.05_res_0.01, file = file.path(resultDir,"gene_module_df_0.05_res_0.01.Rdata"))

load(file = file.path(resultDir,"gene_module_df_0.Rdata"))
load(file = file.path(resultDir,"gene_module_df_0_res_0.01.Rdata"))
load(file = file.path(resultDir,"gene_module_df_0.05.Rdata"))
load(file = file.path(resultDir,"gene_module_df_0.05_res_0.01.Rdata"))


### Plot module heatmap by cluster - to see which module genes are expressed in which clusters

agg_mat_0 <- aggregate_gene_expression(subset_for_DGEA, gene_module_df_0, cell_group_df)
row.names(agg_mat_0) <- stringr::str_c("Module ", row.names(agg_mat_0))

pdf(file =file.path(plotDir, "graph_auto_DEG_clusters_0.pdf"))
pheatmap::pheatmap(agg_mat_0,
                   scale="column", clustering_method="ward.D2")
dev.off()

# check with resolution of 0.01
agg_mat_0_res_0.01 <- aggregate_gene_expression(subset_for_DGEA, gene_module_df_0_res_0.01, cell_group_df)
# change row names to be letters instead to avoid confusion
row.names(agg_mat_0_res_0.01) <- stringr::str_c("Module ", c("A","B","C","D","E"))

pdf(file =file.path(plotDir, "graph_auto_DEG_clusters_0_res_0.01.pdf"), height = 4, width = 4.5)
pheatmap::pheatmap(agg_mat_0_res_0.01,
                   scale="column", clustering_method="ward.D2", fontsize = 16)
dev.off()

# check with altered q_value filtering
agg_mat_0.05 <- aggregate_gene_expression(subset_for_DGEA,gene_module_df_0.05, cell_group_df)
row.names(agg_mat_0.05 ) <- stringr::str_c("Module ", row.names(agg_mat_0.05 ))

pdf(file =file.path(plotDir, "graph_auto_DEG_clusters_0.05.pdf"))
pheatmap::pheatmap(agg_mat_0.05,
                   scale="column", clustering_method="ward.D2")
dev.off()

# check with resolution of 0.01
agg_mat_0.05_res_0.01 <- aggregate_gene_expression(subset_for_DGEA, gene_module_df_0.05_res_0.01, cell_group_df)
row.names(agg_mat_0.05_res_0.01) <- stringr::str_c("Module ", row.names(agg_mat_0.05_res_0.01))

pdf(file =file.path(plotDir, "graph_auto_DEG_clusters_0.05_res_0.01.pdf"))
pheatmap::pheatmap(agg_mat_0.05_res_0.01,
                   scale="column", clustering_method="ward.D2")
dev.off()

# plot expression of gene modules with plot cells on subset data
coregulated_gene_modules_0 <- plot_cells(subset_for_DGEA,
                                         genes=gene_module_df_0[gene_module_df_0$module %in% c(4,10,6,7),],
                                         label_cell_groups=TRUE,
                                         show_trajectory_graph=FALSE,
                                         cell_size=3)
ggsave(coregulated_gene_modules_0, file = file.path(plotDir, "coregulated_gene_module_expression_0.pdf"), height= 8, width = 10)


# plot expression of gene modules with plot cells on subset data - repeat with revised resolution
coregulated_gene_modules_0.01_flipped <- plot_cells(cds_no_MAIT, # plot on flipped UMAP
                                                    genes=gene_module_df_0_res_0.01[gene_module_df_0_res_0.01$module ==1,],
                                                    label_cell_groups=FALSE, # labels are plotting backwards
                                                    show_trajectory_graph=FALSE,
                                                    # group_label_size = 4,
                                                    cell_size=3) + theme(text = element_text(size = 20))
ggsave(coregulated_gene_modules_0.01_flipped, file = file.path(plotDir, "coregulated_gene_modules_0.01_flipped.pdf"), height= 4, width =6)



## Get module gene lists and run GO analysis 
# Example: https://ucdavis-bioinformatics-training.github.io/2019_August_UCD_mRNAseq_Workshop/differential_expression/enrichment

#### Annotate DGEA modules and get GO terms ####

## Get gene lists for each module
module_list_0 <- unique(gene_module_df_0$module)

module_gene_lists_0 <- lapply(module_list_0, get_module_gene_list, gene_module_df_0)
names(module_gene_lists_0) <- paste("module_0_list_",module_list_0, sep= "_")
list2env(module_gene_lists_0, envir = .GlobalEnv)

# get compiled module gene list
module_gene_lists_bind_0 <- purrr::map_df(module_gene_lists_0 , data.frame, .id = 'name')
View(module_gene_lists_bind_0)

## Get gene list for modules with resolution 0.01
module_list_0_res_0.01 <- unique(gene_module_df_0_res_0.01$module)
module_gene_lists_0_res_0.01 <- lapply(module_list_0_res_0.01, get_module_gene_list, gene_module_df_0_res_0.01)
names(module_gene_lists_0_res_0.01) <- paste("module_0_list_0.01",module_list_0_res_0.01, sep= "_")
list2env(module_gene_lists_0_res_0.01, envir = .GlobalEnv)

module_gene_lists_0_res_0.01_bind <- purrr::map_df(module_gene_lists_0_res_0.01 , data.frame, .id = 'name')
module_gene_lists_0_res_0.01_bind %>% dplyr::filter(name == "module_0_list_0.01_1") %>% View()
write.csv(module_gene_lists_0_res_0.01_bind  , file = file.path(resultDir, "module_gene_lists_0_res_0.01.csv"))

module_gene_lists_0_res_0.01_bind_1 <- module_gene_lists_0_res_0.01_bind %>% dplyr::filter(name == "module_0_list_0.01_1") 

# check for exhaustion genes 
figure_ex_gene_list <- c("TOX", "TBX21", "TIGIT", "KLRG1", "EOMES", "PDCD1", "IL7R", "CTLA4", "HAVCR2", "LAG3",
                         "SLAMF6","ZEB2","CX3CR1", "GZMA","GZMB", "CD101","NKG7", "BATF",
                         "KLRD1", "CD160", "LILRB1", "KIR2DL1", "KIR2DS4", "KIR2DL3", "KIR3DL1", "KIR3DL2", "FCGR3A")

module_gene_lists_0_res_0.01_bind_1[module_gene_lists_0_res_0.01_bind_1$id %in% figure_ex_gene_list,]

# get GO enrichment terms for each module in original set
module_0_names <-  list(
  "module_0_list__1"=  module_0_list__1,
  "module_0_list__2"=  module_0_list__2,
  "module_0_list__3"=  module_0_list__3,
  "module_0_list__4"=  module_0_list__4 ,
  "module_0_list__5"=  module_0_list__5,
  "module_0_list__6"=  module_0_list__6,
  "module_0_list__7"=  module_0_list__7,
  "module_0_list__8"=  module_0_list__8,
  "module_0_list__9" = module_0_list__9 ,
  "module_0_list__10" =module_0_list__10)

mod_0_genes_entrez <- lapply(module_0_names, anno_gene_lists, bkgr_genes_entrez)
names(mod_0_genes_entrez) <- paste("module_0_GO_enrichment",module_list_0, sep= "_")
list2env(mod_0_genes_entrez, envir = .GlobalEnv)

# output individual datasets and also output compiled data
mod_0_genes_entrez_bind <- purrr::map_df(mod_0_genes_entrez , data.frame, .id = 'name')

# get GO enrichment terms for each module in original set with res = 0.01
names(module_gene_lists_0_res_0.01) 
module_0_res_0.01_names <-  list(
  "module_0_list_0.01_1" = module_0_list_0.01_1, 
  "module_0_list_0.01_4" = module_0_list_0.01_4,
  "module_0_list_0.01_3" = module_0_list_0.01_3,
  "module_0_list_0.01_5" = module_0_list_0.01_5,
  "module_0_list_0.01_2" = module_0_list_0.01_2)

mod_0_res_0.01_genes_entrez <- lapply(module_0_res_0.01_names, anno_gene_lists, bkgr_genes_entrez)
names(mod_0_res_0.01_genes_entrez) <- paste("module_0_res_0.01_GO_enrichment",module_list_0_res_0.01, sep= "_")
list2env(mod_0_res_0.01_genes_entrez, envir = .GlobalEnv)

# output individual datasets and also output compiled data
mod_0_res_0.01_genes_entrez_bind <- purrr::map_df(mod_0_res_0.01_genes_entrez , data.frame, .id = 'name')

### Remove redundant GO terms 

# Run function on each dataset 
module_0_names_enrichment <- list(
  "module_0_GO_enrichment_1"= module_0_GO_enrichment_1,
  "module_0_GO_enrichment_2"= module_0_GO_enrichment_2,
  "module_0_GO_enrichment_3"= module_0_GO_enrichment_3,
  "module_0_GO_enrichment_4"= module_0_GO_enrichment_4,
  "module_0_GO_enrichment_5"= module_0_GO_enrichment_5 ,
  "module_0_GO_enrichment_6"= module_0_GO_enrichment_6,
  "module_0_GO_enrichment_7"= module_0_GO_enrichment_7,
  "module_0_GO_enrichment_8"= module_0_GO_enrichment_8,
  "module_0_GO_enrichment_9"= module_0_GO_enrichment_9,
  "module_0_GO_enrichment_10"=module_0_GO_enrichment_10
)

# lapply function to reduce the redundant GO terms
mod_0_genes_entrez_bind_reduced_GO <- lapply(module_0_names_enrichment,filter_goana_output )
names(mod_0_genes_entrez_bind_reduced_GO) <- paste("module_0_GO_enrichment_reduced",module_list_0, sep= "_")
# combine the results
mod_0_genes_entrez_bind_reduced_GO_bind <- purrr::map_df(mod_0_genes_entrez_bind_reduced_GO , data.frame, .id = 'name')

### Investigate results

# remove CC terms since I am less interested in these, and only keep the more significant terms
mod_0_genes_entrez_bind_reduced_GO_bind_no_CC <- mod_0_genes_entrez_bind_reduced_GO_bind %>% filter(Ont!="CC") %>% filter(P.DE <= 0.0001)
mod_0_genes_entrez_bind_reduced_GO_bind_no_CC$P.DE <- as.numeric(mod_0_genes_entrez_bind_reduced_GO_bind_no_CC$P.DE)

## Plot results of GO as heatmap to compare and contrast - with also filtering to keep only the top 50 terms
mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide <- mod_0_genes_entrez_bind_reduced_GO_bind_no_CC %>% ungroup() %>% 
  group_by(name) %>% 
  top_n(20) %>%
  dplyr::select(-N, -Term, -Ont, -DE) %>%
  pivot_wider(names_from = GO_ID, values_from = P.DE) %>% column_to_rownames(var = "name")
head(mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide )
mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide[is.na(mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide)] <- 0
head(mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide )

pdf(file.path(plotDir, "module_0_GO_term_heatmap.pdf"), height = 12, width = 7)
ComplexHeatmap::Heatmap(t(as.matrix(mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide)))
dev.off()

# View also with term as the column name
mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term <- mod_0_genes_entrez_bind_reduced_GO_bind_no_CC %>% ungroup() %>% 
  group_by(name) %>% 
  top_n(20) %>%
  dplyr::select(-N, -GO_ID, -Ont, -DE) %>%
  pivot_wider(names_from = Term, values_from = P.DE) %>% column_to_rownames(var = "name")
mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term[is.na(mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term)] <- 0

pdf(file.path(plotDir, "module_0_GO_term_heatmap_term.pdf"), height = 20, width = 7)
ComplexHeatmap::Heatmap(t(as.matrix(mod_0_genes_entrez_bind_reduced_GO_bind_no_CC_wide_term)),
                        row_names_gp = gpar(fontsize = 6))
dev.off()

# export results from tuned moran's I
for (i in names(module_gene_lists_0)) {
  module = module_gene_lists_bind_0 %>% filter(name == i)
  write.csv(module, file = file.path(resultDir, paste0(i, ".csv")))
}

#### Make bubble plot of tuned DGEA KEGG enrichment from string ####


## Load KEGG enrichment results from stringdb for module 1 

module_gene_lists_0_res_0.01_module1_enrichment.KEGG <- read.table(file = file.path(resultDir, "module_gene_lists_0_res_0.01_module1_enrichment.KEGG.tsv"), sep = "\t", header = TRUE)
module_gene_lists_0_res_0.01_module1_enrichment.KEGG$group <- "Module 1"

# Plot bubble plot
module_gene_lists_0_res_0.01_module1_enrichment.KEGG_bubble <- module_gene_lists_0_res_0.01_module1_enrichment.KEGG %>% 
  top_n(-10, false.discovery.rate) %>%
  ggplot( aes(x = group, y = term.description, size =observed.gene.count, color = false.discovery.rate)) +
  geom_point() + 
  theme( text = element_text(size = 12), 
         legend.text = element_text(size = 4)) + labs(x = NULL, y = "KEGG Term") + 
  scale_color_viridis_c(direction = -1) +
  theme_minimal() 

ggsave(module_gene_lists_0_res_0.01_module1_enrichment.KEGG_bubble, file = file.path(plotDir, "module_gene_lists_0_res_0.01_module1_enrichment.KEGG_bubble.pdf"), device = "pdf", width = 5, height = 4)



#### Test overlap of tuned DGEA analysis with ATAC seq data ####

# Load ATAC-seq contrast results for those all genes and not just those that are overlapping the start site

ATAC_file_list = c("ATACseqData_P452_3_norm_db_anno_CD57pos_DN_up.csv",
                   "ATACseqData_P452_3_norm_db_anno_CD57minus_DN_up.csv",
                   "ATACseqData_P452_3_norm_db_anno_up.csv",
                   "ATACseqData_P452_3_norm_db_anno_down.csv")

read_atac <- function(i) {
  
  df <- read.csv(paste0("./P452_3_SAVED_DATA/",i)) %>%
    mutate(atac_list = i) %>% mutate(atac_list = str_remove(atac_list, "ATACseqData_P452_3_norm_db_anno_")) %>%
    mutate(atac_list = str_remove(atac_list, ".csv"))
  df
}
atac_all_up <- lapply(ATAC_file_list, read_atac)  
names(atac_all_up) <- ATAC_file_list
atac_all_up_df <- bind_rows(atac_all_up)

## Find overlaps with between ATAC lists and each module gene list

atac_module_match <- function (i, x) {
  module =x %>% filter(name == i) %>% dplyr::rename(gene_name = id)
  module_atac <- left_join(atac_all_up_df, module) %>% filter(!is.na(module)) %>%
    group_by(atac_list) %>% mutate(atac_list_count = n()) 
  module_atac
}

atac_module_match_df <- lapply(names(module_gene_lists_0),atac_module_match, module_gene_lists_bind_0 )
names(atac_module_match_df)  <- paste0(names(module_gene_lists_0), "atac_match")
atac_module_match_df <- bind_rows(atac_module_match_df)     

## calculate enrichment values for each overlapping list
# get module and atac overlap for each module, total entries in atac list, and total entries in module list
atac_module_match_df_counts <- atac_module_match_df %>% distinct(module, atac_list, atac_list_count) %>%
  left_join((atac_all_up_df %>% group_by(atac_list) %>% summarize(atac_list_total = n()))) %>%
  left_join((module_gene_lists_bind_0 %>% group_by(module) %>% summarize(module_list_total = n())))
summary(atac_module_match_df_counts)
## Calculate total background universe as total annotate RNAseq genes that overlap with annotated ATAC genes
# load full annotated background peaks from ATAC workspace
load( file ="./P452_3_SAVED_DATA/ATAC_background_peak_anno.Rdata")
atac_enrichment_overlap_bg <- nrow(background_peak_anno[background_peak_anno$gene_name %in% rownames(subset_for_DGEA), ]) 
class(atac_enrichment_overlap_bg) # integer

## Calculate enrichment for each overlap in each module

## Run test for enrichment with phyper for each set

# phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

# overlap = number atac list genes in module gene list (atac_list_count)
# group 2 = length of full atac gene list (atac_list_total) 
# group 1 = length total module gene signature (module_list_total)

# perform enrichment in loop 
ATAC_enrich <- NULL;
for(row in 1:nrow(atac_module_match_df_counts)) {
  x <- atac_module_match_df_counts[row, "atac_list"]
  y = atac_module_match_df_counts[row, "module"]
  overlap <- atac_module_match_df_counts[row, "atac_list_count"]
  group2  <- atac_module_match_df_counts[row, "atac_list_total"]
  group1 <- atac_module_match_df_counts[row, "module_list_total"]
  
  tmp <- data.frame(name = x, module= y, 
                    enrich =phyper(as.numeric(overlap-1), as.numeric(group2), as.numeric(atac_enrichment_overlap_bg-group2), as.numeric(group1), lower.tail = FALSE))
  ATAC_enrich <- rbind(ATAC_enrich, tmp)
  
}

atac_module_match_df_counts_enrich <- left_join( atac_module_match_df_counts, ATAC_enrich)

# find enriched lists
atac_module_match_df_counts_enrich %>% filter(enrich <= 0.05)
#atac_list       module atac_list_count atac_list_total module_list_total   enrich
#<chr>           <fct>            <int>           <int>             <int>    <dbl>
#1 CD57pos_DN_up   10                   9            1484                21 1.34e- 4
#2 CD57minus_DN_up 10                   7             849                21 1.49e- 4
#3 CD57pos_DN_up   4                   35            1484                64 2.29e-18
#4 CD57minus_DN_up 4                   31             849                64 1.97e-21
#5 down            3                    6             195                80 8.10e- 4

# investigate the overlaps
atac_module_match_df %>% filter(module==4 & atac_list == "CD57pos_DN_up") %>% distinct(gene_name) %>% View()
atac_module_match_df %>% filter(module==4 & atac_list == "CD57minus_DN_up") %>% distinct(gene_name) %>% View()
atac_module_match_df %>% filter(module==4) %>%  group_by(gene_name) %>% mutate(counts = n()) %>% filter(counts >1 ) %>% distinct(gene_name) %>%
  View()
atac_module_match_df %>% filter(module==10) %>% View() 
atac_module_match_df %>% filter(module==10) %>% distinct(gene_name)  
atac_module_match_df %>% filter(module==7) %>% distinct(gene_name)  

### Repeat this analysis for the new tuned DGEA set with resolution of 0.01

## Find overlaps with between ATAC lists and each module gene list

atac_module_res_0.01_match_df <- lapply( names(module_gene_lists_0_res_0.01) ,atac_module_match, module_gene_lists_0_res_0.01_bind )
names(atac_module_res_0.01_match_df)  <- paste0( names(module_gene_lists_0_res_0.01), "atac_match")
atac_module_res_0.01_match_df <- bind_rows(atac_module_res_0.01_match_df)     

## calculate enrichment values for each overlapping list
# get module and atac overlap for each module, total entries in atac list, and total entries in module list
atac_module_res_0.01_match_df_counts <- atac_module_res_0.01_match_df %>% distinct(module, atac_list, atac_list_count) %>%
  left_join((atac_all_up_df %>% group_by(atac_list) %>% summarize(atac_list_total = n()))) %>%
  left_join((module_gene_lists_0_res_0.01_bind %>% group_by(module) %>% summarize(module_list_total = n())))
summary(atac_module_res_0.01_match_df_counts)

## Run test for enrichment with phyper for each set

# phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

# overlap = number atac list genes in module gene list (atac_list_count)
# group 2 = length of full atac gene list (atac_list_total) 
# group 1 = length total module gene signature (module_list_total)

# perform enrichment in loop 
ATAC_enrich_res_0.01 <- NULL;
for(row in 1:nrow(atac_module_res_0.01_match_df_counts)) {
  x <- atac_module_res_0.01_match_df_counts[row, "atac_list"]
  y = atac_module_res_0.01_match_df_counts[row, "module"]
  overlap <- atac_module_res_0.01_match_df_counts[row, "atac_list_count"]
  group2  <- atac_module_res_0.01_match_df_counts[row, "atac_list_total"]
  group1 <- atac_module_res_0.01_match_df_counts[row, "module_list_total"]
  
  tmp <- data.frame(name = x, module= y, 
                    enrich =phyper(as.numeric(overlap-1), as.numeric(group2), as.numeric(atac_enrichment_overlap_bg-group2), as.numeric(group1), lower.tail = FALSE))
  ATAC_enrich_res_0.01 <- rbind(ATAC_enrich_res_0.01, tmp)
  
}

atac_module_res_0.01_match_df_counts_enrich <- left_join( atac_module_res_0.01_match_df_counts, ATAC_enrich_res_0.01)

# find enriched lists
atac_module_res_0.01_match_df_counts_enrich %>% filter(enrich <= 0.05)
#atac_list       module atac_list_count atac_list_total module_list_total   enrich
#<chr>           <fct>            <int>           <int>             <int>    <dbl>
#1 CD57pos_DN_up   1                   57            1484               171 3.97e-16
#2 CD57minus_DN_up 1                   41             849               171 9.54e-15
#3 down            2                    5             195               135 3.87e- 2

atac_module_res_0.01_match_df  %>% filter(module==1 & atac_list == "CD57pos_DN_up") %>% distinct(gene_name) %>% View()
atac_module_res_0.01_match_df_ATAC_module1_CD57pos_DN_up <- atac_module_res_0.01_match_df  %>% filter(module==1 & atac_list == "CD57pos_DN_up") %>% distinct(gene_name) 
write.csv(atac_module_res_0.01_match_df_ATAC_module1_CD57pos_DN_up, file = file.path(resultDir, "atac_module_res_0.01_match_df_ATAC_module1_CD57pos_DN_up.csv"))

atac_module_res_0.01_match_df  %>% filter(module==1 & atac_list == "CD57minus_DN_up") %>% distinct(gene_name) %>% View()

#### Calculate DGEA overlap with KD gene sets ####

# load kirsten's PD1 and CD57 gene sets
CD57pos_vs_DN_KD_gene_set <- read_csv("./Kirsten_Data/P362-1 T1DAL 10X/input gene sets/CD57pos_vs_DN.csv",
                                      skip=1, col_names = "gene_name")
PD1pos_vs_DN_KD_gene_set <- read_csv("./Kirsten_Data/P362-1 T1DAL 10X/input gene sets/PD1pos_vs_DN.csv", 
                                     skip=1, col_names = "gene_name")

# How many genes overlap with module 1 list?
module_gene_lists_0_res_0.01_bind_module_1 <- module_gene_lists_0_res_0.01_bind %>% filter(name == "module_0_list_0.01_1")

module_gene_lists_0_res_0.01_bind_module_1_CD57_overlap <- module_gene_lists_0_res_0.01_bind_module_1[module_gene_lists_0_res_0.01_bind_module_1$id %in% CD57pos_vs_DN_KD_gene_set$gene_name,]
nrow(module_gene_lists_0_res_0.01_bind_module_1_CD57_overlap) # 19
module_gene_lists_0_res_0.01_bind_module_1_PD1_overlap <- module_gene_lists_0_res_0.01_bind_module_1[module_gene_lists_0_res_0.01_bind_module_1$id %in% PD1pos_vs_DN_KD_gene_set$gene_name,]
nrow(module_gene_lists_0_res_0.01_bind_module_1_PD1_overlap) # 14 

# test enrichment for CD57 gene list overlap
phyper(nrow(module_gene_lists_0_res_0.01_bind_module_1_CD57_overlap)-1, nrow(CD57pos_vs_DN_KD_gene_set),
       nrow(subset_for_DGEA)-nrow(CD57pos_vs_DN_KD_gene_set),  nrow(module_gene_lists_0_res_0.01_bind_module_1),lower.tail= FALSE)
# 7.720915e-20
# test enrichment for PD1 gene list overlap 
phyper(nrow(module_gene_lists_0_res_0.01_bind_module_1_PD1_overlap)-1, nrow(PD1pos_vs_DN_KD_gene_set),
       nrow(subset_for_DGEA)-nrow(PD1pos_vs_DN_KD_gene_set),  nrow(module_gene_lists_0_res_0.01_bind_module_1),lower.tail= FALSE)
#9.982383e-20

#phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

CD57_PD1_overlap
CD57_PD1_overlap[CD57_PD1_overlap$gene_name %in% module_gene_lists_0_res_0.01_bind_module_1_CD57_overlap$id,]$gene_name # 12
CD57_PD1_overlap[CD57_PD1_overlap$gene_name %in% module_gene_lists_0_res_0.01_bind_module_1_PD1_overlap$id,] # 12

# overlap = number KD gene set overlapping with module 1 set
# group 2 = length of full KD set
# group 1 = length total module 1 set 
# total = total genes used to determine significant modules 


#### Upset plot and Venn diagram of fine tuned DGEA showing list overlaps with ATAC lists, PD1 and CD57 gene sets ####

# combine all lists
CD57pos_vs_DN_KD_gene_set_upset <- CD57pos_vs_DN_KD_gene_set %>% mutate(name = "CD57pos_vs_DN_KD")
PD1pos_vs_DN_KD_gene_set_upset <- PD1pos_vs_DN_KD_gene_set %>% mutate(name = "PD1pos_vs_DN_KD")
atac_all_up_df_upset <- atac_all_up_df %>% dplyr::select(gene_name, atac_list) %>% 
  dplyr::rename(name = atac_list) %>% 
  filter(!(name %in% c("up","down") ))
atac_in_all <- atac_all_up_df_upset %>% group_by(gene_name) %>% filter(n() >=2) %>%  distinct(gene_name) %>% mutate(name = "Shared Up")
module_gene_lists_0_res_0.01_bind_module_1_upset <- module_gene_lists_0_res_0.01_bind_module_1 %>% dplyr::rename(gene_name = id) %>% 
  dplyr::select(gene_name, name)

upset_atac_CD57_PD1 <- rbind(#CD57pos_vs_DN_KD_gene_set_upset, 
  #PD1pos_vs_DN_KD_gene_set_upset,
  #atac_all_up_df_upset,
  module_gene_lists_0_res_0.01_bind_module_1_upset,
  atac_in_all) %>% 
  # Add counter column for spreading data
  mutate(Test = 1)

# make A binary matrix/data frame where rows are genes and columns are sets and data is filled in as 1 and zero
upset_atac_CD57_PD1_binary <- upset_atac_CD57_PD1 %>%
  distinct(Test, gene_name, name) %>%
  tidyr::pivot_wider(values_from = Test, names_from = name) %>% 
  replace(is.na(.), 0) %>% column_to_rownames(.,var = "gene_name")

# Plot upset plot using distinct mode (default)
pdf(file.path(plotDir, "upset_atac_CD57_PD1.pdf"), height = 8, width = 10)
UpSetR::upset(upset_atac_CD57_PD1_binary, order.by = "freq")
dev.off()

# try to plot in intersect mode but still ran into same inexplicable error:
# "Error: code for combination set should only contain 0 and 1."
m2 = make_comb_mat(lt, mode = "intersect")
lt = list(a = unique(CD57pos_vs_DN_KD_gene_set_upset$gene_name),
          b = unique(PD1pos_vs_DN_KD_gene_set_upset$gene_name),
          c = unique(module_gene_lists_0_res_0.01_bind_module_1_upset$gene_name))
set = make_comb_mat(lt,mode = "intersect")
UpSetR::upset(set)

## make a venn diagram instead of the two comparisons of interest with scaled venn diagrams

# make list of overlaps
module_1_KD_overlap <- 
  list(CD57pos_vs_DN_KD = unique(CD57pos_vs_DN_KD_gene_set_upset$gene_name),
       PD1pos_vs_DN_KD= unique(PD1pos_vs_DN_KD_gene_set_upset$gene_name),
       module_1 = unique(module_gene_lists_0_res_0.01_bind_module_1_upset$gene_name))
module_1_KD_overlap_venn <- gplots::venn(module_1_KD_overlap, show.plot = FALSE)
lengths(attributes(module_1_KD_overlap_venn)$intersections)

# Input in the form of a named numeric vector
module_1_KD_overlap_fit <- eulerr::euler(c(A = 83, 
                                           B = 9,
                                           C = 150,
                                           "A&B" = 21,
                                           "A&C"=7,
                                           "B&C"= 2,
                                           "A&B&C"= 12))

# plot the resulting Venn diagram
pdf(file = file.path(plotDir, "module_1_KD_overlap_Venn_eulerr.pdf"), height = 4, width = 6)
plot(module_1_KD_overlap_fit,
     quantities = TRUE,
     fill = "transparent",
     lty = 1:3,
     labels = c("CD57", "PD1", "Module 1"))
dev.off()

# plot overlap of ust PD1 and CD57
module_1_KD_overlap_only <- 
  list(CD57pos_vs_DN_KD = unique(CD57pos_vs_DN_KD_gene_set_upset$gene_name),
       PD1pos_vs_DN_KD= unique(PD1pos_vs_DN_KD_gene_set_upset$gene_name))
module_1_KD_overlap_venn_only <- gplots::venn(module_1_KD_overlap_only, show.plot = FALSE)
lengths(attributes(module_1_KD_overlap_venn_only)$intersections)

# Input in the form of a named numeric vector
module_1_KD_overlap_only_fit <- eulerr::euler(c(A = 90, 
                                                B = 11,
                                                "A&B" = 33
))

pdf(file = file.path(plotDir, "module_1_KD_overlap_only_Venn_eulerr.pdf"), height = 4, width = 6)
plot(module_1_KD_overlap_only_fit,
     quantities = TRUE,
     fill = "transparent",
     lty = 1:2,
     labels = c("CD57-like", "PD1-like"))
dev.off()

# export overlapping PD1 and CD57 gene list genes
CD57_PD1_overlap <- CD57pos_vs_DN_KD_gene_set_upset[CD57pos_vs_DN_KD_gene_set_upset$gene_name %in% PD1pos_vs_DN_KD_gene_set_upset$gene_name,] 
write.csv(CD57_PD1_overlap$gene_name, file= file.path(resultDir, "CD57_PD1_overlap.csv"), row.names = FALSE)

# repeat with module 1 vs ATAC
module_1_atac_overlap <- 
  list(shared_up = unique(atac_in_all$gene_name),
       module_1 = unique(module_gene_lists_0_res_0.01_bind_module_1_upset$gene_name))
module_1_atac_overlap_venn <- gplots::venn(module_1_atac_overlap, show.plot = FALSE)
lengths(attributes(module_1_atac_overlap_venn)$intersections)

# confirm enrichment 
phyper(36-1, length(unique(atac_in_all$gene_name)),
       nrow(subset_for_DGEA)-length(unique(atac_in_all$gene_name)),  nrow(module_gene_lists_0_res_0.01_bind_module_1),lower.tail= FALSE)
# 4.996535e-19

#phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

# overlap = number of ATAC set in module 1 set 
# group 2 = length of full atac set used to calculate the overlap
# group 1 = length total module 1 set 
# total = total genes used to determine significant modules 

# Input in the form of a named numeric vector
module_1_atac_overlap_fit <- eulerr::euler(c(A = 722, 
                                             B = 135,
                                             "A&B" = 36))

# plot the resulting Venn diagram
pdf(file = file.path(plotDir, "module_1_atac_overlap_Venn_eulerr.pdf"), height = 3, width = 3)
plot(module_1_atac_overlap_fit,
     quantities = TRUE,
     fill = "transparent",
     lty = 1:3,
     labels = c("Shared Up ATAC", "Module 1"),
     main = list(label = "Enrichment p-value = 5e-19", cex = 1))
dev.off()

# which genes are in this overlap?
module1_shared_up_overlap <- data.frame(gene_name = attributes(module_1_atac_overlap_venn)$intersections$`shared_up:module_1`)
write.csv(module1_shared_up_overlap, file= file.path(resultDir, "module1_shared_up_overlap.csv"))

# which genes unique to shared up
module1_shared_up_unique_atac <- data.frame(gene_name = attributes(module_1_atac_overlap_venn)$intersections$`shared_up`)
write.csv(module1_shared_up_unique_atac, file= file.path(resultDir, "module1_shared_up_unique_atac.csv"))

# which genes unique to module 1
module1_shared_up_unique_module1 <- data.frame(gene_name = attributes(module_1_atac_overlap_venn)$intersections$`module_1`)
write.csv(module1_shared_up_unique_module1, file= file.path(resultDir, "module1_shared_up_unique_module1.csv"))



#### Heatmap of zscore gene expression of module1_shared_up_overlap ####

module1_shared_up_overlap_list <- module1_shared_up_overlap$gene_name

module1_shared_up_overlap_exp <- as.data.frame(t(exprs(subset_for_DGEA[rowData(subset_for_DGEA)$gene_short_name %in% module1_shared_up_overlap_list,])) )
module1_shared_up_overlap_exp$Cell.ID <- row.names(module1_shared_up_overlap_exp)
colData(subset_for_DGEA)$Cluster <- colData(subset_for_DGEA)$Cluster.Name

module1_shared_up_overlap_anno_data <- as.data.frame(colData(subset_for_DGEA)) 
module1_shared_up_overlap_anno_data$Cell.ID <- row.names(colData(subset_for_DGEA))

# Calculate the mean of genes of interest, 
module1_shared_up_overlap_exp_anno <-  module1_shared_up_overlap_exp %>%
  merge(module1_shared_up_overlap_anno_data, by = "Cell.ID") %>%
  filter_at(vars(module1_shared_up_overlap_list), any_vars(. > 0)) %>%
  group_by(Cluster.Name, Timepoint) %>%
  summarise_at(vars(module1_shared_up_overlap_list), ~mean(. > 0)) 

cluster_means <- module1_shared_up_overlap_exp_anno

## Use scale to subtract the mean values and get zscore
# ?scale(): If center is a numeric-alike vector with length equal to the number of columns of x,
#then each column of x has the corresponding value from center subtracted from it.
col_scaled_mean <- scale(cluster_means[,module1_shared_up_overlap_list])
row.names(col_scaled_mean) <- cluster_means$Cluster.Name

cols = c("#FFFF99","#1F78B4","#FB9A99","#FDBF6F","#B2DF8A","#FF7F00","gray","#EF8A62","#CAB2D6","#E31A1C","#E9A3C9","black","#6A3D9A","#A6CEE3","#B15928", "#33A02C")
cluster_cols <- cols
# names(cluster_cols) <- keep_clusts$Cluster 

# heatmap_anno <- rowAnnotation(Cluster = row.names(col_scaled_median_RP), col=list(Cluster=c(cluster_cols)))

pdf(file = file.path(plotDir, "module1_shared_up_overlap_zscore_heatmap.pdf"), height = 7, width =6)
Heatmap(t(col_scaled_mean),row_names_side="right",name="Z-score\nExpression",show_parent_dend_line=F, cluster_rows = T,
        cluster_columns = T,rect_gp = gpar(col = "white", lwd = 1), ##,right_annotation = heatmap_anno)
        column_names_gp = grid::gpar(fontsize = 16), row_names_gp = grid::gpar(fontsize = 16), 
        heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12), title_gp = grid::gpar(fontsize = 16)))

dev.off()

#### Plot all Module 1 gene expression ####
module_1_gene_list <- unique(module_gene_lists_0_res_0.01_bind_module_1_upset$gene_name)

module1_all_exp <- as.data.frame(t(exprs(subset_for_DGEA[rowData(subset_for_DGEA)$gene_short_name %in% module_1_gene_list,])) )
module1_all_exp$Cell.ID <- row.names(module1_all_exp)

module1_all_exp_data <- as.data.frame(colData(subset_for_DGEA)) 
module1_all_exp_data$Cell.ID <- row.names(colData(subset_for_DGEA))

# Calculate the mean of genes of interest across all cells
module1_all_exp_anno <-  module1_all_exp %>%
  merge(module1_all_exp_data, by = "Cell.ID") %>%
  filter_at(vars(module_1_gene_list), any_vars(. > 0)) %>%
  group_by(Cluster.Name, Timepoint) %>%
  summarise_at(vars(module_1_gene_list), ~mean(. > 0)) 

cluster_means_all_module_1 <- module1_all_exp_anno

## Use scale to subtract the mean values and get zscore
# ?scale(): If center is a numeric-alike vector with length equal to the number of columns of x,
#then each column of x has the corresponding value from center subtracted from it.
col_scaled_mean_module_1 <- scale(cluster_means_all_module_1[,module_1_gene_list])
row.names(col_scaled_mean_module_1) <- cluster_means_all_module_1$Cluster.Name

pdf(file = file.path(plotDir, "module1_all_exp_zscore_heatmap.pdf"), height = 30, width =5)
Heatmap(t(col_scaled_mean_module_1),row_names_side="right",name="Col Zscore",show_parent_dend_line=F, cluster_rows = T, cluster_columns = T,rect_gp = gpar(col = "white", lwd = 1)) ##,right_annotation = heatmap_anno)
dev.off()


#### Find genes that change as a function of pseudotime in monocle ####

# again use graph test but pass neighbor_graph="principal_graph"
#cds_no_MAIT_no_90_DE_PT_res <- graph_test(cds_no_MAIT_no_9, neighbor_graph="principal_graph",cores=4)
# kept getting a weird error when running this:
#Error: 'rBind' is defunct.
#Since R version 3.2.0, base's rbind() should work fine with S4 objects# found an online fix to manually edit the function and change rBind to rbing

#Here I find a good solution from https://groups.google.com/g/monocle-3-users/c/tBjYuAxwyEo
#trace('calculateLW', edit = T, where = asNamespace("monocle3"))
#change Matrix::rBind to rbind

# So I manually changed rBind to rbind and saved the function
cds_no_MAIT_no_90_DE_PT_res <- graph_test(cds_no_MAIT_no_9, neighbor_graph="principal_graph",cores=8)

# this running of graph tests now asks if cells at a similar position in the trajectory have similar expression

cds_no_MAIT_no_90_DE_PT_res_ids <- row.names(subset(cds_no_MAIT_no_90_DE_PT_res, q_value < 0.05))

# Collect trajectory variable genes into modules
cds_no_MAIT_no_90_DE_PT_gene_module_df <- find_gene_modules(cds_no_MAIT_no_9[cds_no_MAIT_no_90_DE_PT_res_ids,], resolution=c(10^seq(-6,-1)))

# Then aggregate the module scores and plot as a heatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_no_MAIT_no_9)), 
                                Cluster=colData(cds_no_MAIT_no_9)$Cluster.Name)
cds_no_MAIT_no_90_DE_PT_agg_mat <- aggregate_gene_expression(cds_no_MAIT_no_9, cds_no_MAIT_no_90_DE_PT_gene_module_df, cell_group_df)
row.names(cds_no_MAIT_no_90_DE_PT_agg_mat) <- stringr::str_c("Module ", row.names(cds_no_MAIT_no_90_DE_PT_agg_mat))

pdf(file= file.path(plotDir, "cds_no_MAIT_no_90_DE_pseudotime_modules_heatmap.pdf"), width = 8, height = 12)
pheatmap::pheatmap(cds_no_MAIT_no_90_DE_PT_agg_mat,
                   scale="column", clustering_method="ward.D2")
dev.off()

# also plot on plot_cells
pdf(file= file.path(plotDir, "cds_no_MAIT_no_90_DE_pseudotime_modules_plot_cells.pdf"), width = 20, height = 20)
plot_cells(cds_no_MAIT_no_9,
           genes=cds_no_MAIT_no_90_DE_PT_gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

#### Analyze Differential expression at specific branch points ####

# First choose cells for subset # subset for clusters 5, 6 and 8
cds_no_MAIT_no_9_subset <- choose_cells(cds_no_MAIT_no_9)

subset_pr_test_res <- graph_test(cds_no_MAIT_no_9_subset, neighbor_graph="principal_graph", cores=8)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

# group the genes into modules
gene_module_df <- find_gene_modules(cds_no_MAIT_no_9_subset[pr_deg_ids,], resolution=0.001)

# organize modules by similarity across trajectory
agg_mat <- aggregate_gene_expression(cds_no_MAIT_no_9_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

subset_module_DE_plot_cells <- plot_cells(cds_no_MAIT_no_9_subset,
                                          genes=gene_module_df,
                                          label_cell_groups=FALSE,
                                          show_trajectory_graph=FALSE)
ggsave(subset_module_DE_plot_cells, file = file.path(plotDir, "subset_module_DE_plot_cells_5_6_8.pdf"), height = 12, width = 8 )

# Module 1 is higher in PD1 and could be interesting
# find genes in module 1
subse_module_1_PT <- gene_module_df %>% filter(module == 1)
write.csv(subse_module_1_PT$id, file = file.path(resultDir, "subset_module_1_PT.csv"))


#### Find top marker genes expressed by each cluster with cluster 9 removed ####

# Use top_markers to ask what genes make clusters different from one another
# only re-run top_markers if necessary, it takes a long time to run
marker_test_res <- top_markers(cds_no_MAIT_no_9, 
                               group_cells_by="Cluster.Name")
# marker_test_res contains a number of metrics for how specifically expressed each gene is in each partition. 

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

ggsave(top_markers_expression, file = file.path(plotDir, "top_markers_expression_no9.pdf"), height = 10, width = 8)
# export results
write.csv(top_specific_markers, file = file.path(resultDir, "top_specific_markers_updated_9_26_2023.csv"))
top_specific_markers <- read.csv(file = file.path(resultDir, "top_specific_markers_updated_9_26_2023.csv"))
View(top_specific_markers)

top_specific_markers %>% filter(cell_group %in% c(5,6,8)) %>% View()

## Plotting the expression of some top markers for different clusters


TCM_TEM <-  plot_cells(cds_no_MAIT_no_9, genes = c( "TIGIT","KLRG1", "IL7R", "CCR7", "CD27"), 
                       cell_size=0.5)
ggsave(TCM_TEM, file = file.path(plotDir, "TCM_TEM.pdf"), width = 8, height = 6)

#### Analyze differential expression by cluster using regression analysis ####

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

# test for genes that differ between cluster ex
#gene_fits_ex <- fit_models(cds_no_MAIT_ex , model_formula_str = "~cluster_ex")
#save(gene_fits_ex ,file = file.path(resultDir, "gene_fits_cluster_ex.Rdata"))
load(file = file.path(resultDir, "gene_fits_cluster_ex.Rdata"))

# extract table of coefficients
fit_coefs <- coefficient_table(gene_fits_ex)

# extract the cluster term
cluster_ex_terms <- fit_coefs %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_ex_terms$term)) #  "cluster_exPD1"

# filter for those genes significantly differing between cluster 7 and cluster 5,6,8
cluster_ex_terms_sig <- cluster_ex_terms   %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_ex_terms_sig, file = file.path(resultDir, "cluster_ex_terms_sig.Rdata"))
load(file = file.path(resultDir, "cluster_ex_terms_sig.Rdata"))
rm(gene_fits_ex)

# join category for up or down
cluster_ex_terms_sig_plot <- cluster_ex_terms_sig %>% mutate(up_down = case_when(normalized_effect >=0 ~"PD-1-like",
                                                                                 normalized_effect <0 ~ "CD57-like"))
# add labels for each  - this is the gene list from ATAC figure 3
CD57_PD1_label <- data.frame(gene_short_name = c("IL7R","HAVCR2","CTLA4","KLRG1","FCGR3A",
                                                 "TIGIT","LAG3","PDCD1","CD160","EOMES",
                                                 "TOX", "KIR3DL1","KIR2DS4", "TBX21","KLRD1",
                                                 "LILRB1","KIR3DL2","KIR2DL3","KIR2DL1", "NCR1","KLRF1"),
                             label = c("IL7R","HAVCR2","CTLA4","KLRG1","FCGR3A",
                                       "TIGIT","LAG3","PDCD1","CD160","EOMES",
                                       "TOX", "KIR3DL1","KIR2DS4", "TBX21","KLRD1",
                                       "LILRB1","KIR3DL2","KIR2DL3","KIR2DL1","NCR1","KLRF1"))
cluster_ex_terms_sig_plot <- left_join(cluster_ex_terms_sig_plot,CD57_PD1_label)

# plot as volcano plot
cluster_ex_terms_sig_volcano <- ggplot(cluster_ex_terms_sig, aes(x = normalized_effect, y = -log10(q_value))) + geom_point() 
cluster_ex_terms_sig_volcano_non_log <- ggplot(cluster_ex_terms_sig, aes(x = normalized_effect, y = q_value)) + geom_point() 
ggsave(cluster_ex_terms_sig_volcano, file=file.path(plotDir, "cluster_ex_terms_sig_volcano.pdf"), device = "pdf")

# format the volcano plot to highlight
cluster_ex_terms_sig_plot_highlight <- ggplot(cluster_ex_terms_sig_plot, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                             label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 5, max.overlaps = 10 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#ba4d4cc7", "#fdbf6fff")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16),
        legend.position="bottom") +
  labs(x = "log2 Fold Change",y = "-log10 q-value")

# add annotation under the columns

text_low <- textGrob("Up in CD57+", gp=gpar(fontsize=16, fontface="bold"))
text_high <- textGrob("Up in PD1+", gp=gpar(fontsize=16, fontface="bold"))

cluster_ex_terms_sig_plot_highlight_plot <- cluster_ex_terms_sig_plot_highlight + 
  annotation_custom(text_high,xmin=3,xmax=3,ymin=-50,ymax=-50) + 
  annotation_custom(text_low,xmin=-2,xmax=-2,ymin=-50,ymax=-50) +
  coord_cartesian( clip="off")

# export plot
ggsave(plot = cluster_ex_terms_sig_plot_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_highlight_plot.pdf"),device = "pdf", width = 7, height = 7)
ggsave(plot = cluster_ex_terms_sig_plot_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_highlight_plot.svg"),device = "svg", width = 7, height = 7)


# filter for up and down
cluster_ex_terms_sig_up <- cluster_ex_terms_sig  %>% filter(normalized_effect >=0)
cluster_ex_terms_sig_down <- cluster_ex_terms_sig  %>% filter(normalized_effect <0)

# negative = higher in CD57
# positive = higher in PD1

# write.csv
write.csv(cluster_ex_terms_sig_up, file = file.path(resultDir, "cluster_ex_terms_sig_up.csv"))
write.csv(cluster_ex_terms_sig_down, file = file.path(resultDir, "cluster_ex_terms_sig_down.csv"))

# export all those with higher effect greater than 1
cluster_ex_terms_sig_up_1 <- cluster_ex_terms_sig_up %>% filter(normalized_effect >= 1.0)
cluster_ex_terms_sig_down_1 <- cluster_ex_terms_sig_down %>% filter(normalized_effect <= -1.0)

write.csv(cluster_ex_terms_sig_up_1, file = file.path(resultDir, "cluster_ex_terms_sig_up_1.csv"))
write.csv(cluster_ex_terms_sig_down_1, file = file.path(resultDir, "cluster_ex_terms_sig_down_1.csv"))

# repeat for 0.5
cluster_ex_terms_sig_up_0.5 <- cluster_ex_terms_sig_up %>% filter(normalized_effect >= 0.5)
cluster_ex_terms_sig_down_0.5 <- cluster_ex_terms_sig_down %>% filter(normalized_effect <= -0.5)

write.csv(cluster_ex_terms_sig_up_0.5, file = file.path(resultDir, "cluster_ex_terms_sig_up_0.5.csv"))
write.csv(cluster_ex_terms_sig_down_0.5, file = file.path(resultDir, "cluster_ex_terms_sig_down_0.5.csv"))

# plot a few of the highest and lowest genes by normalized effect to make sure I am interpretting up and down correctly
cds_no_MAIT_ex_subset <- cds_no_MAIT_ex[row.names(subset(rowData(cds_no_MAIT_ex),
                                                         gene_short_name %in% c("TIGIT","KLRG1", "GZMK", "KIR3DL2","KLRF1", "CTLA4","LAG3"))),]
plot_genes_violin(cds_no_MAIT_ex_subset, group_cells_by="cluster_ex", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

### Check overlap with atac lists


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

### Replot volcano plot highlighting genes from Daniel_Wherry_2022_Figure2b_Texterm_TexKLR
Daniel_Wherry_2022_Figure2b_Texterm_TexKLR <- readxl::read_xlsx("/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/RAW_DATA/P452-1/Published_KIR_gene_sets/Daniel_Wherry_2022_Figure2b_Texterm_TexKLR.xlsx")
Daniel_Wherry_2022_Figure2b_Texterm_TexKLR <- Daniel_Wherry_2022_Figure2b_Texterm_TexKLR %>%
  dplyr::rename(gene_short_name = gene_name) %>% mutate(label = gene_short_name) %>% dplyr::select(-Tex_subset)

cluster_ex_terms_sig_plot_Daniel <- left_join(cluster_ex_terms_sig_plot,Daniel_Wherry_2022_Figure2b_Texterm_TexKLR)

# format the volcano plot to highlight
cluster_ex_terms_sig_plot_Daniel_highlight <- ggplot(cluster_ex_terms_sig_plot_Daniel, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                           label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 5, max.overlaps = 10 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#ba4d4cc7", "#fdbf6fff")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16),
        legend.position="bottom") +
  labs(x = "log2 Fold Change",y = "-log10 q-value")

# add annotation under the columns

text_low <- textGrob("Up in CD57+", gp=gpar(fontsize=16, fontface="bold"))
text_high <- textGrob("Up in PD1+", gp=gpar(fontsize=16, fontface="bold"))

cluster_ex_terms_sig_plot_Daniel_highlight_plot <- cluster_ex_terms_sig_plot_Daniel_highlight + 
  annotation_custom(text_high,xmin=3,xmax=3,ymin=-50,ymax=-50) + 
  annotation_custom(text_low,xmin=-2,xmax=-2,ymin=-50,ymax=-50) +
  coord_cartesian( clip="off")

# export plot
ggsave(plot = cluster_ex_terms_sig_plot_Daniel_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_Daniel_highlight_plot.pdf"),device = "pdf", width = 7, height = 7)
ggsave(plot = cluster_ex_terms_sig_plot_Daniel_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_Daniel_highlight_plot.svg"),device = "svg", width = 7, height = 7)


#### Plotting curated gene list

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

# export plot
ggsave(plot = cluster_ex_terms_sig_plot_curated_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_curated_highlight_plot.pdf"),device = "pdf", width = 7, height = 6)


## Now highlight all genes
# format the volcano plot to highlight all genes
cluster_ex_terms_sig_plot_all_highlight <- ggplot(cluster_ex_terms_sig_plot, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                 label = gene_short_name)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 2, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#ba4d4cc7", "#fdbf6fff")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16),
        legend.position="bottom") +
  labs(x = "log2 Fold Change",y = "-log10 q-value")

# add annotation under the columns

text_low <- textGrob("Up in CD57+", gp=gpar(fontsize=16, fontface="bold"))
text_high <- textGrob("Up in PD1+", gp=gpar(fontsize=16, fontface="bold"))

cluster_ex_terms_sig_plot_all_highlight_plot <- cluster_ex_terms_sig_plot_all_highlight + 
  annotation_custom(text_high,xmin=3,xmax=3,ymin=-50,ymax=-50) + 
  annotation_custom(text_low,xmin=-2,xmax=-2,ymin=-50,ymax=-50) +
  coord_cartesian( clip="off")

# export plot
ggsave(plot = cluster_ex_terms_sig_plot_all_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_all_highlight_plot.pdf"),device = "pdf", width = 10, height = 10)


#### Analyze differential expression of all CD57 clusters versus all other clusters ####

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

# test for genes that differ between CD57 and all other
gene_fits <- fit_models(cds_no_MAIT_no_9 , model_formula_str = "~cluster_CD57")

# extract table of coefficients
fit_coefs <- coefficient_table(gene_fits)

# extract the response term
cluster_terms <- fit_coefs %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_terms$term)) #  
#save(cluster_terms, file = file.path(resultDir, "cluster_terms_CD57_other.Rdata"))

# filter for those genes significantly differing between CD57 and all other
cluster_terms_sig <- cluster_terms  %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
save(cluster_terms_sig, file = file.path(resultDir, "cluster_terms_sig_CD57_other.Rdata"))
#load(file = file.path(resultDir, "cluster_terms_sig_CD57_other.Rdata"))
rm(gene_fits)

write.csv(cluster_terms_sig, file = file.path(resultDir, "cluster_terms_sig_CD57_other.csv"))

#### Analyze differential expression in exhausted clusters by response ####

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

# test for genes that differ between response and cluster ex
#gene_fits_response <- fit_models(cds_no_MAIT_ex , model_formula_str = "~Response")
#save(gene_fits_response ,file = file.path(resultDir, "gene_fits_cluster_response.Rdata"))
load(file = file.path(resultDir, "gene_fits_response.Rdata"))

# extract table of coefficients
fit_coefs_response <- coefficient_table(gene_fits_response)

# extract the response term
cluster_response_terms <- fit_coefs_response %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_response_terms$term)) #  "ResponseR"

# filter for those genes significantly differing between cluster 7 and cluster 5,6,8
cluster_response_terms_sig <- cluster_response_terms   %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_response_terms_sig, file = file.path(resultDir, "cluster_response_terms_sig.Rdata"))
load(file = file.path(resultDir, "cluster_response_terms_sig.Rdata"))
rm(gene_fits_ex)

write.csv(cluster_response_terms_sig, file = file.path(resultDir, "cluster_response_terms_sig.csv"))

### Repeat with just the CD57 populations and just the PD1 population

# use the data subset to only include those clusters after cluster 4
cds_no_MAIT_ex_PD1 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(7)))]

# test for genes that differ between response
gene_fits_response_PD1 <- fit_models(cds_no_MAIT_ex_PD1 , model_formula_str = "~Response")

# extract table of coefficients
fit_coefs_response_PD1 <- coefficient_table(gene_fits_response_PD1)

# extract the response term
cluster_response_terms_PD1 <- fit_coefs_response_PD1 %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_response_terms_PD1$term)) #  "ResponseR"

# filter for those genes significantly differing between R and NR in PD1 cluster
cluster_response_terms_PD1_sig <- cluster_response_terms_PD1   %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
save(cluster_response_terms_PD1_sig, file = file.path(resultDir, "cluster_response_terms_PD1_sig.Rdata"))
load(file = file.path(resultDir, "cluster_response_terms_PD1_sig.Rdata"))
rm(gene_fits_response_PD1)

write.csv(cluster_response_terms_PD1_sig, file = file.path(resultDir, "cluster_response_terms_PD1_sig.csv"))
nrow(cluster_response_terms_PD1_sig)

# Repeat for CD57
# use the data subset to only include CD57 clusterrs
cds_no_MAIT_ex_CD57 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,6,8)))]
levels(as.data.frame(colData(cds_no_MAIT_ex_CD57)$Response))
save(cds_no_MAIT_ex_CD57, file = file.path(resultDir, "cds_no_MAIT_ex_CD57.Rdata"))

# test for genes that differ between response
gene_fits_response_CD57 <- fit_models(cds_no_MAIT_ex_CD57 , model_formula_str = "~Response")

# extract table of coefficients
fit_coefs_response_CD57 <- coefficient_table(gene_fits_response_CD57)

# extract the response term
cluster_response_terms_CD57 <- fit_coefs_response_CD57 %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_response_terms_CD57$term)) #  "ResponseR"

# filter for those genes significantly differing between R and NR in CD57 cluster
cluster_response_terms_CD57_sig <- cluster_response_terms_CD57   %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
save(cluster_response_terms_CD57_sig, file = file.path(resultDir, "cluster_response_terms_CD57_sig.Rdata"))
load(file = file.path(resultDir, "cluster_response_terms_CD57_sig.Rdata"))
rm(gene_fits_response_CD57)

write.csv(cluster_response_terms_CD57_sig, file = file.path(resultDir, "cluster_response_terms_CD57_sig.csv"))


### Compare and contrast the genes that are increased in responders
cluster_response_terms_CD57_sig_up <- cluster_response_terms_CD57_sig %>% filter(up_down == "up")
cluster_response_terms_PD1_sig_up <- cluster_response_terms_PD1_sig %>% filter(up_down == "up")

cluster_response_terms_CD57_sig %>% dplyr::count(up_down)
nrow(cluster_response_terms_CD57_sig) # 635
cluster_response_terms_PD1_sig %>% dplyr::count(up_down)
nrow(cluster_response_terms_PD1_sig) # 204

# make list of overlaps
R_NR_overlap <- 
  list(cluster_response_terms_CD57_sig_up = unique(cluster_response_terms_CD57_sig_up$gene_short_name ),
       cluster_response_terms_PD1_sig_up = unique(cluster_response_terms_PD1_sig_up$gene_short_name))
R_NR_overlap_overlap_venn <- gplots::venn(R_NR_overlap, show.plot = FALSE)
lengths(attributes(R_NR_overlap_overlap_venn)$intersections)

length(unique(cluster_response_terms_CD57_sig_up$gene_short_name )) # 475
length(unique(cluster_response_terms_PD1_sig_up$gene_short_name)) # 132

# Input in the form of a named numeric vector
R_NR_overlap_fit <- eulerr::euler(c(A = 73, 
                                    B = 416,
                                    "A&B" = 59))

# plot the resulting Venn diagram
pdf(file = file.path(plotDir, "R_NR_overlap_fit_Venn_eulerr.pdf"), height = 3, width = 3)
plot(R_NR_overlap_fit,
     quantities = TRUE,
     fill = "transparent",
     lty = 1:3,
     labels = c("PD-1+ R", "CD57+ R"))
dev.off()

# which genes are in this overlap?
CD57_PD1_up_shared_up_overlap <- data.frame(gene_name = attributes(R_NR_overlap_overlap_venn)$intersections$`cluster_response_terms_CD57_sig_up:cluster_response_terms_PD1_sig_up`)
write.csv(CD57_PD1_up_shared_up_overlap , file= file.path(resultDir, "CD57_PD1_up_shared_up_overlap.csv"))

# which genes unique to CD57 up 
CD57_R_up_unique <- data.frame(gene_name = attributes(R_NR_overlap_overlap_venn)$intersections$`cluster_response_terms_CD57_sig_up`)
write.csv(CD57_R_up_unique, file= file.path(resultDir, "CD57_R_up_unique.csv"))

# which genes unique to module 1
PD1_R_up_unique <- data.frame(gene_name = attributes(R_NR_overlap_overlap_venn)$intersections$`cluster_response_terms_PD1_sig_up`)
write.csv(PD1_R_up_unique, file= file.path(resultDir, "PD1_R_up_unique.csv"))


### Plot volcano plots of Differential expression by response for PD1 and CD57 clusters
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

# export plot
ggsave(plot = cluster_ex_terms_PD1_sig_plot_curated_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_PD1_sig_plot_curated_highlight_plot.pdf"),device = "pdf", width = 6, height = 6)

## Repeat for CD57 cluster plot
CD57_curated_label <- data.frame(gene_short_name = c( "KIR3DL2", "KIR2DL3",
                                                      "KLRC1", "KLRB1",
))
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

# export plot
ggsave(plot = cluster_ex_terms_CD57_sig_plot_curated_highlight_plot, 
       file = file.path(plotDir, "cluster_ex_terms_CD57_sig_plot_curated_highlight_plot.pdf"),device = "pdf", width = 6, height = 6)


#### Compare whether DEGs between R and NR account for additional heterogeneity ####

load(file = file.path(resultDir, "cluster_response_terms_CD57_sig.Rdata"))

cluster_response_terms_CD57_sig_compare <- cluster_response_terms_CD57_sig %>% dplyr::rename(gene_name = gene_short_name)

length(unique(cluster_response_terms_CD57_sig_compare$gene_name)) # 635

## Compare R vs NR list to CD57 vs PD1 lists # 264 overlaps between these lists

CD57_R_vs_NR_PD1_vs_CD57_overlap <- cluster_response_terms_CD57_sig[cluster_response_terms_CD57_sig$gene_short_name %in% unique(cluster_ex_terms_sig $gene_short_name),]
CD57_R_vs_NR_PD1_vs_CD57_overlap$overlap <- "yes"
cluster_response_terms_CD57_sig[cluster_response_terms_CD57_sig$gene_short_name %in% unique(cluster_ex_terms_sig $gene_short_name),] %>% View()

## Compare with Module 1 genes

# which genes are in this overlap? # 8 genes
cluster_response_terms_CD57_sig_compare[cluster_response_terms_CD57_sig_compare$gene_name %in% module1_shared_up_overlap$gene_name,]


# which genes unique to shared up # 25 genes different in R and NR are in ATAC shared up 
cluster_response_terms_CD57_sig_compare[cluster_response_terms_CD57_sig_compare$gene_name %in% module1_shared_up_unique_atac$gene_name,]

# which genes unique to module 1 # 48 genes that differ between R and NR are in module 1..supporting that some of the patterns of sharing we
# see are in fact based on R and NR differences..I did already show some sharing between PD1 R vs NR and CD57 R vs NR
cluster_response_terms_CD57_sig_compare[cluster_response_terms_CD57_sig_compare$gene_name %in% module1_shared_up_unique_module1$gene_name,]

### Replot volcano plot highlighting those genes differing in PD1 vs CD57 that are also different between R and NR...

cluster_ex_terms_sig_plot_R_NR_overlap <- cluster_ex_terms_sig_plot %>% left_join(.,unique(CD57_R_vs_NR_PD1_vs_CD57_overlap[,c("gene_short_name","overlap")])) %>%
  mutate(up_down_overlap = case_when(overlap == "yes" ~ "R vs. NR\nOverlap",
                                     up_down == "PD-1-like"    ~ "PD-1-like",
                                     up_down == "CD57-like"    ~ "CD57-like"))

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

# export plot
ggsave(plot = cluster_ex_terms_sig_plot_curated_highlight_plot_overlap, 
       file = file.path(plotDir, "cluster_ex_terms_sig_plot_curated_highlight_plot_overlap.pdf"),device = "pdf", width = 9, height = 6)



#### Plot CD57 KIR gene expression between responders and non-responders ####

# calculate mean across all cells per donor
cds_no_MAIT_ex_CD57_KIR_expression <- log10(as.data.frame(t(exprs(cds_no_MAIT_ex_CD57[rowData(cds_no_MAIT_ex_CD57)$gene_short_name %in% 
                                                                                        c("KIR3DL2","KIR2DL3", "KLRC1","KLRB1"),])) )) %>% 
  rownames_to_column("Cell.ID") 

CD57_anno_data <- as.data.frame(colData(cds_no_MAIT_ex_CD57)) #%>% select(-CCR7, -CD28, -IL7R, -IFNG, -MKI67) # this was originally commented out by Kirsten
CD57_anno_data$Cell.ID <- row.names(colData(cds_no_MAIT_ex_CD57))

# Calculate the mean of genes of interest across all cells per donor where expression above 0, removing cluster 9
CD57_anno_data_anno <-  cds_no_MAIT_ex_CD57_KIR_expression %>%
  merge(CD57_anno_data, by = "Cell.ID") %>%
  group_by(Donor.ID, Response) %>%
  summarise_at(vars(c("KIR3DL2","KIR2DL3", "KLRC1","KLRB1")), ~mean(. > 0)) %>%
  pivot_longer(3:6,names_to = "gene_name", values_to = "mean_log10exp")

# calculate mean and sd for boxplot
CD57_anno_data_anno <- CD57_anno_data_anno %>% group_by(Response, gene_name) %>% mutate(mean = mean(mean_log10exp), sd = sd(mean_log10exp))

min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# Plot as boxplot
CD57_anno_data_anno_boxplot <- ggplot(CD57_anno_data_anno, aes(x = Response, y = mean_log10exp, color = Response)) + 
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=2) +
  facet_grid(.~gene_name) +
  theme_minimal() +
  scale_color_manual(values = R_NR_colors) +
  stat_compare_means() +
  theme(text = element_text(size = 20)) +
  labs(y = "Mean log10 Expr.", title = "CD57+ Population Differential Expression")

ggsave(CD57_anno_data_anno_boxplot, file = file.path(plotDir, "CD57_anno_data_anno_boxplot.pdf"), device = "pdf", width = 9,
       height = 5)

#### Explore Markers differentiating individual CD57 clusters ####

cds_CD57 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,6,8)))]


# Use top_markers to ask what genes make clusters different from one another
# only re-run top_markers if necessary, it takes a long time to run

marker_test_res_CD57 <- top_markers(cds_CD57, 
                                    group_cells_by="Cluster.Name")
# marker_test_res contains a number of metrics for how specifically expressed each gene is in each partition. 

# Rank markers according to specificity markers, pseudo_R2 is one of these measures
top_specific_markers_CD57 <- marker_test_res_CD57 %>%
  filter(specificity >= 0.5) %>%
  group_by(cell_group) 

top_specific_markers_ids_CD57 <- unique(top_specific_markers_CD57 %>% pull(gene_short_name))

# how many shared across all three
#top_specific_markers_CD57 %>% ungroup() %>% dplyr::count(gene_id) %>% arrange(desc(n)) %>% View()
# only 1 was a top markers across all three
top_specific_markers_CD57 %>% filter(cell_group == 5) %>% ungroup() %>% dplyr::select(gene_id)

# Now we can plot the expression and fraction of cells that express each 
# marker in each group with the plot_genes_by_group function
top_markers_expression_CD57 <- plot_genes_by_group(cds_CD57,
                                                   top_specific_markers_ids_CD57,
                                                   group_cells_by="Cluster.Name",
                                                   ordering_type="cluster_row_col",
                                                   max.size=5) + theme(text = element_text(size = 16))

ggsave(top_markers_expression_CD57, file = file.path(plotDir, "top_markers_expression_CD57.pdf"),
       height = 8, width = 5)

# export results
write.csv(top_specific_markers_CD57, file = file.path(resultDir, "top_specific_markers_CD57.csv"))

#### Explore differential expression between CD57 clusters ####

# assess differential expression between cluster 5 and 6
cds_no_MAIT_ex_5_6 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,6)))]

# test for genes that differ between cluster ex
#gene_fits_ex_5_6 <- fit_models(cds_no_MAIT_ex_5_6 , model_formula_str = "~Cluster.Name")
#save(gene_fits_ex_5_6,file = file.path(resultDir, "gene_fits_ex_5_6.Rdata"))
#load(file = file.path(resultDir, "gene_fits_ex_5_6.Rdata"))

# extract table of coefficients
fit_coefs_5_6 <- coefficient_table(gene_fits_ex_5_6)

# extract the cluster term
cluster_ex_terms_5_6 <- fit_coefs_5_6 %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_ex_terms_5_6$term)) # "Cluster.Name6"

# filter for those genes significantly differing
cluster_ex_terms_5_6_select <- cluster_ex_terms_5_6 %>% dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
cluster_ex_term_5_6_sig <- cluster_ex_terms_5_6_select  %>% filter (q_value < 0.05) 
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_ex_terms_5_6_select, file = file.path(resultDir, "cluster_ex_terms_5_6.Rdata"))
#save(cluster_ex_term_5_6_sig, file = file.path(resultDir, "cluster_ex_term_5_6_sig.Rdata"))
load(file = file.path(resultDir, "cluster_ex_term_5_6_sig.Rdata"))
write.csv(cluster_ex_term_5_6_sig, file = file.path(resultDir, "cluster_ex_term_6_vs_5_sig.csv"), row.names = FALSE)
rm(gene_fits_ex_5_6)
nrow(cluster_ex_term_5_6_sig) # 1124

## Repeat for 5,8
cds_no_MAIT_ex_5_8 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(5,8)))]

# test for genes that differ between cluster ex
gene_fits_ex_5_8 <- fit_models(cds_no_MAIT_ex_5_8 , model_formula_str = "~Cluster.Name")
#save(gene_fits_ex_5_8,file = file.path(resultDir, "gene_fits_ex_5_8.Rdata"))
#load(file = file.path(resultDir, "gene_fits_ex_5_8.Rdata"))

# extract table of coefficients
fit_coefs_5_8 <- coefficient_table(gene_fits_ex_5_8)

# extract the cluster term
cluster_ex_terms_5_8 <- fit_coefs_5_8 %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_ex_terms_5_8$term)) #  "Cluster.Name8"

# filter for those genes significantly differing 
cluster_ex_terms_5_8_select <- cluster_ex_terms_5_8 %>%
  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
cluster_ex_term_5_8_sig <- cluster_ex_terms_5_8_select %>% filter (q_value < 0.05) 
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_ex_terms_5_8_select, file = file.path(resultDir, "cluster_ex_terms_5_8.Rdata"))
#save(cluster_ex_term_5_8_sig, file = file.path(resultDir, "cluster_ex_term_5_8_sig.Rdata"))
load(file = file.path(resultDir, "cluster_ex_term_5_8_sig.Rdata"))
load(file = file.path(resultDir, "cluster_ex_term_5_8.Rdata"))
write.csv(cluster_ex_term_5_8_sig, file = file.path(resultDir, "cluster_ex_term_8_vs_5_sig.csv"), row.names = FALSE)
rm(gene_fits_ex_5_8)
nrow(cluster_ex_term_5_8_sig) # 458

## Repeat for 6,8

cds_no_MAIT_ex_6_8 <- cds_no_MAIT[,colnames(cds_no_MAIT) %in% row.names(colData(cds_no_MAIT) %>% as.data.frame() %>% filter(Cluster.Name %in% c(6,8)))]

# test for genes that differ between cluster ex
gene_fits_ex_6_8 <- fit_models(cds_no_MAIT_ex_6_8 , model_formula_str = "~Cluster.Name")
#save(gene_fits_ex_6_8,file = file.path(resultDir, "gene_fits_ex_6_8.Rdata"))
#load(file = file.path(resultDir, "gene_fits_ex_6_8.Rdata"))

# extract table of coefficients
fit_coefs_6_8 <- coefficient_table(gene_fits_ex_6_8)

# extract the cluster term
cluster_ex_terms_6_8 <- fit_coefs_6_8 %>% filter(term != "(Intercept)") 
levels(as.factor(cluster_ex_terms_6_8$term)) #  "Cluster.Name8"

# filter for those genes significantly differing 
cluster_ex_terms_6_8_select <- cluster_ex_terms_6_8 %>% 
  dplyr::select(gene_short_name, term, q_value, estimate, normalized_effect)
cluster_ex_term_6_8_sig <- cluster_ex_terms_6_8_select  %>% filter (q_value < 0.05) 
# normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307
#save(cluster_ex_term_6_8_select, file = file.path(resultDir, "cluster_ex_term_6_8.Rdata"))
#save(cluster_ex_term_6_8_sig, file = file.path(resultDir, "cluster_ex_term_6_8_sig.Rdata"))
write.csv(cluster_ex_term_6_8_sig, file = file.path(resultDir, "cluster_ex_term_8_vs_6_sig.csv"), row.names = FALSE)
load(file = file.path(resultDir, "cluster_ex_term_6_8_sig.Rdata"))
rm(gene_fits_ex_6_8)

nrow(cluster_ex_term_6_8_sig)  # 475

#### PLOT VOLCANO PLOT OF CD57 CONTRAST RESULTS ####

### Extract genes in each contrast that are < -1 or > 1
cluster_ex_term_6_8_sig_LFC_1 <- cluster_ex_term_6_8_sig %>% filter(normalized_effect >= 1 | normalized_effect <= -1) %>%
  mutate(contrast = "8_vs_6", label = gene_short_name)
nrow(cluster_ex_term_6_8_sig_LFC_1) # 134
cluster_ex_term_5_8_sig_LFC_1 <- cluster_ex_term_5_8_sig %>% filter(normalized_effect >= 1 | normalized_effect <= -1) %>% 
  mutate(contrast = "8_vs_5", label = gene_short_name)
nrow(cluster_ex_term_5_8_sig_LFC_1) # 145
cluster_ex_term_5_6_sig_LFC_1 <- cluster_ex_term_5_6_sig %>% filter(normalized_effect >= 1 | normalized_effect <= -1) %>% 
  mutate(contrast = "6_vs_5", label = gene_short_name)
nrow(cluster_ex_term_5_6_sig_LFC_1) # 82

### Label all genes with LFC > 1 in volcano plots
# 5 =  "#b84c7d"
# 6 ="#7f63b8"
# 8 = "#56ae6c"

## 6 vs 8
cluster_ex_term_6_8_LFC_label <- left_join(cluster_ex_terms_6_8_select,cluster_ex_term_6_8_sig_LFC_1) %>%
  mutate(up_down = case_when(normalized_effect >=0 & q_value < 0.05 ~"Cluster 8",
                             normalized_effect <0 & q_value < 0.05 ~ "Cluster 6",
                             q_value > 0.05 ~ "N.S."))

cluster_ex_terms_6_8_LFC_label_plot <- ggplot(cluster_ex_term_6_8_LFC_label, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                 label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 3, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#7f63b8", "#56ae6c", "gray70")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16),
        legend.position="bottom") +
  scale_y_continuous(limits = c(0,450)) +
  geom_vline(xintercept = 0) +
  labs(x = "log2 Fold Change",y = "-log10 q-value", title = "DEG Cluster 8 vs Cluster 6")

## 5 vs 8 
cluster_ex_term_5_8_LFC_label <- left_join(cluster_ex_terms_5_8_select,cluster_ex_term_5_8_sig_LFC_1) %>%
  mutate(up_down = case_when(normalized_effect >=0 & q_value < 0.05 ~"Cluster 8",
                             normalized_effect <0 & q_value < 0.05 ~ "Cluster 5",
                             q_value > 0.05 ~ "N.S."))

cluster_ex_term_5_8_LFC_label_plot <- ggplot(cluster_ex_term_5_8_LFC_label, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 3, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#b84c7d", "#56ae6c", "gray70")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16),
        legend.position="bottom") +
  scale_y_continuous(limits = c(0,450)) +
  geom_vline(xintercept = 0) +
  labs(x = "log2 Fold Change",y = "-log10 q-value", title = "DEG Cluster 8 vs Cluster 5")

## 5 vs 6
cluster_ex_term_5_6_LFC_label <- left_join(cluster_ex_terms_5_6_select,cluster_ex_term_5_6_sig_LFC_1) %>%
  mutate(up_down = case_when(normalized_effect >=0 & q_value < 0.05 ~"Cluster 6",
                             normalized_effect <0 & q_value < 0.05 ~ "Cluster 5",
                             q_value > 0.05 ~ "N.S."))
cluster_ex_term_5_6_LFC_label_plot <- ggplot(cluster_ex_term_5_6_LFC_label, aes(x = normalized_effect, y = -log10(q_value), color = up_down,
                                                                                label = label)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 3, max.overlaps = 50 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#b84c7d", "#7f63b8", "gray70")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16),
        legend.position="bottom") +
  scale_y_continuous(limits = c(0,450)) +
  geom_vline(xintercept = 0) +
  labs(x = "log2 Fold Change",y = "-log10 q-value", title = "DEG Cluster 6 vs Cluster 5")

# export all three plots
all_CD57_DEG_LFC_1 <- ggarrange(cluster_ex_terms_6_8_LFC_label_plot, cluster_ex_term_5_8_LFC_label_plot,cluster_ex_term_5_6_LFC_label_plot, nrow = 2, ncol = 2)
ggsave(all_CD57_DEG_LFC_1, file = file.path(plotDir,"all_CD57_DEG_LFC_1.pdf"), width = 15, height = 15)

#### Plot heatmap of top 10 genes that differentiate each CD57 contrast ####


# extract top 10 terms in each contrast
cluster_ex_term_6_8_sig_LFC_1_up <- cluster_ex_term_6_8_sig_LFC_1 %>% top_n(10, normalized_effect)
cluster_ex_term_6_8_sig_LFC_1_down <- cluster_ex_term_6_8_sig_LFC_1 %>% top_n(-10, normalized_effect)
cluster_ex_term_6_8_sig_LFC_1_top <- rbind(cluster_ex_term_6_8_sig_LFC_1_up, cluster_ex_term_6_8_sig_LFC_1_down )

cluster_ex_term_5_8_sig_LFC_1_up <- cluster_ex_term_5_8_sig_LFC_1 %>% top_n(10, normalized_effect)
cluster_ex_term_5_8_sig_LFC_1_down <- cluster_ex_term_5_8_sig_LFC_1 %>% top_n(-10, normalized_effect)
cluster_ex_term_5_8_sig_LFC_1_top <- rbind(cluster_ex_term_5_8_sig_LFC_1_up, cluster_ex_term_5_8_sig_LFC_1_down )

cluster_ex_term_5_6_sig_LFC_1_up <- cluster_ex_term_5_6_sig_LFC_1 %>% top_n(10, normalized_effect)
cluster_ex_term_5_6_sig_LFC_1_down <- cluster_ex_term_5_6_sig_LFC_1 %>% top_n(-10, normalized_effect)
cluster_ex_term_5_6_sig_LFC_1_top <- rbind(cluster_ex_term_5_6_sig_LFC_1_up, cluster_ex_term_5_6_sig_LFC_1_down )



#### FIGURE S5 PLOTING ####

## Add correct public masked ID and replot without cluster 9!
load( file = "./EW_T1DAL_Results/cds_no_MAIT_no_9.Rdata")
cds_no_MAIT_no_9_meta <- as.data.frame(colData(cds_no_MAIT_no_9))
cds_no_MAIT_no_9_meta <- left_join(cds_no_MAIT_no_9_meta, T1DAL_ITN_ID_dictionary)
colData(cds_no_MAIT_no_9)$masked_public_PID <- cds_no_MAIT_no_9_meta$masked_public_PID

# replot
cds_no_MAIT_UMAP_donor_no_9_PID <- plot_cells(cds_no_MAIT_no_9, color_cells_by = "masked_public_PID", label_cell_groups = F) + #+ facet_grid(~ Response) + 
  scale_color_manual(values = donor_colors) + 
  #labs(color = "Donor ID") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size=4),title='Donor ID')) 
ggsave(cds_no_MAIT_UMAP_donor_no_9_PID, file = "./FIGURES/cds_no_MAIT_UMAP_Donor.ID_23_01_23.tiff", device = "tiff",
       width = 6, height = 5)

# replot with R vs NR
cds_no_MAIT_UMAP_response_9_group <- plot_cells(cds_no_MAIT_no_9, color_cells_by = "Response", label_cell_groups = F) + #+ facet_grid(~ Response) + 
  scale_color_manual(values = c("#56ae6c","#b0457b")) + 
  theme(text = element_text(size = 20)) 
ggsave(cds_no_MAIT_UMAP_response_9_group, file = "./FIGURES/cds_no_MAIT_UMAP_response_23_09_22.tiff", device = "tiff",
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

pdf(file = "./FIGURES/combined_cluster_percent_total.pdf", height = 12, width = 7)
egg::ggarrange(total_cells_per_cluster, percent_cells_per_donor_per_cluster_plot, 
               percent_cells_per_response_per_cluster_plot,
               percent_cells_per_cluster_per_donor_plot,
               ncol = 1, heights = c(0.3,0.5, 0.5, 0.5)) 
dev.off()

