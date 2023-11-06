#### P362-1 T1DAL 10X QC Re-analysis ####

# This QC generates the clustering output that is saved and used in 
# the subsequent analysis script 

#### Load libraries ####

library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(org.Hs.eg.db)
library(limma)
library(apird)

# Set wd
setwd("/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES")

# set seed
set.seed(42)

#### QC and batch correction for T1DAL memory CD8 10X analysis ####

# Code from K. Diggins 
# 2020-09-02 P362-1 T1DAL 10X analysis - DGEA analysis.Rmd
# 2020-09-02 P362-1 T1DAL 10X analysis - QC.Rmd
# Note that the DGEA analysis also includes the QC steps 

### Load the data ####

P362_libs <- apird::getProjectLibs("P362-1")
P362_anno <- apird::getAnno(P362_libs) 
lib_fc <- getLibIdsWithFcIds(P362_libs)

# Save annotation data for later:
# Kirsten originally saved the data as "P362-2 annotation.Rdata", but it should be P362-1
save(P362_anno, file="/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/P362-1_annotation.Rdata")
load( file="/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/P362-1_annotation.Rdata")

## Load the 10X data - only available on isilon
# filepath on isilon #/Volumes/Bioinformatics/workspace/mrosasco/200713-10XT1DALRvsNR3Batches_P362-1/P362-1_T1DAL_CRv4NoAdapter_aggrNoNorm/outs
# downloaded data locally from isilon 

all.data <- Read10X("./RAW_DATA/P362_1/P362-1_T1DAL_CRv4NoAdapter_aggrNoNorm/outs/filtered_feature_bc_matrix")

# Create Seurat object to perform QC
aggr <- CreateSeuratObject(counts = all.data, min.features = 100, min.cells = 3, project = "P362-1")
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

## load ITN dictionary matching operational donor IDs to the correct public IDs
T1DAL_ITN_ID_dictionary <- readxl::read_xlsx( "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/RAW_DATA/t1dal_mask_id_EW_formatted.xlsx")
colnames(T1DAL_ITN_ID_dictionary) <- c("operational_PID", "masked_public_PID","Donor.ID")
# remove T1DAL from public_PID column
T1DAL_ITN_ID_dictionary$masked_public_PID <- str_remove(T1DAL_ITN_ID_dictionary$masked_public_PID, "T1DAL_")

#### Format for QC in Seurat ####

## Merge all metadata and annotations
metadata <- read.csv(file="./RAW_DATA/P362_1/P362-1_T1DAL_CRv4NoAdapter_aggrNoNorm/outs/aggregation.csv") %>%
  dplyr::select(-molecule_h5) %>%
  merge(P362_anno,by.x="library_id",by.y="libid") 

gemgroup <- sapply(strsplit(rownames(aggr@meta.data), split="-"), "[[", 2) 

aggr <- AddMetaData(object=aggr,
                    metadata=data.frame(gemgroup=gemgroup,
                                        Lib.ID=metadata[gemgroup,"library_id"],
                                        Donor.ID=metadata[gemgroup,"participantID"],
                                        Pool = metadata[gemgroup,"pool"],
                                        Timepoint = metadata[gemgroup,"timepoint"],
                                        Response = metadata[gemgroup,"response"],
                                        Date.Collected = metadata[gemgroup,"dateCollected"],
                                        row.names=rownames(aggr@meta.data)))

mydata_preQC <- aggr

## Stash QC info as part of the project metadata 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mydata_preQC[["percent.mt"]] <- PercentageFeatureSet(mydata_preQC, pattern = "^MT-")

# plot density 
pdf(file = "./FIGURES/percent.mt_plot.pdf")
plot(density(mydata_preQC$percent.mt))
dev.off()

# preview QC metrics
head(mydata_preQC@meta.data, 5)

#### Quality Filtering with Seurat ####

# Visualize QC metrics, and use these to filter cells.
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(mydata_preQC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
## for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mydata_preQC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata_preQC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Plot the QC data

mydata_preQC_plot <- ggplot(mydata_preQC@meta.data,aes(x = percent.mt, y=nFeature_RNA)) +
  geom_point(alpha=0.5,size=1,color="grey50") +
  theme_classic() +
  geom_hline(yintercept=200, size=1, linetype="dotted",color="red") + geom_hline(yintercept=2500, size=1, linetype="dotted",color="red") +
  geom_vline(xintercept=25,size=0.5,linetype="dotted",color="red") + # fixed the y intercept to 25 from original code
  labs(x="Percent Mitochondrial DNA", y="Features per cell") +
  theme(axis.title = element_text(size=15))
ggsave(mydata_preQC_plot, file = "./FIGURES/mydata_preQC_plot.tiff", device = "tiff", width = 5, height = 4)

## Apply the data filter
mydata_postQC <- subset(mydata_preQC, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# originally the 25% MT was commented out, but her QC results powerpoint said she did it so 
# I added it back in - also noticed that this step was in the 2020-03-05 T1DAL 10X - DGEA with Seurat.Rmd script 
print((ncol(mydata_preQC) - ncol(mydata_postQC))/ncol(mydata_preQC)*100)
# 15.65 # this checks out with what her powerpoint originally said 
plot(density(mydata_postQC@meta.data$percent.mt))

# save post QC seurat data
save(mydata_postQC, file = "./EW_T1DAL_Results/mydata_postQC.RData")

## Log normalize the data FOR PLOTTING
# NOTE THAT THIS LOG NORMALIZATION WAS NOT USED IN THE DOWNSTREAM MONOCLE ANALYSIS
mydata <- NormalizeData(mydata_postQC, normalization.method = "LogNormalize", scale.factor = 10000)
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mydata), 10)
# [1] "S100A8" "S100A9" "LYZ"    "CST3"   "CCL4L2" "FCN1"   "G0S2"   "CCL4"   "CXCL8"  "IFIT2" 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
ggsave(plot2, file = "./FIGURES/standarized_variance_expression.tiff", device = "tiff", width = 7, height = 5)

#### Create Monocle CDS data ####

#subset_metadata <- dplyr::filter(aggr_postQC@meta.data,row.names(aggr_postQC@meta.data) %in% col.names(GetAssayData(aggr_postQC,assay="SCT",slot="scale.data")))
# ^ Kirsten had this line in her code and it was commented out, it refers to postQC data
# choosing to not run it - but preserving code for historical purposes

# Get gemgroup for later? 
gemgroup <- sapply(strsplit(rownames(mydata_postQC@meta.data), split="-"), "[[", 2) 

# Put all metadata together
allmetadata = aggr@meta.data

### Expression matrix
final_matrix = GetAssayData(mydata_postQC,assay="RNA",slot="counts")

## Filter cells and genes according to seurat thresholds
keep_genes = row.names(final_matrix)
keep_cells = colnames(final_matrix)

## Subset Cell metadata
final_meta_data = allmetadata[keep_cells,]
final_meta_data$Percent.MT <- mydata_postQC$percent.mt

## Gene meta data
final_gene_data = data.frame(gene_short_name = keep_genes) #] %>% dplyr::filter(gene_short_name %in% keep_genes)
row.names(final_gene_data) <- row.names(final_matrix)

## Combine into cell data set object
cds <- new_cell_data_set(final_matrix, 
                         cell_metadata = final_meta_data,
                         gene_metadata = final_gene_data)

## Add processing batch ID
colData(cds)$Batch <- colData(cds)$Lib.ID
colData(cds)$Batch <- as.factor(dplyr::recode(colData(cds)$Batch,
                                              "lib48853" = 1,
                                              "lib48854" = 1,
                                              "lib48855" = 1,
                                              "lib48856" = 1,
                                              "lib48857" = 1,
                                              "lib48858" = 1,
                                              "lib48859" = 2,
                                              "lib48860" = 2,
                                              "lib48861" = 2,
                                              "lib48862" = 2,
                                              "lib48863" = 2,
                                              "lib48864" = 2))

### Preprocess and plot 
cds <- preprocess_cds(cds, num_dim = 100) # this was commented out by Kirsten, decided to keep this way
# default is log normalization

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, preprocess_method = "PCA") 

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

# Plot the cells without trajectory learned
plot_cells(cds, color_cells_by = "Donor.ID", label_cell_groups = F) # donor ID relatively evenly mixed

# Output plots both by louvain and normal leiden clustering
percent.MT_UMAP_plot <- plot_cells(cds, color_cells_by = "Percent.MT", label_cell_groups = F)
#ggsave(percent.MT_UMAP_plot, file = "./FIGURES/percent.MT_UMAP_plot.tiff", device = "tiff", width = 7, height = 7)  
#ggsave(percent.MT_UMAP_plot, file = "./FIGURES/percent.MT_UMAP_plot_norm_none_12_16_21.tiff", device = "tiff", width = 7, height = 7)  
#ggsave(percent.MT_UMAP_plot, file = "./FIGURES/percent.MT_UMAP_plot_norm_log_12_16_21.tiff", device = "tiff", width = 7, height = 7)  
#ggsave(percent.MT_UMAP_plot, file = "./FIGURES/percent.MT_UMAP_plot_seurat_norm_12_16_21.tiff", device = "tiff", width = 7, height = 7)  
ggsave(percent.MT_UMAP_plot, file = "./FIGURES/percent.MT_UMAP_plot_monocle_default_norm_12_16_21.tiff", device = "tiff", width = 7, height = 7)  


#there is one cluster on the left that has a lot of high MT cells 

## Manually gate out high percent.MT cells.
cds_post_QC <- choose_cells(cds)
print((ncol(mydata_preQC) - ncol(cds_post_QC))/ncol(mydata_preQC)*100) # 33.95

# Plot cell by date collected and donor ID
post_QC_MT_date_collected_UMAP <- plot_cells(cds_post_QC, color_cells_by = "Date.Collected", label_cell_groups = F) # inherently color_cells_by="cluster"
post_QC_MT_donor_UMAP <- plot_cells(cds_post_QC, color_cells_by = "Donor.ID", label_cell_groups = F) # inherently color_cells_by="cluster"

ggsave(post_QC_MT_date_collected_UMAP + post_QC_MT_donor_UMAP , file = "./FIGURES/post_QC_MT_donor_date_collected_UMAPS.tiff", device = "tiff",
       width = 10, height = 6)
#ggsave(post_QC_MT_date_collected_UMAP + post_QC_MT_donor_UMAP , file = "./FIGURES/post_QC_MT_donor_date_collected_UMAPS_norm_none_12_16_21.tiff", device = "tiff",
#       width = 10, height = 6)

#### Batch correction ####

## Batch correct for DONOR.ID
### Kirsten's code originally said to correct by Donor.ID, but her QC analysis 
# summary powerpoint here: 
# "/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_POWERPOINTS_FROM_KIRSTEN/2020-09-08 T1DAL 10X - QC analysis and preprocessing.pptx"
# suggested she corrected by Date.Collected, so I originally did this
# however the final output in terms of number of clusters did not match her output, so
# we both decided to use Donor.ID for batch correction and NOT Date.Collected
cds_batch_corrected = align_cds(cds_post_QC, alignment_group = "Donor.ID") # residual_model_formula_str = "~ Percent.MT") # this was commented out by Kirsten - kept it this way
cds_batch_corrected <- reduce_dimension(cds_batch_corrected) 
cds_batch_corrected <- cluster_cells(cds_batch_corrected, 
                                     resolution = 1e-5) #smaller number = fewer cluster
# Note from Kirsten: the Batchelor correction doesn't actually change the data..it just changes how it looks..so later in downstream
# analysis or other packages for DGEA, make sure you actually perform normalization 

# save my version of this data
save(cds_batch_corrected,file="./EW_T1DAL_Results/2021-12-06_T1DAL_batch_corrected_cds.Rdata")

# Plot cells after batch correction with and without facet grid
batch_corrected_UMAP_date <- plot_cells(cds_batch_corrected, color_cells_by = "Date.Collected", label_cell_groups = F) #+ facet_grid(~ Donor.ID)
batch_corrected_UMAP_donor_facet <- plot_cells(cds_batch_corrected, color_cells_by = "Donor.ID", label_cell_groups = F)

ggsave(batch_corrected_UMAP_date + batch_corrected_UMAP_donor_facet , file = "./FIGURES/batch_corrected_donor_ID_UMAPS.tiff", device = "tiff",
       width = 12, height = 6)
#ggsave(batch_corrected_UMAP_date + batch_corrected_UMAP_donor_facet , file = "./FIGURES/batch_corrected_donor_ID_UMAPS_norm_none_12_16_21.tiff", device = "tiff",
#       width = 12, height = 6)

#### Exclude MAIT cells ####

# Plot CD161 (KLRB1) - the MAIT cell marker
KLRB1_batch_corrected <- plot_cells(cds_batch_corrected, gene = "KLRB1",label_cell_groups = F)
ggsave(KLRB1_batch_corrected, file = "./FIGURES/KLRB1_batch_corrected.tiff", device = "tiff",
       width = 7, height = 6)
#ggsave(KLRB1_batch_corrected, file = "./FIGURES/KLRB1_batch_corrected_norm_none_12_16_21.tiff", device = "tiff",
#       width = 7, height = 6)

# cds MAIT cells 
cds_no_MAIT <- choose_cells(cds_batch_corrected)

## Load Kirstens data below when comparing our two results 
# load("C:/Users/kdiggins/Box/P362-1 T1DAL 10X/P362-1 T1DAL cds object - postQC no MAIT cells.Rdata")
cds_no_MAIT <- reduce_dimension(cds_no_MAIT)

# Re cluster the cells after the MAIT cells are removed - checked resolution with Kirsten
cds_no_MAIT = cluster_cells(cds_no_MAIT, 
                            resolution = 1e-4) #smaller number = fewer cluster # Louvain here looks very strange!

#### Perform Pseudotime analysis ####

cds_no_MAIT <- learn_graph(cds_no_MAIT,use_partition=F, close_loop=F,
                           learn_graph_control=list(ncenter = 140, minimal_branch_len = 6)) # best balance!

## Rename clusters per trajectory:
colData(cds_no_MAIT)$Cluster.Name <- clusters(cds_no_MAIT)
colData(cds_no_MAIT)$Cluster.Name <- as.factor(dplyr::recode(colData(cds_no_MAIT)$Cluster.Name,
                                                             "1" = 2,
                                                             "2" = 4,
                                                             "3" = 5,
                                                             "4" = 6,
                                                             "5" = 3,
                                                             "6" = 1,
                                                             "7"=7,
                                                             "8"=8,
                                                             "9"=9)) # Kirsten only had 6 here originally?

# Save the cds_no_MAIT file for future analysis 
save(cds_no_MAIT,file = "./EW_T1DAL_Results/P362-1 T1DAL cds object - postQC no MAIT cells.Rdata")
# Later places in Kirsten's code this was also called 
# "P362-1 T1DAL cds object - postQC no MAIT cells - more clusters.Rdata" - Kirsten says this is the same 
# as the object above



