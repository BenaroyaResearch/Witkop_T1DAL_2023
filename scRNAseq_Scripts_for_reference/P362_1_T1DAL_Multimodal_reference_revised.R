#### P361 Seurat Multimodal Reference Mapping ####

# Erin Witkop
# Begun September 22nd, 2022

#### LOAD LIBRARIES ####

library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(monocle3)
library(ggpubr)

#### Set paths and load data ####

setwd("/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/Multimodal_Reference_Mapping")

set.seed(42) # same seed as QC script 

# Load saved annotations 
load("../P362-1_annotation.Rdata")

# Load post QC no MAIT cell data
load("../EW_T1DAL_Results/P362-1 T1DAL cds object - postQC no MAIT cells.Rdata")

# set output directories
plotDir <- "../FIGURES/Multimodal_reference_mapping/"
resultDir <- "../EW_T1DAL_Results/"

# Load Seurat object following filtering step 
# not this data has not been batch corrected, normalized, or MAIT cells removed 
load("../EW_T1DAL_Results/mydata_postQC.RData")
# note this data has not yet been normalized 

# inspect metadata to check it transferred correctly
mydata_postQC@meta.data

#### Preprocess and normalize P362_1 RNA data ####

# Save object for normalization under different name (Multi-Modal Reference = MMR)
mydata_postQC_MMR <- mydata_postQC
mydata_postQC_MMR@meta.data
rm(mydata_postQC)

## RNA
DefaultAssay(mydata_postQC_MMR) <- 'RNA'
mydata_postQC_MMR <- mydata_postQC_MMR %>%
  #NormalizeData(mydata_postQC_MMR, normalization.method = "LogNormalize") %>% # skipping log normalization since in the next step SCT transform is recommended
  FindVariableFeatures() %>% 
  # scale the data but include batch effect correction for Donor.ID, which is what was originally done in monocle 
  ScaleData(vars.to.regress = "Donor.ID") %>% 
  RunPCA()

# save this object for later
save(mydata_postQC_MMR, file = file.path(resultDir, "mydata_postQC_MMR.RData"))

# View initial PCA plot
DimPlot(object = mydata_postQC_MMR, reduction = "pca")

#### Seurat Multimodal Reference Mapping ####

## Load Seurat reference data set and previously saved data
pbmc <- LoadH5Seurat("../RAW_DATA/pbmc_multimodal.h5seurat")
load( file = file.path(resultDir, "mydata_postQC_MMR.RData"))

# plot reference to confirm it downloaded correctly 
DimPlot(object = pbmc, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# Find Anchors between the two Seurat datasets
# this finds anchors between the datasets using the SPCA dimensionality reduction, which is based entirely on the RNAseq data
transferAnchorsData10xMultimodalReference <-
  FindTransferAnchors(
    reference = pbmc,
    query = mydata_postQC_MMR,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50)
# save version of transferAnchors for easy recovery
#Normalizing query using reference SCT model
#Projecting cell embeddings
#Finding neighborhoods
#Finding anchors
#Found 4441 anchors
#save(transferAnchorsData10xMultimodalReference, file = file.path(resultDir, "transferAnchorsData10xMultimodalReference.Rdata"))
load(file.path(resultDir, "transferAnchorsData10xMultimodalReference.Rdata"))

## Map Query
# this maps our 10x dataset into the space of the reference dataset, and pulls in the inferred cell types and inferred ADT values from the reference dataset, 
# based on similarity in expression of genes that are found in both
data10x <- MapQuery(
  anchorset = transferAnchorsData10xMultimodalReference,
  query = mydata_postQC_MMR,
  reference = pbmc,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

data10x@meta.data$predicted.celltype.l1 <-
  factor(data10x@meta.data$predicted.celltype.l1,
         levels = c("B", "CD4 T", "CD8 T", "other T", "NK", "Mono", "DC", "other"))
data10x@meta.data$predicted.celltype.l2 <-
  factor(data10x@meta.data$predicted.celltype.l2,
         levels = str_sort(unique(data10x@meta.data$predicted.celltype.l2), numeric = TRUE))
data10x@meta.data$predicted.celltype.l3 <-
  factor(data10x@meta.data$predicted.celltype.l3,
         levels = str_sort(unique(data10x@meta.data$predicted.celltype.l3), numeric = TRUE))

# Explore reference mapping
seuratMapped_predicted_celltype.l1 <- DimPlot(
  object = data10x, 
  reduction = "ref.umap", 
  group.by = "predicted.celltype.l1",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(seuratMapped_predicted_celltype.l1, file = file.path(plotDir, "seuratMapped_predicted_celltype.l1_UMAP.pdf"), device = "tiff")

seuratMapped_predicted_celltype.l2 <- DimPlot(
  object = data10x, 
  reduction = "ref.umap", 
  group.by = "predicted.celltype.l2",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(seuratMapped_predicted_celltype.l2, file = file.path(plotDir, "seuratMapped_predicted_celltype.l2_UMAP.pdf"), device = "tiff")

# Explore reference mapping
seuratMapped_predicted_celltype.l3 <- DimPlot(
  object = data10x, 
  reduction = "ref.umap", 
  group.by = "predicted.celltype.l3",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(seuratMapped_predicted_celltype.l3, file = file.path(plotDir, "seuratMapped_predicted_celltype.l3_UMAP.pdf"), device = "tiff")

## Plot score for assignments
data10x_CD8_TEM_CD8_TCM <- FeaturePlot(data10x, features = c("CD8 TEM", "CD8 TCM"),  
                                       reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
ggsave(data10x_CD8_TEM_CD8_TCM, file = file.path(plotDir, "data10x_CD8_TEM_CD8_TCM.pdf"), width = 10, height = 5)

data10x_CD4_TEM_CD4_TCM <- FeaturePlot(data10x, features = c("CD4 TEM", "CD4 TCM"),  
                                       reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
ggsave(data10x_CD4_TEM_CD4_TCM, file = file.path(plotDir, "data10x_CD4_TEM_CD4_TCM.pdf"), width = 10, height = 5)

# Export the metadata
Seurat_reference_mapped_metadata <- data10x@meta.data
save(Seurat_reference_mapped_metadata , file = file.path(resultDir, "Seurat_reference_mapped_metadata.RData"))

#### Computing de novo UMAP visualization ####

# This step was recommended by the Seurat package developers
# "We emphasize that if users are attempting to map datasets where the underlying samples are not PBMC, or contain cell types that are not present in the reference, 
# "computing a ‘de novo’ visualization is an important step in interpreting their dataset.
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html#example-1-mapping-human-peripheral-blood-cells-1

#merge reference and query
pbmc$id <- 'reference'
data10x$id <- 'query'
refquery <- merge(pbmc, data10x)
refquery[["spca"]] <- merge(pbmc[["spca"]], data10x[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)

save(refquery, file = file.path(resultDir, "refquery.RData"))

# Plot new de-novo population UMAP
query_vs_reference_UMAP <- DimPlot(refquery, group.by = 'id', shuffle = TRUE)
ggsave(query_vs_reference_UMAP, file = file.path(plotDir, "query_vs_reference_UMAP.pdf"))

# plot predicted cell type de novo
predicted_cell_type_denovo <- DimPlot(
  object = refquery, 
  group.by = "predicted.celltype.l2",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(predicted_cell_type_denovo , file = file.path(plotDir, "predicted_cell_type_denovo.pdf"))

#### Map Seurat cell definitions onto Monocle data ####

# Load Monocle clustering 
load("../EW_T1DAL_Results/P362-1 T1DAL cds object - postQC no MAIT cells.Rdata")
cds_seurat_reference <- cds_no_MAIT

## Put new Seurat metadata in the same order as the monocle metadata
head(colData(cds_seurat_reference))
head(Seurat_reference_mapped_metadata)
Seurat_reference_mapped_metadata %>% dplyr::count(id) # all are query none are reference
#id     n
#1 query 45115

all(row.names(colData(cds_seurat_reference)) %in% row.names(Seurat_reference_mapped_metadata)) # TRUE
Seurat_reference_mapped_metadata_ordered <- Seurat_reference_mapped_metadata[row.names(colData(cds_seurat_reference)),]
head(Seurat_reference_mapped_metadata_ordered)
colnames(Seurat_reference_mapped_metadata_ordered)

# Join metadata columns for the new cell types and their predictions 
colData(cds_seurat_reference)$predicted.celltype.l1.score <- Seurat_reference_mapped_metadata_ordered$predicted.celltype.l1.score
colData(cds_seurat_reference)$predicted.celltype.l1 <- Seurat_reference_mapped_metadata_ordered$predicted.celltype.l1
colData(cds_seurat_reference)$predicted.celltype.l2.score <- Seurat_reference_mapped_metadata_ordered$predicted.celltype.l2.score
colData(cds_seurat_reference)$predicted.celltype.l2 <- Seurat_reference_mapped_metadata_ordered$predicted.celltype.l2
colData(cds_seurat_reference)$predicted.celltype.l3.score <- Seurat_reference_mapped_metadata_ordered$predicted.celltype.l3.score
colData(cds_seurat_reference)$predicted.celltype.l3 <- Seurat_reference_mapped_metadata_ordered$predicted.celltype.l3

colData(cds_seurat_reference)

# set colors for level1
levels(as.factor(colData(cds_seurat_reference)$predicted.celltype.l1))
L1_colors <- c("#6c7ed7","#bc9b3c","#5d3686","#6ca24d","#c26abb","#46c19a","#b94a73","#b8533c")

# Plot cells in monocle with each level of annotation 
plot_cells(cds_seurat_reference, color_cells_by="predicted.celltype.l1",show_trajectory_graph = F, 
           cell_size=0.5, label_cell_groups = FALSE) 
labelled_l1 <- plot_cells(cds_seurat_reference, color_cells_by="predicted.celltype.l1",show_trajectory_graph = F, 
                          cell_size=0.5,  group_label_size = 4) 
ggsave(labelled_l1, file = file.path(plotDir, "labelled_l1.pdf"))


labelled_l2 <-plot_cells(cds_seurat_reference, color_cells_by="predicted.celltype.l2",show_trajectory_graph = F, cell_size=0.5,
                         label_cell_groups = FALSE,  group_label_size = 4) 
ggsave(labelled_l2, file = file.path(plotDir, "labelled_l2.pdf"), width = 6, height = 4)

labelled_l2_plot_label <- plot_cells(cds_seurat_reference, color_cells_by="predicted.celltype.l2",show_trajectory_graph = F, cell_size=0.5,
                                     group_label_size = 4) 
ggsave(labelled_l2_plot_label, file = file.path(plotDir, "labelled_l2_plot_label.pdf"))


plot_cells(cds_seurat_reference, color_cells_by="predicted.celltype.l3",show_trajectory_graph = F, cell_size=1) 

# Plot CD4 and CD8 expression
CD4_CD8_expression <- plot_cells(cds_seurat_reference, genes = c("CD8A", "CD4"),show_trajectory_graph = F)
ggsave(CD4_CD8_expression, file = file.path(plotDir, "CD4_CD8_expression.pdf"), width = 8, height = 4)


# Plot CD4 score
Seurat_reference_mapped_metadata_ordered <- Seurat_reference_mapped_metadata_ordered %>% mutate(CD4_score = case_when(
  predicted.celltype.l2 %in% c( "CD4 CTL","CD4 Naive","CD4 TCM","CD4 TEM" ) ~ predicted.celltype.l2.score,
  !(predicted.celltype.l2 %in% c( "CD4 CTL","CD4 Naive","CD4 TCM","CD4 TEM" )) ~ 0
))
Seurat_reference_mapped_metadata_ordered$CD4_score

# add CD4 score to df
colData(cds_seurat_reference)$CD4_score <- Seurat_reference_mapped_metadata_ordered$CD4_score

# plot CD4 score
CD4_score <- plot_cells(cds_seurat_reference, color_cells_by="CD4_score",show_trajectory_graph = F, cell_size=0.5) 
ggsave(CD4_score , file = file.path(plotDir, "CD4_score.pdf"))

# save the cds object with the cell types mapped onto it
save(cds_seurat_reference, file = file.path(resultDir, "cds_seurat_reference.Rdata"))

#### Map T1DAL data onto to reference CD8 T cells ####

## Load Seurat reference data set and previously saved data
pbmc <- LoadH5Seurat("../RAW_DATA/pbmc_multimodal.h5seurat")
load( file = file.path(resultDir, "mydata_postQC_MMR.RData"))

# plot reference to confirm it downloaded correctly 
DimPlot(object = pbmc, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# Find CD8 categories in the data
levels(as.factor(pbmc@meta.data$celltype.l1)) #  "CD8 T" 

# Subset to keep only CD8 T cells
pbmc_CD8 <- subset(pbmc, subset = celltype.l1 == "CD8 T")

# this finds anchors between the datasets using the SPCA dimensionality reduction, which is based entirely on the RNAseq data
transferAnchorsData10xMultimodalReference_CD8 <-
  FindTransferAnchors(
    reference = pbmc_CD8,
    query = mydata_postQC_MMR,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50)
# save version of transferAnchors for easy recovery
#Normalizing query using reference SCT model
#Projecting cell embeddings
#Finding neighborhoods
#Finding anchors
#Found 3018 anchors
#save(transferAnchorsData10xMultimodalReference_CD8, file = file.path(resultDir, "transferAnchorsData10xMultimodalReference_CD8.Rdata"))
load(file.path(resultDir, "transferAnchorsData10xMultimodalReference_CD8.Rdata"))

## Map Query
# this maps our 10x dataset into the space of the reference dataset, and pulls in the inferred cell types and inferred ADT values from the reference dataset, 
# based on similarity in expression of genes that are found in both
data10x_CD8 <- MapQuery(
  anchorset = transferAnchorsData10xMultimodalReference_CD8,
  query = mydata_postQC_MMR,
  reference = pbmc_CD8,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

data10x_CD8@meta.data$predicted.celltype.l1 <-
  factor(data10x_CD8@meta.data$predicted.celltype.l1,
         levels = c("B", "CD4 T", "CD8 T", "other T", "NK", "Mono", "DC", "other"))
data10x_CD8@meta.data$predicted.celltype.l2 <-
  factor(data10x_CD8@meta.data$predicted.celltype.l2,
         levels = str_sort(unique(data10x_CD8@meta.data$predicted.celltype.l2), numeric = TRUE))
data10x_CD8@meta.data$predicted.celltype.l3 <-
  factor(data10x_CD8@meta.data$predicted.celltype.l3,
         levels = str_sort(unique(data10x_CD8@meta.data$predicted.celltype.l3), numeric = TRUE))

# Explore reference mapping
seuratMapped_predicted_celltype.l1_CD8 <- DimPlot(
  object = data10x_CD8, 
  reduction = "ref.umap", 
  group.by = "predicted.celltype.l1",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(seuratMapped_predicted_celltype.l1_CD8, file = file.path(plotDir, "seuratMapped_predicted_celltype.l1_UMAP_CD8.pdf"), device = "tiff")

seuratMapped_predicted_celltype.l2_CD8 <- DimPlot(
  object = data10x_CD8, 
  reduction = "ref.umap", 
  group.by = "predicted.celltype.l2",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(seuratMapped_predicted_celltype.l2_CD8, file = file.path(plotDir, "seuratMapped_predicted_celltype.l2_UMAP_CD8.pdf"), device = "tiff")

# Explore reference mapping
seuratMapped_predicted_celltype.l3_CD8 <- DimPlot(
  object = data10x_CD8, 
  reduction = "ref.umap", 
  group.by = "predicted.celltype.l3",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(seuratMapped_predicted_celltype.l3_CD8, file = file.path(plotDir, "seuratMapped_predicted_celltype.l3_UMAP_CD8.pdf"), device = "tiff")

# Export the metadata
Seurat_reference_mapped_metadata_CD8 <- data10x_CD8@meta.data
save(Seurat_reference_mapped_metadata_CD8, file = file.path(resultDir, "Seurat_reference_mapped_metadata_CD8.RData"))

#### Computing de novo UMAP visualization with CD8 mapped data ####

# This step was recommended by the Seurat package developers
# "We emphasize that if users are attempting to map datasets where the underlying samples are not PBMC, or contain cell types that are not present in the reference, 
# "computing a ‘de novo’ visualization is an important step in interpreting their dataset.
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html#example-1-mapping-human-peripheral-blood-cells-1

#merge reference and query
pbmc_CD8$id <- 'reference'
data10x_CD8$id <- 'query'
refquery_CD8 <- merge(pbmc_CD8, data10x_CD8)
refquery_CD8[["spca"]] <- merge(pbmc_CD8[["spca"]], data10x_CD8[["ref.spca"]])
refquery_CD8 <- RunUMAP(refquery_CD8, reduction = 'spca', dims = 1:50)

save(refquery_CD8, file = file.path(resultDir, "refquery_CD8.RData"))

# Plot new de-novo population UMAP
query_vs_reference_UMAP_CD8 <- DimPlot(refquery_CD8, group.by = 'id', shuffle = TRUE)
ggsave(query_vs_reference_UMAP_CD8, file = file.path(plotDir, "query_vs_reference_UMAP_CD8.pdf"))

# plot predicted cell type de novo
predicted_cell_type_denovo_CD8 <- DimPlot(
  object = refquery_CD8, 
  group.by = "predicted.celltype.l2",
  label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
ggsave(predicted_cell_type_denovo_CD8 , file = file.path(plotDir, "predicted_cell_type_denovo_CD8.pdf"))

#### Map CD8 Seurat cell definitions onto Monocle data ####

# Load Monocle clustering 
load("../EW_T1DAL_Results/P362-1 T1DAL cds object - postQC no MAIT cells.Rdata")
cds_seurat_reference_CD8 <- cds_no_MAIT
# load final CD8 reference data for use later in this section
load(file = file.path(resultDir, "cds_seurat_reference_CD8.Rdata"))

## Put new Seurat metadata in the same order as the monocle metadata
head(colData(cds_seurat_reference_CD8))
head(Seurat_reference_mapped_metadata_CD8)

all(row.names(colData(cds_seurat_reference_CD8)) %in% row.names(Seurat_reference_mapped_metadata_CD8)) # TRUE
Seurat_reference_mapped_metadata_CD8_ordered <- Seurat_reference_mapped_metadata_CD8[row.names(colData(cds_seurat_reference_CD8)),]
head(Seurat_reference_mapped_metadata_CD8_ordered)
colnames(Seurat_reference_mapped_metadata_CD8_ordered)

# Join metadata columns for the new cell types and their predictions 
colData(cds_seurat_reference_CD8)$predicted.celltype.l1.score <- Seurat_reference_mapped_metadata_CD8_ordered$predicted.celltype.l1.score
colData(cds_seurat_reference_CD8)$predicted.celltype.l1 <-       Seurat_reference_mapped_metadata_CD8_ordered$predicted.celltype.l1
colData(cds_seurat_reference_CD8)$predicted.celltype.l2.score <- Seurat_reference_mapped_metadata_CD8_ordered$predicted.celltype.l2.score
colData(cds_seurat_reference_CD8)$predicted.celltype.l2 <-       Seurat_reference_mapped_metadata_CD8_ordered$predicted.celltype.l2
colData(cds_seurat_reference_CD8)$predicted.celltype.l3.score <- Seurat_reference_mapped_metadata_CD8_ordered$predicted.celltype.l3.score
colData(cds_seurat_reference_CD8)$predicted.celltype.l3 <-       Seurat_reference_mapped_metadata_CD8_ordered$predicted.celltype.l3

colData(cds_seurat_reference_CD8)

# set colors for level1
levels(as.factor(colData(cds_seurat_reference_CD8)$predicted.celltype.l1))
L1_colors <- c("#6c7ed7","#bc9b3c","#5d3686","#6ca24d","#c26abb","#46c19a","#b94a73","#b8533c")

cds_seurat_reference_CD8_figure <- cds_seurat_reference_CD8
colnames(colData(cds_seurat_reference_CD8_figure))[18] <- "Seurat Reference\nMapping L2"
# Plot cells in monocle with each level of annotation 
labelled_l2_CD8 <- plot_cells(cds_seurat_reference_CD8_figure[,colnames(cds_seurat_reference_CD8_figure) %in% row.names(colData(cds_seurat_reference_CD8_figure) %>% as.data.frame() %>% filter(Cluster.Name !="9")) ],
                              color_cells_by="Seurat Reference\nMapping L2",show_trajectory_graph = F, 
                              cell_size=1,
                              label_cell_groups = FALSE,  ) + 
  scale_colour_manual(values=c("#b54f90",
                               "#ad993c",
                               "#7066bc",
                               "#56ae6c"
  )) +
  labs(fill="Seurat Predicted Cell Types L2") + 
  theme(text = element_text(size = 20))

#### FIGURE S5C ####
ggsave(labelled_l2_CD8, file = file.path(plotDir, "labelled_l2_CD8.pdf"), width= 8, height = 6)

labelled_l2_plot_label_CD8 <- plot_cells(cds_seurat_reference_CD8, color_cells_by="predicted.celltype.l2",show_trajectory_graph = F, cell_size=0.5,
                                         group_label_size = 4) 
ggsave(labelled_l2_plot_label_CD8, file = file.path(plotDir, "labelled_l2_plot_label_CD8.pdf"))


labelled_l3_plot_label_CD8 <- plot_cells(cds_seurat_reference_CD8, color_cells_by="predicted.celltype.l3",
                                         show_trajectory_graph = F, cell_size=0.5,  group_label_size = 4) 
ggsave(labelled_l3_plot_label_CD8, file = file.path(plotDir, "labelled_l3_plot_label_CD8.pdf"))

# plot with individual labels for each level
l3_plot_label_CD8 <- plot_cells(cds_seurat_reference_CD8, color_cells_by="predicted.celltype.l3",
                                show_trajectory_graph = F, cell_size=0.5, label_cell_groups = FALSE) 
ggsave(l3_plot_label_CD8 , file = file.path(plotDir, "l3_plot_label_CD8.pdf"))
naive_expression <- plot_cells(cds_seurat_reference_CD8, genes = c("CCR7","SELL","IL7R"),
                               show_trajectory_graph = F, cell_size=0.5, label_cell_groups = FALSE)
ggsave(naive_expression  , file = file.path(plotDir, "naive_expression.pdf"), height = 5, width = 12)


# save the cds object with the cell types mapped onto it
save(cds_seurat_reference_CD8, file = file.path(resultDir, "cds_seurat_reference_CD8.Rdata"))

