#### P452-3 T1DAL ATAC-seq Comparison  ####

## This script is meant to compare the DiffBind peaks between our data and other published datasets


#### LOAD LIBRARIES ####

library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # human transcript
library(EnsDb.Hsapiens.v86)
library(stringr) # string processing
library(readxl)
library(soGGi) # for calcuation of TSSe scores
library(RColorBrewer)
library(igraph)
library(ggplot2) # plotting
library(ATACseqQC)
library(GenomicAlignments)
library(DESeq2)
library(egg)
library(ComplexHeatmap)
library(DiffBind)
library(ChIPpeakAnno)
library(reactome.db) # to get enrichment paths
library(gridExtra)
library(ComplexHeatmap)
library(Repitools)
library(apird)
library(gridGraphics)
library(grid)
library(ggfortify)
library(gplots)
library(devtools)
library(ggrepel)
#install_github("js229/Vennerable")
library(Vennerable)
library(eulerr)
library(KEGGREST)
library(gridExtra)
library(karyoploteR)
library(wiggleplotr)
library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)
library(ensembldb)
library(readxl)
library(IRanges)
library(ggpubr)

# cell pop color codes: # DN non naive= 7497d9ff, Pd1 = ba4d4cc7, CD57= 93a24eff 

setwd("/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup")
set.seed(42) # same seed as QC script 

#### LOAD ALL SAVED DATA AND SET PATHS ####

# Set Paths 
baseDir <- "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup"
plotDir <- file.path(baseDir,'FIGURES')
resultDir <- file.path(baseDir,'SAVED_DATA')
annotationDir <- file.path(baseDir, "RAW_DATA")

# load anno manauly
project <- "P452-3"
getGcqProjectInfo(project)
libids <- getProjectLibs(project, searchType = "regex")
anno <- getAnno(libids)

#### LOAD SIGNIFICANT ATAC PEAKS FROM CD57 POS VS DN AND CD57 MINUS VS DN ####

# note this loads all up and down peaks that were significant at FDR p < 0.05
CD57pos_DN <- read.csv(file = file.path(resultDir, "ATACseqData_P452_3_norm_db_anno_CD57pos_DN.csv"))
CD57minus_DN <- read.csv(file = file.path(resultDir, "ATACseqData_P452_3_norm_db_anno_CD57minus_DN.csv"))
CD57_vs_PD1 <- read.csv(file = file.path(resultDir, "ATACseqData_P452_3_norm_db_anno.csv"))

#### LOAD IN EPIGENETIC DATASETS FOR COMPARISON ####

## LOAD IN GILES ET AL 2022 IMMUNITY EPIGENETIC DATASET
# this data was taken from Giles et al., Supplementary table mmc3 Table S2 from this link: https://europepmc.org/article/pmc/pmc9214622#SM1
# document consists of 9 sheets, read in individually 

# This Giles dataset includes TEMRA cells for comparison, but not NK-like Tex

## Figure Caption: Table S2 (Related to Figure 2 and Figure 5). Differentially accessible peaks between all pairwise
#comparisons for HD CD8 T cell subsets and ACRs used in UMAP construction

multiplesheets <- function(fname) {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # print data frame
  print(data_frame)
}

# specifying the path name
path <- file.path(annotationDir, "Giles_ATAC_NIHMS1781187-supplement-mmc3.xlsx")
giles_Table_S2 <- multiplesheets(path)

## LOAD IN DANIEL ET AL 2022 NATURE IMMUNOLOGY EPIGENETIC DATASET 
# Loading in the Peak Scores section of Supplementary Table 2 from their paper here: https://www.nature.com/articles/s41590-022-01337-5#Sec33 
# description about this figure from the text: Cell-type-specific OCR accessibility was
# correlated with gene expression at marker gene loci that define Tex subsets, including Tcf7, Pdcd1 and Tox (Fig. 1g,h and Supplementary Table 2). Tnaive 

# NOTE: Column A-M refer to gene expression values, Columns O-AB refer to cell specific differential accessibility
Daniel_OCRs <- read_excel(file.path(annotationDir, "Daniel_science_41590_2022_1337_MOESM2_ESM.xlsx"),
                          # susbet to keep Supplementary Table 2 and columns O - AB which has the OCR data
                          sheet =2, range = cell_cols("O:AB"))

#### FIND ALL OVERLAPPING SIGNIFICANT PEAKS IN GILES TEMRA DATA VERSUS CD57 vs PD1 ####

# This time compare only with the TEMRA dataset
giles_Table_S2_EMRA <- giles_Table_S2$EMRA

# Write function to identify any overlaps between list combinations
# each list from here are the results of that cluster vs. other clusters individually, and sometimes multiple datasets in the list have the same significant peak
# based on this, each of my hits may match back to any may match back to multiple original hits 

# get list of TEMRA contrasts
giles_Table_S2_EMRA_list <- unique(giles_Table_S2_EMRA$subset_comparison)

TEMRA_overlaps_giles <- function(subject, query) {
  # subset query to also include only the things that were up in the CD57pos or CD57 minus population
  query_data <- query %>%
    dplyr::filter(FDR <= 0.05) 
  # compare all up and down peaks to the TEMRA set!
  # filter(Fold >= 0) 
  
  # subset giles data
  subject_data <- giles_Table_S2_EMRA %>% dplyr::filter(subset_comparison == subject)
  # filter for sig peaks, those that are up in cell type of interest, then keep only distinct peaks
  subject_data <- subject_data %>% 
    dplyr::filter(padj <= 0.05) %>% 
    #  filter(log2FoldChange >= 0) %>%  # compare all peaks up and down 
    distinct(chr, start, end)
  
  # convert into IRanges objects
  ranges_subject <- split(IRanges(subject_data$start, subject_data$end), subject_data$chr) # note chr column name
  ranges_query <- split(IRanges(query_data$start, query_data$end), query_data$seqnames) # not chr is in seqnames
  
  # subset to keep any overlaps
  overlap <- subsetByOverlaps(ranges_query, ranges_subject )
  
  # export results to dataframe
  overlap_df = as(overlap, "data.frame") 
  overlap_df$dataset <- paste(subject)
  
  # add numbers on for computing jaccard
  overlap_df$giles_length <- nrow(subject_data)
  overlap_df$query_length <- nrow(query_data)
  
  #output
  overlap_df
  
}

CD57_PD1_TEMRA_overlap <- lapply(giles_Table_S2_EMRA_list, TEMRA_overlaps_giles, CD57_vs_PD1) 
names(CD57_PD1_TEMRA_overlap) <- giles_Table_S2_EMRA_list
CD57_PD1_TEMRA_overlap_df <- bind_rows(CD57_PD1_TEMRA_overlap)

# save TEMRA overlap
#save(CD57_PD1_TEMRA_overlap_df, file = file.path(resultDir, "CD57_PD1_TEMRA_overlap_df.Rdata"))

## How many hits were there in each group - this will be the intersection
CD57_PD1_TEMRA_overlap_df_intersect <- CD57_PD1_TEMRA_overlap_df %>% group_by(dataset, giles_length, query_length) %>%
  summarize(intersect = n())

### Calculate Jaccard similarity for each group of overlapping peaks 
# Jaccard similarity is the ratio of the union between two datasets and the intersection of the datasets

# calculate jaccard
CD57_PD1_TEMRA_overlap_df_intersect_jaccard <- CD57_PD1_TEMRA_overlap_df_intersect %>%
  mutate(union = (giles_length + query_length - intersect),
         Jaccard = intersect/union) %>% arrange(desc(Jaccard))
write.csv(CD57_PD1_TEMRA_overlap_df_intersect_jaccard, file.path(resultDir,"CD57_PD1_TEMRA_overlap_df_intersect_jaccard.csv" ), row.names = FALSE)

### Run phyper hypergeometric tests to assess over-representation of CD57+ vs CD57minus list in each TEMRA list
#phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

TEMRA_giles_phyper_results <- data.frame()
for (row in 1:nrow(CD57_PD1_TEMRA_overlap_df_intersect_jaccard)) {
  # group1 = our dataa
  # group2 = giles TEMRA
  total = 18799*2 #which was the number of background genes in the ATAC analysis 
  set <- CD57_PD1_TEMRA_overlap_df_intersect_jaccard[row, ]
  phyper_df <- phyper(set$intersect - 1, set$giles_length, (total - set$giles_length), set$query_length, lower.tail = FALSE )
  TEMRA_giles_phyper_results <- rbind(phyper_df, TEMRA_giles_phyper_results)
  
}

#### WRITE FUNCTION TO ANNOTATE PEAKS ####

# Annotate peaks using function 
# the annotatePeakInBatch function obtains the distance to the nearest TSS, 
# miRNA, exon et al for a list of peak locations leveraging IRanges and biomaRt package

peak_annotate <- function(data, list) {
  
  anno <- annotatePeakInBatch(list, 
                              AnnotationData = data, 
                              output = "overlapping",
                              bindingRegion = c(-5000, 5000),
                              select = "first", 
                              multiple =  FALSE)
  anno <- data.frame(anno) 
  anno
}

## Get anno data

## Format Annotation data using formatted annotation from Matt L/Hannah
annoData <- readRDS(file.path(resultDir,"annoData_human_v86.rds")) 
wantedLevels = paste("chr", c(1:22, "X", "Y"), sep = '')
noChrLevels = c(1:22, "X", "Y")

idxVec = unlist(lapply(annoData@seqnames,FUN = function(y) any(str_detect(y, paste("chr", c(1:22, "X", "Y"), sep = '')))))
annoData = annoData[idxVec]
annoData = keepSeqlevels(annoData, wantedLevels)
seqlevelsStyle(annoData) <- "UCSC" 


#### ANNOTATE ALL PEAKS IN TEMRA IN GILES DATASET AND ASSESS OVERLAP ####

# This time compare only with the TEMRA dataset
giles_Table_S2_EMRA <- giles_Table_S2$EMRA
nrow(giles_Table_S2_EMRA) # 151159

## Annotate ALL GENES IN TEMRA and recombine afterward
# get list of TEMRA contrasts
giles_Table_S2_EMRA_list <- unique(giles_Table_S2_EMRA$subset_comparison)

TEMRA_anno <- function(subject) {
  
  subset_data <- giles_Table_S2_EMRA %>% dplyr::filter(subset_comparison == subject) %>% 
    #filter(log2FoldChange >=0) %>%
    dplyr::filter(padj <= 0.05) 
  
  anno <- makeGRangesFromDataFrame(subset_data[,c("chr","start","end","padj","log2FoldChange","subset_comparison")])
  seqlevelsStyle(anno   ) <-"UCSC"
  anno <- peak_annotate(anno, data = annoData)
  anno$subset <- subject
  anno
}

TEMRA_overlaps_giles_all <- lapply(giles_Table_S2_EMRA_list, TEMRA_anno)
names(TEMRA_overlaps_giles_all) <- giles_Table_S2_EMRA_list
TEMRA_overlaps_giles_all_df <- bind_rows(TEMRA_overlaps_giles_all)

#save(TEMRA_overlaps_giles_all_df, file=file.path(resultDir, "TEMRA_overlaps_giles_all_df.Rdata"))

#### CALCULATE OVERLAP WITH ALL GENES TEMRA AND ALL IN CD57 ####

# load
load(file=file.path(resultDir, "TEMRA_overlaps_giles_all_df.Rdata"))

# This time compare only with the TEMRA dataset
giles_Table_S2_EMRA <- giles_Table_S2$EMRA

# Write function to identify any overlaps between list combinations
# each list from here are the results of that cluster vs. other clusters individually, and sometimes multiple datasets in the list have the same significant peak
# based on this, each of my hits may match back to any may match back to multiple original hits 

# get list of TEMRA contrasts
giles_Table_S2_EMRA_list <- unique(giles_Table_S2_EMRA$subset_comparison)

TEMRA_overlaps_giles <- function(subject, query) {
  # subset query to also include only the things that were up in the CD57pos 
  query_data <- query %>%
    dplyr::filter(FDR <= 0.05) 
  
  # subset giles data - USE DATSET ANNOTATED TO GENES
  subject_data <- TEMRA_overlaps_giles_all_df %>% dplyr::filter(subset == subject) %>%
    #filter(padj <= 0.05)
    distinct(seqnames, start, end)
  
  # convert into IRanges objects
  ranges_subject <- split(IRanges(subject_data$start, subject_data$end), subject_data$seqnames) # note chr column name
  ranges_query <- split(IRanges(query_data$start, query_data$end), query_data$seqnames) # not chr is in seqnames
  
  # subset to keep any overlaps
  overlap <- subsetByOverlaps(ranges_query, ranges_subject )
  
  # export results to dataframe
  overlap_df = as(overlap, "data.frame") 
  overlap_df$dataset <- paste(subject)
  
  # add numbers on for computing jaccard
  overlap_df$giles_length <- nrow(subject_data)
  overlap_df$query_length <- nrow(query_data)
  
  #output
  overlap_df
  
}

CD57_PD1_TEMRA_overlap <- lapply(giles_Table_S2_EMRA_list, TEMRA_overlaps_giles, CD57_vs_PD1) 
names(CD57_PD1_TEMRA_overlap) <- giles_Table_S2_EMRA_list
CD57_PD1_TEMRA_overlap_df <- bind_rows(CD57_PD1_TEMRA_overlap)

## How many hits were there in each group - this will be the intersection
CD57_PD1_TEMRA_overlap_df_intersect <- CD57_PD1_TEMRA_overlap_df %>% group_by(dataset, giles_length, query_length) %>%
  summarize(intersect = n())

### Calculate Jaccard similarity for each group of overlapping peaks 
# Jaccard similarity is the ratio of the union between two datasets and the intersection of the datasets

# calculate jaccard
CD57_PD1_TEMRA_overlap_df_intersect_jaccard <- CD57_PD1_TEMRA_overlap_df_intersect %>%
  mutate(union = (giles_length + query_length - intersect),
         Jaccard = intersect/union) %>% arrange(desc(Jaccard))
write.csv(CD57_PD1_TEMRA_overlap_df_intersect_jaccard, file.path(resultDir,"CD57_PD1_TEMRA_overlap_df_intersect_jaccard.csv" ), row.names = FALSE)

### Run phyper hypergeometric tests to assess over-representation of CD57+list in each TEMRA up list
#phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

TEMRA_giles_phyper_results <- data.frame()
for (row in 1:nrow(CD57_PD1_TEMRA_overlap_df_intersect_jaccard)) {
  # group1 = our data
  # group2 = giles TEMRA
  total = 18799 #which was the number of background genes in the ATAC analysis 
  set <- CD57_PD1_TEMRA_overlap_df_intersect_jaccard[row, ]
  print(set)
  phyper_df <- phyper(set$intersect - 1, set$giles_length, (total - set$giles_length), set$query_length, lower.tail = FALSE )
  phyper_df <- cbind(as.data.frame(phyper_df), set)
  print(phyper_df )
  TEMRA_giles_phyper_results <- rbind(phyper_df, TEMRA_giles_phyper_results)
  
}

TEMRA_giles_phyper_results # very little overlap
write.csv(TEMRA_giles_phyper_results , file.path(resultDir, "TEMRA_giles_phyper_results.csv"), row.names = FALSE)

#### ANNOTATE PEAKS INCREASED IN TEMRA IN GILES DATASET AND ASSESS OVERLAP ####

# This time compare only with the TEMRA dataset
giles_Table_S2_EMRA <- giles_Table_S2$EMRA
nrow(giles_Table_S2_EMRA) # 151159

## Annotate genes that are up in TEMRA and recombine afterward
# get list of TEMRA contrasts
giles_Table_S2_EMRA_list <- unique(giles_Table_S2_EMRA$subset_comparison)

TEMRA_anno_up <- function(subject) {
  
  subset_data <- giles_Table_S2_EMRA %>% dplyr::filter(subset_comparison == subject) %>% 
    dplyr::filter(log2FoldChange >=0) %>%
    dplyr::filter(padj <= 0.05) 
  
  anno <- makeGRangesFromDataFrame(subset_data[,c("chr","start","end","padj","log2FoldChange","subset_comparison")])
  seqlevelsStyle(anno   ) <-"UCSC"
  anno <- peak_annotate(anno, data = annoData)
  anno$subset <- subject
  anno
}

TEMRA_overlaps_giles_all_up <- lapply(giles_Table_S2_EMRA_list, TEMRA_anno_up)
names(TEMRA_overlaps_giles_all_up) <- giles_Table_S2_EMRA_list
TEMRA_overlaps_giles_all_up_df <- bind_rows(TEMRA_overlaps_giles_all_up)

#save(TEMRA_overlaps_giles_all_up_df , file = file.path(resultDir, "TEMRA_overlaps_giles_all_up_df.Rdata"))

#### FIGURE S16E: CALCULATE OVERLAP WITH JUST GENES UP IN TEMRA AND UP IN CD57 ####

# load data
load(file = file.path(resultDir, "TEMRA_overlaps_giles_all_up_df.Rdata"))

# Find overlap between TEMRA up and CD57up and then calculate jacard

# Write function to identify any overlaps between list combinations

# get list of TEMRA contrasts
giles_Table_S2_EMRA_list <- unique(giles_Table_S2_EMRA$subset_comparison)

TEMRA_overlaps_giles_up <- function(subject, query) {
  # subset query to also include only the things that were up in the CD57pos 
  query_data <- query %>%
    dplyr::filter(FDR <= 0.05) %>%
    # compare all up peaks to the TEMRA set!
    dplyr::filter(Fold >= 0) 
  
  # subset giles data - USE DATSET ANNOTATED TO GENES
  subject_data <- TEMRA_overlaps_giles_all_up_df %>% dplyr::filter(subset == subject) %>%
    #filter(padj <= 0.05) %>% # ALREADY DID THIS ABOVE when annotating
    # filter(log2FoldChange >= 0) %>%  # ALREADY DID THIS ABOVE when annotating
    distinct(seqnames, start, end)
  
  # convert into IRanges objects
  ranges_subject <- split(IRanges(subject_data$start, subject_data$end), subject_data$seqnames) # note chr column name
  ranges_query <- split(IRanges(query_data$start, query_data$end), query_data$seqnames) # not chr is in seqnames
  
  # subset to keep any overlaps
  overlap <- subsetByOverlaps(ranges_query, ranges_subject )
  
  # export results to dataframe
  overlap_df = as(overlap, "data.frame") 
  overlap_df$dataset <- paste(subject)
  
  # add numbers on for computing jaccard
  overlap_df$giles_length <- nrow(subject_data)
  overlap_df$query_length <- nrow(query_data)
  
  #output
  overlap_df
  
}

CD57_PD1_TEMRA_overlap_up <- lapply(giles_Table_S2_EMRA_list, TEMRA_overlaps_giles_up, CD57_vs_PD1) 
names(CD57_PD1_TEMRA_overlap_up) <- giles_Table_S2_EMRA_list
CD57_PD1_TEMRA_overlap_df_up <- bind_rows(CD57_PD1_TEMRA_overlap_up)

## How many hits were there in each group - this will be the intersection
CD57_PD1_TEMRA_overlap_df_up_intersect <- CD57_PD1_TEMRA_overlap_df_up %>% group_by(dataset, giles_length, query_length) %>%
  summarize(intersect = n())

### Calculate Jaccard similarity for each group of overlapping peaks 
# Jaccard similarity is the ratio of the union between two datasets and the intersection of the datasets

# calculate jaccard
CD57_PD1_TEMRA_overlap_df_up_intersect_jaccard <- CD57_PD1_TEMRA_overlap_df_up_intersect %>%
  mutate(union = (giles_length + query_length - intersect),
         Jaccard = intersect/union) %>% arrange(desc(Jaccard))
write.csv(CD57_PD1_TEMRA_overlap_df_up_intersect_jaccard, file.path(resultDir,"CD57_PD1_TEMRA_overlap_df_up_intersect_jaccard.csv" ), row.names = FALSE)

### Run phyper hypergeometric tests to assess over-representation of CD57+list in each TEMRA up list
#phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

TEMRA_giles_phyper_results_up <- data.frame()
for (row in 1:nrow(CD57_PD1_TEMRA_overlap_df_up_intersect_jaccard)) {
  # group1 = our data
  # group2 = giles TEMRA
  total = 18799 #which was the number of background genes in the ATAC analysis 
  set <- CD57_PD1_TEMRA_overlap_df_up_intersect_jaccard[row, ]
  print(set)
  phyper_df <- phyper(set$intersect - 1, set$giles_length, (total - set$giles_length), set$query_length, lower.tail = FALSE )
  phyper_df <- cbind(as.data.frame(phyper_df), set)
  print(phyper_df )
  TEMRA_giles_phyper_results_up <- rbind(phyper_df, TEMRA_giles_phyper_results_up)
  
}

TEMRA_giles_phyper_results_up 

# save 
write.csv(TEMRA_giles_phyper_results_up , file.path(resultDir, "TEMRA_giles_phyper_results_up.csv"), row.names = FALSE)
#save(TEMRA_giles_phyper_results_up , file = file.path(resultDir, "TEMRA_giles_phyper_results_up.Rdata"))

### Plot results as a figure rather than a table
load(file = file.path(resultDir, "TEMRA_giles_phyper_results_up.Rdata"))
TEMRA_giles_phyper_results_up 
TEMRA_giles_phyper_results_up_plot <- ggplot(TEMRA_giles_phyper_results_up, aes(x = Jaccard, y = -log10(phyper_df), label = dataset)) + 
  geom_point() +
  ggrepel::geom_text_repel(data = subset(TEMRA_giles_phyper_results_up, dataset != "temra_vs_PD1+CD39+"),
                           color = "black", size = 4) +
  ggrepel::geom_text_repel(data = subset(TEMRA_giles_phyper_results_up, dataset == "temra_vs_PD1+CD39+"),
                           color = "red", fontface = "bold", size = 4) +
  theme_minimal() + 
  labs(x = "Jaccard Similarity Coefficient", y = "-log10 Hypergeometric P-value") + 
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red")

## FIGURE S16E
ggsave(TEMRA_giles_phyper_results_up_plot, file = file.path(plotDir, "TEMRA_giles_phyper_results_up_plot.pdf"), width = 4.5, height = 4.5)

#### CREATE CD57 ATAC SIGNATURE ####

# GOAL: create a signature that I can aggregate and plot with the percent of EMRA cells for each person for both
# CD57 population and the PD-1 population with their percent of EMRA cells 
# going to use the mean log reads for all significant CD57 genes in each group and then sum these values

# Load normalized data and extract peakset
load(file = file.path(resultDir, "ATACseqData_P452_3_norm_final.RData"))
ATACseqData_P452_3_norm 

# count 
ATACseqData_P452_3_norm_counts <- dba.count(ATACseqData_P452_3_norm,  peaks=NULL, score=DBA_SCORE_NORMALIZED) 

# retrieve peaks
ATACseqData_P452_3_norm_counts_table <- dba.peakset(ATACseqData_P452_3_norm_counts, bRetrieve = TRUE)

# annotate peaks
seqlevelsStyle(ATACseqData_P452_3_norm_counts_table ) <- "UCSC"
ATACseqData_P452_3_norm_counts_table_anno <- peak_annotate(ATACseqData_P452_3_norm_counts_table, data = annoData)

# make data long
ATACseqData_P452_3_norm_counts_table_anno_long <- ATACseqData_P452_3_norm_counts_table_anno %>% 
  pivot_longer(names_to = "SampleID", values_to = "Norm_Counts" ,  cols = c(6:29) ) 
ATACseqData_P452_3_norm_counts_table_anno_long$SampleID <- as.numeric(str_remove(ATACseqData_P452_3_norm_counts_table_anno_long$SampleID, "X"))
class(ATACseqData_P452_3_norm_counts_table_anno_long$SampleID)
# join sample metadata
ATACseqData_P452_3_norm_counts_table_anno_long_meta <- left_join(ATACseqData_P452_3_norm_counts_table_anno_long, as.data.frame(ATACseqData_P452_3_norm$samples))
ATACseqData_P452_3_norm_counts_table_anno_long_meta$seqnames <- str_remove(ATACseqData_P452_3_norm_counts_table_anno_long_meta$seqnames, "chr")
ATACseqData_P452_3_norm_counts_table_anno_long_meta$peak <- as.numeric(ATACseqData_P452_3_norm_counts_table_anno_long_meta$peak)

# get list of genes that are in the CD57 increased gene set
# pivot data longer
CD57_sig <- CD57_vs_PD1 %>% 
  # filter for those peaks that were increased in CD57pos
  dplyr::filter(Fold >=0) 
CD57_sig$seqnames <- str_remove(CD57_sig$seqnames, "chr")

# subset for CD57 signature genes
CD57_sig_peaks <- left_join(CD57_sig, 
                            ATACseqData_P452_3_norm_counts_table_anno_long_meta)
CD57_sig %>% dplyr::count(seqnames,start,end,gene_name) # a few genes have lines for the sample peaks..
CD57_sig_peaks  %>% dplyr::count(seqnames,start,end, gene_name)
CD57_sig %>% group_by(seqnames,start,end) %>% filter(n() >1) %>% View

# get composite CD57 score per library
CD57_sig_peaks_score <- CD57_sig_peaks %>%
  group_by(libId, Condition, Factor) %>%
  summarize(CD57_score = sum(Norm_Counts))

# save
#save(CD57_sig_peaks_score, file = file.path(resultDir, "CD57_sig_peaks_score.Rdata"))

#### FIGURE S16D: PLOT CD57 ATAC SIGNATURE BY EMRA ####

## Load Alice W calculate percent EM per donor and the key

sample_key <- read_xlsx(file.path(annotationDir, "T1DAL_Clusters_from_T1DAL_AbATE_Placebo_8_MCs_for_Erin_2023-10-10.xlsx"), sheet=1) %>%
  dplyr::rename(Donor_ID= ParticipantID)
sample_key$Donor_ID <- as.character(sample_key$Donor_ID)
# Join with anno to get libId and donor ID
anno_colnames <- anno %>% dplyr::rename(timepoint = timePoint, libId = libid, Donor_ID=donorId)
sample_key_anno <- left_join(sample_key, anno_colnames[,c("Donor_ID", "libId", "timepoint","sample_name")]) %>% 
  dplyr::filter(!is.na(libId))
unique(sample_key_anno$Donor_ID) # "10213" "10396" "10295" "10458"

# add population key for joining with EMRA 
sample_key_anno <-  sample_key_anno %>%
  mutate(population = case_when(
    grepl("CD57pos_DP_NN_CD8",sample_name) ~ "NN DP CD57+",
    grepl("CD57neg_DP_NN_CD8",sample_name) ~ "NN DP CD57-",
    grepl("DN_NN_CD8",sample_name) ~ "NN DN" ,
    
  )) 

#### Load percent EMRA data from Alice W at week 104
# reorganized the data and put in a new spreadsheet
EMRA_104 <- read_xlsx(file.path(annotationDir, "T1DAL_Clusters_from_T1DAL_AbATE_Placebo_8_MCs_for_Erin_2023-10-10_SHEET_3_REORG.xlsx"))

# separate the title
EMRA_104 <- EMRA_104 %>% separate(title , into = c("title","exp","sort","sample","Donor_ID","timepoint"), "_") %>% dplyr::select(-"title", -"exp", -"sort",-"sample")
EMRA_104$Donor_ID <- str_remove(EMRA_104$Donor_ID, "Subj")

# join with libId 
colnames(sample_key_anno)
colnames(EMRA_104)
EMRA_104  <- left_join(EMRA_104 , sample_key_anno) %>%
  ## remove NN
  mutate(population = case_when(
    population =="NN DP CD57+" ~ "DP CD57+",
    population == "NN DP CD57-" ~ "DP CD57-",
    population == "NN DN" ~ "DN" ,
    
  )) 

## Join with CD57 score!
EMRA_104_CD57 <- left_join(EMRA_104,CD57_sig_peaks_score)

### Plot %EMRA by CD57 score per sample
#DN non naive= 7497d9ff, Pd1 = ba4d4cc7, CD57= 93a24eff 

EMRA_104_CD57_plot <- ggplot(EMRA_104_CD57, aes(x = EMRA, y = CD57_score, color = population)) + 
  geom_point() +
  theme_minimal() + 
  geom_smooth(method='lm', formula= y~x) + 
  stat_cor( label.y.npc = 0.98, label.x.npc  = 0.05, show.legend = FALSE, size = 6) +
  scale_color_manual(values= c( "#7497d9ff", "#ba4d4cc7", "#93a24eff" )) +
  labs(x= "% EMRA Cells per Tex Cluster 104 wk",y="CD57+ Tex ATAC Score at 104 wk") +
  theme(text = element_text(size = 16))

## FIGURE S16D
ggsave(EMRA_104_CD57_plot, file = file.path(plotDir, "EMRA_104_CD57_plot.pdf"), width = 7, height = 6)

#save(EMRA_104_CD57, file=file.path(resultDir, "EMRA_104_CD57.Rdata"))

#### FIGURE S16B, FIGURE S16C:PLOT EMRA FLOW DATA ####

load(file=file.path(resultDir, "EMRA_104_CD57.Rdata"))

# Use data already loaded to plot average percent EMRA, EM, CM and NAIVE per population
# across both timepoints

EM_dist_all_time <- EMRA_104_CD57 %>% group_by(population) %>%
  summarize(Naive = mean(`NAÏVE`),
            EM = mean(EM),
            EMRA = mean(EMRA),
            CM = mean(CM)) %>%
  pivot_longer(cols = c(2:5), values_to = "percent_pheno", names_to = "Population")

EM_dist_all_time_plot <- ggplot(EM_dist_all_time, aes(x = population, y = percent_pheno , fill = Population)) + 
  geom_col(position = "fill") + labs(x = element_blank(), y ="% in Phenotype at 104 wk") + 
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(type = c("#8960b3",
                               "#56ae6c",
                               "#ba495b",
                               "#b0923b")) + 
  theme(text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 75, hjust = 1))

## FIGURE S16B
ggsave(EM_dist_all_time_plot, file = file.path(plotDir, "EM_dist_all_time_plot.pdf"), width = 5, height =5)

### Plot EMRA in CD57+ and CD57- per person
EMRA_percent <- EMRA_104_CD57 %>% dplyr::filter(population !="DN") %>%
  ggplot( aes(x = population, y = EMRA)) + geom_boxplot(outlier.shape  = NA,fill = "#ba495b", alpha = 0.8) + 
  geom_point() + stat_compare_means() + 
  facet_grid(.~timepoint) +
  theme_minimal() + 
  labs(y = "% EMRA in Population", x = element_blank()) + 
  scale_y_continuous(limits = c(0,100))+
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 75, hjust = 1))

## FIGURE S16C
ggsave(EMRA_percent, file = file.path(plotDir,"EMRA_percent_CD57_PD1.pdf"), width = 5, height = 6)  

#### FIGURE S16A: Plot the percent of islet antigen specific CD57 cells at baseline  ####

## Load reorganized data 
# - T1DAL_Clusters_from_T1DAL_AbATE_Placebo_8_MCs_for_Erin_2023-10-10.xlsx
# - Copied and pasted the data from the "Raw data - 0 wk" tab, from the bottom section with totals (rows 89-116)
# - reorganized this data for loading into R and saved the new document as "T1DAL_Clusters_from_T1DAL_AbATE_Placebo_8_MCs_for_Erin_2023-10-10_Raw_Data_0wk_REORG.xlsx"
# the donor IDs were in the same order as donor IDs listed above 

Raw_data_0_wk <- read_excel(file.path(annotationDir, "T1DAL_Clusters_from_T1DAL_AbATE_Placebo_8_MCs_for_Erin_2023-10-10_Raw_Data_0wk_REORG.xlsx"))

# plot percent PD1 exhausted and percent CD57 exhausted for each specificity by R vs NR

Raw_data_0_wk_CD57_Ag_specificity_CD57 <- ggplot(Raw_data_0_wk, aes(y = `CD57hi EXHAUSTED`, x =  Response, fill = Response  )) + geom_boxplot(outlier.fill = NULL) +
  geom_point() + theme_minimal() + stat_compare_means() + facet_grid(.~`Antigen Specificity`)  +
  theme(text= element_text(size = 16)) + labs(y = "% CD57+ Tex at Baseline") +
  scale_y_continuous(limits = c(0,100)) +
  scale_fill_discrete(type = c("#b0913b",
                               "#697cd4"))
Raw_data_0_wk_CD57_Ag_specificity_PD1 <- ggplot(Raw_data_0_wk, aes(y = `PD1hi EXHAUSTED`, x =  Response, fill = Response  )) + geom_boxplot(outlier.fill = NULL) +
  geom_point() + theme_minimal() + stat_compare_means() + facet_grid(.~`Antigen Specificity`)  +
  theme(text= element_text(size = 16)) + labs(y = "% PD-1+ Tex at Baseline") +
  scale_y_continuous(limits = c(0,100)) +
  scale_fill_discrete(type = c("#b0913b",
                               "#697cd4"))

### Plot what percent of total CD8 T cells are islet specific at baseline

# From Alice W who put together the data: Values are % of CD8 T cells.  Visit 1 = 0wk; Visit 3 = 104wk.
percent_Ag_of_CD8 <- read_excel(file.path(annotationDir, "Antigen_specific_percent_CD8_AW_2023_10_20.xlsx"))
percent_Ag_of_CD8_baseline <- percent_Ag_of_CD8 %>% separate(ID, into = c("Subject","Visit"), sep = "_") %>% 
  dplyr::filter(Visit == "Visit1") %>% 
  # filter to only keep islet and virus
  dplyr::select(-INSULIN) %>%
  # pivot the data
  pivot_longer(cols = c(4:5), values_to = "% of CD8 at Baseline", names_to = "Antigen Specificity")

# Make plot
Raw_data_0_wk_CD8_Ag_specificity <- ggplot(percent_Ag_of_CD8_baseline, aes(y = `% of CD8 at Baseline`, x =  Response, fill = Response  )) + 
  geom_boxplot(outlier.fill = NULL) +
  geom_point() + theme_minimal() + stat_compare_means() + facet_grid(.~`Antigen Specificity`)  +
  theme(text= element_text(size = 16))  +
  #scale_y_continuous(limits = c(0,100)) +
  scale_fill_discrete(type = c("#b0913b",
                               "#697cd4"))
Raw_data_0_wk_comb <- egg::ggarrange(Raw_data_0_wk_CD8_Ag_specificity,Raw_data_0_wk_CD57_Ag_specificity_CD57,Raw_data_0_wk_CD57_Ag_specificity_PD1, ncol = 3)

# FIGURE S16A
ggsave(Raw_data_0_wk_comb, file = file.path(plotDir, "Raw_data_0_wk_Ag_specificity_all_Tex_CD8.pdf"), height = 5,width = 15)  



