#### P362-1 T1DAL TCR ReAnalysis ####

# this script assesses patterns of TCR sharing and generated TCR figures 

#### README: 
  # The first part of the script uses clone pairs that were not filtered to remove
    # cells with multiple alphas and betas and were determined by tcrGraph.
    # this data was used to calculate Jaccard similarity. 
  # The second part of the script and onwards from the section titled
    # "CALCULATE TCR SHARING MANUALLY USING V,J CDR3_NT AND REPEAT ANALYES BELOW"
    # uses my manual matching of clone pairs and filtering to only keep pairs from
    # cells that had been retained during scRNA-seq QC so that I could then map patterns back
    # onto the UMAPs and other plots and feel more confident about unique sharing. 

#### LOAD LIBRARIES ####

library(plyr)
library(tidyverse)
# Code to set global ggplot theme if you like
library(ggplot2); theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1),
          axis.text=element_text(colour="black"),
          axis.ticks=element_line(colour="black"),
          legend.key = element_blank(),
          text = element_text(size=20),
          strip.text.x = element_text(size = 14,margin = margin( b = 2, t = 2) ),
          strip.background = element_rect(fill="white", colour="black")))

library(ggthemes)
library(ggbeeswarm)
library(viridis)
library(readxl)
library(kableExtra)
library(RColorBrewer)
library(plotly)
library(gtools)
library(egg)
library(apird) #API for the research database
library(devtools)  #if needed to obtain github packages
# install_github('mjdufort/TCRtools') #if needed to get Matt Dufort's package
library(TCRtools) 
library(annotables)
#devtools::install_github("stephenturner/annotables")
# install_github('benaroyaresearch/RNAseQC')
library(RNAseQC) 
library(data.table)
library(edgeR)
library(ggrepel)
library(ComplexHeatmap)
library(geneSetTools) # For barcode plots
library(egg) #For ggarrange
library(ggpubr) #Also for ggarrange
library(DGETools) #install_github('benaroyaresearch/DGETools')
library(inlmisc) #For colors
library(umap)
library(tcrGraph)
library(igraph)
library(ggalluvial)
library(monocle3) # added loading this so I can get colData() and plot_cells()
library(visNetwork)
library(circlize)
library(reshape)
library(reshape2)
library(gtools)
library(igraph)
library(bayesbio)
library(rstatix)
library(gridExtra)
library(agricolae)
library(grid)
library(viridisLite)
library(UpSetR)
library(vegan)
library(fossil)

# set colors for R and NR 
R_NR_colors <- c("#76a44a","#8d70c9")

options(stringsAsFactors = FALSE)
set.seed(42)

#### LOAD ALL SAVED DATA AND SET PATHS ####

# Set Paths 
baseDir <- "/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/T1DAL_cleanup"
plotDir <- file.path(baseDir,'FIGURES')
resultDir <- file.path(baseDir,'SAVED_DATA')
annotationDir <- file.path(baseDir, "RAW_DATA")

## Load and rename cols for annotations
# Load saved annotations 
load(file= file.path(resultDir, "P362-1_annotation.Rdata"))
P362_anno_IDs <- P362_anno %>% dplyr::rename(Lib.ID = "libid", Donor.ID = "participantID") %>% select(Lib.ID, Donor.ID)

# Load TCR annotations
tcr_anno <- P362_anno %>% dplyr::filter(libraryProtocolId=="next_gem_sc_5prime_10x_genomics_tcr_enriched") 

## Load and rename cds object
# Load post QC no MAIT cell data
load(file= file.path(resultDir,"P362-1 T1DAL cds object - postQC no MAIT cells.Rdata"))
cds_all <- cds_no_MAIT
rm(cds_no_MAIT)

## Load saved data that was created in this script 
load(file.path(resultDir,"tcr_output.Rdata"))
tcr_metrics # metrics on the tcrs
tcrs # the original loaded tcr information
all_cloneID_data
# all_cloneID_data includes barcodes, the actual nucleotide sequences for the different chains, and clonotype info. 
cds_metadata_with_clones # contains merged CDS with TCR info 

## load ITN dictionary matching operational donor IDs to the correct public IDs
T1DAL_ITN_ID_dictionary <- readxl::read_xlsx(file.path(annotationDir, "t1dal_mask_id_EW_formatted.xlsx"))
colnames(T1DAL_ITN_ID_dictionary) <- c("operational_PID", "masked_public_PID","Donor.ID")
# remove T1DAL from public_PID column
T1DAL_ITN_ID_dictionary$masked_public_PID <- str_remove(T1DAL_ITN_ID_dictionary$masked_public_PID, "T1DAL_")

# set cluster pallette
cluster_pal <- c("#a6cee3", "#1f78b4", "#b2df8a" ,"#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")

#### LOAD TCR DATA ####

# Load Per cell TCR data

### Get TCR library IDs ###
libs <- tcr_anno %>% dplyr::select(libid)
libs <- as.vector(libs$libid)

# Load and combine all TCR files 
tcr_files <- file.path(annotationDir, "individualLibraryResults","tcr", paste0(libs, "_all_contig_annotations.csv"))
clonotype_files <- file.path(annotationDir, "individualLibraryResults","tcr", paste0(libs, "_clonotypes.csv"))
tcr_qc_files <- file.path(annotationDir, "individualLibraryResults","tcr", paste0(libs, "_metrics_summary.csv"))

unpack_file <- function(file_in){
  
  data <- read.csv(file_in)
  data$libid <- str_extract(file_in, "lib[0-9]+")
  
  return(data)
  
}

tcrs_list <- lapply(tcr_files, unpack_file)
tcrs <- bind_rows(tcrs_list)

clonotypes_list <- lapply(clonotype_files, unpack_file)
clonotypes <- bind_rows(clonotypes_list)

## Extract and bind all tcr metrics
tcr_metrics <- lapply(tcr_qc_files, unpack_file)

## Fix character conversion of estimated number of cells in files
make_char_numeric <- function(file_in){
  in_file <- file_in
  in_file$Estimated.Number.of.Cells <- as.numeric(str_remove(in_file$Estimated.Number.of.Cells,","))
  in_file$Number.of.Cells.With.Productive.V.J.Spanning.Pair <- as.numeric(str_remove(in_file$Number.of.Cells.With.Productive.V.J.Spanning.Pair,","))
  
  return(in_file)
}

tcr_metrics <- lapply(tcr_metrics,make_char_numeric)

tcr_metrics <- bind_rows(tcr_metrics)

#### Calculate Data quality metrics ####

make_pct_numeric <- function(pct_column){
  
  num_column <- pct_column %>%
    str_remove_all("%") %>%
    as.numeric() 
  num_column <- num_column/100
  
  return(num_column)
}

tcr_metrics <- tcr_metrics %>%
  dplyr::mutate(Reads.Mapped.to.Any.V.D.J.Gene = make_pct_numeric(Reads.Mapped.to.Any.V.D.J.Gene),
                Cells.With.Productive.TRA.Contig = make_pct_numeric(Cells.With.Productive.TRA.Contig),
                Cells.With.Productive.TRB.Contig = make_pct_numeric(Cells.With.Productive.TRB.Contig))

tcr_anno_metrics <- merge(tcr_anno,tcr_metrics,by="libid")

#### ASSESS TCR CLONES USING TCR GRAPH ####

P362_anno_to_merge <- P362_anno %>% dplyr::rename(Donor.ID = participantID) %>% select(libid, Donor.ID)
tcrs$barcode_short <- sapply(strsplit(tcrs$barcode,"-"), `[`, 1)

all_libs_tcrs <-  tcrs %>% 
  dplyr::filter(is_cell=="True", productive=="True", cdr3 != "None") %>%
  merge(P362_anno_to_merge,by="libid") %>%
  dplyr::rename(Lib.ID = libid) %>% 
  # dplyr::rename(libid_new = barcode) %>%
  mutate(libid = paste(barcode_short,Donor.ID,sep="_"))

# make the tcr graph
tcr_graph_output_all = tcrGraph::makeTcrGraph(all_libs_tcrs,link="cdr3_nt")

# count TCR clones and get the number of each clone
tcrGraph_clones_all <- tcrGraph::getClonesFromTcrGraph(tcr_graph_output_all, link = "cdr3_nt")
#save(tcrGraph_clones_all,file=file.path(resultDir, "tcrGraph output.Rdata"))

# load previously saved results
load(file=file.path(resultDir, "tcrGraph output.Rdata"))

#### MERGE TCR GRAPH DATA WITH NETWORK DATA ####

### Reformat output from TCR graph to merge with colData(cds_all) 
all_cloneID_data = data.frame(barcode_donorID = NULL, 
                              Donor.ID = NULL,
                              barcode = NULL, 
                              cloneID = NULL, 
                              clone_count = NULL,
                              vGenes = NULL,
                              jGenes = NULL,
                              cdr3_nt = NULL)

for(i in 1:length(unique(tcrGraph_clones_all$cloneId))){
  single_clone_data <- tcrGraph_clones_all[i,] 
  
  ## Extract IDs
  both_IDs = unlist(str_split(single_clone_data$libs, ", "))
  Donor.IDs =  sapply(strsplit(both_IDs, "_"), `[`, 2)
  barcodes_long = sapply(strsplit(both_IDs, "_"), `[`, 1)
  barcodes = sapply(strsplit(barcodes_long,"-"), `[`, 1)
  
  ## Build single clone data frame
  single_clone <- data.frame(barcode_donorID = both_IDs, 
                             Donor.ID = Donor.IDs,
                             barcode_original = barcodes_long,
                             barcode = barcodes,
                             cloneID = single_clone_data$cloneId, 
                             clone_count = single_clone_data$cloneCounts,
                             vGenes = single_clone_data$vGenes,
                             jGenes = single_clone_data$jGenes,
                             cdr3_nt = single_clone_data$cdr3_nt)
  # Build combined TCR dataframe
  all_cloneID_data <- bind_rows(all_cloneID_data,single_clone)
}

# all_cloneID_data includes barcodes, the actual nucleotide sequences for the different chains, and clonotype info. 
View( all_cloneID_data )

# confirm clone counts
all_cloneID_data %>% dplyr::count(cloneID) %>% arrange(desc(n)) %>% head()

## Reformat colData(cds_all) for merging with TCR data
colData(cds_all)$barcode_original <- row.names(colData(cds_all)) # colData(cds_all) works with monocle3 package loading added 
cds_metadata <- colData(cds_all) %>% as.data.frame() %>% select(Donor.ID,barcode_original)
cds_metadata$barcode <-  sapply(strsplit(cds_metadata$barcode_original,"-"), `[`, 1)
cds_metadata$barcode_donorID <- paste(cds_metadata$barcode,cds_metadata$Donor.ID,sep="_")

## Merge colData(cds_all) with TCR info by barcode and Donor ID; if there's a matching barcode and ID between GEX and TCR lib, it will merge, otherwise those fields will be NA
cds_metadata_with_clones <- left_join(cds_metadata,all_cloneID_data, by = "barcode_donorID") #%>% select(-Donor.ID) ## to avoid redundancy in merged colData

## Check that barcode ID order is the same in the cds metadata object and in colData(cds_all)
all(cds_metadata_with_clones$barcode_original == row.names(SummarizedExperiment::colData(cds_all))) # TRUE
## True = TCR and clonotype data is in the same order as the rows in colData(cds)  

### Add new TCR coldata to the cds object ##
## Since barcodes in metadata match exactly with colData(cds), no need to merge (removes row names); can just append desired columns to colData(cds)

colData(cds_all)$cloneID <- cds_metadata_with_clones$cloneID
colData(cds_all)$clone_count <- cds_metadata_with_clones$clone_count
colData(cds_all)$vGenes <- cds_metadata_with_clones$vGenes
colData(cds_all)$cdr3_nt <- cds_metadata_with_clones$cdr3_nt
colData(cds_all)$barcode_donorID <- cds_metadata$barcode_donorID

head(colData(cds_all))

#### Calculate TCR Clone Counts ####

# Calculate clone counts and add to data frame
clone_counts <- colData(cds_all)$clone_count
clone_counts[is.na(clone_counts)] <- 0
colData(cds_all)$clone_count_numeric <- clone_counts
cds_metadata_ordered <- colData(cds_all) %>% as.data.frame() %>% arrange(clone_count_numeric,)

# do these counts match with the calculated clone counts tcrgraph? NO because some of the cells were filtered out in creating the cds object
colData(cds_all) %>% as.data.frame()  %>% dplyr::count(cloneID) %>% arrange(desc(n)) %>% head()

# Full cds metadata with TCR info, cluster number, sequence
cds_metadata_ordered


#### SAVE ANALYSIS OUTPUT DATAFRAMES WITH TCRGRAPH CLONEID ####

#save(tcr_metrics, tcrs, all_cloneID_data, cds_metadata_with_clones, cds_metadata_ordered, file = file.path(resultDir,"tcr_output.Rdata"))
# tcr_metrics # metrics on the tcrs
# tcrs # the original loaded tcr information
# all_cloneID_data
# # all_cloneID_data includes barcodes, the actual nucleotide sequences for the different chains, and clonotype info. 
# cds_metadata_with_clones # contains merged CDS with TCR info from tcrgraph
# cds_metadata_ordered # Full cds metadata with TCR info, cluster number, sequence

# save all_libs_tcrs separately
#save(all_libs_tcrs, file = file.path(resultDir,"all_libs_tcrs.Rdata"))

#### FIGURE 5B: Calculate Jaccard similarity index between clusters using TCRgraph clones ####

# Load previously saved data
load(file = file.path(resultDir, "tcr_output.Rdata"))
load(file = file.path(resultDir,"all_libs_tcrs.Rdata"))

# The Jaccard similarity index = number observations in both sets/ number observations in either set
# The index ranges between 0 and 1 and closer to 1 is a higher Jaccard similarity index

# Write loop to calculate the Jaccard similarity index over each cluster 
cluster_list1 <- levels(cds_metadata_ordered$Cluster.Name)

# get vector of clusters for each cluster
cluster1 <- cds_metadata_ordered %>% filter(Cluster.Name == 1 ) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster1 <- cluster1[,1]
cluster2 <- cds_metadata_ordered %>% filter(Cluster.Name == 2 ) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster2 <- cluster2[,1]
cluster3 <- cds_metadata_ordered %>% filter(Cluster.Name == 3) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster3 <- cluster3[,1]
cluster4 <- cds_metadata_ordered %>% filter(Cluster.Name == 4 ) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster4 <- cluster4[,1]
cluster5 <- cds_metadata_ordered %>% filter(Cluster.Name == 5 ) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster5 <- cluster5[,1]
cluster6 <- cds_metadata_ordered %>% filter(Cluster.Name == 6 ) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster6 <- cluster6[,1]
cluster7 <- cds_metadata_ordered %>% filter(Cluster.Name == 7 ) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster7 <- cluster7[,1]
cluster8 <- cds_metadata_ordered %>% filter(Cluster.Name == 8 ) %>% select(cloneID) %>% filter(!is.na(cloneID))
cluster8 <- cluster8[,1]

## Calculate the Jaccard similarity overall between clusters
cluster_list <- list(cluster1 = cluster1,
                     cluster2 = cluster2,
                     cluster3 = cluster3,
                     cluster4 = cluster4,
                     cluster5 = cluster5,
                     cluster6 = cluster6,
                     cluster7 = cluster7,
                     cluster8 = cluster8)


## Create empty dataframe to house jaccard output
compare_clones_jaccard_df <- data.frame(Cluster1 = character(),
                                        Cluster2 = character(),
                                        Jaccard = numeric())

# calculate for cluster 1 and add to df
for (name in names(cluster_list)) {
  c <- name
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster1",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster1, cluster_list[[name]]))
  # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

# calculate for cluster 2 and add to df
for (name in names(cluster_list)) {
  c <- name
  
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster2",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster2, cluster_list[[name]] ) ) # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

# calculate for cluster 3 and add to df
for (name in names(cluster_list)) {
  c <- name
  
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster3",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster3, cluster_list[[name]] ) ) # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

# calculate for cluster 4 and add to df
for (name in names(cluster_list)) {
  c <- name
  
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster4",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster4, cluster_list[[name]] ) ) # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

# calculate for cluster 5 and add to df
for (name in names(cluster_list)) {
  c <- name
  
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster5",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster5, cluster_list[[name]] ) ) # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

# calculate for cluster 6 and add to df
for (name in names(cluster_list)) {
  c <- name
  
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster6",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster6, cluster_list[[name]] ) ) # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

# calculate for cluster 7 and add to df
for (name in names(cluster_list)) {
  c <- name
  
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster7",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster7, cluster_list[[name]] ) ) # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

# calculate for cluster 8 and add to df
for (name in names(cluster_list)) {
  c <- name
  
  # Calculate the Jaccard similarity index for each set of clones
  new <- data.frame( Cluster1 = "cluster8",
                     Cluster2 = c,
                     Jaccard = jaccardSets(cluster8, cluster_list[[name]] ) ) # Create new row
  compare_clones_jaccard_df[nrow(compare_clones_jaccard_df) + 1, ] <- new 
}

compare_clones_jaccard_df
#save(compare_clones_jaccard_df, file = file.path(resultDir ,"compare_clones_jaccard_df.Rdata"))

### Create jaccard heatmap ###

load(file = file.path(resultDir ,"compare_clones_jaccard_df.Rdata"))

# set the cluster1 as column and cluster 2 as rows with values in the middle
compare_clones_jaccard_df_spread <- spread(compare_clones_jaccard_df, Cluster1, Jaccard, fill = 0)
rownames(compare_clones_jaccard_df_spread) <- compare_clones_jaccard_df_spread$Cluster2 
compare_clones_jaccard_df_spread <- compare_clones_jaccard_df_spread[,-1]

# make into matrix 
compare_clones_jaccard_df_mat <- as.matrix(compare_clones_jaccard_df_spread)

# Plot heatmap with redundant values removed
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
lower_jaccard <- get_lower_tri(compare_clones_jaccard_df_mat)
upper_jaccard <- get_upper_tri(compare_clones_jaccard_df_mat)

## Plot FIGURE 5B
pdf(file = file.path(plotDir,"jaccard_similarity_heatmap_upper_purple.pdf"), width = 4.5, height = 4)
ComplexHeatmap::Heatmap(upper_jaccard, 
                        heatmap_legend_param =  list(title= "TCR\nRepertoire\nOverlap"),
                        col = c("#F7F7F7"  , "#b44dc6c7","#762A83"))
dev.off()

#### CALCULATE TCR SHARING MANUALLY USING V,J CDR3_NT  ####

# This code determines TCR clones using my manual method rather than the tcrgraph method

load(file = file.path(resultDir,"all_libs_tcrs.Rdata"))

#View(all_libs_tcrs)
nrow(all_libs_tcrs) # total rows 26404
# are there multiple tcrs per cell?
length(unique(all_libs_tcrs$barcode)) # 18200, most have an alpha and a beta chain

# Find the unique cdr3nt sequences
all_libs_tcrs_unique_cdr3nt <- all_libs_tcrs %>% distinct(cdr3_nt) 
nrow(all_libs_tcrs_unique_cdr3nt) # 15898 unique cdr3nt sequences! 

# Get all the possible V, J, and CDR3 combinations and how many times they occur total across the dataset
tcr_v_j_cdr3_comb <- ddply(all_libs_tcrs,.(v_gene, j_gene,cdr3_nt), plyr::summarize, sum = length(cdr3_nt)) %>%
  arrange(desc(sum))
# note: this is regardless of whether the sequence is alpha or beta chain 
nrow(tcr_v_j_cdr3_comb) # 15925
nrow(all_libs_tcrs) # TCRs
length(unique(all_libs_tcrs$barcode)) #  18200

# add arbitrary ID to the TCR clones - but the lowest numbers have the highest number of clones 
tcr_v_j_cdr3_comb_id <- tcr_v_j_cdr3_comb %>% mutate(tcr_id = paste0("VJcdr3nt_", row.names(.)))

# join back with list of all tcrs - with each one in a different row
tcr_v_j_cdr3_comb_all <- left_join(all_libs_tcrs, tcr_v_j_cdr3_comb_id)
# Joining, by = c("v_gene", "j_gene", "cdr3_nt")

## Create barcode_donorID column for joining with metadata 
# the libid column of tcr_v_j_cdr3_comb_all is the barcode_donorID! change the name
colnames(tcr_v_j_cdr3_comb_all)[22] <- "barcode_donorID" 
nrow(tcr_v_j_cdr3_comb_all) # 26404
length(unique(tcr_v_j_cdr3_comb_all$barcode_donorID)) # 18413 # 
# there are multiple TCRs from the same cell- alpha and beta chains! 

## Spread the alpha and beta lists so you can recombined and determine combined clonotype
## Split into alpha and beta lists
# save one version for later
tcr_v_j_cdr3_comb_all_TRA_full_seq <- tcr_v_j_cdr3_comb_all %>% filter(chain == "TRA")
nrow(tcr_v_j_cdr3_comb_all_TRA_full_seq) # 7573
tcr_v_j_cdr3_comb_all_TRB_full_seq <- tcr_v_j_cdr3_comb_all %>% filter(chain == "TRB")
nrow(tcr_v_j_cdr3_comb_all_TRB_full_seq) # 18821

# generate version for further editing
tcr_v_j_cdr3_comb_all_TRA <- tcr_v_j_cdr3_comb_all %>% filter(chain == "TRA")
nrow(tcr_v_j_cdr3_comb_all_TRA) # 7573
tcr_v_j_cdr3_comb_all_TRB <- tcr_v_j_cdr3_comb_all %>% filter(chain == "TRB")
nrow(tcr_v_j_cdr3_comb_all_TRB) # 18821
# 7573+18821 = 26394 = we are missing ten cells now?

## Subset and rename columns to label TRA and TRB
tcr_v_j_cdr3_comb_all_TRA <- tcr_v_j_cdr3_comb_all_TRA %>% select(Lib.ID, barcode_donorID, Donor.ID, cdr3_nt, sum, tcr_id)
colnames(tcr_v_j_cdr3_comb_all_TRA)[4:6] <- paste(colnames(tcr_v_j_cdr3_comb_all_TRA)[4:6], "TRA", sep = "_")

# how many cell have double alpha?
double_alpha <- tcr_v_j_cdr3_comb_all_TRA %>% dplyr::count(barcode_donorID) %>% filter(n >=2) # 356 cells have multiple alpha chains
nrow(double_alpha)

tcr_v_j_cdr3_comb_all_TRB <- tcr_v_j_cdr3_comb_all_TRB %>% select(Lib.ID, barcode_donorID, Donor.ID, cdr3_nt, sum, tcr_id)
colnames(tcr_v_j_cdr3_comb_all_TRB)[4:6] <- paste(colnames(tcr_v_j_cdr3_comb_all_TRB)[4:6], "TRB", sep = "_")

# how many cells have double beta?
double_beta <- tcr_v_j_cdr3_comb_all_TRB %>% dplyr::count(barcode_donorID) %>% filter(n >=2) # 1042 cells have multiple beta chains - kind of a lot!
nrow(double_beta)

# Are there any cells with muliple alpha and multiple beta chains?
nrow(double_alpha[double_alpha$barcode_donorID %in% double_beta$barcode_donorID,]) #50 
nrow(double_beta[double_beta$barcode_donorID %in% double_alpha$barcode_donorID,]) #50 

# cells with multiple beta should be removed because those are likely doublet cells!

# rbind these 
tcr_v_j_cdr3_comb_all_chain <- full_join(tcr_v_j_cdr3_comb_all_TRA , tcr_v_j_cdr3_comb_all_TRB) # Joining, by = c("Lib.ID", "barcode_donorID", "Donor.ID")
nrow(tcr_v_j_cdr3_comb_all_chain) # 19855 - ther are only 18200 cells, so there must be TCRs with multiple alpha or beta sequences
colnames(tcr_v_j_cdr3_comb_all_chain)

## Remove cells with more than one beta chain since these are likely due to doublet cells
tcr_v_j_cdr3_comb_all_chain_db_beta_rm <- tcr_v_j_cdr3_comb_all_chain[!(tcr_v_j_cdr3_comb_all_chain$barcode_donorID %in% double_beta$barcode_donorID),] 
nrow(tcr_v_j_cdr3_comb_all_chain_db_beta_rm) # 17671

## Remove cells with more than one alpha chain to increase ability to interpret the results and track clonal expansion
tcr_v_j_cdr3_comb_all_chain_db_a_b_rm <- tcr_v_j_cdr3_comb_all_chain_db_beta_rm[!(tcr_v_j_cdr3_comb_all_chain_db_beta_rm$barcode_donorID %in% double_alpha$barcode_donorID),]
nrow(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm ) # 17059

## Keep cells that only have an alpha beta pair 
tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_pair <- tcr_v_j_cdr3_comb_all_chain_db_a_b_rm %>% filter(!is.na(tcr_id_TRA) & !is.na(tcr_id_TRB))
nrow(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_pair) # 5920
any(duplicated(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_pair$barcode_donorID)) # FALSE

## How many cells only had one alpha with no beta, or only one beta with no alpha
nrow(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm ) # 17059
tcr_v_j_cdr3_comb_all_chain_db_a_b_rm %>% filter(is.na(tcr_id_TRA)) %>% nrow() #10557 cells had one beta but no alpha
tcr_v_j_cdr3_comb_all_chain_db_a_b_rm %>% filter(is.na(tcr_id_TRB)) %>% nrow() #582 cells had one alpha but no beta

## Determine clonotype id based on unique alpha beta pairs
tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype <- tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_pair %>% distinct(tcr_id_TRA, tcr_id_TRB) %>%
  mutate(clonotype_id = paste(tcr_id_TRA, tcr_id_TRB, sep = ":"))
nrow(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype) # 4092 unique clonotypes

## Add back clonotype id to main dataframe to get final clonotype df for analysis below
tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype <- left_join(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_pair,  tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype)
nrow(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype) # 5920
any(duplicated(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype$barcode_donorID)) # FALSE

## Join back clonotype info and filter out cells with multiple alphas and beta in the dfs with the full sequence info 
#will use this later when looking at sequences and specificity
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered <- left_join(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype,  
                                                         tcr_v_j_cdr3_comb_all_TRA_full_seq[!(tcr_v_j_cdr3_comb_all_TRA_full_seq$barcode_donorID %in% double_alpha$barcode_donorID),])
# Joining, by = c("Lib.ID", "barcode_donorID", "Donor.ID")
nrow(tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered)# 5920

tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered <- left_join(tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype,
                                                         tcr_v_j_cdr3_comb_all_TRB_full_seq[!(tcr_v_j_cdr3_comb_all_TRB_full_seq$barcode_donorID %in% double_beta$barcode_donorID),])
# Joining, by = c("Lib.ID", "barcode_donorID", "Donor.ID")
nrow(tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered) # 5920

## resave final df with less complicated name
tcr_clones <- tcr_v_j_cdr3_comb_all_chain_db_a_b_rm_clonotype
#save(tcr_clones,tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered,tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered , file = file.path(resultDir,"tcr_clones.RData"))

#### Analysis of Expanded and Private vs Public TCRs - using manual assignment ####

# load previously saved data
load(file.path(resultDir, "tcr_clones.RData"))

## How many tcrs are public vs private. Public = full alpha beta chain combination found in multiple donors
tcr_clones_public_private <- tcr_clones %>% ungroup() %>% distinct(clonotype_id, Donor.ID) %>% 
  group_by(clonotype_id) %>% dplyr::count() %>% 
  # any with more than 1 count are public
  mutate(clone_sharing = case_when(
    n >= 2 ~ "public",
    n == 1 ~ "private")) %>%
  dplyr::rename(donor_count = n)

# How many are public vs private?
tcr_clones_public_private %>% ungroup() %>% dplyr::count(clone_sharing)
#clone_sharing     n
#<chr>         <int>
#  1 private        4089
# 2 public            3

## How many tcrs are public vs private when public defined as shared alpha chains between donors
tcr_clones_public_private_alpha_only  <- tcr_clones %>% ungroup() %>% distinct(tcr_id_TRA, Donor.ID) %>% 
  group_by(tcr_id_TRA) %>% dplyr::count() %>% 
  # any with more than 1 count are public
  mutate(clone_sharing_alpha = case_when(
    n >= 2 ~ "public",
    n == 1 ~ "private")) %>%
  dplyr::rename(donor_count = n)

# rejoin public and private status to df - using definition with only shared alpha
tcr_clones_sharing <- left_join(tcr_clones, tcr_clones_public_private_alpha_only) # Joining, by = "tcr_id_TRA"
nrow(tcr_clones_sharing) # 5920

# How many TCRs are expanded? (n >=2) (including both public and private!) - defined by full alpha beta clonotype
tcr_clones_sharing_expanded <-  tcr_clones_sharing %>% group_by(clonotype_id) %>% mutate(clone_count= n()) %>%
  mutate(expanded = case_when(clone_count > 1 ~ "expanded", clone_count == 1 ~ "not_expanded")) # 2472 TCRs have expanded
tcr_clones_sharing_expanded %>% filter(clone_count > 1) %>% distinct(clonotype_id) %>% nrow() # 485 are expanded!
tcr_clones_sharing_expanded %>% ungroup() %>% distinct(clonotype_id,expanded) %>% dplyr::count(expanded)
#expanded         n
#<chr>        <int>
#  1 expanded       485
#2 not_expanded  3607

# How many TCRs are expanded (n>=2) (including public and private) when defined by shared alpha chain?
tcr_clones_sharing_expanded_shared_alpha <- tcr_clones_sharing %>% group_by(tcr_id_TRA) %>% mutate(clone_count_alpha = n()) %>%
  mutate(expanded_alpha = case_when(clone_count_alpha > 1 ~ "expanded", clone_count_alpha == 1 ~ "not_expanded"))

# How many expanded are public vs. private
tcr_clones_sharing_expanded_shared_alpha   %>% distinct(tcr_id_TRA, expanded_alpha, clone_sharing_alpha) %>% ungroup() %>% 
  dplyr::count(expanded_alpha, clone_sharing_alpha) 
#expanded_alpha clone_sharing_alpha     n
#<chr>          <chr>               <int>
#  1 expanded       private               433
#2 expanded       public                 81
#3 not_expanded   private              3055

## final df with expansion and public private both defined by shared alpha
tcr_clones_sharing_expanded_shared_alpha
#save(tcr_clones_sharing_expanded_shared_alpha, file = file.path(resultDir,"tcr_clones_sharing_expanded_shared_alpha.RData"))

#### Join manual clone data with network data  ####

## Reload cds data
# Load and rename cds object
load(file.path(resultDir,"P362-1 T1DAL cds object - postQC no MAIT cells.Rdata"))
cds_all_tcr_manual <- cds_no_MAIT
rm(cds_no_MAIT)

colData(cds_all_tcr_manual) # DataFrame with 27050 rows and 14 columns

## Don't need to reformat the tcr data further since it already has the donor ID and the barcode for joining

## Format cds_metadata for joining with UMAP
colData(cds_all_tcr_manual)$barcode_original <- row.names(colData(cds_all_tcr_manual)) # set rownames to be a column for joining
cds_metadata_tcr_df <- colData(cds_all_tcr_manual) %>% as.data.frame() %>% select(Donor.ID,barcode_original)
cds_metadata_tcr_df$barcode <-  sapply(strsplit(cds_metadata_tcr_df$barcode_original,"-"), `[`, 1)
cds_metadata_tcr_df$barcode_donorID <- paste(cds_metadata_tcr_df$barcode, cds_metadata_tcr_df$Donor.ID, sep="_")
nrow(cds_metadata_tcr_df) # 27050
length(unique(cds_metadata_tcr_df$barcode_donorID)) # 27050
# You should be able to merge with the row.names (barcodes), but KD explicitly created a new column with that barcode ID to make it easier to access via dplyr.

## Merge colData(cds_all) with TCR info (by alpha sharing!_) by barcode and Donor ID; if there's a matching barcode and ID between GEX and TCR lib, it will merge, otherwise those fields will be NA
cds_metadata_with_tcr_manual <- left_join(cds_metadata_tcr_df,tcr_clones_sharing_expanded_shared_alpha, by = c("Donor.ID","barcode_donorID")) #%>% select(-Donor.ID) ## to avoid redundancy in merged colData
nrow(cds_metadata_with_tcr_manual) # 27050

## Check that barcode ID order is the same in the cds metadata object and in colData(cds_all)
all(cds_metadata_with_tcr_manual$barcode_original == row.names(SummarizedExperiment::colData(cds_all_tcr_manual))) # TRUE
## True = TCR and clonotype data is in the same order as the rows in colData(cds)  - and also TRUE because no duplicate joining

# replace NAs for clone count with 0 for plotting in UMAP
cds_metadata_with_tcr_manual$clone_count_alpha <- cds_metadata_with_tcr_manual$clone_count_alpha %>% replace_na(0) 
class(cds_metadata_with_tcr_manual$clone_count_alpha) # integer
cds_metadata_with_tcr_manual$clone_count_alpha <- as.numeric(cds_metadata_with_tcr_manual$clone_count_alpha)

### Add new TCR coldata to the cds object ##
## Since barcodes in metadata match exactly with colData(cds), no need to merge (removes row names); can just append desired columns to colData(cds)
colnames(cds_metadata_with_tcr_manual)

colData(cds_all_tcr_manual)$cdr3_nt_TRA <- cds_metadata_with_tcr_manual$cdr3_nt_TRA
#colData(cds_all_tcr_manual)$sum_TRA <- cds_metadata_with_tcr_manual$sum_TRA
colData(cds_all_tcr_manual)$tcr_id_TRA <- cds_metadata_with_tcr_manual$tcr_id_TRA
colData(cds_all_tcr_manual)$cdr3_nt_TRB <- cds_metadata_with_tcr_manual$cdr3_nt_TRB
#colData(cds_all_tcr_manual)$sum_TRB <- cds_metadata_with_tcr_manual$sum_TRB
colData(cds_all_tcr_manual)$tcr_id_TRB <- cds_metadata_with_tcr_manual$tcr_id_TRB
colData(cds_all_tcr_manual)$clonotype_id <- cds_metadata_with_tcr_manual$clonotype_id
colData(cds_all_tcr_manual)$donor_count <- cds_metadata_with_tcr_manual$donor_count
colData(cds_all_tcr_manual)$clone_sharing_alpha <- cds_metadata_with_tcr_manual$clone_sharing_alpha
colData(cds_all_tcr_manual)$clone_count_alpha <- cds_metadata_with_tcr_manual$clone_count_alpha
colData(cds_all_tcr_manual)$expanded_alpha <- cds_metadata_with_tcr_manual$expanded_alpha
colData(cds_all_tcr_manual)$barcode_donorID <- cds_metadata_with_tcr_manual$barcode_donorID

head(colData(cds_all_tcr_manual))

# save
#save(cds_all_tcr_manual, file = file.path(resultDir,"cds_all_tcr_manual.Rdata"))

#### Format MANUAL TCR Clone Counts for plotting ####

# load previously saved data
load( file = file.path(resultDir, "cds_all_tcr_manual.Rdata"))
# make df that has all the cds info, metadata
cds_tcr_manual_metadata <- colData(cds_all_tcr_manual) %>% as.data.frame()

# Re-calculate clone count so that it only includes those cells that have TCR AND are in the UMAP
cds_tcr_manual_metadata_clone_count <- cds_tcr_manual_metadata %>% filter(!is.na(clonotype_id)) %>% group_by(tcr_id_TRA) %>% mutate(clone_count_alpha_UMAP = n()) %>%
  mutate(expanded_alpha_UMAP = case_when(clone_count_alpha > 1 ~ "expanded", clone_count_alpha_UMAP == 1 ~ "not_expanded"))

# rejoin with full metadata
cds_tcr_manual_metadata <- left_join(cds_tcr_manual_metadata, cds_tcr_manual_metadata_clone_count)
cds_tcr_manual_metadata$clone_count_alpha_UMAP <- replace_na(cds_tcr_manual_metadata$clone_count_alpha_UMAP, 0) # replace Nas with zero for plotting
cds_tcr_manual_metadata$clone_count_alpha_UMAP  <- as.numeric(cds_tcr_manual_metadata$clone_count_alpha_UMAP )

# order by clone count alpha UMAP
cds_tcr_manual_metadata_ordered <- cds_tcr_manual_metadata %>% arrange(clone_count_alpha_UMAP)

# add clone_count_alpha_UMAP to TCR manual
all(cds_tcr_manual_metadata$barcode_original == row.names(SummarizedExperiment::colData(cds_all_tcr_manual))) # TRUE
colData(cds_all_tcr_manual)$clone_count_alpha_UMAP <- cds_tcr_manual_metadata$clone_count_alpha_UMAP

# order cds by UMAP clone count
cds_all_tcr_manual_ordered <- cds_all_tcr_manual[,cds_tcr_manual_metadata_ordered$barcode_original]

#save(cds_all_tcr_manual_ordered, file = file.path(resultDir,"cds_all_tcr_manual_ordered.Rdata"))
#save(cds_tcr_manual_metadata, cds_tcr_manual_metadata_ordered, file =file.path(resultDir,"cds_tcr_manual_metadata.RData"))

#### Clone count re-analysis with manual TCR counts ####

# load previously saved data
load(file.path(resultDir, "cds_all_tcr_manual_ordered.Rdata"))

# Which cells in each cluster are part of any expanded clone based on alpha chain sharing - excluding cluster 9
cells_in_clones_tcr_manual <- cds_tcr_manual_metadata_ordered   %>%
  # keep only those with more than 1 clone defined by shared alpha
  filter(clone_count_alpha_UMAP > 1) %>% 
  # remove cluster 9
  filter(Cluster.Name != "9")
nrow(cells_in_clones_tcr_manual) # WAS 1919 before using defintion based on expansion after joining with UMAP, is now is 1877
#save(cells_in_clones_tcr_manual, file = file.path(resultDir,"cells_in_clones_tcr_manual.Rdata"))

# for use later - Which cells in each cluster are part of any expanded clone based on clonotype id sharing
cells_in_clones_tcr_manual_clonotype <- cds_tcr_manual_metadata_ordered   %>%
  filter(!is.na(clonotype_id)) %>% group_by(clonotype_id) %>% mutate(clonotype_count = n()) %>%
  filter(clonotype_count > 1) %>%
  # remove cluster 9
  filter(Cluster.Name != "9")

nrow(cells_in_clones_tcr_manual_clonotype) # 1767

# load previously saved object if necessary
load(file = file.path(resultDir, "cells_in_clones_tcr_manual.Rdata"))

# Count the number of cells in an expanded clone per donor, response, and cluster (since each cell is a row in the metadata)
cells_in_clones_tcr_manual_grouped <- cells_in_clones_tcr_manual %>% group_by(Donor.ID, Response, Cluster.Name) %>% dplyr::count() 
#save(cells_in_clones_tcr_manual_grouped, file = file.path(resultDir,"cells_in_clones_tcr_manual_grouped.Rdata"))

#### FIGURE S10: Plot distribution of TCRs across patients ####

load(file.path(resultDir, "cds_all_tcr_manual_ordered.Rdata"))

# total tcrs per person from those cells that are in the UMAP
total_tcrs_per_person_in_UMAP <- cds_tcr_manual_metadata_ordered %>%  filter(!is.na(clonotype_id)) %>% dplyr::count(Donor.ID)

# add public donor ID for plotting
total_tcrs_per_person_in_UMAP <- left_join(total_tcrs_per_person_in_UMAP , T1DAL_ITN_ID_dictionary)

total_tcrs_per_person_in_UMAP %>% summarize(mean = mean(n), sd = sd(n), min = min(n),max = max(n))
#mean       sd min max
#1 368.3333 235.3552  56 808

total_UMAP_TCRs <- ggplot(total_tcrs_per_person_in_UMAP, aes(x = masked_public_PID, y = n)) + geom_col() +
  labs(x = "Donor", y = "Total UMAP Cells with TCR") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text= element_text(size = 20)) 

## FIGURE S10A
ggsave(total_UMAP_TCRs, file = file.path(plotDir, "total_UMAP_TCRs.pdf"), height = 5, width = 7)

# total tcrs with expanded alpha per donor
TCRs_per_donor_expanded_alpha <- cells_in_clones_tcr_manual %>% dplyr::count(Donor.ID) 

# set donor colors
donor_colors <- c("#bb5595", "#8fb03d", "#6677db", "#c99f3b", "#583687", "#66a453", "#c26ac6", "#45c097", "#b64468", "#848bd3", "#ab7434", "#bb4c41")

# plotting percent of cell per donor per cluster in the UMAP as a stacked bar plot
Donor_ID_cell_number <- cds_tcr_manual_metadata_ordered %>%  filter(!is.na(clonotype_id)) %>% filter(Cluster.Name != 9)
Donor_ID_cell_number$Cluster.Name <- factor(Donor_ID_cell_number$Cluster.Name, levels = c("8","7","6","5","4","3","2","1"))

percent_cells_per_donor_per_cluster <- Donor_ID_cell_number %>%
  group_by(Cluster.Name) %>%
  mutate(total_cells_cluster = n()) %>%
  ungroup() %>%
  group_by(Cluster.Name, Donor.ID) %>%
  mutate(total_per_cluster = n()) %>%
  distinct(Donor.ID, Cluster.Name, total_cells_cluster, total_per_cluster) %>%
  mutate(percent_per_cluster = total_per_cluster/total_cells_cluster*100) 

percent_cells_per_donor_per_cluster <- left_join(percent_cells_per_donor_per_cluster, T1DAL_ITN_ID_dictionary)

# plot total cells per cluster per donor
total_cells_per_donor_per_cluster <- percent_cells_per_donor_per_cluster %>%
  ggplot(aes(x = Cluster.Name, y = total_per_cluster, color =masked_public_PID)) + geom_point(size =2) + 
  scale_color_manual(values = donor_colors) + labs(x = "Cluster",y = "UMAP Cells per Donor with TCR", color = "Donor ID") + 
  theme_bw() + theme(text= element_text(size = 24))

percent_cells_per_donor_per_cluster_plot <- percent_cells_per_donor_per_cluster %>%
  ggplot(aes(x = Cluster.Name, y = percent_per_cluster, fill =masked_public_PID)) + geom_col() + 
  scale_fill_manual(values = donor_colors) + labs(x = "Cluster",y = "UMAP Cells with TCR per Donor (%)", fill = "Donor ID") + 
  theme_bw() + theme(text= element_text(size = 24))

## FIGURE S10B
# export plots of total cells with TCR 
ggsave(total_cells_per_donor_per_cluster, file = file.path(plotDir,"total_cells_per_donor_per_cluster.pdf"), height = 7, width = 7, device  = "pdf")

## FIGURE S10C
ggsave(percent_cells_per_donor_per_cluster_plot, file = file.path(plotDir, "percent_cells_per_donor_per_cluster_plot.pdf"), device  = "pdf",
       width = 7, height = 7)


#### FIGURE 5A: Are there more cells overall in the exhausted clusters than the precursor clusters? #### 

# load previously saved data
load(file = file.path(resultDir,"cells_in_clones_tcr_manual_grouped.Rdata"))

cells_in_clones_tcr_manual_grouped_ex_precursor <- cells_in_clones_tcr_manual_grouped %>% mutate(cluster_type = case_when(
  Cluster.Name %in% c(5,6,7,8) ~ "Partially-Exhausted",
  Cluster.Name %in% c(1,2,3,4) ~ "Non-Exhausted")) %>%
  group_by(Donor.ID, cluster_type) %>%
  summarize(cluster_type_sum = sum(n))

# run T_test
cells_in_clones_tcr_manual_grouped_ex_precursor_t_test <- t.test(cluster_type_sum ~ cluster_type, data= cells_in_clones_tcr_manual_grouped_ex_precursor )
# p = 0.037

cells_in_clones_tcr_manual_grouped_ex_precursor$cluster_type <- factor(cells_in_clones_tcr_manual_grouped_ex_precursor$cluster_type,
                                                                       levels = c("Non-Exhausted"   ,
                                                                                  "Partially-Exhausted"),
                                                                       labels = c("Non-Exhausted\nClusters 1-4", "Exhausted\nClusters\n5-8"))

# normalize values by total number of TCRs
cells_in_clones_tcr_manual_grouped_ex_precursor_donor <- left_join(cells_in_clones_tcr_manual_grouped_ex_precursor, total_tcrs_per_person_in_UMAP) %>%
  mutate(percent_total_expanded =cluster_type_sum/n*100 )

# plot
cells_in_clones_tcr_manual_grouped_ex_precursor_barplot_percent <- ggplot(cells_in_clones_tcr_manual_grouped_ex_precursor_donor, aes(x = cluster_type, y = percent_total_expanded, fill = cluster_type)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point()  +
  labs(x = NULL, y = "% Cells with Expanded TRA", fill = "Cluster Group") + 
  stat_compare_means(method = "t.test") +
  scale_fill_manual(values = c("#7f63b8","#ba6437")) + 
  scale_y_continuous(limits = c(0,100)) +
  theme(text = element_text(size = 16),  legend.position = "none")

## FIGURE 5A
ggsave(cells_in_clones_tcr_manual_grouped_ex_precursor_barplot_percent, file = file.path(plotDir,"cells_in_clones_tcr_manual_grouped_ex_precursor_barplot_percent.pdf"),
       width =3.5, height = 4)

#### FIGURE 5C, FIGURE 5F, FIGURE S12B: Create airline plot using manual TCR clonotype defined by alpha chain sharing  ####
load(file.path(resultDir, "cds_all_tcr_manual.Rdata"))

# Use the UMAP coldata where the TCR information has been joined

# make a copy of the annotation plus UMAP coordinates, for easier manipulation
data.tmp_manual <-
  as.data.frame(colData(cds_all_tcr_manual)) %>%
  # Combine the colData with the reduced dimension representation coordinated for these clones
  cbind(
    as.data.frame(reducedDims(cds_all_tcr_manual)$UMAP) %>%
      magrittr::set_colnames(c("V1", "V2"))) %>% 
  dplyr::rename("clone_id_tcr_graph_clonal_expansion" = "tcr_id_TRA")

# create data frame to store links
curves.tmp_manual <-
  data.frame(
    clone_id_tcr_graph_clonal_expansion = character(),
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric())

# loop over each clone, and extract UMAP coordinates for cells from the same clone
for (clone_id.tmp_manual in na.omit(unique(data.tmp_manual$clone_id_tcr_graph_clonal_expansion))) {
  clone_id_curves.tmp_manual <- curves.tmp_manual[0,]
  data_for_curves.tmp_manual <-
    data.tmp_manual %>%
    dplyr::filter(clone_id_tcr_graph_clonal_expansion %in% clone_id.tmp_manual)
  if (nrow(data_for_curves.tmp_manual) > 1) {
    for (i in 1:(nrow(data_for_curves.tmp_manual)-1)) {
      for (j in (i+1):nrow(data_for_curves.tmp_manual)) {
        clone_id_curves.tmp_manual <-
          rbind(
            clone_id_curves.tmp_manual,
            list(
              clone_id_tcr_graph_clonal_expansion =
                data_for_curves.tmp_manual$clone_id_tcr_graph_clonal_expansion[i],
              x = data_for_curves.tmp_manual$V1[i],
              y = data_for_curves.tmp_manual$V2[i],
              cluster1 = data_for_curves.tmp_manual$Cluster.Name[i],
              xend = data_for_curves.tmp_manual$V1[j],
              yend = data_for_curves.tmp_manual$V2[j],
              cluster2 = data_for_curves.tmp_manual$Cluster.Name[j]))
      }
    }
  }
  curves.tmp_manual <-
    rbind(curves.tmp_manual, clone_id_curves.tmp_manual)
}

# code to color the curves by averaging the colors of the clusters
# not used; does not work if coloring the points, as ggplot does not easily allow two separate color palettes
# curves.tmp$curve_color <- as.character(NA)
# for (i in 1:nrow(curves.tmp)) {
#   curves.tmp$curve_color[i] <-
#     average_colors(
#       c(pal.cluster_renumbered[as.character(curves.tmp$cluster1[i])],
#         pal.cluster_renumbered[as.character(curves.tmp$cluster2[i])]))
# }

# set clor pallette
pal.cluster_renumbered = toupper(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffed6f','#b15928', "gray","black","blue","red")) 

## FIGURE 5C make airline plot of the full TCR sharing cluster 9 removed
full_TCR_airline_plot_manual_no_9 <- plot_cells(
  cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% filter(Cluster.Name != 9))],
  color_cells_by = "Cluster.Name", cell_size=1, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  geom_curve(
    data = curves.tmp_manual,
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.07,
    alpha=0.5) +
  theme(text = element_text(size = 24))

## plot Figure 5C
ggsave(full_TCR_airline_plot_manual_no_9, file = file.path(plotDir,"full_TCR_airline_plot_manual_ALPHA_no_9.pdf"), device = "pdf",height = 8, width = 8)

## FIGURE S12B: plot connections to each cluster individually
subset_clusts_curves = c(4,5)
cluster_4_5_TCR_manual <- plot_cells(
  cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% filter(Cluster.Name != 9))],
  color_cells_by = "Cluster.Name", cell_size=1,alpha=0.8, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(subset_clusts_curves) & cluster2 %in% c(subset_clusts_curves)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  labs(title = "Precursor 4-5 Sharing")

subset_clusts_curves = c(4,6)
cluster_4_6_TCR_manual <- plot_cells(
  cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% filter(Cluster.Name != 9))],
  color_cells_by = "Cluster.Name", cell_size=1,alpha=0.8, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(subset_clusts_curves) & cluster2 %in% c(subset_clusts_curves)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  labs(title = "Precursor 4-6 Sharing")

subset_clusts_curves = c(4,7)
cluster_4_7_TCR_manual <- plot_cells(
  cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% filter(Cluster.Name != 9))],
  color_cells_by = "Cluster.Name", cell_size=1,alpha=0.8, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(subset_clusts_curves) & cluster2 %in% c(subset_clusts_curves)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1)+
  labs(title = "Precursor 4-7 Sharing")

subset_clusts_curves = c(4,8)
cluster_4_8_TCR_manual <- plot_cells(
  cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% filter(Cluster.Name != 9))],
  color_cells_by = "Cluster.Name", cell_size=1,alpha=0.8, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(subset_clusts_curves) & cluster2 %in% c(subset_clusts_curves)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1)+
  labs(title = "Precursor 4-8 Sharing")

precursor_clusters <- egg::ggarrange(cluster_4_5_TCR_manual, cluster_4_6_TCR_manual,cluster_4_7_TCR_manual,cluster_4_8_TCR_manual)

# FIGURE S12B
ggsave(precursor_clusters , file = file.path(plotDir,"precursor_clusters_all.tiff"), device = "tiff",height = 10, width = 10)

### FIGURE 5F: Plot all connections that start or end in 4 or one Ex cluster

cluster_4_to_Tex_TCR_manual <- plot_cells(
  cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% filter(Cluster.Name != 9))],
  color_cells_by = "Cluster.Name", cell_size=1,alpha=0.8, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  # plot 4 to 5
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(4) & cluster2 %in% c(5)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1)+
  # plot 4 to 6
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(4) & cluster2 %in% c(6)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  # plot 4 to 7
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(4) & cluster2 %in% c(7)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  # plot 4 to 8
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(4) & cluster2 %in% c(8)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  ### repeat with inverse
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(5) & cluster2 %in% c(4)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1)+
  # plot 4 to 6
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(6) & cluster2 %in% c(4)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  # plot 4 to 7
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(7) & cluster2 %in% c(4)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  # plot 4 to 8
  geom_curve(
    data = curves.tmp_manual %>% dplyr::filter(cluster1 %in% c(8) & cluster2 %in% c(4)),
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.2,
    alpha=0.1) +
  labs(title = "Precursor 4 to Any Tex Sharing")

# FIGURE 5F
ggsave(cluster_4_to_Tex_TCR_manual, file = file.path(plotDir,"cluster_4_to_Tex_TCR_manual.tiff"), device = "tiff",height = 10, width = 10)

#### Set up data for downstream plotting ####

# initially this code was created to make Circos plots but
# I decided not to use the Circos plots in downstream analysis
# so I have only kept the code that created data objects that were used downstream

# load previously saved data
load(file.path(resultDir, "tcr_clones.RData"))

# start with filtered TRA and TRB list objects that have the clonotype info and had cells with multiple alphas and betas removed  

# filter for only those cells in the UMAP
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered %>%
  filter(barcode_donorID %in% cds_tcr_manual_metadata[!is.na(cds_tcr_manual_metadata$clonotype_id),]$barcode_donorID)
tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP <- tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered %>%
  filter(barcode_donorID %in% cds_tcr_manual_metadata[!is.na(cds_tcr_manual_metadata$clonotype_id),]$barcode_donorID)

# Join the cluster information for plotting colors
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP <- left_join(tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP, cds_tcr_manual_metadata[,c("barcode_donorID", "Cluster.Name")])
tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP <- left_join(tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP, cds_tcr_manual_metadata[,c("barcode_donorID", "Cluster.Name")])

# change barcode column name
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP %>% dplyr::rename(tcr1 = barcode_donorID)
tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP_circos <- tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP %>% dplyr::rename(tcr1 = barcode_donorID)
class(tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos$Cluster.Name)

# Make cluster a character so it can be plotted
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos$Cluster.Name <- as.character(tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos$Cluster.Name)
tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP_circos$Cluster.Name <- as.character(tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP_circos$Cluster.Name)

# Make df housing cluster names and colors of choice 
cluster_col <- data.frame(Cluster.Name = c("1","2","3","4","5","6","7","8"),
                          cluster_col = c("#b74560","#4bc490","#644196","#c1a13c","#7087de","#729b45","#c060a7","#b75b36"))

# Make named vector of cluster and df names for producing legend
cluster_col_vec <- structure(cluster_col$cluster_col, names=cluster_col$Cluster.Name, class="character")

# Join df with colors onto cell dataframe
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos <- left_join(tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos, cluster_col)
tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP_circos <- left_join(tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP_circos, cluster_col)

# save output here for TRA and TRB
#save(tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos ,tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP_circos , file= file.path(resultDir, "TCR_UMAP_circos_data.Rdata") )

#### FIGURE 5E, FIGURE S12A: Alpha sharing: Common precursor populations? Count number of cells shared between clusters by donor ####

# load previously saved data
load(file= file.path(resultDir, "TCR_UMAP_circos_data.Rdata"))
load(file =file.path(resultDir,"cds_tcr_manual_metadata.RData"))

donor_list <- unique(tcr_anno$participantID)

# Match alpha chains in a loop by donor
match_TCR_donor <- function(donor, data) {
  data_filter <- data %>% filter(Donor.ID == donor)
  match <- TCRtools::match_TCR_chains(data_filter, id_col= "tcr1", junction_col = "cdr3_nt")
  match$Donor.ID <- donor 
  match
}

match_TCR_alpha_donor <- lapply(donor_list, match_TCR_donor, tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos)
match_TCR_alpha_donor_all <- bind_rows(match_TCR_alpha_donor)

# Join cluster info for both tcr1 and tcr2
match_TCR_alpha_donor_all_cluster <- match_TCR_alpha_donor_all %>% dplyr::rename(barcode_donorID = tcr1) %>% 
  left_join(., cds_tcr_manual_metadata[,c("Cluster.Name", "barcode_donorID", "Response")]) %>% 
  dplyr::rename(Cluster.Name.From = Cluster.Name, tcr1 = barcode_donorID, barcode_donorID = tcr2) %>%
  left_join(., cds_tcr_manual_metadata[,c("Cluster.Name", "barcode_donorID")]) %>% 
  dplyr::rename(Cluster.Name.To = Cluster.Name, tcr2 = barcode_donorID)

# change cluster name
match_TCR_alpha_donor_all_cluster_filter <- match_TCR_alpha_donor_all_cluster %>% filter(Cluster.Name.From != Cluster.Name.To) %>%
  mutate(from_to = paste(Cluster.Name.From, Cluster.Name.To, sep = "-")) 

# save for use downstream
#save(match_TCR_alpha_donor_all_cluster, file = file.path(resultDir, "match_TCR_alpha_donor_all_cluster.RData"))
#save(match_TCR_alpha_donor_all_cluster_filter, file= file.path(resultDir, "match_TCR_alpha_donor_all_cluster_filter.Rdata"))

# calculate total number of links per donor, after aggregating to count the total matches in each category
match_TCR_alpha_donor_all_cluster_filter_total_per_donor <- match_TCR_alpha_donor_all_cluster_filter %>% group_by(Donor.ID) %>%
  summarize(total_matches_all_connections = n())

# Keep only specific links of interest to ask the question about which initial population (1,2,3,4) are the shared precursor
match_TCR_alpha_donor_all_cluster_filter_precursor <- match_TCR_alpha_donor_all_cluster_filter  %>%
  # Keep links with following patterns of sharing to set up for running statistical test
  filter(from_to %in% c("1-5","1-6","1-7","1-8","2-5","2-6","2-7","2-8","3-5", "3-6", "3-7","3-8", "4-5","4-6","4-7", "4-8"))

## Count matches per category for each donor
match_TCR_alpha_donor_all_cluster_filter_match <- match_TCR_alpha_donor_all_cluster_filter_precursor %>%
  group_by(Donor.ID, Response, from_to) %>%
  summarize(number_matches = n())

## Combine values 1 to exhausted, 2 to exhausted, 3 to exhausted, and 4 to exhausted and run ANOVA
match_TCR_alpha_donor_all_cluster_filter_match_precursor <- match_TCR_alpha_donor_all_cluster_filter_match %>%
  mutate(precursor_test = case_when(from_to %in% c("1-5","1-6","1-7","1-8") ~ "1",
                                    from_to %in% c("2-5","2-6","2-7","2-8") ~ "2",
                                    from_to %in% c("3-5", "3-6", "3-7","3-8") ~ "3",
                                    from_to %in% c("4-5","4-6","4-7", "4-8") ~ "4")) %>%
  group_by(Donor.ID, precursor_test) %>%
  summarize(number_matches = sum(number_matches))

# run ANOVA to see if there are differences in the groups
match_TCR_alpha_donor_all_cluster_filter_match_precursor$precursor_test <- factor(match_TCR_alpha_donor_all_cluster_filter_match_precursor$precursor_test ,
                                                                                  levels= c("1","2","3","4"))

match_TCR_alpha_donor_all_cluster_filter_match_precursor_fit <- lm(number_matches~precursor_test, data = match_TCR_alpha_donor_all_cluster_filter_match_precursor)
match_TCR_alpha_donor_all_cluster_filter_match_precursor_aov <- car::Anova(match_TCR_alpha_donor_all_cluster_filter_match_precursor_fit, type = "III")
precursor_HSD <- agricolae::HSD.test(match_TCR_alpha_donor_all_cluster_filter_match_precursor_fit, "precursor_test", group = TRUE)

# AOV and tukey are not significant, likely because of large difference in sample size

# normalize plot by the total amount of between cluster sharing for each donor
match_TCR_alpha_donor_all_cluster_filter_total_per_donor 
match_TCR_alpha_donor_all_cluster_filter_match_precursor_percent_total <- 
  left_join(match_TCR_alpha_donor_all_cluster_filter_match_precursor,match_TCR_alpha_donor_all_cluster_filter_total_per_donor ) %>%
  group_by(Donor.ID, precursor_test) %>%
  mutate(percent_total = number_matches/total_matches_all_connections *100)

### FIGURE 5E
precursor_test_plot_stats_percent_total <- ggplot(match_TCR_alpha_donor_all_cluster_filter_match_precursor_percent_total, aes(x = precursor_test , y = percent_total, colour = precursor_test)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point() + 
  scale_colour_manual(values = cluster_pal[1:4]) +
  theme(legend.position = "none", text = element_text(size = 16)) +
  stat_compare_means()+
  labs(x = "Potential Shared\nLineage Population", y = "% Exhausted Cluster Sharing per Donor")

## FIGURE 5E
ggsave(precursor_test_plot_stats_percent_total , file = file.path(plotDir,"precursor_test_plot_stats_percent_total.pdf"), device = "pdf", height = 5, width = 4)

### Repeat but split out the exhausted cluster matches into CD57 and PD1 and also normalize by the total number per person
## Combine values 1 to exhausted, 2 to exhausted, 3 to exhausted, and 4 to exhausted and run ANOVA
match_TCR_alpha_donor_all_cluster_filter_match_precursor_CD57 <- match_TCR_alpha_donor_all_cluster_filter_match %>%
  mutate(precursor_test = case_when(from_to %in% c("1-5","1-6","1-8") ~ "1",
                                    from_to %in% c("2-5","2-6","2-8") ~ "2",
                                    from_to %in% c("3-5", "3-6","3-8") ~ "3",
                                    from_to %in% c("4-5","4-6", "4-8") ~ "4")) %>%
  filter(!is.na(precursor_test)) %>%
  group_by(Donor.ID, precursor_test) %>%
  summarize(number_matches = sum(number_matches)) %>%
  ungroup() %>%
  left_join(match_TCR_alpha_donor_all_cluster_filter_total_per_donor ) %>%
  group_by(Donor.ID, precursor_test) %>%
  mutate(percent_total = number_matches/total_matches_all_connections *100)


match_TCR_alpha_donor_all_cluster_filter_match_precursor_PD1 <- match_TCR_alpha_donor_all_cluster_filter_match %>%
  mutate(precursor_test = case_when(from_to =="1-7" ~ "1",
                                    from_to == "2-7" ~ "2",
                                    from_to == "3-7" ~ "3",
                                    from_to == "4-7" ~ "4")) %>%
  filter(!is.na(precursor_test)) %>%
  group_by(Donor.ID, precursor_test) %>%
  summarize(number_matches = sum(number_matches)) %>%
  ungroup() %>%
  left_join(match_TCR_alpha_donor_all_cluster_filter_total_per_donor ) %>%
  group_by(Donor.ID, precursor_test) %>%
  mutate(percent_total = number_matches/total_matches_all_connections *100)

# FIGURE S12A: plot for CD57 and PD1 separately
precursor_test_plot_stats_percent_total_CD57 <- ggplot(match_TCR_alpha_donor_all_cluster_filter_match_precursor_CD57 , aes(x = precursor_test , y = percent_total, colour = precursor_test)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point() + 
  scale_colour_manual(values = cluster_pal[1:4]) +
  theme(legend.position = "none", text = element_text(size = 16)) +
  stat_compare_means()+
  labs(x = "Potential Shared\nLineage Population", y = "% CD57+ Tex Cluster Sharing per Donor")

# FIGURE S12A:
ggsave(precursor_test_plot_stats_percent_total_CD57, file = file.path(plotDir,"precursor_test_plot_stats_percent_total_CD57.pdf"), device = "pdf", height = 5, width = 4)

precursor_test_plot_stats_percent_total_PD1 <- ggplot(match_TCR_alpha_donor_all_cluster_filter_match_precursor_PD1 , aes(x = precursor_test , y = percent_total, colour = precursor_test)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point() + 
  scale_colour_manual(values = cluster_pal[1:4]) +
  theme(legend.position = "none", text = element_text(size = 16)) +
  stat_compare_means()+
  labs(x = "Potential Shared\nLineage Population", y = "% PD-1+ Tex Cluster Sharing per Donor")

# FIGURE S12A:
ggsave(precursor_test_plot_stats_percent_total_PD1, file = file.path(plotDir,"precursor_test_plot_stats_percent_total_PD1.pdf"), device = "pdf", height = 5, width = 4)

#### Beta sharing: Common precursor populations? Count number of cells shared between clusters by donor ####

# load previously saved data
load(file= file.path(resultDir, "TCR_UMAP_circos_data.Rdata"))

# Match beta chains in a loop by donor - only repeat the code necessary to generate the downstream network plots
match_TCR_donor <- function(donor, data) {
  data_filter <- data %>% filter(Donor.ID == donor)
  match <- TCRtools::match_TCR_chains(data_filter, id_col= "tcr1", junction_col = "cdr3_nt")
  match$Donor.ID <- donor 
  match
}

match_TCR_beta_donor <- lapply(donor_list, match_TCR_donor, tcr_v_j_cdr3_comb_all_TRB_full_seq_filtered_UMAP_circos)
match_TCR_beta_donor_all <- bind_rows(match_TCR_beta_donor)

# Join cluster info for both tcr1 and tcr2
match_TCR_beta_donor_all_cluster <- match_TCR_beta_donor_all %>% dplyr::rename(barcode_donorID = tcr1) %>% 
  left_join(., cds_tcr_manual_metadata[,c("Cluster.Name", "barcode_donorID", "Response")]) %>% 
  dplyr::rename(Cluster.Name.From = Cluster.Name, tcr1 = barcode_donorID, barcode_donorID = tcr2) %>%
  left_join(., cds_tcr_manual_metadata[,c("Cluster.Name", "barcode_donorID")]) %>% 
  dplyr::rename(Cluster.Name.To = Cluster.Name, tcr2 = barcode_donorID)

# save for future
#save(match_TCR_beta_donor_all_cluster, file = file.path(resultDir, "match_TCR_beta_donor_all_cluster.Rdata"))

#### Make network layout ####

# Create handmade positions for the clusters
x    <- c(5,4.5,4.5,1.7,-3,-3,-3,-5)
y    <- c(4.5,3,1,3,3,1,4.5,3)
grid <- matrix(c(x, y), ncol=2)

x    <- c(1,1,1, 0.7, 0.2, 0.2, 0.2, -0.2)
y    <- c(2.5, 1, -1, 1,1,-1,2.5,1)

#test layout
g <- make_graph(~ 1-2-3-4-5-6-7-8)
plot(g, layout = grid)

#### FIGURE 5D: Network plot of any TRA sharing between or within clusters across all donors ####

load(file = file.path(resultDir, "match_TCR_alpha_donor_all_cluster.RData"))

# remove reciprocal cell matches so those aren't duplicated
match_TCR_alpha_donor_all_cluster_rm_dup <- match_TCR_alpha_donor_all_cluster %>% 
  mutate(key = paste0(pmin(tcr1, tcr2), pmax(tcr1, tcr2), sep = "")) %>% 
  distinct( key, .keep_all = TRUE)

## Create edge df: count up all the connections across all the donors
match_TCR_alpha_donor_all_cluster_count_all <- match_TCR_alpha_donor_all_cluster_rm_dup %>% group_by(Cluster.Name.From, Cluster.Name.To) %>% dplyr::count() %>%
  dplyr::rename(from = Cluster.Name.From, to = Cluster.Name.To, number_sharing = n) %>% mutate(from_to = paste(from, to, sep = "-"))
nrow(match_TCR_alpha_donor_all_cluster_count_all) # 50

# sort to remove reciprocal cluster matches
match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed <- unique(data.frame(t(apply(match_TCR_alpha_donor_all_cluster_count_all[,c(1,2)], 1, sort))))
colnames(match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed) <- c("from","to") 
match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed <- match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed %>% mutate(from_to = paste(from, to, sep = "-"))

# filter original edge df now that reciprocal combinations have been removed
match_TCR_alpha_donor_all_cluster_count_all <- match_TCR_alpha_donor_all_cluster_count_all %>% 
  filter(from_to %in% match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed$from_to)
nrow(match_TCR_alpha_donor_all_cluster_count_all) # 29

# Create node df: of how many cells in each cluster have a tcr - regardless of donor
match_TCR_alpha_donor_all_cluster_donor_tcr_count <- match_TCR_alpha_donor_all_cluster %>% distinct( Cluster.Name.From, tcr1) %>% 
  dplyr::count( Cluster.Name.From) %>% dplyr::rename(from = Cluster.Name.From, number_cells_with_tcr = n)

# create network
all_sharing_net_donor <- graph_from_data_frame(d=match_TCR_alpha_donor_all_cluster_count_all, vertices=match_TCR_alpha_donor_all_cluster_donor_tcr_count , directed=F) 
all_sharing_net_donor_cluster_col <- graph_from_data_frame(d=match_TCR_alpha_donor_all_cluster_count_all, vertices=match_TCR_alpha_donor_all_cluster_donor_tcr_count , directed=F) 

# plot network
# set colors of clusters
network_colors  <- c("#a2539b","#a2539b","#a2539b", "#6973ca", "#93a24e", "#93a24e", "#ba4d4c", "#93a24e")
network_pheno  <- c("#a2539b", "#6973ca", "#93a24e", "#ba4d4c")

V(all_sharing_net_donor)$color <- network_colors
V(all_sharing_net_donor_cluster_col)$color <- cluster_pal

# Set node size based on cell number:
V(all_sharing_net_donor)$size <- V(all_sharing_net_donor)$number_cells_with_tcr/10
V(all_sharing_net_donor_cluster_col)$size <- V(all_sharing_net_donor_cluster_col)$number_cells_with_tcr/10

# Set edge width based on weight:
E(all_sharing_net_donor)$width <- E(all_sharing_net_donor)$number_sharing/100
E(all_sharing_net_donor_cluster_col)$width <- E(all_sharing_net_donor_cluster_col)$number_sharing/100


# FIGURE 5D: plot with new grid layout and colors as cluster colors
pdf(file = file.path(plotDir,"Network_ALL_TRA_combos_GRID.pdf"), width = 8, height = 7.5)
plot(all_sharing_net_donor_cluster_col, layout = grid, vertex.label.font  = 2, vertex.label.color="black", 
     vertex.label.cex = 1.5)
legend(title = "Cluster", x=-2, y=1, legend = c("1","2","3","4","5","6","7","8"), pch=21,
       col="#777777", pt.bg=cluster_pal, pt.cex=2.5, cex=1.5, bty="n", ncol=1)
title(main = "TRA Pairs Shared Within\nand Between Clusters",cex.main = 1.5)
dev.off()

####  FIGURE S11A: Network plot of any TRB sharing between or within clusters across all donors ####

# load previously saved data
load(file = file.path(resultDir, "match_TCR_beta_donor_all_cluster.RData"))
match_TCR_beta_donor_all_cluster

# remove reciprocal cell matches so those aren't duplicated
match_TCR_beta_donor_all_cluster_rm_dup <- match_TCR_beta_donor_all_cluster %>% 
  mutate(key = paste0(pmin(tcr1, tcr2), pmax(tcr1, tcr2), sep = "")) %>% 
  distinct( key, .keep_all = TRUE)
#save(match_TCR_beta_donor_all_cluster_rm_dup, file= file.path(resultDir, "match_TCR_beta_donor_all_cluster_rm_dup.RData"))

## Create edge df: count up all the connections across all the donors
match_TCR_beta_donor_all_cluster_count_all <- match_TCR_beta_donor_all_cluster_rm_dup %>% group_by(Cluster.Name.From, Cluster.Name.To) %>% dplyr::count() %>%
  dplyr::rename(from = Cluster.Name.From, to = Cluster.Name.To, number_sharing = n) %>% mutate(from_to = paste(from, to, sep = "-"))
nrow(match_TCR_beta_donor_all_cluster_count_all) # 49

# sort to remove reciprocal matches
match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed <- unique(data.frame(t(apply(match_TCR_beta_donor_all_cluster_count_all[,c(1,2)], 1, sort))))
colnames(match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed) <- c("from","to") 
match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed <- match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed %>% mutate(from_to = paste(from, to, sep = "-"))

# filter original edge df now that reciprocal combinations have been removed
match_TCR_beta_donor_all_cluster_count_all <- match_TCR_beta_donor_all_cluster_count_all %>% 
  filter(from_to %in% match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed$from_to)
nrow(match_TCR_beta_donor_all_cluster_count_all) # 29

# Create node df: of how many cells in each cluster have a tcr - regardless of donor
match_TCR_beta_donor_all_cluster_donor_tcr_count <- match_TCR_beta_donor_all_cluster %>% distinct( Cluster.Name.From, tcr1) %>% 
  dplyr::count( Cluster.Name.From) %>% dplyr::rename(from = Cluster.Name.From, number_cells_with_tcr = n)

# create network
all_sharing_net_donor_beta <- graph_from_data_frame(d=match_TCR_beta_donor_all_cluster_count_all, vertices=match_TCR_beta_donor_all_cluster_donor_tcr_count , directed=F) 
all_sharing_net_donor_beta_cluster_col <- graph_from_data_frame(d=match_TCR_beta_donor_all_cluster_count_all, vertices=match_TCR_beta_donor_all_cluster_donor_tcr_count , directed=F) 

# plot network
# set colors of clusters
network_colors  <- c("#a2539b","#a2539b","#a2539b", "#6973ca", "#93a24e", "#93a24e", "#ba4d4c", "#93a24e")
network_pheno  <- c("#a2539b", "#6973ca", "#93a24e", "#ba4d4c")

V(all_sharing_net_donor_beta)$color <- network_colors
V(all_sharing_net_donor_beta_cluster_col)$color <- cluster_pal

# Set node size based on cell number:
V(all_sharing_net_donor_beta)$size <- V(all_sharing_net_donor_beta)$number_cells_with_tcr/10
V(all_sharing_net_donor_beta_cluster_col)$size <- V(all_sharing_net_donor_beta_cluster_col)$number_cells_with_tcr/10

# Set edge width based on weight:
E(all_sharing_net_donor_beta)$width <- E(all_sharing_net_donor_beta)$number_sharing/100
E(all_sharing_net_donor_beta_cluster_col)$width <- E(all_sharing_net_donor_beta_cluster_col)$number_sharing/100


#  FIGURE S11A: plot with new grid layout and colors as cluster colors
pdf(file = file.path(plotDir,"Network_ALL_TRB_combos_GRID.pdf"), width = 8.5, height = 7.5)
plot(all_sharing_net_donor_beta_cluster_col, layout = grid, vertex.label.font  = 2, vertex.label.color="black", vertex.label.cex = 1.5)
legend(title = "Cluster", x=-2, y=1, legend = c("1","2","3","4","5","6","7","8"), pch=21,
       col="#777777", pt.bg=cluster_pal, pt.cex=2.5, cex=1.5, bty="n", ncol=1)
title(main = "TRB Pairs Shared Within\nand Between Clusters",cex.main = 1.5)
dev.off()

#### FIGURE S11B: Heatmap of any TRA and TRB sharing between clusters by donor ####

# get within and between cluster sharing per donor

load(file = file.path(resultDir, "match_TCR_alpha_donor_all_cluster.RData"))

match_TCR_alpha_donor_all_cluster

donor_list <- unique(tcr_anno$participantID)

donor_response_match <- match_TCR_alpha_donor_all_cluster %>% distinct(Donor.ID, Response)

TRA_all_donor_sharing <- function(donor, data) {
  
  ## Create edge df: count up all the connections in each donor individually
  ## Create edge df: count up all the connections across all the donors
  match_TCR_alpha_donor_all_cluster_count_all_DONOR <- data %>% 
    filter(Donor.ID == donor) %>%
    group_by(Cluster.Name.From, Cluster.Name.To) %>% dplyr::count() %>%
    dplyr::rename(from = Cluster.Name.From, to = Cluster.Name.To, number_sharing = n) %>%
    mutate(from_to = paste(from, to, sep = "-"), Donor.ID = donor)
  
  # sort to remove reciprocal matches
  match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed_DONOR <- unique(data.frame(t(apply(match_TCR_alpha_donor_all_cluster_count_all_DONOR[,c(1,2)], 1, sort))))
  colnames(match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed_DONOR) <- c("from","to") 
  match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed_DONOR <- match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed_DONOR %>% mutate(from_to = paste(from, to, sep = "-"))
  
  # filter original edge df now that reciprocal combinations have been removed
  match_TCR_alpha_donor_all_cluster_count_all_DONOR <- match_TCR_alpha_donor_all_cluster_count_all_DONOR %>% 
    filter(from_to %in% match_TCR_alpha_donor_all_cluster_count_all_reciprocal_removed_DONOR$from_to)
  
  # Create node df: of how many cells in each cluster have a tcr - in each donor individually
  match_TCR_alpha_donor_all_cluster_donor_tcr_count_DONOR <- data %>% 
    # filter for donor of interest
    filter(Donor.ID == donor) %>% 
    distinct( Cluster.Name.From, tcr1) %>% 
    dplyr::count( Cluster.Name.From) %>% dplyr::rename(from = Cluster.Name.From, number_cells_with_tcr = n)
  
  # return output by donor so I can concatenate
  match_TCR_alpha_donor_all_cluster_count_all_DONOR
}

TRA_all_sharing_donor <- lapply(donor_list, TRA_all_donor_sharing, match_TCR_alpha_donor_all_cluster_rm_dup)
names(TRA_all_sharing_donor ) <- for (donor in donor_list) {
  print(paste0("TRA_all_sharing_", donor))
}
TRA_all_sharing_donor <- bind_rows(TRA_all_sharing_donor)

# save
#save(TRA_all_sharing_donor , file=file.path(resultDir, "TRA_all_sharing_donor.Rdata"))

# TRA
# load previous data calculating sharing per donor
load(file=file.path(resultDir, "TRA_all_sharing_donor.Rdata"))
TRA_all_sharing_donor
class(TRA_all_sharing_donor$number_sharing) # integer

# spread the donor and make the from to as the rownames
TRA_all_sharing_donor_heatmap <- TRA_all_sharing_donor %>% ungroup() %>%
  dplyr::select(-c(from,to)) %>%
  pivot_wider(names_from = Donor.ID, values_from = number_sharing) %>% column_to_rownames(.,var = "from_to")

TRA_all_sharing_donor_heatmap[is.na(TRA_all_sharing_donor_heatmap)] <- 0

# set up heatmap annotation
hc <- columnAnnotation(Response = donor_response_match$Response,
                       col = list(Response = c("R" = "#353a33",
                                               "NR" = "#c6cec4")))
#### Repeat for TRB

load( file= file.path(resultDir, "match_TCR_beta_donor_all_cluster_rm_dup.RData"))

match_TCR_beta_donor_all_cluster_rm_dup

TRB_all_donor_sharing <- function(donor, data) {
  
  ## Create edge df: count up all the connections in each donor individually
  match_TCR_beta_donor_all_cluster_count_all_DONOR <- data %>% 
    filter(Donor.ID == donor) %>%
    group_by(Cluster.Name.From, Cluster.Name.To) %>% dplyr::count() %>%
    dplyr::rename(from = Cluster.Name.From, to = Cluster.Name.To, number_sharing = n) %>% mutate(from_to = paste(from, to, sep = "-"), Donor.ID = donor)
  
  # sort to remove reciprocal matches
  match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed_DONOR <- unique(data.frame(t(apply(match_TCR_beta_donor_all_cluster_count_all_DONOR[,c(1,2)], 1, sort))))
  colnames(match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed_DONOR) <- c("from","to") 
  match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed_DONOR <- match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed_DONOR %>% mutate(from_to = paste(from, to, sep = "-"))
  
  # filter original edge df now that reciprocal combinations have been removed
  match_TCR_beta_donor_all_cluster_count_all_DONOR <- match_TCR_beta_donor_all_cluster_count_all_DONOR %>% 
    filter(from_to %in% match_TCR_beta_donor_all_cluster_count_all_reciprocal_removed_DONOR$from_to)
  
  # Create node df: of how many cells in each cluster have a tcr - in each donor individually
  match_TCR_beta_donor_all_cluster_donor_tcr_count_DONOR <- data %>% 
    # filter for donor of interest
    filter(Donor.ID == donor) %>% 
    distinct( Cluster.Name.From, tcr1) %>% 
    dplyr::count( Cluster.Name.From) %>% dplyr::rename(from = Cluster.Name.From, number_cells_with_tcr = n)
  
    ## return output by donor so I can concatenate
  match_TCR_beta_donor_all_cluster_count_all_DONOR
}

TRB_all_sharing_donor <- lapply(donor_list, TRB_all_donor_sharing, match_TCR_beta_donor_all_cluster_rm_dup)
names(TRB_all_sharing_donor ) <- for (donor in donor_list) {
  print(paste0("TRB_all_sharing_", donor))
}
TRB_all_sharing_donor <- bind_rows(TRB_all_sharing_donor)
#save(TRB_all_sharing_donor, file = file.path(resultDir,"TRB_all_sharing_donor.Rdata"))

# TRB
# load previously saved data
load(file = file.path(resultDir,"TRB_all_sharing_donor.Rdata"))
TRB_all_sharing_donor
class(TRB_all_sharing_donor$number_sharing) # integer

# spread the donor and make the from to as the rownames
TRB_all_sharing_donor_heatmap <- TRB_all_sharing_donor %>% ungroup() %>%
  dplyr::select(-c(from,to)) %>%
  pivot_wider(names_from = Donor.ID, values_from = number_sharing) %>% column_to_rownames(.,var = "from_to")

TRB_all_sharing_donor_heatmap[is.na(TRB_all_sharing_donor_heatmap)] <- 0
TRB_all_sharing_donor_heatmap <- as.matrix(TRB_all_sharing_donor_heatmap)

# get heatmap annotation
hc <- columnAnnotation(Response = donor_response_match$Response,
                       col = list(Response = c("R" = "#353a33",
                                               "NR" = "#c6cec4")))

### Separate heatmap into between cluster

# TRA
within_cluster <- c("1-1","2-2","3-3", "4-4","5-5","6-6","7-7","8-8")
TRA_all_sharing_donor_heatmap_between_cluster <- TRA_all_sharing_donor_heatmap[!(rownames(TRA_all_sharing_donor_heatmap) %in% within_cluster),]

# change donors for paper figure to use operational PID and not donor ID
TRA_all_sharing_donor_heatmap_between_cluster_PID <- as.data.frame(t(TRA_all_sharing_donor_heatmap_between_cluster)) %>% 
  rownames_to_column(var = "Donor.ID") %>% 
  left_join(.,T1DAL_ITN_ID_dictionary[,c("Donor.ID","masked_public_PID")]) %>%
  column_to_rownames(var = "masked_public_PID") %>% select(-Donor.ID) 
TRA_all_sharing_donor_heatmap_between_cluster_PID <- t(TRA_all_sharing_donor_heatmap_between_cluster_PID)  

# FIGURE S11B TRA: repeat plot with larger text size and using public donor ID
pdf(file = file.path(plotDir,"TRA_all_sharing_donor_heatmap_between_cluster_PID.pdf"),height = 8, width = 8)
ComplexHeatmap::Heatmap(TRA_all_sharing_donor_heatmap_between_cluster_PID, 
                        cluster_rows = FALSE, 
                        name = "Total TRA Sharing\nBetween Cluster",
                        top_annotation = hc,
                        col =plasma(8, alpha  = 1, begin = 0, end = 1, direction = 1),
                        # increase text size
                        column_names_gp = grid::gpar(fontsize = 16), row_names_gp = grid::gpar(fontsize = 16), 
                        heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12), title_gp = grid::gpar(fontsize = 16)))
dev.off()


# TRB
TRB_all_sharing_donor_heatmap_between_cluster <- TRB_all_sharing_donor_heatmap[!(rownames(TRB_all_sharing_donor_heatmap) %in% within_cluster),]

# change donors for paper figure to use operational PID and not donor ID
TRB_all_sharing_donor_heatmap_between_cluster_PID <- as.data.frame(t(TRB_all_sharing_donor_heatmap_between_cluster)) %>% 
  rownames_to_column(var = "Donor.ID") %>% 
  left_join(.,T1DAL_ITN_ID_dictionary[,c("Donor.ID","masked_public_PID")]) %>%
  column_to_rownames(var = "masked_public_PID") %>% select(-Donor.ID) 
TRB_all_sharing_donor_heatmap_between_cluster_PID <- t(TRB_all_sharing_donor_heatmap_between_cluster_PID)  

# FIGURE S11B TRB: repeat plot with larger text size and using public donor ID
pdf(file = file.path(plotDir,"TRB_all_sharing_donor_heatmap_between_cluster_PID.pdf"),height = 8, width = 8)
ComplexHeatmap::Heatmap(TRB_all_sharing_donor_heatmap_between_cluster_PID, 
                        cluster_rows = FALSE, 
                        name = "Total TRB Sharing\nBetween Cluster",
                        top_annotation = hc,
                        #col=brewer.pal(9,"Blues"))
                        col =plasma(8, alpha  = 1, begin = 0, end = 1, direction = 1) ,
                        # increase text size
                        column_names_gp = grid::gpar(fontsize = 16), row_names_gp = grid::gpar(fontsize = 16), 
                        heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12), title_gp = grid::gpar(fontsize = 16)))

dev.off()

#### Recategorize TRA sharing as different differentiation trajectories ####

# load previously saved data
load(file= file.path(resultDir, "TCR_UMAP_circos_data.Rdata"))

### Set up binary matrix that will be used to recategorize the cluster sharing patterns below!

# start with original data used to look at the TCR connections between clusters
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos

# filter the data just to make it unique by the TCRs found in each cluster - regardless of donor
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos %>%
  distinct(cdr3_nt, Cluster.Name)

# make A binary matrix/data frame where rows are elements and columns are sets and data is filled in as 1 and zero
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster%>%
  # add counter
  mutate(Test = 1) %>%
  tidyr::pivot_wider(values_from = Test, names_from = Cluster.Name) %>% 
  replace(is.na(.), 0) %>% column_to_rownames(.,var = "cdr3_nt")

tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary[,c("1","2","3","4","5","6","7","8")]
head(tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary )

### Recategorize the clusters
# start at DN cluster 4 since there is little sharing with the first few clusters

## Tex PD1 trajectory = sharing with cluster 4 - 7 ONLY
## Tex CD7 trajectory = sharing with 4-6, 4-8, 4-5 ONLY OR a combination of any two or three of these 
## Tex Branching = sharing with both 4-7, and any or some of 4-6, 4-8, 4-5
## Tex Fluid = sharing of UNIQUE TCRs not shared with 4 between 5-7,6-7, 8-7 or multiple of these 

# start with the binary matrix generated to create the upset plot above 
tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary

# get all combinations from 4 onward
tra_combos_4_after <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary %>% rownames_to_column("cdr3_nt_TRA") %>%
  dplyr::select(-1,-2,-3,-4) %>%
  dplyr::distinct()

write.csv(tra_combos_4_after, file = file.path(resultDir, "tra_combos_4_after.csv"))

# get all combinations
tra_combos_all <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary %>% rownames_to_column("cdr3_nt_TRA") %>%
  dplyr::select(-1) %>%
  dplyr::distinct()

write.csv(tra_combos_all, file = file.path(resultDir,"tra_combos_all.csv"))

# load binary trajectory matrix I generated manually that categorizes trajectories
binary_trajectories = readxl::read_xlsx(file.path(annotationDir, "Binary_trajectory_options.xlsx"))

# join with binary matrix for 4 onward
tra_cluster_binary_trajectories <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary %>% 
  rownames_to_column("cdr3_nt_TRA") %>% dplyr::select("cdr3_nt_TRA", "4","5","6","7","8") %>%
  left_join(., binary_trajectories) 

# find distinct categories after cluster 4 that I missed
tra_cluster_binary_trajectories_distinct <- tra_cluster_binary_trajectories %>%
  dplyr::select(2:7) %>%
  dplyr::distinct()

# export version so I can fill in my own personal gaps in the categories
write.csv(tra_cluster_binary_trajectories_distinct , file = file.path(resultDir,"tra_cluster_binary_trajectories_distinct.csv"))

# import edited version with the remaining categories filled in
tra_full_trajectory <- read_csv(file = file.path(annotationDir,"tra_cluster_binary_trajectories_distinct_edited.csv")) %>% 
  # remove rownames
  dplyr::select(-1)

# join back with cells 
tra_cluster_binary_trajectories_4_onward <- tcr_v_j_cdr3_comb_all_TRA_full_seq_filtered_UMAP_circos_unique_cluster_binary %>% 
  rownames_to_column("cdr3_nt_TRA") %>% dplyr::select("cdr3_nt_TRA", "4","5","6","7","8") %>%
  left_join(., tra_full_trajectory) 
tra_cluster_binary_trajectories_4_onward %>% filter(is.na(Trajectory)) # none, all TRA clones have been categorized

# save for future reference
#save(tra_cluster_binary_trajectories_4_onward, file = file.path(resultDir,"tra_cluster_binary_trajectories_4_onward.Rdata"))
nrow(tra_cluster_binary_trajectories_4_onward) # 2906

#### FIGURE S13B: Recategorized TRA statistics ####

load(file.path(resultDir, "tra_cluster_binary_trajectories_4_onward.Rdata"))
load(file =file.path(resultDir,"cds_tcr_manual_metadata.RData"))

# Add up the number of TRAs in each category
tra_cluster_binary_trajectories_4_onward_stats <- tra_cluster_binary_trajectories_4_onward %>% dplyr::count(Trajectory)

## how many per donor? 
tra_cluster_binary_trajectories_4_onward_stats_donor <- left_join(tra_cluster_binary_trajectories_4_onward, 
                                                                  unique(cds_tcr_manual_metadata_ordered[ ,c("Donor.ID","cdr3_nt_TRA")] )) %>%
  dplyr::count(Trajectory, Donor.ID) %>%
  filter(Trajectory %in% c("Tex_Branching", "Tex_CD57" ,"Tex_Fluid","Tex_PD1" ))

# repeat but don't condense by counts
tra_cluster_binary_trajectories_4_onward_stats_all_donor <- left_join(tra_cluster_binary_trajectories_4_onward, 
                                                                      unique(cds_tcr_manual_metadata_ordered[ ,c("Donor.ID","cdr3_nt_TRA")] )) %>% 
  filter(Trajectory %in% c("Tex_Branching", "Tex_CD57" ,"Tex_Fluid","Tex_PD1" ))

tra_cluster_binary_trajectories_4_onward_stats_donor$Trajectory <- factor(tra_cluster_binary_trajectories_4_onward_stats_donor$Trajectory,
                                                                          levels = c("Tex_Branching", "Tex_CD57","Tex_Fluid","Tex_PD1"  ),
                                                                          labels = c("Tex-Branching", "Tex-CD57+","Tex-Fluid","Tex-PD1")) 

# save for future use
#save(tra_cluster_binary_trajectories_4_onward_stats_all_donor, file = file.path(resultDir,"tra_cluster_binary_trajectories_4_onward_stats_all_donor.Rdata"))

# Join with UMAP data to find where these cells are!
cds_tcr_manual_metadata_ordered_TRA_trajectory <- left_join(cds_tcr_manual_metadata_ordered, tra_cluster_binary_trajectories_4_onward)
# Joining, by = "cdr3_nt_TRA"

# save for use downstream 
#save(cds_tcr_manual_metadata_ordered_TRA_trajectory, file = file.path(resultDir, "cds_tcr_manual_metadata_ordered_TRA_trajectory.RData"))

# do all cells with TCR have a trajectory label?
cds_tcr_manual_metadata_ordered_TRA_trajectory  %>% filter(!is.na(clonotype_id)) %>% filter(is.na(Trajectory)) # 0 rows, all have a trajectory

# keep only cells with TCR 
cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR <- cds_tcr_manual_metadata_ordered_TRA_trajectory %>% filter(!is.na(clonotype_id))
nrow(cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR ) # 4420

# plot number of cells across all donor trajectory types of interest - each dot is number of cells per clone
cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj <- cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR %>%
  dplyr::count(cdr3_nt_TRA,Trajectory) %>%
  filter(Trajectory %in% c("Tex_Branching", "Tex_CD57" ,"Tex_Fluid","Tex_PD1" )) 

cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj$Trajectory <- factor(cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj$Trajectory,
                                                                                   levels = c("Tex_Branching", "Tex_CD57","Tex_Fluid","Tex_PD1"  ),
                                                                                   labels = c("Tex-Branching", "Tex-CD57+","Tex-Fluid","Tex-PD1")) 
# get stats to plot
cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj_wilcox <- 
  cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj %>%
  pairwise_wilcox_test(n ~ Trajectory) %>%
  add_y_position(step.increase = 0.02) %>% 
  filter(p.adj <= 0.05)

## FIGURE S13B
cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj_plot <- cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj %>%
  ggplot(aes(x = Trajectory, y = n)) +
  geom_jitter(width = 0.2) + 
  geom_violin(alpha = 0.1) +
  stat_pvalue_manual(cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj_wilcox,
                     tip.length = 0,step.increase = 0.02) +
  labs(x = "Trajectory", y = "Total Cells per Clone")
# FIGURE S13B 
ggsave(cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj_plot, file = file.path(plotDir,"cds_tcr_manual_metadata_ordered_TRA_trajectory_TCR_clone_traj_plot.pdf"),
       width = 7, height = 5)

#### FIGURE S13A: Calculate TRA diversity by donor and trajectory ####

load(file = file.path(resultDir,"tra_cluster_binary_trajectories_4_onward_stats_all_donor.Rdata"))

tra_cluster_binary_trajectories_4_onward_stats_all_donor
donor_list <- unique(tra_cluster_binary_trajectories_4_onward_stats_all_donor$Donor.ID)

shannon_diversity <- function(x, i) {
  tra_cluster_binary_trajectories_4_onward_stats_all_donor %>% filter(Donor.ID == x & Trajectory == i) %>% 
    group_by(cdr3_nt_TRA) %>%
    mutate(number = cur_group_id()) %>%
    ungroup() %>%
    summarize(diversity = vegan::diversity(number)) 
}

Tex_CD57_donor <- lapply(donor_list,shannon_diversity, "Tex_CD57" )
Tex_CD57_donor <- bind_rows(Tex_CD57_donor)
Tex_CD57_donor$Donor.ID <- donor_list
Tex_CD57_donor$Trajectory <- "Tex_CD57"

Tex_PD1_donor <- lapply(donor_list,shannon_diversity, "Tex_PD1" )
Tex_PD1_donor <- bind_rows(Tex_PD1_donor)
Tex_PD1_donor$Donor.ID <- donor_list
Tex_PD1_donor$Trajectory <- "Tex_PD1"

Tex_Branching_donor <- lapply(donor_list,shannon_diversity, "Tex_Branching" )
Tex_Branching_donor <- bind_rows(Tex_Branching_donor)
Tex_Branching_donor$Donor.ID <- donor_list
Tex_Branching_donor$Trajectory <- "Tex_Branching"

Tex_Fluid_donor <- lapply(donor_list,shannon_diversity, "Tex_Fluid" )
Tex_Fluid_donor <- bind_rows(Tex_Fluid_donor)
Tex_Fluid_donor$Donor.ID <- donor_list
Tex_Fluid_donor$Trajectory <- "Tex_Fluid"

Tex_shannon <- rbind(Tex_CD57_donor, Tex_PD1_donor,Tex_Branching_donor, Tex_Fluid_donor) %>% 
  filter(diversity !=0)

Tex_shannon$Trajectory <- factor(Tex_shannon$Trajectory,
                                 levels = c("Tex_Branching", "Tex_CD57","Tex_Fluid","Tex_PD1"  ),
                                 labels = c("Tex-Branching", "Tex-CD57+","Tex-Fluid","Tex-PD1")) 

# convert operational ID to public ID
Tex_shannon_PID <- left_join(Tex_shannon , T1DAL_ITN_ID_dictionary)

# plot diversity
Tex_shannon_plot <- 
  ggplot(Tex_shannon_PID, aes(x = Trajectory, y = diversity)) + 
  geom_jitter(width = 0.2, aes(color = masked_public_PID),   size = 3) + 
  geom_violin(alpha = 0.1) +
  labs(x = "Trajectory", y = "Shannon-Wiener Index", colour= "Donor", title = "Shannon-Wiener Index") +
  stat_compare_means(method = "anova", size = 6) + theme(text = element_text(size = 18),
                                                         axis.text.x = element_text(angle = 70, hjust = 1)
  )
## FIGURE S13A
ggsave(Tex_shannon_plot , file = file.path(plotDir,"Tex_shannon_plot.pdf"),height = 7, width = 7)

#### FIGURE 6A: Upset plot with recategorized TRA ####

load(file.path(resultDir, "tra_cluster_binary_trajectories_4_onward.Rdata"))

# Plot upset plot using distinct mode (default)
levels(as.factor(tra_cluster_binary_trajectories_4_onward$Trajectory))

tra_cluster_binary_trajectories_4_onward_upset <- tra_cluster_binary_trajectories_4_onward %>% 
  filter(Trajectory %in% c("Tex_Branching", "Tex_CD57" ,"Tex_Fluid","Tex_PD1" )) %>%
  column_to_rownames("cdr3_nt_TRA") %>% select(-Trajectory)
colnames(tra_cluster_binary_trajectories_4_onward_upset)
tra_cluster_binary_trajectories_4_onward_upset_m <- as.matrix(tra_cluster_binary_trajectories_4_onward_upset)
m <- make_comb_mat(tra_cluster_binary_trajectories_4_onward_upset_m)

# get correct grouping for row order
degree_join <- data.frame(size = comb_degree(m)) %>% 
  rownames_to_column(var = "comb_degree")

tra_full_trajectory_order  <- 
  unite(tra_full_trajectory, comb_degree, -c("Trajectory"), sep = "") 

tra_full_trajectory_order_join <- left_join(tra_full_trajectory_order, degree_join) %>%
  filter(!is.na(size)) %>% arrange(Trajectory)

# set order to look like format of comb_degree(m)
traj_order <- setNames( tra_full_trajectory_order_join $size,tra_full_trajectory_order_join $comb_degree)

# additional annotations were added after the fact in inkscape
## FIGURE 6A
pdf(file.path(plotDir,"TRA_recategorize_trajectory_upset.pdf"), height = 5, width = 5)
UpSet(t(m),
      comb_order = c(9,1,2,7,4,5,10,6,12,11,14,8,3,15,17,16,13),
      top_annotation = upset_top_annotation(t(m), 
                                            show_annotation_name = FALSE, axis_param=list(gp=gpar(fontsize = 14)),
                                            numbers_param =list(gp=gpar(fontsize = 14))), 
      right_annotation = upset_right_annotation(t(m),  
                                                show_annotation_name = FALSE, axis_param=list(gp=gpar(fontsize = 14))),
      column_names_gp = grid::gpar(fontsize = 16),
      row_names_gp = grid::gpar(fontsize = 16),
      row_title_gp = gpar(fontsize = 16))


dev.off()


#### FIGURE 6B: plot UMAP with different clone trajectories ####

load(file =file.path(resultDir,"cds_tcr_manual_metadata.RData"))
load(file = file.path(resultDir, "cds_tcr_manual_metadata_ordered_TRA_trajectory.RData"))

## add trajectory info to metadata - first check order
all(cds_tcr_manual_metadata_ordered$barcode_original == cds_tcr_manual_metadata_ordered_TRA_trajectory$barcode_original) # TRUE in the same order

# create separate column for each trajectory to filter
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP <- cds_tcr_manual_metadata_ordered_TRA_trajectory %>% 
  mutate(Tex_branching = case_when(Trajectory == "Tex_Branching"~"2",Trajectory != "Tex_branching"~"1", is.na(Trajectory) ~ "0")) %>% 
  mutate(Tex_fluid = case_when(Trajectory == "Tex_Fluid"~"2",Trajectory != "Tex_Fluid"~"1", is.na(Trajectory) ~ "0")) %>%
  mutate(Tex_CD57 = case_when(Trajectory == "Tex_CD57"~"2",Trajectory != "Tex_CD57"~"1", is.na(Trajectory) ~ "0")) %>%
  mutate(Tex_PD1 = case_when(Trajectory == "Tex_PD1"~"2",Trajectory != "Tex_PD1"~"1", is.na(Trajectory) ~ "0")) 

# add trajectory info
cds_all_tcr_manual_ordered_trajectory <- cds_all_tcr_manual_ordered
colData(cds_all_tcr_manual_ordered_trajectory)$Trajectory <-cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP$ Trajectory
colData(cds_all_tcr_manual_ordered_trajectory)$Tex_branching <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP$ Tex_branching
colData(cds_all_tcr_manual_ordered_trajectory)$Tex_fluid <-cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP$ Tex_fluid
colData(cds_all_tcr_manual_ordered_trajectory)$Tex_CD57 <-cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP$ Tex_CD57
colData(cds_all_tcr_manual_ordered_trajectory)$Tex_PD1 <-cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP$ Tex_PD1

# plot individual trajectories
cds_tcr_manual_metadata_ordered_TRA_trajectory_branching_order <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP %>%
  arrange(Tex_branching)
cds_all_tcr_manual_ordered_trajectory_tex_branching <- plot_cells(cds_all_tcr_manual_ordered_trajectory[,cds_tcr_manual_metadata_ordered_TRA_trajectory_branching_order$barcode_original] ,color_cells_by="Tex_branching",show_trajectory_graph = F, 
                                                                  label_cell_groups = FALSE, cell_size=1) +
  scale_color_manual(values  = c("gray80", "gray80","#ab62c0")) +
  theme(text = element_text(size = 10), strip.text.x = element_text(size = 10)) +
  guides(colour = "none") +
  ggtitle("Tex-Branching Cells")

cds_tcr_manual_metadata_ordered_TRA_trajectory_fluid_order <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP %>%
  arrange(Tex_fluid)
cds_all_tcr_manual_ordered_trajectory_tex_fluid <- plot_cells(cds_all_tcr_manual_ordered_trajectory[,cds_tcr_manual_metadata_ordered_TRA_trajectory_fluid_order$barcode_original] ,
                                                              color_cells_by="Tex_fluid",show_trajectory_graph = F, cell_size=1, label_cell_groups = FALSE) +
  scale_color_manual(values  = c("gray80", "gray80","#c2843c")) +
  theme(text = element_text(size = 10), strip.text.x = element_text(size = 10)) +
  guides(colour = "none") +
  ggtitle("Tex-Fluid Cells")

cds_tcr_manual_metadata_ordered_TRA_trajectory_PD1_order <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP %>%
  arrange(Tex_PD1)
cds_all_tcr_manual_ordered_trajectory_tex_pd1 <- plot_cells(cds_all_tcr_manual_ordered_trajectory[,cds_tcr_manual_metadata_ordered_TRA_trajectory_PD1_order$barcode_original] ,color_cells_by="Tex_PD1",
                                                            show_trajectory_graph = F, cell_size=1, label_cell_groups = FALSE) +
  scale_color_manual(values  = c("gray80","gray80" ,"#ba4d4cc7")) +
  theme(text = element_text(size = 10), strip.text.x = element_text(size = 10)) +
  guides(colour = "none") +
  ggtitle("Tex-PD-1+ Cells")

cds_tcr_manual_metadata_ordered_TRA_trajectory_CD57_order <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP %>%
  arrange(Tex_CD57)
cds_all_tcr_manual_ordered_trajectory_tex_cd57 <- plot_cells(cds_all_tcr_manual_ordered_trajectory[,cds_tcr_manual_metadata_ordered_TRA_trajectory_CD57_order$barcode_original] ,color_cells_by="Tex_CD57",
                                                             show_trajectory_graph = F, cell_size=1, label_cell_groups = FALSE) +
  scale_color_manual(values  = c("gray80","gray80", "#93a24eff"),) +
  theme(text = element_text(size = 10), strip.text.x = element_text(size = 10)) +
  guides(colour = "none") +
  ggtitle("Tex-CD57+ Cells")

# FIGURE 6B: export all umaps
all_trajectory_umaps <- ggarrange(cds_all_tcr_manual_ordered_trajectory_tex_branching,cds_all_tcr_manual_ordered_trajectory_tex_cd57 ,
                                  cds_all_tcr_manual_ordered_trajectory_tex_fluid,
                                  cds_all_tcr_manual_ordered_trajectory_tex_pd1 )

ggsave(all_trajectory_umaps, file = file.path(plotDir, "all_trajectory_umaps.pdf"), height = 5, width = 5)

#### Assess tcr Ag specificity of shared TRAs with any cluster ####

load(file= file.path(resultDir, "match_TCR_alpha_donor_all_cluster_filter.Rdata"))
load(file = file.path(resultDir,"all_libs_tcrs.Rdata"))

# get df with matched tcr alpha chains shared between any cluster
match_TCR_alpha_donor_all_cluster_filter

# get unique shared alpha chains per donor
match_TCR_alpha_donor_all_cluster_filter_donor_distinct <- match_TCR_alpha_donor_all_cluster_filter %>%
  distinct(Donor.ID, tcr1, v_gene,j_gene,cdr3_nt) %>% dplyr::rename(libid = tcr1)
nrow(match_TCR_alpha_donor_all_cluster_filter_donor_distinct) # 1479

# match tcr1 and cell barcode with the full sequence information for TRB and TRA
all_libs_tcrs_cluster_TRA_TRB <- all_libs_tcrs %>% filter(libid  %in% match_TCR_alpha_donor_all_cluster_filter_donor_distinct$libid)
nrow(all_libs_tcrs_cluster_TRA_TRB) # 2958

# Load VDJdb database from txt file
vdjdb <- read_table(file = file.path(annotationDir,"vdjdb-2022-03-30/vdjdb.txt"))

# Remove allele info from v and j gene name for matching with our table, and rename gene to chain
vdjdb_v_gene_j_gene <- vdjdb %>% separate(v.segm, into =  c("v_gene","v_allele_group"), sep = "\\*") %>% 
  separate(j.segm, into =  c("j_gene","j_allele_group"), sep = "\\*") %>% dplyr::rename(chain = gene)

# Match TRA and TRB with their specificities
all_libs_tcrs_cluster_TRA_TRB_specificity <- left_join(all_libs_tcrs_cluster_TRA_TRB, vdjdb_v_gene_j_gene, by = c("v_gene","j_gene","cdr3","chain"))
nrow(all_libs_tcrs_cluster_TRA_TRB_specificity) # 4863

## Count the number of matches per TRA or TRB
all_libs_tcrs_cluster_TRA_TRB_specificity %>% filter(is.na(antigen.species)) %>% nrow() # 2767 of 4863 have no match if you search for v, j and cdr3
all_libs_tcrs_cluster_TRA_TRB_specificity %>% filter(chain == "TRA" & !is.na(antigen.species)) %>% dplyr::count(libid) %>% arrange(desc(n)) 
# 129 cells have a TRA hit,  60 cells have 16 hits, the rest have 1 or 2
all_libs_tcrs_cluster_TRA_TRB_specificity %>% filter(chain == "TRB"& !is.na(antigen.species)) %>% dplyr::count(libid) %>% arrange(desc(n)) 
# 62 total TRB, 59 cells have 18 hits, the rest have 1

# keep df of those with hits 
all_libs_tcrs_cluster_TRA_TRB_specificity_filtered <- all_libs_tcrs_cluster_TRA_TRB_specificity %>% filter(!is.na(antigen.species))

#save(all_libs_tcrs_cluster_TRA_TRB_specificity_filtered , file = file.path(resultDir, "all_libs_tcrs_cluster_TRA_TRB_specificity_filtered.Rdata"))

### Categorize distinct TCRs into public and private and cross reactive

# find all distinct matches regardless of how many times they show up in the database
all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct <- all_libs_tcrs_cluster_TRA_TRB_specificity_filtered  %>% 
  dplyr::distinct(Donor.ID, cdr3, chain, antigen.species, antigen.epitope, antigen.gene, antigen.species, cdr3_nt ) 

# find public and private 
all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct_pub_priv <- all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct %>%
  distinct(cdr3, Donor.ID) %>%  group_by(cdr3) %>% mutate(public_private = case_when(n() >1 ~ "public", n() == 1 ~ "private"))

# find cross reactive
all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct_ag_species <- all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct %>%
  distinct(cdr3, antigen.species) %>% dplyr::count(cdr3) %>% dplyr::rename(number_ag_species = n) 

# rejoin grouping information 
all_libs_tcrs_cluster_TRA_TRB_specificity_pub_priv_ag_species <- left_join(all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct, 
                                                                           all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct_pub_priv ) %>%
  left_join(all_libs_tcrs_cluster_TRA_TRB_specificity_filtered_distinct_ag_species)

# split into TRA and TRB
all_libs_tcrs_cluster_pub_priv_ag_species_TRA <- all_libs_tcrs_cluster_TRA_TRB_specificity_pub_priv_ag_species %>% filter(chain == "TRA") %>%
  # rename cdr3_nt column for rejoining with cds object
  dplyr::rename(cdr3_nt_TRA = cdr3_nt)

all_libs_tcrs_cluster_pub_priv_ag_species_TRB <- all_libs_tcrs_cluster_TRA_TRB_specificity_pub_priv_ag_species %>% filter(chain == "TRB") %>%
  # rename cdr3_nt column for rejoining with cds object
  dplyr::rename(cdr3_nt_TRB = cdr3_nt)

# save these objects
#save(all_libs_tcrs_cluster_pub_priv_ag_species_TRA, file = file.path(resultDir, "all_libs_tcrs_cluster_pub_priv_ag_species_TRA.Rdata"))
#save(all_libs_tcrs_cluster_pub_priv_ag_species_TRB, file = file.path(resultDir, "all_libs_tcrs_cluster_pub_priv_ag_species_TRB.Rdata"))

#### Determine trajectories on a donor level for CMV specific TRAs ####

# Load previous TRA object
load(file = file.path(resultDir,"all_libs_tcrs_cluster_pub_priv_ag_species_TRA.Rdata"))
load(file.path(resultDir, "tra_cluster_binary_trajectories_4_onward.Rdata"))
load(file.path(resultDir, "cds_all_tcr_manual_ordered.Rdata"))

# subset for CMV specific cdr3
all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV <- all_libs_tcrs_cluster_pub_priv_ag_species_TRA %>% filter(antigen.species == "CMV")

# Join with differentiation trajectory for each cell
all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV_traj_4_onward <- left_join(all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV, tra_cluster_binary_trajectories_4_onward)

# join with umap to get number of cells with each
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP %>% 
  filter(cdr3_nt_TRA %in% all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV$cdr3_nt_TRA)
# save for plotting airline plot later
#save(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv, file = file.path(resultDir,"cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv.Rdata"))

# calculate the number of cells with each trajectory per donor
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_donor <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv  %>% filter(Trajectory != "Only_4") %>% dplyr::count(Donor.ID, Trajectory, Response)

# total tcrs per person from those cells that are in the UMAP
total_tcrs_per_person_in_UMAP <- cds_tcr_manual_metadata_ordered %>%  filter(!is.na(clonotype_id)) %>% dplyr::count(Donor.ID)

# add public donor ID for plotting
total_tcrs_per_person_in_UMAP <- left_join(total_tcrs_per_person_in_UMAP , T1DAL_ITN_ID_dictionary)

# Normalize values by percent of the total cells with TRA for each donor
total_tcrs_per_person_in_UMAP_cmv <- total_tcrs_per_person_in_UMAP %>% dplyr::rename(total_UMAP_TCR = n)
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_donor <- left_join(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_donor, total_tcrs_per_person_in_UMAP_cmv) %>%
  mutate(percent_of_UMAP_TCR = n/total_UMAP_TCR *100)

#### Determine donor level trajectories for CMV AND EBV specific TRAs ####

# Load previous TRA object
load(file = file.path(resultDir, "all_libs_tcrs_cluster_pub_priv_ag_species_TRA.Rdata"))

# subset for CMV AND EBV specific cdr3
all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV_EBV <- all_libs_tcrs_cluster_pub_priv_ag_species_TRA %>% 
  filter(antigen.species == "CMV" | antigen.species == "EBV")

# Join with differentiation trajectory for each cell
all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV_EBV_traj_4_onward <- left_join(all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV_EBV, tra_cluster_binary_trajectories_4_onward)

# join with umap to get number of cells with each
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP %>% 
  filter(cdr3_nt_TRA %in% all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV_EBV$cdr3_nt_TRA)

# save for use in airline plot
#save(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv, file = file.path(resultDir,"cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv.Rdata"))

cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv[!(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv$cdr3_nt_TRA %in% unique(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv$cdr3_nt_TRA )),]
# this last code revealed that only two additional cdr3_nt_TRAs were identified that were only EBV specific, and these cells had the CD57_no_4 trajectory and I remove them in 
# the next piece of code because we are not interested in cells with this trajectory! 

# calculate the number of cells with each trajectory per donor - EBV and CMV combined 
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv  %>%
  filter(!(Trajectory %in% c("Only_4", "CD57_no_4"))) %>%
  dplyr::count(Donor.ID, Trajectory, Response)
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_donor

#Normalize values by percent of the total cells with TRA for each donor
total_tcrs_per_person_in_UMAP_cmv_ebv <- total_tcrs_per_person_in_UMAP %>% dplyr::rename(total_UMAP_TCR = n)
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor <- left_join(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor, total_tcrs_per_person_in_UMAP_cmv_ebv) %>%
  mutate(percent_of_UMAP_TCR = n/total_UMAP_TCR *100)

#save(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor, file = file.path(resultDir, "cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor.Rdata"))

#### FIGURE S14C: Combine CMV/EBV specific scRNAseq data with CYTOF panel data from Alice ####

load(file = file.path(resultDir, "cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor.Rdata"))

# Load reformatted version of Alice Wiedeman CYTOF data with the percent of TMR specific cells for each donor across TN10 and T1DAL
TMR_specific_CYTOF_total_exhausted <- readxl::read_xlsx(file.path(annotationDir,"AAI_2022_Oral_Presentation_2022-05-07_Plus_Virus_for_Erin_AEW_2023-01-06_RAW_DATA.pptx_EW_reformatted.xlsx"),
                                                        sheet = 3)
TMR_specific_CYTOF_PD1_exhausted <- readxl::read_xlsx(file.path(annotationDir,"AAI_2022_Oral_Presentation_2022-05-07_Plus_Virus_for_Erin_AEW_2023-01-06_RAW_DATA.pptx_EW_reformatted.xlsx"),
                                                      sheet = 4)
TMR_specific_CYTOF_CD57_exhausted <- readxl::read_xlsx(file.path(annotationDir,"AAI_2022_Oral_Presentation_2022-05-07_Plus_Virus_for_Erin_AEW_2023-01-06_RAW_DATA.pptx_EW_reformatted.xlsx"),
                                                       sheet = 5)
# load and format CYTOF key
CYTOF_donor_ID_key <-  readxl::read_xlsx(file.path(annotationDir,"AAI_2022_Oral_Presentation_2022-05-07_Plus_Virus_for_Erin_AEW_2023-01-06_RAW_DATA.pptx_EW_reformatted.xlsx"),
                                         sheet = 2) %>% dplyr::rename(Subject_ID =`Subject ID`, Donor.ID = ParticipantID)
CYTOF_donor_ID_key$Subject_ID <- as.character(CYTOF_donor_ID_key$Subject_ID)
CYTOF_donor_ID_key$Subject_ID <- paste0("SubjT", CYTOF_donor_ID_key$Subject_ID) # add SubjT header to match with other spreadsheets

# Turn turn all sheets into long format and combine
TMR_specific_CYTOF_total_exhausted <- TMR_specific_CYTOF_total_exhausted %>% 
  pivot_longer(cols = c(4:15), names_to = c("Specificity", "Visit"), values_to = "Percent", names_sep = "_") 
TMR_specific_CYTOF_PD1_exhausted <- TMR_specific_CYTOF_PD1_exhausted %>% 
  pivot_longer(cols = c(4:15), names_to = c("Specificity", "Visit"), values_to = "Percent", names_sep = "_") 
TMR_specific_CYTOF_CD57_exhausted <- TMR_specific_CYTOF_CD57_exhausted %>% 
  pivot_longer(cols = c(4:15), names_to = c("Specificity", "Visit"), values_to = "Percent", names_sep = "_") 

TMR_specific_CYTOF_all_pop <- rbind(TMR_specific_CYTOF_total_exhausted, TMR_specific_CYTOF_PD1_exhausted, TMR_specific_CYTOF_CD57_exhausted)

# decode all subject IDs to match donor IDs
TMR_specific_CYTOF_all_pop <- left_join(TMR_specific_CYTOF_all_pop, CYTOF_donor_ID_key) 

# Filter data for CMV/EBV specific ("Virus")
TMR_specific_CYTOF_all_pop_virus <- TMR_specific_CYTOF_all_pop %>% filter(Specificity == "Virus") %>%
  # and remove TN10 patients that don't have a listed donor ID
  filter(!is.na(Donor.ID))
TMR_specific_CYTOF_all_pop_virus$Donor.ID <- as.character(TMR_specific_CYTOF_all_pop_virus$Donor.ID)

# Join with the percent of cells that are EBV or CMV specific with the different trajectories 
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus <- left_join(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor, TMR_specific_CYTOF_all_pop_virus)

# add weeks as column for plotting rather than visit for CYTOF
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus %>%
  mutate(CYTOF_weeks = case_when(Visit == "Visit1" ~ "0", Visit  == "Visit2" ~ "35", Visit == "Visit3"~"104"))
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus$CYTOF_weeks <- as.numeric(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus$CYTOF_weeks)

# make donor factor
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus$Donor.ID <- factor(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus$Donor.ID,
                                                                                                 levels = c( "10289","10469","10748","10256","10295","10396","10458"))
# set trajectory colors
#PD1: "#ba4d4cc7" CD57: "#93a24eff" Tex branching "#ab62c0", tex fluid 
trajectory_colors <- c("#ab62c0", "#93a24eff", "#c2843c","#ba4d4cc7")

# Plot percent EMV and EBV specific PD1 and CD57 data stacked for each donor for scRNAseq then CYTOF data
CMV_EBV_traj_TCR <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus %>%
  distinct(Donor.ID, Response, Trajectory, percent_of_UMAP_TCR)

# Fix the donor IDs to be public IDs
CMV_EBV_traj_TCR <- left_join(CMV_EBV_traj_TCR, T1DAL_ITN_ID_dictionary)

# FIGURE S14C
CMV_EBV_TMR_CYTOF <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv_ebv_donor_CYTOF_virus %>%
  filter(Phenotype_group != "Total_exhausted") %>%
  # join with masked donor ID
  left_join(., T1DAL_ITN_ID_dictionary) %>%
  ggplot( aes(x = CYTOF_weeks, y = Percent,  fill = Phenotype_group)) + 
  #geom_jitter(alpha = 1) + 
  #geom_boxplot(outlier.shape = NA) +
  geom_col(position = "dodge") +
  #geom_point() + 
  #geom_line() +
  theme(text = element_text(size = 20)) +
  labs(x = "Visit (Wk)", y = "% CMV/EBV in CYTOF Tex Cluster", fill = "Trajectory") +  facet_grid( .~masked_public_PID + Response) +
  #scale_color_manual(values = c("#93a24eff","#ba4d4cc7")) +
  scale_fill_manual(values = c("#93a24eff","#ba4d4cc7"), breaks = c("CD57_exhausted","PD1_exhausted"), labels = c("Tex-CD57+","Tex-PD-1+")) +
  scale_x_continuous( breaks = c(0,35,104)) + scale_y_continuous(limits = c(0,100))

## FIGURE S14C
ggsave(CMV_EBV_TMR_CYTOF, file = file.path(plotDir,"CMV_EBV_TMR_CYTOF.pdf"), device = "pdf", height = 6, width = 10)

#### FIGURE 6D: Airline plot of CMV Specific TRAs on UMAP ####

# start with object without cluster 9
load(file.path(resultDir,"cds_no_MAIT_no_9.Rdata"))

# load subset cds metadata with the cmv specific cells
load(file = file.path(resultDir,"cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv.Rdata"))

# Use the UMAP coldata where the TCR information has been joined

# make a copy of the annotation plus UMAP coordinates, for easier manipulation
data.tmp_manual <-
  as.data.frame(colData(cds_all_tcr_manual)) %>%
  # Combine the colData with the reduced dimension representation coordinated for these clones
  cbind(
    as.data.frame(reducedDims(cds_all_tcr_manual)$UMAP) %>%
      magrittr::set_colnames(c("V1", "V2"))) %>% 
  dplyr::rename("clone_id_tcr_graph_clonal_expansion" = "tcr_id_TRA") %>%
  filter(barcode_original %in% cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv$barcode_original)

# create data frame to store links
curves.tmp_manual <-
  data.frame(
    clone_id_tcr_graph_clonal_expansion = character(),
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric())

# loop over each clone, and extract UMAP coordinates for cells from the same clone
for (clone_id.tmp_manual in na.omit(unique(data.tmp_manual$clone_id_tcr_graph_clonal_expansion))) {
  clone_id_curves.tmp_manual <- curves.tmp_manual[0,]
  data_for_curves.tmp_manual <-
    data.tmp_manual %>%
    dplyr::filter(clone_id_tcr_graph_clonal_expansion %in% clone_id.tmp_manual)
  if (nrow(data_for_curves.tmp_manual) > 1) {
    for (i in 1:(nrow(data_for_curves.tmp_manual)-1)) {
      for (j in (i+1):nrow(data_for_curves.tmp_manual)) {
        clone_id_curves.tmp_manual <-
          rbind(
            clone_id_curves.tmp_manual,
            list(
              clone_id_tcr_graph_clonal_expansion =
                data_for_curves.tmp_manual$clone_id_tcr_graph_clonal_expansion[i],
              x = data_for_curves.tmp_manual$V1[i],
              y = data_for_curves.tmp_manual$V2[i],
              cluster1 = data_for_curves.tmp_manual$Cluster.Name[i],
              xend = data_for_curves.tmp_manual$V1[j],
              yend = data_for_curves.tmp_manual$V2[j],
              cluster2 = data_for_curves.tmp_manual$Cluster.Name[j]))
      }
    }
  }
  curves.tmp_manual <-
    rbind(curves.tmp_manual, clone_id_curves.tmp_manual)
}

# set clor pallette
pal.cluster_renumbered = toupper(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffed6f','#b15928', "gray","black","blue","red")) 

# make airline plot of the full TCR sharing
CMV_TCR_airline_plot_manual <- plot_cells(
  cds_no_MAIT_no_9,
  color_cells_by = "Cluster.Name", cell_size=1, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  geom_curve(
    data = curves.tmp_manual,
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.07,
    alpha=0.5) +
  theme(text = element_text(size = 16)) +
  labs(title = "CMV-Specific TRAs")

## FIGURE 6D:
ggsave(CMV_TCR_airline_plot_manual, file = file.path(plotDir,"CMV_TCR_airline_plot_manual_ALPHA.pdf"), device = "pdf",height = 5, width = 5)

#### FIGURE S14A: Airline plot of EBV specific TRAs on the UMAP ####

# Load previous necessary objects
load(file = file.path(resultDir, "all_libs_tcrs_cluster_pub_priv_ag_species_TRA.Rdata"))
load(file.path(resultDir, "tra_cluster_binary_trajectories_4_onward.Rdata"))
# start with object without cluster 9
load(file.path(resultDir,"cds_no_MAIT_no_9.Rdata"))

# subset for CMV AND EBV specific cdr3
all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV <- all_libs_tcrs_cluster_pub_priv_ag_species_TRA %>% 
  filter( antigen.species == "EBV")
#save(all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV , file = file.path(resultDir,"all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV.Rdata"))

# Join with differentiation trajectory for each cell
all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV_traj_4_onward <- left_join(all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV, tra_cluster_binary_trajectories_4_onward)

# join with umap to get number of cells with each
cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_ebv <- cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP %>% 
  filter(cdr3_nt_TRA %in% all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV$cdr3_nt_TRA)
# save for plotting airline plot 
#save(cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_ebv, file = file.path(resultDir,"cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_ebv.Rdata"))

# Use the UMAP coldata where the TCR information has been joined

# make a copy of the annotation plus UMAP coordinates, for easier manipulation
data.tmp_manual <-
  as.data.frame(colData(cds_all_tcr_manual)) %>%
  # Combine the colData with the reduced dimension representation coordinated for these clones
  cbind(
    as.data.frame(reducedDims(cds_all_tcr_manual)$UMAP) %>%
      magrittr::set_colnames(c("V1", "V2"))) %>% 
  dplyr::rename("clone_id_tcr_graph_clonal_expansion" = "tcr_id_TRA") %>%
  filter(barcode_original %in% cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_ebv$barcode_original)

# create data frame to store links
curves.tmp_manual <-
  data.frame(
    clone_id_tcr_graph_clonal_expansion = character(),
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric())

# loop over each clone, and extract UMAP coordinates for cells from the same clone
for (clone_id.tmp_manual in na.omit(unique(data.tmp_manual$clone_id_tcr_graph_clonal_expansion))) {
  clone_id_curves.tmp_manual <- curves.tmp_manual[0,]
  data_for_curves.tmp_manual <-
    data.tmp_manual %>%
    dplyr::filter(clone_id_tcr_graph_clonal_expansion %in% clone_id.tmp_manual)
  if (nrow(data_for_curves.tmp_manual) > 1) {
    for (i in 1:(nrow(data_for_curves.tmp_manual)-1)) {
      for (j in (i+1):nrow(data_for_curves.tmp_manual)) {
        clone_id_curves.tmp_manual <-
          rbind(
            clone_id_curves.tmp_manual,
            list(
              clone_id_tcr_graph_clonal_expansion =
                data_for_curves.tmp_manual$clone_id_tcr_graph_clonal_expansion[i],
              x = data_for_curves.tmp_manual$V1[i],
              y = data_for_curves.tmp_manual$V2[i],
              cluster1 = data_for_curves.tmp_manual$Cluster.Name[i],
              xend = data_for_curves.tmp_manual$V1[j],
              yend = data_for_curves.tmp_manual$V2[j],
              cluster2 = data_for_curves.tmp_manual$Cluster.Name[j]))
      }
    }
  }
  curves.tmp_manual <-
    rbind(curves.tmp_manual, clone_id_curves.tmp_manual)
}

# set clor pallette
pal.cluster_renumbered = toupper(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffed6f','#b15928', "gray","black","blue","red")) 

# make airline plot of the full TCR sharing
EBV_TCR_airline_plot_manual <- plot_cells(
  cds_no_MAIT_no_9,
  color_cells_by = "Cluster.Name", cell_size=1, show_trajectory_graph = FALSE,
  group_label_size=8) +
  scale_color_manual(values=pal.cluster_renumbered) +
  geom_curve(
    data = curves.tmp_manual,
    mapping = aes(x=x, y=y, xend=xend, yend=yend),
    size = 0.07,
    alpha=0.5) +
  theme(text = element_text(size = 16)) +
  labs(title = "EBV-Specific TRA")

# FIGURE S14A
ggsave(EBV_TCR_airline_plot_manual, file = file.path(plotDir, "EBV_TCR_airline_plot_manual_ALPHA.pdf"), device = "pdf",height = 5, width = 5)

#### FIGURE S14B: Plot CMV specificity for each donor separately ####

# load previously saved data
load(file.path(resultDir,"cds_no_MAIT_no_9.Rdata"))
load(file = file.path(resultDir, "all_libs_tcrs_cluster_pub_priv_ag_species_TRA.Rdata"))

# subset for CMV specific cdr3
all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV <- all_libs_tcrs_cluster_pub_priv_ag_species_TRA %>% filter(antigen.species == "CMV")

donor_cmv_list <- unique(all_libs_tcrs_cluster_pub_priv_ag_species_TRA_CMV$Donor.ID)

for(donor in donor_cmv_list) {
  
  cds_all_tcr_manual_donor <- cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% 
                                                   row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% 
                                                               filter( Donor.ID == donor))]
  cds_all_tcr_manual_donor <- cds_all_tcr_manual_donor[,colnames(cds_all_tcr_manual_donor) %in% row.names(colData(cds_all_tcr_manual_donor) %>% as.data.frame() %>% filter(Cluster.Name != 9))]
  
  # make a copy of the annotation plus UMAP coordinates, for easier manipulation
  data.tmp_manual_donor <-
    as.data.frame(colData(cds_all_tcr_manual_donor)) %>%
    # Combine the colData with the reduced dimension representation coordinated for these clones
    cbind(
      as.data.frame(reducedDims(cds_all_tcr_manual_donor)$UMAP) %>%
        magrittr::set_colnames(c("V1", "V2"))) %>% 
    dplyr::rename("clone_id_tcr_graph_clonal_expansion" = "tcr_id_TRA") %>%
    filter(barcode_original %in% cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_cmv$barcode_original)
  
  # create data frame to store links
  # create data frame to store links
  curves.tmp_manual_donor <-
    data.frame(
      clone_id_tcr_graph_clonal_expansion = character(),
      x = numeric(),
      y = numeric(),
      xend = numeric(),
      yend = numeric())
  
  # loop over each clone, and extract UMAP coordinates for cells from the same clone
  for (clone_id.tmp_manual_donor in na.omit(unique(data.tmp_manual_donor$clone_id_tcr_graph_clonal_expansion))) {
    clone_id_curves.tmp_manual_donor <- curves.tmp_manual_donor[0,]
    data_for_curves.tmp_manual_donor <-
      data.tmp_manual_donor %>%
      dplyr::filter(clone_id_tcr_graph_clonal_expansion %in% clone_id.tmp_manual_donor)
    if (nrow(data_for_curves.tmp_manual_donor) > 1) {
      for (i in 1:(nrow(data_for_curves.tmp_manual_donor)-1)) {
        for (j in (i+1):nrow(data_for_curves.tmp_manual_donor)) {
          clone_id_curves.tmp_manual_donor <-
            rbind(
              clone_id_curves.tmp_manual_donor,
              list(
                clone_id_tcr_graph_clonal_expansion =
                  data_for_curves.tmp_manual_donor$clone_id_tcr_graph_clonal_expansion[i],
                x = data_for_curves.tmp_manual_donor$V1[i],
                y = data_for_curves.tmp_manual_donor$V2[i],
                cluster1 = data_for_curves.tmp_manual_donor$Cluster.Name[i],
                xend = data_for_curves.tmp_manual_donor$V1[j],
                yend = data_for_curves.tmp_manual_donor$V2[j],
                cluster2 = data_for_curves.tmp_manual_donor$Cluster.Name[j]))
        }
      }
    }
    curves.tmp_manual_donor <-
      rbind(curves.tmp_manual_donor, clone_id_curves.tmp_manual_donor)
  }
  
  # set clor pallette
  pal.cluster_renumbered = toupper(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffed6f','#b15928', "gray","black","blue","red")) 
  
  donor_plot <- T1DAL_ITN_ID_dictionary %>% filter(Donor.ID == donor)
  
  # make airline plots of the full TCR sharing
  full_TCR_airline_plot_manual_donor <- plot_cells(
    cds_all_tcr_manual_donor,
    color_cells_by = "Cluster.Name", cell_size=1, show_trajectory_graph = FALSE,
    group_label_size=8) +
    scale_color_manual(values=pal.cluster_renumbered) +
    geom_curve(
      data = curves.tmp_manual_donor,
      mapping = aes(x=x, y=y, xend=xend, yend=yend),
      size = 0.5,
      alpha=0.1) + 
    ggtitle(paste("CMV-specific TRA Sharing\nBetween All Clusters in",  donor_plot$masked_public_PID)) +
    theme(text = element_text(size = 16))
  
  plots = list(full_TCR_airline_plot_manual_donor)
  
  # PLOTTED AS SEPARATE FIGURES FOR EACH DONOR FOR FIGURE S14B
  library(gridExtra)
  pdf(paste0(plotDir,"/",donor,"_CMV_specificity_airline.pdf"), onefile = TRUE, width = 5, height = 6)
  do.call("grid.arrange", plots)  
  dev.off()
  
}

#### FIGURE S14A: Plot EBV specificity for each donor separately ####
load(file = file.path(resultDir, "all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV.Rdata"))
load(file = file.path(resultDir, "cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_ebv.Rdata"))

donor_ebv_list <- unique(all_libs_tcrs_cluster_pub_priv_ag_species_TRA_EBV$Donor.ID)

for(donor in donor_ebv_list) {
  
  cds_all_tcr_manual_donor <- cds_all_tcr_manual[,colnames(cds_all_tcr_manual) %in% 
                                                   row.names(colData(cds_all_tcr_manual) %>% as.data.frame() %>% 
                                                               filter( Donor.ID == donor))]
  cds_all_tcr_manual_donor <- cds_all_tcr_manual_donor[,colnames(cds_all_tcr_manual_donor) %in% row.names(colData(cds_all_tcr_manual_donor) %>% as.data.frame() %>% filter(Cluster.Name != 9))]
  
  # make a copy of the annotation plus UMAP coordinates, for easier manipulation
  data.tmp_manual_donor <-
    as.data.frame(colData(cds_all_tcr_manual_donor)) %>%
    # Combine the colData with the reduced dimension representation coordinated for these clones
    cbind(
      as.data.frame(reducedDims(cds_all_tcr_manual_donor)$UMAP) %>%
        magrittr::set_colnames(c("V1", "V2"))) %>% 
    dplyr::rename("clone_id_tcr_graph_clonal_expansion" = "tcr_id_TRA") %>%
    filter(barcode_original %in% cds_tcr_manual_metadata_ordered_TRA_trajectory_UMAP_ebv$barcode_original)
  
  # create data frame to store links
  # create data frame to store links
  curves.tmp_manual_donor <-
    data.frame(
      clone_id_tcr_graph_clonal_expansion = character(),
      x = numeric(),
      y = numeric(),
      xend = numeric(),
      yend = numeric())
  
  # loop over each clone, and extract UMAP coordinates for cells from the same clone
  for (clone_id.tmp_manual_donor in na.omit(unique(data.tmp_manual_donor$clone_id_tcr_graph_clonal_expansion))) {
    clone_id_curves.tmp_manual_donor <- curves.tmp_manual_donor[0,]
    data_for_curves.tmp_manual_donor <-
      data.tmp_manual_donor %>%
      dplyr::filter(clone_id_tcr_graph_clonal_expansion %in% clone_id.tmp_manual_donor)
    if (nrow(data_for_curves.tmp_manual_donor) > 1) {
      for (i in 1:(nrow(data_for_curves.tmp_manual_donor)-1)) {
        for (j in (i+1):nrow(data_for_curves.tmp_manual_donor)) {
          clone_id_curves.tmp_manual_donor <-
            rbind(
              clone_id_curves.tmp_manual_donor,
              list(
                clone_id_tcr_graph_clonal_expansion =
                  data_for_curves.tmp_manual_donor$clone_id_tcr_graph_clonal_expansion[i],
                x = data_for_curves.tmp_manual_donor$V1[i],
                y = data_for_curves.tmp_manual_donor$V2[i],
                cluster1 = data_for_curves.tmp_manual_donor$Cluster.Name[i],
                xend = data_for_curves.tmp_manual_donor$V1[j],
                yend = data_for_curves.tmp_manual_donor$V2[j],
                cluster2 = data_for_curves.tmp_manual_donor$Cluster.Name[j]))
        }
      }
    }
    curves.tmp_manual_donor <-
      rbind(curves.tmp_manual_donor, clone_id_curves.tmp_manual_donor)
  }
  
  # set clor pallette
  pal.cluster_renumbered = toupper(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffed6f','#b15928', "gray","black","blue","red")) 
  
  donor_plot <- T1DAL_ITN_ID_dictionary %>% filter(Donor.ID == donor)
  
  # make airline plots of the full TCR sharing
  full_TCR_airline_plot_manual_donor <- plot_cells(
    cds_all_tcr_manual_donor,
    color_cells_by = "Cluster.Name", cell_size=1, show_trajectory_graph = FALSE,
    group_label_size=8) +
    scale_color_manual(values=pal.cluster_renumbered) +
    geom_curve(
      data = curves.tmp_manual_donor,
      mapping = aes(x=x, y=y, xend=xend, yend=yend),
      size = 0.5,
      alpha=0.1) + 
    ggtitle(paste("EBV-specific TRA Sharing\nBetween All Clusters in",  donor_plot$masked_public_PID)) +
    theme(text = element_text(size = 16))
  
  plots = list(full_TCR_airline_plot_manual_donor)
  
  # PLOTTED AS SEPARATE PLOT FOR EACH DONOR FOR FIGURE S14A AND ONLY HIGHLIGHTED TWO SPECIFIC DONORS
  library(gridExtra)
  pdf(paste0(plotDir,"/",donor,"_EBV_specificity_airline.pdf"), onefile = TRUE, width = 5, height = 6)
  do.call("grid.arrange", plots)  
  dev.off()
  
}

