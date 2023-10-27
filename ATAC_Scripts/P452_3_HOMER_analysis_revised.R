#### HOMER Transcription factor Motif Enrichment Analysis ####

# This script exports differentially accessible peaks as BED files to 
# run HOMER, and also analyzes the results 


#### LOAD LIBRARIES ####

library(ComplexHeatmap)
library(DiffBind)
library(tidyverse)
library(biomaRt)
library(cowplot)
library(egg)
library(biomartr)
library(stringr)
library(universalmotif)
library(Biostrings)
library(igraph)

# installed marge package for parsing homer results
devtools::install_github('robertamezquita/marge', ref = 'master')
library(marge)

#### LOAD DATA and SET PATHS ####

# setwd
setwd("/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/ATAC_Analysis")

# Paths 
baseDir <- "/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES"
plotDir <- file.path(baseDir,'/T1DAL_cleanup/FIGURES/ATAC_REPEAT_FIGURES')
resultDir <- file.path(baseDir,'P452_3_SAVED_DATA')
annotationDir <- file.path(baseDir, "RAW_DATA/P452_3/")
homerDir <- file.path(baseDir, "ATAC_Analysis/HOMER_Analysis")

# Load R data files with the lists of DA peaks from each contrast of interest: Cd57pos vs PD1, CD57 vs DN, PD1 vs DN
load(file.path(resultDir, "ATACseqData_P452_3_norm_db_anno_CD57_vs_PD1.RData"))
load(file.path(resultDir, "ATACseqData_P452_3_norm_db_anno_CD57pos_DN.Rdata"))
load(file.path(resultDir, "ATACseqData_P452_3_norm_db_anno_CD57minus_DN.RData"))

# Load normalized data for export of contrasts
load(file.path(resultDir, "ATACseqData_P452_3_norm_final.RData"))

## Load data on differentially expressed genes between PD1 and CD57
# load all the data
load(file = file.path(baseDir, "EW_T1DAL_Results/cluster_ex_terms_sig.Rdata"))

# load tables split up
# negative = higher in CD57
# positive = higher in PD1
cluster_ex_terms_sig_up_PD1 <- read.csv(file = file.path(baseDir, "EW_T1DAL_Results/cluster_ex_terms_sig_up.csv"))
cluster_ex_terms_sig_down_CD57 <- read.csv( file = file.path(baseDir, "EW_T1DAL_Results/cluster_ex_terms_sig_down.csv"))

#### Export contrast results as BED files ####

### Function to export three contrasts
export_contrast <- function(x, data) {
  x <- dba.report(data, bAll= TRUE,
                  method=DBA_EDGER, 
                  contrast = x,  # contrast 1 = CD57_pos_DP       8 CD57_minus_DP 8
                  bUsePval = FALSE, # tells it to use the FDR threshold 
                  th=0.05,DataType = DBA_DATA_FRAME)
  # add PeakID for HOMER
  x$PeakID <- as.numeric(rownames(x))
  x
}
contrast_list <- c(1,2,3)
contrast_names <- c("CD57_vs_PD1", "CD57pos_vs_DN","CD57minus_vs_DN")

contrasts <- lapply(contrast_list,export_contrast, ATACseqData_P452_3_norm)
names(contrasts) <- contrast_names
list2env(contrasts, .GlobalEnv)

### Export each full contrast as BED 
for (contrast in names(contrasts)) {
  x <- contrasts[[contrast]]
  BED <- x %>% select(Chr, Start, End, PeakID) %>% mutate(Chr = paste0("chr",Chr)) %>% 
    # HOMER wants a unique PEAK identifier for each and chromosome in UCSC format
    distinct()
  write.table(BED, file = file.path(homerDir, paste(contrast, "BED.txt", sep = "_") ), sep="\t", quote=F, row.names=F, col.names=F)
}

# Note HOMER also asks for + or - strand information which I'm not including right now
export_contrast <- function(x, data) {
  x <- dba.report(data, bAll= TRUE,
                  method=DBA_EDGER, 
                  contrast = x,  # contrast 1 = CD57_pos_DP       8 CD57_minus_DP 8
                  bUsePval = FALSE, # tells it to use the FDR threshold 
                  th=0.05,DataType = DBA_DATA_FRAME)
  # add PeakID for HOMER
  x$PeakID <- as.numeric(rownames(x))
  x
}
contrast_list <- c(1,2,3)
contrast_names <- c("CD57_vs_PD1", "CD57pos_vs_DN","CD57minus_vs_DN")

contrasts <- lapply(contrast_list,export_contrast, ATACseqData_P452_3_norm)
names(contrasts) <- contrast_names
list2env(contrasts, .GlobalEnv)

### Export each full contrast as BED 
for (contrast in names(contrasts)) {
  x <- contrasts[[contrast]]
  BED <- x %>% select(Chr, Start, End, PeakID) %>% mutate(Chr = paste0("chr",Chr)) %>% 
    # HOMER wants a unique PEAK identifier for each and chromosome in UCSC format
    distinct()
  write.table(BED, file = file.path(homerDir, paste(contrast, "BED.txt", sep = "_") ), sep="\t", quote=F, row.names=F, col.names=F)
}

### Export contrasts separated into positive and negative LFC

for (contrast in names(contrasts)) {
  x <- contrasts[[contrast]]
  BED_UP <- x %>% filter(Fold >=0) %>% 
    select(Chr, Start, End, PeakID) %>% mutate(Chr = paste0("chr",Chr)) %>% 
    # HOMER wants a unique PEAK identifier for each and chromosome in UCSC format
    distinct()
  BED_DOWN <- x %>% filter(Fold < 0) %>%
    select(Chr, Start, End, PeakID) %>% mutate(Chr = paste0("chr",Chr)) %>% 
    # HOMER wants a unique PEAK identifier for each and chromosome in UCSC format
    distinct()
  write.table(BED_UP, file = file.path(homerDir, paste(contrast, "BED_UP.txt", sep = "_") ), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(BED_DOWN, file = file.path(homerDir, paste(contrast, "BED_DOWN.txt", sep = "_") ), sep="\t", quote=F, row.names=F, col.names=F)
  
}


#### RUN HOMER VIA COMMAND LINE #### 

# I locally installed HOMER. I wrote a script to run the procedure, called run_HOMER_P452_3_edited.sh

# The results for the analysis are here:
# /Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/ATAC_Analysis/HOMER_Analysis


#### IMPORT HOMER RESULTS from CD57 vs PD1 contrast ####

# Import HOMER results using marge API package
denovo_PD1_tbl <- marge::read_denovo_results(path = file.path(homerDir, "CD57_vs_PD1_DOWN"), homer_dir = TRUE)
denovo_motif_PD1 <- marge::read_denovo_html_results(path = file.path(homerDir, "CD57_vs_PD1_DOWN"), homer_dir = TRUE)

denovo_CD57_tbl <- marge::read_denovo_results(path = file.path(homerDir, "CD57_vs_PD1_UP"), homer_dir = TRUE)
denovo_motif_CD57 <- marge::read_denovo_html_results(path = file.path(homerDir, "CD57_vs_PD1_UP"), homer_dir = TRUE)

# for each motif you have the de novo motif info and then all matches to known motifs in a table 

### concatenate_motif_results fpr PD1
denovo_PD1_motifs <- as.data.frame(do.call(rbind, lapply(denovo_motif_PD1, as.vector)))
denovo_PD1_motifs <- cbind(motif=rownames(denovo_PD1_motifs), denovo_PD1_motifs)
PD1_motif_information <- as.data.frame(do.call(rbind, lapply(denovo_PD1_motifs$Motif_information, as.vector)))
PD1_motif_information$p_value <- as.numeric(PD1_motif_information$p_value)
PD1_motif_information$motif_number <- rownames(PD1_motif_information)

# missing motif family column for element 6, add
denovo_PD1_motifs$Matches_to_known_motifs$`6`$motif_family <- NA
# make dataframe with the concatenated list of motif matches
PD1_matches_to_known_motifs <- as.data.frame(do.call(rbind, lapply(denovo_PD1_motifs$Matches_to_known_motifs, as.vector)))
PD1_matches_to_known_motifs <- cbind(motif_match=rownames(PD1_matches_to_known_motifs), PD1_matches_to_known_motifs)

### Concatenate motif results for CD57
denovo_CD57_motifs <- as.data.frame(do.call(rbind, lapply(denovo_motif_CD57, as.vector)))
denovo_CD57_motifs <- cbind(motif=rownames(denovo_CD57_motifs), denovo_CD57_motifs)
CD57_motif_information <- as.data.frame(do.call(rbind, lapply(denovo_CD57_motifs$Motif_information, as.vector)))
CD57_motif_information$p_value <- as.numeric(CD57_motif_information$p_value)
CD57_motif_information$motif_number <- rownames(CD57_motif_information)

#  make dataframe with the concatenated list of motif matches
CD57_matches_to_known_motifs <- as.data.frame(do.call(rbind, lapply(denovo_CD57_motifs$Matches_to_known_motifs, as.vector)))
CD57_matches_to_known_motifs <- cbind(motif_match=rownames(CD57_matches_to_known_motifs), CD57_matches_to_known_motifs)

## filter results to keep only de novo motifs with p value < 1e-12 - which are likely to be false positives based on HOMER documentation
class(PD1_motif_information$p_value) # character
PD1_motif_information_keep <- PD1_motif_information %>% filter(p_value <=1e-12)
CD57_motif_information_keep <- CD57_motif_information %>% filter(p_value <=1e-12) 

PD1_matches_to_known_motifs_keep <- PD1_matches_to_known_motifs %>% filter(original_alignment %in% PD1_motif_information_keep$consensus)
CD57_matches_to_known_motifs_keep <- CD57_matches_to_known_motifs %>% filter(original_alignment %in% CD57_motif_information_keep$consensus)

# export motif information and matches to known motifs
write.csv(PD1_matches_to_known_motifs_keep, file = file.path(resultDir, "PD1_matches_to_known_motifs_keep.csv"), row.names = FALSE)
write.csv(CD57_matches_to_known_motifs_keep, file = file.path(resultDir, "CD57_matches_to_known_motifs_keep.csv"), row.names = FALSE)

write.csv(PD1_motif_information_keep, file = file.path(resultDir, "PD1_motif_information_keep.csv"), row.names = FALSE)
write.csv(CD57_motif_information_keep, file = file.path(resultDir, "CD57_motif_information_keep.csv"), row.names = FALSE)

#### Find TFs in list of CD57 vs PD1 DEGs using biomart ####

# use the getGO biomartr function to get GOs and identify transcription factos

# get GO terms for the PD1 list
cluster_ex_terms_sig_up_PD1_GO_tbl <- biomartr::getGO(organism = "Homo sapiens", 
                                                      genes    = cluster_ex_terms_sig_up_PD1$gene_short_name,
                                                      filters  = "hgnc_symbol")

# get GO terms for the CD57 list
cluster_ex_terms_sig_down_CD57_GO_tbl <- biomartr::getGO(organism = "Homo sapiens", 
                                                         genes    = cluster_ex_terms_sig_down_CD57$gene_short_name,
                                                         filters  = "hgnc_symbol")

### Identify GO terms related to transcription ###

## Investigate PD1
cluster_ex_terms_sig_up_PD1_GO_tbl_transcription <- cluster_ex_terms_sig_up_PD1_GO_tbl %>% 
  filter(grepl("transcription", goslim_goa_description)) %>% distinct()
unique(cluster_ex_terms_sig_up_PD1_GO_tbl_transcription[,c("goslim_goa_description","goslim_goa_accession")])

# 1                       DNA-templated transcription           GO:0006351 : The synthesis of an RNA transcript from a DNA template.
# 2         regulation of DNA-templated transcription           GO:0006355 : Any process that modulates the frequency, rate or extent of cellular DNA-templated transcription.
# 10                 transcription regulator activity           GO:0140110 : A molecular function that controls the rate, timing and/or magnitude of gene transcription.
# 57 general transcription initiation factor activity           GO:0140223 : General transcription factors (GTFs) bind to and open promoter DNA, initiate RNA synthesis and stimulate the escape of the polymerase from the promoter

PD1_GO_gene_list <- unique(cluster_ex_terms_sig_up_PD1_GO_tbl_transcription$hgnc_symbol)

## Investigate CD57
cluster_ex_terms_sig_down_CD57_GO_tbl_transcription <- cluster_ex_terms_sig_down_CD57_GO_tbl %>% 
  filter(grepl("transcription", goslim_goa_description)) %>% distinct()
unique(cluster_ex_terms_sig_down_CD57_GO_tbl_transcription[,c("goslim_goa_description","goslim_goa_accession")])

#goslim_goa_description goslim_goa_accession
#1                        DNA-templated transcription           GO:0006351
#2          regulation of DNA-templated transcription           GO:0006355
#3                   transcription regulator activity           GO:0140110
#125 general transcription initiation factor activity           GO:0140223


#### FIND OVERLAP BETWEEN ATAC TF ENRICHMENT AND DEGS AND PLOT VOLCANO ####

PD1_matches_to_known_motifs_keep 
CD57_matches_to_known_motifs_keep

# find DEGs that are in the list of transcription factor motifs
PD1_matches_to_known_motifs_keep_DEG <- PD1_matches_to_known_motifs_keep[str_to_upper(PD1_matches_to_known_motifs_keep$motif_name) %in% cluster_ex_terms_sig_up_PD1_GO_tbl_transcription$hgnc_symbol,] 
CD57_matches_to_known_motifs_keep_DEG <- CD57_matches_to_known_motifs_keep[str_to_upper(CD57_matches_to_known_motifs_keep$motif_name) %in% cluster_ex_terms_sig_down_CD57_GO_tbl_transcription$hgnc_symbol,] 

# Join HOMER p-value for motif as a whole and DEG log expression
# remember: normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307

PD1_matches_to_known_motifs_keep_DEG_heatmap <- PD1_matches_to_known_motifs_keep_DEG %>% 
  separate(motif_match, into = c("motif_number","motif_match_number")) %>%
  # join motif information for that entire motif to get enrichment p value
  left_join(.,PD1_motif_information, by = "motif_number") %>%
  dplyr::rename(gene_short_name = "motif_name.x") %>%
  mutate(gene_short_name = str_to_upper(gene_short_name)) %>%
  # join DEG info 
  left_join(., cluster_ex_terms_sig_up_PD1) %>%
  mutate(group = "PD1")

CD57_matches_to_known_motifs_keep_DEG_heatmap <- CD57_matches_to_known_motifs_keep_DEG %>% 
  separate(motif_match, into = c("motif_number","motif_match_number")) %>%
  # join motif information for that entire motif to get enrichment p value
  left_join(.,CD57_motif_information, by = "motif_number") %>%
  dplyr::rename(gene_short_name = "motif_name.x") %>%
  mutate(gene_short_name = str_to_upper(gene_short_name)) %>%
  # join DEG info 
  left_join(., cluster_ex_terms_sig_down_CD57) %>%
  mutate(group = "CD57")

## Join together tables for joint volcano plotting
PD1_CD57_volcano <- rbind(PD1_matches_to_known_motifs_keep_DEG_heatmap,CD57_matches_to_known_motifs_keep_DEG_heatmap)
# remove duplicate rows for hits to multiple databases
PD1_CD57_volcano <- PD1_CD57_volcano %>% distinct(gene_short_name, p_value, normalized_effect, group, q_value)


## Plot volcano
PD1_CD57_volcano_plot <- ggplot(PD1_CD57_volcano, aes(y = -log10(p_value), x = normalized_effect, 
                                                      color = group, label = gene_short_name)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 5, max.overlaps = 10 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#93a24eff", "#ba4d4cc7")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 12), text = element_text(size = 12),
        legend.position="bottom") +
  labs(x = "log2 Fold Change of TF DEGs",y = "Differential Peak Motif\nEnrichment p-value", title = "Differential Expression of\nEnriched Transcription Factors") +
  geom_vline(xintercept = 0, linetype = "dotted") 

# export plot
ggsave(plot = PD1_CD57_volcano_plot, 
       file = file.path(plotDir, "PD1_CD57_volcano_plot.pdf"),device = "pdf", width = 4, height = 4)

## Replot the volcano plot plotting the DEG q value instead of p value

PD1_CD57_volcano_plot_log10q <- ggplot(PD1_CD57_volcano, aes(y = -log10(p_value), x = -log10(q_value), 
                                                             color = group, label = gene_short_name)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 5, max.overlaps = 10 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#ba4d4cc7", "#fdbf6fff")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 12), text = element_text(size = 12),
        legend.position="bottom") +
  labs(x = "-log10 q-value TF DEGs",y = "Differential Peak Motif\nEnrichment -log10 p-value", title = "Differential Expression of\nEnriched Transcription Factors") 
#geom_vline(xintercept = 0, linetype = "dotted") 

# export plot
ggsave(plot = PD1_CD57_volcano_plot_log10q, 
       file = file.path(plotDir, "PD1_CD57_volcano_plot_log10q.pdf"),device = "pdf", width = 4, height = 4)


#### Plot R vs NR log2 fold change with TF motif enrichment volcano plot ####

# Load R vs NR significant transcripts within PD1 cluster and CD57 cluster
# positive values mean increased expression in responders
cluster_response_terms_CD57_sig <- read.csv(file =  "../EW_T1DAL_Results/cluster_response_terms_CD57_sig.csv")
cluster_response_terms_PD1_sig <- read.csv(file =  "../EW_T1DAL_Results/cluster_response_terms_PD1_sig.csv")

# make list of all enriched motifs
PD1_matches_to_known_motifs_keep
CD57_matches_to_known_motifs_keep

# find R vs NR DEGs that are in the list of all enriched transcription factor motifs
# only join with the PD1 or CD57 motif enrichment that is applicable to their gene list of interest
cluster_response_terms_CD57_sig_motifs_CD57 <- CD57_matches_to_known_motifs_keep[str_to_upper(CD57_matches_to_known_motifs_keep$motif_name) %in% cluster_response_terms_CD57_sig$gene_short_name,] 
# RUNX3 and TGIF1 overlapped here
cluster_response_terms_PD1_sig_motifs_PD1 <- PD1_matches_to_known_motifs_keep[str_to_upper(PD1_matches_to_known_motifs_keep$motif_name) %in% cluster_response_terms_PD1_sig$gene_short_name,] 
# no overlaps 

# Join HOMER p-value for motif as a whole and DEG log expression
# remember: normalized effect is the log2 fold change normalized by the intercept https://github.com/cole-trapnell-lab/monocle3/issues/307

cluster_response_terms_CD57_sig_motifs_CD57_heatmap <- cluster_response_terms_CD57_sig_motifs_CD57 %>% 
  separate(motif_match, into = c("motif_number","motif_match_number")) %>%
  # join motif information for that entire motif to get enrichment p value
  left_join(.,CD57_motif_information, by = "motif_number") %>%
  dplyr::rename(gene_short_name = "motif_name.x") %>%
  mutate(gene_short_name = str_to_upper(gene_short_name)) %>%
  # join DEG info 
  left_join(., cluster_response_terms_CD57_sig) %>%
  mutate(group = "Responders")

## Plot volcano
CD57_R_vs_NR_TF_volcano_plot <- ggplot(cluster_response_terms_CD57_sig_motifs_CD57_heatmap, 
                                       aes(y = -log10(p_value), x = normalized_effect, color = group,
                                           label = gene_short_name)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 5, max.overlaps = 10 ) +
  scale_color_manual(values = c("#5ba965")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 12), text = element_text(size = 12),
        legend.position="bottom") +
  labs(x = "log2 Fold Change of CD57+ TF DEGs",y = "Differential Peak Motif\nEnrichment p-value", title = "CD57 Cluster R vs NR\nDifferential Expression of Enriched TFs") +
  geom_vline(xintercept = 0, linetype = "dotted") 

# export plot
ggsave(plot = CD57_R_vs_NR_TF_volcano_plot, 
       file = file.path(plotDir, "CD57_R_vs_NR_TF_volcano_plot.pdf"),device = "pdf", width = 4, height = 4)

### Replot PD1 vs CD57 volcano adding label for RUNX3 and TGIF1

## Join together tables for joint volcano plotting
PD1_CD57_volcano_labels <- PD1_CD57_volcano 
PD1_CD57_volcano_labels[PD1_CD57_volcano_labels$gene_short_name == "RUNX3",31] <- "R vs. NR\nOverlap"
PD1_CD57_volcano_labels[PD1_CD57_volcano_labels$gene_short_name == "TGIF1",31] <- "R vs. NR\nOverlap"

## Plot volcano
PD1_CD57_volcano_labels_plot <- ggplot(PD1_CD57_volcano_labels, aes(y = -log10(p_value), x = normalized_effect, 
                                                                    color = group, label = gene_short_name)) +
  theme_minimal() +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(color = "black", size = 5, max.overlaps = 10 ) +
  #annotate("text", x = -1, y = 6.5, label = "PPI enrichment\np-value 2.74e-07", size = 3) +
  scale_color_manual(values = c("#ba4d4cc7", "#fdbf6fff", "gray")) +
  theme(legend.title= element_blank(), legend.text = element_text(size = 12), text = element_text(size = 12),
        legend.position="bottom") +
  labs(x = "log2 Fold Change of TF DEGs",y = "Differential Peak Motif\nEnrichment p-value", title = "Differential Expression of\nEnriched Transcription Factors") +
  geom_vline(xintercept = 0, linetype = "dotted") 

# export plot
ggsave(plot = PD1_CD57_volcano_labels_plot, 
       file = file.path(plotDir, "PD1_CD57_volcano_labels_plot.pdf"),device = "pdf", width = 4, height = 4)


#### FIND DEGs within chromatin regions defined by GREAT ####

# load GREAT all genes files
CD57_GREAT_genes <- read.table(file = "./GREAT_Results/CD57pos_vs_DN_BED_UP_gene_to_genomic_region-hg38-all-gene.txt",fill = TRUE,
                               sep = "\t", col.names = c("gene_short_name","region_og")) 
PD1_GREAT_genes <- read.table(file = "./GREAT_Results/CD57pos_vs_DN_BED_DOWN_gene_to_genomic_region-hg38-all-gene.txt",fill = TRUE,
                              sep = "\t", col.names = c("gene_short_name","region_og"))

# Load ATAC bed files exported above
CD57_vs_PD1_UP_BED <- read.table(file = file.path(homerDir, "CD57_vs_PD1_BED_UP.txt"),
                                 sep = "\t", col.names = c("chr","start","end","region"))
CD57_vs_PD1_DOWN_BED <- read.table(file = file.path(homerDir, "CD57_vs_PD1_BED_DOWN.txt"),
                                   sep = "\t", col.names = c("chr","start","end","region"))

# Identify overlap between DEG list and GREAT list of genes associated with open chromatin regions
CD57_GREAT_genes_DEG <- CD57_GREAT_genes %>% filter(gene_short_name %in% cluster_ex_terms_sig_down_CD57$gene_short_name)
PD1_GREAT_genes_DEG <- PD1_GREAT_genes %>% filter(gene_short_name %in% cluster_ex_terms_sig_up_PD1$gene_short_name)

# Separate region column and join with ATAC
process_HOMER_region <- function(x, bed) {
  sep <- HOMER_object_list[[x]]
  sep <- sep %>% tidyr::separate_rows(region_og, sep = ",") 
  sep$region_og <- stringr::str_trim(sep$region_og, side= "left")
  sep <- sep %>%
    separate(region_og, into = c("region","distance"), sep = " ", remove = FALSE)
  sep$distance <- str_remove(sep$distance, "\\(")
  sep$distance <- str_remove(sep$distance, "\\)")
  sep <- sep %>% mutate(direction = case_when(
    grepl("\\+", distance) ~"downstream",
    grepl("\\-", distance)  ~ "upstream"
  ))
  sep$region <- as.numeric(sep$region)
  sep$distance_numeric <- str_remove(sep$distance,"\\+")
  sep$distance_numeric <- str_remove(sep$distance_numeric,"\\-")
  sep$distance_numeric <- as.numeric(sep$distance_numeric)
  
  # output
  sep
}

HOMER_object_list <- list(A = CD57_GREAT_genes_DEG, B = PD1_GREAT_genes_DEG)

HOMER_DEG <- lapply(names(HOMER_object_list), process_HOMER_region)
names(HOMER_DEG) <- c("CD57_GREAT_genes_DEG_sep", "PD1_GREAT_genes_DEG_sep")
list2env(HOMER_DEG, envir = .GlobalEnv)

# join ATAC regions 
CD57_GREAT_genes_DEG_sep_ATAC <- left_join(CD57_GREAT_genes_DEG_sep, CD57_vs_PD1_UP_BED) 
PD1_GREAT_genes_DEG_sep_ATAC <- left_join(PD1_GREAT_genes_DEG_sep, CD57_vs_PD1_DOWN_BED) 

# subtract window based on up or down region for that DEG
window_adjust <- function(x) {
  # + means downstream of the TSS
  # - mean upstream of the TSS https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655462/Output 
  # this means that I really only need to extend the start for those that go upstream, and extend the end for those that go further downstream
  
  upstream <- window_adjust_list[[x]] %>% filter(direction == "upstream") %>%
    mutate(start = start - distance_numeric)
  downstream <- window_adjust_list[[x]] %>% filter(direction == "downstream") %>%
    mutate(new_end = start+distance_numeric, # only some rows now have a larger end! for those rows keep new end, for other rows keep original end
           end = case_when(new_end > end ~ new_end,
                           new_end <= end ~ end)) %>% dplyr::select(-new_end)
  # recombine
  window_full <- rbind(upstream,downstream)
  
  # output
  window_full
  
}

window_adjust_list <- list(A = CD57_GREAT_genes_DEG_sep_ATAC, B = PD1_GREAT_genes_DEG_sep_ATAC)
window_ATAC_adjust <- lapply(names(window_adjust_list), window_adjust)
names(window_ATAC_adjust) <- c("CD57_GREAT_genes_DEG_sep_ATAC_window_full","PD1_GREAT_genes_DEG_sep_ATAC_window_full")
list2env(window_ATAC_adjust, .GlobalEnv)

# export new extended regions as BED file
write.table(CD57_GREAT_genes_DEG_sep_ATAC_window_full[,c("chr","start","end","region")],file = "./TFMotifView/CD57_GREAT_genes_DEG_sep_ATAC_window_full_BED.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.table(PD1_GREAT_genes_DEG_sep_ATAC_window_full[,c("chr","start","end","region")],file = "./TFMotifView/PD1_GREAT_genes_DEG_sep_ATAC_window_full_BED.txt", sep="\t", quote=F, row.names=F, col.names=F)

# export full file as CSV for viewing in excel
write.csv(CD57_GREAT_genes_DEG_sep_ATAC_window_full,file = file.path(homerDir, "CD57_GREAT_genes_DEG_sep_ATAC_window_full.csv"))
write.csv(PD1_GREAT_genes_DEG_sep_ATAC_window_full,file = file.path(homerDir, "PD1_GREAT_genes_DEG_sep_ATAC_window_full.csv"))

#### SCAN GREAT REGIONS FOR KEY TF MOTIFS ####

## First extract sequence ranges identified by GREAT with DEGS in them using biostrings
CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges <- CD57_GREAT_genes_DEG_sep_ATAC_window_full[,c("chr","start","end","gene_short_name")]
PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges <- PD1_GREAT_genes_DEG_sep_ATAC_window_full[,c("chr","start","end","gene_short_name","region_og")]

# Write for loop to get peptide sequence for each window and concatenate
CD57_GREAT_genes_DEG_sep_ATAC_seq <- data.frame()
for (row in 1:nrow(CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges)) {
  chr <- CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "chr"]
  chr_rm <-  str_remove(chr,"chr")
  str <-CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "start"]
  stop <- CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "end"]
  
  seq = getSequence(chr = chr_rm, start = str , end= stop, 
                    type = "hgnc_symbol", 
                    seqType = "cdna", 
                    mart = mart)
  
  # join on metadata about the original GREAT window if there were results
  if (nrow(seq) > 0){
    seq_meta <- data.frame(chr = CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "chr"],
                           start = CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "start"],
                           end = CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "end"],
                           gene_short_name = CD57_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "gene_short_name"])
    seq_meta <- do.call("rbind", replicate( 
      nrow(seq), seq_meta, simplify = FALSE))
    print(seq_meta)
    seq_total <- cbind(seq,seq_meta)
    seq_total
    print(seq_total)
    
    CD57_GREAT_genes_DEG_sep_ATAC_seq <- rbind(CD57_GREAT_genes_DEG_sep_ATAC_seq, seq_total)
  }
  
}

### Repeat for the PD1 ranges
# Write for loop to get peptide sequence for each window and concatenate
PD1_GREAT_genes_DEG_sep_ATAC_seq <- data.frame()
for (row in 1:nrow(PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges)) {
  chr <- PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "chr"]
  chr_rm <-  str_remove(chr,"chr")
  str <-PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "start"]
  stop <- PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "end"]
  
  seq = getSequence(chr = chr_rm, start = str , end= stop, 
                    type = "hgnc_symbol", 
                    seqType = "cdna", 
                    mart = mart)
  
  # join on metadata about the original GREAT window if there were results
  if (nrow(seq) > 0){
    seq_meta <- data.frame(chr = PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "chr"],
                           start = PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "start"],
                           end = PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "end"],
                           gene_short_name = PD1_GREAT_genes_DEG_sep_ATAC_window_full_ranges[row, "gene_short_name"])
    seq_meta <- do.call("rbind", replicate( 
      nrow(seq), seq_meta, simplify = FALSE))
    print(seq_meta)
    seq_total <- cbind(seq,seq_meta)
    seq_total
    print(seq_total)
    
    PD1_GREAT_genes_DEG_sep_ATAC_seq <- rbind(PD1_GREAT_genes_DEG_sep_ATAC_seq, seq_total)
  }
  
}

### Filter list to keep those where the query DEG gene and the main hgnc region gene match
# this means I'm focusing on genes with more proximal regulatory regions identified by GREAT regions!
CD57_GREAT_genes_DEG_sep_ATAC_seq_filter <- CD57_GREAT_genes_DEG_sep_ATAC_seq  %>% filter(hgnc_symbol == gene_short_name)
PD1_GREAT_genes_DEG_sep_ATAC_seq_filter <- PD1_GREAT_genes_DEG_sep_ATAC_seq  %>% filter(hgnc_symbol == gene_short_name)

#### Scan each set of sequences for TF motifs of interest ####

# get list of motif sequences identified
PD1_matches_to_known_motifs_keep
CD57_matches_to_known_motifs_keep

### CD57 motifs
# create motif sequences for TBX21
CD57_matches_to_known_motifs_keep_TBX21 <- CD57_matches_to_known_motifs_keep %>% filter(motif_name == "TBX21")
# use the original alignment sequence
CD57_TBX21_motif <- universalmotif::create_motif(CD57_matches_to_known_motifs_keep_TBX21$original_alignment, 
                                                 alphabet ="DNA", name = "TBX21")

# scan each DEG sequence ranges for TBX21 motif
CD57_TBX21_hits <- data.frame()
for (row in 1:nrow(CD57_GREAT_genes_DEG_sep_ATAC_seq_filter)) {
  
  seq <- CD57_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "cdna"]
  DNA_seq <- DNAStringSet(seq)
  scan <- scan_sequences(CD57_TBX21_motif, DNA_seq, threshold.type = "pvalue", threshold = 0.001)
  print(scan)
  
  if (nrow(scan) > 0){
    # create metadata to label the sequence
    seq_meta <- data.frame(chr = CD57_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "chr"],
                           start = CD57_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "start"],
                           end = CD57_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "end"],
                           gene_short_name = CD57_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "gene_short_name"])
    seq_meta <- do.call("rbind", replicate( 
      nrow(scan), seq_meta, simplify = FALSE))
    print(seq_meta)
    scan_total <- cbind(as.data.frame(scan),seq_meta)
    # keep only one hit per sequence, sometimes there are multiple domains present inside
    scan_unique <- scan_total %>% distinct(motif, chr, start, end,gene_short_name, pvalue)
    print(scan_unique)
    # add to results
    CD57_TBX21_hits <- rbind(CD57_TBX21_hits, scan_unique)
    
  }  
}

CD57_TBX21_hits_unique <- CD57_TBX21_hits %>% distinct(motif, gene_short_name, pvalue)

### Repeat for PD1 TCF7
# create motif sequences for TCF7
PD1_matches_to_known_motifs_keep_TCF7 <- PD1_matches_to_known_motifs_keep %>% filter(motif_name == "TCF7")
# use the original alignment sequence
PD1_TCF7_motif <- universalmotif::create_motif(PD1_matches_to_known_motifs_keep_TCF7$original_alignment, 
                                               alphabet ="DNA", name = "TCF7")

# scan each DEG sequence ranges for TCF7 motif
PD1_TCF7_hits <- data.frame()
for (row in 1:nrow(PD1_GREAT_genes_DEG_sep_ATAC_seq_filter)) {
  
  seq <- PD1_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "cdna"]
  DNA_seq <- DNAStringSet(seq)
  scan <- scan_sequences(PD1_TCF7_motif, DNA_seq, threshold.type = "pvalue", threshold = 0.001)
  print(scan)
  
  if (nrow(scan) > 0){
    # create metadata to label the sequence
    seq_meta <- data.frame(chr = PD1_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "chr"],
                           start = PD1_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "start"],
                           end = PD1_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "end"],
                           gene_short_name = PD1_GREAT_genes_DEG_sep_ATAC_seq_filter[row, "gene_short_name"])
    seq_meta <- do.call("rbind", replicate( 
      nrow(scan), seq_meta, simplify = FALSE))
    print(seq_meta)
    scan_total <- cbind(as.data.frame(scan),seq_meta)
    # keep only one hit per sequence, sometimes there are multiple domains present inside
    
    scan_unique <- scan_total %>% distinct(motif, chr, start, end,gene_short_name, pvalue)
    print(scan_unique)
    # add to results
    PD1_TCF7_hits <- rbind(PD1_TCF7_hits, scan_unique)
    
  }  
}

PD1_TCF7_hits_unique <- PD1_TCF7_hits %>% distinct(motif, gene_short_name, pvalue)

#### PLOT VOLCANO PLOT OF INDIVIDUAL GENE RANGE MOTIF ENRICHMENT ####

# create matrix for each list of TFs
CD57_TBX21_hits_unique_matrix <- CD57_TBX21_hits_unique %>% 
  # keep lowest p-value hit per gene name
  group_by(gene_short_name) %>%
  filter( pvalue== min(pvalue)) %>%
  dplyr::select(gene_short_name, pvalue) %>%
  column_to_rownames("gene_short_name")
CD57_TBX21_hits_unique_matrix$pvalue <- as.numeric(CD57_TBX21_hits_unique_matrix$pvalue)

PD1_TCF7_hits_unique_matrix <- PD1_TCF7_hits_unique %>% 
  # keep lowest p-value hit per gene name
  group_by(gene_short_name) %>%
  filter( pvalue== min(pvalue)) %>%
  dplyr::select(gene_short_name, pvalue) %>%
  column_to_rownames("gene_short_name")
PD1_TCF7_hits_unique_matrix$pvalue <- as.numeric(PD1_TCF7_hits_unique_matrix$pvalue)

# plot heatmaps
ComplexHeatmap::Heatmap(CD57_TBX21_hits_unique_matrix)
ComplexHeatmap::Heatmap(PD1_TCF7_hits_unique_matrix)

#### PLOT AS NETWORK PLOT ####

# create edge df
CD57_TBX21_hits_network_edge <- CD57_TBX21_hits_unique_matrix %>% rownames_to_column("node2") %>%
  mutate(node1 = "TBX21", score = -log10(pvalue)) %>% 
  dplyr::select(node1, node2, pvalue, score)
PD1_TCF7_hits_network_edge <- PD1_TCF7_hits_unique_matrix %>% rownames_to_column("node2") %>%
  mutate(node1 = "TCF7", score = -log10(pvalue)) %>% 
  dplyr::select(node1, node2, pvalue, score)

# create node df 
CD57_TBX21_hits_network_node <- data.frame(gene_short_name = unique(c(CD57_TBX21_hits_network_edge$node1, CD57_TBX21_hits_network_edge$node2))) %>% 
  left_join( cluster_ex_terms_sig_down_CD57) %>%
  dplyr::rename(node = gene_short_name) %>%
  mutate(type = case_when(node == "TBX21" ~"TBX21",
                          node!="TBX21"~"DEG")) %>%
  mutate(normalized_effect = normalized_effect *-1)

PD1_TCF7_hits_network_node <- data.frame(gene_short_name = unique(c(PD1_TCF7_hits_network_edge $node1, PD1_TCF7_hits_network_edge$node2))) %>% 
  left_join( cluster_ex_terms_sig_up_PD1) %>%
  dplyr::rename(node = gene_short_name) %>%
  mutate(type = case_when(node == "TCF7" ~"TCF7",
                          node!="TCF7"~"DEG")) 
# create network
CD57_TBX21_hits_network <- graph_from_data_frame(d=CD57_TBX21_hits_network_edge , 
                                                 vertices=CD57_TBX21_hits_network_node, directed = F) 

PD1_TCF7_hits_network <- graph_from_data_frame(d=PD1_TCF7_hits_network_edge, 
                                               vertices=PD1_TCF7_hits_network_node, directed = F) 

# plot networks
# view igraph plotting parameters
#?igraph.plotting
# Set node size
V(CD57_TBX21_hits_network)$size <- V(CD57_TBX21_hits_network)$normalized_effect*5
V(PD1_TCF7_hits_network)$size <- V(PD1_TCF7_hits_network)$normalized_effect*5
V(CD57_TBX21_hits_network)$color <- ifelse(V(CD57_TBX21_hits_network)$type == "TBX21", "#ba543d", "#93a24eff")
V(PD1_TCF7_hits_network)$color <- ifelse(V(PD1_TCF7_hits_network)$type == "TCF7", "#6780d8", "#ba4d4cc7")

# Set edge color 
#E(upset_pairwise_net)$color = E(upset_pairwise_net)$edge_color

l <- layout.reingold.tilford(CD57_TBX21_hits_network)
# plot
pdf(file = file.path(plotDir, "CD57_TBX21_hits_network.pdf") )
plot(CD57_TBX21_hits_network, layout=-l[, 2:1], main = "CD57+ DEGs with TBX21 Motifs",  
     label.font  = 2, vertex.label.color="black", 
     edge.arrow.size = 0.5, 
     vertex.label.cex = 1, 
     vertex.label.dist=2,
     xlim = c(-1,1),ylim=c(-1,1))
dev.off()


l_PD1 <- layout_as_tree(PD1_TCF7_hits_network)
pdf(file = file.path(plotDir, "PD1_TCF7_hits_network.pdf"), width = 6, height = 6 )
plot(PD1_TCF7_hits_network, layout = l_PD1, main = "PD-1+ DEGs with TCF7 Motifs",  
     label.font  = 2, vertex.label.color="black", 
     edge.arrow.size = 0.5, 
     vertex.label.cex = 1, 
     vertex.label.dist=2)
dev.off()

#### PLOT INSTEAD AS BUBBLE PLOT ####

CD57_TBX21_hits_network_bubble <- CD57_TBX21_hits_network_edge %>% dplyr::rename(node = node2) %>% 
  left_join(CD57_TBX21_hits_network_node) %>% mutate(group = "CD57+ Tex", `-log10(p-value)` = -log10(pvalue)) %>%
  dplyr::rename(LFC = normalized_effect)

PD1_TCF7_hits_network_bubble <- PD1_TCF7_hits_network_edge %>% dplyr::rename(node = node2) %>% 
  left_join(PD1_TCF7_hits_network_node) %>%  mutate(group = "PD-1+ Tex", `-log10(p-value)` = -log10(pvalue)) %>%
  dplyr::rename(LFC = normalized_effect)

# plot as bubble plot
CD57_TBX21_hits_network_bubble_plot <- ggplot(CD57_TBX21_hits_network_bubble, aes(y = group, x = reorder(node, LFC), size = LFC, color = `-log10(p-value)`)) +
  geom_point() + theme_minimal() + labs(x = element_blank(), y = "DEG Gene Name", title = "CD57+ Tex DEGs\nwith TBX21 Motifs") + scale_color_viridis_c() + 
  theme(text = element_text(size = 12, face = "bold"),  
        #legend.position = "bottom", 
        legend.box = "vertical", 
        plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90)) 

PD1_TCF7_hits_network_bubble_plot <- ggplot(PD1_TCF7_hits_network_bubble, aes(y = group, x = reorder(node, LFC), size = LFC, color = `-log10(p-value)`)) +
  geom_point() + theme_minimal() + labs(x = element_blank(), y = "DEG Gene Name", title = "PD-1+ Tex DEGs\nwith TCF7 Motifs") + scale_color_viridis_c()+ 
  theme(text = element_text(size = 12, face = "bold") ,
        #legend.position = "bottom", 
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90)) 

bubble_plots <- egg::ggarrange(CD57_TBX21_hits_network_bubble_plot,PD1_TCF7_hits_network_bubble_plot , ncol = 2, nrow = 1)
ggsave(bubble_plots, file = file.path(plotDir, "TF_DEG_bubble_plots.pdf"), height = 5, width = 17)

#### PLOT RESULTS AS PROTEIN-PROTEIN INTERACTION NETWORK FOR GENES WITH TF MOTIF ####

# export tables with the TF motif data and DEG data combined

write.csv(CD57_TBX21_hits_network_bubble, file = file.path(resultDir, "CD57_TBX21_hits_network_bubble.csv"), row.names = FALSE)
write.csv(PD1_TCF7_hits_network_bubble, file = file.path(resultDir, "PD1_TCF7_hits_network_bubble.csv"), row.names = FALSE)

# Plugging gene list into stringDB to get protein protein interaction network

### Load short networks with only 1-1 connections included
CD57_TBX21_string_interactions <- read.table(file = file.path(resultDir, "CD57_TBX21_motif_stringDB_string_interactions_string_interactions_short.tsv"), sep = "\t", header = TRUE)
PD1_TCF7_string_interactions <- read.table(file = file.path(resultDir, "PD1_TCF7_motif_stringDB_string_interactions_short.tsv"), sep = "\t", header = TRUE)

# create edge df 
CD57_TBX21_string_interactions_edge <- CD57_TBX21_string_interactions  %>%
  dplyr::select(node1,node2, combined_score) 

PD1_TCF7_string_interactions_edge <- PD1_TCF7_string_interactions  %>%
  dplyr::select(node1,node2, combined_score) 

# Create node files
CD57_TBX21_string_interactions_nodes <- 
  data.frame(node = unique(c(CD57_TBX21_string_interactions$node1, CD57_TBX21_string_interactions$node2)))

PD1_TCF7_string_interactions_nodes <- 
  data.frame(node = unique(c(PD1_TCF7_string_interactions$node1, PD1_TCF7_string_interactions$node2)))


# create network
CD57_TBX21_string_interactions_net <- graph_from_data_frame(d=CD57_TBX21_string_interactions_edge , 
                                                            vertices=CD57_TBX21_string_interactions_nodes, directed=F) 

PD1_TCF7_string_interactions_net <- graph_from_data_frame(d=PD1_TCF7_string_interactions_edge , 
                                                          vertices=PD1_TCF7_string_interactions_nodes, directed=F) 

# plot networks
# view igraph plotting parameters
#?igraph.plotting
# Set node size
V(CD57_TBX21_string_interactions_net)$size <- 7
V(PD1_TCF7_string_interactions_net)$size <- 7

# Set edge width based on weight:
E(CD57_TBX21_string_interactions_net)$width <- E(CD57_TBX21_string_interactions_net)$combined_score*5
E(PD1_TCF7_string_interactions_net)$width <- E(PD1_TCF7_string_interactions_net)$combined_score*5

# Set edge color 
#E(upset_pairwise_net)$color = E(upset_pairwise_net)$edge_color

# plot
pdf(file = file.path(plotDir, "CD57_TBX21_string_interactions_network.pdf"), width = 5, height =5 )
plot(CD57_TBX21_string_interactions_net, main = "CD57+ Tex DEGs with TBX21 Motifs\nPPI enrichment P-value = 4.38e-05",  
     label.font  = 2, vertex.label.color="black", vertex.label.cex = 1, 
     vertex.label.dist=1.5,
     vertex.color = "#93a24eff")
dev.off()

pdf(file = file.path(plotDir, "PD1_TCF7_string_interactions_network.pdf"), width = 6, height = 6 )
plot(PD1_TCF7_string_interactions_net, main = "PD-1+ Tex DEGs with TCF7 Motifs\nPPI enrichment P-value = 1.11e-06",  
     label.font  = 2, vertex.label.color="black", vertex.label.cex = 0.8, 
     vertex.label.dist=1.5,
     vertex.color = "#ba4d4cc7")
dev.off()


