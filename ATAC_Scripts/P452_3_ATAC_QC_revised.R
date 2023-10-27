#### P452-3 T1DAL CD57 Exhausted Cell ATAC Re-analysis ####

# P452-3 used the same samples as P452-1, but re-pooled and re-quantitated
# the samples prior to sequencing in order to decrease the variability
# in final sequencing read depth 

### Erin Witkop
### Bioinformatics Postdoctoral Research Associate

#### LOAD LIBRARIES ####

library(TxDb.Hsapiens.UCSC.hg38.knownGene) # human transcript
library(stringr) # string processing
library(tidyverse)
library(readxl)
library(soGGi) # for calcuation of TSSe scores
library(RColorBrewer)
library(ggplot2); # plotting
# set theme
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1),
          axis.text=element_text(colour="black"),
          axis.ticks=element_line(colour="black"),
          legend.key = element_blank(),
          text = element_text(size=10),
          strip.text.x = element_text(size = 10,margin = margin( b = 2, t = 2) ),
          strip.background = element_rect(fill="white", colour="black"))) +
  scale_color_brewer(palette = "Dark2")
library(ATACseqQC)
library(GenomicAlignments)

Palette = c('#DAE6F2', "#90F0B0", "#B7F794", "#E0DD91", "#FAE3A0", "#F0C69C", "#F5B9BA")

# setwd
setwd("/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES")

#### SET PATHS ####
projectname = "P452-3"
annoFileName = "P452-3 Repeat Final Annotation.xlsx"

# Paths 
baseDir <- "/Users/ewitkop/Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES"
plotDir <- file.path(baseDir,'FIGURES/ATAC_REPEAT_FIGURES')
resultDir <- file.path(baseDir,'P452_3_SAVED_DATA')
annotationDir <- file.path(baseDir, "RAW_DATA/P452_3")

fragmentDir <- file.path("/Volumes/Bioinformatics/pipeline/Illumina/220802_VH00126_203_AAAMVKNHV/Project_P452-3Processed_globus_220807/insertSizes")
alignmentDir <- file.path("/Volumes/Bioinformatics/pipeline/Illumina/220802_VH00126_203_AAAMVKNHV/Project_P452-3Processed_globus_220807/shiftedAlignments")

annoFile <- file.path(annotationDir,annoFileName)

#### LOAD ANNOTATION AND METRICS DATA ####

# Load annotation compiled by the Genomics Core
P452_3_anno <- read_excel(annoFile) %>% 
  # remove duplicated columns in the Final_annotation sheet
  select(-c( "Date collected...16", "index_name...21","index_name...23", "Date collected...26" ))

# reformat columns names
colnames(P452_3_anno) <- str_replace_all(colnames(P452_3_anno), "-","_") 
colnames(P452_3_anno) <- str_replace_all(colnames(P452_3_anno), " ","_") 
colnames(P452_3_anno) <- str_remove_all(colnames(P452_3_anno), "\\...*") 
colnames(P452_3_anno)
head(P452_3_anno)
# rename New_Lib_ID to match metrics data
colnames(P452_3_anno)[19] <- "libId"

# Save anno data
save(P452_3_anno, file =  file.path(resultDir, "P452_3_anno.Rdata"))

### Load Compiled Metrics data from Stephan 

post_alignment_metrics <- read_csv(file.path(annotationDir, "P452-3_AAAMVKNHV_220807_combined_metrics.csv"))
overrep_seq <- read_csv(file.path(annotationDir, "P452-3_AAAMVKNHV_220807_combined_overrep_seqs.csv"))
combined_qc <- read_csv(file.path(annotationDir, "P452-3_AAAMVKNHV_220807_combined_qc.csv"))
summary_data <- read_csv(file.path(annotationDir, "P452-3_AAAMVKNHV_220807_combined_summary-data.csv")) %>%
  # Separate libId and flow cell
  separate(libId, into= c("libId","Flow_Cell"), sep = "_")

## Join annotation data to summary_data
all(summary_data$libId %in% P452_3_anno$libId) # TRUE

summary_data_anno <- left_join(summary_data, P452_3_anno[,c("Sample_Name","Sort","Timepoint","libId","Donor_ID", "Treatment", "Cell_Number")])
summary_data_anno$Sample_Name <- factor(summary_data_anno$Sample_Name, levels = c("1_10398_BL_CD57pos_DP_NN_CD8","2_10398_BL_CD57neg_DP_NN_CD8","3_10398_BL_DN_NN_CD8",
                                                                                  "4_10398_V30_CD57pos_DP_NN_CD8","5_10398_V30_CD57neg_DP_NN_CD8","6_10398_V30_DN_NN_CD8",
                                                                                  "7_10213_BL_CD57pos_DP_NN_CD8","8_10213_BL_CD57neg_DP_NN_CD8","9_10213_BL_DN_NN_CD8",
                                                                                  "10_10213_V30_CD57pos_DP_NN_CD8","11_10213_V30_CD57neg_DP_NN_CD8","12_10213_V30_DN_NN_CD8",
                                                                                  "13_10458_BL_CD57pos_DP_NN_CD8","14_10458_BL_CD57neg_DP_NN_CD8","15_10458_BL_DN_NN_CD8",
                                                                                  "16_10458_V30_CD57pos_DP_NN_CD8","17_10458_V30_CD57neg_DP_NN_CD8","18_10458_V30_DN_NN_CD8",
                                                                                  "19_10295_BL_CD57pos_DP_NN_CD8","20_10295_BL_CD57neg_DP_NN_CD8","21_10295_BL_DN_NN_CD8",
                                                                                  "22_10295_V30_CD57pos_DP_NN_CD8","23_10295_V30_CD57neg_DP_NN_CD8","24_10295_V30_DN_NN_CD8"))
# Shorten the sort name
summary_data_anno <- summary_data_anno %>%
  mutate(Sort_short = case_when(
    Sort == "CD57- DP Non-naïve CD8 T cells: singlet, live, CD14-CD19-CD56-CD3+CD8+CD4-, not CD45RA+CCR7+, KLRG1+TIGIT+, CD57-" ~ "CD57-_DP",
    Sort == "CD57+ DP Non-naïve CD8 T cells: singlet, live, CD14-CD19-CD56-CD3+CD8+CD4-, not CD45RA+CCR7+, KLRG1+TIGIT+, CD57+" ~ "CD57+_DP",
    Sort == "DN Non-naïve CD8 T cells: singlet, live, CD14-CD19-CD56-CD3+CD8+CD4-, not CD45RA+CCR7+, KLRG1-TIGIT-" ~ "DN_non_naive"
  ))

summary_data_anno$Sort_short <- factor(summary_data_anno$Sort_short, levels = c("CD57-_DP","CD57+_DP","DN_non_naive"),
                                       labels = c("CD57- DP","CD57+ DP", "DN Non-naive"))

# rename %gc for downstream plotting and make into decimal, similar to other metircs
summary_data_anno <- mutate(summary_data_anno, percent_gc = `%gc`*0.01) 

# Make DONOR ID a factor
unique(summary_data_anno$Donor_ID)
summary_data_anno$Donor_ID  <- factor(summary_data_anno$Donor_ID, levels = c("10213" ,"10295", "10398" ,"10458"))

# Join library ID to overrep_seq
overrep_seq_anno <- overrep_seq %>%
  # Separate libId and flow cell
  separate(libId, into= c("libId","Flow_Cell"), sep = "_") %>% 
  left_join(., unique(summary_data_anno[,c("Sample_Name","Sort_short","Timepoint","libId","Donor_ID")]))

# Make DONOR ID a factor
overrep_seq_anno$Donor_ID  <- factor(overrep_seq_anno$Donor_ID, levels = c("10213" ,"10295", "10398" ,"10458"))

# save summary annotation data
 save(summary_data_anno, file = file.path(resultDir, "summary_data_anno.Rdata"))
load( file = file.path(resultDir, "summary_data_anno.Rdata"))

#### PLOT PRE-ALIGNMENT METRICS ####

# Ensure metrics of interest are numeric
class(summary_data_anno$fastq_total_reads)
class(summary_data_anno$median_cv_coverage)
class(summary_data_anno$median_3prime_bias)
class(summary_data_anno$median_5prime_bias)
class(summary_data_anno$percent_gc)
class(summary_data_anno$pct_correct_strand_reads)
class(summary_data_anno$percent_duplication)

### Make list of pre-QC metrics I'm interested in
pre_align_metrics_list <- c("fastq_total_reads",  "mean_read_length")
pre_align_perc_metrics <- c("percent_gc","pct_correct_strand_reads","percent_duplication")                  

### Write function for plotting pre-alignment metrics
# Design functions for plotting metrics and percentages separately
pre_align_ggplot <- function(data, metric){
  ggplot(data, aes_string(x = "Sample_Name", y = metric, fill = "Donor_ID")) + 
    geom_col(position = "dodge") + 
    facet_grid(.~Sort_short+Timepoint, scales = "free", space = "free", drop = TRUE) + 
    theme(axis.text.x = element_text(angle  = 90)) + 
    labs(y = metric, x = "Sample Name") + 
    scale_color_manual(values = Palette) + 
    scale_fill_discrete(name = "Donor ID")
}

pre_align_ggplot_perc <- function(data, metric){
  ggplot(data, aes_string(x = "Sample_Name", y = metric, fill = "Donor_ID")) + 
    geom_col(position = "dodge") + 
    facet_grid(.~Sort_short+Timepoint, scales = "free", space = "free", drop = TRUE) + 
    theme(axis.text.x = element_text(angle  = 90)) + 
    labs(y = metric, x = "Sample Name") + 
    scale_y_continuous(limits = c(0,1), labels = function(x) paste0(x*100, "%")) +
    scale_color_manual(values = Palette) + 
    scale_fill_discrete(name = "Donor ID")
}

# Apply the functions
pre_align_plot_list <- lapply(pre_align_metrics_list, pre_align_ggplot, data = summary_data_anno)
pre_align_perc_plot_list <- lapply(pre_align_perc_metrics, pre_align_ggplot_perc, data = summary_data_anno)

# Output each list into a single pdf
pdf(file.path(plotDir, "pre_align_qc_plots.pdf"), height = 8, width = 10)
invisible(lapply(pre_align_plot_list, print))
dev.off()

pdf(file.path(plotDir, "pre_align_qc_plots_perc.pdf"), height = 8, width = 10)
invisible(lapply(pre_align_perc_plot_list, print))
dev.off()

### Plot the percentage of over-represented sequences 
over_rep_perc_plot <-
  ggplot(overrep_seq_anno, aes(x = Sample_Name, y = percentage, fill = Donor_ID)) + 
  geom_col(position = "dodge") + 
  facet_grid(.~Sort_short+Timepoint, scales = "free", space = "free", drop = TRUE) + 
  theme(axis.text.x = element_text(angle  = 90)) + 
  labs(y = "Percentage Over-Represented Sequences", x = "Sample Name") + 
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = Palette) + 
  scale_fill_discrete(name = "Donor ID")

ggsave(over_rep_perc_plot, file = file.path(plotDir, "over_rep_perc_plot.jpg"),device = "jpg",height = 8, width = 10)

#### PLOT THE PERCENT OF MITOCHONDRIAL READS IN EACH ALIGNMENT ####

percent_mt_reads <- read.table(file.path(annotationDir, "percent_mt_reads_shifted_sorted_full.txt"), col.names = c("path","percent_mt_reads"))

# join libID and Sample_Name
percent_mt_reads$libId <- str_extract(percent_mt_reads$path, "lib[0-9]+")
percent_mt_reads <- left_join(percent_mt_reads, anno_unique)

# Plot 
percent_mt_reads_plot <- ggplot(percent_mt_reads, aes(x = Sample_Name, y =  percent_mt_reads, fill = Donor_ID)) + 
  geom_col(position = "dodge") + 
  facet_grid(.~Sort_short+Timepoint, scales = "free", space = "free", drop = TRUE) + 
  theme(axis.text.x = element_text(angle  = 90)) + 
  labs(y = "Percent MT Reads", x = "Sample Name") + 
  scale_y_continuous(limits = c(0,100))+
  scale_color_manual(values = Palette) + 
  scale_fill_discrete(name = "Donor ID")

ggsave(percent_mt_reads_plot, file= file.path(plotDir, "percent_mt_reads_plot.jpg"), device = "jpg", height = 5, width = 7)


#### PLOT POST-ALIGNMENT METRICS ####

post_align_metrics_list <- c("pf_hq_aligned_reads", "pf_reads_aligned", "median_cv_coverage", "median_3prime_bias","median_5prime_bias")
post_align_perc_metrics <- c("pct_pf_reads_aligned")                  

### Write function for plotting post-alignment metrics
# Design functions for plotting metrics and percentages separately
post_align_ggplot <- function(data, metric){
  ggplot(data, aes_string(x = "Sample_Name", y = metric, fill = "Donor_ID")) + 
    geom_col(position = "dodge") + 
    facet_grid(.~Sort_short+Timepoint, scales = "free", space = "free", drop = TRUE) + 
    theme(axis.text.x = element_text(angle  = 90)) + 
    labs(y = metric, x = "Sample Name") + 
    scale_color_manual(values = Palette) + 
    scale_fill_discrete(name = "Donor ID")
}

post_align_ggplot_perc <- function(data, metric){
  ggplot(data, aes_string(x = "Sample_Name", y = metric, fill = "Donor_ID")) + 
    geom_col(position = "dodge") + 
    facet_grid(.~Sort_short+Timepoint, scales = "free", space = "free", drop = TRUE) + 
    theme(axis.text.x = element_text(angle  = 90)) + 
    labs(y = metric, x = "Sample Name") + 
    scale_y_continuous(limits = c(0,1), labels = function(x) paste0(x*100, "%")) +
    scale_color_manual(values = Palette) + 
    scale_fill_discrete(name = "Donor ID")
}

# Apply the functions
post_align_plot_list <- lapply(post_align_metrics_list, post_align_ggplot, data = summary_data_anno)
post_align_perc_plot_list <- lapply(post_align_perc_metrics, post_align_ggplot_perc, data = summary_data_anno)

# Output each list into a single pdf
pdf(file.path(plotDir, "post_align_qc_plots.pdf"), height = 8, width = 10)
invisible(lapply(post_align_plot_list, print))
dev.off()

pdf(file.path(plotDir, "post_align_qc_plots_perc.pdf"), height = 8, width = 10)
invisible(lapply(post_align_perc_plot_list, print))
dev.off()

#### PLOT READ FRAGMENT LENGTH DISTRIBUTION ####

#Get a list of frag size files
fragFiles <- list.files(fragmentDir, 
                        full.names = T)

# make unique df of anno data
anno_unique <- summary_data_anno %>% dplyr::distinct(libId, Sample_Name, Donor_ID, Timepoint, Sort_short, Treatment)

#Function to read and process files
unpackFragData <- function(fragFileIn, annotation){
  libId <- str_extract(fragFileIn, "lib[0-9]+")
  Sample_Name <- annotation$Sample_Name[annotation$libId == libId]
  Treatment <- annotation$Treatment[annotation$libId == libId]
  Donor_ID <- annotation$Donor_ID[annotation$libId == libId]
  Sort_short <- annotation$Sort_short[annotation$libId == libId]
  Timepoint <- annotation$Timepoint[annotation$libId == libId]
  
  
  fragDist <- read.table(fragFileIn)
  colnames(fragDist) <- c("nReads", "fragLen")
  
  fragDist$fragLen <- abs(fragDist$fragLen)
  fragDist$normReads <- fragDist$nReads/sum(fragDist$nReads) *10^3
  fragDist$libId <- libId
  fragDist$Donor_ID <- Donor_ID
  fragDist$Sample_Name <-Sample_Name
  fragDist$Treatment <- Treatment
  fragDist$Sort_short <- Sort_short
  fragDist$Timepoint <- Timepoint
  
  
  return(fragDist)
  
}

#  This code takes a long time!
#fragList <- lapply(fragFiles, function(x) unpackFragData(x, anno_unique))
#fragComb <- bind_rows(fragList)
head(fragComb)

# save
save(fragComb, file = file.path(resultDir, "fragComb.Rdata"))
load(file = file.path(resultDir, "fragComb.Rdata"))

# Plot by Timepoint and cell sort 
gFragComb  <- fragComb %>%
  dplyr::filter(fragLen > 0 & fragLen < 1000) %>%
  ggplot(aes(x = fragLen,
             y = normReads/max(normReads),
             color = Donor_ID))+
  facet_grid(.~Sort_short + Timepoint) +
  geom_line()+
  xlim(c(0,1000)) +
  labs(x = "Fragment length, bp", y = expression(paste("Normalized read density (a.u.)")), color = "") + 
  theme(text = element_text(size = 12)) +
  theme(legend.position="top") + 
  theme(legend.direction='horizontal')

ggsave(gFragComb, file= file.path(plotDir, "FragmentSizeDistCombined_sort_time.jpg"), device = "jpg", height = 5, width = 10)

# Plot by Individual lib to check each sample
gFraglib <- fragComb %>%
  dplyr::filter(fragLen > 0 & fragLen < 1000) %>%
  ggplot(aes(x = fragLen,
             y = normReads/max(normReads),
             color = Donor_ID))+
  facet_grid(.~Sample_Name) +
  geom_line()+
  xlim(c(0,1000)) +
  labs(x = "Fragment length, bp", y = expression(paste("Normalized read density (a.u.)")), color = "") + 
  theme(text = element_text(size = 12)) +
  theme(legend.position="top") + 
  theme(legend.direction='horizontal')

ggsave(gFraglib, file= file.path(plotDir, "FragmentSizeDistCombined_sample.jpg"), device = "jpg", height = 5, width = 20)

#### PLOT ENRICHMENT AROUND THE TF START SITE ####

# Run P452_3_TSSe_biometal.R on a server to calculate transcription factor start site enrichment 

## Load the data once it has finished
load(file.path(resultDir, "TSSe_output.Rdata"))

TSSEComb

# Join with annotation data
TSSEComb_anno <- left_join(TSSEComb, anno_unique)
colnames(TSSEComb_anno)

# Plot TSSe 
TSSEplot_P452_3  <-
  ggplot(TSSEComb_anno, aes(x = distanceToTSS-25,y = aggregateTSSscore, color = Sample_Name))+
  geom_line()+
  labs(x = "distance to TSS (bp)", y = "aggregate TSS score", color = "") + 
  theme(text = element_text(size = 13)) +
  theme(legend.position="right") + 
  theme(legend.direction='vertical') +
  facet_wrap(. ~ Sort_short + Timepoint) +
  scale_x_continuous(limits = c(-900, 900))

ggsave(TSSEplot_P452_3 , file = file.path(plotDir, "TSSEplot_P452_3_ALL_LIBRARIES.jpg"), device = "jpg", height = 5, width = 15)



