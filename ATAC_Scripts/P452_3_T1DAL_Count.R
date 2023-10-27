
#### P452-3 Run Counting for ATAC data in srvbiometal ####


#### LOAD LIBRARIES 
library(tidyverse)
library(DiffBind)

resultDir <- "/mnt/bioinformatics/workspace/ewitkop/P452_3_T1DAL_Repeat_ATACseq"

#### Pre-process ATAC data with correct paths ####

## Loading - loading in the peak set using dba - with all samples with under 30 million reads removed 
samplesheet <- read.csv(file=file.path(resultDir, "samplesheet_P452_3.csv"))
samplesheet$bamReads <- gsub("/Volumes/", "/mnt/", samplesheet$bamReads)
samplesheet$bamReads <- gsub("/Bioinformatics/", "/bioinformatics/", samplesheet$bamReads)
samplesheet$Peaks <- gsub("/Volumes/", "/mnt/", samplesheet$Peaks)
samplesheet$Peaks <- gsub("/Bioinformatics/", "/bioinformatics/", samplesheet$Peaks)

# Load data
ATACseqData_P452_3 <- dba(sampleSheet=samplesheet)
# this object has the initial peak set and gives metadata about the number of peaks 

## Counting (slow)
ATACseqData_P452_3 <- dba.count(ATACseqData_P452_3, 
                            bRemoveDuplicates = TRUE, 
                            bParallel = TRUE,
                            summit = 1000,
)
# computing summits 
# recentering peaks

# Save output - run rest of the pipeline in R locally 
save(ATACseqData_P452_3, file = file.path(resultDir, "ATACseqData_P452_3.RData"))
