#### P452-3 CALCULATE TRANSCRIPTION START SITE ENRICHMENT ON BIOMETAL R STUDIO SERVER

# I previously ran this code on my local computer, but it takes forever to run over the 
# internet!

#### LOAD LIBRARIES ####

library(TxDb.Hsapiens.UCSC.hg38.knownGene) # human transcript
library(stringr) # string processing
library(tidyverse)
library(readxl)
library(soGGi) # for calcuation of TSSe scores
library(ATACseqQC)
library(GenomicAlignments)

#### SET PATHS ####
resultDir <- "/mnt/bioinformatics/workspace/ewitkop/P452_3_T1DAL_Repeat_ATACseq"

fragmentDir <- "/mnt/bioinformatics/pipeline/Illumina/220802_VH00126_203_AAAMVKNHV/Project_P452-3Processed_globus_220807/insertSizes"
alignmentDir <- "/mnt/bioinformatics/pipeline/Illumina/220802_VH00126_203_AAAMVKNHV/Project_P452-3Processed_globus_220807/shiftedAlignments"


#### Calculate enrichment around the transcription start site ####

# get a list of alignment files
alignmentFiles <- list.files(alignmentDir, full.names = T, pattern = "^.*\\_shifted_sorted.bam$")

# Get transcript and sequence information from human genome to get transcription start site info
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
seqlevelsStyle(txs) <- "NCBI" 
seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
which <- as(seqinformation, "GRanges")

TSSdistance_axis = seq(-950, 1000, by = 50)

chromosomes = paste0(seq(1,23))

# Function to get the TSSe score
getTSSEscore <- function(alignmentFile){
  
  libId <- str_extract(alignmentFile, "lib[0-9]+")
  print(sprintf("Processing %s ...", libId))
  bamfile <- alignmentFile
  # reading the GAlignments takes a long time
  gln <- readGAlignments(bamfile)
  tsse <- TSSEscore(gln, txs, seqlev = chromosomes, width = 50, step = 50)
  # tsse <- TSSEscore(gln, txs, seqlev = c("chr12"), width = 50, step = 50)
  rm(gln)
  
  TSSElib = data.frame(TSSdistance_axis, tsse$values)
  colnames(TSSElib) <- c("distanceToTSS", "aggregateTSSscore")
  TSSElib$libId = libId
  TSSElib$TSSEscore <- tsse$TSSEscore
  
  return(TSSElib)
}

## Run on full list of libraries 

TSSEList <- lapply(alignmentFiles, function(x) getTSSEscore(x)) 
TSSEComb <- bind_rows(TSSEList)

# save output  
save(TSSEList, TSSEComb,
     file = file.path(resultDir, "TSSe_output.Rdata"))