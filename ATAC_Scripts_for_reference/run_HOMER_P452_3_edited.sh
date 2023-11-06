
#!/bin/bash

# Script to Run HOMER transcription factor motif enrichment analysis
# see this helpful tutorial for running analysis: http://homer.ucsd.edu/homer/ngs/peakMotifs.html

# Note I already loaded hg38 using the following command
#loadGenome.pl -name hg38 -org null -fasta ~/Documents/Genomes/hg38.fa -gtf ~/Documents/Genomes/hg38.ensGene.gtf

# Set paths
DIR=/Users/ewitkop/Library/CloudStorage/Box-Box/EW_Bioinformatics_Postdoc_Research/T1DAL_10X_Project/ANALYSIS_FILES/ATAC_Analysis/HOMER_Analysis
TF=/Applications/homer/bin

# NOTES on parameters
# -size is mandatory and 200 is the default window
# removed -mask which was to allow use of the masked genome. I loaded the soft masked genome in the first place for hg38, so we are already set there

# run Find motifs on peaks separated out into the contrasts of interest
#findMotifsGenome.pl $DIR/CD57_vs_PD1_BED_UP.txt hg38 $DIR/CD57_vs_PD1_UP -size 200 
$TF/findMotifsGenome.pl $DIR/CD57_vs_PD1_BED_DOWN.txt hg38 $DIR/CD57_vs_PD1_DOWN -size 200 

$TF/findMotifsGenome.pl $DIR/CD57minus_vs_DN_BED_UP.txt hg38 $DIR/CD57minus_vs_DN_UP -size 200 
$TF/findMotifsGenome.pl $DIR/CD57pos_vs_DN_BED_UP.txt hg38 $DIR/CD57pos_vs_DN_UP -size 200 