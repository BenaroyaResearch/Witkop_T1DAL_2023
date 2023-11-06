#!/bin/bash

# THIS CODE PUTS BAM FILES INTO SCALED BIGWIG FORMAT FOR VIEWING IN UCSC GENOME BROWSER
INPUT_DIR=$1 #
# run this command from srvbiometal
OUTPUT_DIR=${INPUT_DIR}/bigWigs
mkdir -p ${OUTPUT_DIR}

eval "$(conda shell.bash hook)"
conda activate deeptools

INDEX=1

for currbam in `ls ${INPUT_DIR}/*shifted_sorted.bam`
do
    currLib=`echo ${currbam} | grep -Eo lib[0-9]+`

    echo "Generating scaled bigwig"
    bamCoverage -b ${currbam} -o $OUTPUT_DIR/${currLib}_scaled.bigWig --scaleFactor 1 --binSize 10 --extendReads --effectiveGenomeSize 2913022398
   # bamCoverage -b ${currbam} -o $OUTPUT_DIR/${currLib}_unscaled.bigWig --binSize 10 --extendReads --effectiveGenomeSize 2913022398

    let INDEX=${INDEX}+1
done