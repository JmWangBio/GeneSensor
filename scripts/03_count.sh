#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Specify paths to experiment name, input, output, and genome
EXP_NAME='gene_sensor'
INPUT_DIR='/path/to/data/processed/aligned'
OUTPUT_DIR='/path/to/data/processed/count'
GTF='/path/to/gtf/CriGri.gtf.gz'

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Collect all BAM files into a single variable
BAM_FILES=$(ls ${INPUT_DIR}/${EXP_NAME}-*_CriGri_sorted.bam)

# Echo progress message
echo "Processing featureCounts for all samples..."

# Run featureCounts on all BAM files simultaneously
featureCounts -T 4 -p -t exon -g gene_id -a $GTF \
  -o ${OUTPUT_DIR}/${EXP_NAME}_CriGri_counts.txt $BAM_FILES

# Clean up the output to retain only counts
tail -n +2 ${OUTPUT_DIR}/${EXP_NAME}_CriGri_counts.txt | \
  cut -f1,7- > ${OUTPUT_DIR}/${EXP_NAME}_CriGri_final_counts.txt

# Echo completion message
echo "FeatureCounts completed successfully!"
