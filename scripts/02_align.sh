#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Specify paths to experiment name, input, output, and genome
EXP_NAME='gene_sensor'
INPUT_DIR='/path/to/data/processed/trimmed'
OUTPUT_DIR='/path/to/data/processed/aligned'
GENOME_DB='/path/to/fasta/CriGri'

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Loop over all samples
for i in {1..9}
do
  echo "Aligning sample ${EXP_NAME}-${i}..."
  
  # Align reads with HISAT2 and sort the output
  hisat2 -t -p 4 -x $GENOME_DB \
    -1 ${INPUT_DIR}/${EXP_NAME}-${i}_trimmed_1.fq.gz \
    -2 ${INPUT_DIR}/${EXP_NAME}-${i}_trimmed_2.fq.gz | \
    samtools sort -O BAM -o ${OUTPUT_DIR}/${EXP_NAME}-${i}_CriGri_sorted.bam
    
  # Collect alignment statistics 
  samtools flagstat ${OUTPUT_DIR}/${EXP_NAME}-${i}_CriGri_sorted.bam > \
    ${OUTPUT_DIR}/${EXP_NAME}-${i}_CriGri.flagstat
  
  # Index the BAM file
  samtools index ${OUTPUT_DIR}/${EXP_NAME}-${i}_CriGri_sorted.bam
  
  echo "Sample ${EXP_NAME}-${i} aligned successfully!"
done
