#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Specify paths to experiment name, input, and output
EXP_NAME='gene_sensor'
INPUT_DIR='/path/to/data/raw'
OUTPUT_DIR='/path/to/data/processed/trimmed'

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Loop over all samples
for i in {1..9}
do
  echo "Trimming sample ${EXP_NAME}-${i}..."
  trim_galore -q 25 --phred33 --length 20 --stringency 3 --paired -o ${OUTPUT_DIR} \
    ${INPUT_DIR}/${EXP_NAME}-${i}_1.fastq.gz ${INPUT_DIR}/${EXP_NAME}-${i}_2.fastq.gz
  
  # Rename output files for clarity
  mv ${OUTPUT_DIR}/${EXP_NAME}-${i}_1_val_1.fq.gz ${OUTPUT_DIR}/${EXP_NAME}-${i}_trimmed_1.fq.gz
  mv ${OUTPUT_DIR}/${EXP_NAME}-${i}_2_val_2.fq.gz ${OUTPUT_DIR}/${EXP_NAME}-${i}_trimmed_2.fq.gz
    
  echo "Sample ${EXP_NAME}-${i} trimmed successfully!" 
done
