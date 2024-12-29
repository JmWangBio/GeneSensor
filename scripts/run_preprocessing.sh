#!/bin/bash
# Master preprocessing script

echo "Starting RNA-seq preprocessing pipeline..."

# Step 1: Trimming
echo "Trimming reads..."
bash 01_trim.sh

# Step 2: Alignment
echo "Aligning reads..."
bash 02_align.sh

# Step 3: Quantification
echo "Counting reads..."
bash 03_count.sh

echo "Preprocessing complete!"
