#!/bin/bash

# Exit immediately if a command fails
set -e

# Define the paths to the Snakefile and config file
SNAKEFILE="Snakefile"
CONFIG="config.yaml"
LOG_DIR="logs"

# Create the log directory if it does not exist
mkdir -p ${LOG_DIR}

# Execute Snakemake
snakemake --snakefile ${SNAKEFILE} \
	--configfile ${CONFIG} \
	--cores 4 \
	--keep-going \
	--rerun-incomplete \
	--printshellcmds \
	> ${LOG_DIR}/workflow.log 2>&1

# Print success message
echo "Workflow executed successfully. Check ${LOG_DIR}/workflow.log for details."
