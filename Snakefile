
# Configurations
configfile: "config.yaml"

# Define the samples
SAMPLES= config["samples"]
GENOME_DB = config["genome"]
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
GTF = config["gtf"]

rule all:
	input:
		expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_1_trimmed.fq.gz", sample=SAMPLES)
		expand(f"{OUTPUT_DIR}/aligned/{{sample}}_CriGri_sorted.bam", sample=SAMPLES)
		expand(f"{OUTPUT_DIR}/count/{{sample}}_CriGri_count.txt", sample=SAMPLES)
		f"{OUTPUT_DIR}/count/all_counts.txt"
		f"{OUTPUT_DIR}/count/final_all_counts.txt"
		f"{OUTPUT_DIR}/../output/dea_results.rds"
		f"{OUTPUT_DIR}/../output/vst_transformed_counts_long.rds"
		f"{OUTPUT_DIR}/../output/volcano.png"
		f"{OUTPUT_DIR}/../output/jitter.png"

# Rule: Trim reads using Trim Galore
rule trim:
	input:
		r1=f"{INPUT_DIR}/{{sample}}_1.fastq.gz", # Raw read forward
		r2=f"{INPUT_DIR}/{{sample}}_2.fastq.gz"  # Raw read reverse
	output:
		dir=f"{OUTPUT_DIR}/trimmed"
		r1_trimmed=f"{OUTPUT_DIR}/trimmed/{{sample}}_1_trimmed.fq.gz", # Trimmed read forward
		r2_trimmed=f"{OUTPUT_DIR}/trimmed/{{sample}}_2_trimmed.fq.gz"  # Trimmed read reverse
	log:
		f"logs/{{sample}}_trim.log"
	shell:
		"""
		trim_galore -q 25 --phred33 --length 20 --stringency 3 --paired \
		--output_dir {output.dir} {input.r1} {input.r2} > {log} 2>&1
		mv data/processed/trimmed/{wildcards.sample}_1_val_1.fq.gz {output.r1_trimmed}
		mv data/processed/trimmed/{wildcards.sample}_2_val_2.fq.gz {output.r2_trimmed}
		"""

# Rule: Align reads using HISAT2 and sort using samtools
rule align:
	input:
		r1_trimmed=f"{OUTPUT_DIR}/trimmed/{{sample}}_1_trimmed.fq.gz", # Trimmed read forward
		r2_trimmed=f"{OUTPUT_DIR}/trimmed/{{sample}}_2_trimmed.fq.gz"  # Trimmed read reverse
	output:
		bam=f"{OUTPUT_DIR}/aligned/{{sample}}.bam"	
	log:
		f"logs/{{sample}}_align.log"
	threads: 4
	shell:
		"""
		hisat2 -p {threads} -x {GENOME_DB} -1 {input.r1_trimmed} -2 {input.r2_trimmed} | \
		samtools sort -O BAM -o {output.bam} > {log} 2>&1
		"""

# Rule: Index BAM files using samtools
rule index:
	input:
		bam=f"{OUTPUT_DIR}/aligned/{{sample}}.bam"
	output:
		bai=f"{OUTPUT_DIR}/aligned/{{sample}}.bam.bai"
	log:
		f"logs/{{sample}}_index.log"
	shell:
		"samtools index {input.bam} > {log} 2>&1"

# Rule: Count reads using featureCounts; merge into one summary file
rule count:
	input:
		bam=expand(f"{OUTPUT_DIR}/aligned/{{sample}}.bam", sample=SAMPLES)
	output:
		count=f"{OUTPUT_DIR}/count/all_counts.txt"
	log:
		f"logs/all_counts.log"
	shell:
		"""
		featureCounts -T 4 -p -t exon -g gene_id -a {GTF} \
		-o {output.count} {input.bam} > {log} 2>&1
		"""
		
# Rule: Convert the summary file into the correct format
rule process_counts:
	input:
		counts=f"{OUTPUT_DIR}/count/all_counts.txt"
	output:
		final_counts=f"{OUTPUT_DIR}/count/final_all_counts.txt"
	shell:
		"tail -n +2 {input.counts} | cut -f1,7- > {output.final_counts}"

# Rule: Differential expression analysis
rule run_analysis:
	input:
		f"{OUTPUT_DIR}/count/final_all_counts.txt"
	output:
		f"{OUTPUT_DIR}/../output/dea_results.rds"
		f"{OUTPUT_DIR}/../output/vst_transformed_counts_long.rds"
	script:
		f"{OUTPUT_DIR}/../../scripts/04_analysis.R"
	
# Rule: Visualization
rule run_visualization:
	input:
		f"{OUTPUT_DIR}/../output/dea_results.rds"
		f"{OUTPUT_DIR}/../output/vst_transformed_counts_long.rds"
	output:
		f"{OUTPUT_DIR}/../output/volcano.png"
		f"{OUTPUT_DIR}/../output/jitter.png"
	script:
		f"{OUTPUT_DIR}/../../scripts/05_visualization.R"
