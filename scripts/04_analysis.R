# RNA-seq Data Analysis Script
# Author: Junmin Wang
# Description: This script performs normalization, differential expression analysis,
# gene annotation, and variance stabilization for RNA-seq data.

# ===============================================
# Load Required Libraries
# ===============================================
required_packages <- c("tidyverse", "biomaRt", "edgeR", "limma", "DESeq2")
lapply(required_packages, library, character.only = TRUE)

# ===============================================
# Load Raw Count Data
# ===============================================
cat("Loading raw count data...\n")
data_wide <- read.delim("/path/to/data/processed/count/gene_sensor_CriGri_final_counts.txt", 
                        sep = "\t")

# Rename columns and set rownames to Gene IDs
data_wide <- data_wide %>% 
  rename_with(~ gsub("_CriGri_.*.", "", .)) %>%
  column_to_rownames(var = "Geneid")

# Convert to matrix for edgeR compatibility
data_wide_mat <- as.matrix(data_wide)

# ===============================================
# Gene Annotation: Convert Ensembl IDs to Gene Symbols
# ===============================================
cat("Fetching gene annotations from Ensembl...\n")
hamster <- useMart(
  'ensembl', 
  dataset = 'cgcrigri_gene_ensembl',
  host = "nov2020.archive.ensembl.org"
)

# Retrieve Ensembl ID to gene symbol mapping
gene_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  mart = hamster
)

# ===============================================
# Create DGEList Object and Normalize Counts
# ===============================================
cat("Creating DGEList object and normalizing counts...\n")
d0 <- DGEList(counts = data_wide_mat)

# Normalize with TMM method (default)
d <- calcNormFactors(d0, method = "TMM")

# ===============================================
# Define Metadata: Experimental Groups
# ===============================================
group <- rep(c("WT", "Bcl2.const", "Bcl2.gs"), each = 3)

# ===============================================
# Define Custom Functions for Limma-Voom Analysis
# ===============================================
ContrastsFit <- function(Data, contrasts.pair, Group) {
  if (!inherits(Data, "DGEList")) stop("Data must be a DGEList object.")
  ct <- factor(Group)
  design <- model.matrix(~ 0 + ct)
  colnames(design) <- levels(ct)
  keep <- filterByExpr(Data, design)
  y <- voom(Data[keep, ], design, plot = FALSE)
  fit <- lmFit(y, design)
  contrasts <- makeContrasts(contrasts = contrasts.pair, levels = design)
  contrasts.fit <- eBayes(contrasts.fit(fit, contrasts))
  return(list(trData = y, fit = contrasts.fit))
}

LimmaResult <- function(contrasts.fit, contrasts.pair, contrasts.name) {
  results <- purrr::map2(
    seq_along(contrasts.pair),
    contrasts.pair,
    ~ {
      r <- topTable(contrasts.fit, coef = .x, number = Inf, adjust.method = "BH") %>%
        filter(!is.na(P.Value)) %>%
        rownames_to_column(var = "ID") %>%
        select(ID, logFC, t, P.Value, adj.P.Val) %>%
        rename_with(~ paste(contrasts.name[.x], ., sep = "_"), -ID)
      return(r)
    }
  )
  r_combine <- purrr::reduce(results, full_join, by = "ID")  # Combine results
  return(r_combine)
}

# ===============================================
# Perform Limma-Voom Differential Expression Analysis
# ===============================================
cat("Performing differential expression analysis...\n")

# Define contrasts for comparisons
contrasts.pair <- c("Bcl2.const-WT", "Bcl2.gs-WT", "Bcl2.gs-Bcl2.const")
contrasts.name <- gsub("-", "vs", contrasts.pair)

# Fit model and compute contrasts
contrasts.fit.res <- ContrastsFit(d, contrasts.pair, group)
contrasts.fit <- contrasts.fit.res[["fit"]]

# Extract results
LIMMAoutput <- LimmaResult(contrasts.fit, contrasts.pair, contrasts.name)

# ===============================================
# Combine Results with Gene Annotations
# ===============================================
cat("Merging differential expression results with gene annotations...\n")
Result <- LIMMAoutput %>%
  left_join(gene_annotations, by = c("ID" = "ensembl_gene_id")) %>%
  rename(Gene = external_gene_name) %>%
  relocate(Gene, .after = ID) %>%
  mutate(Gene = ifelse(is.na(Gene), "Unknown", Gene))  # Replace missing symbols with "Unknown"

# ===============================================
# Save Results
# ===============================================
cat("Saving results...\n")
saveRDS(Result, file = "/path/to/data/output/dea_results.rds")

# ===============================================
# Perform VST Transformation with DESeq2
# ===============================================
# Create DESeq2 object
cat("Creating DESeq2 object and performing transformation...\n")
col_data <- data.frame(
  condition = group,
  row.names = colnames(data_wide_mat)
)
dds <- DESeqDataSetFromMatrix(
  countData = data_wide_mat,
  colData = col_data,
  design = ~condition
)

# Apply variance-stabilizing transformation (VST)
vst_data <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_data)

# Save transformed data for downstream visualization
vst_df <- as.data.frame(vst_mat) %>%
  rownames_to_column(var = "GeneID")

# Add gene symbols using the `table` from earlier
vst_df <- vst_df %>%
  left_join(table, by = c("GeneID" = "ensembl_gene_id")) %>%
  rename("external_gene_name" = "GeneSymbol") %>%
  relocate(GeneSymbol, .after = GeneID)

# Add group metadata
group <- rep(c("WT", "Bcl2.const", "Bcl2.gs"), each = 3)  # Define group info
colnames(vst_df)[3:ncol(vst_df)] <- paste0(colnames(vst_df)[3:ncol(vst_df)], "_", group)

# Convert to long format
vst_long <- vst_df %>%
  pivot_longer(
    cols = -c(GeneID, GeneSymbol),
    names_to = c("Sample", "Group"),
    names_sep = "_",  # Split sample and group info
    values_to = "Expression"
  )

saveRDS(vst_long, file = "/path/to/data/output/vst_transformed_counts_long.rds")
