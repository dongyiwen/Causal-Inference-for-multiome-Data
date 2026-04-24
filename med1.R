
suppressPackageStartupMessages({
  library(dplyr)
  library(pbapply)
  library(CMAverse)
})

# -----------------------------
# MEDIATION FUNCTION 
# -----------------------------

run_mediation_for_gene <- function(gene_id, seu, peak_to_gene, 
                                   exposure_var = "Sample.Age", 
                                   exposed_level = "pcw20",
                                   verbose = TRUE) {
  
  if (verbose) cat(sprintf("  Analyzing %s...\n", gene_id))
  
  if (!gene_id %in% rownames(seu[["RNA"]])) {
    if (verbose) warning(sprintf("Gene %s not found in RNA data", gene_id))
    return(NULL)
  }
  
  tryCatch({
    # Try Seurat v5 method first
    expr_vec <- as.numeric(LayerData(seu, assay = "RNA", layer = "data")[gene_id, ])
  }, error = function(e) {
    # Fall back to Seurat v4 method
    expr_vec <- as.numeric(seu[["RNA"]]@data[gene_id, ])
  })
  
  # Get peaks linked to this gene 
  gene_peaks <- peak_to_gene$peak[peak_to_gene$gene_id == gene_id]
  
  if (length(gene_peaks) == 0) {
    if (verbose) warning(sprintf("No promoter peaks found for %s", gene_id))
    return(NULL)
  }
  
  # Check which peaks actually exist in ATAC assay
  available_peaks <- rownames(seu[["ATAC"]])
  gene_peaks <- intersect(gene_peaks, available_peaks)
  
  if (length(gene_peaks) == 0) {
    if (verbose) warning(sprintf("No matching peaks in ATAC data for %s", gene_id))
    return(NULL)
  }
  
  # Get chromatin accessibility (Seurat v5 compatible)
  tryCatch({
    # Try Seurat v5 method
    atac_counts <- LayerData(seu, assay = "ATAC", layer = "counts")[gene_peaks, , drop = FALSE]
  }, error = function(e) {
    # Fall back to Seurat v4 method
    atac_counts <- seu[["ATAC"]]@counts[gene_peaks, , drop = FALSE]
  })
  
  chrom_vec <- colMeans(as.matrix(atac_counts))
  
  # developmental stage
  stage_vec <- ifelse(seu[[exposure_var]][, 1] %in% exposed_level, 1, 0)
  
  # Create mediation dataset
  med_data <- data.frame(
    stage = stage_vec,
    chrom_access = chrom_vec,
    expression = expr_vec
  )
  
  # Remove NAs
  med_data <- na.omit(med_data)
  
  # Check sample size
  if (nrow(med_data) < 30) {
    if (verbose) warning(sprintf("Insufficient samples for %s (n=%d)", gene_id, nrow(med_data)))
    return(NULL)
  }
  
  # Check for sufficient variation
  if (sd(med_data$chrom_access) < 1e-6 || sd(med_data$expression) < 1e-6) {
    if (verbose) warning(sprintf("No variation in data for %s", gene_id))
    return(NULL)
  }
  
  # Run causal mediation analysis
  result <- tryCatch({
    fit <- cmest(
      data = med_data,
      model = "rb",
      casecontrol = FALSE,
      outcome = "expression",
      exposure = "stage",
      mediator = "chrom_access",
      EMint = TRUE,
      mreg = list("linear"),
      yreg = "linear",
      a = 1,
      astar = 0,
      mval = list(0),
      estimation = "paramfunc",
      inference = "delta"
    )
    
    summary_fit <- summary(fit)
    
    list(
      gene_id = gene_id,
      gene_name = peak_to_gene$gene_name[peak_to_gene$gene_id == gene_id][1],
      n_peaks = length(gene_peaks),
      n_cells = nrow(med_data),
      fit = fit,
      summary = summary_fit
    )
    
  }, error = function(e) {
    if (verbose) warning(sprintf("Mediation failed for %s: %s", gene_id, e$message))
    NULL
  })
  
  return(result)
}

# -----------------------------
# FIX PEAK FORMAT MISMATCH
# -----------------------------

# Convert peak names from "chr1:123-456" to "chr1-123-456"
peak_to_gene$peak_fixed <- gsub(":", "-", peak_to_gene$peak)

# Verify match
n_match <- sum(peak_to_gene$peak_fixed %in% rownames(seu[["ATAC"]]))
cat(sprintf("Peaks matching after format fix: %d / %d (%.1f%%)\n",
            n_match, nrow(peak_to_gene),
            100 * n_match / nrow(peak_to_gene)))

# Replace peak column with fixed version
peak_to_gene$peak <- peak_to_gene$peak_fixed
peak_to_gene$peak_fixed <- NULL

# -----------------------------
# SELECT 20 GENES FOR ANALYSIS
# -----------------------------

cat("\n=== Selecting Genes for Mediation Analysis ===\n")

# Get genes with promoter peaks
genes_with_peaks_id <- unique(peak_to_gene$gene_id)
valid_genes <- intersect(genes_with_peaks_id, rownames(seu[["RNA"]]))

cat(sprintf("Found %d genes with promoter peaks in RNA data\n", length(valid_genes)))

# Select 20 genes with most peaks 
gene_peak_counts <- peak_to_gene %>%
  dplyr::filter(gene_id %in% valid_genes) %>%
  group_by(gene_id) %>%
  summarise(n_peaks = n(), .groups = "drop") %>%
  arrange(desc(n_peaks))

# Take top 20
gene_list <- head(gene_peak_counts$gene_id, 100)

# Show gene info
gene_info <- peak_to_gene %>%
  dplyr::filter(gene_id %in% gene_list) %>%
  group_by(gene_id, gene_name) %>%
  summarise(n_peaks = n(), .groups = "drop") %>%
  arrange(desc(n_peaks))

print(gene_info)

# -----------------------------
# RUN MEDIATION ANALYSIS
# -----------------------------

cat("\n=== Running Mediation Analysis ===\n")

# Run with progress bar
med_results <- pblapply(gene_list, function(g) {
  run_mediation_for_gene(g, seu, peak_to_gene, verbose = FALSE)
})

names(med_results) <- gene_list

# Remove failed analyses
med_results <- Filter(Negate(is.null), med_results)

cat(sprintf("Completed: %d/%d genes successfully analyzed\n", 
            length(med_results), length(gene_list)))

# -----------------------------
# EXTRACT RESULTS
# -----------------------------

extract_mediation_stats <- function(result) {
  if (is.null(result)) return(NULL)
  
  effects <- result$summary$effect.pe
  
  # returns NA if the named effect is missing
  get_eff <- function(nm) {
    if (!is.null(effects) && nm %in% names(effects)) {
      return(effects[[nm]])
    } else {
      return(NA_real_)
    }
  }
  
  data.frame(
    gene_id        = result$gene_id,
    gene_name      = result$gene_name,
    n_peaks        = result$n_peaks,
    n_cells        = result$n_cells,
    # te   = total effect,
    # pnde = pure natural direct effect,
    # pnie = pure natural indirect effect,
    # pm   = overall proportion mediated,
    total_effect   = get_eff("te"),
    direct_effect  = get_eff("pnde"),
    indirect_effect = get_eff("pnie"),
    prop_mediated  = get_eff("pm"),
    stringsAsFactors = FALSE
  )
}

# Apply extraction to each result
results_list <- lapply(med_results, extract_mediation_stats)

# Drop any NULLs that may still be present
results_list <- Filter(Negate(is.null), results_list)
names(results_list) <- NULL

results_table <- do.call(rbind, results_list)


head(results_table)

# -----------------------------
# SAVE RESULTS
# -----------------------------

saveRDS(med_results, "mediation_results_20genes.rds")
write.csv(results_table, "mediation_results_20genes.csv", row.names = FALSE)

