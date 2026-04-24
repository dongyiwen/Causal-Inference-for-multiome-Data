suppressPackageStartupMessages({
  library(dplyr)
  library(pbapply)
  library(CMAverse)
  library(mediation)
  library(ggplot2)
  library(gridExtra)
})

# =============================================================================
# PART 1: RUN CMAVERSE MEDIATION ANALYSIS
# =============================================================================


# Function: CMAverse mediation for one gene
run_mediation_cmaverse <- function(gene_id, seu, peak_to_gene, 
                                   exposure_var = "Sample.Age", 
                                   exposed_level = "pcw20",
                                   verbose = TRUE) {
  
  if (verbose) cat(sprintf("  Analyzing %s...\n", gene_id))
  
  if (!gene_id %in% rownames(seu[["RNA"]])) {
    if (verbose) warning(sprintf("Gene %s not found", gene_id))
    return(NULL)
  }
  
  expr_vec <- tryCatch({
    as.numeric(LayerData(seu, assay = "RNA", layer = "data")[gene_id, ])
  }, error = function(e) {
    as.numeric(seu[["RNA"]]@data[gene_id, ])
  })
  
  gene_peaks <- peak_to_gene$peak[peak_to_gene$gene_id == gene_id]
  if (length(gene_peaks) == 0) return(NULL)
  
  gene_peaks <- intersect(gene_peaks, rownames(seu[["ATAC"]]))
  if (length(gene_peaks) == 0) return(NULL)
  
  atac_counts <- tryCatch({
    LayerData(seu, assay = "ATAC", layer = "counts")[gene_peaks, , drop = FALSE]
  }, error = function(e) {
    seu[["ATAC"]]@counts[gene_peaks, , drop = FALSE]
  })
  
  chrom_vec <- colMeans(as.matrix(atac_counts))
  stage_vec <- ifelse(seu[[exposure_var]][, 1] %in% exposed_level, 1, 0)
  
  med_data <- data.frame(
    stage = stage_vec,
    chrom_access = chrom_vec,
    expression = expr_vec
  )
  
  med_data <- na.omit(med_data)
  
  if (nrow(med_data) < 30 || sd(med_data$chrom_access) < 1e-6 || 
      sd(med_data$expression) < 1e-6) {
    return(NULL)
  }
  
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
    
    list(
      gene_id = gene_id,
      gene_name = peak_to_gene$gene_name[peak_to_gene$gene_id == gene_id][1],
      n_peaks = length(gene_peaks),
      n_cells = nrow(med_data),
      fit = fit,
      summary = summary(fit)
    )
  }, error = function(e) NULL)
  
  return(result)
}

# Prepare peak data
if (any(grepl(":", peak_to_gene$peak))) {
  peak_to_gene$peak <- gsub(":", "-", peak_to_gene$peak)
}

# Select 100 genes
genes_with_peaks <- unique(peak_to_gene$gene_id)
valid_genes <- intersect(genes_with_peaks, rownames(seu[["RNA"]]))

gene_peak_counts <- peak_to_gene %>%
  dplyr::filter(gene_id %in% valid_genes) %>%
  group_by(gene_id) %>%
  summarise(n_peaks = n(), .groups = "drop") %>%
  arrange(desc(n_peaks))

gene_list <- head(gene_peak_counts$gene_id, 100)

# Run CMAverse
cmaverse_results <- pblapply(gene_list, function(g) {
  run_mediation_cmaverse(g, seu, peak_to_gene, verbose = FALSE)
})
names(cmaverse_results) <- gene_list
cmaverse_results <- Filter(Negate(is.null), cmaverse_results)

# Extract results
extract_cmaverse <- function(result) {
  if (is.null(result)) return(NULL)
  effects <- result$summary$effect.pe
  get_eff <- function(nm) {
    if (!is.null(effects) && nm %in% names(effects)) effects[[nm]] else NA_real_
  }
  data.frame(
    gene_id = result$gene_id,
    gene_name = result$gene_name,
    n_peaks = result$n_peaks,
    n_cells = result$n_cells,
    total_effect = get_eff("te"),
    direct_effect = get_eff("pnde"),
    indirect_effect = get_eff("pnie"),
    prop_mediated = get_eff("pm"),
    stringsAsFactors = FALSE
  )
}

cmaverse_list <- lapply(cmaverse_results, extract_cmaverse)
cmaverse_list <- Filter(Negate(is.null), cmaverse_list)
names(cmaverse_list) <- NULL
cmaverse_table <- do.call(rbind, cmaverse_list)

# =============================================================================
# PART 2: RUN MEDIATION PACKAGE ANALYSIS
# =============================================================================


run_mediation_package <- function(gene_id, seu, peak_to_gene, 
                                  exposure_var = "Sample.Age", 
                                  exposed_level = "pcw20",
                                  n_sims = 100,
                                  verbose = TRUE) {
  
  gene_name <- gene_id
  if (gene_id %in% peak_to_gene$gene_id) {
    tmp <- peak_to_gene$gene_name[peak_to_gene$gene_id == gene_id][1]
    if (!is.na(tmp) && tmp != "") gene_name <- tmp
  }
  
  if (!gene_id %in% rownames(seu[["RNA"]])) return(NULL)
  
  expr_vec <- tryCatch({
    as.numeric(LayerData(seu, assay = "RNA", layer = "data")[gene_id, ])
  }, error = function(e) {
    as.numeric(seu[["RNA"]]@data[gene_id, ])
  })
  
  gene_peaks <- peak_to_gene$peak[peak_to_gene$gene_id == gene_id]
  if (length(gene_peaks) == 0) return(NULL)
  
  atac_mat <- tryCatch({
    LayerData(seu, assay = "ATAC", layer = "counts")[gene_peaks, , drop = FALSE]
  }, error = function(e) {
    seu[["ATAC"]]@counts[gene_peaks, , drop = FALSE]
  })
  chrom_vec <- colMeans(as.matrix(atac_mat))
  
  stage_vec <- ifelse(seu[[exposure_var]][, 1] == exposed_level, 1, 0)
  
  med_data <- data.frame(
    treat = stage_vec,
    mediator = chrom_vec,
    outcome = expr_vec,
    stringsAsFactors = FALSE
  )
  
  med_data <- na.omit(med_data)
  
  if (sd(med_data$mediator) < 1e-6 || length(unique(med_data$treat)) < 2) {
    return(NULL)
  }
  
  tryCatch({
    med_model <- lm(mediator ~ treat, data = med_data)
    out_model <- lm(outcome ~ mediator + treat, data = med_data)
    
    med_out <- mediate(
      model.m = med_model,
      model.y = out_model,
      treat = "treat",
      mediator = "mediator",
      boot = TRUE,
      sims = n_sims
    )
    
    med_summary <- summary(med_out)
    
    list(
      gene_id = gene_id,
      gene_name = gene_name,
      n_peaks = length(gene_peaks),
      n_cells = nrow(med_data),
      
      # NIE (treated): Y(1, M(1)) - Y(1, M(0))
      NIE_treated = med_out$d1,
      NIE_treated_pval = med_out$d1.p,
      NIE_treated_ci_lower = med_out$d1.ci[1],
      NIE_treated_ci_upper = med_out$d1.ci[2],
      
      # NDE (control): Y(1, M(0)) - Y(0, M(0))
      NDE_control = med_out$z0,
      NDE_control_pval = med_out$z0.p,
      NDE_control_ci_lower = med_out$z0.ci[1],
      NDE_control_ci_upper = med_out$z0.ci[2],
      
      # Average effects (for comparison)
      ACME_average = med_out$d.avg,
      ADE_average = med_out$z.avg,
      total_effect = med_out$tau.coef,
      total_pval = med_out$tau.p,
      prop_mediated = med_out$n.avg,
      prop_mediated_pval = med_out$n.avg.p,
      
      # Full objects
      mediation_object = med_out,
      mediator_model = med_model,
      outcome_model = out_model,
      summary = med_summary
    )
  }, error = function(e) NULL)
}

cat("\n=== Running Mediation Package ===\n")
mediation_results <- pblapply(gene_list, function(g) {
  run_mediation_package(g, seu, peak_to_gene, n_sims = 100, verbose = FALSE)
})
names(mediation_results) <- gene_list
mediation_results <- Filter(Negate(is.null), mediation_results)

cat(sprintf("✓ Mediation: %d genes analyzed\n", length(mediation_results)))

extract_mediation <- function(result) {
  if (is.null(result)) return(NULL)
  data.frame(
    gene_id = result$gene_id,
    gene_name = result$gene_name,
    n_peaks = result$n_peaks,
    n_cells = result$n_cells,
    
    # NIE (treated): Y(1, M(1)) - Y(1, M(0))
    NIE_treated = result$NIE_treated,
    NIE_treated_pval = result$NIE_treated_pval,
    NIE_treated_ci_lower = result$NIE_treated_ci_lower,
    NIE_treated_ci_upper = result$NIE_treated_ci_upper,
    
    # NDE (control): Y(1, M(0)) - Y(0, M(0))
    NDE_control = result$NDE_control,
    NDE_control_pval = result$NDE_control_pval,
    NDE_control_ci_lower = result$NDE_control_ci_lower,
    NDE_control_ci_upper = result$NDE_control_ci_upper,
    
    # Average effects
    ACME_avg = result$ACME_average,
    ADE_avg = result$ADE_average,
    total_effect = result$total_effect,
    total_pval = result$total_pval,
    prop_mediated = result$prop_mediated,
    prop_mediated_pval = result$prop_mediated_pval,
    
    # Significance flags
    significant_NIE = result$NIE_treated_pval < 0.05,
    significant_NDE = result$NDE_control_pval < 0.05,
    significant_mediation = result$NIE_treated_pval < 0.05,
    significant_direct = result$NDE_control_pval < 0.05,
    
    stringsAsFactors = FALSE
  )
}

mediation_table <- do.call(rbind, lapply(mediation_results, extract_mediation))

# =============================================================================
# PART 3: INDIVIDUAL METHOD SUMMARIES
# =============================================================================


cat("\n=== CMAverse Summary ===\n")
cat(sprintf("Genes: %d | Cells: %.0f | Peaks: %.1f\n", 
            nrow(cmaverse_table), mean(cmaverse_table$n_cells), 
            mean(cmaverse_table$n_peaks)))
cat(sprintf("Indirect: %.4f ± %.4f\n", mean(cmaverse_table$indirect_effect), 
            sd(cmaverse_table$indirect_effect)))
cat(sprintf("Direct:   %.4f ± %.4f\n", mean(cmaverse_table$direct_effect), 
            sd(cmaverse_table$direct_effect)))
cat(sprintf("Total:    %.4f ± %.4f\n", mean(cmaverse_table$total_effect), 
            sd(cmaverse_table$total_effect)))
cat(sprintf("Prop Med: %.4f ± %.4f\n", mean(cmaverse_table$prop_mediated, na.rm=T), 
            sd(cmaverse_table$prop_mediated, na.rm=T)))

cat("\n=== Mediation Package Summary ===\n")
cat(sprintf("Genes: %d | Cells: %.0f | Peaks: %.1f\n", 
            nrow(mediation_table), mean(mediation_table$n_cells), 
            mean(mediation_table$n_peaks)))
cat("\nNIE (treated): Y(1,M(1)) - Y(1,M(0))\n")
cat(sprintf("  Estimate: %.4f ± %.4f\n", mean(mediation_table$NIE_treated), 
            sd(mediation_table$NIE_treated)))
cat(sprintf("  Significant: %d (%.1f%%)\n", sum(mediation_table$significant_NIE),
            100*mean(mediation_table$significant_NIE)))
cat("\nNDE (control): Y(1,M(0)) - Y(0,M(0))\n")
cat(sprintf("  Estimate: %.4f ± %.4f\n", mean(mediation_table$NDE_control), 
            sd(mediation_table$NDE_control)))
cat(sprintf("  Significant: %d (%.1f%%)\n", sum(mediation_table$significant_NDE),
            100*mean(mediation_table$significant_NDE)))
cat("\nAverage Effects:\n")
cat(sprintf("  Indirect (ACME avg): %.4f ± %.4f\n", mean(mediation_table$ACME_avg), 
            sd(mediation_table$ACME_avg)))
cat(sprintf("  Direct (ADE avg):    %.4f ± %.4f\n", mean(mediation_table$ADE_avg), 
            sd(mediation_table$ADE_avg)))
cat(sprintf("  Total:               %.4f ± %.4f\n", mean(mediation_table$total_effect), 
            sd(mediation_table$total_effect)))
cat(sprintf("  Prop Med:            %.4f ± %.4f\n", mean(mediation_table$prop_mediated, na.rm=T), 
            sd(mediation_table$prop_mediated, na.rm=T)))
cat(sprintf("Significant: Indirect=%d (%.1f%%), Direct=%d (%.1f%%)\n",
            sum(mediation_table$significant_mediation),
            100*mean(mediation_table$significant_mediation),
            sum(mediation_table$significant_direct),
            100*mean(mediation_table$significant_direct)))

# =============================================================================
# PART 4: METHOD COMPARISON
# =============================================================================

comparison_df <- inner_join(
  cmaverse_table %>% rename_with(~paste0("CMA_", .), -c(gene_id, gene_name)),
  mediation_table %>% rename_with(~paste0("MED_", .), -c(gene_id, gene_name)),
  by = "gene_id", suffix = c("_CMA", "_MED")
)

comparison_df <- comparison_df %>%
  mutate(
    diff_indirect = MED_NIE_treated - CMA_indirect_effect,
    diff_direct = MED_NDE_control - CMA_direct_effect,
    diff_total = MED_total_effect - CMA_total_effect,
    indirect_match = sign(MED_NIE_treated) == sign(CMA_indirect_effect),
    direct_match = sign(MED_NDE_control) == sign(CMA_direct_effect)
  )

cat(sprintf("\nMerged: %d genes\n", nrow(comparison_df)))
cat("\nCorrelations:\n")
cat(sprintf("  Indirect: r=%.3f\n", cor(comparison_df$CMA_indirect_effect, 
                                        comparison_df$MED_NIE_treated)))
cat(sprintf("  Direct:   r=%.3f\n", cor(comparison_df$CMA_direct_effect, 
                                        comparison_df$MED_NDE_control)))
cat(sprintf("  Total:    r=%.3f\n", cor(comparison_df$CMA_total_effect, 
                                        comparison_df$MED_total_effect)))

cat("\nDirection Agreement:\n")
cat(sprintf("  Indirect: %.1f%%\n", 100*mean(comparison_df$indirect_match)))
cat(sprintf("  Direct:   %.1f%%\n", 100*mean(comparison_df$direct_match)))

# =============================================================================
# PART 5: FOREST PLOTS
# =============================================================================


# Mediation package forest plot
plot_med_forest <- function(data_table, effect = "indirect") {
  if (effect == "indirect") {
    df <- data_table %>% arrange(desc(abs(NIE_treated))) %>%
      mutate(label = ifelse(is.na(gene_name)|gene_name=="", gene_id, gene_name),
             label = factor(label, levels=label),
             sig = NIE_treated_pval < 0.05,
             val = NIE_treated, 
             lo = NIE_treated_ci_lower, 
             hi = NIE_treated_ci_upper)
    title <- "Indirect Effect - NIE (treated): Y(1,M(1)) - Y(1,M(0))"
    xlab <- "NIE (treated)"
  } else if (effect == "direct") {
    df <- data_table %>% arrange(desc(abs(NDE_control))) %>%
      mutate(label = ifelse(is.na(gene_name)|gene_name=="", gene_id, gene_name),
             label = factor(label, levels=label),
             sig = NDE_control_pval < 0.05,
             val = NDE_control, 
             lo = NDE_control_ci_lower, 
             hi = NDE_control_ci_upper)
    title <- "Direct Effect - NDE (control): Y(1,M(0)) - Y(0,M(0))"
    xlab <- "NDE (control)"
  } else {
    df <- data_table %>% arrange(desc(abs(total_effect))) %>%
      mutate(label = ifelse(is.na(gene_name)|gene_name=="", gene_id, gene_name),
             label = factor(label, levels=label),
             sig = total_pval < 0.05,
             val = total_effect, 
             lo = val-1.96*abs(val)/3, 
             hi = val+1.96*abs(val)/3)
    title <- "Total Effect: NIE + NDE"
    xlab <- "Total Effect"
  }
  
  ggplot(df, aes(x=val, y=label)) +
    geom_vline(xintercept=0, linetype="dashed", color="gray50") +
    geom_errorbarh(aes(xmin=lo, xmax=hi, color=sig), height=0.3, linewidth=0.8) +
    geom_point(aes(color=sig, size=abs(val))) +
    scale_color_manual(values=c("FALSE"="gray60","TRUE"="#E41A1C"),
                       labels=c("FALSE"="NS","TRUE"="p<0.05")) +
    scale_size_continuous(range=c(2,5)) +
    labs(title=paste("Mediation Pkg:", title), x=xlab, y="Gene",
         color="Significance", size="|Effect|") +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold"), legend.position="right",
          panel.grid.minor=element_blank(), axis.text.y=element_text(size=7))
}

# CMAverse forest plot
plot_cma_forest <- function(data_table, effect = "indirect") {
  if (effect == "indirect") {
    df <- data_table %>% arrange(desc(abs(indirect_effect))) %>%
      mutate(label = ifelse(is.na(gene_name)|gene_name=="", gene_id, gene_name),
             label = factor(label, levels=label),
             val = indirect_effect)
    title <- "Indirect Effect (NIE)"
    xlab <- "Indirect Effect"
  } else if (effect == "direct") {
    df <- data_table %>% arrange(desc(abs(direct_effect))) %>%
      mutate(label = ifelse(is.na(gene_name)|gene_name=="", gene_id, gene_name),
             label = factor(label, levels=label),
             val = direct_effect)
    title <- "Direct Effect (NDE)"
    xlab <- "Direct Effect"
  } else {
    df <- data_table %>% arrange(desc(abs(total_effect))) %>%
      mutate(label = ifelse(is.na(gene_name)|gene_name=="", gene_id, gene_name),
             label = factor(label, levels=label),
             val = total_effect)
    title <- "Total Effect"
    xlab <- "Total Effect"
  }
  
  df <- df %>% mutate(dir = ifelse(val > 0, "Positive", "Negative"))
  
  ggplot(df, aes(x=val, y=label)) +
    geom_vline(xintercept=0, linetype="dashed", color="gray50") +
    geom_point(aes(color=dir, size=abs(val))) +
    scale_color_manual(values=c("Positive"="#E41A1C","Negative"="#377EB8")) +
    scale_size_continuous(range=c(2,5)) +
    labs(title=paste("CMAverse:", title), x=xlab, y="Gene",
         color="Direction", size="|Effect|") +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold"), legend.position="right",
          panel.grid.minor=element_blank(), axis.text.y=element_text(size=7))
}

# Create plots
med_ind <- plot_med_forest(mediation_table, "indirect")
med_dir <- plot_med_forest(mediation_table, "direct")
med_tot <- plot_med_forest(mediation_table, "total")

cma_ind <- plot_cma_forest(cmaverse_table, "indirect")
cma_dir <- plot_cma_forest(cmaverse_table, "direct")
cma_tot <- plot_cma_forest(cmaverse_table, "total")

# Save plots
ggsave("mediation_indirect.pdf", med_ind, width=10, height=12)
ggsave("mediation_direct.pdf", med_dir, width=10, height=12)
ggsave("mediation_total.pdf", med_tot, width=10, height=12)
ggsave("cmaverse_indirect.pdf", cma_ind, width=10, height=12)
ggsave("cmaverse_direct.pdf", cma_dir, width=10, height=12)
ggsave("cmaverse_total.pdf", cma_tot, width=10, height=12)

comparison_plot <- grid.arrange(cma_ind, med_ind, ncol=2)
ggsave("comparison_indirect.pdf", comparison_plot, width=18, height=12)

# Save tables
write.csv(cmaverse_table, "cmaverse_results.csv", row.names=FALSE)
write.csv(mediation_table, "mediation_results.csv", row.names=FALSE)
write.csv(comparison_df, "method_comparison.csv", row.names=FALSE)
