
# ---------------------------
# 0. Libraries
# ---------------------------
rm(list=ls())
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2020)   
  library(TFBSTools)    
  library(dplyr)
  library(pbapply)
  library(CMAverse)  
})





# -----------------------------
# 1. Load metadata
# -----------------------------
cat("Loading metadata...\n")
meta <- fread("/Users/yiwendong/Downloads/GSE162170_multiome_cell_metadata.txt.gz")
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell
true_cell_barcodes <- meta$Cell.Barcode

cat(sprintf("  Loaded %d cells\n", nrow(meta)))

summary(factor(meta$Sample.Age))

set.seed(1) 
n <- nrow(meta)
idx <- sample(seq_len(n), size = floor(n/2))
meta$Sample.Age <- as.character(meta$Sample.Age)
meta$Sample.Age[idx] <- "pcw20"
meta$Sample.Age <- factor(meta$Sample.Age)





# -----------------------------
# 2. Load RNA counts
# -----------------------------
cat("Loading RNA counts...\n")
rna_dt <- fread("/Users/yiwendong/Downloads/GSE162170_multiome_rna_counts.tsv.gz")
rna_mat <- as.matrix(rna_dt[, -1])
rownames(rna_mat) <- rna_dt[[1]]
colnames(rna_dt)[1] <- "gene"

# Convert to sparse matrix (saves memory)
rna_mat <- as(rna_mat, "dgCMatrix")

cat(sprintf("  RNA matrix: %d genes x %d cells\n", nrow(rna_mat), ncol(rna_mat)))




# -----------------------------
# 3. Select 3000 most variable genes
# -----------------------------

cat("Selecting top 3000 variable genes...\n")
gene_vars <- apply(rna_mat, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:3000]
rna_mat <- rna_mat[top_genes, ]

cat(sprintf("  Retained %d genes\n", length(top_genes)))


# -----------------------------
# 4. Load peak coordinates
# -----------------------------
# Load consensus peaks
cat("Loading ATAC peak coordinates...\n")
peaks <- fread("/Users/yiwendong/Downloads/GSE162170_multiome_atac_consensus_peaks.txt.gz",
               skip = 1, 
               header = FALSE)

# Assign peak names
peaks$name <- peaks$V6

# Create GRanges object for genomic operations
gr.peaks <- makeGRangesFromDataFrame(
  peaks,
  keep.extra.columns = TRUE,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3",
  ignore.strand = TRUE
)
names(gr.peaks) <- peaks$name

cat(sprintf("  Loaded %d ATAC peaks\n", length(gr.peaks)))



# -----------------------------
# 5. Load ATAC peak matrix
# -----------------------------
cat("Loading ATAC counts (matched cells only)...\n")

# Standardize barcode format: remove prefix before last underscore
rna_bc_clean <- gsub(".*_", "", colnames(rna_mat))
atac_bc_clean <- true_cell_barcodes

# Find cells present in BOTH datasets
subset_cells <- intersect(rna_bc_clean, atac_bc_clean)
cat(sprintf("  Found %d cells in both RNA and ATAC\n", length(subset_cells)))


cell_order_atac <- true_cell_barcodes
col_indices <- match(subset_cells, cell_order_atac)
cols_to_read <- sort(unique(c(1, col_indices)))  ## <-- CHANGED LINE

# Read only relevant columns from ATAC counts
atac_dt <- fread(
  "/Users/yiwendong/Downloads/GSE162170_multiome_atac_counts.tsv.gz",
  header = FALSE,
  select = cols_to_read
)

colnames(atac_dt)[1] <- "peak"



## ------------------- BUILD ATAC MATRIX ------------------- ##
cat("Converting ATAC data.table to sparse matrix...\n")

n_rows <- nrow(atac_dt)          # number of peaks
n_cols <- ncol(atac_dt) - 1L     # number of cells (selected)

# Lists to collect indices and values for each column
i_list <- vector("list", n_cols)
x_list <- vector("list", n_cols)

for (j in seq_len(n_cols)) {
  col_vec <- as.numeric(atac_dt[[j + 1L]])
  nz_idx <- which(col_vec != 0)
  
  if (length(nz_idx) > 0L) {
    i_list[[j]] <- nz_idx
    x_list[[j]] <- col_vec[nz_idx]
  } else {
    i_list[[j]] <- integer(0)
    x_list[[j]] <- numeric(0)
  }
}

# Concatenate all column-wise indices/values
i <- unlist(i_list, use.names = FALSE)
x <- unlist(x_list, use.names = FALSE)

# Column pointers: cumulative counts of non-zeros per column
p <- c(0L, cumsum(vapply(i_list, length, integer(1L))))

# Safety: ensure integer types and consistent sizes
i <- as.integer(i - 1L)                   
p <- as.integer(p)
Dim <- as.integer(c(n_rows, n_cols))      

# Optional sanity checks (you can leave these in or remove once it works)
stopifnot(length(p) == (n_cols + 1L))
stopifnot(length(x) == p[length(p)])

atac_mat <- new(
  "dgCMatrix",
  i = i,
  p = p,
  x = as.numeric(x),
  Dim = Dim,
  Dimnames = list(
    atac_dt[["peak"]],           
    colnames(atac_dt)[-1]        
  )
)

cat(sprintf("  ATAC matrix: %d peaks x %d cells (sparse)\n",
            nrow(atac_mat), ncol(atac_mat)))


peak_ids <- peaks$V6
if (nrow(atac_mat) != length(peak_ids)) {
  warning(sprintf("Peak count mismatch: %d peaks vs %d rows. Using first %d peaks.",
                  length(peak_ids), nrow(atac_mat), nrow(atac_mat)))
  peak_ids <- peak_ids[1:nrow(atac_mat)]
  gr.peaks <- gr.peaks[1:nrow(atac_mat)]
}

if (ncol(atac_mat) != length(subset_cells)) {
  warning(sprintf(
    "Cell count mismatch: %d subset_cells vs %d ATAC columns. Using first %d cells.",
    length(subset_cells), ncol(atac_mat), ncol(atac_mat)
  ))
  
  subset_cells <- subset_cells[1:ncol(atac_mat)]
  
  if (ncol(rna_mat) != ncol(atac_mat)) {
    rna_mat <- rna_mat[, 1:ncol(atac_mat)]
  }
}

colnames(atac_mat) <- subset_cells
colnames(rna_mat) <- subset_cells

cat(sprintf("  ATAC matrix: %d peaks x %d cells\n", nrow(atac_mat), ncol(atac_mat)))

cat("Aligning RNA, ATAC, and metadata...\n")

# Standardize RNA column names
colnames(rna_mat) <- gsub(".*_", "", colnames(rna_mat))

# Find shared cells
shared_cells <- Reduce(intersect, list(
  colnames(rna_mat),
  colnames(atac_mat),
  meta$Cell.Barcode
))

cat(sprintf("  Final cell count: %d\n", length(shared_cells)))

# Subset all to shared cells in same order
rna_mat <- rna_mat[, shared_cells]
atac_mat <- atac_mat[, shared_cells]
meta <- meta[match(shared_cells, meta$Cell.Barcode), ]

# Verify alignment
stopifnot(
  all(colnames(rna_mat) == colnames(atac_mat)),
  all(colnames(rna_mat) == meta$Cell.Barcode),
  all(!is.na(meta$Cell.Barcode))
)

cat("All datasets aligned\n")


cat("Building Seurat object...\n")

# Create ChromatinAssay for ATAC data
chrom_assay <- CreateChromatinAssay(
  counts = atac_mat,
  ranges = gr.peaks,
  genome = "hg38"
)

# Create Seurat object with RNA data
seu <- CreateSeuratObject(
  counts = rna_mat,
  assay = "RNA",
  meta.data = meta
)

# Add ATAC as second assay
seu[["ATAC"]] <- chrom_assay
DefaultAssay(seu) <- "RNA"

cat(" Multiome object created\n")


# -----------------------------
# 8. Basic normalization
# -----------------------------

# RNA

cat("Normalizing data...\n")

# RNA normalization
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, assay = "RNA")
seu <- FindVariableFeatures(seu, assay = "RNA", nfeatures = 3000)
seu <- ScaleData(seu, assay = "RNA")

# ATAC normalization (TF-IDF: term frequency-inverse document frequency)
DefaultAssay(seu) <- "ATAC"
seu <- RunTFIDF(seu)
seu <- FindTopFeatures(seu, min.cutoff = "q75")
seu <- RunSVD(seu)

DefaultAssay(seu) <- "RNA"


cat("  ✓ Saved: multiome_object_normalized.rds\n")



# -----------------------------
# 9. Link peaks to genes (promoter peaks)
# -----------------------------
# Use the Seurat multiome object

cat("Linking peaks to genes...\n")


library(GenomeInfoDb)
library(ensembldb)

DefaultAssay(seu) <- "ATAC"

## 0. Convert peak names (rownames of ATAC assay) to GRanges
peak_ranges <- StringToGRanges(
  rownames(seu[["ATAC"]]),
  sep = c(":", "-")
)

## 1. Build gene annotation GRanges
library(EnsDb.Hsapiens.v86)
gene_ann <- ensembldb::genes(EnsDb.Hsapiens.v86)

# Make seqlevels style consistent between peaks and genes
seqlevelsStyle(gene_ann) <- seqlevelsStyle(peak_ranges)

## 2. Find nearest gene for each peak
hits <- distanceToNearest(
  peak_ranges,
  gene_ann,
  ignore.strand = TRUE
)

peak_idx <- queryHits(hits)
gene_idx <- subjectHits(hits)

## 3. Build peak strings explicitly (chr:start-end)
peak_strings <- paste0(
  as.character(seqnames(peak_ranges)[peak_idx]), ":",
  start(peak_ranges)[peak_idx], "-",
  end(peak_ranges)[peak_idx]
)

## 4. Extract gene metadata
gene_meta <- mcols(gene_ann)

# Raw Ensembl IDs
gene_id_vec_raw <- gene_meta$gene_id
gene_id_vec <- sub("\\.\\d+$", "", gene_id_vec_raw)

gene_name_vec <- if ("gene_name" %in% colnames(gene_meta)) {
  gene_meta$gene_name
} else if ("symbol" %in% colnames(gene_meta)) {
  gene_meta$symbol
} else {
  rep(NA_character_, length(gene_ann))
}

## 5. Assemble peak to gene table
peak_to_gene <- data.frame(
  peak      = peak_strings,
  gene_id   = gene_id_vec[gene_idx],    # cleaned Ensembl IDs
  gene_name = gene_name_vec[gene_idx],  # symbols if available
  distance  = mcols(hits)$distance,
  stringsAsFactors = FALSE
)

rownames(peak_to_gene) <- peak_to_gene$peak

## 6. Filter to promoter peaks (within ±2 kb of TSS)
peak_to_gene <- peak_to_gene %>%
  dplyr::filter(abs(distance) <= 2000)

cat(sprintf("  Linked %d promoter peaks to genes\n", nrow(peak_to_gene)))

## 7. Coverage: how many RNA genes (Ensembl IDs) have promoter peaks

# RNA feature names are Ensembl IDs
rna_genes <- rownames(seu[["RNA"]])

# Genes that have at least one promoter peak (by Ensembl ID)
genes_with_peaks <- unique(peak_to_gene$gene_id)

# Intersect with RNA features
covered_rna_genes <- intersect(genes_with_peaks, rna_genes)

cat(sprintf("  %d RNA genes (of %d) have promoter peaks\n",
            length(covered_rna_genes), length(rna_genes)))

cat("Finished linking peaks to genes.\n")
