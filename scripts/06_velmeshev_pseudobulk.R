# 06_velmeshev_pseudobulk.R
# Build "pseudobulk" from the Velmeshev 2019 single-nucleus data:
# add up the raw UMI counts of all cells from the same sample and cell type.
# Why: single cells are noisy and are not independent (many cells per person), so
# the statistical unit must be the PERSON/sample, not the cell.
#
# The count file is 3.3 GB (238 million non-zero entries) and my Mac has ~8.6 GB RAM,
# so the file is read in chunks and only the per-group sums are kept.
#
# Outputs (data_clean/):
#   velmeshev_pb_counts.rds   sparse genes x (sample|celltype) matrix of summed counts
#   velmeshev_pb_meta.csv     one row per column of that matrix
#   velmeshev_pb_total_*.csv  all cells of a sample summed (acts like "bulk" tissue)
#
# Run:  Rscript scripts/06_velmeshev_pseudobulk.R

suppressMessages({
  library(data.table)
  library(Matrix)
})

in_dir    <- "data_raw/velmeshev2019/extracted"
clean_dir <- "data_clean"

# ---- 1. which group (sample|celltype) does each cell belong to? ----
meta <- fread(file.path(in_dir, "meta.txt"))
barcodes <- readLines(file.path(in_dir, "barcodes.tsv"))
stopifnot(length(barcodes) == nrow(meta))
stopifnot(all(barcodes == meta$cell))   # same cell order as the matrix columns
meta$group <- paste(meta$sample, meta$cluster, sep = "|")
group_levels <- sort(unique(meta$group))
cell_group <- match(meta$group, group_levels)
cat("cells:", nrow(meta), "| sample x celltype groups:", length(group_levels), "\n")

genes <- fread(file.path(in_dir, "genes.tsv"), header = FALSE)
n_genes <- nrow(genes)

# ---- 2. read the matrix in chunks and sum counts per (gene, group) ----
mtx <- file.path(in_dir, "matrix.mtx")
header <- scan(mtx, what = "", nlines = 3, sep = "\n", quiet = TRUE)
dims <- as.numeric(strsplit(header[3], " ")[[1]])
n_entries <- dims[3]
cat("matrix:", dims[1], "genes x", dims[2], "cells,", n_entries, "entries\n")

chunk_size <- 25e6
n_chunks <- ceiling(n_entries / chunk_size)
sums_list <- list()
for (k in 1:n_chunks) {
  skip_n <- 3 + (k - 1) * chunk_size
  take_n <- min(chunk_size, n_entries - (k - 1) * chunk_size)
  d <- fread(mtx, skip = skip_n, nrows = take_n, header = FALSE,
             col.names = c("gene", "cell", "count"), colClasses = "integer")
  d[, grp := cell_group[cell]]
  agg <- d[, .(count = sum(count)), by = .(gene, grp)]   # sum over cells in a group
  sums_list[[k]] <- agg
  cat(" chunk", k, "of", n_chunks, "done;", nrow(agg), "gene-group sums\n")
  rm(d, agg); gc(verbose = FALSE)
}
all_sums <- rbindlist(sums_list)
all_sums <- all_sums[, .(count = sum(count)), by = .(gene, grp)]   # groups split across chunks
pb <- sparseMatrix(i = all_sums$gene, j = all_sums$grp, x = all_sums$count,
                   dims = c(n_genes, length(group_levels)))
rownames(pb) <- genes$V2          # gene symbols
colnames(pb) <- group_levels

# duplicated symbols (a few genes share one) -> add them together
if (any(duplicated(rownames(pb)))) {
  cat("duplicated gene symbols:", sum(duplicated(rownames(pb))), "-> summed\n")
  pb <- rowsum(as.matrix(pb), group = rownames(pb))   # dense after this, fine (see below)
  pb <- Matrix(pb, sparse = TRUE)
}
cat("pseudobulk matrix:", nrow(pb), "genes x", ncol(pb), "groups\n")
saveRDS(pb, file.path(clean_dir, "velmeshev_pb_counts.rds"))

# ---- 3. metadata for each pseudobulk column ----
first <- meta[!duplicated(group), ]
first <- first[match(group_levels, first$group), ]
n_cells <- as.numeric(table(meta$group)[group_levels])
pbm <- data.frame(group_id = group_levels, sample = first$sample, individual = first$individual,
                  celltype = first$cluster, region = first$region, age = first$age,
                  sex = first$sex, diagnosis = first$diagnosis, seqbatch = first$Seqbatch,
                  capbatch = first$Capbatch, pmi = first$`post-mortem interval (hours)`,
                  rin = first$`RNA Integrity Number`, n_cells = n_cells,
                  total_umis = colSums(pb), stringsAsFactors = FALSE)
write.csv(pbm, file.path(clean_dir, "velmeshev_pb_meta.csv"), row.names = FALSE)

# ---- 4. total pseudobulk per sample (all cell types added) ----
samples <- sort(unique(pbm$sample))
tot <- matrix(0, nrow = nrow(pb), ncol = length(samples), dimnames = list(rownames(pb), samples))
for (s in samples) {
  cols <- which(pbm$sample == s)
  tot[, s] <- as.numeric(rowSums(pb[, cols, drop = FALSE]))
}
sm <- pbm[!duplicated(pbm$sample), c("sample", "individual", "region", "age", "sex", "diagnosis", "seqbatch", "capbatch", "pmi", "rin")]
sm <- sm[match(samples, sm$sample), ]
sm$n_cells <- as.numeric(table(meta$sample)[samples])
write.csv(tot, file.path(clean_dir, "velmeshev_pb_total_counts.csv"))
write.csv(sm, file.path(clean_dir, "velmeshev_pb_total_meta.csv"), row.names = FALSE)

cat("\n=========== PLAIN-ENGLISH SUMMARY ===========\n")
cat("Samples:", length(samples), "(", sum(sm$diagnosis == "ASD"), "ASD /", sum(sm$diagnosis == "Control"), "control samples)\n")
cat("Cell types:", length(unique(pbm$celltype)), "\n")
cat("CHECK: groups with very few cells (<20) are unreliable and are filtered before testing.\n")
