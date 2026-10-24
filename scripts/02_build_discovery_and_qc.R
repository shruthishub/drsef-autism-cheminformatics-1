# 02_build_discovery_and_qc.R
# Step: build the human-only discovery matrix (autism from GSE59288 + human
# controls from GSE51264) and run basic QC (PCA, sample correlation, batch check).
#
# Why base R + limma only: the data are already normalized-ish expression
# tables (RPKM-like), so we only need to read, filter, log-transform, and plot.
#
# Run:  Rscript scripts/02_build_discovery_and_qc.R

library(GEOquery)
library(limma)   # used here for plotMDS-style checks and later for DE

raw_dir   <- "data_raw"
clean_dir <- "data_clean"
res_dir   <- "results"
fig_dir   <- "figures"

# ------------------------------------------------------------------
# 1. Read the two expression tables (genes in rows, samples in columns)
# ------------------------------------------------------------------
aut <- read.table(gzfile(file.path(raw_dir, "GSE59288/GSE59288_exp_mRNA.txt.gz")),
                  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
ctl <- read.table(gzfile(file.path(raw_dir, "GSE51264/GSE51264_rpkm_human.tsv.gz")),
                  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

cat("GSE59288 table:", nrow(aut), "genes x", ncol(aut), "samples\n")
cat("GSE51264 human table:", nrow(ctl), "genes x", ncol(ctl), "samples\n")

# ------------------------------------------------------------------
# 2. Keep HUMAN samples only: drop chimpanzee (ptr) columns from GSE59288
# ------------------------------------------------------------------
keep_aut <- grep("^aut", colnames(aut))
aut <- aut[, keep_aut]
cat("Autism columns kept:", ncol(aut), "(chimpanzee 'ptr' columns dropped)\n")
# GSE51264 'human' file should only have human_XX columns
stopifnot(all(grepl("^human_", colnames(ctl))))

# ------------------------------------------------------------------
# 3. Keep only genes with a proper Ensembl ID that appear in BOTH tables
#    (GSE51264 also has CUFF.* novel lincRNAs and ENSGR chrY copies; skip those)
# ------------------------------------------------------------------
# GSE51264 IDs carry a version suffix (ENSG00000000003.10): remove it so IDs match
ctl <- ctl[grepl("^ENSG[0-9]+\\.[0-9]+$", rownames(ctl)), ]  # drops CUFF.* and ENSGR*
rownames(ctl) <- sub("\\.[0-9]+$", "", rownames(ctl))
aut <- aut[grepl("^ENSG[0-9]+$", rownames(aut)), ]
ctl <- ctl[!duplicated(rownames(ctl)), ]
shared <- intersect(rownames(aut), rownames(ctl))
cat("Genes shared by both tables:", length(shared), "\n")
aut <- aut[shared, ]
ctl <- ctl[shared, ]
expr_raw <- cbind(aut, ctl)

# ------------------------------------------------------------------
# 4. Sample table: age, sex, group
#    Age and sex come from the GEO sample tables (see script 01 outputs)
# ------------------------------------------------------------------
parse_age_years <- function(txt) {
  # GSE59288 ages are written like "age: 730 days" or "age: 12 years"
  out <- rep(NA_real_, length(txt))
  for (i in 1:length(txt)) {
    num <- as.numeric(gsub("[^0-9.]", "", sub("^age[^:]*:", "", txt[i])))
    if (grepl("day", txt[i])) {
      out[i] <- num / 365.25
    } else {
      out[i] <- num
    }
  }
  return(out)
}

p1 <- read.csv(file.path(res_dir, "GSE59288_part1_sample_table.csv"),
               row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
# sample names in the matrix are aut1..aut34; the GEO title ends with "[autN]"
p1$sample <- sub(".*\\[(.*)\\]", "\\1", p1$title)
p1$age_years <- parse_age_years(p1[["characteristics_ch1.2"]])
p1$sex <- sub("gender: ", "", p1[["characteristics_ch1.1"]])
meta_aut <- data.frame(sample = p1$sample, group = "autism",
                       age_years = p1$age_years, sex = p1$sex,
                       source = "GSE59288", stringsAsFactors = FALSE)

# controls: pull sample table from GEO (human part = GPL11154)
g <- getGEO("GSE51264", GSEMatrix = TRUE, getGPL = FALSE, destdir = raw_dir)
pc <- NULL
for (j in 1:length(g)) {
  if (annotation(g[[j]]) == "GPL11154") pc <- pData(g[[j]])
}
stopifnot(!is.null(pc))
write.csv(pc, file.path(res_dir, "GSE51264_human_sample_table.csv"))
# age is split over two columns: "age (year): X" and "age (day): Y"
age_cols <- grep("^characteristics", colnames(pc), value = TRUE)
ctl_age <- rep(NA_real_, nrow(pc))
ctl_sex <- rep(NA_character_, nrow(pc))
for (i in 1:nrow(pc)) {
  yrs <- 0; days <- 0
  for (cc in age_cols) {
    v <- as.character(pc[i, cc])
    if (grepl("^age \\(year\\)", v)) yrs  <- as.numeric(sub(".*: ", "", v))
    if (grepl("^age \\(day\\)", v))  days <- as.numeric(sub(".*: ", "", v))
    if (grepl("^Sex", v)) ctl_sex[i] <- sub(".*: ", "", v)
  }
  ctl_age[i] <- yrs + days / 365.25
}
meta_ctl <- data.frame(sample = pc$title, group = "control",
                       age_years = ctl_age, sex = ctl_sex,
                       source = "GSE51264", stringsAsFactors = FALSE)

meta <- rbind(meta_aut, meta_ctl)
rownames(meta) <- meta$sample
# make sure matrix columns and meta rows line up
missing <- setdiff(colnames(expr_raw), rownames(meta))
if (length(missing) > 0) cat("!! samples in matrix but not in meta:", missing, "\n")
meta <- meta[colnames(expr_raw), ]

cat("\nGroup counts (human only):\n"); print(table(meta$group))
cat("\nAge (years) by group:\n")
print(tapply(meta$age_years, meta$group, summary))
cat("\nSex by group:\n"); print(table(meta$group, meta$sex))

# ------------------------------------------------------------------
# 5. Filter low genes, log2 transform
# ------------------------------------------------------------------
# keep genes with RPKM > 1 in at least 20% of samples
keep <- rep(FALSE, nrow(expr_raw))
min_samples <- ceiling(0.2 * ncol(expr_raw))
for (i in 1:nrow(expr_raw)) {
  if (sum(expr_raw[i, ] > 1) >= min_samples) keep[i] <- TRUE
}
expr_f <- expr_raw[keep, ]
cat("\nGenes kept after low-expression filter:", nrow(expr_f), "of", nrow(expr_raw), "\n")
logexpr <- log2(as.matrix(expr_f) + 1)

# ------------------------------------------------------------------
# 6. QC: PCA, sample correlation, batch (source) check
# ------------------------------------------------------------------
pca <- prcomp(t(logexpr), scale. = TRUE)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cat("\nPC1 and PC2 explain", var_expl[1], "% and", var_expl[2], "% of variance\n")

colors <- ifelse(meta$group == "autism", "red", "blue")
png(file.path(fig_dir, "02_pca_group.png"), width = 900, height = 700)
plot(pca$x[, 1], pca$x[, 2], col = colors, pch = 19,
     xlab = paste0("PC1 (", var_expl[1], "%)"), ylab = paste0("PC2 (", var_expl[2], "%)"),
     main = "PCA of discovery samples (red = autism, blue = control)")
legend("topright", legend = c("autism", "control"), col = c("red", "blue"), pch = 19)
dev.off()

# PC1 vs age: are samples ordered by development?
png(file.path(fig_dir, "02_pca_age.png"), width = 900, height = 700)
plot(log10(meta$age_years + 0.01), pca$x[, 1], col = colors, pch = 19,
     xlab = "log10(age in years + 0.01)", ylab = "PC1",
     main = "PC1 vs age (red = autism, blue = control)")
dev.off()

# how much does each PC correlate with group, age, sex?
cat("\nCorrelation of top PCs with group / age / sex:\n")
grp_num <- as.numeric(meta$group == "autism")
sex_num <- as.numeric(meta$sex == "male")
for (k in 1:5) {
  cat(" PC", k, ": group r =", round(cor(pca$x[, k], grp_num), 2),
      "| age r =", round(cor(pca$x[, k], meta$age_years, use = "complete.obs"), 2),
      "| sex r =", round(cor(pca$x[, k], sex_num, use = "complete.obs"), 2), "\n")
}

# sample-to-sample correlation outliers
cm <- cor(logexpr)
mean_cor <- (rowSums(cm) - 1) / (ncol(cm) - 1)
cat("\nLowest mean sample correlations (possible outliers):\n")
print(round(sort(mean_cor)[1:5], 3))

# ------------------------------------------------------------------
# 7. Save cleaned data (never touch data_raw)
# ------------------------------------------------------------------
write.csv(meta, file.path(clean_dir, "discovery_sample_meta.csv"), row.names = FALSE)
write.csv(logexpr, file.path(clean_dir, "discovery_log2rpkm_filtered.csv"))
write.csv(pca$x[, 1:5], file.path(res_dir, "02_pca_scores.csv"))

cat("\n=========== PLAIN-ENGLISH SUMMARY ===========\n")
cat("Samples:", sum(meta$group == "autism"), "autism,", sum(meta$group == "control"), "control, all human PFC.\n")
cat("Check: if PC1 correlates ~1 with 'group', batch (different GEO series) is\n")
cat("confounded with diagnosis and DE results cannot be trusted as-is.\n")
