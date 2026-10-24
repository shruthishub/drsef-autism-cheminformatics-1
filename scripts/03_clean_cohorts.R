# 03_clean_cohorts.R
# Step: turn each bulk cohort into the same simple format so each can be
# analysed separately and later combined by meta-analysis (script 05).
#   data_clean/<cohort>_expr.csv  genes (HGNC symbol) x samples
#   data_clean/<cohort>_meta.csv  one row per sample
#
# Why separate cohorts: each study is processed in its own lab/pipeline, so we
# compare autism vs control INSIDE a study, never across studies.
#
# Run:  Rscript scripts/03_clean_cohorts.R

suppressMessages({
  library(GEOquery)
  library(readxl)
  library(org.Hs.eg.db)   # Ensembl ID -> gene symbol for the RNA-seq sets
  library(AnnotationDbi)
})

raw_dir   <- "data_raw"
clean_dir <- "data_clean"
res_dir   <- "results"

# ---------------------------------------------------------------
# helper: collapse many rows (probes/Ensembl IDs) to one row per gene symbol,
# keeping the row with the highest mean expression for each symbol
# ---------------------------------------------------------------
collapse_to_symbol <- function(mat, symbols) {
  ok <- !is.na(symbols) & symbols != "" & !grepl("^LOC[0-9]+$", symbols)
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]
  row_mean <- rowMeans(mat, na.rm = TRUE)
  ord <- order(row_mean, decreasing = TRUE)   # best row first
  mat <- mat[ord, , drop = FALSE]
  symbols <- symbols[ord]
  keep <- !duplicated(symbols)                # first = highest mean
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- symbols[keep]
  return(mat)
}

ensembl_to_symbol <- function(ids) {
  s <- mapIds(org.Hs.eg.db, keys = ids, keytype = "ENSEMBL", column = "SYMBOL",
              multiVals = "first")
  return(as.character(s))
}

summary_lines <- c()

# ===============================================================
# 1. GSE28521  (Illumina array, already log2 + quantile normalized by the authors,
#    cortex and cerebellum normalized separately -> keep them separate)
# ===============================================================
cat("\n--- GSE28521 ---\n")
g <- getGEO("GSE28521", GSEMatrix = TRUE, getGPL = FALSE, destdir = raw_dir)[[1]]
p <- pData(g)
e <- exprs(g)

gpl <- Table(getGEO("GPL6883", destdir = raw_dir))
probe_symbol <- gpl$Symbol[match(rownames(e), gpl$ID)]
cat("probes with a symbol:", sum(!is.na(probe_symbol)), "of", nrow(e), "\n")

meta <- data.frame(
  sample    = colnames(e),
  title     = p$title,
  donor     = sub("^[A-Z]_(.*)_[A-Z]$", "\\1", p$title),
  group     = ifelse(grepl("autism", p[["disease status:ch1"]]), "autism", "control"),
  region    = tolower(p[["tissue (brain region):ch1"]]),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample
cat("GEO has NO age/sex/RIN/PMI for this series (they are in the paper's supplement).\n")

for (reg in c("frontal cortex", "temporal cortex", "cerebellum")) {
  cat(" ", reg, ":", table(meta$group[meta$region == reg]), "(", paste(names(table(meta$group[meta$region == reg])), collapse = "/"), ")\n")
}

# cortex = frontal + temporal; a person can appear twice (one sample per region)
is_cx <- meta$region %in% c("frontal cortex", "temporal cortex")
for (set_name in c("cortex", "cerebellum")) {
  if (set_name == "cortex") use <- is_cx else use <- meta$region == "cerebellum"
  m <- meta[use, ]
  x <- e[, m$sample]
  # keep probes measured in every sample of this tissue (the authors' own
  # detection filter makes probe lists differ between cortex and cerebellum)
  full <- rowSums(is.na(x)) == 0
  x <- x[full, ]
  sym <- gpl$Symbol[match(rownames(x), gpl$ID)]
  x <- collapse_to_symbol(x, sym)
  cat(" GSE28521", set_name, ":", ncol(x), "samples,", nrow(x), "genes,",
      length(unique(m$donor)), "donors; autism/control samples:",
      sum(m$group == "autism"), "/", sum(m$group == "control"), "\n")
  write.csv(x, file.path(clean_dir, paste0("GSE28521_", set_name, "_expr.csv")))
  write.csv(m, file.path(clean_dir, paste0("GSE28521_", set_name, "_meta.csv")), row.names = FALSE)
}

# ===============================================================
# 2. GSE102741  (RNA-seq, DLPFC, log2 RPKM given by the authors)
#    xlsx columns are sample1..sample52; the GEO "description" column links
#    each GEO sample to its xlsx column name (checked above).
# ===============================================================
cat("\n--- GSE102741 ---\n")
x <- as.data.frame(read_excel(file.path(raw_dir, "GSE102741/GSE102741_log2RPKMcounts.xlsx")))
rownames(x) <- x[[1]]
x <- x[, -1]
p <- read.csv(file.path(res_dir, "GSE102741_part1_sample_table.csv"), check.names = FALSE,
              stringsAsFactors = FALSE)
meta <- data.frame(
  sample = p$description,                              # "sample1" ... "sample52"
  geo    = p[[1]],
  donor  = as.character(p[["brainid:ch1"]]),
  group  = ifelse(grepl("Autism", p[["disease status:ch1"]]), "autism", "control"),
  age    = as.numeric(p[["age:ch1"]]),
  sex    = p[["Sex:ch1"]],
  rin    = as.numeric(p[["rin:ch1"]]),
  ph     = as.numeric(p[["ph:ch1"]]),
  race   = p[["race:ch1"]],
  region = "dlpfc",
  stringsAsFactors = FALSE
)
stopifnot(all(meta$sample %in% colnames(x)))
x <- as.matrix(x[, meta$sample])
sym <- ensembl_to_symbol(rownames(x))
x <- collapse_to_symbol(x, sym)
# drop genes that are essentially off: need log2RPKM > 0.5 in at least 20% of samples
keep <- rowSums(x > 0.5) >= ceiling(0.2 * ncol(x))
x <- x[keep, ]
cat(" samples:", ncol(x), "genes kept:", nrow(x), " autism/control:",
    sum(meta$group == "autism"), "/", sum(meta$group == "control"), "\n")
cat(" xlsx value range (should look like log2):", round(range(x), 2), "\n")
write.csv(x, file.path(clean_dir, "GSE102741_expr.csv"))
write.csv(meta, file.path(clean_dir, "GSE102741_meta.csv"), row.names = FALSE)

# ===============================================================
# 3. GSE64018  (RNA-seq counts, superior temporal gyrus, 12 vs 12)
#    SENSITIVITY cohort only: samples were chosen for the largest expression
#    difference, which can exaggerate effects.
# ===============================================================
cat("\n--- GSE64018 ---\n")
cnt <- read.table(gzfile(file.path(raw_dir, "GSE64018/GSE64018_countlevel_12asd_12ctl.txt.gz")),
                  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
cat(" count table:", nrow(cnt), "genes x", ncol(cnt), "samples\n")
p <- read.csv(file.path(res_dir, "GSE64018_sample_table.csv"), row.names = 1, check.names = FALSE,
              stringsAsFactors = FALSE)
get_val <- function(txt) sub("^[^:]*: *", "", txt)
meta <- data.frame(
  donor = get_val(p[["characteristics_ch1"]]),
  group = p$group,
  age   = as.numeric(get_val(p[["characteristics_ch1.3"]])),
  sex   = get_val(p[["characteristics_ch1.4"]]),
  rin   = as.numeric(get_val(p[["characteristics_ch1.5"]])),
  bank  = get_val(p[["characteristics_ch1.6"]]),
  pmi   = as.numeric(get_val(p[["characteristics_ch1.7"]])),
  region = "superior temporal gyrus",
  stringsAsFactors = FALSE
)
# matrix columns look like AN02987_ba41-42-22_8.6 -> donor is the first piece
col_donor <- sub("_.*$", "", colnames(cnt))
idx <- match(col_donor, meta$donor)
cat(" columns matched to metadata:", sum(!is.na(idx)), "of", ncol(cnt), "\n")
stopifnot(all(!is.na(idx)))
meta <- meta[idx, ]
meta$sample <- colnames(cnt)
rownames(meta) <- meta$sample
cat(" sex coding:", unique(meta$sex), " bank:", unique(meta$bank), "\n")
sym <- ensembl_to_symbol(sub("\\..*$", "", rownames(cnt)))
cnt_sym <- collapse_to_symbol(as.matrix(cnt), sym)
# keep genes with at least 10 counts in at least 6 samples
keep <- rowSums(cnt_sym >= 10) >= 6
cnt_sym <- cnt_sym[keep, ]
cat(" samples:", ncol(cnt_sym), "genes kept:", nrow(cnt_sym), " autism/control:",
    sum(meta$group == "autism"), "/", sum(meta$group == "control"), "\n")
write.csv(cnt_sym, file.path(clean_dir, "GSE64018_counts.csv"))
write.csv(meta, file.path(clean_dir, "GSE64018_meta.csv"), row.names = FALSE)

# ===============================================================
# 4. Donor overlap, with IDs normalised (strip AN / UMB prefixes, leading zeros)
# ===============================================================
cat("\n--- Donor overlap (normalised IDs) ---\n")
norm_id <- function(x) {
  x <- toupper(as.character(x))
  x <- sub("^(AN|UMB|UM)", "", x)
  x <- sub("^0+", "", x)
  return(x)
}
m28521 <- read.csv(file.path(clean_dir, "GSE28521_cortex_meta.csv"))
m28521all <- unique(c(m28521$donor, read.csv(file.path(clean_dir, "GSE28521_cerebellum_meta.csv"))$donor))
vm <- read.delim(file.path(raw_dir, "velmeshev2019/meta.tsv"), check.names = FALSE)
p298 <- read.csv(file.path(res_dir, "GSE298243_sample_table.csv"), row.names = 1, check.names = FALSE)
d298 <- unique(c(sub("donor: ", "", p298[["characteristics_ch1.1"]]), sub("donor alt: ", "", p298[["characteristics_ch1.2"]])))
p38322 <- read.csv(file.path(res_dir, "GSE38322_sample_table.csv"), row.names = 1, check.names = FALSE)
d38322 <- sub("atp_id: ", "", p38322[["characteristics_ch1.1"]])

lists <- list(GSE28521 = norm_id(m28521all),
              GSE38322 = norm_id(d38322),
              GSE64018 = norm_id(meta$donor),
              GSE102741 = norm_id(read.csv(file.path(clean_dir, "GSE102741_meta.csv"))$donor),
              GSE298243 = norm_id(d298),
              Velmeshev = norm_id(unique(vm$individual)))
rows <- list()
nm <- names(lists)
for (a in 1:(length(lists) - 1)) {
  for (b in (a + 1):length(lists)) {
    sh <- intersect(lists[[a]], lists[[b]])
    cat(" ", nm[a], "vs", nm[b], ":", length(sh), "shared", if (length(sh) > 0) paste0("(", paste(head(sh, 10), collapse = ","), ")") else "", "\n")
    rows[[length(rows) + 1]] <- data.frame(a = nm[a], b = nm[b], n_shared = length(sh), ids = paste(sh, collapse = ";"))
  }
}
write.csv(do.call(rbind, rows), file.path(res_dir, "03_donor_overlap_normalised.csv"), row.names = FALSE)

cat("\n=========== PLAIN-ENGLISH SUMMARY ===========\n")
cat("Each bulk cohort now has its own expression table (gene symbols) and metadata in data_clean/.\n")
cat("CHECK: GSE28521 has no age/sex in GEO, so it will be adjusted with surrogate variables.\n")
cat("CHECK: overlapping donors listed above must not be counted twice in the meta-analysis.\n")
