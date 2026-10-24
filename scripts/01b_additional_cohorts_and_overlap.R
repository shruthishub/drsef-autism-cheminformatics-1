# 01b_additional_cohorts_and_overlap.R
# Step: collect sample tables for the extra human brain cohorts found in the GEO
# search, count autism/control, and check whether the SAME PEOPLE appear in more
# than one dataset (donors reused across studies are not independent replication).
#
# Run:  Rscript scripts/01b_additional_cohorts_and_overlap.R

suppressMessages(library(GEOquery))

raw_dir <- "data_raw"
res_dir <- "results"

# helper: get the sample table of one GEO series (first platform / part j)
get_pheno <- function(acc, destdir) {
  g <- getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE, destdir = destdir)
  out <- list()
  for (j in 1:length(g)) {
    out[[j]] <- pData(g[[j]])
  }
  return(out)
}

# helper: paste all text columns of a sample into one string for regex search
row_text <- function(p) {
  keep <- grep("characteristics|title|source_name|description", colnames(p))
  txt <- character(nrow(p))
  for (i in 1:nrow(p)) {
    txt[i] <- paste(as.character(unlist(p[i, keep])), collapse = " | ")
  }
  return(txt)
}

# ---------------------------------------------------------------
# 1. GSE38322 (Illumina HT12 microarray, ATP brain bank) and GSE64018 (RNA-seq)
# ---------------------------------------------------------------
p38322 <- get_pheno("GSE38322", file.path(raw_dir, "GSE38322"))[[1]]
t38322 <- row_text(p38322)
p38322$group <- ifelse(grepl("^Autism|autism", p38322$title) | grepl("Autism_", p38322$source_name_ch1), "autism", "control")
p38322$region <- sub("brain region: ", "", p38322[["characteristics_ch1"]])
write.csv(p38322, file.path(res_dir, "GSE38322_sample_table.csv"))
cat("\nGSE38322: group x region\n")
print(table(p38322$group, p38322$region))
p38322$donor <- sub("atp_id: ", "", p38322[["characteristics_ch1.1"]])

p64018 <- get_pheno("GSE64018", file.path(raw_dir, "GSE64018"))[[1]]
p64018$group <- ifelse(grepl("ASD", p64018[["characteristics_ch1.2"]]), "autism", "control")
p64018$donor <- sub("brainid: ", "", p64018[["characteristics_ch1"]])
write.csv(p64018, file.path(res_dir, "GSE64018_sample_table.csv"))
cat("\nGSE64018 group counts (superior temporal gyrus):\n")
print(table(p64018$group))

# ---------------------------------------------------------------
# 2. GSE298243 (single-nucleus methylation + RNA, BA9, males)
# ---------------------------------------------------------------
p298243 <- get_pheno("GSE298243", file.path(raw_dir, "candidates"))
cat("\nGSE298243: number of series-matrix parts:", length(p298243), "\n")
p298 <- p298243[[1]]
write.csv(p298, file.path(res_dir, "GSE298243_sample_table.csv"))
cat("columns:\n"); print(colnames(p298)[grepl("characteristics|title", colnames(p298))])
cat("samples:", nrow(p298), "\n")
for (cc in grep("characteristics", colnames(p298), value = TRUE)) {
  v <- unique(as.character(p298[[cc]]))
  cat(cc, ":", head(v, 4), "\n")
}

# ---------------------------------------------------------------
# 3. Collect donor IDs from every cohort and compare
# ---------------------------------------------------------------
# GSE28521: titles look like "A_AN11989_C"  -> donor AN11989
p28521 <- read.csv(file.path(res_dir, "GSE28521_part1_sample_table.csv"), row.names = 1, check.names = FALSE)
d28521 <- unique(sub("^[A-Z]_([A-Za-z0-9]+)_.*$", "\\1", p28521$title))
# GSE102741: brainid column (numbers, from the Mount Sinai / other bank)
p102741 <- read.csv(file.path(res_dir, "GSE102741_part1_sample_table.csv"), row.names = 1, check.names = FALSE)
d102741 <- unique(as.character(p102741[["brainid:ch1"]]))
# Velmeshev: individual column (numbers like 5278 -> UMB/ATP style codes)
vm <- read.delim(file.path(raw_dir, "velmeshev2019/meta.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
dvel <- unique(as.character(vm$individual))
# GSE28475 titles do not carry donor IDs (e.g. "Frozen Tissue_DASL_57B"), so skip

cat("\nDonor ID examples:\n")
cat(" GSE28521  :", head(d28521, 4), "\n")
cat(" GSE38322  :", head(unique(p38322$donor), 4), "\n")
cat(" GSE64018  :", head(unique(p64018$donor), 4), "\n")
cat(" GSE102741 :", head(d102741, 4), "\n")
cat(" Velmeshev :", head(dvel, 4), "\n")

lists <- list(GSE28521 = d28521, GSE38322 = unique(p38322$donor),
              GSE64018 = unique(p64018$donor), GSE102741 = d102741, Velmeshev = dvel)
cat("\nShared donor IDs between cohorts (exact string match):\n")
names_l <- names(lists)
overlap_rows <- list()
for (a in 1:(length(lists) - 1)) {
  for (b in (a + 1):length(lists)) {
    shared <- intersect(lists[[a]], lists[[b]])
    cat(" ", names_l[a], "vs", names_l[b], ":", length(shared), "shared",
        if (length(shared) > 0) paste0("(", paste(head(shared, 8), collapse = ", "), ")") else "", "\n")
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      a = names_l[a], b = names_l[b], n_shared = length(shared),
      shared_ids = paste(shared, collapse = ";"))
  }
}
write.csv(do.call(rbind, overlap_rows), file.path(res_dir, "01b_donor_overlap.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# 4. Velmeshev donor summary (autism vs control individuals)
# ---------------------------------------------------------------
vm_donor <- unique(vm[, c("individual", "diagnosis", "age", "sex", "region", "post-mortem interval (hours)", "RNA Integrity Number")])
cat("\nVelmeshev individuals by diagnosis:\n")
print(table(unique(vm[, c("individual", "diagnosis")])$diagnosis))
cat("Velmeshev age range (years) by diagnosis:\n")
print(tapply(vm$age, vm$diagnosis, range))
cat("Cells per diagnosis:\n"); print(table(vm$diagnosis))
cat("Cell types:\n"); print(sort(table(vm$cluster), decreasing = TRUE))
write.csv(vm_donor, file.path(res_dir, "01b_velmeshev_donor_table.csv"), row.names = FALSE)
