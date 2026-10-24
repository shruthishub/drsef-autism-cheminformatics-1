# 05_replication_and_meta.R
# Aim 1, part 2: do the differential-expression results replicate, and what does
# a combined (meta-analysis) result look like?
#
# Design decisions (see notes):
#   PRIMARY cohorts (no shared donors): GSE28521 cortex (29/29 samples, 32 donors)
#                                       GSE102741 DLPFC (13/39)
#   SENSITIVITY cohort: GSE64018 (shares 11 donors with GSE28521 and was selected
#                       for large expression differences) -> never used as evidence
#                       of independent replication.
#
# Meta-analysis method: weighted Stouffer Z. Each cohort uses a different platform
# (array vs RNA-seq) so log fold changes are on different scales; z-scores from
# each cohort's moderated t-test are comparable. Weight = sqrt(sample size).
#
# Packages: fgsea (pathway enrichment on ranked genes), msigdbr (gene set lists)
#
# Run:  Rscript scripts/05_replication_and_meta.R

suppressMessages({
  library(fgsea)
  library(msigdbr)
})

res_dir <- "results"
fig_dir <- "figures"
set.seed(1)

# t-statistic -> signed z-score using the p-value (works for any df)
signed_z <- function(tab) {
  z <- qnorm(1 - tab$P.Value / 2) * sign(tab$t)
  z[is.infinite(z)] <- sign(z[is.infinite(z)]) * 8
  return(z)
}

a <- read.csv(file.path(res_dir, "04_deg_GSE28521_cortex.csv"), stringsAsFactors = FALSE)
b <- read.csv(file.path(res_dir, "04_deg_GSE102741.csv"), stringsAsFactors = FALSE)
c3 <- read.csv(file.path(res_dir, "04_deg_GSE64018.csv"), stringsAsFactors = FALSE)
a$z <- signed_z(a); b$z <- signed_z(b); c3$z <- signed_z(c3)
n_a <- 58; n_b <- 52; n_c <- 24    # samples per cohort (donors for GSE28521 = 32)

# ===============================================================
# 1. REPLICATION of GSE28521 hits in GSE102741
# ===============================================================
cat("\n=== 1. Replication: GSE28521 cortex (discovery) -> GSE102741 (replication) ===\n")
shared <- intersect(a$gene, b$gene)
A <- a[match(shared, a$gene), ]
B <- b[match(shared, b$gene), ]
cat("Genes measured in both:", length(shared), "\n")

hit <- A$adj.P.Val < 0.05
cat("Discovery hits (FDR<0.05 in GSE28521):", sum(hit), "\n")
same_dir <- sign(A$z[hit]) == sign(B$z[hit])
cat(" Same direction in GSE102741:", sum(same_dir), "of", sum(hit),
    "(", round(100 * mean(same_dir), 1), "%; 50% expected by chance)\n")
bt <- binom.test(sum(same_dir), sum(hit), p = 0.5)
cat(" Binomial test p =", signif(bt$p.value, 3), "\n")
# nominal replication: same direction AND p<0.05 in the replication cohort
rep_nom <- same_dir & B$P.Value[hit] < 0.05
cat(" Same direction AND p<0.05 in GSE102741:", sum(rep_nom), "\n")
# how many would be expected by chance? permute the replication cohort's gene labels
n_perm <- 2000
perm_counts <- numeric(n_perm)
hit_idx <- which(hit)
for (i in 1:n_perm) {
  shuf <- sample(length(shared))
  zB <- B$z[shuf][hit_idx]
  pB <- B$P.Value[shuf][hit_idx]
  perm_counts[i] <- sum(sign(A$z[hit_idx]) == sign(zB) & pB < 0.05)
}
perm_p <- (sum(perm_counts >= sum(rep_nom)) + 1) / (n_perm + 1)
cat(" Expected by chance (permutation mean):", round(mean(perm_counts), 1),
    "| permutation p =", signif(perm_p, 3), "\n")

# whole-transcriptome correlation of effects
r_all <- cor(A$z, B$z, method = "spearman")
r_hit <- cor(A$z[hit], B$z[hit], method = "spearman")
cat(" Spearman r of z across all shared genes =", round(r_all, 3), "| among discovery hits =", round(r_hit, 3), "\n")

# the same check but using ONLY the top 200 discovery genes by p-value
top <- order(A$P.Value)[1:200]
cat(" Top-200 discovery genes: same direction in replication =",
    round(100 * mean(sign(A$z[top]) == sign(B$z[top])), 1), "%\n")

png(file.path(fig_dir, "05_replication_scatter.png"), width = 800, height = 800)
plot(A$z, B$z, pch = 19, cex = 0.4, col = ifelse(hit, "red", "grey70"),
     xlab = "GSE28521 cortex z (autism vs control)", ylab = "GSE102741 DLPFC z",
     main = paste0("Replication of effects (Spearman r = ", round(r_all, 2), ")"))
abline(h = 0, v = 0, lty = 2)
legend("topleft", legend = "FDR<0.05 in GSE28521", col = "red", pch = 19)
dev.off()

# ===============================================================
# 2. META-ANALYSIS (weighted Stouffer)
# ===============================================================
cat("\n=== 2. Meta-analysis ===\n")
w_a <- sqrt(n_a); w_b <- sqrt(n_b); w_c <- sqrt(n_c)

meta_two <- data.frame(gene = shared, z_28521 = A$z, z_102741 = B$z, stringsAsFactors = FALSE)
meta_two$z_meta <- (w_a * A$z + w_b * B$z) / sqrt(w_a^2 + w_b^2)
meta_two$P <- 2 * pnorm(-abs(meta_two$z_meta))
meta_two$FDR <- p.adjust(meta_two$P, method = "BH")
meta_two$logFC_28521 <- A$logFC
meta_two$logFC_102741 <- B$logFC
meta_two$same_direction <- sign(meta_two$z_28521) == sign(meta_two$z_102741)
meta_two$p_28521 <- A$P.Value
meta_two$p_102741 <- B$P.Value
# "replicated gene" = meta FDR<0.05 AND same direction in both AND p<0.05 in both
meta_two$replicated <- meta_two$FDR < 0.05 & meta_two$same_direction &
                       meta_two$p_28521 < 0.05 & meta_two$p_102741 < 0.05
meta_two <- meta_two[order(meta_two$P), ]
write.csv(meta_two, file.path(res_dir, "05_meta_primary_2cohorts.csv"), row.names = FALSE)
cat("Primary meta-analysis (GSE28521 + GSE102741):", nrow(meta_two), "genes\n")
cat(" meta FDR<0.05:", sum(meta_two$FDR < 0.05),
    "| of those same direction in both:", sum(meta_two$FDR < 0.05 & meta_two$same_direction),
    "| fully replicated (p<0.05 in both):", sum(meta_two$replicated), "\n")
cat(" replicated, up in autism:", sum(meta_two$replicated & meta_two$z_meta > 0),
    "| down in autism:", sum(meta_two$replicated & meta_two$z_meta < 0), "\n")
cat("Top 15 replicated genes:\n")
print(head(meta_two[meta_two$replicated, c("gene", "z_28521", "z_102741", "z_meta", "FDR")], 15), row.names = FALSE)

# leave-one-out: how well does each single cohort alone agree with the meta z?
# (with 2 cohorts this is just the two cohorts themselves; reported for clarity)

# sensitivity: add GSE64018 as a third cohort
sh3 <- Reduce(intersect, list(a$gene, b$gene, c3$gene))
A3 <- a[match(sh3, a$gene), ]; B3 <- b[match(sh3, b$gene), ]; C3 <- c3[match(sh3, c3$gene), ]
m3 <- data.frame(gene = sh3, z_28521 = A3$z, z_102741 = B3$z, z_64018 = C3$z)
m3$z_meta <- (w_a * A3$z + w_b * B3$z + w_c * C3$z) / sqrt(w_a^2 + w_b^2 + w_c^2)
m3$P <- 2 * pnorm(-abs(m3$z_meta))
m3$FDR <- p.adjust(m3$P, method = "BH")
m3 <- m3[order(m3$P), ]
write.csv(m3, file.path(res_dir, "05_meta_sensitivity_3cohorts.csv"), row.names = FALSE)
cat("\nSensitivity meta (adds GSE64018, NOT independent of GSE28521):", nrow(m3), "genes, FDR<0.05 =", sum(m3$FDR < 0.05), "\n")
ov <- intersect(meta_two$gene[meta_two$FDR < 0.05], m3$gene[m3$FDR < 0.05])
cat(" overlap with primary meta FDR<0.05 genes:", length(ov), "\n")

# ===============================================================
# 3. Pathways: which biological programs are shifted in autism?
#    fgsea on the meta z-scores (ranked list), Hallmark + GO biological process
# ===============================================================
cat("\n=== 3. Pathway enrichment on the ranked meta z-scores ===\n")
ranks <- meta_two$z_meta
names(ranks) <- meta_two$gene
ranks <- ranks[!duplicated(names(ranks))]

get_sets <- function(coll, sub = NULL) {
  df <- if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
        else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  sets <- split(df$gene_symbol, df$gs_name)
  return(sets)
}
hall <- get_sets("H")
gobp <- get_sets("C5", "GO:BP")
cat("Hallmark sets:", length(hall), "| GO:BP sets:", length(gobp), "\n")

run_fgsea <- function(sets, label) {
  fg <- fgsea(pathways = sets, stats = ranks, minSize = 15, maxSize = 500, eps = 0)
  fg <- fg[order(fg$padj), ]
  out <- data.frame(pathway = fg$pathway, size = fg$size, NES = fg$NES, pval = fg$pval, padj = fg$padj)
  write.csv(out, file.path(res_dir, paste0("05_fgsea_", label, ".csv")), row.names = FALSE)
  return(out)
}
fh <- run_fgsea(hall, "hallmark")
fg <- run_fgsea(gobp, "GOBP")
cat("\nTop Hallmark sets (NES>0 = up in autism):\n")
print(head(fh[, c("pathway", "NES", "padj")], 10), row.names = FALSE)
cat("\nTop GO:BP sets, up in autism:\n")
print(head(fg[fg$NES > 0, c("pathway", "NES", "padj")], 8), row.names = FALSE)
cat("\nTop GO:BP sets, down in autism:\n")
print(head(fg[fg$NES < 0, c("pathway", "NES", "padj")], 8), row.names = FALSE)

# dedicated checks that match the project hypotheses
cat("\nHypothesis-specific sets:\n")
for (nm in c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")) {
  row <- fh[fh$pathway == nm, ]
  cat(" ", nm, ": NES =", round(row$NES, 2), ", padj =", signif(row$padj, 3), "\n")
}
syn <- fg[grepl("SYNAP", fg$pathway), ]
cat(" synapse-related GO:BP sets tested:", nrow(syn), "| significantly DOWN (padj<0.05):", sum(syn$padj < 0.05 & syn$NES < 0), "| UP:", sum(syn$padj < 0.05 & syn$NES > 0), "\n")

cat("\n=========== PLAIN-ENGLISH SUMMARY ===========\n")
cat("See the counts above. If the same-direction rate is near 50% and the permutation\n")
cat("p-value is not small, GSE28521 effects do NOT replicate in GSE102741, and that is a result.\n")
