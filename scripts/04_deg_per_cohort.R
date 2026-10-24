# 04_deg_per_cohort.R
# Aim 1, part 1: differential expression (autism vs control) INSIDE each cohort.
#
# Packages and why:
#   limma  - linear models for expression; works for arrays and (with voom) RNA-seq
#   edgeR  - library-size normalisation (TMM) for the count data (GSE64018)
#   sva    - surrogate variable analysis; GSE28521 has no age/sex/RIN in GEO, so we
#            estimate hidden technical/biological variation from the data instead
#
# Model per cohort (group = autism vs control is the effect we test):
#   GSE28521 cortex : ~ group + region + SV1..SVk   (+ donor correlation, because the
#                     same person can give a frontal and a temporal sample)
#   GSE102741       : ~ group + age + sex + RIN + race
#   GSE64018        : ~ group + age + sex + RIN + PMI + brain bank   (voom on counts)
#
# Run:  Rscript scripts/04_deg_per_cohort.R

suppressMessages({
  library(limma)
  library(edgeR)
  library(sva)
})

clean_dir <- "data_clean"
res_dir   <- "results"
fig_dir   <- "figures"
set.seed(1)

# ---------------------------------------------------------------
# helper: pull a tidy table out of a fitted limma model
# se = standard error of the log fold change (needed for meta-analysis)
# ---------------------------------------------------------------
tidy_fit <- function(fit, coef_name) {
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  se <- fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post)
  out <- data.frame(gene = rownames(tt), logFC = tt$logFC, se = se, t = tt$t,
                    P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
                    AveExpr = tt$AveExpr, stringsAsFactors = FALSE)
  return(out)
}

# helper: volcano plot (base R)
make_volcano <- function(tab, title, file) {
  png(file, width = 900, height = 700)
  sig <- tab$adj.P.Val < 0.05
  col <- rep("grey70", nrow(tab))
  col[sig & tab$logFC > 0] <- "red"
  col[sig & tab$logFC < 0] <- "blue"
  plot(tab$logFC, -log10(tab$P.Value), pch = 19, cex = 0.5, col = col,
       xlab = "log2 fold change (autism - control)", ylab = "-log10 p-value", main = title)
  legend("topright", legend = c("up in autism (FDR<0.05)", "down in autism (FDR<0.05)"),
         col = c("red", "blue"), pch = 19)
  dev.off()
}

summary_rows <- list()

report <- function(name, tab, n_aut, n_ctl) {
  n_up <- sum(tab$adj.P.Val < 0.05 & tab$logFC > 0)
  n_dn <- sum(tab$adj.P.Val < 0.05 & tab$logFC < 0)
  cat(" ", name, ": genes tested =", nrow(tab), "| FDR<0.05: up =", n_up, ", down =", n_dn,
      "| smallest FDR =", signif(min(tab$adj.P.Val), 3), "\n")
  # lambda: genomic-inflation-like check on the p-values. Near 1 is healthy.
  lambda <- median(qchisq(1 - tab$P.Value, 1)) / qchisq(0.5, 1)
  cat("     inflation lambda =", round(lambda, 2), "(1 = no inflation)\n")
  summary_rows[[length(summary_rows) + 1]] <<- data.frame(
    cohort = name, n_autism = n_aut, n_control = n_ctl, genes = nrow(tab),
    up_fdr05 = n_up, down_fdr05 = n_dn, lambda = round(lambda, 2))
}

# ===============================================================
# 1. GSE28521 cortex (microarray)
# ===============================================================
cat("\n--- GSE28521 cortex ---\n")
x <- as.matrix(read.csv(file.path(clean_dir, "GSE28521_cortex_expr.csv"), row.names = 1, check.names = FALSE))
m <- read.csv(file.path(clean_dir, "GSE28521_cortex_meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(m$sample == colnames(x)))
m$group  <- factor(m$group, levels = c("control", "autism"))
m$region <- factor(m$region)

# surrogate variables: protect the group effect (mod) and region, estimate the rest
mod  <- model.matrix(~ group + region, data = m)
mod0 <- model.matrix(~ region, data = m)
n_sv <- num.sv(x, mod, method = "be", B = 20, seed = 1)
n_sv <- min(n_sv, 5)       # cap so we do not eat all degrees of freedom
cat(" surrogate variables used:", n_sv, "\n")
if (n_sv > 0) {
  sv <- sva(x, mod, mod0, n.sv = n_sv)$sv
  colnames(sv) <- paste0("SV", 1:n_sv)
  design <- cbind(mod, sv)
} else {
  design <- mod
}
# do the surrogate variables just copy the diagnosis? check, to be safe
if (n_sv > 0) {
  for (k in 1:n_sv) {
    cat("   cor(SV", k, ", group) = ", round(cor(sv[, k], as.numeric(m$group == "autism")), 2), "\n", sep = "")
  }
}
# same person can have 2 cortical samples -> estimate within-donor correlation
dc <- duplicateCorrelation(x, design, block = m$donor)
cat(" within-donor correlation:", round(dc$consensus.correlation, 3), "\n")
fit <- lmFit(x, design, block = m$donor, correlation = dc$consensus.correlation)
fit <- eBayes(fit, robust = TRUE)
tab28521 <- tidy_fit(fit, "groupautism")
write.csv(tab28521, file.path(res_dir, "04_deg_GSE28521_cortex.csv"), row.names = FALSE)
make_volcano(tab28521, "GSE28521 cortex: autism vs control", file.path(fig_dir, "04_volcano_GSE28521_cortex.png"))
report("GSE28521_cortex", tab28521, sum(m$group == "autism"), sum(m$group == "control"))

# ===============================================================
# 2. GSE102741 (RNA-seq, log2 RPKM)
# ===============================================================
cat("\n--- GSE102741 ---\n")
x <- as.matrix(read.csv(file.path(clean_dir, "GSE102741_expr.csv"), row.names = 1, check.names = FALSE))
m <- read.csv(file.path(clean_dir, "GSE102741_meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(m$sample == colnames(x)))
m$group <- factor(m$group, levels = c("control", "autism"))
m$sex   <- factor(m$sex)
m$race  <- factor(m$race)
cat(" covariate balance - age mean autism/control:",
    round(tapply(m$age, m$group, mean), 1), "| RIN:", round(tapply(m$rin, m$group, mean), 2), "\n")
design <- model.matrix(~ group + age + sex + rin + race, data = m)
fit <- lmFit(x, design)
# trend=TRUE because log2 RPKM variance depends on the mean expression
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
tab102741 <- tidy_fit(fit, "groupautism")
write.csv(tab102741, file.path(res_dir, "04_deg_GSE102741.csv"), row.names = FALSE)
make_volcano(tab102741, "GSE102741 DLPFC: autism vs control", file.path(fig_dir, "04_volcano_GSE102741.png"))
report("GSE102741", tab102741, sum(m$group == "autism"), sum(m$group == "control"))

# ===============================================================
# 3. GSE64018 (RNA-seq counts) - SENSITIVITY cohort
# ===============================================================
cat("\n--- GSE64018 (sensitivity only) ---\n")
cnt <- as.matrix(read.csv(file.path(clean_dir, "GSE64018_counts.csv"), row.names = 1, check.names = FALSE))
m <- read.csv(file.path(clean_dir, "GSE64018_meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(m$sample == colnames(cnt)))
m$group <- factor(m$group, levels = c("control", "autism"))
m$sex   <- factor(m$sex)
m$bank  <- factor(m$bank)
dge <- DGEList(counts = cnt)
dge <- calcNormFactors(dge)                 # TMM normalisation
design <- model.matrix(~ group + age + sex + rin + pmi + bank, data = m)
v <- voom(dge, design)                      # turns counts into log-CPM with precision weights
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)
tab64018 <- tidy_fit(fit, "groupautism")
write.csv(tab64018, file.path(res_dir, "04_deg_GSE64018.csv"), row.names = FALSE)
make_volcano(tab64018, "GSE64018 STG: autism vs control", file.path(fig_dir, "04_volcano_GSE64018.png"))
report("GSE64018", tab64018, sum(m$group == "autism"), sum(m$group == "control"))

# ===============================================================
# 4. Quick cross-cohort look (full meta-analysis is script 05)
# ===============================================================
summ <- do.call(rbind, summary_rows)
write.csv(summ, file.path(res_dir, "04_deg_summary.csv"), row.names = FALSE)

cat("\n--- Direction/rank agreement between cohorts (genes in both) ---\n")
tabs <- list(GSE28521_cortex = tab28521, GSE102741 = tab102741, GSE64018 = tab64018)
nm <- names(tabs)
for (a in 1:(length(nm) - 1)) {
  for (b in (a + 1):length(nm)) {
    shared <- intersect(tabs[[a]]$gene, tabs[[b]]$gene)
    ta <- tabs[[a]][match(shared, tabs[[a]]$gene), ]
    tb <- tabs[[b]][match(shared, tabs[[b]]$gene), ]
    r_s <- cor(ta$t, tb$t, method = "spearman")
    same_dir <- mean(sign(ta$t) == sign(tb$t))
    cat(" ", nm[a], "vs", nm[b], ": shared genes =", length(shared),
        "| Spearman r of t-stats =", round(r_s, 3), "| same direction =", round(100 * same_dir, 1), "%\n")
  }
}
cat("\n=========== PLAIN-ENGLISH SUMMARY ===========\n")
print(summ)
cat("CHECK: lambda far above 1 means leftover confounding; r near 0 between cohorts means no replication.\n")
