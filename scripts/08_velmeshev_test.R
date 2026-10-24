# 08_velmeshev_test.R
# Third, technologically different cohort (single-nucleus RNA-seq, 15 ASD / 16 control
# people, PFC + anterior cingulate, ages 4-22). Questions:
#   (a) does the total-pseudobulk autism-vs-control signal agree with GSE28521?
#   (b) do the GSE28521 WGCNA modules differ by diagnosis here too?
#   (c) are cell-type proportions different between autism and control?
#
# Packages: limma + edgeR (voom for counts), same donor-correlation model as before
# because one person can contribute both a PFC and an ACC sample.
#
# Run:  Rscript scripts/08_velmeshev_test.R

suppressMessages({
  library(limma)
  library(edgeR)
  library(WGCNA)
})
set.seed(1)
clean_dir <- "data_clean"
res_dir   <- "results"
fig_dir   <- "figures"

signed_z <- function(P, t) {
  z <- qnorm(1 - P / 2) * sign(t)
  z[is.infinite(z)] <- sign(z[is.infinite(z)]) * 8
  return(z)
}

# ===============================================================
# (a) total pseudobulk differential expression
# ===============================================================
cnt <- as.matrix(read.csv(file.path(clean_dir, "velmeshev_pb_total_counts.csv"), row.names = 1, check.names = FALSE))
m <- read.csv(file.path(clean_dir, "velmeshev_pb_total_meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(m$sample == colnames(cnt)))
m$group  <- factor(ifelse(m$diagnosis == "ASD", "autism", "control"), levels = c("control", "autism"))
m$region <- factor(m$region)
m$sex    <- factor(m$sex)
m$seqbatch <- factor(m$seqbatch)
cat("samples:", nrow(m), "| ASD/control samples:", sum(m$group == "autism"), "/", sum(m$group == "control"),
    "| individuals:", length(unique(m$individual)), "\n")
cat("cells per sample: min", min(m$n_cells), "median", median(m$n_cells), "max", max(m$n_cells), "\n")
cat("RIN missing for", sum(is.na(m$rin)), "samples; PMI missing for", sum(is.na(m$pmi)), "\n")
cat("age mean ASD/control:", round(tapply(m$age, m$group, mean), 1), "| sex counts:\n"); print(table(m$group, m$sex))
cat("sequencing batches by group (confounding check):\n"); print(table(m$group, m$seqbatch))

keep <- rowSums(cnt >= 10) >= 8
cnt <- cnt[keep, ]
cat("genes kept:", nrow(cnt), "\n")
dge <- DGEList(counts = cnt)
dge <- calcNormFactors(dge)
# covariates: region, age, sex, sequencing batch (RIN/PMI are partly missing so left out)
design <- model.matrix(~ group + region + age + sex + seqbatch, data = m)
v <- voom(dge, design)
dc <- duplicateCorrelation(v, design, block = m$individual)
cat("within-person correlation:", round(dc$consensus.correlation, 3), "\n")
v <- voom(dge, design, block = m$individual, correlation = dc$consensus.correlation)
fit <- lmFit(v, design, block = m$individual, correlation = dc$consensus.correlation)
fit <- eBayes(fit, robust = TRUE)
tt <- topTable(fit, coef = "groupautism", number = Inf, sort.by = "none")
tabV <- data.frame(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
tabV$z <- signed_z(tabV$P.Value, tabV$t)
write.csv(tabV, file.path(res_dir, "08_deg_velmeshev_total.csv"), row.names = FALSE)
lambda <- median(qchisq(1 - tabV$P.Value, 1)) / qchisq(0.5, 1)
cat("Velmeshev total pseudobulk: FDR<0.05 up =", sum(tabV$adj.P.Val < 0.05 & tabV$logFC > 0),
    "down =", sum(tabV$adj.P.Val < 0.05 & tabV$logFC < 0), "| lambda =", round(lambda, 2), "\n")

# compare with the bulk cohorts
a <- read.csv(file.path(res_dir, "04_deg_GSE28521_cortex.csv"))
b <- read.csv(file.path(res_dir, "04_deg_GSE102741.csv"))
a$z <- signed_z(a$P.Value, a$t); b$z <- signed_z(b$P.Value, b$t)
for (nm in c("GSE28521", "GSE102741")) {
  other <- if (nm == "GSE28521") a else b
  sh <- intersect(tabV$gene, other$gene)
  zv <- tabV$z[match(sh, tabV$gene)]
  zo <- other$z[match(sh, other$gene)]
  cat(" Velmeshev vs", nm, ": shared genes =", length(sh), "| Spearman r of z =", round(cor(zv, zo, method = "spearman"), 3),
      "| same direction =", round(100 * mean(sign(zv) == sign(zo)), 1), "%\n")
  hit <- other$adj.P.Val[match(sh, other$gene)] < 0.05
  if (sum(hit) > 0) {
    sd <- sign(zv[hit]) == sign(zo[hit])
    bt <- binom.test(sum(sd), sum(hit), 0.5)
    cat("   among", sum(hit), nm, "hits: same direction in Velmeshev =", sum(sd), "(", round(100 * mean(sd), 1), "%), binomial p =", signif(bt$p.value, 3), "\n")
  }
}

# ===============================================================
# (b) GSE28521 WGCNA modules, tested inside Velmeshev total pseudobulk
# ===============================================================
w <- readRDS(file.path(clean_dir, "wgcna_gse28521.rds"))
colors <- w$colors; names(colors) <- w$genes
logcpm <- v$E                                       # voom log2-CPM
shared <- intersect(names(colors), rownames(logcpm))
cat("\nnetwork genes found in Velmeshev:", length(shared), "of", length(colors), "\n")
# remove known covariates, keep diagnosis
x_adj <- removeBatchEffect(logcpm[shared, ], batch = m$seqbatch,
                           covariates = cbind(age = m$age, male = as.numeric(m$sex == "M"), acc = as.numeric(m$region == "ACC")),
                           design = model.matrix(~ group, data = m))
me <- moduleEigengenes(t(x_adj), colors = colors[shared], excludeGrey = TRUE)$eigengenes
colnames(me) <- sub("^ME", "", colnames(me))
dm <- model.matrix(~ group, data = m)
dcm <- duplicateCorrelation(t(me), dm, block = m$individual)
fm <- lmFit(t(me), dm, block = m$individual, correlation = dcm$consensus.correlation)
fm <- eBayes(fm)
mod_v <- data.frame(module = rownames(fm$coefficients), diff_ME_velm = fm$coefficients[, "groupautism"],
                    t_velm = fm$t[, "groupautism"], P_velm = fm$p.value[, "groupautism"])
mod_v$FDR_velm <- p.adjust(mod_v$P_velm, method = "BH")
prev <- read.csv(file.path(res_dir, "07_module_trait_and_replication.csv"))
mod_all <- merge(prev[, c("module", "n_genes", "diff_ME", "t", "P", "FDR", "P_102741", "Zsummary")], mod_v, by = "module")
mod_all$same_dir_velm <- sign(mod_all$diff_ME) == sign(mod_all$diff_ME_velm)
mod_all <- mod_all[order(mod_all$P), ]
write.csv(mod_all, file.path(res_dir, "08_module_replication_velmeshev.csv"), row.names = FALSE)
cat("\nModule effects: GSE28521 vs GSE102741 vs Velmeshev\n")
print(mod_all[, c("module", "n_genes", "t", "FDR", "P_102741", "t_velm", "P_velm", "same_dir_velm", "Zsummary")], row.names = FALSE, digits = 3)
cat("Correlation of module t-statistics, GSE28521 vs Velmeshev: r =", round(cor(mod_all$t, mod_all$t_velm), 2), "\n")
sig <- mod_all$FDR < 0.05
cat("Of", sum(sig), "modules significant in GSE28521: same direction in Velmeshev =", sum(mod_all$same_dir_velm[sig]),
    "| p<0.05 in Velmeshev =", sum(mod_all$P_velm[sig] < 0.05 & mod_all$same_dir_velm[sig]), "\n")

# ===============================================================
# (c) cell-type composition (fraction of each sample's cells), autism vs control
# ===============================================================
pbm <- read.csv(file.path(clean_dir, "velmeshev_pb_meta.csv"), stringsAsFactors = FALSE)
types <- sort(unique(pbm$celltype))
frac <- matrix(0, nrow = length(types), ncol = nrow(m), dimnames = list(types, m$sample))
for (i in 1:nrow(pbm)) {
  frac[pbm$celltype[i], pbm$sample[i]] <- pbm$n_cells[i]
}
frac <- sweep(frac, 2, colSums(frac), "/")
eps <- 0.001
logit_frac <- log((frac + eps) / (1 - frac + eps))
dcc <- duplicateCorrelation(logit_frac, dm, block = m$individual)
fc <- lmFit(logit_frac, design = model.matrix(~ group + region, data = m), block = m$individual, correlation = dcc$consensus.correlation)
fc <- eBayes(fc)
comp <- data.frame(celltype = rownames(fc$coefficients),
                   mean_frac_control = rowMeans(frac[, m$group == "control"]),
                   mean_frac_autism = rowMeans(frac[, m$group == "autism"]),
                   t = fc$t[, "groupautism"], P = fc$p.value[, "groupautism"])
comp$FDR <- p.adjust(comp$P, method = "BH")
comp <- comp[order(comp$P), ]
write.csv(comp, file.path(res_dir, "08_celltype_composition.csv"), row.names = FALSE)
cat("\nCell-type proportion differences (autism vs control):\n")
print(comp, row.names = FALSE, digits = 3)

cat("\n=========== PLAIN-ENGLISH SUMMARY ===========\n")
cat("Compare r and same-direction % above with the GSE28521 vs GSE102741 values (r=0.18, 55%).\n")
cat("CHECK: Velmeshev donors are only 4-22 years old, so agreement is partial at best.\n")
