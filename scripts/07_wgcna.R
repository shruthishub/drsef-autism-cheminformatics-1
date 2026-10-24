# 07_wgcna.R
# Aim 2: co-expression modules (WGCNA) in the biggest clean cortex cohort
# (GSE28521, 58 cortical samples), then:
#   - which modules differ between autism and control (module eigengene test)
#   - SFARI autism-risk gene enrichment (permutation test)
#   - interferon / inflammation enrichment (MSigDB Hallmark)
#   - do the modules replicate in the independent cohort GSE102741 (module preservation
#     + eigengene test there)
#
# Packages and why:
#   WGCNA  - builds gene co-expression networks and modules, calculates preservation
#   limma  - tests module eigengenes with the same within-donor correlation model as script 04
#   sva    - same surrogate variables as script 04 (the data have no age/sex in GEO)
#
# Run:  Rscript scripts/07_wgcna.R

suppressMessages({
  library(WGCNA)
  library(limma)
  library(sva)
  library(msigdbr)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)
cor <- WGCNA::cor
set.seed(1)

clean_dir <- "data_clean"
res_dir   <- "results"
fig_dir   <- "figures"

# ===============================================================
# 1. Load GSE28521 cortex and remove nuisance variation
# ===============================================================
x <- as.matrix(read.csv(file.path(clean_dir, "GSE28521_cortex_expr.csv"), row.names = 1, check.names = FALSE))
m <- read.csv(file.path(clean_dir, "GSE28521_cortex_meta.csv"), stringsAsFactors = FALSE)
stopifnot(all(m$sample == colnames(x)))
m$group  <- factor(m$group, levels = c("control", "autism"))
m$region <- factor(m$region)

# surrogate variables (same call as script 04, so same result with seed 1)
mod  <- model.matrix(~ group + region, data = m)
mod0 <- model.matrix(~ region, data = m)
n_sv <- min(num.sv(x, mod, method = "be", B = 20, seed = 1), 5)
sv <- sva(x, mod, mod0, n.sv = n_sv)$sv
colnames(sv) <- paste0("SV", 1:n_sv)

# remove region and SVs but KEEP the diagnosis difference
x_adj <- removeBatchEffect(x, covariates = cbind(region_temporal = as.numeric(m$region == "temporal cortex"), sv),
                           design = model.matrix(~ group, data = m))
cat("expression matrix after adjustment:", nrow(x_adj), "genes x", ncol(x_adj), "samples\n")

# keep the 6000 most variable genes (WGCNA works best on varying genes, and it saves memory)
gene_var <- apply(x_adj, 1, var)
keep <- order(gene_var, decreasing = TRUE)[1:6000]
datExpr <- t(x_adj[keep, ])      # WGCNA wants samples in rows, genes in columns
gsg <- goodSamplesGenes(datExpr, verbose = 0)
cat("goodSamplesGenes all OK:", gsg$allOK, "\n")

# sample clustering to look for outliers
sampleTree <- hclust(dist(datExpr), method = "average")
png(file.path(fig_dir, "07_sample_tree.png"), width = 1200, height = 600)
plot(sampleTree, main = "GSE28521 cortex sample clustering (outlier check)", cex = 0.5)
dev.off()

# ===============================================================
# 2. Soft-threshold power
# ===============================================================
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
fit_tab <- sft$fitIndices
fit_tab$signed_R2 <- -sign(fit_tab[, 3]) * fit_tab[, 2]
print(fit_tab[, c("Power", "signed_R2", "mean.k.")], row.names = FALSE)
ok <- fit_tab$Power[fit_tab$signed_R2 > 0.8]
if (length(ok) > 0) power <- min(ok) else power <- 12   # WGCNA's own default fallback for signed networks, n<20 would use 18
cat("chosen soft power:", power, if (length(ok) == 0) "(no power reached R^2>0.8, fallback used)" else "", "\n")
png(file.path(fig_dir, "07_soft_threshold.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))
plot(fit_tab$Power, fit_tab$signed_R2, type = "b", xlab = "soft power", ylab = "scale-free fit (signed R^2)")
abline(h = 0.8, col = "red")
plot(fit_tab$Power, fit_tab$mean.k., type = "b", xlab = "soft power", ylab = "mean connectivity")
dev.off()

# ===============================================================
# 3. Build modules
# ===============================================================
net <- blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2,
                        numericLabels = FALSE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
                        maxBlockSize = 7000, verbose = 0)
colors <- net$colors
mods <- sort(setdiff(unique(colors), "grey"))
cat("\nmodules found:", length(mods), "| genes in grey (unassigned):", sum(colors == "grey"), "of", length(colors), "\n")
print(table(colors))
png(file.path(fig_dir, "07_dendrogram.png"), width = 1200, height = 600)
plotDendroAndColors(net$dendrograms[[1]], colors[net$blockGenes[[1]]], "module", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "GSE28521 cortex gene dendrogram")
dev.off()

# ===============================================================
# 4. Module eigengenes vs diagnosis (same donor-correlation model as script 04)
# ===============================================================
MEs <- net$MEs                       # samples x modules
colnames(MEs) <- sub("^ME", "", colnames(MEs))
MEs <- MEs[, colnames(MEs) != "grey", drop = FALSE]
design_me <- model.matrix(~ group, data = m)          # nuisance already removed from the data
dc <- duplicateCorrelation(t(MEs), design_me, block = m$donor)
fit_me <- lmFit(t(MEs), design_me, block = m$donor, correlation = dc$consensus.correlation)
fit_me <- eBayes(fit_me)
me_tab <- data.frame(module = rownames(fit_me$coefficients),
                     n_genes = as.numeric(table(colors)[rownames(fit_me$coefficients)]),
                     diff_ME = fit_me$coefficients[, "groupautism"],
                     t = fit_me$t[, "groupautism"],
                     P = fit_me$p.value[, "groupautism"])
me_tab$FDR <- p.adjust(me_tab$P, method = "BH")
me_tab <- me_tab[order(me_tab$P), ]
cat("\nModule eigengene vs diagnosis (positive diff_ME = higher in autism):\n")
print(me_tab, row.names = FALSE)

# ===============================================================
# 5. Module gene lists, DE direction, SFARI + interferon enrichment
# ===============================================================
deg <- read.csv(file.path(res_dir, "04_deg_GSE28521_cortex.csv"))
net_genes <- colnames(datExpr)
gene_tab <- data.frame(gene = net_genes, module = colors, stringsAsFactors = FALSE)
gene_tab$t_28521 <- deg$t[match(gene_tab$gene, deg$gene)]
write.csv(gene_tab, file.path(res_dir, "07_module_genes.csv"), row.names = FALSE)
saveRDS(list(colors = colors, MEs = MEs, genes = net_genes, power = power), file.path(clean_dir, "wgcna_gse28521.rds"))

sfari <- read.csv("data_raw/reference/sfari_genes.csv", check.names = FALSE)
sfari_all <- unique(sfari[["gene-symbol"]])
sfari_12  <- unique(sfari[["gene-symbol"]][sfari[["gene-score"]] %in% c(1, 2)])
sfari_1   <- unique(sfari[["gene-symbol"]][sfari[["gene-score"]] %in% 1])
cat("\nSFARI genes: all =", length(sfari_all), "| score 1-2 =", length(sfari_12), "| score 1 =", length(sfari_1), "\n")
cat("SFARI genes present in the 6000-gene network: all =", sum(net_genes %in% sfari_all), "| score 1-2 =", sum(net_genes %in% sfari_12), "\n")

hall <- split(msigdbr(species = "Homo sapiens", collection = "H")$gene_symbol,
              msigdbr(species = "Homo sapiens", collection = "H")$gs_name)
ifn_sets <- list(IFN_alpha = hall[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
                 IFN_gamma = hall[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
                 Inflammatory = hall[["HALLMARK_INFLAMMATORY_RESPONSE"]])

# permutation test: draw random gene sets of the same size from the network universe
perm_enrich <- function(module_genes, target, universe, n_perm = 10000) {
  target <- intersect(target, universe)
  obs <- sum(module_genes %in% target)
  rnd <- numeric(n_perm)
  for (i in 1:n_perm) {
    rnd[i] <- sum(sample(universe, length(module_genes)) %in% target)
  }
  p <- (sum(rnd >= obs) + 1) / (n_perm + 1)
  return(c(observed = obs, expected = mean(rnd), fold = obs / max(mean(rnd), 1e-9), p = p))
}
enrich_rows <- list()
for (mod_name in mods) {
  mg <- net_genes[colors == mod_name]
  for (set_name in c("SFARI_all", "SFARI_1_2", "SFARI_1", names(ifn_sets))) {
    target <- switch(set_name, SFARI_all = sfari_all, SFARI_1_2 = sfari_12, SFARI_1 = sfari_1, ifn_sets[[set_name]])
    r <- perm_enrich(mg, target, net_genes, n_perm = 5000)
    enrich_rows[[length(enrich_rows) + 1]] <- data.frame(module = mod_name, n_genes = length(mg), set = set_name,
                                                          observed = r["observed"], expected = round(r["expected"], 2),
                                                          fold = round(r["fold"], 2), P = r["p"])
  }
}
enr <- do.call(rbind, enrich_rows)
enr$FDR <- p.adjust(enr$P, method = "BH")
rownames(enr) <- NULL
write.csv(enr, file.path(res_dir, "07_module_enrichment_sfari_ifn.csv"), row.names = FALSE)
cat("\nModules enriched (permutation FDR<0.1) for SFARI or immune sets:\n")
print(enr[enr$FDR < 0.1, ], row.names = FALSE)
cat("Best raw p per set:\n")
for (s in unique(enr$set)) {
  r <- enr[enr$set == s, ]
  r <- r[which.min(r$P), ]
  cat("  ", s, ": module", r$module, "observed", r$observed, "vs expected", r$expected, "(p =", signif(r$P, 3), ", FDR =", signif(r$FDR, 3), ")\n")
}

# ===============================================================
# 6. Replication: module preservation + eigengene test in GSE102741
# ===============================================================
cat("\n=== Module replication in GSE102741 ===\n")
x2 <- as.matrix(read.csv(file.path(clean_dir, "GSE102741_expr.csv"), row.names = 1, check.names = FALSE))
m2 <- read.csv(file.path(clean_dir, "GSE102741_meta.csv"), stringsAsFactors = FALSE)
m2$group <- factor(m2$group, levels = c("control", "autism"))
shared <- intersect(net_genes, rownames(x2))
cat("network genes also measured in GSE102741:", length(shared), "of", length(net_genes), "\n")
# adjust GSE102741 for its known covariates, keep diagnosis
x2_adj <- removeBatchEffect(x2[shared, ], covariates = cbind(age = m2$age, rin = m2$rin, male = as.numeric(m2$sex == "Male")),
                            batch = m2$race, design = model.matrix(~ group, data = m2))
datExpr2 <- t(x2_adj)
colors_sh <- colors[match(shared, net_genes)]

# eigengenes of the discovery modules, computed inside GSE102741
colors_sh_noGrey <- colors_sh
me2 <- moduleEigengenes(datExpr2, colors = colors_sh_noGrey, excludeGrey = TRUE)$eigengenes
colnames(me2) <- sub("^ME", "", colnames(me2))
fit2 <- lmFit(t(me2), model.matrix(~ group, data = m2))
fit2 <- eBayes(fit2)
rep_tab <- data.frame(module = rownames(fit2$coefficients), diff_ME_102741 = fit2$coefficients[, "groupautism"],
                      t_102741 = fit2$t[, "groupautism"], P_102741 = fit2$p.value[, "groupautism"])
rep_tab <- merge(me_tab, rep_tab, by = "module")
rep_tab$same_direction <- sign(rep_tab$diff_ME) == sign(rep_tab$diff_ME_102741)
rep_tab <- rep_tab[order(rep_tab$P), ]
cat("\nDiscovery vs replication module-level diagnosis effect:\n")
print(rep_tab[, c("module", "n_genes", "diff_ME", "P", "FDR", "diff_ME_102741", "P_102741", "same_direction")], row.names = FALSE)
cat("module-level effect correlation across modules (GSE28521 t vs GSE102741 t): r =",
    round(cor(rep_tab$t, rep_tab$t_102741), 2), "\n")

# preservation of module structure (density / connectivity), 100 permutations
multiExpr <- list(disc = list(data = datExpr[, shared]), rep = list(data = datExpr2))
multiColor <- list(disc = colors_sh)
mp <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1, nPermutations = 100,
                         networkType = "signed", randomSeed = 1, quickCor = 0, verbose = 0, maxGoldModuleSize = 100, maxModuleSize = 400)
pres <- mp$preservation$Z$ref.disc$inColumnsAlsoPresentIn.rep
pres_tab <- data.frame(module = rownames(pres), size = pres$moduleSize, Zsummary = round(pres$Zsummary.pres, 1))
pres_tab <- pres_tab[!pres_tab$module %in% c("gold"), ]
pres_tab$verdict <- ifelse(pres_tab$Zsummary > 10, "strong", ifelse(pres_tab$Zsummary > 2, "weak/moderate", "not preserved"))
cat("\nModule preservation Zsummary (>10 strong, 2-10 weak/moderate, <2 none):\n")
print(pres_tab[order(-pres_tab$Zsummary), ], row.names = FALSE)

final <- merge(rep_tab, pres_tab[, c("module", "Zsummary", "verdict")], by = "module", all.x = TRUE)
final <- final[order(final$P), ]
write.csv(final, file.path(res_dir, "07_module_trait_and_replication.csv"), row.names = FALSE)

# eigengene heatmap-like plot of the diagnosis effect
png(file.path(fig_dir, "07_module_diagnosis_effects.png"), width = 1000, height = 600)
barplot(final$t[order(final$t)], names.arg = final$module[order(final$t)], col = final$module[order(final$t)],
        las = 2, ylab = "t-statistic (autism vs control)", main = "Module eigengene differences, GSE28521 cortex")
abline(h = c(-2, 2), lty = 2)
dev.off()

cat("\n=========== PLAIN-ENGLISH SUMMARY ===========\n")
cat("Modules:", length(mods), "| modules with FDR<0.05 for diagnosis:", sum(me_tab$FDR < 0.05), "\n")
cat("Look for modules that (1) differ by diagnosis, (2) point the same way in GSE102741, (3) have Zsummary>2,\n")
cat("and (4) are enriched for SFARI or interferon genes. Only modules passing (1)-(3) count as replicated.\n")
