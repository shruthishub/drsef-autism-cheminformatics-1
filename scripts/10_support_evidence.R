# 10_support_evidence.R
# Step 2a: extra evidence for the UP candidates from script 09.
#  (1) Cerebellum direction: limma on GSE28521 cerebellum (21 samples, 1 per donor) + 2 surrogate variables.
#      GSE28521 is the SAME cohort as the cortex data, so this is NOT independent replication; it only shows
#      whether the gene also points up in another region of the same brains.
#  (2) WGCNA module membership (GSE28521 cortex) and the module's diagnosis effect.
# Run from the project folder:  Rscript scripts/10_support_evidence.R
suppressMessages({library(limma); library(sva)})
res_dir <- "results"

cand <- read.csv(file.path(res_dir, "09_up_candidates.csv"), check.names = FALSE)

# ---- (1) cerebellum DE ----
meta <- read.csv("data_clean/GSE28521_cerebellum_meta.csv", check.names = FALSE)
ex   <- read.csv("data_clean/GSE28521_cerebellum_expr.csv", check.names = FALSE, row.names = 1)
ex   <- as.matrix(ex[, meta$sample])
stopifnot(!anyDuplicated(meta$donor))
cat("Cerebellum:", table(meta$group), "\n")
grp  <- factor(meta$group, levels = c("control", "autism"))
mod  <- model.matrix(~ grp); mod0 <- model.matrix(~ 1, data = meta)
sv   <- sva(ex, mod, mod0, n.sv = 2)$sv
design <- cbind(mod, sv)
fit  <- eBayes(lmFit(ex, design))
tt   <- topTable(fit, coef = "grpautism", number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
lam <- median(tt$t^2) / qchisq(0.5, 1)
cat("Cerebellum lambda:", round(lam, 2), "| FDR<0.05:", sum(tt$adj.P.Val < 0.05), "\n")
write.csv(tt[, c("gene","logFC","t","P.Value","adj.P.Val")], file.path(res_dir, "10_deg_GSE28521_cerebellum.csv"), row.names = FALSE)

# ---- (2) WGCNA modules ----
mods <- read.csv(file.path(res_dir, "07_module_genes.csv"), check.names = FALSE)
mt   <- read.csv(file.path(res_dir, "07_module_trait_and_replication.csv"), check.names = FALSE)

cand$cereb_logFC <- tt$logFC[match(cand$gene, tt$gene)]
cand$cereb_P     <- tt$P.Value[match(cand$gene, tt$gene)]
cand$cereb_up    <- cand$cereb_logFC > 0
cand$module      <- mods$module[match(cand$gene, mods$gene)]
cand$module_diffME <- mt$diff_ME[match(cand$module, mt$module)]
cand$module_FDR    <- mt$FDR[match(cand$module, mt$module)]
cand$module_up_in_asd <- !is.na(cand$module_diffME) & cand$module_diffME > 0 & cand$module_FDR < 0.05

cat("\nOf", nrow(cand), "candidates (up in all 3 cohorts, P_meta<0.01):\n")
cat("  also up in cerebellum:", sum(cand$cereb_up, na.rm = TRUE), "| and cereb P<0.05:", sum(cand$cereb_up & cand$cereb_P < 0.05, na.rm = TRUE), "\n")
cat("  in a WGCNA module that is UP in autism (FDR<0.05):", sum(cand$module_up_in_asd), "\n")
print(table(cand$module, useNA = "ifany"))
cat("Chance for cerebellum 'up': 50% of ", sum(!is.na(cand$cereb_up)), " = ", round(0.5 * sum(!is.na(cand$cereb_up))), "\n")
write.csv(cand, file.path(res_dir, "10_candidates_with_support.csv"), row.names = FALSE)
