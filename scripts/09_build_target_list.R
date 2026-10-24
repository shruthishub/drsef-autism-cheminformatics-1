# 09_build_target_list.R
# Step 1 of the target -> drug plan: which genes are UP in autism brain and agree across cohorts?
# Uses only existing result tables (no new data). Independent cohorts = GSE28521 cortex,
# GSE102741, Velmeshev pseudobulk. GSE64018 shares 11 donors with GSE28521 so it is shown
# but NOT counted as independent.
#
# Run from the project folder:  Rscript scripts/09_build_target_list.R

res_dir <- "results"
rd <- function(f) read.csv(file.path(res_dir, f), check.names = FALSE)

a <- rd("04_deg_GSE28521_cortex.csv");  a$z <- a$t
b <- rd("04_deg_GSE102741.csv");        b$z <- b$t
c <- rd("04_deg_GSE64018.csv");         c$z <- c$t
v <- rd("08_deg_velmeshev_total.csv")   # already has z

keep <- function(d, tag) {
  d <- d[!duplicated(d$gene), c("gene", "logFC", "P.Value", "z")]
  names(d)[-1] <- paste0(names(d)[-1], "_", tag)
  d
}
m <- Reduce(function(x, y) merge(x, y, by = "gene", all = FALSE),
            list(keep(a, "28521"), keep(b, "102741"), keep(v, "velm")))
m <- merge(m, keep(c, "64018"), by = "gene", all.x = TRUE)   # sensitivity only
cat("Genes measured in all 3 independent cohorts:", nrow(m), "\n")

# Equal-weight Stouffer across the 3 independent cohorts (sample sizes differ a lot;
# equal weights stop the biggest cohort dominating).
m$z_meta <- (m$z_28521 + m$z_102741 + m$z_velm) / sqrt(3)
m$P_meta <- 2 * pnorm(-abs(m$z_meta))
m$FDR_meta <- p.adjust(m$P_meta, "BH")
m$n_up <- (m$logFC_28521 > 0) + (m$logFC_102741 > 0) + (m$logFC_velm > 0)
m$n_p05_up <- (m$logFC_28521 > 0 & m$P.Value_28521 < 0.05) +
              (m$logFC_102741 > 0 & m$P.Value_102741 < 0.05) +
              (m$logFC_velm > 0 & m$P.Value_velm < 0.05)

# Known vs novel: SFARI only for now. Literature cross-check is a separate step.
sfari <- read.csv("data_raw/reference/sfari_genes.csv", check.names = FALSE)
m$sfari_score <- sfari[["gene-score"]][match(m$gene, sfari[["gene-symbol"]])]
m$known_sfari <- !is.na(m$sfari_score)

up <- m[m$z_meta > 0, ]
cat("\nUP candidates under different rules (z_meta > 0):\n")
rules <- list(
  "up in all 3 cohorts"                       = up$n_up == 3,
  "up in all 3, FDR_meta < 0.10"              = up$n_up == 3 & up$FDR_meta < 0.10,
  "up in all 3, FDR_meta < 0.05"              = up$n_up == 3 & up$FDR_meta < 0.05,
  "up in all 3, P_meta < 0.01"                = up$n_up == 3 & up$P_meta < 0.01,
  "up in all 3, nominal p<0.05 in >= 2"       = up$n_up == 3 & up$n_p05_up >= 2,
  "up in >= 2, nominal p<0.05 in >= 2"        = up$n_up >= 2 & up$n_p05_up >= 2)
for (nm in names(rules)) {
  k <- rules[[nm]]
  cat(sprintf("  %-40s %5d genes (%d SFARI)\n", nm, sum(k), sum(up$known_sfari[k])))
}

# Expected by chance: if no real signal, P(up in all 3) = 1/8 of genes.
cat(sprintf("\nChance check: %d of %d genes are up in all 3 cohorts (%.1f%%); pure noise gives 12.5%%.\n",
            sum(m$n_up == 3), nrow(m), 100 * mean(m$n_up == 3)))
# Permutation of direction would be stronger; this is a quick screen only.

cand <- up[up$n_up == 3 & up$P_meta < 0.01, ]
cand <- cand[order(cand$P_meta), ]
cand$tier <- ifelse(cand$known_sfari, "known(SFARI)", "novel(vs SFARI only)")
write.csv(m,    file.path(res_dir, "09_cross_cohort_all_genes.csv"), row.names = FALSE)
write.csv(cand, file.path(res_dir, "09_up_candidates.csv"), row.names = FALSE)
cat("\nTop candidates:\n")
print(head(cand[, c("gene","logFC_28521","logFC_102741","logFC_velm","z_meta","P_meta","FDR_meta","sfari_score","tier")], 25), row.names = FALSE)
