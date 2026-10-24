# 19_export_genesets_and_oligo.R
# Exports (a) GO:BP + GO:CC gene sets (msigdbr, human) as long csv and (b) Velmeshev oligodendrocyte pseudobulk (log2 CPM) + meta,
# so the Python pathway test (20_myelin_pathway.py) needs no R. Run: Rscript scripts/19_export_genesets_and_oligo.R
suppressMessages({library(msigdbr); library(Matrix)})
dir.create("data_clean", showWarnings = FALSE)
gs <- rbind(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
            msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC"))
write.csv(unique(gs[, c("gs_name", "gene_symbol", "ncbi_gene")]), "data_clean/go_bp_cc_sets.csv", row.names = FALSE)
cnt <- readRDS("data_clean/velmeshev_pb_counts.rds"); meta <- read.csv("data_clean/velmeshev_pb_meta.csv")
meta <- meta[match(colnames(cnt), meta$group_id), ]
ol <- which(meta$celltype == "Oligodendrocytes" & meta$n_cells >= 10)
x <- as.matrix(cnt[, ol]); lib <- colSums(x); lc <- log2(t(t(x) / lib * 1e6) + 1)
lc <- lc[rowSums(x >= 1) >= 10, ]
write.csv(lc, "data_clean/velmeshev_oligo_logcpm.csv"); write.csv(meta[ol, ], "data_clean/velmeshev_oligo_meta.csv", row.names = FALSE)
cat("oligodendrocyte pseudobulks:", length(ol), "| genes kept:", nrow(lc), "| by diagnosis:", paste(names(table(meta$diagnosis[ol])), table(meta$diagnosis[ol]), collapse = " "), "\n")
