# 11_druggability.R
# Step 2b: which candidate genes can be drugged / docked? Free public APIs, no login.
#  UniProt  -> reviewed human accession, protein name, number of PDB structures, AlphaFold model
#  ChEMBL   -> target, drugs with a recorded mechanism (max_phase), number of potent compounds (pChEMBL >= 6)
# Records the query date because both databases change over time.
# Run from the project folder:  Rscript scripts/11_druggability.R
suppressMessages(library(jsonlite))
res_dir <- "results"
cand <- read.csv(file.path(res_dir, "10_candidates_with_support.csv"), check.names = FALSE)
genes <- cand$gene

get_txt <- function(url) { Sys.sleep(0.15); tryCatch(paste(system2("curl", c("-s", "-m", "30", shQuote(url)), stdout = TRUE), collapse = "\n"), error = function(e) "") }
get_json <- function(url) { x <- get_txt(url); if (!nzchar(x)) return(NULL); tryCatch(fromJSON(x, simplifyVector = FALSE), error = function(e) NULL) }

out <- lapply(genes, function(g) {
  u <- read.delim(text = get_txt(sprintf(
    "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:%s+AND+organism_id:9606+AND+reviewed:true&fields=accession,protein_name,xref_pdb,xref_alphafolddb&format=tsv", g)),
    sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(u) == 0) return(data.frame(gene = g, uniprot = NA, protein = NA, n_pdb = NA, alphafold = NA,
                                      chembl_target = NA, n_mech = NA, max_phase = NA, mech_drugs = NA, n_potent = NA))
  acc <- u[[1]][1]
  pdb <- as.character(u[["PDB"]][1]); af <- as.character(u[["AlphaFoldDB"]][1])
  pdb <- if (is.na(pdb)) "" else pdb; af <- if (is.na(af)) "" else af
  n_pdb <- sum(nzchar(strsplit(pdb, ";")[[1]]))
  tg <- get_json(sprintf("https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=%s&target_type=SINGLE%%20PROTEIN&limit=1", acc))
  tid <- if (!is.null(tg) && length(tg$targets)) tg$targets[[1]]$target_chembl_id else NA
  n_mech <- NA; maxp <- NA; drugs <- NA; n_pot <- NA
  if (!is.na(tid)) {
    me <- get_json(sprintf("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id=%s&limit=100", tid))
    if (!is.null(me)) {
      n_mech <- length(me$mechanisms)
      ph <- sapply(me$mechanisms, function(x) ifelse(is.null(x$max_phase), NA, as.numeric(x$max_phase)))
      maxp <- if (length(ph)) suppressWarnings(max(ph, na.rm = TRUE)) else NA
      if (is.infinite(maxp)) maxp <- NA
      drugs <- paste(unique(na.omit(sapply(me$mechanisms, function(x) x$molecule_chembl_id))), collapse = ";")
    }
    ac <- get_json(sprintf("https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=%s&pchembl_value__gte=6&limit=1", tid))
    n_pot <- if (!is.null(ac)) ac$page_meta$total_count else NA
  }
  data.frame(gene = g, uniprot = acc, protein = u[[2]][1], n_pdb = n_pdb, alphafold = nzchar(af),
             chembl_target = tid, n_mech = n_mech, max_phase = maxp, mech_drugs = drugs, n_potent = n_pot)
})
d <- do.call(rbind, out)
d$query_date <- as.character(Sys.Date())
all <- merge(cand, d, by = "gene", sort = FALSE)
write.csv(all, file.path(res_dir, "11_candidates_druggability.csv"), row.names = FALSE)

cat("Queried", nrow(d), "genes on", as.character(Sys.Date()), "\n")
cat("  with reviewed UniProt entry:", sum(!is.na(d$uniprot)), "| with PDB structure:", sum(d$n_pdb > 0, na.rm = TRUE), "| with AlphaFold:", sum(d$alphafold, na.rm = TRUE), "\n")
cat("  with ChEMBL single-protein target:", sum(!is.na(d$chembl_target)), "| with drug/clinical mechanism:", sum(d$n_mech > 0, na.rm = TRUE),
    "| >=100 potent compounds:", sum(d$n_potent >= 100, na.rm = TRUE), "\n")
