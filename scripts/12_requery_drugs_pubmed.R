# 12_requery_drugs_pubmed.R
# Step 3 of target plan:
#  (1) re-query candidates that had no ChEMBL result in script 11 (some were API timeouts), with retries
#  (2) drug names, mechanism, max phase and indications for every candidate target that has a mechanism
#  (3) PubMed novelty check: how many papers mention the gene together with autism
# Run from the project folder:  Rscript scripts/12_requery_drugs_pubmed.R
suppressMessages(library(jsonlite))
res_dir <- "results"
d <- read.csv(file.path(res_dir, "11_candidates_druggability.csv"), check.names = FALSE)

fetch <- function(url, tries = 4, wait = 0.2) {
  for (i in seq_len(tries)) {
    Sys.sleep(wait * i)
    x <- suppressWarnings(system2("curl", c("-s", "-m", "60", shQuote(url)), stdout = TRUE, stderr = FALSE))
    x <- paste(x, collapse = "\n")
    j <- tryCatch(fromJSON(x, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(j)) return(j)
  }
  NULL
}

# ---- (1) re-query missing ChEMBL targets ----
redo <- which(is.na(d$chembl_target) & !is.na(d$uniprot))
cat("Re-querying", length(redo), "genes without a ChEMBL target...\n")
n_fixed <- 0; failed <- character()
for (i in redo) {
  tg <- fetch(sprintf("https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=%s&target_type=SINGLE%%20PROTEIN&limit=1", d$uniprot[i]))
  if (is.null(tg)) { failed <- c(failed, d$gene[i]); next }       # still a failure, NOT a true absence
  if (length(tg$targets)) {
    tid <- tg$targets[[1]]$target_chembl_id; d$chembl_target[i] <- tid; n_fixed <- n_fixed + 1
    me <- fetch(sprintf("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id=%s&limit=100", tid))
    ac <- fetch(sprintf("https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=%s&pchembl_value__gte=6&limit=1", tid))
    if (!is.null(me)) {
      d$n_mech[i] <- length(me$mechanisms)
      ph <- suppressWarnings(as.numeric(unlist(lapply(me$mechanisms, function(x) x$max_phase))))
      d$max_phase[i] <- if (length(ph) && any(!is.na(ph))) max(ph, na.rm = TRUE) else NA
      d$mech_drugs[i] <- paste(unique(unlist(lapply(me$mechanisms, function(x) x$molecule_chembl_id))), collapse = ";")
    }
    if (!is.null(ac)) d$n_potent[i] <- ac$page_meta$total_count
  }
}
cat("  now found:", n_fixed, "| still failed after retries:", length(failed), if (length(failed)) paste0("(", paste(failed, collapse = ","), ")"), "\n")
cat("  total with ChEMBL target:", sum(!is.na(d$chembl_target)), "| with mechanism:", sum(d$n_mech > 0, na.rm = TRUE), "\n")

# ---- (2) drug names + indications ----
brain_re <- "schizo|bipolar|mania|depress|epilep|seizure|parkinson|alzheimer|dementia|autis|down syndrome|anxiety|attention deficit|hyperactiv|huntington|psychosis|migraine|neuro"
rows <- list()
for (i in which(d$n_mech > 0)) {
  me <- fetch(sprintf("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id=%s&limit=100", d$chembl_target[i]))
  for (m in me$mechanisms) {
    mol <- m$molecule_chembl_id; if (is.null(mol)) next
    rows[[length(rows) + 1]] <- data.frame(gene = d$gene[i], molecule = mol,
      action = ifelse(is.null(m$action_type), NA, m$action_type),
      mechanism = ifelse(is.null(m$mechanism_of_action), NA, m$mechanism_of_action), stringsAsFactors = FALSE)
  }
}
dr <- unique(do.call(rbind, rows))
cat("Mechanism records:", nrow(dr), "for", length(unique(dr$gene)), "genes,", length(unique(dr$molecule)), "molecules\n")
info <- lapply(unique(dr$molecule), function(mol) {
  mj <- fetch(sprintf("https://www.ebi.ac.uk/chembl/api/data/molecule/%s.json", mol))
  ind <- fetch(sprintf("https://www.ebi.ac.uk/chembl/api/data/drug_indication.json?molecule_chembl_id=%s&limit=200", mol))
  hd <- if (!is.null(ind)) unique(na.omit(unlist(lapply(ind$drug_indications, function(x) x$mesh_heading)))) else character()
  data.frame(molecule = mol, drug_name = ifelse(is.null(mj$pref_name), NA, mj$pref_name),
             drug_max_phase = ifelse(is.null(mj$max_phase), NA, mj$max_phase),
             indications = paste(head(hd, 12), collapse = "; "),
             brain_indication = any(grepl(brain_re, hd, ignore.case = TRUE)), stringsAsFactors = FALSE)
})
info <- do.call(rbind, info)
dr <- merge(dr, info, by = "molecule", all.x = TRUE)
dr <- dr[order(dr$gene, -as.numeric(dr$drug_max_phase)), ]
write.csv(dr, file.path(res_dir, "12_target_drug_pairs.csv"), row.names = FALSE)

# ---- (3) PubMed novelty ----
pm <- function(term) {
  x <- fetch(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=0&term=%s", URLencode(term, reserved = TRUE)), wait = 0.35)
  if (is.null(x)) NA else as.integer(x$esearchresult$count)
}
d$pubmed_asd <- sapply(d$gene, function(g) pm(sprintf('%s[tiab] AND (autism[tiab] OR "autism spectrum disorder"[tiab] OR ASD[tiab])', g)))
d$pubmed_total <- sapply(d$gene, function(g) pm(sprintf("%s[tiab]", g)))
d$pubmed_date <- as.character(Sys.Date())
write.csv(d, file.path(res_dir, "12_candidates_full.csv"), row.names = FALSE)
cat("PubMed: autism papers per gene, median", median(d$pubmed_asd, na.rm = TRUE), "| 0 papers:", sum(d$pubmed_asd == 0, na.rm = TRUE),
    "| 1-4:", sum(d$pubmed_asd %in% 1:4), "| >=5:", sum(d$pubmed_asd >= 5, na.rm = TRUE), "| failed:", sum(is.na(d$pubmed_asd)), "\n")
