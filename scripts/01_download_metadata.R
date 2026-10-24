# 01_download_metadata.R
# Step: download the sample tables (metadata) for the 4 GEO datasets and count
# autism vs control per dataset and per brain region.
#
# Why GEOquery: it is the standard R package for reading GEO "series matrix"
# files into R objects, so we get the sample labels without scraping web pages.
#
# Run from the project folder:  Rscript scripts/01_download_metadata.R

library(GEOquery)

# ---- settings ----
accessions <- c("GSE59288", "GSE28521", "GSE28475", "GSE102741")
raw_dir <- "data_raw"
res_dir <- "results"
options(timeout = 600)  # downloads can be slow

# ---- helper: join all "characteristics" columns of a sample into one text ----
# GEO stores things like "diagnosis: autism" or "brain region: frontal cortex"
# in several columns. Pasting them together lets us search them with a regex.
get_all_text <- function(pheno) {
  # pick columns that hold sample descriptions
  keep <- grep("characteristics|title|source_name|description", colnames(pheno))
  all_text <- character(nrow(pheno))
  for (i in 1:nrow(pheno)) {
    pieces <- c()
    for (k in keep) {
      pieces <- c(pieces, as.character(pheno[i, k]))
    }
    all_text[i] <- paste(pieces, collapse = " | ")
  }
  return(all_text)
}

# ---- helper: guess diagnosis from the text (CHECK the printed table by eye!) ----
guess_diagnosis <- function(txt) {
  txt <- tolower(txt)
  out <- rep("unknown", length(txt))
  for (i in 1:length(txt)) {
    if (grepl("control|ctl|normal|healthy", txt[i])) {
      out[i] <- "control"
    } else if (grepl("autis|asd", txt[i])) {
      out[i] <- "autism"
    }
  }
  return(out)
}

# ---- helper: guess brain region ----
guess_region <- function(txt) {
  txt <- tolower(txt)
  out <- rep("unknown", length(txt))
  for (i in 1:length(txt)) {
    if (grepl("cerebell", txt[i])) {
      out[i] <- "cerebellum"
    } else if (grepl("temporal|BA41|BA21|BA22|ba41|ba21|ba22", txt[i])) {
      out[i] <- "temporal cortex"
    } else if (grepl("frontal|prefrontal|dlpfc|ba9|ba10|ba44", txt[i])) {
      out[i] <- "frontal cortex"
    }
  }
  return(out)
}

summary_rows <- list()

for (acc in accessions) {
  cat("\n=====================================================\n")
  cat("Downloading", acc, "\n")
  cat("=====================================================\n")

  # getGPL = FALSE: we only need sample info now, not the platform annotation
  gse_list <- tryCatch(
    getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE, destdir = raw_dir),
    error = function(e) {
      cat("!! FAILED to download", acc, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(gse_list)) next

  cat("Number of series-matrix objects returned:", length(gse_list), "\n")

  # some GEO series have more than one platform -> more than one object
  for (j in 1:length(gse_list)) {
    gse <- gse_list[[j]]
    pheno <- pData(gse)
    tag <- paste0(acc, "_part", j)
    cat("\n--", tag, ": platform", annotation(gse), ",", nrow(pheno), "samples,",
        nrow(gse) , "features in matrix\n")

    all_text <- get_all_text(pheno)
    pheno$guess_diagnosis <- guess_diagnosis(all_text)
    pheno$guess_region <- guess_region(all_text)

    # save the full table, never overwrite raw downloads
    write.csv(pheno, file.path(res_dir, paste0(tag, "_sample_table.csv")),
              row.names = TRUE)

    # show what the characteristic columns look like, so I can check the guesses
    cat("\nColumns available:\n")
    print(colnames(pheno))
    char_cols <- grep("characteristics", colnames(pheno), value = TRUE)
    for (cc in char_cols) {
      cat("\nValues in", cc, "(first 8 unique):\n")
      print(head(unique(as.character(pheno[[cc]])), 8))
    }

    cat("\nGuessed diagnosis counts:\n")
    print(table(pheno$guess_diagnosis))
    cat("\nGuessed diagnosis x region:\n")
    print(table(pheno$guess_diagnosis, pheno$guess_region))

    # one summary row per series
    n_aut <- sum(pheno$guess_diagnosis == "autism")
    n_ctl <- sum(pheno$guess_diagnosis == "control")
    n_unk <- sum(pheno$guess_diagnosis == "unknown")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      series = tag, platform = annotation(gse), n_samples = nrow(pheno),
      n_autism = n_aut, n_control = n_ctl, n_unknown = n_unk
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(res_dir, "01_sample_counts_summary.csv"), row.names = FALSE)

cat("\n\n================ SUMMARY (guessed labels, VERIFY) ================\n")
print(summary_df)
cat("\nPlan numbers to compare against (from the paper, NOT verified):\n")
cat(" GSE59288: 34 autism / 40 control (GEO lists 73 samples)\n")
cat(" GSE28521: 58 cortex samples (29/29) + 21 cerebellum\n")
cat(" GSE28475: ~20 autism cortex (or up to 52 arrays)\n")
cat(" GSE102741: 13 autism / 39 control\n")
cat("\nWHAT TO CHECK: open results/*_sample_table.csv and confirm that\n")
cat("'unknown' is 0 and that the guessed labels match the real GEO labels.\n")
