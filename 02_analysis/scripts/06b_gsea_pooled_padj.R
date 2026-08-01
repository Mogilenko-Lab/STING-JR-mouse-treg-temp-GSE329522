# 06b_gsea_pooled_padj.R — COMPUTE
## Compute pooled Benjamini-Hochberg FDR correction across all GSEA database families
## for focal contrasts (WT_heat and Interaction).
##
## Reads:
##   03_results/06_gsea/tables/by_contrast/<contrast>/gsea_msigdb.csv
##   03_results/06_gsea/tables/by_contrast/<contrast>/gsea_custom.csv
##
## Writes:
##   03_results/06_gsea/tables/_overview/gsea_pooled_overview_WT_heat.csv
##   03_results/06_gsea/tables/_overview/gsea_pooled_overview_Interaction.csv
##   03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv
##   03_results/06_gsea/tables/_overview/gsea_pooled_summary_by_db.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# 1. Setup + Paths
STAGE <- "06_gsea"
FOCAL_CONTRASTS <- c("WT_heat", "Interaction")

base_dir <- file.path("03_results", STAGE, "tables")
by_contrast_dir <- file.path(base_dir, "by_contrast")
overview_dir <- file.path(base_dir, "_overview")
dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)

all_pooled_rows <- list()
all_summary_rows <- list()

for (co in FOCAL_CONTRASTS) {
  fp_m <- file.path(by_contrast_dir, co, "gsea_msigdb.csv")
  fp_c <- file.path(by_contrast_dir, co, "gsea_custom.csv")

  if (!file.exists(fp_m) || !file.exists(fp_c)) {
    stop(sprintf("Source GSEA tables missing for contrast %s at %s or %s", co, fp_m, fp_c))
  }

  df_m <- readr::read_csv(fp_m, show_col_types = FALSE)
  df_c <- readr::read_csv(fp_c, show_col_types = FALSE)
  df_combined <- dplyr::bind_rows(df_m, df_c)

  # Check tie floor
  min_p <- min(df_combined$pvalue, na.rm = TRUE)
  n_min_p <- sum(df_combined$pvalue == min_p, na.rm = TRUE)

  # Calculate n_tests_in_db and published sig per database
  db_counts <- df_combined %>%
    dplyr::group_by(database) %>%
    dplyr::summarise(
      n_tests_in_db = dplyr::n(),
      sig_before_in_db = sum(padj < 0.05, na.rm = TRUE),
      .groups = "drop"
    )

  # Pool pvalues across all databases for this contrast and adjust with BH
  n_pooled <- nrow(df_combined)
  df_pooled <- df_combined %>%
    dplyr::mutate(
      padj_pooled = stats::p.adjust(pvalue, method = "BH"),
      n_tests_pooled = n_pooled
    ) %>%
    dplyr::left_join(db_counts, by = "database")

  # Add per-database post-pooling count
  db_post <- df_pooled %>%
    dplyr::group_by(database) %>%
    dplyr::summarise(
      sig_after_in_db = sum(padj_pooled < 0.05, na.rm = TRUE),
      .groups = "drop"
    )

  db_summary <- db_counts %>%
    dplyr::left_join(db_post, by = "database") %>%
    dplyr::mutate(
      sig_lost_in_db = sig_before_in_db - sig_after_in_db
    )

  df_pooled <- df_pooled %>%
    dplyr::left_join(db_summary %>% dplyr::select(database, sig_after_in_db, sig_lost_in_db), by = "database")

  # Build per-database summary dataframe
  summary_by_db <- db_summary %>%
    dplyr::mutate(
      contrast = co,
      n_tests_pooled = n_pooled,
      pct_sig_before = round(100 * sig_before_in_db / n_tests_in_db, 1),
      pct_sig_after = round(100 * sig_after_in_db / n_tests_in_db, 1),
      min_pvalue = min_p,
      n_at_min_pvalue = n_min_p
    ) %>%
    dplyr::select(contrast, database, n_tests_in_db, n_tests_pooled, sig_before_in_db, sig_after_in_db, sig_lost_in_db, pct_sig_before, pct_sig_after, min_pvalue, n_at_min_pvalue)

  # Reorder columns according to spec
  req_cols <- c("contrast", "database", "pathway_id", "pathway_name", "direction",
                "nes", "pvalue", "padj", "padj_pooled", "set_size",
                "n_tests_in_db", "n_tests_pooled", "sig_before_in_db", "sig_after_in_db", "sig_lost_in_db")

  df_out <- df_pooled %>% dplyr::select(dplyr::all_of(req_cols), dplyr::everything())

  # Save contrast-specific table
  out_co_fp <- file.path(overview_dir, sprintf("gsea_pooled_overview_%s.csv", co))
  readr::write_csv(df_out, out_co_fp)

  all_pooled_rows[[co]] <- df_out
  all_summary_rows[[co]] <- summary_by_db

  cat(sprintf("[06b] %s: %d tests pooled across %d databases. Min pvalue: %g (n=%d).\n",
              co, n_pooled, dplyr::n_distinct(df_out$database), min_p, n_min_p))
  cat(sprintf("      Sig before: %d, Sig after: %d, Net change: %d\n",
              sum(summary_by_db$sig_before_in_db), sum(summary_by_db$sig_after_in_db),
              sum(summary_by_db$sig_before_in_db) - sum(summary_by_db$sig_after_in_db)))
}

# Combine and write master pooled table and summary table
master_pooled <- dplyr::bind_rows(all_pooled_rows)
master_summary <- dplyr::bind_rows(all_summary_rows)

readr::write_csv(master_pooled, file.path(overview_dir, "gsea_pooled_overview.csv"))
readr::write_csv(master_summary, file.path(overview_dir, "gsea_pooled_summary_by_db.csv"))

cat(sprintf("[06b] Master pooled GSEA overview tables written to %s/\n", overview_dir))
