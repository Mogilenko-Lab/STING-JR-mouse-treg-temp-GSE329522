#!/usr/bin/env Rscript
# 00_setup_metadata.R - Build the sample-metadata object for GSE329522
# Project: STING-cGAS-GSE329522 (2x2 genotype x temperature iTreg bulk RNA-seq)
# Inputs:
#   - 02_analysis/config/analysis_config.yaml (via config.R)
#   - GSE329522 deposited CPM column order (12630-RS-021 .. 040) [knowledge, not read here]
#   - Owner sample sheet (2026-07-22) -- see 00_data/processed/PROVENANCE.md
# Outputs:
#   - 03_results/objects/sample_metadata.rds
#   - 03_results/master/sample_metadata.csv
# Dependencies: yaml (via config.R). No heavy packages required.
#
# ----------------------------------------------------------------------------
# TWO ORDERINGS (the historical mapping risk, now resolved):
#
#  (A) CPM column order = TEMPERATURE-MAJOR (this is what the analysis uses).
#      OWNER-CONFIRMED by the 2026-07-22 sample sheet; the Hspa1b/Hsph1 heat-shock
#      thermometer + Cgas (WT>KO) inference matched it exactly, 20/20:
#        021-025 = WT_37     026-030 = cGASKO_37
#        031-035 = WT_39     036-040 = cGASKO_39
#
#  (B) GEO GSM accession order = GENOTYPE-MAJOR (stored only to document the
#      accession-vs-column discrepancy):
#        GSM9705690..694 = WT 37      GSM9705695..699 = WT 39
#        GSM9705700..704 = cGAS KO 37 GSM9705705..709 = cGAS KO 39
#
#  TRUE GSM per library = match each library's confirmed genotype+temp to its GSM:
#        WT_37     021-025 -> GSM690-694
#        cGASKO_37 026-030 -> GSM700-704
#        WT_39     031-035 -> GSM695-699
#        cGASKO_39 036-040 -> GSM705-709
#
#  gsm_id_positional = the GSM accession list (690..709) lined up POSITIONALLY
#  to columns 021..040 -- what a naive positional join would (wrongly) produce.
#  Because (A) and (B) differ, that naive join mislabels libraries 026-035.
# ----------------------------------------------------------------------------

source("02_analysis/config/config.R")

# ---- Library IDs in deposited CPM column order (temperature-major) ----------
library_id <- sprintf("12630-RS-%03d", 21:40)

# ---- Owner-confirmed group per library (TEMPERATURE-MAJOR; analysis truth) ---
group <- rep(c("WT_37", "cGASKO_37", "WT_39", "cGASKO_39"), each = 5)

genotype <- ifelse(grepl("^WT_", group), "WT", "cGASKO")
temp     <- sub(".*_", "", group)

# ---- TRUE GSM per library (match inferred genotype+temp to GSM) -------------
# WT_37 -> 690-694, cGASKO_37 -> 700-704, WT_39 -> 695-699, cGASKO_39 -> 705-709
gsm_true_by_group <- list(
  WT_37     = sprintf("GSM%d", 9705690:9705694),
  cGASKO_37 = sprintf("GSM%d", 9705700:9705704),
  WT_39     = sprintf("GSM%d", 9705695:9705699),
  cGASKO_39 = sprintf("GSM%d", 9705705:9705709)
)
gsm_id <- unlist(lapply(c("WT_37", "cGASKO_37", "WT_39", "cGASKO_39"),
                        function(g) gsm_true_by_group[[g]]))

# ---- Naive positional GSM (accession order lined up to columns 021..040) ----
# This is the WRONG join, kept only to document the naive-positional-join trap.
gsm_id_positional <- sprintf("GSM%d", 9705690:9705709)

# ---- Assemble data.frame with proper factor reference levels ---------------
metadata <- data.frame(
  library_id        = library_id,
  genotype          = factor(genotype, levels = c("WT", "cGASKO")),  # ref = WT
  temp              = factor(temp,     levels = c("37", "39")),       # ref = 37
  group             = factor(group, levels = c("WT_37", "cGASKO_37", "WT_39", "cGASKO_39")),
  gsm_id            = gsm_id,
  gsm_id_positional = gsm_id_positional,
  mapping_status    = "CONFIRMED",   # owner sample sheet 2026-07-22 (was INFERRED; matched 20/20)
  stringsAsFactors  = FALSE
)
rownames(metadata) <- metadata$library_id

# ---- Save -------------------------------------------------------------------
dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MASTER,  recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(DIR_OBJECTS, "sample_metadata.rds")
csv_path <- file.path(DIR_MASTER,  "sample_metadata.csv")
saveRDS(metadata, rds_path)
write.csv(metadata, csv_path, row.names = FALSE)

# ---- Summary ----------------------------------------------------------------
cat("\n=== Sample metadata summary ===\n")
cat("Rows:", nrow(metadata), "libraries\n\n")
cat("Group counts (expect 5/5/5/5):\n")
print(table(metadata$group))
cat("\nGenotype reference level:", levels(metadata$genotype)[1], "\n")
cat("Temp reference level:    ", levels(metadata$temp)[1], "\n\n")
cat("Crosstab genotype x temp:\n")
print(table(metadata$genotype, metadata$temp))

stopifnot(all(table(metadata$group) == 5))
stopifnot(levels(metadata$genotype)[1] == "WT")
stopifnot(levels(metadata$temp)[1] == "37")

cat("\n[OK] Wrote:", rds_path, "\n")
cat("[OK] Wrote:", csv_path, "\n")
cat("[OK] 5/5/5/5 balanced; refs = WT & 37; mapping_status = CONFIRMED (owner sample sheet 2026-07-22).\n")
