# 08b_coresh_annotate.R — COMPUTE / DATA (annotation ingest)
## Turn the frozen web-research evidence for the CoReSh-recovered GEO datasets into a
## tracked, tabular annotation table (stage 08_coresh). This is the "what are these
## datasets actually" layer: the mmu compendium ships only gseId + gplId (variance
## structure, no titles), so dataset identity is researched externally and frozen here.
##
## Inputs (frozen, tracked — so the CSV is reproducible from committed files):
##   03_results/08_coresh/tables/_overview/coresh_dataset_annotation.json
##       — structured research evidence, one object per seeded GSE. Produced by the
##         agy CLI (Google Antigravity / Gemini) from a fixed research prompt, then
##         validated and frozen by the orchestrator. The prompt, raw output, and a
##         chain-of-custody log (tool, model, command, timestamp) are kept with the
##         project's reasoning notes; the frozen JSON here is the committed input.
##   03_results/08_coresh/tables/coresh_provenance.csv
##       — the 14 seeded GSEs (cross-check: annotation set must equal seeded set).
##
## Output (compute-only; no plots):
##   03_results/08_coresh/tables/_overview/coresh_dataset_annotation.csv
##       — one row per seeded GSE: identity + IFN/viral/STING/thermal flags + evidence,
##         PMID, source_url, and provenance columns (source_tool/prompt/date).
##
## Run from project root:
##   Rscript 02_analysis/scripts/08b_coresh_annotate.R
##
## NOTE: annotation is descriptive metadata only — it never enters the confirmatory
## spine and cannot change any NES/pctVar. It records what each recovered dataset is, so
## the reader can see what public biology each query signature surfaced. The read is
## exploratory: the context_class and boolean columns describe the datasets themselves.

# ============================================================================
# 0. Environment (config.R FIRST — paths, %||%)
# ============================================================================

source("02_analysis/config/config.R")   # PROJECT_ROOT, YAML_CONFIG, DIR_RESULTS, stage_dir(), %||%

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(readr)
})
options(stringsAsFactors = FALSE)

STAGE       <- "08_coresh"
SCRIPT_PATH <- "02_analysis/scripts/08b_coresh_annotate.R"
tbl_dir     <- stage_dir(STAGE, "tables")
OVERVIEW_DIR <- YAML_CONFIG$figures$overview_dir %||% "_overview"
ov_tbl_dir  <- file.path(tbl_dir, OVERVIEW_DIR)

## Provenance constants (see reasoning-dir RUN.md for the full command + timestamp).
SOURCE_TOOL   <- "agy (Google Antigravity CLI / Gemini, --dangerously-skip-permissions, print-mode)"
SOURCE_PROMPT <- "docs/_internal/reasoning/2026-07-23_coresh_gse_annotation/PROMPT.md"
RETRIEVED_DATE <- "2026-07-23"

json_fp <- file.path(ov_tbl_dir, "coresh_dataset_annotation.json")
prov_fp <- file.path(tbl_dir, "coresh_provenance.csv")
out_fp  <- file.path(ov_tbl_dir, "coresh_dataset_annotation.csv")

# ============================================================================
# 1. GUARD — the frozen research evidence must exist (NEVER fabricate)
# ============================================================================

if (!file.exists(json_fp)) {
  stop(
    "Frozen dataset-annotation evidence not found at ", json_fp, ".\n",
    "  It is produced by external web research (agy) from ", SOURCE_PROMPT,
    " and frozen by the orchestrator; without it there is nothing to ingest."
  )
}

# ============================================================================
# 2. Read + validate the evidence (14 objects, required keys, GSE set matches seed)
# ============================================================================

ann <- jsonlite::fromJSON(json_fp, simplifyDataFrame = TRUE)
if (!is.data.frame(ann)) stop("08b: annotation JSON did not parse to a data.frame.")

req_keys <- c("gse", "recovering_query", "title", "organism", "tissue_cell",
              "perturbation", "context_class", "interferon_related", "viral",
              "cgas_sting", "thermal_heat_shock", "incidental_ifn_note",
              "evidence", "pubmed_id", "source_url", "confidence")
missing_keys <- setdiff(req_keys, colnames(ann))
if (length(missing_keys) > 0L)
  stop("08b: annotation JSON missing required key(s): ", paste(missing_keys, collapse = ", "))

## Cross-check: the annotated GSE set must equal the seeded GSE set (coresh_provenance.csv).
if (file.exists(prov_fp)) {
  seeded <- sort(unique(readr::read_csv(prov_fp, show_col_types = FALSE)$gse))
  annotated <- sort(unique(ann$gse))
  if (!identical(seeded, annotated)) {
    stop("08b: annotated GSE set does not match seeded GSE set.\n",
         "  seeded-only:    ", paste(setdiff(seeded, annotated), collapse = ", "), "\n",
         "  annotated-only: ", paste(setdiff(annotated, seeded), collapse = ", "))
  }
  message(sprintf("[08b] GSE set cross-check OK: %d datasets match coresh_provenance.csv.",
                  length(seeded)))
} else {
  message("[08b] NOTE: coresh_provenance.csv absent — skipping seed cross-check.")
}

## Coerce the four flags to logical (tolerate string "true"/"false" from odd JSON).
to_lgl <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "1")
}
for (b in c("interferon_related", "viral", "cgas_sting", "thermal_heat_shock"))
  ann[[b]] <- to_lgl(ann[[b]])

# ============================================================================
# 3. Derive query family + assemble the tracked table
# ============================================================================
## query_family collapses the gate variants: the two Interaction gates share the
## cGAS-dependent ISG expectation; the two WT_heat gates share the generic-thermal one.

ann <- ann %>%
  mutate(
    query_family = ifelse(grepl("Interaction", recovering_query, fixed = TRUE),
                          "Interaction_up (cGAS-dependent ISG)",
                          "WT_heat_up (generic thermal)"),
    incidental_ifn_note = ifelse(is.na(incidental_ifn_note), "", incidental_ifn_note),
    pubmed_id           = ifelse(is.na(pubmed_id), "", as.character(pubmed_id)),
    source_tool         = SOURCE_TOOL,
    source_prompt       = SOURCE_PROMPT,
    retrieved_date      = RETRIEVED_DATE
  ) %>%
  arrange(query_family, recovering_query, gse) %>%
  select(gse, recovering_query, query_family, title, organism, tissue_cell,
         perturbation, context_class, interferon_related, viral, cgas_sting,
         thermal_heat_shock, incidental_ifn_note, evidence, pubmed_id,
         source_url, confidence, source_tool, source_prompt, retrieved_date)

dir.create(ov_tbl_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(ann, out_fp)
message(sprintf("[08b] wrote %s (%d rows).", out_fp, nrow(ann)))

# ============================================================================
# 4. Console landscape — what each query signature recovered
# ============================================================================
## For each query signature, tabulate the context_class of the datasets CoReSh surfaced
## and list them, so the reader can see what biology each signature co-regulates with.

cat("\n=== CoReSh recovered-dataset landscape (context_class by query family) ===\n")
land <- ann %>%
  count(query_family, context_class, name = "n") %>%
  arrange(query_family, desc(n))
print(as.data.frame(land))

for (qf in unique(ann$query_family)) {
  sub <- ann[ann$query_family == qf, ]
  cat(sprintf("\n%s  (%d datasets)\n", qf, nrow(sub)))
  for (i in seq_len(nrow(sub))) {
    flags <- c(if (sub$interferon_related[i]) "IFN",
               if (sub$viral[i]) "viral",
               if (sub$cgas_sting[i]) "cGAS/STING",
               if (sub$thermal_heat_shock[i]) "thermal")
    cat(sprintf("  %-10s %-22s [%s] %s\n",
                sub$gse[i], sub$context_class[i],
                paste(flags, collapse = ","),
                substr(sub$title[i], 1, 60)))
  }
}

## Flag counts per query family.
cat("\n--- flag counts per query family ---\n")
tally <- ann %>%
  group_by(query_family) %>%
  summarise(n = n(),
            interferon_related = sum(interferon_related),
            viral = sum(viral),
            cgas_sting = sum(cgas_sting),
            thermal = sum(thermal_heat_shock),
            .groups = "drop")
print(as.data.frame(tally))

cat(sprintf("\n[08b] COMPLETE — annotation table at %s\n", out_fp))
