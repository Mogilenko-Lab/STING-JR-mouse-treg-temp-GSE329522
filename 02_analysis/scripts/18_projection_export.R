# 18_projection_export.R — COMPUTE
# =============================================================================
# Human projection export (stage 11_projection) — FREEZE the mouse→human signature
# contract that every human compartment (Phases 1–4) reads.
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   11_projection   (script 18)
#
# ROLE: COMPUTE ONLY. No new statistics. It ORTHOLOG-MAPS the mouse-symbol sets that
#   17_signature_derive.R froze (17_signature_sets.rds) into human HGNC symbols using
#   the pinned/offline babelgene (via ortholog_utils.R) and writes the published,
#   byte-stable contract under 03_results/human_projection/. There is NO ggplot/ggsave
#   here; figures live in 18_projection_export_viz.R.
#
# ┏━ BREAKPOINT-10 GATE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃ This script REFUSES TO RUN until the human signs off on the signature          ┃
# ┃ definitions. It stops unless analysis_config.yaml::decisions.projection.status ┃
# ┃ == "APPROVED". The freeze is downstream-load-bearing (Phases 1–4 depend on     ┃
# ┃ these exact paths/columns), so it must not happen before the review notebook   ┃
# ┃ (notebooks/17_signature_review/) is signed off. See phase-00 §BREAKPOINT 10.   ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#
# Inputs:
#   03_results/objects/17_signature_sets.rds          (mouse-symbol sets + ranked lists)
#   analysis_config.yaml::decisions.projection.*       (the human's BREAKPOINT-10 choices)
#
# Collision policy (from decisions.projection.ortholog_ambiguity; ortholog_utils encodes it):
#   one mouse → many human : binary sets UNION; ranked list assigns the mouse t to EACH human.
#   many mouse → one human : ranked keeps MAX |t| per human; binary dedupe by union.
#   no human ortholog      : dropped and logged (auditable).
#   Interaction (1 df, underpowered): if its post-mapping human set is trivially small
#     (< trivial_min_genes) and drop_interaction_if_trivial, it is DEMOTED from the export
#     and flagged in SIGNATURES.md rather than shipped hollow.
#
# Outputs — the FROZEN CONTRACT (stable API; Phases 1–4 read these paths verbatim):
#   03_results/human_projection/SIGNATURES.md                         (human-readable provenance)
#   03_results/human_projection/manifest.csv                          (machine index: 1 row / (contrast,direction))
#   03_results/human_projection/ortholog_map.tsv                      (the applied map)
#   03_results/human_projection/signatures/<contrast>/<contrast>_up.txt|_down.txt|_ranked.rnk
#   + stage tables for the viz sibling:
#   03_results/11_projection/tables/_overview/{human_signature_sizes,mapping_loss}.csv
#
# IDEMPOTENT + BYTE-STABLE: pure read->map->round->write. Re-running yields identical bytes.
#
# Run from project root (only after sign-off):
#   Rscript 02_analysis/scripts/18_projection_export.R
# =============================================================================

source("02_analysis/config/config.R")            # PROJECT_ROOT, YAML_CONFIG, DIR_*, %||%
source("02_analysis/helpers/de_gsea_helpers.R")   # round_numeric_cols
source("02_analysis/helpers/ortholog_utils.R")    # build_ortholog_map, map_binary_set, map_ranked_list, babelgene_provenance

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

STAGE  <- "11_projection"
SCRIPT <- "02_analysis/scripts/18_projection_export.R"

# ============================================================================
# 0. BREAKPOINT-10 GATE — refuse to freeze until the human has signed off.
# ============================================================================

dcn <- YAML_CONFIG$decisions$projection
if (is.null(dcn))
  stop("[18] decisions.projection missing from analysis_config.yaml — register it and sign off at BREAKPOINT 10 first.")
status <- as.character(dcn$status %||% "PROPOSED_DEFAULT")
if (!identical(status, "APPROVED"))
  stop("[18] BREAKPOINT 10 NOT signed off: decisions.projection.status = '", status,
       "' (need 'APPROVED'). Render notebooks/17_signature_review/, get human sign-off on the ",
       "contrast set / gate / ortholog-ambiguity policy, set status: APPROVED, then re-run.")

# ============================================================================
# 1. LOAD the frozen mouse-symbol sets + the human's decisions.
# ============================================================================

rds <- file.path(DIR_OBJECTS, "17_signature_sets.rds")
if (!file.exists(rds))
  stop("[18] 17_signature_sets.rds not found — run 17_signature_derive.R first.")
sig <- readRDS(rds)

CONTRASTS_EXPORT <- unlist(dcn$contrasts_primary %||% list("WT_heat"))
ROLE_PRIMARY     <- unlist(dcn$role_primary %||% list("WT_heat"))
GATE             <- as.character(dcn$gate %||% "fdr_logfc")
amb              <- dcn$ortholog_ambiguity %||% list()
MIN_SUPPORT      <- as.integer(amb$min_support %||% 3L)
DROP_INT_TRIVIAL <- isTRUE(amb$drop_interaction_if_trivial %||% TRUE)
TRIVIAL_MIN      <- as.integer(amb$trivial_min_genes %||% 10L)

if (!GATE %in% sig$gates)
  stop("[18] decisions.projection.gate='", GATE, "' not among derived gates: ",
       paste(sig$gates, collapse = ", "))
missing_co <- setdiff(CONTRASTS_EXPORT, sig$contrasts)
if (length(missing_co))
  stop("[18] contrasts_primary names absent from 17_signature_sets.rds: ",
       paste(missing_co, collapse = ", "))

role_of <- function(co) ifelse(co %in% ROLE_PRIMARY, "primary", "comparator")

message("=================================================================")
message("18_projection_export: FREEZE mouse->human contract (11_projection)")
message("  status=APPROVED  gate=", GATE, "  export=", paste(CONTRASTS_EXPORT, collapse = ", "))
message("  min_support=", MIN_SUPPORT, "  drop_interaction_if_trivial=", DROP_INT_TRIVIAL,
        " (< ", TRIVIAL_MIN, ")")
message("=================================================================")

# ============================================================================
# 2. BUILD the applied ortholog map ONCE over every mouse gene that appears in any
#    exported artifact (the ranked lists span the DE universe), then reuse it for
#    binary sets + ranked lists. This is the ortholog_map.tsv the contract publishes.
# ============================================================================

mouse_universe <- sort(unique(unlist(lapply(CONTRASTS_EXPORT, function(co)
  sig$sets[[co]]$ranked$gene_symbol), use.names = FALSE)))
message("[map] building applied ortholog map over ", length(mouse_universe),
        " mouse symbols (babelgene min_support=", MIN_SUPPORT, ") ...")
omap <- build_ortholog_map(mouse_universe, min_support = MIN_SUPPORT)
prov <- babelgene_provenance(min_support = MIN_SUPPORT)
message(sprintf("  [map] edges=%d  mapped mouse genes=%d  babelgene=%s (data %s)",
                nrow(omap), length(unique(omap$mouse_symbol)), prov$version, prov$data_date))

# helper: count how many mouse genes in a set are unmapped / one→many (for the manifest).
n_unmapped   <- function(ms) sum(!(unique(ms) %in% omap$mouse_symbol))
n_many_mapped <- function(ms) {
  ms <- unique(ms); tt <- table(omap$mouse_symbol[omap$mouse_symbol %in% ms])
  sum(ms %in% names(tt)[tt > 1L])
}

# ============================================================================
# 3. MAP each exported contrast; assemble manifest rows + mapping-loss rows; decide
#    the trivial-Interaction demotion. Files are written after the demotion check.
# ============================================================================

hp_dir  <- file.path(DIR_RESULTS, "human_projection")
sig_dir <- file.path(hp_dir, "signatures")

per_contrast <- list()      # holds the mapped payload per exported contrast
manifest_rows <- list()
mapping_loss_rows <- list()
demoted <- character(0)

for (co in CONTRASTS_EXPORT) {
  up_mouse   <- sig$sets[[co]]$up[[GATE]]
  down_mouse <- sig$sets[[co]]$down[[GATE]]
  ranked_mouse <- sig$sets[[co]]$ranked            # data.frame(gene_symbol, t)

  up_human   <- map_binary_set(up_mouse, omap)
  down_human <- map_binary_set(down_mouse, omap)
  ranked_human <- map_ranked_list(ranked_mouse, omap, mouse_col = "gene_symbol", stat_col = "t")

  # trivial-Interaction demotion: judged on the mapped up+down human set size.
  human_updown <- length(unique(c(up_human, down_human)))
  if (identical(co, "Interaction") && DROP_INT_TRIVIAL && human_updown < TRIVIAL_MIN) {
    demoted <- c(demoted, co)
    message(sprintf("  [DEMOTE] Interaction post-mapping up+down human set = %d (< %d) — ",
                    human_updown, TRIVIAL_MIN),
            "demoted from the export (flagged in SIGNATURES.md), not shipped hollow.")
    next
  }

  per_contrast[[co]] <- list(
    role = role_of(co),
    up_human = up_human, down_human = down_human, ranked_human = ranked_human,
    up_mouse = up_mouse, down_mouse = down_mouse)

  # manifest rows (one per direction: up, down, ranked)
  add_manifest <- function(direction, n_mouse, n_human, ms, file_rel) {
    manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
      contrast = co, role = role_of(co), direction = direction, gate = GATE,
      n_mouse = n_mouse, n_human = n_human,
      n_dropped_no_ortholog = n_unmapped(ms), n_many_mapped = n_many_mapped(ms),
      file = file_rel, stringsAsFactors = FALSE)
  }
  add_manifest("up",     length(up_mouse),           length(up_human),
               up_mouse,   file.path("signatures", co, paste0(co, "_up.txt")))
  add_manifest("down",   length(down_mouse),         length(down_human),
               down_mouse, file.path("signatures", co, paste0(co, "_down.txt")))
  add_manifest("ranked", nrow(ranked_mouse),         nrow(ranked_human),
               ranked_mouse$gene_symbol, file.path("signatures", co, paste0(co, "_ranked.rnk")))

  # mapping-loss rows for the viz (per contrast x direction: mouse in -> human out breakdown)
  for (dir in c("up", "down")) {
    ms <- if (dir == "up") up_mouse else down_mouse
    hs <- if (dir == "up") up_human else down_human
    mapping_loss_rows[[length(mapping_loss_rows) + 1L]] <- data.frame(
      contrast = co, role = role_of(co), direction = dir, gate = GATE,
      n_mouse = length(ms),
      n_unmapped = n_unmapped(ms), n_many_mapped = n_many_mapped(ms),
      n_human = length(hs), stringsAsFactors = FALSE)
  }
}

if (length(per_contrast) == 0L)
  stop("[18] every exported contrast was demoted/empty — nothing to freeze. Revisit the gate/contrast set.")

# ============================================================================
# 4. WRITE the frozen contract (byte-stable). Own the namespace: clear a stale
#    signatures/ tree first so a re-run with a changed contrast set leaves nothing behind.
# ============================================================================

if (dir.exists(sig_dir)) unlink(sig_dir, recursive = TRUE)
dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

write_lines_sorted <- function(x, path) writeLines(x, path)   # already sorted upstream

for (co in names(per_contrast)) {
  pc  <- per_contrast[[co]]
  cdir <- file.path(sig_dir, co)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(pc$up_human,   file.path(cdir, paste0(co, "_up.txt")))
  writeLines(pc$down_human, file.path(cdir, paste0(co, "_down.txt")))
  # .rnk: 2-col TSV, NO header (GSEA convention), signed t (9 sig figs), descending.
  rnk <- pc$ranked_human
  rnk$t <- signif(rnk$t, 9L)
  readr::write_tsv(rnk, file.path(cdir, paste0(co, "_ranked.rnk")), col_names = FALSE)
}

# manifest.csv (machine-readable index)
manifest <- dplyr::bind_rows(manifest_rows)
manifest_path <- file.path(hp_dir, "manifest.csv")
readr::write_csv(round_numeric_cols(manifest), manifest_path)

# ortholog_map.tsv (the applied map)
omap_out <- omap[, c("mouse_symbol", "mouse_ensembl", "human_symbol",
                     "mapping_type", "babelgene_version")]
readr::write_tsv(omap_out, file.path(hp_dir, "ortholog_map.tsv"))

# stage tables for the viz sibling
tbl_dir <- stage_dir(STAGE, "tables")
ov_dir  <- file.path(tbl_dir, YAML_CONFIG$figures$overview_dir %||% "_overview")
dir.create(ov_dir, recursive = TRUE, showWarnings = FALSE)

human_sizes <- manifest %>%
  dplyr::filter(direction %in% c("up", "down")) %>%
  dplyr::transmute(contrast, role, direction, gate, n_human)
readr::write_csv(round_numeric_cols(human_sizes), file.path(ov_dir, "human_signature_sizes.csv"))

mapping_loss <- dplyr::bind_rows(mapping_loss_rows) %>%
  dplyr::mutate(n_mapped_1to1 = n_mouse - n_unmapped - n_many_mapped)
readr::write_csv(round_numeric_cols(mapping_loss), file.path(ov_dir, "mapping_loss.csv"))

# ============================================================================
# 5. SIGNATURES.md — human-readable provenance manifest.
# ============================================================================

git_sha <- tryCatch(
  suppressWarnings(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE))[1],
  error = function(e) NA_character_) %||% NA_character_

md <- c(
  "# Mouse → human projection signatures (frozen contract)",
  "",
  "Frozen mouse 39 °C Treg-stress signature export from the internal bulk RNA-seq anchor",
  "(GSE329522 iTreg, 2×2 genotype × temperature). Every human compartment (Phases 1–4)",
  "scores against THIS artifact. It is a reshape + ortholog map of the finished DE master —",
  "**no new biology**. Language is correlative: these are the INPUT the human phases test",
  "for presence, never a claim that fever/HIF/STING drives human disease.",
  "",
  "## Provenance chain",
  "",
  "`03_results/master/master_de_genes.csv` → `03_results/objects/17_signature_sets.rds`",
  "(17_signature_derive.R) → this contract (18_projection_export.R).",
  sprintf("- git SHA: `%s`", git_sha),
  sprintf("- built_at: %s UTC", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")),
  "",
  "## Thresholds / gate applied",
  "",
  sprintf("- gate: **%s** (binary up/down sets)", GATE),
  sprintf("- de_fdr = %s ; de_logfc = %s", sig$thresholds$de_fdr, sig$thresholds$de_logfc),
  sprintf("- rank metric: signed **%s** (ranked .rnk lists; NEVER logFC)", sig$thresholds$rank_metric),
  "",
  "## Ortholog source + collision policy",
  "",
  sprintf("- %s %s, bundled data %s; direction: %s; min_support = %d.",
          prov$package, prov$version, prov$data_date, prov$mapping_dir, prov$min_support),
  "- one mouse → many human: binary sets UNION; ranked assigns the mouse t to each human ortholog.",
  "- many mouse → one human: ranked keeps MAX |t| per human symbol; binary dedupe by union.",
  "- no human ortholog: dropped (see `ortholog_map.tsv` for every applied edge; unmapped genes absent).",
  "",
  "## Exported contrasts",
  "")

# per-contrast provenance block, from the manifest
for (co in names(per_contrast)) {
  m_up   <- manifest %>% dplyr::filter(contrast == co, direction == "up")
  m_dn   <- manifest %>% dplyr::filter(contrast == co, direction == "down")
  m_rk   <- manifest %>% dplyr::filter(contrast == co, direction == "ranked")
  lab    <- tryCatch(YAML_CONFIG$design$contrast_labels[[co]], error = function(e) co) %||% co
  note   <- if (identical(co, "Geno_at_39"))
              " (genotype proxy for the death/vulnerability-at-39°C axis; the inhibitor arm has no GEO matrix — NOT an inhibitor signature)"
            else if (identical(co, "Interaction"))
              " (tests cGAS-dependence of the heat response; 1 df, underpowered at n=5 — labelled by what it TESTS, not a result)"
            else ""
  md <- c(md,
    sprintf("### %s — role: %s", co, role_of(co)),
    sprintf("- definition: %s%s", lab, note),
    sprintf("- up:   %d mouse → %d human (dropped %d, many-mapped %d)",
            m_up$n_mouse, m_up$n_human, m_up$n_dropped_no_ortholog, m_up$n_many_mapped),
    sprintf("- down: %d mouse → %d human (dropped %d, many-mapped %d)",
            m_dn$n_mouse, m_dn$n_human, m_dn$n_dropped_no_ortholog, m_dn$n_many_mapped),
    sprintf("- ranked: %d mouse → %d human genes (`%s_ranked.rnk`, signed t, descending)",
            m_rk$n_mouse, m_rk$n_human, co),
    "")
}

if (length(demoted)) {
  md <- c(md, "## Demoted (not shipped)", "",
    sprintf("- **%s**: post-mapping human up+down set < %d genes (trivially small); demoted per decisions.projection.ortholog_ambiguity.drop_interaction_if_trivial. Not shipped hollow.",
            paste(demoted, collapse = ", "), TRIVIAL_MIN),
    "")
}

md <- c(md,
  "## Files",
  "",
  "- `manifest.csv` — one row per (contrast, direction): role, gate, n_mouse, n_human, n_dropped_no_ortholog, n_many_mapped, file.",
  "- `ortholog_map.tsv` — the applied map: mouse_symbol, mouse_ensembl, human_symbol, mapping_type, babelgene_version.",
  "- `signatures/<contrast>/<contrast>_up.txt` / `_down.txt` — human HGNC symbols, one per line (AUCell/UCell).",
  "- `signatures/<contrast>/<contrast>_ranked.rnk` — 2-col TSV (human_symbol⇥t), signed, descending (fgsea/decoupleR).",
  "")

writeLines(md, file.path(hp_dir, "SIGNATURES.md"))

# ============================================================================
# 6. ACCEPTANCE CHECKS (structural) + console summary.
# ============================================================================

stopifnot(
  "SIGNATURES.md missing"   = file.exists(file.path(hp_dir, "SIGNATURES.md")),
  "manifest.csv missing"    = file.exists(manifest_path),
  "ortholog_map.tsv missing"= file.exists(file.path(hp_dir, "ortholog_map.tsv")))
for (co in names(per_contrast)) stopifnot(
  file.exists(file.path(sig_dir, co, paste0(co, "_up.txt"))),
  file.exists(file.path(sig_dir, co, paste0(co, "_down.txt"))),
  file.exists(file.path(sig_dir, co, paste0(co, "_ranked.rnk"))))

# WT_heat sanity: the primary up+down human set must be non-trivial (not decimated).
if ("WT_heat" %in% names(per_contrast)) {
  wh <- length(unique(c(per_contrast$WT_heat$up_human, per_contrast$WT_heat$down_human)))
  message(sprintf("  [check] WT_heat primary up+down human set = %d genes.", wh))
  if (wh < TRIVIAL_MIN)
    warning("[18] WT_heat human set is trivially small (", wh, ") — inspect mapping_loss before shipping.")
}

message("=================================================================")
message("18_projection_export COMPLETE — contract frozen at ", hp_dir)
message("  Exported: ", paste(names(per_contrast), collapse = ", "),
        if (length(demoted)) paste0("  (demoted: ", paste(demoted, collapse = ", "), ")") else "")
message("  Run 18_projection_export_viz.R for the human-space size + mapping-loss figures.")
message("  Next: add a DATA_MANIFEST.md pointer + note completion in docs/_internal.")
message("=================================================================")
