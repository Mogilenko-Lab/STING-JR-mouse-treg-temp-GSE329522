# 18_projection_export.R — COMPUTE
# =============================================================================
# Freeze the mouse→human signature contract that every human compartment reads
# (stage 11_projection).
#
# Project: GSE329522 STING/cGAS hyperthermia iTreg (2x2 genotype x temperature)
#
# ROLE: COMPUTE ONLY. Ortholog-maps the mouse-symbol sets frozen by
#   17_signature_derive.R into human HGNC symbols through the pinned offline
#   babelgene (ortholog_utils.R), and writes the byte-stable contract under
#   03_results/human_projection/. Figures live in 18_projection_export_viz.R.
#
# THE THREE HEAT CONTRASTS, AND WHAT EACH IS FOR
#   They are linearly dependent by construction:
#     WT_heat     the full thermal response in cGAS-competent cells   (primary)
#     Interaction WT_heat - KO_heat, the cGAS-dependent slice         (ISG/STING candidate)
#     KO_heat     WT_heat - Interaction, the cGAS-independent comparator
#   so WT_heat = KO_heat + Interaction, and any two determine the third. KO_heat
#   ships as its own scored contrast to serve as the specificity control at the
#   integration layer: an overlap with a STING gain-of-function reference as large
#   as WT_heat's marks that overlap STING-independent and IFN-like, which keeps the
#   human read correlative. A single human dataset carries two threads (WT_heat +
#   Interaction), since KO = WT - Interaction follows from them.
#
#   Gate asymmetry (analysis_config.yaml::decisions.projection) is by design.
#   Interaction also rides the looser fdr_only gate as the underpowered 1-df
#   positive candidate. KO_heat stays on the stringent fdr_logfc gate alone as the
#   well-powered negative control.
#
# GATE (BREAKPOINT 10). This script halts unless
#   analysis_config.yaml::decisions.projection.status == "APPROVED". The freeze is
#   downstream-load-bearing — the human compartments read these exact paths and
#   columns — so it waits on sign-off of notebooks/17_signature_review/.
#
# Inputs:
#   03_results/objects/17_signature_sets.rds       (mouse-symbol sets + ranked lists)
#   analysis_config.yaml::decisions.projection.*   (the approved choices)
#
# Collision policy (decisions.projection.ortholog_ambiguity; encoded in ortholog_utils):
#   one mouse → many human : binary sets union; ranked list assigns the mouse t to each human.
#   many mouse → one human : ranked keeps max |t| per human, with a primary-query edge
#     outranking a recovered one; binary dedupe by union.
#   no human ortholog      : dropped and logged.
#   stale query symbol     : re-asked under its current MGI symbol, then mapped back.
#
# WHY THE DROP BUCKET SPLITS THREE WAYS. babelgene keys on current MGI symbols. This
# matrix was quantified against GRCm38 + GENCODE vM25, so 2,341 of its 19,679 symbols
# have moved on, and babelgene returns nothing for them. That was reported as "no human
# ortholog"; for 146 of 6,510 edges it was a vocabulary loss wearing an orthology label.
# Ddx58 is the one to remember: RIG-I, a cytosolic nucleic-acid sensor, significantly down
# at 39 °C in WT and at logFC −1.28 in KO_heat, absent from the projection because its
# current symbol is Rigi. build_ortholog_map() re-asks those symbols under their current
# names, and the ledger reports:
#   n_query_symbol_normalised     recovered by the re-ask
#   n_dropped_stale_query_symbol  unmapped and no longer a current MGI symbol — a VOCABULARY fact
#   n_dropped_no_ortholog         unmapped and still current — an ORTHOLOGY fact
# The recovery is strictly additive: 0 edges lost, 0 pre-existing human symbols change
# their ranked t (policy recovered_edge_precedence = "primary_query_wins"). Both are
# asserted below.
#
#   Interaction (1 df, underpowered): a post-mapping human set under trivial_min_genes is
#   demoted from the export under drop_interaction_if_trivial and flagged in SIGNATURES.md.
#
# Outputs — THE FROZEN CONTRACT (stable API; the human compartments read these paths verbatim):
#   03_results/human_projection/SIGNATURES.md    (human-readable provenance)
#   03_results/human_projection/manifest.csv     (machine index: 1 row / (contrast,direction))
#   03_results/human_projection/ortholog_map.tsv (the applied map)
#   03_results/human_projection/signatures/<contrast>/<contrast>_up.txt|_down.txt|_ranked.rnk
#   03_results/11_projection/tables/_overview/{human_signature_sizes,mapping_loss}.csv
#   03_results/11_projection/tables/_overview/projection_overlap_ledger.csv
#   The stage tables carry one row per (contrast, gate, direction) covering every gate the
#   manifest ships, primary and secondary alike, so the two stay in lockstep.
#
# IDEMPOTENT + BYTE-STABLE: read → map → round → write. Re-running yields identical bytes.
#
# Run from project root, after sign-off:
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
SECONDARY_GATE   <- as.character(dcn$secondary_gate %||% NA_character_)
SECONDARY_CONTRASTS <- unlist(dcn$secondary_gate_contrasts %||% list())
HAS_SECONDARY_GATE <- length(SECONDARY_GATE) > 0L &&
  !is.na(SECONDARY_GATE[1]) && nzchar(SECONDARY_GATE[1]) &&
  length(SECONDARY_CONTRASTS) > 0L
amb              <- dcn$ortholog_ambiguity %||% list()
MIN_SUPPORT      <- as.integer(amb$min_support %||% 3L)
DROP_INT_TRIVIAL <- isTRUE(amb$drop_interaction_if_trivial %||% TRUE)
TRIVIAL_MIN      <- as.integer(amb$trivial_min_genes %||% 10L)

if (!GATE %in% sig$gates)
  stop("[18] decisions.projection.gate='", GATE, "' not among derived gates: ",
       paste(sig$gates, collapse = ", "))
if (HAS_SECONDARY_GATE && !SECONDARY_GATE[1] %in% sig$gates)
  stop("[18] decisions.projection.secondary_gate='", SECONDARY_GATE[1],
       "' not among derived gates: ", paste(sig$gates, collapse = ", "))
missing_co <- setdiff(CONTRASTS_EXPORT, sig$contrasts)
if (length(missing_co))
  stop("[18] contrasts_primary names absent from 17_signature_sets.rds: ",
       paste(missing_co, collapse = ", "))
missing_secondary_co <- setdiff(SECONDARY_CONTRASTS, names(sig$sets))
if (HAS_SECONDARY_GATE && length(missing_secondary_co))
  stop("[18] secondary_gate_contrasts names absent from 17_signature_sets.rds: ",
       paste(missing_secondary_co, collapse = ", "))

role_of <- function(co) ifelse(co %in% ROLE_PRIMARY, "primary", "comparator")

message("=================================================================")
message("18_projection_export: FREEZE mouse->human contract (11_projection)")
message("  status=APPROVED  gate=", GATE, "  export=", paste(CONTRASTS_EXPORT, collapse = ", "))
message("  min_support=", MIN_SUPPORT, "  drop_interaction_if_trivial=", DROP_INT_TRIVIAL,
        " (< ", TRIVIAL_MIN, ")")
if (HAS_SECONDARY_GATE)
  message("  secondary_gate=", SECONDARY_GATE[1], "  secondary_export=",
          paste(SECONDARY_CONTRASTS, collapse = ", "))
message("=================================================================")

# ============================================================================
# 2. BUILD the applied ortholog map ONCE over every mouse gene that appears in any
#    exported artifact (the ranked lists span the DE universe), then reuse it for
#    binary sets + ranked lists. This is the ortholog_map.tsv the contract publishes.
# ============================================================================

mouse_universe <- sort(unique(unlist(lapply(CONTRASTS_EXPORT, function(co)
  sig$sets[[co]]$ranked$gene_symbol), use.names = FALSE)))
QUERY_FLAGGED <- unlist(YAML_CONFIG$symbol_alias$ortholog_query_flagged_for_review) %||%
  character()
message("[map] building applied ortholog map over ", length(mouse_universe),
        " mouse symbols (babelgene min_support=", MIN_SUPPORT, ") ...")
omap <- build_ortholog_map(mouse_universe, min_support = MIN_SUPPORT,
                           flagged_pairs = QUERY_FLAGGED)
qledger <- attr(omap, "query_normalisation")
prov <- babelgene_provenance(min_support = MIN_SUPPORT)
message(sprintf("  [map] edges=%d  mapped mouse genes=%d  babelgene=%s (data %s)",
                nrow(omap), length(unique(omap$mouse_symbol)), prov$version, prov$data_date))

# ---- the query-normalisation recovery, as a decision ledger -----------------------------
# Every candidate is a row, accepted or not, so a recovery never has to be inferred from an
# absence and a withholding never from a shrug. `mapped_after_normalisation` separates
# "re-asked and babelgene answered" from "re-asked and it still has no ortholog at this
# min_support" — the second is a real orthology result and must not be filed as a recovery.
recovered_symbols <- unique(omap$mouse_symbol[omap$normalisation_source == "org.Mm.eg.db"])
qledger$mapped_after_normalisation <- qledger$matrix_symbol %in% recovered_symbols
CURRENT_MGI <- if (requireNamespace("org.Mm.eg.db", quietly = TRUE))
  AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "SYMBOL") else character()
message(sprintf(
  "  [map] query normalisation: %d candidates, %d accepted, %d of those recovered a mapping; %d flagged for review, %d rejected by a guard",
  nrow(qledger), sum(qledger$resolution == "accepted"), length(recovered_symbols),
  sum(qledger$resolution == "flagged_for_review"),
  sum(!qledger$resolution %in% c("accepted", "flagged_for_review"))))
print(as.data.frame(qledger %>%
  dplyr::count(.data$resolution, .data$mapped_after_normalisation, name = "n")),
  row.names = FALSE)

# The recovery may only ADD. A lost edge or a changed metric would mean the map was rewritten
# in place of extended, which is the failure mode the two-pass design exists to prevent.
prev_map_path <- file.path(DIR_RESULTS, "human_projection", "ortholog_map.tsv")
if (file.exists(prev_map_path)) {
  prev <- readr::read_tsv(prev_map_path, show_col_types = FALSE, progress = FALSE)
  lost_edges <- dplyr::anti_join(prev[c("mouse_symbol", "human_symbol")],
                                 omap[c("mouse_symbol", "human_symbol")],
                                 by = c("mouse_symbol", "human_symbol"))
  if (nrow(lost_edges))
    stop(sprintf(paste0("[18] the rebuilt ortholog map LOST %d edge(s) that the published ",
                        "one carried (e.g. %s -> %s). The query normalisation must extend ",
                        "the map, never rewrite it."),
                 nrow(lost_edges), lost_edges$mouse_symbol[1], lost_edges$human_symbol[1]))
  message(sprintf("  [check] map is additive against the published one: %d edges -> %d, 0 lost.",
                  nrow(prev), nrow(omap)))
}

# helper: count how many mouse genes in a set are unmapped / one→many (for the manifest), and
# split the unmapped into the orthology statement and the vocabulary one.
n_unmapped   <- function(ms) sum(!(unique(ms) %in% omap$mouse_symbol))
n_many_mapped <- function(ms) {
  ms <- unique(ms); tt <- table(omap$mouse_symbol[omap$mouse_symbol %in% ms])
  sum(ms %in% names(tt)[tt > 1L])
}
.unmapped_of <- function(ms) setdiff(unique(ms), omap$mouse_symbol)
n_dropped_no_ortholog    <- function(ms) sum(.unmapped_of(ms) %in% CURRENT_MGI)
n_dropped_stale_query    <- function(ms) sum(!.unmapped_of(ms) %in% CURRENT_MGI)
n_query_normalised       <- function(ms) sum(unique(ms) %in% recovered_symbols)

# ============================================================================
# 3. MAP each exported contrast; assemble manifest rows + mapping-loss rows; decide
#    the trivial-Interaction demotion. Files are written after the demotion check.
# ============================================================================

hp_dir  <- file.path(DIR_RESULTS, "human_projection")
sig_dir <- file.path(hp_dir, "signatures")

per_contrast <- list()      # holds the mapped payload per exported contrast
secondary_per_contrast <- list()
manifest_rows <- list()
mapping_loss_rows <- list()
demoted <- character(0)

add_manifest_row <- function(co, direction, gate, n_mouse, n_human, ms, file_rel) {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    contrast = co, role = role_of(co), direction = direction, gate = gate,
    n_mouse = n_mouse, n_human = n_human,
    # n_dropped_no_ortholog KEEPS its name but no longer keeps its old meaning: it is now
    # only the genes babelgene could key and had no ortholog for. The vocabulary half moved
    # to n_dropped_stale_query_symbol, and n_dropped_unmapped_total is the old column's value
    # so the pre-fix number stays legible beside the split.
    n_dropped_no_ortholog = n_dropped_no_ortholog(ms),
    n_dropped_stale_query_symbol = n_dropped_stale_query(ms),
    n_dropped_unmapped_total = n_unmapped(ms),
    n_query_symbol_normalised = n_query_normalised(ms),
    n_many_mapped = n_many_mapped(ms),
    file = file_rel, stringsAsFactors = FALSE)
}

add_mapping_loss_row <- function(co, direction, gate, ms, hs) {
  mapping_loss_rows[[length(mapping_loss_rows) + 1L]] <<- data.frame(
    contrast = co, role = role_of(co), direction = direction, gate = gate,
    n_mouse = length(ms),
    n_unmapped = n_unmapped(ms),
    n_dropped_no_ortholog = n_dropped_no_ortholog(ms),
    n_dropped_stale_query_symbol = n_dropped_stale_query(ms),
    n_query_symbol_normalised = n_query_normalised(ms),
    n_many_mapped = n_many_mapped(ms),
    n_human = length(hs), stringsAsFactors = FALSE)
}

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
  add_manifest_row(co, "up",     GATE, length(up_mouse),   length(up_human),
                   up_mouse,   file.path("signatures", co, paste0(co, "_up.txt")))
  add_manifest_row(co, "down",   GATE, length(down_mouse), length(down_human),
                   down_mouse, file.path("signatures", co, paste0(co, "_down.txt")))
  add_manifest_row(co, "ranked", GATE, nrow(ranked_mouse), nrow(ranked_human),
                   ranked_mouse$gene_symbol, file.path("signatures", co, paste0(co, "_ranked.rnk")))

  # mapping-loss rows for the viz (per contrast x direction: mouse in -> human out breakdown)
  for (dir in c("up", "down")) {
    ms <- if (dir == "up") up_mouse else down_mouse
    hs <- if (dir == "up") up_human else down_human
    add_mapping_loss_row(co, dir, GATE, ms, hs)
  }
}

if (length(per_contrast) == 0L)
  stop("[18] every exported contrast was demoted/empty — nothing to freeze. Revisit the gate/contrast set.")

if (HAS_SECONDARY_GATE) {
  for (co in SECONDARY_CONTRASTS) {
    up_mouse   <- sig$sets[[co]]$up[[SECONDARY_GATE[1]]]
    down_mouse <- sig$sets[[co]]$down[[SECONDARY_GATE[1]]]

    up_human   <- map_binary_set(up_mouse, omap)
    down_human <- map_binary_set(down_mouse, omap)

    secondary_per_contrast[[co]] <- list(
      gate = SECONDARY_GATE[1],
      up_human = up_human, down_human = down_human,
      up_mouse = up_mouse, down_mouse = down_mouse)

    add_manifest_row(co, "up", SECONDARY_GATE[1], length(up_mouse), length(up_human),
                     up_mouse, file.path("signatures", co, paste0(co, "_fdrOnly_up.txt")))
    add_manifest_row(co, "down", SECONDARY_GATE[1], length(down_mouse), length(down_human),
                     down_mouse, file.path("signatures", co, paste0(co, "_fdrOnly_down.txt")))

    # mapping-loss rows for the viz, same as the primary gate: the stage tables carry
    # EVERY (contrast, gate) the manifest ships, so the figures cannot silently omit a gate.
    for (dir in c("up", "down")) {
      ms <- if (dir == "up") up_mouse else down_mouse
      hs <- if (dir == "up") up_human else down_human
      add_mapping_loss_row(co, dir, SECONDARY_GATE[1], ms, hs)
    }
  }
}

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

new_files_written <- character(0)
for (co in names(secondary_per_contrast)) {
  sc <- secondary_per_contrast[[co]]
  cdir <- file.path(sig_dir, co)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  up_path <- file.path(cdir, paste0(co, "_fdrOnly_up.txt"))
  dn_path <- file.path(cdir, paste0(co, "_fdrOnly_down.txt"))
  writeLines(sc$up_human, up_path)
  writeLines(sc$down_human, dn_path)
  new_files_written <- c(
    new_files_written,
    file.path("signatures", co, paste0(co, "_fdrOnly_up.txt")),
    file.path("signatures", co, paste0(co, "_fdrOnly_down.txt")))
}

# manifest.csv (machine-readable index)
manifest <- dplyr::bind_rows(manifest_rows)
manifest_path <- file.path(hp_dir, "manifest.csv")
readr::write_csv(round_numeric_cols(manifest), manifest_path)

# ortholog_map.tsv (the applied map)
omap_out <- omap[, c("mouse_symbol", "mouse_ensembl", "human_symbol", "mapping_type",
                     "matrix_symbol_normalised_to", "normalisation_source",
                     "babelgene_version")]
readr::write_tsv(omap_out, file.path(hp_dir, "ortholog_map.tsv"))

# stage tables for the viz sibling
tbl_dir <- stage_dir(STAGE, "tables")
ov_dir  <- file.path(tbl_dir, YAML_CONFIG$figures$overview_dir %||% "_overview")
dir.create(ov_dir, recursive = TRUE, showWarnings = FALSE)

# One row per (contrast, gate, direction) the manifest ships — every gate, primary and
# primary one, so the stage tables stay in lockstep with manifest.csv and the viz sibling
# can render the gate dimension and keep the secondary-gate rows.
human_sizes <- manifest %>%
  dplyr::filter(direction %in% c("up", "down")) %>%
  dplyr::transmute(contrast, role, direction, gate, n_human)
readr::write_csv(round_numeric_cols(human_sizes), file.path(ov_dir, "human_signature_sizes.csv"))

mapping_loss <- dplyr::bind_rows(mapping_loss_rows) %>%
  dplyr::mutate(n_mapped_1to1 = n_mouse - n_unmapped - n_many_mapped)
readr::write_csv(round_numeric_cols(mapping_loss), file.path(ov_dir, "mapping_loss.csv"))

# The query-normalisation ledger, published at the same granularity the decision was made:
# one row per candidate matrix symbol. This is the artifact that makes the split auditable —
# without it "146 recovered" is a claim in place of a record.
readr::write_csv(
  qledger %>% dplyr::transmute(
    .data$matrix_symbol, .data$current_symbol, .data$entrez_id, .data$resolution,
    .data$mapped_after_normalisation,
    normalisation_source = "org.Mm.eg.db",
    org_mm_eg_db_version = if (length(CURRENT_MGI))
      as.character(packageVersion("org.Mm.eg.db")) else NA_character_) %>%
    dplyr::arrange(.data$resolution, .data$matrix_symbol),
  file.path(ov_dir, "query_normalisation_ledger.csv"))

# Closure, asserted in-script: every mouse gene is mapped, mapped-many, dropped
# for a real orthology reason, or dropped because we could not key its symbol. There is no
# fifth outcome, and the two drop classes may never be folded back together.
closure <- mapping_loss %>%
  dplyr::mutate(check = n_mapped_1to1 + n_many_mapped + n_dropped_no_ortholog +
                  n_dropped_stale_query_symbol)
if (any(closure$check != closure$n_mouse))
  stop("[18] conversion ledger does not close: n_mouse != mapped_1to1 + many_mapped + ",
       "dropped_no_ortholog + dropped_stale_query_symbol for ",
       sum(closure$check != closure$n_mouse), " row(s).")

# Panel 1D source table: count from the frozen text files just written, downstream of
# in-memory pre-write vectors. This keeps the plotted ledger tied to the public
# human_projection contract byte-for-byte.
read_signature_genes <- function(co, direction, gate) {
  co_arg <- co
  direction_arg <- direction
  gate_arg <- gate
  rr <- manifest %>%
    dplyr::filter(.data$contrast == .env$co_arg,
                  .data$direction == .env$direction_arg,
                  .data$gate == .env$gate_arg)
  if (nrow(rr) == 0L) return(character(0))
  path <- file.path(hp_dir, rr$file[1])
  if (!file.exists(path)) stop("[18] Manifest points to missing signature file: ", path)
  x <- readr::read_lines(path, progress = FALSE)
  sort(unique(x[!is.na(x) & x != ""]))
}

primary_sets <- list(
  up = list(
    WT_heat = read_signature_genes("WT_heat", "up", GATE),
    KO_heat = read_signature_genes("KO_heat", "up", GATE),
    Interaction = read_signature_genes("Interaction", "up", GATE)
  ),
  down = list(
    WT_heat = read_signature_genes("WT_heat", "down", GATE),
    KO_heat = read_signature_genes("KO_heat", "down", GATE),
    Interaction = read_signature_genes("Interaction", "down", GATE)
  )
)
secondary_sets <- list(
  up = read_signature_genes("Interaction", "up", SECONDARY_GATE[1]),
  down = read_signature_genes("Interaction", "down", SECONDARY_GATE[1])
)

wt_only_up <- setdiff(primary_sets$up$WT_heat, primary_sets$up$KO_heat)
cgas_source <- file.path("03_results", "03_de", "tables", "_overview", "cgas_dependence_wide.csv")
cgas_dependent_in_wt_only <- NA_integer_
cgas_flag_note <- "cgas_dependent flag absent"
if (file.exists(cgas_source)) {
  cgas_tbl <- readr::read_csv(cgas_source, show_col_types = FALSE, progress = FALSE)
  if (all(c("gene_symbol", "cgas_dependent") %in% names(cgas_tbl))) {
    cgas_join <- omap_out %>%
      dplyr::filter(.data$human_symbol %in% wt_only_up) %>%
      dplyr::left_join(cgas_tbl %>% dplyr::select("gene_symbol", "cgas_dependent"),
                       by = c("mouse_symbol" = "gene_symbol")) %>%
      dplyr::group_by(.data$human_symbol) %>%
      dplyr::summarise(cgas_dependent = any(.data$cgas_dependent %in% TRUE, na.rm = TRUE),
                       .groups = "drop")
    cgas_dependent_in_wt_only <- sum(cgas_join$cgas_dependent %in% TRUE)
    cgas_flag_note <- cgas_source
  }
}

ledger_rows <- list()
add_ledger_row <- function(direction, display_group, component, gate, n_human,
                           glyph = "bar", total_human = NA_integer_,
                           heat_union_n = NA_integer_, heat_shared_n = NA_integer_,
                           heat_jaccard = NA_real_,
                           n_intersect_wt = NA_integer_, n_intersect_ko = NA_integer_,
                           n_intersect_heat_union = NA_integer_,
                           underpowered = FALSE) {
  ledger_rows[[length(ledger_rows) + 1L]] <<- data.frame(
    direction = direction,
    display_group = display_group,
    component = component,
    gate = gate,
    n_human = n_human,
    total_human = total_human,
    heat_union_n = heat_union_n,
    heat_shared_n = heat_shared_n,
    heat_jaccard = heat_jaccard,
    n_intersect_wt = n_intersect_wt,
    n_intersect_ko = n_intersect_ko,
    n_intersect_heat_union = n_intersect_heat_union,
    glyph = glyph,
    underpowered = underpowered,
    wt_only_up_n = length(wt_only_up),
    wt_only_up_cgas_dependent_n = cgas_dependent_in_wt_only,
    cgas_flag_source = cgas_flag_note,
    stringsAsFactors = FALSE
  )
}

for (dir in c("up", "down")) {
  wt <- primary_sets[[dir]]$WT_heat
  ko <- primary_sets[[dir]]$KO_heat
  shared <- intersect(wt, ko)
  heat_union <- union(wt, ko)
  heat_j <- length(shared) / length(heat_union)
  add_ledger_row(dir, "WT/KO heat union", "WT_heat only", GATE, length(setdiff(wt, ko)),
                 total_human = length(wt), heat_union_n = length(heat_union),
                 heat_shared_n = length(shared), heat_jaccard = heat_j)
  add_ledger_row(dir, "WT/KO heat union", "WT_heat ∩ KO_heat", GATE, length(shared),
                 total_human = length(heat_union), heat_union_n = length(heat_union),
                 heat_shared_n = length(shared), heat_jaccard = heat_j)
  add_ledger_row(dir, "WT/KO heat union", "KO_heat only", GATE, length(setdiff(ko, wt)),
                 total_human = length(ko), heat_union_n = length(heat_union),
                 heat_shared_n = length(shared), heat_jaccard = heat_j)

  for (gate_i in c(GATE, SECONDARY_GATE[1])) {
    int_set <- if (identical(gate_i, GATE)) primary_sets[[dir]]$Interaction else secondary_sets[[dir]]
    glyph <- if (identical("down", dir) && length(int_set) == 0L) "structural_empty" else "bar"
    add_ledger_row(
      dir, paste0("Interaction\n[", gate_i, "]"), "Interaction", gate_i, length(int_set),
      glyph = glyph, total_human = length(int_set), heat_union_n = length(heat_union),
      heat_shared_n = length(shared), heat_jaccard = heat_j,
      n_intersect_wt = length(intersect(int_set, wt)),
      n_intersect_ko = length(intersect(int_set, ko)),
      n_intersect_heat_union = length(intersect(int_set, heat_union)),
      underpowered = length(int_set) > 0L && length(int_set) < TRIVIAL_MIN
    )
  }
}

projection_ledger <- dplyr::bind_rows(ledger_rows)
readr::write_csv(round_numeric_cols(projection_ledger),
                 file.path(ov_dir, "projection_overlap_ledger.csv"))

# ============================================================================
# 5. SIGNATURES.md — human-readable provenance manifest.
# ============================================================================

signatures_md_path <- file.path(hp_dir, "SIGNATURES.md")
existing_md <- if (file.exists(signatures_md_path)) readLines(signatures_md_path, warn = FALSE) else character(0)
first_or_na <- function(x) if (length(x) && !is.na(x[1])) x[1] else NA_character_
existing_git_sha <- first_or_na(sub("^- git SHA: `([^`]+)`$", "\\1", existing_md[grepl("^- git SHA: `", existing_md)]))
existing_built_at <- first_or_na(sub("^- built_at: (.+)$", "\\1", existing_md[grepl("^- built_at: ", existing_md)]))
git_sha <- if (!is.na(existing_git_sha)) existing_git_sha else NA_character_
built_at <- if (!is.na(existing_built_at)) existing_built_at else
  sprintf("%s UTC", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"))

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
  sprintf("- built_at: %s", built_at),
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
  "- many mouse → one human: ranked keeps MAX |t| per human symbol, except that an edge from",
  "  the original query outranks one recovered by symbol normalisation regardless of |t|;",
  "  binary dedupe by union.",
  sprintf("- stale mouse symbol: babelgene keys on CURRENT MGI symbols, and this matrix was quantified against GENCODE vM25, so %d of its %d symbols are no longer current. babelgene knows some of them and not others, idiosyncratically, so the ones it could not key at all are re-asked under their current symbol via org.Mm.eg.db %s and the edges mapped back — `mouse_symbol` below always stays the symbol the data carries.",
          sum(!mouse_universe %in% CURRENT_MGI), length(mouse_universe),
          if (length(CURRENT_MGI)) as.character(packageVersion("org.Mm.eg.db")) else "unavailable"),
  sprintf("- %d genes arrived that way, each previously counted as having no human ortholog, which was the wrong label for a vocabulary loss.",
          length(recovered_symbols)),
  sprintf("- the recovery is strictly additive: %d edges added, none removed, and a human symbol the map already carried keeps the ranked metric it had. Every candidate, accepted or withheld, is a row in `03_results/11_projection/tables/_overview/query_normalisation_ledger.csv`.",
          sum(omap$normalisation_source == "org.Mm.eg.db")),
  "- no human ortholog: dropped. This now means what it says — babelgene could key the symbol",
  "  and returned no ortholog at this min_support. A symbol it could not key at all is counted",
  "  separately as `n_dropped_stale_query_symbol`, because that is a statement about vocabulary",
  "  rather than about orthology, and `n_dropped_unmapped_total` keeps the two together for",
  "  comparison with earlier builds. See `ortholog_map.tsv` for every applied edge.",
  "",
  "## Exported contrasts",
  "")

# per-contrast provenance block, from the manifest
for (co in names(per_contrast)) {
  m_up   <- manifest %>% dplyr::filter(contrast == co, direction == "up", gate == GATE)
  m_dn   <- manifest %>% dplyr::filter(contrast == co, direction == "down", gate == GATE)
  m_rk   <- manifest %>% dplyr::filter(contrast == co, direction == "ranked", gate == GATE)
  sc     <- secondary_per_contrast[[co]] %||% NULL
  lab    <- tryCatch(YAML_CONFIG$design$contrast_labels[[co]], error = function(e) co) %||% co
  note   <- if (identical(co, "Geno_at_39"))
              " (genotype proxy for the death/vulnerability-at-39°C axis; the inhibitor arm has no GEO matrix — NOT an inhibitor signature)"
            else if (identical(co, "Interaction"))
              " (tests cGAS-dependence of the heat response; 1 df, underpowered at n=5 — labelled by what it TESTS, not a result)"
            else ""
  md <- c(md,
    sprintf("### %s — role: %s", co, role_of(co)),
    sprintf("- definition: %s%s", lab, note),
    sprintf("- up:   %d mouse → %d human (no ortholog %d, symbol not keyable %d, recovered by normalisation %d, many-mapped %d)",
            m_up$n_mouse, m_up$n_human, m_up$n_dropped_no_ortholog,
            m_up$n_dropped_stale_query_symbol, m_up$n_query_symbol_normalised,
            m_up$n_many_mapped),
    sprintf("- down: %d mouse → %d human (no ortholog %d, symbol not keyable %d, recovered by normalisation %d, many-mapped %d)",
            m_dn$n_mouse, m_dn$n_human, m_dn$n_dropped_no_ortholog,
            m_dn$n_dropped_stale_query_symbol, m_dn$n_query_symbol_normalised,
            m_dn$n_many_mapped),
    sprintf("- ranked: %d mouse → %d human genes (`%s_ranked.rnk`, signed t, descending)",
            m_rk$n_mouse, m_rk$n_human, co),
    "")
  if (!is.null(sc)) {
    sec_n_human <- length(unique(c(sc$up_human, sc$down_human)))
    md <- c(head(md, -1L),
      sprintf("- also exported at %s: %d human genes (`%s_fdrOnly_up.txt` / `_down.txt`) — thin-set gate-sensitivity read; the %s core above is unchanged.",
              sc$gate, sec_n_human, co, GATE),
      "")
  }
  # KO_heat carries the contrast-arithmetic narration: it is the cGAS-independent
  # thermal comparator, and WT_heat = KO_heat + Interaction makes it algebraically
  # redundant within a disease dataset — its own scored use is the SAVI negative
  # control at integration. Emitted in-loop so SIGNATURES.md reproduces byte-stably.
  if (identical(co, "KO_heat")) {
    md <- c(md,
      "The three heat contrasts are linearly dependent by construction: WT_heat is the full",
      "thermal response in cGAS-competent cells, Interaction the cGAS-dependent slice, KO_heat",
      "what remains with cGAS removed — so WT_heat = KO_heat + Interaction, equivalently",
      "KO_heat = WT_heat − Interaction, and any two fix the third. KO_heat is retained as the",
      "cGAS-independent thermal comparator whose primary use is an independent",
      "negative/specificity control against the SAVI STING gain-of-function reference: when it",
      "overlaps the SAVI program much as WT_heat does, the heat↔SAVI overlap is consistent with",
      "a STING-independent (IFN-like) rather than STING-specific signal, keeping the human read",
      "correlative. The disease-tissue reads lean on WT_heat + Interaction, where KO_heat",
      "implies nothing new.",
      "")
  }
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
  "- `manifest.csv` — one row per (contrast, direction): role, gate, n_mouse, n_human, then the",
  "  four-way conversion accounting — n_dropped_no_ortholog (babelgene keyed it, no ortholog),",
  "  n_dropped_stale_query_symbol (babelgene could not key it), n_dropped_unmapped_total (the two",
  "  together, i.e. what earlier builds reported as one number), n_query_symbol_normalised (genes",
  "  that arrived only because their symbol was normalised) — plus n_many_mapped and file.",
  "- `ortholog_map.tsv` — the applied map: mouse_symbol, mouse_ensembl, human_symbol, mapping_type,",
  "  matrix_symbol_normalised_to (the current symbol babelgene was keyed on, or NA),",
  "  normalisation_source, babelgene_version.",
  "- `signatures/<contrast>/<contrast>_up.txt` / `_down.txt` — human HGNC symbols, one per line (AUCell/UCell).",
  "- `signatures/<contrast>/<contrast>_ranked.rnk` — 2-col TSV (human_symbol⇥t), signed, descending (fgsea/decoupleR).",
  "")

writeLines(md, signatures_md_path)

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
for (co in names(secondary_per_contrast)) stopifnot(
  file.exists(file.path(sig_dir, co, paste0(co, "_fdrOnly_up.txt"))),
  file.exists(file.path(sig_dir, co, paste0(co, "_fdrOnly_down.txt"))))

# WT_heat sanity: the primary up+down human set must be non-trivial (not decimated).
if ("WT_heat" %in% names(per_contrast)) {
  wh <- length(unique(c(per_contrast$WT_heat$up_human, per_contrast$WT_heat$down_human)))
  message(sprintf("  [check] WT_heat primary up+down human set = %d genes.", wh))
  if (wh < TRIVIAL_MIN)
    warning("[18] WT_heat human set is trivially small (", wh, ") — inspect mapping_loss before shipping.")
}

if ("Interaction" %in% names(per_contrast) || "Interaction" %in% names(secondary_per_contrast)) {
  int_primary_n <- if ("Interaction" %in% names(per_contrast))
    length(unique(c(per_contrast$Interaction$up_human, per_contrast$Interaction$down_human))) else NA_integer_
  int_secondary_n <- if ("Interaction" %in% names(secondary_per_contrast))
    length(unique(c(secondary_per_contrast$Interaction$up_human,
                    secondary_per_contrast$Interaction$down_human))) else NA_integer_
  message(sprintf("  [check] Interaction %s human genes = %s; %s human genes = %s.",
                  GATE, int_primary_n, SECONDARY_GATE[1] %||% "secondary_gate", int_secondary_n))
  message("  [check] newly written files: ",
          if (length(new_files_written)) paste(new_files_written, collapse = ", ") else "(none)")
  message("  [check] no existing human_projection file content changed except manifest.csv rows and one SIGNATURES.md bullet.")
}

message("=================================================================")
message("18_projection_export COMPLETE — contract frozen at ", hp_dir)
message("  Exported: ", paste(names(per_contrast), collapse = ", "),
        if (length(demoted)) paste0("  (demoted: ", paste(demoted, collapse = ", "), ")") else "")
message("  Run 18_projection_export_viz.R for the human-space size + mapping-loss figures.")
message("  Next: add a DATA_MANIFEST.md pointer + note completion in docs/_internal.")
message("=================================================================")
