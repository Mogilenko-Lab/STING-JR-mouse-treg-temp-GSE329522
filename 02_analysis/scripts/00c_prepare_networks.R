#!/usr/bin/env Rscript
# 00c_prepare_networks.R - Build decoupleR regulatory networks LOCALLY (CollecTRI, DoRothEA, PROGENy)
# Project: STING-cGAS-GSE329522
# Phase: 0 (infrastructure / pinned reference networks)
# Dependencies: 02_analysis/config/config.R; progeny, dorothea (Bioconductor data pkgs); babelgene
# Outputs:
#   03_results/objects/net_collectri_mouse.rds       (primary TF net; source/target/mor)
#   03_results/objects/net_dorothea_mouse_ABC.rds    (forensics + robustness net; source/target/mor/confidence)
#   03_results/objects/net_progeny_mouse.rds          (Phase 4 pathway net; source/target/weight)
#   03_results/objects/networks_prep_log.txt          (provenance + coverage audit)
#
# WHY LOCAL (not decoupleR::get_collectri / get_progeny):
#   decoupleR 2.16.0 + OmnipathR 3.18.4 (latest Bioconductor-release versions) cannot fetch CollecTRI
#   or PROGENy: OmnipathR's web wrapper falls back to a static table whose loader errors
#   ("argument is of length zero"; OmnipathR::collectri() dies on an internal ncbi_tax_id join).
#   omnipathdb.org is reachable (HTTP 200), so this is a package-version bug, not a network outage,
#   and there is no newer release to upgrade into. We therefore build the SAME networks from their
#   authoritative primary sources and pin them as RDS. This is more reproducible (no live-API
#   dependency) and provenance-transparent. Re-point Phases 3/4 to these cached RDS.

source("02_analysis/config/config.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(babelgene) })

LOG <- file.path(DIR_OBJECTS, "networks_prep_log.txt")
con <- file(LOG, open = "wt")
say <- function(...) { msg <- sprintf(...); cat(msg, "\n"); writeLines(msg, con) }

dataset_symbols <- unique(read.csv(FILE_CPM, check.names = FALSE)$gene_name)

# Standard human->mouse symbol fallback (matches the toolkit MitoPathways loader convention)
mouse_titlecase <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))

say("== 00c network prep | %s ==", as.character(Sys.time()))
say("decoupleR %s | OmnipathR %s | progeny %s | dorothea %s | babelgene %s",
    packageVersion("decoupleR"), packageVersion("OmnipathR"),
    packageVersion("progeny"), packageVersion("dorothea"), packageVersion("babelgene"))

# ===========================================================================
# 1) PROGENy (Mouse, top=500) via the bundled progeny model (NO OmniPath)
# ===========================================================================
prog_mat <- progeny::getModel("Mouse", top = 500)                  # genes x 14 pathways
net_progeny <- as.data.frame(prog_mat) %>%
  tibble::rownames_to_column("target") %>%
  tidyr::pivot_longer(-target, names_to = "source", values_to = "weight") %>%
  dplyr::filter(weight != 0) %>%
  dplyr::select(source, target, weight)
attr(net_progeny, "provenance") <- list(
  source = "progeny::getModel('Mouse', top=500)", pkg = paste0("progeny ", packageVersion("progeny")),
  note = "Bypasses OmniPath; identical model matrix decoupleR::get_progeny would return.")
saveRDS(net_progeny, file.path(DIR_OBJECTS, "net_progeny_mouse.rds"))
say("\n[PROGENy] %d edges | %d pathways | %d unique targets | %d in dataset",
    nrow(net_progeny), dplyr::n_distinct(net_progeny$source),
    dplyr::n_distinct(net_progeny$target), sum(unique(net_progeny$target) %in% dataset_symbols))
stopifnot(dplyr::n_distinct(net_progeny$source) == 14)

# ===========================================================================
# 2) DoRothEA mouse (confidence A,B,C) via bundled dorothea_mm (NO OmniPath)
#    This is the network class Phylo used (DoRothEA A/B/C) -> primary forensics net.
# ===========================================================================
data("dorothea_mm", package = "dorothea")
net_dorothea <- dorothea_mm %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::transmute(source = tf, target = target, mor = mor, confidence = confidence) %>%
  dplyr::distinct()
attr(net_dorothea, "provenance") <- list(
  source = "dorothea::dorothea_mm (confidence A,B,C)", pkg = paste0("dorothea ", packageVersion("dorothea")),
  note = "Native mouse regulons; the network class the Phylo agent used for its HIF1a ranking.")
saveRDS(net_dorothea, file.path(DIR_OBJECTS, "net_dorothea_mouse_ABC.rds"))
say("[DoRothEA A/B/C] %d edges | %d TFs | %d unique targets | %d in dataset",
    nrow(net_dorothea), dplyr::n_distinct(net_dorothea$source),
    dplyr::n_distinct(net_dorothea$target), sum(unique(net_dorothea$target) %in% dataset_symbols))

# ===========================================================================
# 3) CollecTRI (house-default primary TF net): Zenodo human regulons -> mouse
#    Source: Mueller-Dott et al., CollecTRI, Zenodo record 8192729.
#    Human symbols mapped to mouse via babelgene orthology + title-case fallback.
# ===========================================================================
CTRI_CSV <- file.path(PROJECT_ROOT, "00_data/references/networks/CollecTRI_regulons_human.csv")
stopifnot(file.exists(CTRI_CSV))   # downloaded in the network-prep step
ctri_h <- read.csv(CTRI_CSV, stringsAsFactors = FALSE)   # source,target,weight(+/-1),resources,references,sign_decision

genes_h <- unique(c(ctri_h$source, ctri_h$target))
orth <- babelgene::orthologs(genes = genes_h, species = "mouse", human = TRUE)
map <- setNames(orth$symbol, orth$human_symbol)          # human -> mouse (babelgene)
to_mouse <- function(x) { m <- unname(map[x]); ifelse(is.na(m), mouse_titlecase(x), m) }

net_collectri <- ctri_h %>%
  dplyr::transmute(source = to_mouse(source), target = to_mouse(target), mor = weight) %>%
  dplyr::filter(!is.na(source), !is.na(target), nzchar(source), nzchar(target)) %>%
  dplyr::distinct(source, target, .keep_all = TRUE)
n_babel <- sum(genes_h %in% names(map)); n_fallback <- length(genes_h) - n_babel
attr(net_collectri, "provenance") <- list(
  source = "CollecTRI (Zenodo 8192729, human) -> mouse",
  doi = "10.1093/bioadv/vbad026 (Mueller-Dott et al., CollecTRI)",
  mapping = sprintf("babelgene %s orthology (%d/%d genes) + title-case fallback (%d)",
                    packageVersion("babelgene"), n_babel, length(genes_h), n_fallback),
  note = "Local reconstruction because decoupleR::get_collectri is broken at installed versions.")
saveRDS(net_collectri, file.path(DIR_OBJECTS, "net_collectri_mouse.rds"))
say("[CollecTRI->mouse] %d edges | %d TFs | %d unique targets | %d in dataset",
    nrow(net_collectri), dplyr::n_distinct(net_collectri$source),
    dplyr::n_distinct(net_collectri$target), sum(unique(net_collectri$target) %in% dataset_symbols))
say("  human->mouse: %d via babelgene, %d via title-case fallback (of %d unique symbols)",
    n_babel, n_fallback, length(genes_h))

# ---- sanity: key TFs present as sources; ISG targets reachable ----
key_present <- intersect(KEY_TFS, unique(net_collectri$source))
say("\n[sanity] KEY_TFS present in CollecTRI source: %s", paste(key_present, collapse = ", "))
say("[sanity] Hif1a in DoRothEA source: %s | Hif1a in CollecTRI source: %s",
    "Hif1a" %in% net_dorothea$source, "Hif1a" %in% net_collectri$source)

close(con)
cat("\n[OK] wrote net_collectri_mouse.rds, net_dorothea_mouse_ABC.rds, net_progeny_mouse.rds + log\n")
cat("Log:", LOG, "\n")
