## install_deps.R — install any missing R / Bioconductor packages for the
## STING-cGAS-GSE329522 standard-sweep backbone (04-16).
##
## Usage (once, from the project root):
##   Rscript 02_analysis/helpers/install_deps.R
##
## Adapted from 01_modules/.ref/14839-DM/02_analysis/helpers/install_deps.R.
## Logic: BiocManager::install() for everything (handles CRAN + Bioconductor,
## respects the Bioc release tied to R 4.5). Skips anything already loadable so
## re-runs are cheap. The base DE/GSEA/decoupleR stack is already provisioned in
## this container; this list covers the arms whose deps were still missing
## (GATOM: gatom/mwcsr + solver/network deps; CoReSh IO: qs2; entrez map:
## org.Mm.eg.db) plus the GATOM-viz network renderers.

suppressPackageStartupMessages(library(BiocManager))

## WORKAROUND: /tmp is mounted noexec in this devcontainer, which makes R's
## internal file_test("-x", "configure") check fail for packages shipping a
## configure script (e.g. stringfish, qs2). Redirect tempdir under $HOME.
if (!nzchar(Sys.getenv("TMPDIR")) || startsWith(normalizePath(Sys.getenv("TMPDIR", "/tmp")), "/tmp")) {
  exec_tmpdir <- path.expand("~/tmp_exec")
  dir.create(exec_tmpdir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TMPDIR = exec_tmpdir)
}

## Packages required by the standard-sweep arms not guaranteed in the base stack.
REQUIRED_PKGS <- c(
  ## GATOM metabolic modules (scripts 10/14) — gatom + the MWCS solver + graph deps
  "gatom", "mwcsr", "igraph", "BioNet", "sna",
  ## GATOM network viz (script 14) — static network PDFs
  "ggraph", "tidygraph",
  ## CoReSh chunk IO + entrez mapping (scripts 07/08; gated on the compendium)
  "qs2", "org.Mm.eg.db", "BiocParallel"
)

missing_pkgs <- REQUIRED_PKGS[
  !vapply(REQUIRED_PKGS, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) == 0L) {
  message("[install_deps] All packages already installed.")
} else {
  message("[install_deps] Installing missing packages: ",
          paste(missing_pkgs, collapse = ", "))
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  still_missing <- missing_pkgs[
    !vapply(missing_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(still_missing))
    warning("[install_deps] These packages could not be installed: ",
            paste(still_missing, collapse = ", "))
  else
    message("[install_deps] Done. All packages now available.")
}
