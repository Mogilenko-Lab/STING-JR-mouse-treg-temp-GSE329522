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
## org.Mm.eg.db) plus the GATOM-viz network renderers and the set-membership
## renderers used by the stage 12_hsr_decomp viz (eulerr / ggvenn / ComplexUpset).
##
## Idempotent by construction: every package is probed with requireNamespace() first,
## anything already loadable is left alone, and the run ends with an explicit
## INSTALLED / SKIPPED / FAILED ledger, making a re-run self-documenting.

suppressPackageStartupMessages(library(BiocManager))

## WORKAROUND: /tmp is mounted noexec in this devcontainer, which makes R's
## internal file_test("-x", "configure") check fail for any package shipping a
## configure script (stringfish, qs2, eulerr, RcppArmadillo, ...). R fixes
## tempdir() at STARTUP, so Sys.setenv(TMPDIR=) inside a running session is too
## late — install.packages() still unpacks under the old, noexec tempdir. The
## only reliable fix is to re-exec this script in a fresh R process that already
## has an exec-capable TMPDIR in its environment. The guard env var makes that
## re-exec happen at most once.
EXEC_TMPDIR <- path.expand("~/tmp_exec")
.tempdir_is_exec <- function() {
  probe <- file.path(tempdir(), "sciagent_exec_probe.sh")
  ok <- tryCatch({
    writeLines("#!/bin/sh\nexit 0", probe)
    Sys.chmod(probe, "0755")
    file_test("-x", probe) && system2(probe, stdout = FALSE, stderr = FALSE) == 0L
  }, error = function(e) FALSE, warning = function(w) FALSE)
  unlink(probe)
  isTRUE(ok)
}
if (!.tempdir_is_exec() && !nzchar(Sys.getenv("SCIAGENT_INSTALL_DEPS_REEXEC"))) {
  dir.create(EXEC_TMPDIR, recursive = TRUE, showWarnings = FALSE)
  message("[install_deps] tempdir() is on a noexec mount (", tempdir(), ").")
  message("[install_deps] Re-executing with TMPDIR=", EXEC_TMPDIR, " so source builds can run.")
  script <- sub("^--file=", "",
                grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                    env = c(paste0("TMPDIR=", EXEC_TMPDIR),
                            "SCIAGENT_INSTALL_DEPS_REEXEC=1"))
  quit(save = "no", status = status)
}

## Packages required by the standard-sweep arms not guaranteed in the base stack.
REQUIRED_PKGS <- c(
  ## GATOM metabolic modules (scripts 10/14) — gatom + the MWCS solver + graph deps
  "gatom", "mwcsr", "igraph", "BioNet", "sna",
  ## GATOM network viz (script 14) — static network PDFs
  "ggraph", "tidygraph",
  ## CoReSh chunk IO + entrez mapping (scripts 07/08; gated on the compendium)
  "qs2", "org.Mm.eg.db", "BiocParallel",
  ## Set-membership renderers for stage 12_hsr_decomp (script 19 viz):
  ##   eulerr       - area-proportional Euler fit + its own diagError/stress residual
  ##   ggvenn       - conventional fixed-layout Venn, pure ggplot2 (no sf/GDAL stack)
  ##   ComplexUpset - UpSet plot as a ggplot2/patchwork object, so it composes with
  ##                  project_theme()/save_figure() like every other panel here
  "eulerr", "ggvenn", "ComplexUpset"
)

## PINNED SOURCE TARBALLS — packages whose CRAN HEAD cannot be built here.
##   eulerr >= 8.0.0 moved its optimiser to Rust and its configure step aborts with
##   [RUST NOT FOUND]; no rustc/cargo is provisioned in this container. 7.0.2 is the
##   last pure C++/Rcpp release and exposes the same euler() API and the same
##   stress / diagError fit diagnostics the membership panels report. Its extra
##   hard dependency GenSA is listed so a clean container resolves it first.
PINNED_SRC <- c(
  eulerr = "https://cloud.r-project.org/src/contrib/Archive/eulerr/eulerr_7.0.2.tar.gz"
)
PINNED_DEPS <- list(eulerr = c("GenSA", "polyclip", "polylabelr", "RcppArmadillo"))

## ---------------------------------------------------------------------------
## Probe -> install only what is missing -> re-probe -> print the ledger.
## ---------------------------------------------------------------------------
have <- function(p) vapply(p, requireNamespace, logical(1), quietly = TRUE)

present_before <- REQUIRED_PKGS[have(REQUIRED_PKGS)]
missing_pkgs   <- setdiff(REQUIRED_PKGS, present_before)

installed_now <- character(0)
still_missing <- character(0)

if (length(missing_pkgs) == 0L) {
  message("[install_deps] Nothing to do; every declared package is already loadable.")
} else {
  message("[install_deps] Installing missing packages: ",
          paste(missing_pkgs, collapse = ", "))

  ## Pinned-source packages first (with their own hard deps), then the rest via
  ## BiocManager so CRAN and Bioconductor names resolve from one call.
  pinned_missing <- intersect(missing_pkgs, names(PINNED_SRC))
  for (p in pinned_missing) {
    dep_missing <- PINNED_DEPS[[p]][!have(PINNED_DEPS[[p]])]
    if (length(dep_missing)) {
      message("[install_deps]   ", p, ": resolving pinned deps ",
              paste(dep_missing, collapse = ", "))
      BiocManager::install(dep_missing, ask = FALSE, update = FALSE)
    }
    message("[install_deps]   ", p, ": installing pinned source ", PINNED_SRC[[p]])
    tryCatch(utils::install.packages(PINNED_SRC[[p]], repos = NULL, type = "source"),
             error = function(e) message("[install_deps]   ", p, ": ", conditionMessage(e)),
             warning = function(w) message("[install_deps]   ", p, ": ", conditionMessage(w)))
  }

  rest <- setdiff(missing_pkgs, pinned_missing)
  if (length(rest)) BiocManager::install(rest, ask = FALSE, update = FALSE)

  ok            <- have(missing_pkgs)
  installed_now <- missing_pkgs[ok]
  still_missing <- missing_pkgs[!ok]
}

## Ledger: what this run changed vs. what it left alone.
report <- function(label, pkgs) {
  if (length(pkgs) == 0L) {
    message(sprintf("[install_deps] %-9s (0): -", label))
  } else {
    vers <- vapply(pkgs, function(p)
      tryCatch(as.character(utils::packageVersion(p)), error = function(e) "?"),
      character(1))
    message(sprintf("[install_deps] %-9s (%d): %s", label, length(pkgs),
                    paste(sprintf("%s %s", pkgs, vers), collapse = ", ")))
  }
}
report("INSTALLED", installed_now)
report("SKIPPED",   present_before)

if (length(still_missing)) {
  message(sprintf("[install_deps] FAILED    (%d): %s", length(still_missing),
                  paste(still_missing, collapse = ", ")))
  warning("[install_deps] These packages could not be installed: ",
          paste(still_missing, collapse = ", "))
} else {
  message("[install_deps] FAILED    (0): -")
  message("[install_deps] Done. All declared packages are available.")
}
