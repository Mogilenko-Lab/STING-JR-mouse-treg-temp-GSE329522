# Project-local R library — `01_modules/Rlib/`

The container's system library at `/usr/local/lib/R/library` stays as the image built
it. When a script needs a build of a package that differs from the image's, that build
goes into `01_modules/Rlib/` and the script opts in with `use_project_rlib()` from
`02_analysis/config/config.R`. The path is configured once, at `paths.r_library` in
`analysis_config.yaml`, so no script hardcodes it.

Sourcing `config.R` does not touch `.libPaths()`. A script that never calls
`use_project_rlib()` sees the system library exactly as before, so one pinned package
cannot leak into the rest of the pipeline.

The directory holds compiled binaries and is gitignored. What is tracked is this file
and the `require_pkgs` argument in each script that depends on a pinned build, so a
missing install stops the run instead of silently resolving to the system version.

## Opting in

```r
source("02_analysis/config/config.R")
use_project_rlib(require_pkgs = "GOSemSim")   # stops if GOSemSim is not in the project library
library(GOSemSim)
```

`PROJECT_RLIB` overrides the configured path for a single run, which is how a stage is
re-run against the system build for comparison:

```bash
PROJECT_RLIB=none Rscript 02_analysis/scripts/20_semantic_decomposition.R   # system library only
PROJECT_RLIB=/some/other/lib Rscript ...                                    # a different library
```

Setting `PROJECT_RLIB` also suspends the `require_pkgs` check, since an explicit
override is a deliberate choice.

## Contents

### GOSemSim 2.39.2

Installed from the vendored source at `../../../01_modules/.ref/GOSemSim` (umbrella
path `01_modules/.ref/GOSemSim`, upstream `YuLab-SMU/GOSemSim`).

```bash
R CMD INSTALL -l 01_modules/Rlib /workspaces/STING-JR/01_modules/.ref/GOSemSim
```

This is the Bioconductor **devel** branch, where the fix (upstream issue #51) landed
first. A release-branch build carrying the same fix would be numbered differently and
is what would normally be pinned; the behavioural guard below is what makes that swap a
one-line config change rather than an edit to an assertion.

The system library carries 2.36.0. In that build, `GOSemSim:::getSV` keys the Wang
edge-weight table on `is_a` / `part_of` and matches it against the relationship column
of the packaged `gotbl` without normalising the spelling, while `gotbl` writes those
relationships as `isa` and `part of`. Every edge fails both matches and falls through to
the `other` weight, so the whole DAG is traversed at a uniform 0.7 and the is_a 0.8 /
part_of 0.6 distinction the Wang measure is built on never reaches the calculation. On
Biological Process the effect is total: all 66,947 edges take 0.7, and none receives 0.8
or 0.6.

| relationship as `gotbl` spells it | BP edges |
|---|---|
| `isa` | 53,016 |
| `part of` | 5,191 |
| `regulates` | 3,216 |
| `positively regulates` | 2,756 |
| `negatively regulates` | 2,768 |

Count `gotbl`, not GO.db's `GOBPPARENTS`. `getSV` walks `gotbl`, and `GOBPPARENTS` omits
the root edges, so it reports 57,462 for the same graph.

Consumed by `02_analysis/scripts/20_semantic_decomposition.R`, which does not trust this
version string. It recovers the effective is_a weight from the measure itself: for a BP
term whose only `gotbl` parent is an is_a edge to the BP root, Wang similarity to that
root is `(1 + w) / (2 + w)`, so `w = (2s - 1) / (1 - s)`. That reads 0.8011 on 2.39.2 and
0.7027 on 2.36.0, and the stage stops unless it matches `semantic.expect_isa_weight`.
The measured value is written to `semantic_provenance.csv` and is part of the term-matrix
cache key. Upgrading GOSemSim therefore cannot silently change the numbers, and a
regressed build cannot pass by carrying the right version.

The size of the difference between the two builds is measured by
`02_analysis/scripts/20b_semantic_engine_validation.R`, which writes
`03_results/13_semantic_decomp/tables/semantic_engine_validation.csv`.
