#!/usr/bin/env bash
## run_pathway_explorer.sh — pathway-explorer adapter driver (emit universe -> consolidate -> validate -> emit HTML).
## Re-runnable; degrades gracefully on missing masters. Run from project root.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"   # 02_analysis/scripts/pathway_explorer_adapter -> project root
cd "$ROOT"
MST="03_results/master"; UNI="$MST/master_unified.csv"; DE="$MST/master_de_table.csv"; OUT="03_results/interactive"
ADP="02_analysis/scripts/pathway_explorer_adapter"
mkdir -p "$OUT"

## 1) pin the VENDORED copy as the imported package, deterministically, then verify its version.
VENDOR_DIR="$ROOT/01_modules/.ref/pathway-explorer"
VENDOR_SRC="$VENDOR_DIR/src"
VENDOR_FROM="$VENDOR_DIR/VENDORED_FROM.txt"

if [ -d "$VENDOR_SRC/pathway_explorer" ]; then
  ## (a) Shadow any installed copy with the vendored src for every python call below.
  export PYTHONPATH="$VENDOR_SRC:${PYTHONPATH:-}"

  ## (b) Parse the expected version from VENDORED_FROM.txt.
  EXPECTED_VER="$(sed -n 's/^[[:space:]]*__version__[[:space:]]*:[[:space:]]*\([0-9][0-9.]*\).*/\1/p' "$VENDOR_FROM" 2>/dev/null | head -n1)"
  [ -z "$EXPECTED_VER" ] && EXPECTED_VER="$(sed -n 's/^[[:space:]]*pinned tag[[:space:]]*:[[:space:]]*v\?\([0-9][0-9.]*\).*/\1/p' "$VENDOR_FROM" 2>/dev/null | head -n1)"
  if [ -z "$EXPECTED_VER" ]; then
    echo "[adapter] ERROR: could not parse a version from $VENDOR_FROM."; exit 1
  fi

  ## (c) Assert the package python actually resolves is the vendored one AND matches the pin.
  read -r ACT_FILE ACT_VER < <(python -c "import pathway_explorer as p; print(f'ADAPTER_CONFIRMATION:{p.__file__}:{getattr(p, \"__version__\", \"MISSING\")}')" 2>/dev/null | sed -n 's/^ADAPTER_CONFIRMATION:\([^:]*\):\([^:]*\)/\1 \2/p') || true
  if [ -z "${ACT_FILE:-}" ]; then
    echo "[adapter] ERROR: pathway_explorer not importable even with the vendored src on PYTHONPATH ($VENDOR_SRC)."
    exit 1
  fi
  ACT_PKG_DIR="$(cd "$(dirname "$ACT_FILE")" && pwd -P)"
  VENDOR_PKG_DIR="$(cd "$VENDOR_SRC/pathway_explorer" && pwd -P)"
  if [ "$ACT_PKG_DIR" != "$VENDOR_PKG_DIR" ]; then
    echo "[adapter] ERROR: the imported pathway_explorer is NOT the vendored copy."
    echo "[adapter]   expected package dir : $VENDOR_PKG_DIR"
    echo "[adapter]   actual   package dir : $ACT_PKG_DIR  ($ACT_FILE)"
    exit 1
  fi
  if [ "$ACT_VER" != "$EXPECTED_VER" ]; then
    echo "[adapter] ERROR: vendored pathway_explorer version mismatch."
    echo "[adapter]   expected (VENDORED_FROM.txt) : $EXPECTED_VER"
    echo "[adapter]   actual   (package __version__): $ACT_VER"
    exit 1
  fi
  echo "[adapter] using VENDORED pathway_explorer v$ACT_VER from $VENDOR_PKG_DIR (PYTHONPATH-pinned)."
else
  echo "[adapter] WARNING: vendored src not found at $VENDOR_SRC — falling back to install-based import logic."
  python -c "import pathway_explorer" 2>/dev/null || {
    echo "[adapter] pathway_explorer not importable — attempting editable install from 01_modules/.ref/pathway-explorer"
    pip install -e "01_modules/.ref/pathway-explorer" >/dev/null 2>&1
    python -c "import pathway_explorer" 2>/dev/null || {
      echo "[adapter] ERROR: pathway_explorer still not importable. Cannot emit HTML."; exit 1
    }
  }
fi

## 2) emit the set-level explorer universe (custom for STING project)
Rscript "$ADP/emit_explorer_universe.R" || { echo "[adapter] universe emission failed"; exit 1; }

## 3) consolidate the bundle (idempotent; exits 0 even if masters absent)
Rscript "$ADP/consolidate_explorer_bundle.R" || { echo "[adapter] consolidation failed"; exit 1; }

## 4) graceful skip if no unified master
[ -f "$UNI" ] || { echo "[adapter] $UNI absent — nothing to visualize. Exit 0."; exit 0; }

## 5) fail-fast validate the DE bundle
python "$ADP/validate_pathway_explorer_input.py" "$DE" || { echo "[adapter] DE validation FAILED"; exit 1; }

## 6) emit per-contrast HTML + index.html
DE_FLAG=(); [ -f "$DE" ] && DE_FLAG=(--de-data "$DE")

UNIV="$MST/explorer_universe.csv"
UNIV_FLAG=(); [ -f "$UNIV" ] && UNIV_FLAG=(--universe "$UNIV")

if [ -f "$UNIV" ]; then
  HELP_TXT="$(python -m pathway_explorer --help 2>&1)"
  if ! grep -q -- '--universe' <<<"$HELP_TXT"; then
    echo "[adapter] ERROR: $UNIV exists but the resolved pathway_explorer does NOT support --universe."
    exit 1
  fi
  if ! grep -q -- '--force-refit' <<<"$HELP_TXT"; then
    echo "[adapter] ERROR: $UNIV exists but the resolved pathway_explorer lacks --force-refit."
    exit 1
  fi
fi

REF_PARQUET="$OUT/reference_embedding.parquet"
REFIT_FLAG=()
if [ "${EXPLORER_FORCE_REFIT:-auto}" = "1" ]; then
  REFIT_FLAG=(--force-refit)
  echo "[adapter] --force-refit forced via EXPLORER_FORCE_REFIT=1."
elif [ "${EXPLORER_FORCE_REFIT:-auto}" = "auto" ] && [ -f "$UNIV" ]; then
  if [ ! -f "$REF_PARQUET" ] || [ "$UNIV" -nt "$REF_PARQUET" ]; then
    REFIT_FLAG=(--force-refit)
    echo "[adapter] explorer_universe.csv newer than reference_embedding.parquet — adding --force-refit."
  fi
fi

python -m pathway_explorer --data "$UNI" "${DE_FLAG[@]}" "${UNIV_FLAG[@]}" "${REFIT_FLAG[@]}" --all --output "$OUT" --design none
render_rc=$?
if [ "$render_rc" -ne 0 ]; then
  echo "[adapter] ERROR: pathway_explorer render failed (exit $render_rc)."
  exit "$render_rc"
fi
echo "[adapter] HTML written under $OUT:"; ls -1 "$OUT"/*.html 2>/dev/null

## 6a) Post-render assertion: at least 1 HTML must embed basis_contrast_independent: true.
n_bci=$(grep -rlE '"basis_contrast_independent"[[:space:]]*:[[:space:]]*true' "$OUT"/*.html 2>/dev/null | wc -l)
if [ "$n_bci" -lt 1 ]; then
  echo "[adapter] FAIL: no rendered HTML embeds basis_contrast_independent=true."
  exit 1
fi
echo "[adapter] basis_contrast_independent=true confirmed in $n_bci HTML file(s)."

## 7) append per-contrast + index HTML captions to interactive/README.md (idempotent)
python - <<'PYEOF'
import os, glob, pathlib

readme_path = pathlib.Path("03_results/interactive/README.md")
existing = readme_path.read_text() if readme_path.exists() else ""

# per-contrast caption template
per_contrast_caption = """## pathway_explorer_<contrast>.html
**Per-contrast interactive pathway-explorer dashboard (UMAP of GSEA pathways + CollecTRI TF + PROGENy activities, NES-colored, FDR-sliderable, with the t-ranked running-sum panel): this dashboard enables interactive cross-entity similarity mapping and detailed exploration of regulatory cascades.**
| | |
|---|---|
| Script   | `02_analysis/scripts/pathway_explorer_adapter/run_pathway_explorer.sh` |
| Function | `generate_all_dashboards()` (pathway_explorer) |
| Config   | `paths.interactive = 03_results/interactive/; paths.master = 03_results/master/` |
| Input    | `03_results/master/master_unified.csv; 03_results/master/master_de_table.csv` |

"""

index_caption = """## index.html
**Landing page linking every per-contrast pathway-explorer dashboard: the entry point for interactive exploration of MSigDB/custom pathways, TF activities, and PROGENy activities.**
| | |
|---|---|
| Script   | `02_analysis/scripts/pathway_explorer_adapter/run_pathway_explorer.sh` |
| Function | `generate_index_page()` (pathway_explorer) |
| Config   | `paths.interactive = 03_results/interactive/; paths.master = 03_results/master/` |
| Input    | `03_results/master/master_unified.csv` |

"""

to_append = ""
if "## pathway_explorer_<contrast>.html" not in existing:
    to_append += per_contrast_caption
if "## index.html" not in existing:
    to_append += index_caption

if to_append:
    with open(readme_path, "a") as fh:
        fh.write(to_append)
    print("[adapter] README.md HTML captions written.")
else:
    print("[adapter] README.md HTML captions already present (idempotent).")
PYEOF

echo "[adapter] run_pathway_explorer.sh complete."
