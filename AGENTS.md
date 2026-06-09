<!-- BEGIN SCIAGENT:ROLES v1 hash=20c3f484422f0f7e278a1dde1b8e9776e7e599fb -->
# Active roles

Stack (in order, last-wins on name collisions):
1. **base** — `roles/base.yaml` — Default bioinformatics analysis role with full agent suite
2. **pathway-signature** — `roles/pathway-signature.yaml` — Pathway/TF/signature functional interpretation — GSEA + decoupleR + CoReSh + pathway-explorer  *(overlay)*

## Skills (effective)
- `anndata` — (pathway-signature) [shadows base]
- `scanpy` — (pathway-signature) [shadows base]
- `single-cell-rna-qc` — (base)
- `anndatar-seurat-scanpy-conversion` — (pathway-signature) [shadows base]
- `bulk-rnaseq-gsea` — (pathway-signature) [shadows base]
- `bulk-rnaseq-activity-inference` — (pathway-signature) [shadows base]
- `bulk-rnaseq-pathway-explorer` — (pathway-signature) [shadows base]
- `gatom-metabolomic-predictions` — (pathway-signature) [shadows base]
- `coresh-signature-search` — (pathway-signature) [shadows base]
- `starsolo-spliced-unspliced` — (base)
- `rna-velocity-trajectory` — (base)
- `genenmf-metaprogram-discovery` — (base)
- `pycistopic-atac-topic-modeling` — (base)
- `crescendo-scatac-cre-analysis` — (base)
- `chromvar-motif-accessibility` — (base)
- `tf-footprint-differential-analysis` — (base)
- `scenic-grn-inference` — (base)
- `scvi-framework` — (base)
- `scvi-basic` — (base)
- `scvi-scanvi` — (base)
- `scvi-multivi` — (base)
- `scvi-peakvi` — (base)
- `scvi-mrvi` — (base)
- `scvi-contrastivevi` — (base)
- `scvi-linearscvi` — (base)
- `scvi-lda` — (base)
- `scvi-hub-models` — (base)
- `scvi-scarches-reference-mapping` — (base)
- `treearches-hierarchy-learning` — (base)
- `scglue-unpaired-multiomics-integration` — (base)
- `scrna-pipeline-conventions` — (base)
- `cellranger-multi-to-anndata` — (base)
- `scrna-cxg-host` — (base)
- `shinymultiome-uio-host` — (base)
- `consensus-nmf-multirun` — (base)
- `skill-creator` — (base)

## Skills (inherited via requires:)
- `tobias-footprint-bindetect` — (via `tf-footprint-differential-analysis`)
- `hint-atac-differential-footprint` — (via `tf-footprint-differential-analysis`)
- `signac-footprint-visualization` — (via `tf-footprint-differential-analysis`)

## Sub-agents (effective, Claude-only)
- `docs-librarian` — (pathway-signature) [shadows base]
- `bio-interpreter` — (pathway-signature) [shadows base]
- `insight-explorer` — (pathway-signature) [shadows base]
- `captions` — (pathway-signature) [shadows base]
- `doc-curator` — (pathway-signature) [shadows base]
- `code-reviewer` — (pathway-signature) [shadows base]
- `handoff` — (pathway-signature) [shadows base]

## Slash commands (effective, Claude-only)
- `/commit` — (pathway-signature) [shadows base]

---

# Analysis conventions (STING-cGAS-GSE329522)

## Figure & conclusion discipline
- **One claim = one dedicated, captioned figure.** No conclusions that live only in chat/handoff — every claim we produce must be persisted as a figure + written caption.
- **Multi-panel ONLY when sub-panels share the same axis and are directly comparable.** Never stack distinct statements into one panel — that makes them impossible to audit/disambiguate (the original `fig3b_hif1a_robustness.pdf` was the anti-pattern: network-swap ranks + target decomposition + signature comparison crammed onto one illegible panel).
- **Every figure documents how it was generated** (which inputs, which computation) — auditable, never a black box.
- **Each `03_results/<phase>/` carries a `README.md`** captioning, laconically: (1) the STATEMENT the artifacts make, and (2) the MECHANISM behind that statement.
- **Stop to explain surprising mechanics** rather than reporting a number and moving on (e.g. a transcription factor that ranks differently across inference methods/networks deserves its own figure explaining the math, not a footnote).

## Normalize then visualize (compute / viz split)
- Each analysis phase = a COMPUTE script `NN_<name>.R` (all statistics; writes checkpoints to `03_results/objects/*.rds`, master tables to `03_results/master/*.csv`, and plot-ready tidy tables to `03_results/<stage>/tables/*.csv`; contains **no** `ggplot`/`ggsave`) **plus** a VIZ script `NN_<name>_viz.R` (reads those normalized tables and renders figures; contains **no** statistical computation — no `lmFit`/`eBayes`/`run_ulm`/`p.adjust`/`prcomp`/…).
- Viz must run standalone after compute and must never recompute statistics.

## decoupleR networks are pre-built locally (OmniPath is broken here)
`decoupleR::get_collectri()` / `get_progeny()` FAIL in this environment — do NOT call them. Use the cached RDS built by `02_analysis/scripts/00c_prepare_networks.R`: `03_results/objects/{net_collectri_mouse,net_dorothea_mouse_ABC,net_progeny_mouse}.rds`, then `run_ulm(.mor="mor")` / `run_mlm(.mor="weight")` directly. Full root-cause + the generic local-build recipe live in the skill: `01_modules/SciAgent-toolkit/skills/bulk-rnaseq-activity-inference/references/known-issues.md`. Env snapshot: `02_analysis/config/env/`.
<!-- END SCIAGENT:ROLES -->
