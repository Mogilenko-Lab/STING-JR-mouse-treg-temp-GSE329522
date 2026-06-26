<!-- BEGIN SCIAGENT:ROLES v1 hash=e3769f72bc396457e4505b9d87a85dc7b87a8748 -->
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
- `figure-style` — (base)
- `scrna-pipeline-conventions` — (base)
- `cellranger-multi-to-anndata` — (base)
- `scrna-cxg-host` — (base)
- `consensus-nmf-multirun` — (base)
- `skill-creator` — (base)
- `reasoning-trace` — (base)

## Skills (inherited via requires:)
- `tobias-footprint-bindetect` — (via `tf-footprint-differential-analysis`)
- `hint-atac-differential-footprint` — (via `tf-footprint-differential-analysis`)
- `signac-footprint-visualization` — (via `tf-footprint-differential-analysis`)

## Sub-agents (effective, Claude-only)
- `docs-librarian` — (pathway-signature) [shadows base]
- `bio-interpreter` — (pathway-signature) [shadows base]
- `insight-explorer` — (pathway-signature) [shadows base]
- `captions` — (pathway-signature) [shadows base]
- `figure-audit` — (base)
- `doc-curator` — (pathway-signature) [shadows base]
- `code-reviewer` — (pathway-signature) [shadows base]
- `handoff` — (pathway-signature) [shadows base]

## Slash commands (effective, Claude-only)
- `/commit` — (pathway-signature) [shadows base]
<!-- END SCIAGENT:ROLES -->

<!-- BEGIN SCIAGENT:CRAFT v1 hash=1e98ccc98cac4ef01baff5fbbcea627f4f3e2654 -->
# Craft standards

Toolkit-managed standing conventions for this repo — do not hand-edit.

- **Figures** — legible both shrunk in a journal column and projected to the back of a room: bigger, fewer, bolder (one legible tier, base >= 14pt; emit a vector PDF + raster PNG from one plot object). Style only via the project theme entry point (no inline `theme()`/`ggsave(width=)`/raw hex); cap to top-N; never truncate axis labels; disambiguate glyphs; prefer the residualized channel.
- **Results** — every artifact under `03_results/<stage>/{figures,tables}/` with `by_contrast/<c>/` + `_overview/`; a figure's source table is its same-stem neighbor. Compute never plots; viz never computes.
- **README** — a task is unfinished until the sibling `README.md` captions every `03_results/` file you create/edit/delete, including *how to read* it (glyphs, sign convention, Δρ, claim tier).
- **Planning** — plans in `docs/_internal/plans/{date-slug}/` as `00_INDEX.md` + `NN_<slug>.md`; one phase == one script == one implementer (~35% of context); review every 3 phases that code runs and produces its artifacts.
- **Reproducibility** — no ephemeral scripts: every `03_results/` artifact reproducible from a committed `02_analysis/scripts/NN_*`; log non-trivial decisions to `docs/_internal/reasoning/` before proceeding (`_scratch/` is the only sanctioned throwaway zone).
<!-- END SCIAGENT:CRAFT -->
