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
<!-- END SCIAGENT:ROLES -->
