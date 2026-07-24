# Mouse → human projection signatures (frozen contract)

Frozen mouse 39 °C Treg-stress signature export from the internal bulk RNA-seq anchor
(GSE329522 iTreg, 2×2 genotype × temperature). Every human compartment (Phases 1–4)
scores against THIS artifact. It is a reshape + ortholog map of the finished DE master —
**no new biology**. Language is correlative: these are the INPUT the human phases test
for presence, never a claim that fever/HIF/STING drives human disease.

## Provenance chain

`03_results/master/master_de_genes.csv` → `03_results/objects/17_signature_sets.rds`
(17_signature_derive.R) → this contract (18_projection_export.R).
- git SHA: `ea09c85a6afd38bba290435a4b432f7d3fed8cc0`
- built_at: 2026-07-03 22:21:41 UTC

## Thresholds / gate applied

- gate: **fdr_logfc** (binary up/down sets)
- de_fdr = 0.05 ; de_logfc = 1
- rank metric: signed **t** (ranked .rnk lists; NEVER logFC)

## Ortholog source + collision policy

- babelgene 22.9, bundled data 2022-09-29 09:40:17 UTC; direction: mouse->human (species='mouse', human=FALSE); min_support = 3.
- one mouse → many human: binary sets UNION; ranked assigns the mouse t to each human ortholog.
- many mouse → one human: ranked keeps MAX |t| per human symbol; binary dedupe by union.
- no human ortholog: dropped (see `ortholog_map.tsv` for every applied edge; unmapped genes absent).

## Exported contrasts

### WT_heat — role: primary
- definition: WT: heat (39 vs 37 °C)
- up:   213 mouse → 199 human (dropped 12, many-mapped 0)
- down: 126 mouse → 94 human (dropped 28, many-mapped 1)
- ranked: 19679 mouse → 12986 human genes (`WT_heat_ranked.rnk`, signed t, descending)

### KO_heat — role: comparator
- definition: cGAS-KO: heat (39 vs 37 °C)
- up:   239 mouse → 218 human (dropped 19, many-mapped 0)
- down: 153 mouse → 113 human (dropped 34, many-mapped 2)
- ranked: 19679 mouse → 12986 human genes (`KO_heat_ranked.rnk`, signed t, descending)

### Interaction — role: comparator
- definition: Heat × genotype interaction (cGAS-dependence of the heat response) (tests cGAS-dependence of the heat response; 1 df, underpowered at n=5 — labelled by what it TESTS, not a result)
- up:   9 mouse → 7 human (dropped 1, many-mapped 0)
- down: 0 mouse → 0 human (dropped 0, many-mapped 0)
- ranked: 19679 mouse → 12986 human genes (`Interaction_ranked.rnk`, signed t, descending)
- also exported at fdr_only: 18 human genes (`Interaction_fdrOnly_up.txt` / `_down.txt`) — thin-set gate-sensitivity read; the fdr_logfc core above is unchanged.

## Files

- `manifest.csv` — one row per (contrast, direction): role, gate, n_mouse, n_human, n_dropped_no_ortholog, n_many_mapped, file.
- `ortholog_map.tsv` — the applied map: mouse_symbol, mouse_ensembl, human_symbol, mapping_type, babelgene_version.
- `signatures/<contrast>/<contrast>_up.txt` / `_down.txt` — human HGNC symbols, one per line (AUCell/UCell).
- `signatures/<contrast>/<contrast>_ranked.rnk` — 2-col TSV (human_symbol⇥t), signed, descending (fgsea/decoupleR).

