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
- many mouse → one human: ranked keeps MAX |t| per human symbol, except that an edge from
  the original query outranks one recovered by symbol normalisation regardless of |t|;
  binary dedupe by union.
- stale mouse symbol: babelgene keys on CURRENT MGI symbols, and this matrix was quantified against GENCODE vM25, so 2341 of its 19679 symbols are no longer current. babelgene knows some of them and not others, idiosyncratically, so the ones it could not key at all are re-asked under their current symbol via org.Mm.eg.db 3.22.0 and the edges mapped back — `mouse_symbol` below always stays the symbol the data carries.
- 144 genes arrived that way, each previously counted as having no human ortholog, which was the wrong label for a vocabulary loss.
- the recovery is strictly additive: 155 edges added, none removed, and a human symbol the map already carried keeps the ranked metric it had. Every candidate, accepted or withheld, is a row in `03_results/11_projection/tables/_overview/query_normalisation_ledger.csv`.
- no human ortholog: dropped. This now means what it says — babelgene could key the symbol
  and returned no ortholog at this min_support. A symbol it could not key at all is counted
  separately as `n_dropped_stale_query_symbol`, because that is a statement about vocabulary
  rather than about orthology, and `n_dropped_unmapped_total` keeps the two together for
  comparison with earlier builds. See `ortholog_map.tsv` for every applied edge.

## Exported contrasts

### WT_heat — role: primary
- definition: WT: heat (39 vs 37 °C)
- up:   213 mouse → 202 human (no ortholog 8, symbol not keyable 1, recovered by normalisation 3, many-mapped 0)
- down: 126 mouse → 96 human (no ortholog 19, symbol not keyable 7, recovered by normalisation 2, many-mapped 1)
- ranked: 19679 mouse → 13128 human genes (`WT_heat_ranked.rnk`, signed t, descending)

### KO_heat — role: comparator
- definition: cGAS-KO: heat (39 vs 37 °C)
- up:   239 mouse → 221 human (no ortholog 11, symbol not keyable 5, recovered by normalisation 3, many-mapped 0)
- down: 153 mouse → 115 human (no ortholog 24, symbol not keyable 8, recovered by normalisation 2, many-mapped 2)
- ranked: 19679 mouse → 13128 human genes (`KO_heat_ranked.rnk`, signed t, descending)

The three heat contrasts are linearly dependent by construction: WT_heat is the full
thermal response in cGAS-competent cells, Interaction the cGAS-dependent slice, KO_heat
what remains with cGAS removed — so WT_heat = KO_heat + Interaction, equivalently
KO_heat = WT_heat − Interaction, and any two fix the third. KO_heat is retained as the
cGAS-independent thermal comparator whose primary use is an independent
negative/specificity control against the SAVI STING gain-of-function reference: when it
overlaps the SAVI program much as WT_heat does, the heat↔SAVI overlap is consistent with
a STING-independent (IFN-like) rather than STING-specific signal, keeping the human read
correlative. The disease-tissue reads lean on WT_heat + Interaction, where KO_heat
implies nothing new.

### Interaction — role: comparator
- definition: Heat × genotype interaction (cGAS-dependence of the heat response) (tests cGAS-dependence of the heat response; 1 df, underpowered at n=5 — labelled by what it TESTS, not a result)
- up:   9 mouse → 7 human (no ortholog 1, symbol not keyable 0, recovered by normalisation 0, many-mapped 0)
- down: 0 mouse → 0 human (no ortholog 0, symbol not keyable 0, recovered by normalisation 0, many-mapped 0)
- ranked: 19679 mouse → 13128 human genes (`Interaction_ranked.rnk`, signed t, descending)
- also exported at fdr_only: 19 human genes (`Interaction_fdrOnly_up.txt` / `_down.txt`) — thin-set gate-sensitivity read; the fdr_logfc core above is unchanged.

## Files

- `manifest.csv` — one row per (contrast, direction): role, gate, n_mouse, n_human, then the
  four-way conversion accounting — n_dropped_no_ortholog (babelgene keyed it, no ortholog),
  n_dropped_stale_query_symbol (babelgene could not key it), n_dropped_unmapped_total (the two
  together, i.e. what earlier builds reported as one number), n_query_symbol_normalised (genes
  that arrived only because their symbol was normalised) — plus n_many_mapped and file.
- `ortholog_map.tsv` — the applied map: mouse_symbol, mouse_ensembl, human_symbol, mapping_type,
  matrix_symbol_normalised_to (the current symbol babelgene was keyed on, or NA),
  normalisation_source, babelgene_version.
- `signatures/<contrast>/<contrast>_up.txt` / `_down.txt` — human HGNC symbols, one per line (AUCell/UCell).
- `signatures/<contrast>/<contrast>_ranked.rnk` — 2-col TSV (human_symbol⇥t), signed, descending (fgsea/decoupleR).

