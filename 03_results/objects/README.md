# objects: reusable state, and the three tables that travel with it

`03_results/objects/` holds recomputable checkpoints rather than deliverables, so the `.rds`
files here are not tracked. Three CSVs are, because they are accounting rather than state and a
reader has to be able to check them without re-running anything.

## geneset_manifest.csv

The authoritative index of every gene-set collection available to `05_gsea_msigdb_run.R`,
`06_gsea_custom_run.R` and the viz arm. Written by `02_analysis/scripts/04_gsea_set_prep.R`.
One row per collection: `database`, `type`, `n_sets`, `source`, `cache_path`, then the
vocabulary ledger rolled up over the collection's unique reference symbols.

**How to read the ledger columns.** They close, exactly:

```
n_unique_set_genes == n_matched + n_matched_via_alias + n_alias_flagged_for_review +
                      n_alias_rejected_ambiguous + n_expression_filtered +
                      n_below_detection + n_absent_from_reference
```

- `n_matched` — exact string match. This is the pre-fix number and it must never move; every
  set's copy is asserted against a fresh recount on each run.
- `n_matched_via_alias` — resolved through `00_data/references/symbol_alias/symbol_alias_map.csv`
  with the ownership guard passed. This is what the fix recovered.
- `n_alias_rejected_ambiguous` — the ownership guard firing: the candidate is the official symbol
  of a *different* gene. `Gck→Gk`, `Tacr1→Spr`, `Slc6a4→Htt`, `Nkx3-1→Bax` are the kind of trap
  it catches, and accepting one would attach one gene's expression to another gene's set
  membership.
- `n_alias_flagged_for_review` — the guard would accept it; a human declined. Two pairs:
  `Ndufb1→Ndufb1-ps`, because handing a curated NDUFB1 set member to a pseudogene-named row is a
  judgement about the vM25 annotation rather than about nomenclature, and `H2ac25→H2aw`, because
  histone cluster renaming is many-to-many and the one-to-one condition passing there is close to
  accidental.
- `n_expression_filtered` — **structurally 0 here.** The delivered CPM table's 19,679 unique
  `gene_name` values are exactly the 19,679 modelled symbols; the duplicate-symbol collapse at
  `01_mapping_qc.R` drops Ensembl ids, never symbols. There is one vocabulary layer, so there is
  no expression filter to blame a loss on. The column is kept so the schema matches the human
  compartment's.
- `n_below_detection` — not in the delivered quantification. This is the terminal bucket.
- `n_absent_from_reference` — **always 0, with `reference_vocabulary_available = FALSE`.** It
  would need the GENCODE vM25 feature list, and the collaborators delivered a CPM table rather
  than a quantification against a tracked GTF. Separating "not in the annotation" from "in the
  annotation and not detected" is the one statement this project cannot make, so it is reported
  as unavailable rather than guessed.

There is deliberately no recovery fraction and no pass/fail floor. A floor on `n_matched` would
reward exactly the conflation these columns exist to undo.

## geneset_symbol_ledger.csv

The same seven counts at **gene-set** granularity — one row per set the reference ships,
including the sets the size filter then drops, which is precisely where a vocabulary loss would
otherwise become invisible. Adds `database`, `n_overlap_before_alias` (the pre-fix matched count,
kept beside the new one so both stay legible), `set_size_resolved` (unique matrix symbols the set
resolves to), `n_alias_collapsed` (pairs that landed on a symbol the set already carried), and
the applied / flagged / rejected pairs as `/`-joined `REF->MATRIX` strings.

Set granularity is not a convenience. Every claim in this family is per-set, and the two that
matter most cannot be seen in a roll-up: `MITOPATHWAYS_OXPHOS.Complex_V` matched **7 of its 22**
genes and `.CV_subunits` **2 of 17**, both fell under `gsea_min_size = 15`, and both were
therefore absent from `03_results/master/master_gsea_table.csv` entirely — no contrast, no row,
no NES. Resolved they are 22 and 17 and testable in all seven contrasts. The matrix carries the
legacy ATP-synthase block `Atp5a1 … Atpif1` while MitoCarta 3.0 ships `Atp5f1a … Atp5po`, and
that alone is 15 of the 16 genes `MITOPATHWAYS_OXPHOS` (129 → 145) was missing.

## geneset_alias_applied.csv

Every `(database, gene_set, reference_symbol, matrix_symbol)` pair actually applied. The audit
trail behind `n_matched_via_alias`: a recovery should never have to be inferred from a count.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/04_gsea_set_prep.R` | `resolve_and_filter` → `helpers/symbol_alias.R::{resolve_sets,symbol_ledger}` | `symbol_alias.map_path; symbol_alias.matrix_vocabulary; symbol_alias.flagged_for_review; thresholds.gsea_min_size=15; thresholds.gsea_max_size=500` | `00_data/references/symbol_alias/symbol_alias_map.csv` + `03_results/objects/gene_universe.txt` + the 8 MSigDB collections and 5 custom databases named under `databases:` |
