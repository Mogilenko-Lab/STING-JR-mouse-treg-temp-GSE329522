# objects — Recomputable state, and the three tables that travel with it

This directory holds checkpoints, so the `.rds` and `.h5ad` files here stay untracked and
regenerate from the scripts that wrote them. Three CSVs are tracked, because they are
accounting: a reader has to be able to check them without re-running anything.

## What the checkpoints hold

| Group | Files | Written by |
|---|---|---|
| Expression and design | `01_eda.rds`, `sample_metadata.rds`, `gene_universe.txt`, `gene_universe_vec.rds`, `universe_frame_obj.rds` | `01_mapping_qc.R`, `00_setup_metadata.R`, `11_emit_universe.R` |
| Differential expression | `02_de_results.rds` | `02_de_limma_trend.R` |
| Networks and activity | `net_collectri_mouse.rds`, `net_dorothea_mouse_ABC.rds`, `net_progeny_mouse.rds`, `03_tf_collectri.rds`, `03_tf_forensics.rds`, `09_progeny_activity.rds`, `networks_prep_log.txt` | `00c_prepare_networks.R`, `03_decoupler_tf.R`, `09_activity_progeny.R` |
| Gene-set collections | `geneset_msigdb_*.rds` (8), `geneset_custom_*.rds` (5) | `04_gsea_set_prep.R` |
| Enrichment results | `gsea_msigdb_*.rds`, `gsea_custom_*.rds`, `gsea_coresh_*.rds` (7 contrasts each) | `05_gsea_msigdb_run.R`, `06_gsea_custom_run.R`, `08_coresh_derived_gsea.R` |
| CoReSh | `coresh_ranked.rds`, `coresh_derived_sets.rds`, `coresh_gmt_lines.rds`, `coresh_query_entrez.rds` | `07_coresh_search.R`, `08_coresh_derived_gsea.R` |
| Metabolic modules | `10_gatom_all.rds`, `10_gatom_<contrast>.rds` (4) | `10_gatom_modules.R` |
| Signatures and decompositions | `16_synthesis.rds`, `17_signature_sets.rds`, `19_hsr_decomp_gsea.rds`, `20_semantic_decomp.rds` (plus one per GOSemSim build), `24_go_decomposition.rds`, `26_arm_composition.rds` | The stage that owns each |

`20_semantic_decomp__GOSemSim-2.36.0.rds` and `…-2.39.2.rds` sit side by side on purpose: the
filename carries the build, so the two can be compared without either evicting the other. See
[`../13_semantic_decomp/`](../13_semantic_decomp/).

---

## The vocabulary ledger

Three tracked CSVs record what happened when reference gene sets met this matrix's symbol
vocabulary. They exist because **a gene symbol that fails to match is a vocabulary result until
proven a biological one**, and a bare match count conflates the two.

### `geneset_manifest.csv`

The authoritative index of every gene-set collection available to `05_gsea_msigdb_run.R`,
`06_gsea_custom_run.R` and the viz arm. One row per collection: `database`, `type`, `n_sets`,
`source`, `cache_path`, then the vocabulary ledger rolled up over the collection's unique
reference symbols. Written by `04_gsea_set_prep.R`.

**The ledger columns close, exactly:**

```
n_unique_set_genes == n_matched + n_matched_via_alias + n_alias_flagged_for_review +
                      n_alias_rejected_ambiguous + n_expression_filtered +
                      n_below_detection + n_absent_from_reference
```

- `n_matched` — exact string match. This is the pre-fix number and it holds still, with every
  set's copy asserted against a fresh recount on each run.
- `n_matched_via_alias` — resolved through
  `00_data/references/symbol_alias/symbol_alias_map.csv` with the ownership guard passed. This
  is what the fix recovered.
- `n_alias_rejected_ambiguous` — the ownership guard firing: the candidate is the official symbol
  of a *different* gene. `Gck→Gk`, `Tacr1→Spr`, `Slc6a4→Htt` and `Nkx3-1→Bax` are the kind of
  trap it catches, and accepting one would attach one gene's expression to another gene's set
  membership.
- `n_alias_flagged_for_review` — the guard would accept it and a human declined. Two pairs:
  `Ndufb1→Ndufb1-ps`, because handing a curated NDUFB1 set member to a pseudogene-named row
  passes judgment on the vM25 annotation itself, and `H2ac25→H2aw`, because histone cluster
  renaming is many-to-many and the one-to-one condition passing there is close to accidental.
- `n_expression_filtered` — **structurally 0 here.** The delivered CPM table's 19,679 unique
  `gene_name` values are exactly the 19,679 modelled symbols, and the duplicate-symbol collapse
  at `01_mapping_qc.R` drops Ensembl ids while keeping every symbol. One vocabulary layer means
  no expression filter to attribute a loss to. The column is kept so the schema matches the
  human compartment's.
- `n_below_detection` — outside the delivered quantification. This is the terminal bucket.
- `n_absent_from_reference` — **always 0, with `reference_vocabulary_available = FALSE`.** It
  would take the GENCODE vM25 feature list, and what the collaborators delivered is a CPM table.
  Separating "outside the annotation" from "in the annotation and undetected" is the one
  statement this compartment leaves unmade, so it is reported as unavailable.

There is deliberately no recovery fraction and no pass/fail floor. A floor on `n_matched` would
reward exactly the conflation these columns exist to undo.

### `geneset_symbol_ledger.csv`

The same seven counts at **gene-set** granularity — one row per set the reference ships,
including the sets the size filter then drops, which is precisely where a vocabulary loss would
otherwise become invisible. It adds `database`, `n_overlap_before_alias` (the pre-fix matched
count, kept beside the new one so both stay legible), `set_size_resolved` (unique matrix symbols
the set resolves to), `n_alias_collapsed` (pairs landing on a symbol the set already carried),
and the applied, flagged and rejected pairs as `/`-joined `REF->MATRIX` strings.

Set granularity is load-bearing. Every claim in this family is per-set, and the two that matter
most are invisible in a roll-up: `MITOPATHWAYS_OXPHOS.Complex_V` matched **7 of its 22** genes
and `.CV_subunits` **2 of 17**, both fell under `gsea_min_size = 15`, and both were therefore
absent from `../master/master_gsea_table.csv` entirely — no contrast, no row, no NES. Resolved
they are 22 and 17 and testable in all seven contrasts. The matrix carries the legacy
ATP-synthase block `Atp5a1 … Atpif1` while MitoCarta 3.0 ships `Atp5f1a … Atp5po`, and that alone
is 15 of the 16 genes `MITOPATHWAYS_OXPHOS` (129 → 145) was missing.

### `geneset_alias_applied.csv`

Every `(database, gene_set, reference_symbol, matrix_symbol)` pair actually applied. The audit
trail behind `n_matched_via_alias`, so a recovery never has to be inferred from a count.

All three are written by `04_gsea_set_prep.R` through
`helpers/symbol_alias.R::{resolve_sets, symbol_ledger}`, against
`00_data/references/symbol_alias/symbol_alias_map.csv` and `gene_universe.txt`, under
`thresholds.gsea_min_size = 15` and `gsea_max_size = 500`.

---

## Signature provenance

The checkpoints cache resources built elsewhere. Each one's source:

| Checkpoint | Provenance |
|---|---|
| `geneset_msigdb_*.rds` (8) | MSigDB collections H, C2:CP:KEGG, C2:CP:REACTOME, C2:CP:WIKIPATHWAYS, C5:GO:BP, C5:GO:MF, C5:GO:CC and C3:TFT:GTRD, retrieved through msigdbr 26.1.0 for *Mus musculus*. |
| `geneset_custom_*.rds` (5) | TransportDB, the MitoCarta 3.0 MitoPathways build and MitoXplorer, all prebuilt mouse objects from the RNAseq-toolkit reference tree under `01_modules/RNAseq-toolkit/data/references/`, plus this compartment's own `HSR_lens` and `TCR_activation`, frozen under `00_data/references/gene_sets/`. `geneset_manifest.csv` carries the exact source path of each. |
| `net_collectri_mouse.rds` | The human TF→target regulon table pinned at `00_data/references/networks/CollecTRI_regulons_human.csv`, mapped to mouse symbols. Built locally because `decoupleR::get_collectri()` fails against OmnipathR 3.18.4. |
| `net_dorothea_mouse_ABC.rds` | The bundled mouse DoRothEA regulons at confidence classes A, B and C. |
| `net_progeny_mouse.rds` | `progeny::getModel("Mouse", top = 500)`, fourteen pathway footprints. |
| `coresh_*.rds` | Sweeps and derived sets over the public mouse GEO compendium on Synapse, **syn66227307** (mmu), read from the shared read-only reference cache. |
| `10_gatom_*.rds` | Modules found on the KEGG atom-transition network under `00_data/references/gatom/`. |
| Everything else | Derived in this compartment from **GSE329522** — induced regulatory T cells from primary murine splenic CD4⁺ T cells, genotype × temperature, five libraries per cell of the design. |
