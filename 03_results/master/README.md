# master — The accumulator tables

Cross-stage tables that gather one row per result from the individual stages. A stage computes;
these files collect. Read them for every effect in one place rather than opening each stage.

Twelve files in three groups. **Differential expression** (`master_de_genes.csv`,
`master_de_table.csv`, `universe_frame.csv`, `atlas_gene_universe.txt`) is per gene.
**Enrichment and activity** (`master_gsea_table.csv`, `master_tf_activities.csv`,
`master_progeny_activities.csv`, `master_gatom_modules.csv`) is per gene set, and the four share
one schema so they stack. **Explorer substrate** (`master_unified.csv`, `explorer_universe.csv`,
`explorer_manifest.json`) is that stack plus gene membership, built for
[`../interactive/`](../interactive/). `sample_metadata.csv` sits on its own and keys the design.

Seven contrasts run through every table: `WT_heat` and `KO_heat` (39 versus 37 °C within a
genotype), `Temp_main` and `Geno_main` (the two marginal effects), `Geno_at_37` and `Geno_at_39`
(wild-type versus cGAS-knockout within a temperature), and `Interaction` (the heat response's
cGAS-dependence). Positive statistics point to the numerator — 39 °C, or wild-type.

**Symbol vocabulary.** Every gene symbol here is frozen to the build the counts were quantified
against, GRCm38 + GENCODE vM25. Reference gene sets ship current symbols, so
`00_symbol_alias_map.R` resolves reference symbols into this vocabulary before matching and
`04_gsea_set_prep.R` applies the map. Match a symbol from an outside source through that map;
[`../objects/README.md`](../objects/) carries the ledger of what it recovered.

---

## Differential expression

### `master_de_genes.csv` · `master_de_table.csv`

**137,753 rows — 19,679 genes × 7 contrasts. The limma-trend gene-level result.**

One row per (gene × contrast). `t` is the moderated t-statistic and is what every ranked list in
this compartment is built from; `logFC` is the effect on the log2 scale. `significant` and
`direction` are cut at the config FDR threshold, so a reader filters without re-deciding the
cutoff.

| Column | Meaning |
|---|---|
| `gene_symbol` | MGI symbol in the GENCODE vM25 vocabulary. |
| `ensembl` | Ensembl gene id, the stable key across symbol vintages. |
| `logFC`, `t` | Effect size, and the moderated t the rankings use. |
| `P.Value`, `adj.P.Val` | Raw and BH-adjusted p. |
| `contrast` | One of the seven. |
| `significant`, `direction` | The FDR call and its sign. |

`master_de_table.csv` holds the same 137,753 rows under the name the pathway explorer expects.
It is written through a rounding pass, so the two carry the same numbers to nine significant
digits and differ in bytes. Read `master_de_genes.csv`.

### `universe_frame.csv`

**19,679 rows — one per modelled gene, all seven contrasts side by side.**

The wide form: `<contrast>_logFC`, `<contrast>_t`, `<contrast>_padj`, `<contrast>_sig` and
`<contrast>_direction` for each of the seven, on one row per gene. Use it to ask whether a gene
moves in two contrasts at once without self-joining the long table. `detected`, `in_universe` and
`expressed` mark the gene's standing in the ranked universe; all 19,679 rows carry all three,
which makes this frame the compartment's background list.

### `atlas_gene_universe.txt`

**19,679 lines — the background gene list, one symbol per line, sorted.**

The symbol universe every enrichment test in this compartment ran against. A gene absent from
this file was never in the denominator, so its absence from a gene set's matched count is a
vocabulary fact. `explorer_manifest.json` records this file's SHA-256, so a bundle can be checked
against the universe it was built on.

---

## Enrichment and activity

### `master_gsea_table.csv`

**33,740 rows — gene-set enrichment across 15 databases × 7 contrasts.**

One row per (gene set × contrast). `nes` is the normalised enrichment score, positive when the
set sits toward the numerator end of the ranking. `set_size` counts the set's genes that reached
the ranked universe, which runs below the set's published size. `core_enrichment` is the leading
edge as a `/`-separated symbol list, and its length is the leading-edge count other tables quote.

The fifteen databases are Hallmark, KEGG, Reactome, WikiPathways, GO_BP, GO_CC, GO_MF,
MitoPathways, MitoXplorer, TransportDB, TF_Targets, HSR_lens, TCR_activation, PROGENy and
CoReSh_derived. **Filter on `database` before reading a p-value**: adjustment happens within a
database, so a `padj` from one is off the scale of a `padj` from another.
`06b_gsea_pooled_padj.R` supplies the pooled alternative for a cross-database comparison.

| Column | Meaning |
|---|---|
| `pathway_id` | The set's own identifier, e.g. `HALLMARK_HYPOXIA`. |
| `pathway_name` | Display name. |
| `database` | The collection, and the family the `padj` was adjusted within. |
| `nes`, `pvalue`, `padj` | Enrichment statistic and its p-values. |
| `set_size` | The set's genes present in the ranked universe. |
| `core_enrichment` | Leading-edge symbols, `/`-separated. |
| `contrast`, `direction` | One of the seven, and the sign of `nes`. |

### `master_tf_activities.csv`

**4,869 rows — inferred transcription-factor activity from two regulatory networks.**

Same columns as the enrichment table, so the two stack. `nes` is the univariate linear-model
activity score, positive when the factor's targets sit toward the numerator end. `set_size` is
the regulon's genes present in the ranked universe.

`database` separates the two networks and their coverage differs. `CollecTRI` is the primary run
— 4,606 rows, 658 factors across all seven contrasts. `DoRothEA_ABC` is the comparator network
held for the method comparison in `03b_decoupler_method_comparison.R` and scored on `WT_heat`
alone — 263 rows, 263 factors. A factor therefore appears twice in `WT_heat` with two scores;
filter on one network before ranking.

### `master_progeny_activities.csv`

**98 rows — 14 PROGENy pathways × 7 contrasts.**

Signalling-pathway activity from the PROGENy consensus signatures, in the shared schema. The
table is complete by construction, so a missing row is an error. `09_activity_progeny.R` also
appends these rows to `master_gsea_table.csv` under `database = PROGENy`; count PROGENy once
when pooling.

### `master_gatom_modules.csv`

**3 rows — the metabolic modules GATOM returned, under `database = GATOM_KEGG`.**

Maximally-regulated metabolic subnetworks found on the KEGG atom-transition graph, carried in the
shared schema so they join the rest. Three rows is the whole result; the module topology lives in
[`../09_gatom/`](../09_gatom/).

---

## Explorer substrate

### `master_unified.csv`

**38,307 rows — the four enrichment and activity tables stacked, plus gene membership.**

The shared schema with two columns added: `entity_type` names which table a row came from
(33,600 `Pathway`, 4,606 `TF`, 98 `PROGENy`, 3 `GATOM`) and `genes_full_set` carries the set's
membership intersected with `atlas_gene_universe.txt`. Row counts run below the source tables
because a row with no members inside the universe cannot be embedded. This is the pathway
explorer's input; read the source tables for analysis.

### `explorer_universe.csv`

**6,427 rows — one per distinct gene set, contrast collapsed away.**

`master_unified.csv` de-duplicated to (`pathway_id`, `entity_type`, `genes_full_set`), which is
what the explorer's Jaccard neighbour graph and gene-set UMAP are fitted on: 5,304 pathways,
1,108 factors, 14 PROGENy pathways, 1 GATOM module. The table carries membership, so it has no
`nes` or `padj` column.

### `explorer_manifest.json`

The provenance record for the explorer bundle. `produced_by`, `produced_at` and `git_commit` say
which run built it; `atlas_file`, `atlas_sha256` and `atlas_size` pin the gene universe it was
built against; `namespace` is `mouse_symbol_GRCm39` and `genes_full_set_basis` is
`atlas_intersected`, so a reader knows the membership lists were cut to the universe. `n_rows`,
`contrasts` and `entity_types` restate the bundle's shape. Check `atlas_sha256` against
`atlas_gene_universe.txt` before trusting a stale bundle.

---

## The design

### `sample_metadata.csv`

**20 rows — one per library. The design, and the evidence for it.**

`library_id`, `genotype`, `temp` and `group` are the 2×2 design: wild-type and cGAS-knockout at
37 °C and 39 °C, five libraries per cell. `gsm_id` is the GEO accession under the marker-derived
mapping; `gsm_id_positional` is the accession a naive positional join would assign. The two
differ for 10 of the 20 libraries, because the deposited CPM columns run temperature-major while
the accessions run genotype-major. **Join on `gsm_id`.**

`mapping_status` reads `CONFIRMED` for all 20: each library's data-derived label — the heat-shock
thermometer for temperature, `Cgas` expression for genotype — agrees with the owner's sample
sheet. [`../01_qc/`](../01_qc/) carries that check.
