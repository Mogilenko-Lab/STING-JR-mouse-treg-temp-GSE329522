# master — The accumulator tables

Cross-stage tables that gather one row per result from the individual stages. A stage computes,
and these files collect. Read them for every effect in one place.

Twelve files in three groups. **Differential expression** (`master_de_genes.csv`,
`master_de_table.csv`, `universe_frame.csv`, `atlas_gene_universe.txt`) is per gene.
**Enrichment and activity** (`master_gsea_table.csv`, `master_tf_activities.csv`,
`master_progeny_activities.csv`, `master_gatom_modules.csv`) is per gene set, and the four share
one schema so they stack. `sample_metadata.csv` sits on its own and keys the design.

Seven contrasts run through every table: `WT_heat` and `KO_heat` (39 versus 37 °C within a
genotype), `Temp_main` and `Geno_main` (the two marginal effects), `Geno_at_37` and `Geno_at_39`
(wild-type versus cGAS-knockout within a temperature), and `Interaction` (the heat response's
cGAS-dependence). Positive statistics point to the numerator — 39 °C, or wild-type.

**Symbol vocabulary.** Every gene symbol here is frozen to the build the counts were quantified
against, GRCm38 + GENCODE vM25. Reference gene sets ship current symbols, so
`00_symbol_alias_map.R` resolves reference symbols into this vocabulary before matching and
`04_gsea_set_prep.R` applies the map. Match a symbol from an outside source through that map.
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

`master_de_table.csv` holds the same 137,753 rows under an alternate name kept for a
downstream consumer. It is written through a rounding pass, so the two carry the same numbers to
nine significant digits and differ in bytes. Read `master_de_genes.csv`.

### `universe_frame.csv`

**19,679 rows — one per modelled gene, all seven contrasts side by side.**

The wide form: `<contrast>_logFC`, `<contrast>_t`, `<contrast>_padj`, `<contrast>_sig` and
`<contrast>_direction` for each of the seven, on one row per gene. Use it to ask whether a gene
moves in two contrasts at once without self-joining the long table. `detected`, `in_universe` and
`expressed` mark the gene's standing in the ranked universe. All 19,679 rows carry all three,
which makes this frame the compartment's background list.

### `atlas_gene_universe.txt`

**19,679 lines — the background gene list, one symbol per line, sorted.**

The symbol universe every enrichment test in this compartment ran against. A gene absent from
this file was never in the denominator, so its absence from a gene set's matched count is a
vocabulary fact.

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
alone — 263 rows, 263 factors. A factor therefore appears twice in `WT_heat` with two scores, so
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

---

## Signature provenance

Every row in these tables keys to one of the sources below.

| Source | Provenance |
|---|---|
| The libraries and every DE row | **GSE329522** — bulk RNA-seq of induced regulatory T cells differentiated from primary murine splenic CD4⁺ T cells, genotype (wild-type, cGAS-knockout) × temperature (37 °C, 39 °C), five biological replicates per group, twenty libraries. The owner's sample sheet dated 2026-07-22 is the authoritative library-to-condition mapping. |
| `database` = Hallmark · KEGG · Reactome · WikiPathways · GO_BP · GO_MF · GO_CC · TF_Targets | MSigDB collections H, C2:CP:KEGG, C2:CP:REACTOME, C2:CP:WIKIPATHWAYS, C5:GO:BP, C5:GO:MF, C5:GO:CC and C3:TFT:GTRD, retrieved through msigdbr 26.1.0 for *Mus musculus* and frozen by `04_gsea_set_prep.R`. |
| `database` = TransportDB · MitoPathways · MitoXplorer | Prebuilt mouse gene-set objects from the RNAseq-toolkit reference tree: `transportdb/processed/Mus_musculus/transportdb_genesets.rds`, the MitoCarta 3.0 build `mitocarta3.0/processed/Mus_musculus/mito_mitopathways.rds`, and `mitoxplorer3.0/processed/Mus_musculus/mito_mitoxplorer.rds`. |
| `database` = HSR_lens | `HSR_core` and `HSR_sensitivity`, the union of `REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS`, `REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE` and `GOBP_RESPONSE_TO_HEAT` from MSigDB v2026.1.Hs through msigdbr 26.1.0, built by `00d_curate_temp_hsr.R` and `00e_curate_temp_hsr_lens.R`. |
| `database` = TCR_activation | A frozen 66-symbol human T-cell activation panel curated in this repository, mouse symbols through `msigdbr(species = "Mus musculus")`, built by `00f_curate_tcr_activation.R`. |
| `database` = CoReSh_derived | Modules mined from the public mouse GEO compendium on Synapse, **syn66227307** (mmu), read from the shared read-only reference cache. [`../08_coresh/`](../08_coresh/) names the seeding query and the source accession of each module. |
| `database` = PROGENy | `progeny::getModel("Mouse", top = 500)`, fourteen footprints, cached by `00c_prepare_networks.R`. |
| `database` = CollecTRI · DoRothEA_ABC | CollecTRI from the human regulon table pinned at `00_data/references/networks/CollecTRI_regulons_human.csv`, mapped to mouse. DoRothEA from the bundled mouse regulons at confidence A, B and C. Both built by `00c_prepare_networks.R`, because `decoupleR::get_collectri()` fails against OmnipathR 3.18.4. |
| `database` = GATOM_KEGG | Modules found on the KEGG atom-transition network under `00_data/references/gatom/`. |
