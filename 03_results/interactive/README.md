# GSEA Bump Interaction Dashboard

Standalone interactive HTML dashboard and companion CSV for exploring pathway-level Genotype × Temperature interaction. The dashboard exposes two views — "Interaction (single panel)" and "Heat trajectory (WT vs cGASKO)" — toggled via a radio selector at the top of the sidebar; all filter controls (database checkboxes, pattern checkboxes, |NES| range sliders, interaction 3-state filter, keyword highlight) are shared across both views and switching views preserves every active filter.

## interactive/gsea_bump_interaction.html

### View 1 — Interaction (single panel)

**Finding:** The interactive bump chart reveals pathway-level responses to heat in WT and cGAS-KO cells, showing that while most pathways have similar heat responses in both genotypes, a select subset (including interferon-stimulated genes) shows significant cGAS-dependence of the heat response.

**How to read:**
- **x-axis:** The two genotype heat contrasts, `WT heat` and `cGASKO heat`.
- **y-axis:** Plotly y-value displays GSEA NES (Normalized Enrichment Score) or ordinal Rank (1 = highest NES), toggled via the sidebar.
- **Lines:** Each pathway is represented by a line connecting its WT heat and cGASKO heat NES/Rank.
- **Gap:** The vertical gap between the endpoints represents the visual heat-response difference between the genotypes. The significance of this difference is driven by the Genotype × Temperature `Interaction` contrast.
- **Flat Lines:** A flat or non-significant line indicates **no detectable cGAS-dependence (n=5)**.
- **Curvature:** Curved lines bow based on the Genotype × Temperature `Interaction` contrast: curvature magnitude corresponds to `|Interaction NES|` and direction corresponds to the sign of `Interaction NES` (upward bow for positive, downward bow for negative). Curves are only drawn for pathways with significant interaction.
- **Color modes:** Color by pathway pattern, WT heat NES, cGASKO heat NES, or Interaction NES. Pattern color uses project-canonical palette per trajectory-classification category. NES color uses a blue–white–orange diverging scale (blue = negative/down-regulated, white = near zero, orange = positive/up-regulated; capped at |NES| = 3.5).
- **Highlight search:** Type any fragment of a pathway description; matching pathways are drawn on top in bold black with all other pathways dimmed.
- **Interaction filter:** Filter down to "Interaction-affected only (cGAS-dependent)" or "No detectable cGAS-dependence (n=5)" to isolate or exclude cGAS-dependent responses.
- **Claim Tier:** **L3** (enrichment statistics derived from ranked gene list; causal mechanism excluded).
- **Sample mapping** is owner-confirmed, 20 of 20 libraries concordant with the label-blind marker call. Every comparison here still awaits validation in downstream cohorts.
- **Caveat:** GSEA pathway-level interaction effects represent collective transcriptional changes and cannot be uniquely attributed to a single transcription factor (non-identifiability caveat, see fig3p/q/r).

### View 2 — Heat trajectory (WT vs cGASKO)

**Finding:** The heat-trajectory view renders each pathway as a slope from a shared 37 °C reference origin to its heat-contrast NES in two side-by-side Plotly panels (left = WT, right = cGASKO), enabling direct visual comparison of the magnitude and direction of heat responses between genotypes. Most pathways show matched slopes in both panels, consistent with **no detectable cGAS-dependence (n=5)**; pathways with a significant Genotype × Temperature Interaction NES that are steep in WT but flat in cGASKO represent candidates for cGAS-dependent heat responses.

**How to read:**
- **Layout:** Two side-by-side Plotly panels sharing a locked y-axis range. Left panel = WT; right panel = cGASKO. Panel titles ("WT" and "cGASKO") appear above each panel.
- **x-axis (both panels):** Temperature in °C; only two tick positions are shown — 37 °C (reference) and 39 °C (heat condition). The x-range is fixed to [36.5, 39.5].
- **y-axis (shared range):** GSEA NES labeled "NES (slope from 37 °C reference)". The y-range is symmetric around zero and set to encompass the maximum |NES| of visible pathways (scaled by 1.08), enforced identically on both panels so slopes are directly comparable.
- **Lines:** Each pathway is a straight line from (37 °C, 0) to (39 °C, heat-contrast NES). Left panel uses WT_heat NES; right panel uses KO_heat NES. A positive slope = pathway up-regulated at 39 °C vs 37 °C in that genotype; a negative slope = down-regulated.
- **CRITICAL CAVEAT — slope-from-reference construction:** The 37 °C anchor sits at 0 by construction because this view plots a single heat contrast (39 vs 37 °C) per genotype, not absolute per-temperature pathway activity. All lines share a common origin regardless of the pathway's true 37 °C baseline activity. This is a slope chart, NOT a trajectory of absolute pathway activity; the y-metric is NES only (the Rank toggle is hidden in this view because rank at the 37 °C reference is undefined). A true absolute 37/39 trajectory would require per-group ssGSEA scores, which were not computed for this dataset.
- **Interpreting genotype comparison:** Compare a pathway's slope in the WT panel to its slope in the cGASKO panel. A pathway steep in WT but flat in cGASKO (with a significant Interaction padj < 0.05) indicates cGAS-dependent heat response. Matched slopes in both panels = **"no detectable cGAS-dependence (n=5)"** — this phrasing is required; matched slopes are not proof of independence, only absence of detectable signal at this sample size.
- **Color modes:** Identical to View 1 — pattern color or NES diverging scale (blue–white–orange, |NES| capped at 3.5). Color identity is preserved when switching between views.
- **Highlight search:** Matching pathways are drawn with a thick black outline and coloured inner line on top of dimmed non-matching pathways in both panels simultaneously.
- **Y-metric toggle and Curved-lines controls:** Both are hidden in this view. Trajectory mode always uses NES; Bézier curvature encoding does not apply to slope-from-reference lines.
- **Shared filter pipeline:** All sidebar filters (database, pattern, |NES| range sliders, interaction filter, keyword highlight) apply identically to this view. Switching from View 1 to View 2 preserves every active filter.
- **Claim Tier:** **L3** (enrichment statistics; mechanistic attribution excluded). **Sample mapping** is owner-confirmed, 20 of 20 libraries concordant with the label-blind marker call. Non-identifiability caveat applies: a pathway-level interaction signal cannot be uniquely attributed to a single transcription factor (see fig3p/q/r).

**Provenance (both views):**

| Script | Function | Config | Input |
| :--- | :--- | :--- | :--- |
| `02_analysis/scripts/bump_dashboard.py` | `DashboardPipeline.run()` via `bump_dashboard.presentation.html_fragments` (`renderTrajectoryChart` JS in `SCRIPT`) | `02_analysis/config/analysis_config.yaml` | `03_results/master/master_gsea_table.csv` |

---

## interactive/gsea_bump_interaction.csv

**Finding:** Flat wide table containing the source statistics (NES, padj, ranks, and computed response pattern category) for all 4,798 analyzed pathways across WT heat, cGASKO heat, and G×T Interaction contrasts.

**Provenance Table:**

| Script | Function | Config | Input |
| :--- | :--- | :--- | :--- |
| [bump_dashboard.py](file:///workspaces/STING-cGAS-GSE329522/02_analysis/scripts/bump_dashboard.py) | Renders the interactive pathway-trajectory bump-chart dashboard and companion CSV | [analysis_config.yaml](file:///workspaces/STING-cGAS-GSE329522/02_analysis/config/analysis_config.yaml) | [master_gsea_table.csv](file:///workspaces/STING-cGAS-GSE329522/03_results/master/master_gsea_table.csv) |
## master_unified.csv (explorer bundle view)
**The consolidated pathway-explorer input accumulator (entity_type-tagged Pathway/TF/PROGENy/GATOM rows, lowercase nes, genes_full_set = membership ∩ atlas, contrast-invariant): the cross-stratum enrichment/activity universe the dashboards lay out.**
| | |
|---|---|
| Script   | `02_analysis/scripts/pathway_explorer_adapter/consolidate_explorer_bundle.R` |
| Function | (validation + manifest) |
| Config   | `paths.master = 03_results/master/` |
| Input    | `03_results/master/master_unified.csv; 03_results/master/atlas_gene_universe.txt` |

## master_de_table.csv (explorer schema view)
**The explorer-schema DE view (gene_symbol, t, logFC, adj.P.Val, contrast): the t-ranked per-contrast gene list that drives the dashboards' running-sum panel.**
| | |
|---|---|
| Script   | `02_analysis/scripts/pathway_explorer_adapter/consolidate_explorer_bundle.R` |
| Function | (validation) |
| Config   | `paths.master = 03_results/master/` |
| Input    | `03_results/master/master_de_table.csv` |
## pathway_explorer_<contrast>.html
**Per-contrast interactive pathway-explorer dashboard (UMAP of GSEA pathways + CollecTRI TF + PROGENy activities, NES-colored, FDR-sliderable, with the t-ranked running-sum panel): this dashboard enables interactive cross-entity similarity mapping and detailed exploration of regulatory cascades.**
| | |
|---|---|
| Script   | `02_analysis/scripts/pathway_explorer_adapter/run_pathway_explorer.sh` |
| Function | `generate_all_dashboards()` (pathway_explorer) |
| Config   | `paths.interactive = 03_results/interactive/; paths.master = 03_results/master/` |
| Input    | `03_results/master/master_unified.csv; 03_results/master/master_de_table.csv` |

## index.html
**Landing page linking every per-contrast pathway-explorer dashboard: the entry point for interactive exploration of MSigDB/custom pathways, TF activities, and PROGENy activities.**
| | |
|---|---|
| Script   | `02_analysis/scripts/pathway_explorer_adapter/run_pathway_explorer.sh` |
| Function | `generate_index_page()` (pathway_explorer) |
| Config   | `paths.interactive = 03_results/interactive/; paths.master = 03_results/master/` |
| Input    | `03_results/master/master_unified.csv` |

