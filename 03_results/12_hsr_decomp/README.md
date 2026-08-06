# 12_hsr_decomp — What the frozen 39 °C set is made of

`WT_heat_up` is 213 mouse genes selected by the export gate. This stage asks what those genes
are, by direct membership against two curated lenses: a heat-shock core and a TCR /
immediate-early activation panel. **198 of the 213 fall in neither.**

That answer raises a second question, and the stage answers it in the same directory. Is the
heat-shock response absent from the 39 °C response? No — warming induces the curated `HSR_core`
in both genotypes (`WT_heat` NES +2.049 at adjusted p 1.9e-05; `KO_heat` +2.065 at 3.4e-05).

The third question resolves the two. **The gate thresholds effect size, not pathway.**
`WT_heat_up` genes sit at median rank percentile 1.4% with median t 10.80; `HSR_core` genes sit
at 24.4% with median t 2.05; the deepest gene the gate admitted sits at 10.21%. The heat-shock
response moved moderately and the gate kept only extreme movers.

The last panel is the handoff. The mouse gate output is 213 genes, the human set the disease
compartments receive is 202, and 17 of that 202 are present in every one of the eighteen human
ranked lists that consume it.

**A standing distinction this whole stage rests on.** A lens enriching in a ranking and a lens
being contained in a thresholded set are two different measurements. Membership here is counted
over the whole arm; classifying leading-edge genes and re-testing those subsets would test genes
selected because they enriched.

**The ceiling on the reading.** `HSR_core` is a curated proteotoxic-stress reference and covers
more than fever. The 37/39 °C contrast is an experimental perturbation in mouse iTregs, so it
supports the response interpretation and leaves human fever untouched.

**Reading order:** `wtheatup_attribution` → `lens_nes_by_contrast` → `hsr_rank_position_panel` →
`gate_projection_bridge`, with the three `hsr_lens_membership_*` views corroborating the first.

**A third lens on the euler panels.** `Lombardi2022_HIF` is the Lombardi 2022 conserved HIF
signature, re-derived from that paper's supplement and mapped to 100 mouse symbols. It enters
here as a membership lens over gene symbols and is scored as a gene-set database by no stage in
these results; a pan-cancer consensus derived under hypoxia in cancer cells is the wrong
reference for a 39 °C in-vitro iTreg contrast. It meets `WT_heat_up` at 7 of 213 genes, 3.3%.
[`../06_gsea/README.md`](../06_gsea/) records where it was removed from and why.

**Re-render order for the bridge census.** Refresh the human compartments that publish the
ranked lists first, review `tables/source_hash_manifest.csv`, then re-run
`19_hsr_decomposition.R`. The stage writes `tables/source_reads_observed.csv` and stops if a
present sibling artifact differs from its pinned hash.

---

## Figures

### `figures/_overview/wtheatup_attribution.png`

**Of 213 genes in the frozen 39 °C set, 198 sit in neither curated lens.**
One stacked bar over the 213 genes, partitioned by direct set membership into `HSR_core_only` (3
genes, 1.4%), `TCR_activation_only` (12, 5.6%), `shared_both` (0) and `neither` (198). Segment
labels give count and percent, with the thin segments called out to the right of the bar.

This partitions membership. Rank enrichment is the separate measurement, and
`lens_nes_by_contrast` is where it is made. The 93% remainder is reported as a remainder: it is
unnamed and it supports no mechanism.
*Source* `tables/_overview/wtheatup_attribution.csv` ·
`02_analysis/scripts/19_hsr_decomposition_viz.R`.

### `figures/_overview/lens_nes_by_contrast.png`

**Warming induces the curated heat-shock core in both genotypes.**
Horizontal grouped bars. x, enrichment score on each contrast's signed-t ranking, positive
toward genes higher at 39 °C; y, lens, grouped by contrast. The NES prints at the bar end and
the right-hand column gives the effective set size and FDR as *n present of n nominal*.

`HSR_core` is the stronger of the two lenses in all three contrasts: `WT_heat` +2.05 (FDR
1.9e-05), heat main effect +2.01 (4.5e-05), `KO_heat` +2.07 (3.4e-05). Both lenses are wholly
recovered here — 47 of 47 and 66 of 66 — so the three contrasts are comparable on that count.
The cGAS-knockout row is what places the induction outside detectable cGAS-dependence at this
design's power.
*Source* `tables/_overview/lens_nes_by_contrast.csv` ·
`02_analysis/scripts/19_hsr_decomposition_viz.R`.

### `figures/_overview/hsr_rank_position_panel.png`

**The heat-shock core is a moderate responder and the gate kept extreme ones.**
Each point is one gene; ridgelines give each set's density along the same `WT_heat`
rank-percentile axis, with 0% the most-up end. The shaded band is the empirical gate span,
0.0051% to 10.21% — rank 1 to 2,010 of 19,679 — and the dashed line is its deep edge.

`HSR_core` genes centre at 24.4% with median t 2.05, roughly 2.4-fold deeper than the gate's
deepest admission. The gate's own output centres at 1.4% with median t 10.80. Only 213 of the
2,010 genes inside the band passed, because the gate also required log2FC ≥ 1; the deepest gene
it admitted is `Cpne6` at t 4.55, log2FC 1.03. The `WT_heat_up` row is left-shifted by
construction and is labelled as the gate's output.
*Source* `tables/_overview/hsr_rank_position_panel.csv` ·
`02_analysis/scripts/19_hsr_decomposition_viz.R`.

### `figures/_overview/gate_projection_bridge.png`

**The 213-gene mouse set and the 202-gene human set are one object either side of projection.**
Two halves. The top is a funnel of counts on a log axis, because the first and last steps differ
by three orders of magnitude: 19,679 measured genes, 2,010 inside the gate's rank span, 213
through the gate, 202 after ortholog projection, 17 present in every downstream ranked list.
Colour marks which species' symbols each row is counted in.

Curated-lens membership prints beside the two set rows: `HSR_core` 3 of 213 in mouse against a
47-gene lens and 2 of 202 in human against a 56-gene lens; `TCR_activation` 12 of 213 and 12 of
202 against a 66-gene lens both sides.

The two lenses behave differently under projection because they were built differently.
`TCR_activation` asserts a strictly one-to-one map — 66 human symbols to 66 unique mouse symbols,
the build stopping on any duplicate — so its twelve are the same twelve genes (`ATF3 CD40LG CD83
EGR1 EGR2 EGR3 IER3 IL2 LTA NR4A1 NR4A2 NR4A3`). `HSR_core` is assembled paralog-complete from
Reactome and GO under a one-to-many map, so it is 47 mouse against 56 human, and its
`Hspa1a`/`Hspa1b` pair collapses to one human symbol.

The bottom half spans the fewest to the most of the 202 present in any single ranked list of
each compartment, with the list count and length range printed; the dashed rule at 17 marks the
genes present in all of them. Recovery is a property of the analysis that produced each ranked
list, and the shallowest compartment sets the 17.

The human count read 199 until the ortholog step stopped reporting a vocabulary loss as a
biological one. Three genes of this arm arrive now (`DYNLT2B`, `FHIP1A`, `GARIN3`, from
`Tctex1d2`, `Fam160a1`, `Fam71b`), and neither lens count moved with them.
*Source* `tables/_overview/gate_projection_bridge.csv` ·
`02_analysis/scripts/19_hsr_decomposition_viz.R`.

### The three membership views

All three read one table, `tables/hsr_lens_membership.csv`, so they cannot disagree about a
count — only about how a count is drawn.

**`hsr_lens_membership_euler.png`** — area-proportional. `WT_heat_up` (213) against `HSR_core`
(47), `TCR_activation` (66) and `Lombardi2022_HIF` (100). Overlaps with the arm are 3, 12 and 7;
191 of 213 genes (89.7%) sit in no named lens. Areas are fitted to exact counts by `eulerr`
(stress 2.4e-17, diagError 1.4e-09).

**`hsr_lens_provenance_euler.png`** — the same treatment asking where `Interaction` originates.
`WT_heat_up` (213) against `Interaction_up` (7 human / 9 mouse), `Lombardi2022_HIF` (100) and
`HSR_core` (47). `WT_heat_up` overlaps `Interaction_up` at 0 genes. The unassigned remainder is
203 of 213 (95.3%). Fitted at stress 7.5e-19, diagError 4.3e-10.

**`hsr_lens_membership_venn.png`** — the same seven counts where geometry carries nothing. All
three circles are equal and equally overlapping by convention, so only the printed counts vary.
Every region is drawn at conventional size whatever its count, so a printed 0 marks a region
that exists and is empty: 198 in `WT_heat_up` alone, 3 and 12 shared with the two curated
lenses, 0 in the `HSR_core` ∩ `TCR_activation` region and 0 three-way.

**`hsr_lens_membership_upset.png`** — the same counts ranked as intersections. Each bar is one
intersection and the dot matrix below marks which sets it belongs to; the horizontal bars on the
left give the whole-set sizes (213, 47, 66). Every combination gets a column whether or not any
gene falls in it, so the two empty regions appear as zero-height bars labelled "empty". This view
stays readable past three sets, where circle layouts stop being drawable.
*Source* `tables/hsr_lens_membership.csv` and the four same-stem `_overview` tables ·
`02_analysis/scripts/19_hsr_decomposition_viz.R`.

---

## Tables

| File | What it holds | How to read it |
|---|---|---|
| `hsr_decomp_lens_nes.csv` | `HSR_core`, `HSR_sensitivity` and `TCR_activation` enrichment across `WT_heat`, `Temp_main` and `KO_heat`. | `nes > 0` means enrichment toward genes higher at 39 °C. `set_size` is the effective size, the members present in that ranking. `leading_edge_genes` is reported and never used to define a subset for re-testing. |
| `hsr_decomp_wtheatup_attribution.csv` | One row per `WT_heat_up` member, labelled by direct lens membership. | `HSR_core_only` means in the heat-shock core and outside the activation panel; `shared_both` means both; `wt_heat_t > 0` means higher at 39 °C. |
| `hsr_decomp_summary.csv` | The four attribution buckets over the 213-gene denominator, for `HSR_core` and the wider `HSR_sensitivity`. | `fraction` is `n / denominator`: the `HSR_core_only` row gives 3/213 = 1.4%, and the wider lens moves that to 5/213, which is the sensitivity read on the same count. |
| `hsr_decomp_overlap.csv` | Pairwise set overlaps. | `HSR_core` and `TCR_activation` share 0 genes; `WT_heat_up` shares 3 and 12 with them. `jaccard` is `n_intersect / union`. |
| `hsr_decomp_rank_concordance.csv` | Where each set sits in the full `WT_heat` signed-t ranking. | `median_rank_pct = 0` is the most-up end. A Wilcoxon diagnostic is carried; the rank-position panel displays none of it. |
| `hsr_decomp_conditional.csv` | Each lens's NES after the other lens is purged from the ranking. | `delta_nes = nes_cond − nes_uncond`, about +0.005, which supports treating the two as separate lenses for this audit. |
| `hsr_lens_membership.csv` | The shared source of all three membership panels: seven mouse region counts, the external human `WT_heat_up` / `HSR_core` overlap, and the conditional NES deltas. | The seven rows partition the union of the three sets. `human_wt_hsr_intersect` is carried when the external table is present and well formed. |
| `source_hash_manifest.csv` · `source_reads_observed.csv` | The pinned sibling artifacts the bridge census may read, and what it actually read. | `status = read` means the file existed, matched its pinned SHA-256, and was consumed. The census read 19 artifacts and every present file matched. A standalone clone records missing optional sources rather than failing. |

### `tables/_overview/`

Same-stem sources of the eight figures. `wtheatup_attribution.csv` carries the plotted counts and
labels over denominator 213. `lens_nes_by_contrast.csv` carries the NES values with `n_nominal`,
`n_effective` and the printed size/FDR string. `hsr_rank_position_panel.csv` carries every
plotted gene position plus the gate-span guard — `gate_pass_n` must be 213 and
`gate_deepest_gene` must be `Cpne6` at rank 2,010, which is what proves the plotted gate and the
frozen set are the same object. `gate_projection_bridge.csv` reconciles the two set sizes,
carries each side's lens membership, and records the downstream recovery census; read `block`
first to select the funnel, the lens rows or the per-compartment brackets. The four
`hsr_lens_membership_*` and `hsr_lens_provenance_euler` tables carry the fitted areas, their
residuals and the unassigned remainder for each drawing.
