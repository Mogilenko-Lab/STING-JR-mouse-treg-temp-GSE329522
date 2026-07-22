# GSE329522 — sample provenance (iTreg 2×2 genotype × temperature)

Bulk RNA-seq of induced regulatory T cells (iTreg) differentiated from primary murine
splenic CD4⁺ T cells, in a 2×2 design: **genotype** (WT vs cGAS-KO) × **temperature**
(37 °C vs 39 °C), **5 biological replicates per group**, 20 libraries total.

Processed matrix: `GSE329522_normalized_counts_CPM_iTreg.csv` (genes × 20 libraries,
columns `12630-RS-021 … 12630-RS-040`, deposited in **temperature-major** order).

## Sample assignment — owner-confirmed 2026-07-22

The library → condition mapping below is the **authoritative** assignment, provided by
the data owner (sample sheet, 2026-07-22). It matched the prior Hspa1b/Hsph1 heat-shock
thermometer + Cgas-expression inference **exactly (20/20 concordant)**, so `mapping_status`
is now `CONFIRMED` (previously `INFERRED`). Built by `02_analysis/scripts/00_setup_metadata.R`.

| Library | Genotype | Treatment | Group |
|---|---|---|---|
| 12630-RS-021 | WT | 37 °C | WT_37 |
| 12630-RS-022 | WT | 37 °C | WT_37 |
| 12630-RS-023 | WT | 37 °C | WT_37 |
| 12630-RS-024 | WT | 37 °C | WT_37 |
| 12630-RS-025 | WT | 37 °C | WT_37 |
| 12630-RS-026 | cGAS-KO | 37 °C | cGASKO_37 |
| 12630-RS-027 | cGAS-KO | 37 °C | cGASKO_37 |
| 12630-RS-028 | cGAS-KO | 37 °C | cGASKO_37 |
| 12630-RS-029 | cGAS-KO | 37 °C | cGASKO_37 |
| 12630-RS-030 | cGAS-KO | 37 °C | cGASKO_37 |
| 12630-RS-031 | WT | 39 °C | WT_39 |
| 12630-RS-032 | WT | 39 °C | WT_39 |
| 12630-RS-033 | WT | 39 °C | WT_39 |
| 12630-RS-034 | WT | 39 °C | WT_39 |
| 12630-RS-035 | WT | 39 °C | WT_39 |
| 12630-RS-036 | cGAS-KO | 39 °C | cGASKO_39 |
| 12630-RS-037 | cGAS-KO | 39 °C | cGASKO_39 |
| 12630-RS-038 | cGAS-KO | 39 °C | cGASKO_39 |
| 12630-RS-039 | cGAS-KO | 39 °C | cGASKO_39 |
| 12630-RS-040 | cGAS-KO | 39 °C | cGASKO_39 |

Reference levels for modelling: genotype **WT**, temperature **37 °C**. The DE model
(`02_de_limma_trend.R`) is `~0 + group`; contrasts derive from these groups, e.g.
`WT_heat = WT_39 − WT_37`, `Interaction = (WT_39 − cGASKO_39) − (WT_37 − cGASKO_37)`.

## GSM accession caveat

GEO deposits the GSM accessions in **genotype-major** order, whereas the CPM columns are
**temperature-major**. The true per-library GSM therefore differs from a naive positional
join for libraries 026–035. `00_setup_metadata.R` carries both `gsm_id` (correct, matched
by genotype+temperature) and `gsm_id_positional` (the wrong naive join, kept to document
the trap). The analysis is keyed on `library_id`, so this accession ordering does not
affect any result — it matters only for GEO cross-referencing.

| Group | Libraries | True GSM |
|---|---|---|
| WT_37 | 021–025 | GSM9705690–694 |
| cGASKO_37 | 026–030 | GSM9705700–704 |
| WT_39 | 031–035 | GSM9705695–699 |
| cGASKO_39 | 036–040 | GSM9705705–709 |

## Wet-lab protocols (owner-provided)

**Cell type:** purified splenic CD4⁺ T cells, differentiated to iTreg.

**Growth / differentiation.** Primary murine CD4 T cells were isolated from mouse spleens
by negative selection (STEMCELL Technologies or Miltenyi Biotec) and cultured (5% CO₂,
RPMI 1640 with 10 mM HEPES, 50 µM 2-ME, penicillin-streptomycin 100 U/ml, 2 mM glutamine).
Cells were activated on plate-bound anti-CD3 (3 µg/ml) + anti-CD28 (2 µg/ml) at
0.75–1.0 × 10⁶ cells/well (24-well) and cultured 3–4 days with iTreg-polarizing factors:
TGF-β (1.5 ng/ml), IL-2 (100 U/ml), anti-IL-4 (10 µg/ml), anti-IFN-γ (10 µg/ml).

**Temperature treatment.** All cells were cultured at **37 °C or 39 °C for 4 days** prior
to analysis (the sole experimental temperature variable).

**Extraction / library prep.** Cells washed in cold PBS, pellets flash-frozen; RNA
extracted on the Autogen Xtract 16+ (triXact RNA Kit) with on-column DNase. Libraries
from 120–200 ng total RNA with the NEBNext rRNA-Depletion Kit (RNaseH-based, cytoplasmic +
mitochondrial rRNA), poly-A enrichment (oligo-dT), thermal fragmentation, cDNA synthesis,
adaptor ligation, PCR. Extracted molecule: poly-A RNA.

## Downstream validity

Because the confirmed mapping is identical to the mapping the pipeline already used, all
prior results (DE contrasts, the WT_heat / KO_heat / Interaction / genotype signatures, the
CoReSh compendium search) are wired on the correct sample→condition assignment and require
no recomputation. The change is one of status only: `INFERRED → CONFIRMED`.
