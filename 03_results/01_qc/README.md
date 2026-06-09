# Phase 1 -- Sample-mapping QC verdict

**Generated:** 2026-06-08 19:35 by `02_analysis/scripts/01_mapping_qc.R`

## Verdict (PROVISIONAL -- pending collaborator sample sheet)

The inferred temperature-major mapping (021-025 WT_37, 026-030 cGASKO_37, 031-035 WT_39, 036-040 cGASKO_39) is **SUPPORTED (gate PASS)**: **20/20** libraries have a data-derived label (thermometer-predicted temperature + Cgas-predicted genotype) that matches the inferred label. No discordant libraries.

**Evidence.** The heat-shock thermometer (Hspa1b/Hsph1/Hspa1a/Dnajb1) is monotone with the inferred temperature: mean log2CPM = 4.91 in the hot half (031-040) vs 3.79 in the cool half (021-030), a +1.12 log2 shift. Cgas (Cgas) is higher in WT than cGAS-KO within both temperature halves (37C: WT-KO = 1.02; 39C: WT-KO = 1.24 log2CPM). PCA on the top 2000 variable genes places 78.9%/3.2% of variance on PC1/PC2; temperature and genotype track the leading axes, forming the expected 2x2.

**Scramble caveat.** The deposited CPM column order is temperature-major; the GEO GSM accessions are genotype-major. A naive positional GSM->column join therefore mislabels **10** libraries (026, 027, 028, 029, 030, 031, 032, 033, 034, 035) -- see `figures/fig1d_scramble.pdf`. This is why the mapping must be marker-derived, not accession-positional.

**Bottom line:** inferred mapping **SUPPORTED**, PROVISIONAL pending the collaborator sample sheet.

See `tables/mapping_verdict.csv` for the per-library machine-checkable verdict.

