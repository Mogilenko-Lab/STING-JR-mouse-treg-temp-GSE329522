#!/usr/bin/env python
"""
02_analysis/scripts/bump_dashboard.py
======================================

Thin producing script to run the interactive bump-chart dashboard generation pipeline.
"""

import sys
from pathlib import Path

# Set PYTHONPATH to include 01_modules/.ref
project_root = Path(__file__).resolve().parents[2]
ref_dir = project_root / "01_modules" / ".ref"
if str(ref_dir) not in sys.path:
    sys.path.insert(0, str(ref_dir))

from bump_dashboard import DashboardPipeline, DashboardConfig


def main() -> None:
    import pandas as pd
    print("Initializing DashboardConfig and DashboardPipeline...")
    cfg = DashboardConfig(
        config_yaml_path="02_analysis/config/analysis_config.yaml"
    )
    pipeline = DashboardPipeline(config=cfg)
    
    print("Running bump dashboard pipeline...")
    output_path = pipeline.run()
    print(f"Pipeline completed successfully. Dashboard generated at: {output_path}")

    print("Generating companion CSV...")
    from bump_dashboard.application.data_service import DashboardDataService
    ds = DashboardDataService(config=cfg)
    records, _ = ds.build()
    
    rows = []
    for r in records:
        rows.append({
            "pathway_id": r.pathway_id,
            "pathway_name": r.description,
            "database": r.database,
            "nes_WT_heat": r.wt.nes_heat,
            "padj_WT_heat": r.wt.padj_heat,
            "nes_KO_heat": r.ko.nes_heat,
            "padj_KO_heat": r.ko.padj_heat,
            "nes_interaction": r.nes_interaction,
            "padj_interaction": r.padj_interaction,
            "sig_interaction": r.sig_interaction,
            "pattern": r.pattern,
            "nes_temp_main": r.nes_temp_main,
            "nes_geno_main": r.nes_geno_main,
        })
    df = pd.DataFrame(rows)
    df = df.sort_values(by=["pathway_id", "database"])
    csv_path = Path(output_path).parent / "gsea_bump_interaction.csv"
    df.to_csv(csv_path, index=False)
    print(f"Companion CSV written to: {csv_path}")


if __name__ == "__main__":
    main()
