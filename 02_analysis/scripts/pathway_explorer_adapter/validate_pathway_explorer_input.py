#!/usr/bin/env python3
"""validate_pathway_explorer_input.py — fail-fast schema check for master_de_table.csv.
Adapted from reference. Exit 0 = pass, 1 = any violation."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

REQUIRED_COLS = ["gene_symbol", "t", "logFC", "adj.P.Val", "contrast"]


def check_schema(df: pd.DataFrame) -> list[str]:
    errs = []
    for col in REQUIRED_COLS:
        if col not in df.columns:
            errs.append(f"missing column: {col}")
    if "gene_symbol" in df.columns and df["gene_symbol"].isna().any():
        errs.append(f"gene_symbol has {int(df['gene_symbol'].isna().sum())} NaN rows")
    if "t" in df.columns and df["t"].isna().all():
        errs.append("column t is entirely NaN")
    return errs


def main() -> int:
    master_de = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("03_results/master/master_de_table.csv")
    if not master_de.exists():
        print(f"SKIP: {master_de} missing — running-sum panel will be empty (not fatal)."); return 0
    df = pd.read_csv(master_de)
    errs = check_schema(df)
    if errs:
        print("Schema errors in master_de_table.csv:"); [print(f"  - {e}") for e in errs]; return 1
    print(f"OK schema {REQUIRED_COLS}: {len(df)} rows, {df['contrast'].nunique()} contrasts")
    try:
        from pathway_explorer.data_loader import load_gene_rankings
    except ImportError as e:
        print(f"FAIL: cannot import pathway_explorer.data_loader ({e}); pip install -e 01_modules/.ref/pathway-explorer"); return 1
    for contrast in sorted(df["contrast"].unique()):
        try:
            ranked = load_gene_rankings(master_de, contrast=contrast)
        except Exception as e:
            print(f"  FAIL {contrast}: {type(e).__name__}: {e}"); return 1
        if len(ranked) == 0:
            print(f"  FAIL {contrast}: empty ranking after filter"); return 1
        print(f"  OK {contrast}: {len(ranked)} genes, t-range [{ranked['t'].min():.2f}, {ranked['t'].max():.2f}]")
    print("All checks passed."); return 0


if __name__ == "__main__":
    sys.exit(main())
