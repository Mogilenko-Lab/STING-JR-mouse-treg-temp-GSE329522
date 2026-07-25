import marimo

__generated_with = "0.23.14"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _():
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    import yaml
    from plotly.subplots import make_subplots

    return Path, go, make_subplots, pd, yaml


@app.cell(hide_code=True)
def _(Path):
    # Resolve the project root by walking up until the config sentinel appears, so
    # the app loads its tables whether launched from the repo root, the notebook
    # dir, or under the umbrella container.
    def _find_root():
        bases = []
        try:
            bases.append(Path(__file__).resolve())
        except NameError:
            pass
        bases.append(Path.cwd().resolve())
        for b in bases:
            p = b
            for _ in range(8):
                if (p / "02_analysis" / "config" / "analysis_config.yaml").exists():
                    return p
                p = p.parent
        raise FileNotFoundError("could not locate 02_analysis/config/analysis_config.yaml")

    PROJECT_ROOT = _find_root()
    RESULTS = PROJECT_ROOT / "03_results"
    CONFIG = PROJECT_ROOT / "02_analysis" / "config" / "analysis_config.yaml"
    return CONFIG, RESULTS


@app.cell(hide_code=True)
def _(CONFIG, RESULTS, pd, yaml):
    overview = RESULTS / "10_signature" / "tables" / "_overview"
    projection_overview = RESULTS / "11_projection" / "tables" / "_overview"
    decomp = RESULTS / "12_hsr_decomp" / "tables"

    signature_sizes = pd.read_csv(overview / "signature_sizes.csv")
    updown_overlap = pd.read_csv(overview / "updown_overlap.csv")
    ortholog_coverage = pd.read_csv(overview / "ortholog_coverage_preview.csv")
    human_signature_sizes = pd.read_csv(
        projection_overview / "human_signature_sizes.csv"
    )
    master_de_genes = pd.read_csv(RESULTS / "master" / "master_de_genes.csv")

    hsr_decomp_summary = pd.read_csv(decomp / "hsr_decomp_summary.csv")
    hsr_decomp_lens_nes = pd.read_csv(decomp / "hsr_decomp_lens_nes.csv")
    lens_nes_by_contrast = pd.read_csv(decomp / "lens_nes_by_contrast.csv")
    wtheatup_attribution = pd.read_csv(decomp / "wtheatup_attribution.csv")
    hsr_decomp_overlap = pd.read_csv(decomp / "hsr_decomp_overlap.csv")
    hsr_decomp_conditional = pd.read_csv(decomp / "hsr_decomp_conditional.csv")
    hsr_decomp_rank_concordance = pd.read_csv(
        decomp / "hsr_decomp_rank_concordance.csv"
    )

    with open(CONFIG) as _fh:
        cfg = yaml.safe_load(_fh)
    projection_decision = cfg.get("decisions", {}).get("projection", {})
    return (
        hsr_decomp_conditional,
        hsr_decomp_overlap,
        human_signature_sizes,
        lens_nes_by_contrast,
        master_de_genes,
        ortholog_coverage,
        projection_decision,
        signature_sizes,
        updown_overlap,
        wtheatup_attribution,
    )


@app.cell(hide_code=True)
def _():
    palette = {
        "green": "#009E73",
        "orange": "#E69F00",
        "burnt_orange": "#B35806",
        "sky": "#56B4E9",
        "blue": "#0072B2",
        "grey": "#8C8C8C",
        "light_grey": "#E6E6E6",
    }
    contrast_order = [
        "WT_heat",
        "KO_heat",
        "Temp_main",
        "Interaction",
        "Geno_at_39",
        "Geno_at_37",
        "Geno_main",
    ]
    attribution_order = ["thermal_HSR", "activation", "shared_both", "neither"]
    attribution_colors = {
        "thermal_HSR": palette["green"],
        "activation": palette["burnt_orange"],
        "shared_both": palette["sky"],
        "neither": palette["light_grey"],
    }
    return attribution_colors, attribution_order, contrast_order, palette


@app.cell(hide_code=True)
def _(mo, projection_decision):
    _export = ", ".join(projection_decision.get("contrasts_primary", []))
    _gate = projection_decision.get("gate", "fdr_logfc")
    mo.md(
        rf"""
        # Mouse 39 °C iTreg signature review

        Exploration for GSE329522 iTregs. 
        The pipline for GSE329522 includes standard DE, battery of GSEA, TF, PROGENy, GATOM, and CoReSh summaries for all seven contrasts. 
    
        The human compartments need one frozen, ortholog-mapped signature they can all score against.

       Here I double check what I am to freeze and make three calls:

        - which contrasts tell the story,
        - which significance gate to keep,
        - how to handle ortholog ambiguity.

        The current projection decision ships **`{_export}`** at gate **`{_gate}`**.
        """
    )
    return


@app.cell(hide_code=True)
def _(projection_decision):
    _gate = projection_decision.get("gate", "fdr_logfc")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## How big are the DE sets resulting from different contrasts run?

    How many genes actually survive each threshold. The two gates give very different answers.

    The design of the experiment is 2x2 factorial: genotype WT, cGASKO x temperature 37, 39, n=5/group.

    ### Factor reference levels (set in 00_setup_metadata.R): genotype=WT, temp=37.
    **groups:**
    * "WT_37"
    * "cGASKO_37"
    * "WT_39"
    * "cGASKO_39"

    **contrasts:**
    * Effect of temp in WT and KO
      * "WT_heat": WT_39 - WT_37
      * "KO_heat" KO_39 - KO_37
    * Effect of genotype in temp group
      * 'Geno_at_37': WT_37 - cGAS-KO_37
      * 'Geno_at_39': WT_39 - cGAS-KO_39
    * What genes are in the heat response are cGAS dependent?
      * 'Interaction': (WT_39 - cGASKO_39) - (WT_37 - cGASKO_37),
      i.e. genes upregulated are cGAS-dependent
    * What is the averaged effect of heat across genotypes
      * 0.5*(WT_39 + cGASKO_39) - 0.5*(WT_37 + cGASKO_37)

    The freeze keeps **`{_gate}`**. That keeps `WT_heat` lean enough to score
    as a focused set, while the full signed-`t` ranked lists remain available
    for ranking-based enrichment.
    """)
    return


@app.cell(hide_code=True)
def _(mo, projection_decision):
    _gate = projection_decision.get("gate", "fdr_logfc")

    mo.md(rf"""
    ## How big are the DE sets resulting from different contrasts run?

    How many genes actually survive each threshold. The two gates give very
    different answers.

    The design of the experiment is 2x2 factorial: genotype WT, cGASKO x temperature 37, 39, n=5/group.

    ### Factor reference levels (set in 00_setup_metadata.R): genotype=WT, temp=37.
      groups:
        * "WT_37"
        * "cGASKO_37"
        * "WT_39"
        * "cGASKO_39"
      contrasts:
          * Effect of temp in WT and KO
              * "WT_heat": WT_39 - WT_37
              * "KO_heat" KO_39 - KO_37
          * Effect of genotype in temp group 
              * 'Geno_at_37': WT_37 - cGAS-KO_37
              * 'Geno_at_39': WT_39 - cGAS-KO_39
          * What genes are in the heat response are cGAS dependent?
              * 'Interaction': (WT_39 - cGASKO_39) - (WT_37 - cGASKO_37),
              i.e. genes upregulated are cGAS-dependent 
          * What is the averaged effect of heat across genotypes 
              * 0.5*(WT_39 + cGASKO_39) - 0.5*(WT_37 + cGASKO_37)

    The freeze keeps **`{_gate}`**. That keeps `WT_heat` lean enough to score
    as a focused set, while the full signed-`t` ranked lists remain available
    for ranking-based enrichment.
    """)
    return


@app.cell(hide_code=True)
def _(go, make_subplots, mo, palette, signature_sizes):
    _gates = ["fdr_logfc", "fdr_only"]
    _fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=False,
        subplot_titles=["gate = fdr_logfc", "gate = fdr_only"],
        vertical_spacing=0.15,
    )
    for _row, _gate in enumerate(_gates, start=1):
        _d = signature_sizes[signature_sizes["gate"] == _gate].copy()
        for _direction, _color in [("up", palette["burnt_orange"]), ("down", palette["blue"])]:
            _dd = _d[_d["direction"] == _direction]
            _fig.add_trace(
                go.Bar(
                    x=_dd["contrast"],
                    y=_dd["n_genes"],
                    name=_direction if _row == 1 else None,
                    marker_color=_color,
                    text=_dd["n_genes"],
                    textposition="outside",
                    showlegend=_row == 1,
                ),
                row=_row,
                col=1,
            )
    _fig.update_layout(
        barmode="group",
        height=680,
        template="plotly_white",
        title=dict(text="Set sizes per contrast (mouse symbols, pre-ortholog)", x=0.5),
        font=dict(size=14),
        legend=dict(orientation="h", yanchor="bottom", y=1.03, x=0.5, xanchor="center"),
        margin=dict(l=40, r=20, t=100, b=80),
    )
    _fig.update_yaxes(title_text="genes in set")
    _fig.update_xaxes(tickangle=30)

    _wt = signature_sizes[
        (signature_sizes["contrast"] == "WT_heat")
        & (signature_sizes["gate"] == "fdr_logfc")
    ]
    _wt_up = int(_wt[_wt["direction"] == "up"]["n_genes"].iloc[0])
    _wt_down = int(_wt[_wt["direction"] == "down"]["n_genes"].iloc[0])
    _cap = mo.md(
        f"`WT_heat` is 4,207 up and 4,516 down at FDR-only. Adding |log2FC| >= 1 "
        f"cuts it to **{_wt_up} up / {_wt_down} down**. The smaller set is still "
        "large enough to score, and it focuses on genes that moved strongly."
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Is the Temp_main averaged temp effect same as the  the contrasts saying different things?

    If `Temp_main` is basically a copy of `WT_heat`, there is little point in
    shipping both. Jaccard overlap tells me how much the sets share.

    The design is a 2 × 2 factorial: genotype {WT, cGASKO} by temperature
    {37, 39}, n = 5 per group. `WT_heat` is the heat effect in WT. `KO_heat` is
    the heat effect in cGAS-KO samples. `Temp_main` is the pooled temperature
    effect across genotypes. `Interaction` is the cGAS-dependent heat effect.
    """)
    return


@app.cell(hide_code=True)
def _(contrast_order, go, make_subplots, mo, palette, updown_overlap):
    _gate = "fdr_logfc"
    _dirs = ["up", "down"]
    _fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=["UP sets", "DOWN sets"],
        horizontal_spacing=0.12,
    )
    for _col, _direction in enumerate(_dirs, start=1):
        _d = updown_overlap[
            (updown_overlap["gate"] == _gate)
            & (updown_overlap["direction"] == _direction)
        ].copy()
        _order = [c for c in contrast_order if c in set(_d["contrast_a"])]
        _z, _text = [], []
        for _a in _order:
            _row, _text_row = [], []
            for _b in _order:
                _m = _d[(_d["contrast_a"] == _a) & (_d["contrast_b"] == _b)]
                if len(_m):
                    _j = float(_m["jaccard"].iloc[0])
                    _row.append(_j)
                    _text_row.append(f"{_j:.2f}")
                else:
                    _row.append(None)
                    _text_row.append("")
            _z.append(_row)
            _text.append(_text_row)
        _fig.add_trace(
            go.Heatmap(
                z=_z,
                x=_order,
                y=_order,
                text=_text,
                texttemplate="%{text}",
                textfont=dict(size=12),
                colorscale=[[0, "#F7F7F7"], [1, palette["blue"]]],
                zmin=0,
                zmax=1,
                colorbar=dict(title="Jaccard") if _col == 2 else None,
                showscale=_col == 2,
                xgap=2,
                ygap=2,
            ),
            row=1,
            col=_col,
        )
    _fig.update_layout(
        height=520,
        template="plotly_white",
        title=dict(text="Set overlap across contrasts (gate = fdr_logfc)", x=0.5),
        font=dict(size=14),
        margin=dict(l=40, r=40, t=90, b=100),
    )
    _fig.update_xaxes(tickangle=90)
    _fig.update_yaxes(autorange="reversed")

    _cap = mo.md(
        "`WT_heat`, `KO_heat`, and `Temp_main` share a substantial heat-response "
        "core. `Interaction` is mostly separate and small. The committed overlap "
        "table is directional, so I read the up and down limbs separately here."
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Does it survive the trip to human?

    I need to count how many genes survive orthology mapping. The exported
    sets are mouse-derived, but the human compartments score human symbols.
    The mapping rule keeps one-to-one orthologs, unions one-mouse-to-many-human
    mappings, drops unmapped genes, and logs every loss.
    """)
    return


@app.cell(hide_code=True)
def _(
    go,
    human_signature_sizes,
    make_subplots,
    mo,
    ortholog_coverage,
    palette,
):
    _gates = ["fdr_logfc", "fdr_only"]
    _cats = [
        ("mapped_1to1", "1:1", palette["green"]),
        ("one_to_many", "one-to-many", palette["orange"]),
        ("unmapped", "no ortholog", palette["grey"]),
    ]
    _fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=False,
        subplot_titles=["gate = fdr_logfc", "gate = fdr_only"],
        vertical_spacing=0.15,
    )
    for _row, _gate in enumerate(_gates, start=1):
        _d = ortholog_coverage[ortholog_coverage["gate"] == _gate].copy()
        for _cat, _label, _color in _cats:
            _fig.add_trace(
                go.Bar(
                    x=_d["contrast"],
                    y=_d[_cat] / _d["n_input"],
                    name=_label if _row == 1 else None,
                    marker_color=_color,
                    showlegend=_row == 1,
                    customdata=_d[_cat],
                    hovertemplate="%{x}<br>"
                    + _label
                    + ": %{customdata}<br>fraction: %{y:.1%}<extra></extra>",
                ),
                row=_row,
                col=1,
            )
    _fig.update_layout(
        barmode="stack",
        height=640,
        template="plotly_white",
        title=dict(text="Ortholog coverage of the up+down sets", x=0.5),
        font=dict(size=14),
        legend=dict(orientation="h", yanchor="bottom", y=1.03, x=0.5, xanchor="center"),
        margin=dict(l=40, r=20, t=100, b=80),
    )
    _fig.update_yaxes(title_text="fraction of input genes", tickformat=".0%")
    _fig.update_xaxes(tickangle=30)

    _h = human_signature_sizes[human_signature_sizes["gate"] == "fdr_logfc"].copy()
    _rows = []
    for _, _r in _h.iterrows():
        _rows.append(
            f"| {_r['contrast']} | {_r['role']} | {_r['direction']} | "
            f"{int(_r['n_human'])} |"
        )
    _tbl = (
        "| contrast | role | direction | human genes |\n"
        "|---|---|---|---|\n" + "\n".join(_rows)
    )
    _cap = mo.md(
        "The major contrasts keep enough one-to-one human orthologs to score. "
        "`Interaction` remains a small set after mapping, so I carry it as an "
        "exploratory contrast rather than the headline.\n\n" + _tbl
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Does `WT_heat` look like heat-stressed Tregs?

    What are the genes sorted by signed `t`: noise or biology?

    This is the hinge. The top genes tell me whether `WT_heat` reads as a
    clean heat-shock list or as a broader iTreg response at 39 °C.
    """)
    return


@app.cell(hide_code=True)
def _(go, master_de_genes, mo):
    _wt = master_de_genes[master_de_genes["contrast"] == "WT_heat"].copy()
    _top_up = _wt.sort_values("t", ascending=False).head(30).reset_index(drop=True)
    _top_down = _wt.sort_values("t", ascending=True).head(30).reset_index(drop=True)
    _fig = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=["rank", "up gene", "up t", "down gene", "down t"],
                    fill_color="#F2F2F2",
                    align="left",
                    font=dict(size=14),
                ),
                cells=dict(
                    values=[
                        list(range(1, 31)),
                        _top_up["gene_symbol"],
                        [f"{x:.2f}" for x in _top_up["t"]],
                        _top_down["gene_symbol"],
                        [f"{x:.2f}" for x in _top_down["t"]],
                    ],
                    align="left",
                    font=dict(size=13),
                    height=26,
                ),
            )
        ]
    )
    _fig.update_layout(
        height=880,
        template="plotly_white",
        title=dict(text="WT_heat top-30 up and down by signed t", x=0.5),
        margin=dict(l=10, r=10, t=60, b=10),
    )
    _cap = mo.md(
        "The list reads as a general stress-protective response with migratory and "
        "activation hints. The up side includes chaperone/proteostasis genes such as "
        "`Hsph1`, `Hspa4l`, and `Hspa1a`, but it also carries complement, "
        "immediate-early, matrix, cytoskeletal, and differentiation signals. That mix "
        "sets up the next question: how much of `WT_heat_up` is heat, and how much is "
        "activation at 39 °C?"
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Is `WT_heat_up` heat, or just activation at 39 °C?

    `WT_heat_up` is iTregs activated at 39 °C — its top genes are part chaperone (Hspa1a, Hsph1),
    part activation (complement C3, immediate-early). The list can't separate the two. So I read it
    through two lenses that can: a clean heat-shock core (HSF1, chaperones) and a separate activation
    lens, both scored on the same ranking.

    They disagree, and the disagreement is the answer.

    Count the genes: it looks like activation. Of 213 `WT_heat_up` genes, 3 are chaperones, 12 are
    activation. The gate keeps the sharpest movers, and activation moves sharpest.

    Weigh the whole ranking: heat wins. `HSR_core` enriches at NES 2.05, activation at 1.72. The
    chaperone program is broad and coordinated — it runs the length of the ranking, top to bottom.

    So the list is activation; the program underneath is heat. A handful of activation genes clear the
    bar loudly, while a wide chaperone program hums below it.

    Two checks confirm two independent signals. They share no genes. And dropping every activation
    gene leaves `HSR_core` exactly where it was (2.05 → 2.05).

    The heat program doesn't need cGAS either — it's just as strong in `KO_heat` (2.07) as in
    `WT_heat`. That points to HSF1 as the driver, the chaperone response, working without the DNA sensor.

    The 39 °C response carries a real heat-shock program, stronger than the activation it's tangled
    with, even though the gene list alone reads as activated. That's what the human data now tests.
    The core answers proteotoxic stress of any kind — oxidative, proteasomal, thermal — so "heat" here
    means the program 39 °C induces, measured by the 37/39 contrast.
    """)
    return


@app.cell(hide_code=True)
def _(attribution_colors, attribution_order, go, mo, wtheatup_attribution):
    _d = (
        wtheatup_attribution[
            wtheatup_attribution["attribution"].isin(attribution_order)
        ]
        .set_index("attribution")
        .loc[attribution_order]
        .reset_index()
    )
    _fig = go.Figure()
    for _, _r in _d.iterrows():
        _fig.add_trace(
            go.Bar(
                x=["WT_heat_up"],
                y=[int(_r["n"])],
                name=str(_r["attribution_label"]),
                marker_color=attribution_colors[_r["attribution"]],
                text=[str(_r["label"])],
                textposition="inside" if int(_r["n"]) >= 20 else "outside",
                hovertemplate="%{fullData.name}: %{y} genes<extra></extra>",
            )
        )
    _fig.update_layout(
        barmode="stack",
        height=430,
        template="plotly_white",
        title=dict(text="WT_heat_up membership attribution", x=0.5),
        yaxis_title="genes",
        xaxis_title="",
        font=dict(size=14),
        legend=dict(orientation="h", yanchor="bottom", y=1.03, x=0.5, xanchor="center"),
        margin=dict(l=40, r=20, t=90, b=40),
    )
    _fig.update_yaxes(range=[0, int(_d["denominator"].iloc[0]) * 1.12])
    _cap = mo.md(
        "Counting genes makes the list look activation-heavy: 12 activation genes, "
        "3 thermal-HSR genes, no shared genes, and 198 genes in neither curated lens."
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell(hide_code=True)
def _(lens_nes_by_contrast, mo):
    _contrasts = [
        c
        for c in ["WT_heat", "Temp_main", "KO_heat"]
        if c in set(lens_nes_by_contrast["contrast"])
    ]
    if not _contrasts:
        _contrasts = sorted(lens_nes_by_contrast["contrast"].unique().tolist())
    heat_contrast = mo.ui.dropdown(
        options=_contrasts,
        value="WT_heat" if "WT_heat" in _contrasts else _contrasts[0],
        label="Heat contrast",
    )
    heat_contrast
    return (heat_contrast,)


@app.cell(hide_code=True)
def _(go, heat_contrast, lens_nes_by_contrast, mo, palette):
    _selected = heat_contrast.value
    _terms = ["HSR_core", "TCR_activation"]
    _labels = {"HSR_core": "HSR core", "TCR_activation": "activation"}
    _colors = {"HSR_core": palette["green"], "TCR_activation": palette["burnt_orange"]}
    _d = (
        lens_nes_by_contrast[
            (lens_nes_by_contrast["contrast"] == _selected)
            & (lens_nes_by_contrast["term"].isin(_terms))
        ]
        .set_index("term")
        .loc[_terms]
        .reset_index()
    )

    _fig = go.Figure()
    for _, _r in _d.iterrows():
        _fig.add_trace(
            go.Bar(
                x=[_selected],
                y=[float(_r["nes"])],
                name=_labels[_r["term"]],
                marker_color=_colors[_r["term"]],
                text=[f"{float(_r['nes']):.2f}"],
                textposition="outside",
                customdata=[[float(_r["padj"]), int(_r["set_size"])]],
                hovertemplate="%{fullData.name}<br>NES %{y:.2f}<br>FDR %{customdata[0]:.2g}<br>set size %{customdata[1]}<extra></extra>",
            )
        )
    _fig.add_hline(y=0, line_width=1, line_color="#444444")
    _fig.update_layout(
        barmode="group",
        height=420,
        template="plotly_white",
        title=dict(text=f"Lens NES for {_selected}", x=0.5),
        yaxis_title="NES",
        xaxis_title="",
        font=dict(size=14),
        legend=dict(orientation="h", yanchor="bottom", y=1.03, x=0.5, xanchor="center"),
        margin=dict(l=40, r=20, t=90, b=40),
    )
    _fig.update_yaxes(range=[0, max(2.4, float(_d["nes"].max()) + 0.25)])

    _rows = []
    for _, _r in _d.iterrows():
        _rows.append(
            f"| {_labels[_r['term']]} | {float(_r['nes']):.2f} | "
            f"{float(_r['padj']):.2g} | {int(_r['set_size'])} | "
            f"{int(_r['leading_edge_n'])} |"
        )
    _tbl = (
        f"**{_selected}.** The whole ranking favors the heat-shock lens over the "
        "activation lens.\n\n"
        "| lens | NES | FDR | set size | leading edge |\n"
        "|---|---|---|---|---|\n"
        + "\n".join(_rows)
    )
    mo.hstack([mo.ui.plotly(_fig), mo.md(_tbl)])
    return


@app.cell(hide_code=True)
def _(hsr_decomp_conditional, hsr_decomp_overlap, mo):
    _ov = hsr_decomp_overlap[
        (hsr_decomp_overlap["set_a"] == "HSR_core")
        & (hsr_decomp_overlap["set_b"] == "TCR_activation")
    ].iloc[0]
    _cond = hsr_decomp_conditional[
        (hsr_decomp_conditional["term"] == "HSR_core")
        & (hsr_decomp_conditional["conditioned_on"] == "TCR_activation")
    ].iloc[0]
    mo.md(
        f"""
        **Reading the panels.** The membership panel and the ranking panel answer
        different questions. Counting genes says the thresholded `WT_heat_up` list is
        activation-heavy. Ranking all genes says `HSR_core` is the stronger coordinated
        program. The independence checks match that reading: `HSR_core` and
        activation share {int(_ov['n_intersect'])} genes, and conditioning on
        activation leaves `HSR_core` at {float(_cond['nes_uncond']):.2f} →
        {float(_cond['nes_cond']):.2f}.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo, projection_decision):
    _orth = projection_decision.get("ortholog_ambiguity", {})
    mo.md(
        rf"""
        ## Decision — what freezes

        Pulling it together, this is what the signature freeze commits to.

        1. **Contrasts.** I keep `WT_heat` as the primary contrast and carry `KO_heat`
        plus `Interaction` for the temperature and cGAS-dependence questions.
        2. **Gate.** I keep `{projection_decision.get("gate", "fdr_logfc")}`.
        3. **Ortholog ambiguity.** I keep one-to-many as
        `{_orth.get("one_mouse_to_many_human")}`, many-to-one as
        `{_orth.get("many_mouse_to_one_human")}`, and genes with no human ortholog as
        `{_orth.get("no_human_ortholog")}`. Every dropped gene stays auditable in the
        mapping table.

        The projection status is **`{projection_decision.get("status")}`**. The ranked
        files feed human donor-level pseudobulk enrichment. The up/down text sets feed
        per-cell scoring where that is the available evidence tier. I keep `WT_heat`
        as the headline and use the other sets to ask narrower follow-up questions.
        """
    )
    return


if __name__ == "__main__":
    app.run()
