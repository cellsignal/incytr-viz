import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc, html
from plotly.subplots import make_subplots

from incytr_viz.util import *


def create_hist_figure(paths, has_tpds, has_ppds, has_p_value):

    plot_order = [(1, 1), (1, 2), (2, 1), (2, 2)]
    curr_idx = 0

    common_hist_params = dict(
        nbinsx=100,
    )

    fig = make_subplots(2, 2)

    fig.add_trace(
        go.Histogram(x=paths["sigprob"], name="SigProb", **common_hist_params),
        row=plot_order[curr_idx][0],
        col=plot_order[curr_idx][1],
    )
    curr_idx += 1

    if has_tpds:
        fig.add_trace(
            go.Histogram(
                x=paths["tpds"],
                name="TPDS",
                **common_hist_params,
            ),
            row=plot_order[curr_idx][0],
            col=plot_order[curr_idx][1],
        )
        curr_idx += 1
    if has_ppds:
        fig.add_trace(
            go.Histogram(
                x=paths["ppds"],
                name="PPDS",
                **common_hist_params,
            ),
            row=plot_order[curr_idx][0],
            col=plot_order[curr_idx][1],
        )
        curr_idx += 1

    if has_p_value:
        fig.add_trace(
            go.Histogram(x=paths["p_value"], name="P-Value", **common_hist_params),
            row=plot_order[curr_idx][0],
            col=plot_order[curr_idx][1],
        )
        curr_idx += 1

    # Update layout for subplots
    fig.update_xaxes(title_text="Value")
    fig.update_yaxes(title_text="Count")

    # Create subplots grid
    fig.update_layout(
        # xaxis=dict(domain=[0, 0.5]),  # Adjust x-axis domain for the first two subplots
        # xaxis2=dict(
        #     domain=[0.5, 1]
        # ),  # Adjust x-axis domain for the second two subplots
        # yaxis=dict(domain=[0, 0.5]),  # Adjust y-axis domain for the first two subplots
        # yaxis2=dict(
        #     domain=[0.5, 1]
        # ),  # Adjust y-axis domain for the second two subplots
        showlegend=True,
    )

    return fig


def umap_graph(group_id, has_umap, all_pathways):

    if not has_umap:
        return None

    fig = px.scatter(
        all_pathways,
        x="umap1",
        y="umap2",
        color="afc",
        custom_data=["path"],
        color_continuous_scale=px.colors.diverging.Spectral[::-1],
    )
    scatter = dcc.Graph(
        id=f"umap-graph-{group_id}",
        figure=fig,
    )

    return scatter


def cytoscape_container(
    id,
    title,
    elements=[],
    show_network_weights=False,
    layout_name="circle",
):

    return html.Div(
        [
            cyto.Cytoscape(
                id=id,
                elements=elements,
                layout={"name": layout_name},
                minZoom=0.1,
                maxZoom=10,
                stylesheet=[
                    {
                        "selector": "node",
                        "style": {
                            "label": "data(label)",
                            "text-wrap": "ellipsis",
                            "text-valign": "top",
                            "text-halign": "right",
                            "font-size": "28px",
                            "height": "data(height)",
                            "width": "data(width)",
                            "backgroundColor": "data(background_color)",
                        },
                    },
                    {
                        "selector": "edge",
                        "style": {
                            "curve-style": "unbundled-bezier",
                            "target-arrow-shape": "vee",
                            "font-size": "28px",
                            "arrow-scale": ".75",
                            "label": "data(label)" if show_network_weights else "",
                            "loop-sweep": "30deg",
                            "width": "data(width)",
                            "line-color": "data(line_color)",
                            "target-arrow-color": "data(line_color)",
                        },
                    },
                ],
                style={"width": "100%", "height": "900px"},
            ),
        ],
        className="cytoscapeContainer",
    )


def sankey_container(
    clusters,
    ids,
    labels,
    source,
    target,
    value,
    color,
    group_id,
    color_flow,
):

    num_links = len(source)
    is_empty = num_links == 0
    num_targets = len([x for x in ids if "_target" in x])
    num_effectors = len([x for x in ids if "_em" in x])
    no_targets = num_targets == 0

    def get_sankey_height(is_empty, no_targets, num_targets, num_effectors):

        num_terminal_nodes = num_effectors if no_targets else num_targets
        if is_empty:
            out = 250
        elif num_terminal_nodes < 50:
            out = num_terminal_nodes * 25
        else:
            out = num_terminal_nodes * 15
        return f"{max(out, 250)}px"

    sankey_style = {
        "height": get_sankey_height(is_empty, no_targets, num_targets, num_effectors)
    }

    warn_style = {"display": "none"} if not no_targets else {}

    # Drop duplicate rows based on row index
    unique_clusters = clusters.loc[~clusters.index.duplicated(keep="first")]

    if num_links > 2000:
        return html.Div(
            [
                sankey_legend_container(),
                html.Div(
                    "Too many links to display. Please apply additional filters.",
                    className="linksWarning",
                ),
            ],
            className="sankeyTitleAndLegend",
        )

    else:

        return html.Div(
            [
                html.Div(
                    [
                        sankey_legend_container(),
                        html.Div(
                            [
                                html.H4("Cell type"),
                                html.Div(
                                    [
                                        dbc.Table(
                                            [
                                                html.Tbody(
                                                    [
                                                        html.Tr(
                                                            [
                                                                html.Td(r[0]),
                                                                html.Td(
                                                                    [],
                                                                    style={
                                                                        "backgroundColor": r[
                                                                            1
                                                                        ][
                                                                            "color"
                                                                        ],
                                                                        "width": "20px",
                                                                    },
                                                                ),
                                                            ],
                                                            className="sankeyLinkColorLegendRow",
                                                        )
                                                        for r in unique_clusters.iterrows()
                                                    ]
                                                )
                                            ]
                                        ),
                                    ],
                                ),
                            ],
                            className="sankeyLinkColorLegend",
                            style=(
                                {}
                                if (color_flow in ["sender", "receiver"])
                                else {"display": "none"}
                            ),
                        ),
                        html.Div(
                            [
                                html.H4("Kinase-Substrate Relationship"),
                                html.Div(
                                    [
                                        dbc.Table(
                                            [
                                                html.Tbody(
                                                    [
                                                        html.Tr(
                                                            [
                                                                html.Td(k),
                                                                html.Td(
                                                                    [],
                                                                    style={
                                                                        "backgroundColor": v,
                                                                        "width": "20px",
                                                                    },
                                                                ),
                                                            ],
                                                            className="sankeyLinkColorLegendRow",
                                                        )
                                                        for k, v in kinase_color_map().items()
                                                    ]
                                                )
                                            ]
                                        ),
                                    ],
                                ),
                            ],
                            className="sankeyLinkColorLegend",
                            style=(
                                {} if (color_flow == "kinase") else {"display": "none"}
                            ),
                        ),
                    ],
                    className="sankeyTitleAndLegend",
                ),
                dcc.Graph(
                    figure=go.Figure(
                        go.Sankey(
                            arrangement="fixed",
                            node=dict(
                                pad=15,
                                thickness=20,
                                line=dict(color="black", width=0.5),
                                label=labels,
                                customdata=ids,
                                hovertemplate="%{label}: %{value:.0f} pathways<extra></extra>",
                                color=get_node_colors(ids),
                            ),
                            link=dict(
                                source=source,
                                target=target,
                                value=value,
                                color=color,
                                customdata=color,
                                hovertemplate="%{source.customdata} --> %{target.customdata}: %{value:.0f} pathways<extra></extra>",
                            ),
                        ),
                    ),
                    id=f"sankey-{group_id}",
                    className="sankey",
                    style=sankey_style,
                ),
                html.Div(
                    [
                        dbc.Button(
                            "!",
                            id=f"sankey-warning-{group_id}",
                            color="white",
                            style=warn_style,
                            className="sankeyWarning",
                        ),
                        dbc.Popover(
                            dbc.PopoverBody(
                                "Too many target genes to display. To display targets, please apply additional filters"
                            ),
                            trigger="hover",
                            body=True,
                            target=f"sankey-warning-{group_id}",
                        ),
                    ],
                    className="sankeyWarningAndLegendContainer",
                ),
            ],
            className="sankeyContainer",
        )


def _sankey_legend(label, color):
    return html.Span(
        [
            label,
            html.Div(
                [],
                style={
                    "backgroundColor": color,
                },
                className="sankeyColorLegendBox",
            ),
        ],
        className="sankeyColorLegend",
    )


def sankey_legend_container() -> html.Div:
    return html.Div(
        [
            _sankey_legend("Ligand", "red"),
            _sankey_legend("Receptor", "blue"),
            _sankey_legend("EM", "green"),
            _sankey_legend("Target", "purple"),
        ],
        className="sankeyColorLegendsContainer",
    )


def filter_container(sender, receiver, em, target, ligand, receptor):

    all_molecules = list(set(em + target + ligand + receptor))
    return html.Div(
        children=[
            html.Div(
                [
                    dcc.Dropdown(
                        id="sender-select",
                        placeholder="Filter Senders",
                        multi=True,
                        clearable=True,
                        options=sender,
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="receiver-select",
                        placeholder="Filter Receivers",
                        multi=True,
                        clearable=True,
                        options=receiver,
                        className="filter",
                    ),
                ],
                className="filterColumn",
            ),
            html.Div(
                [
                    dcc.Dropdown(
                        id="ligand-select",
                        placeholder="Filter Ligands",
                        multi=True,
                        clearable=True,
                        options=ligand,
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="receptor-select",
                        placeholder="Filter Receptors",
                        multi=True,
                        clearable=True,
                        options=receptor,
                        className="filter",
                    ),
                ],
                className="filterColumn",
            ),
            html.Div(
                [
                    dcc.Dropdown(
                        id="em-select",
                        placeholder="Filter Effectors",
                        multi=True,
                        clearable=True,
                        options=em,
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="target-select",
                        placeholder="Filter Target Genes",
                        multi=True,
                        clearable=True,
                        options=target,
                        className="filter",
                    ),
                ],
                className="filterColumn",
            ),
            html.Div(
                [
                    dcc.Dropdown(
                        id="any-role-select",
                        placeholder="Filter Any Role",
                        multi=True,
                        clearable=True,
                        options=all_molecules,
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="kinase-select",
                        placeholder="Filter Kinase Interaction",
                        multi=False,
                        clearable=True,
                        options=["r_em", "r_t", "em_t", "em_r", "t_r", "t_em"],
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="umap-select-a",
                        disabled=False,
                        options=[],
                        className="filter",
                        style={"display": "none"},
                    ),
                    dcc.Dropdown(
                        id="umap-select-b",
                        disabled=False,
                        options=[],
                        className="filter",
                        style={"display": "none"},
                    ),
                ],
                className="filterColumn",
            ),
        ],
        className="filterContainer sidebarElement",
    )


def slider(
    label: str,
    range: bool = False,
    **slider_kwargs,
):

    component = dcc.Slider if not range else dcc.RangeSlider

    tooltip_format = slider_kwargs.pop("tooltip", default_slider_tooltip)

    return html.Div(
        [
            component(
                tooltip=tooltip_format,
                **slider_kwargs,
            ),
            html.Span(label),
        ],
        className="sliderContainer",
    )


def slider_container(
    has_tpds,
    has_ppds,
    has_p_value,
):

    pval_map = p_value_slider_map()

    sliders = [
        slider(
            "Sigprob",
            id={"type": "numerical-filter", "index": "sigprob"},
            min=0,
            max=1,
            step=0.01,
            value=filter_defaults()["sigprob"],
            marks={0: "0", 1: "1"},
            className="slider invertedSlider",
        ),
        slider(
            "P-Value",
            id={"type": "numerical-filter", "index": "p-value"},
            min=pval_map[0][0],
            max=pval_map[-1][0],
            step=1,
            value=filter_defaults()["p_value"],
            marks={str(x[0]): x[1] for x in pval_map},
            disabled=not has_p_value,
            tooltip={"style": {"display": "none"}},
            className="slider" if has_p_value else "slider disabledSlider",
        ),
        slider(
            "TPDS",
            range=True,
            id={"type": "numerical-filter", "index": "tpds"},
            min=-1.1,
            max=1.1,
            step=0.01,
            value=filter_defaults()["tpds"],
            disabled=not has_tpds,
            marks={-1.1: "-1.1", 1.1: "1.1"},
            allowCross=False,
            tooltip=(
                {"style": {"display": "none"}}
                if not has_tpds
                else default_slider_tooltip
            ),
            className="slider invertedSlider" if has_tpds else "slider disabledSlider",
        ),
        slider(
            "PPDS",
            range=True,
            id={"type": "numerical-filter", "index": "ppds"},
            min=-1.1,
            max=1.1,
            step=0.01,
            value=filter_defaults()["ppds"],
            disabled=not has_ppds,
            marks={-1.1: "-1.1", 1.1: "1.1"},
            allowCross=False,
            tooltip=(
                {"style": {"display": "none"}}
                if not has_tpds
                else default_slider_tooltip
            ),
            className="slider invertedSlider" if has_ppds else "slider disabledSlider",
        ),
    ]
    return html.Div(
        [
            html.Div(sliders[0:2], className="sliderColumn"),
            html.Div(sliders[2:4], className="sliderColumn"),
        ],
        className="sidebarElement allSlidersContainer",
        id="allSlidersContainer",
    )
