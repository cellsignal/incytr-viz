import pandas as pd
import numpy as np
from dash import html, dcc
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


from incytr_viz.util import *


def create_hist_figure(paths, has_tprs, has_prs, has_p_value):

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

    if has_tprs:
        fig.add_trace(
            go.Histogram(
                x=paths["tprs"],
                name="TPRS",
                **common_hist_params,
            ),
            row=plot_order[curr_idx][0],
            col=plot_order[curr_idx][1],
        )
        curr_idx += 1
    if has_prs:
        fig.add_trace(
            go.Histogram(
                x=paths["prs"],
                name="PRS",
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
                            "label": "data(id)",
                            "text-wrap": "ellipsis",
                            "text-valign": "top",
                            "text-halign": "right",
                            "font-size": "20px",
                            "height": "data(height)",
                            "width": "data(width)",
                            "background-color": "data(background_color)",
                        },
                    },
                    {
                        "selector": "edge",
                        "style": {
                            "curve-style": "unbundled-bezier",
                            "target-arrow-shape": "vee",
                            "arrow-scale": ".75",
                            "label": "data(label)" if show_network_weights else "",
                            "loop-sweep": "30deg",
                            "width": "data(width)",
                            "line-color": "data(line_color)",
                            "target-arrow-color": "data(line_color)",
                        },
                    },
                ],
                style={"width": "800px", "height": "800px"},
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
    title,
    group_id,
    warn,
    color_flow,
):

    def get_sankey_height(num_targets, num_links):
        if num_links == 0:
            out = 250
        elif num_targets == 0:
            out = 400
        elif num_targets < 50:
            out = num_targets * 25
        else:
            out = num_targets * 15
        return f"{max(out, 250)}px"

    num_targets = len([x for x in ids if "_target" in x])
    num_links = len(ids)
    sankey_style = {"height": get_sankey_height(num_targets, num_links)}

    warn_style = {"display": "none"} if not warn else {}
    celltype_legend_style = {"display": "none"} if not color_flow else {}

    # Drop duplicate rows based on row index
    unique_clusters = clusters.loc[~clusters.index.duplicated(keep="first")]

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
                                            html.Tr(
                                                [
                                                    html.Td(r[0]),
                                                    html.Td(
                                                        [],
                                                        style={
                                                            "backgroundColor": r[1][
                                                                "color"
                                                            ],
                                                            "width": "20px",
                                                        },
                                                    ),
                                                ],
                                                className="sankeyCellTypeLegendRow",
                                            )
                                            for r in unique_clusters.iterrows()
                                        ]
                                    ),
                                ],
                            ),
                        ],
                        className="sankeyCellTypeLegend",
                        style=celltype_legend_style,
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
                            hovertemplate="Node %{customdata} has total value %{value}<extra></extra>",
                            color=get_node_colors(ids),
                        ),
                        link=dict(
                            source=source,
                            target=target,
                            value=value,
                            color=color,
                            customdata=color,
                            hovertemplate="Link has value %{value} $%{customdata}<extra></extra>",
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


def sankey_legend_container() -> html.Div:
    return html.Div(
        [
            html.Span(
                [
                    "Ligand",
                    html.Div(
                        [],
                        style={
                            "background-color": "red",
                        },
                        className="sankeyColorLegendBox",
                    ),
                ],
                className="sankeyColorLegend",
            ),
            html.Span(
                [
                    "Receptor",
                    html.Div(
                        [],
                        style={
                            "background-color": "blue",
                        },
                        className="sankeyColorLegendBox",
                    ),
                ],
                className="sankeyColorLegend",
            ),
            html.Span(
                [
                    "EM",
                    html.Div(
                        [],
                        style={
                            "background-color": "green",
                        },
                        className="sankeyColorLegendBox",
                    ),
                ],
                className="sankeyColorLegend",
            ),
            html.Span(
                [
                    "Target",
                    html.Div(
                        [],
                        style={
                            "background-color": "purple",
                        },
                        className="sankeyColorLegendBox",
                    ),
                ],
                className="sankeyColorLegend",
            ),
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
                        options=["r_em", "r_t", "em_t"],
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
    id: str,
    minval: int,
    maxval: int,
    step: int,
    value: int,
    label: str,
    **slider_kwargs,
):

    tooltip_format = {
        "placement": "left",
        "always_visible": True,
    }

    class_name = "slider invertedSlider" if label.lower() == "sigprob" else "slider"

    marks = {0: "0", 1: "1"} if label.lower() in ["sigprob", "p-value"] else None

    return html.Div(
        [
            dcc.Slider(
                min=minval,
                max=maxval,
                step=step,
                value=value,
                marks=marks,
                tooltip=tooltip_format,
                id=id,
                className=class_name,
                **slider_kwargs,
            ),
            html.Span(label),
        ],
        className="sliderContainer",
    )


def range_slider(
    id: str,
    minval: int,
    maxval: int,
    step: int,
    value: list,
    label: str,
    disabled: bool,
):

    tooltip_format = {
        "placement": "left",
        "always_visible": True,
    }

    if label.lower() in ["tprs", "prs"]:
        marks = {-2: "-2", 2: "2"}
    else:
        marks = None

    return html.Div(
        [
            dcc.RangeSlider(
                min=minval,
                max=maxval,
                step=step,
                value=value,
                marks=marks,
                tooltip=tooltip_format,
                id=id,
                className="slider",
                disabled=disabled,
            ),
            html.Span(label),
        ],
        className="sliderContainer",
    )


def slider_container(
    has_tprs,
    has_prs,
    has_p_value,
):

    sliders = [
        slider(
            id={"type": "numerical-filter", "index": "sigprob"},
            minval=0,
            maxval=1,
            step=0.01,
            value=0.9,
            label="Sigprob",
            disabled=False,
        )
    ]
    sliders.append(
        slider(
            {"type": "numerical-filter", "index": "p-value"},
            minval=0,
            maxval=1,
            step=0.01,
            value=0.05,
            label="P-Value",
            disabled=not has_p_value,
        )
    )
    sliders.append(
        range_slider(
            {"type": "numerical-filter", "index": "tprs"},
            minval=-2,
            maxval=2,
            step=0.01,
            value=[-2, 2],
            label="TPRS",
            disabled=not has_tprs,
        )
    )
    sliders.append(
        range_slider(
            {"type": "numerical-filter", "index": "prs"},
            minval=-2,
            maxval=2,
            step=0.01,
            value=[-2, 2],
            label="PRS",
            disabled=not has_prs,
        )
    )
    return html.Div(
        [
            html.Div(sliders[0:2], className="sliderColumn"),
            html.Div(sliders[2:4], className="sliderColumn"),
        ],
        className="sidebarElement allSlidersContainer",
        id="allSlidersContainer",
    )
