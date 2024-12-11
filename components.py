import pandas as pd
from dash import html, dcc
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import plotly.graph_objects as go
import plotly.express as px

from util import *
from stylesheet import cytoscape_styles


def hist(df, col, title, num_bins=20):
    try:
        return dcc.Graph(
            figure=px.histogram(
                df,
                x=col,
                title=title,
                nbins=num_bins,
                histfunc="count",
            ),
            className="hist",
        )

    except:
        return dcc.Graph(
            figure=px.histogram(
                pd.DataFrame({title: []}),
                title=title,
            ),
            className="hist",
        )


def hist_container(group_id, *histograms):
    return html.Div(
        children=list(histograms),
        id=f"hist-{group_id}-container",
        className="histContainer",
    )


def umap_container(group_id, group_name, has_umap, all_pathways):
    if has_umap:
        fig = px.scatter(
            all_pathways,
            x="umap1",
            y="umap2",
            color=f"sigweight",
            custom_data=["path"],
        )
        scatter = dcc.Graph(
            id=f"scatter-plot-{group_id}",
            figure=fig,
            style={"width": "600px", "height": "600px"},
        )

        return html.Div([scatter], className="umapContainer")

    else:
        return html.Div([], style={"display": "none"})


def cytoscape_container(
    id,
    title,
    elements=[],
    layout_name="circle",
):

    return html.Div(
        [
            html.H3(title),
            cyto.Cytoscape(
                id=id,
                elements=elements,
                layout={"name": layout_name},
                stylesheet=cytoscape_styles,
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
    show_celltype_legend,
):

    warn_style = {"display": "none"} if not warn else {}
    celltype_legend_style = {"display": "none"} if not show_celltype_legend else {}

    # pdb.set_trace()

    # Drop duplicate rows based on row index
    unique_clusters = clusters.loc[~clusters.index.duplicated(keep="first")]
    return html.Div(
        [
            html.Div(
                [html.H2(title), sankey_legend_container()],
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
                        # legend="legend2",
                    ),
                    # layout=dict(
                    #     legend2=dict(
                    #         title=dict(text="By continent"),
                    #         xref="container",
                    #         yref="container",
                    #         y=0.85,
                    #         bgcolor="Gold",
                    #     )
                    # ),
                ),
                id=f"sankey-{group_id}",
                className="sankey",
            ),
            html.Div(
                [
                    html.H4("cell type"),
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
                                                    "backgroundColor": r[1]["color"],
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
                        style=celltype_legend_style,
                        className="sankeyCellTypeLegend",
                    ),
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


def filter_container(pathways):

    all_molecules = pd.concat(
        [
            pathways["ligand"],
            pathways["receptor"],
            pathways["em"],
            pathways["target"],
        ],
        axis=0,
    ).unique()

    return html.Div(
        children=[
            html.Div(
                [
                    dcc.Dropdown(
                        id="sender-select",
                        placeholder="Filter Senders",
                        multi=True,
                        clearable=True,
                        options=pathways["sender"].unique(),
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="receiver-select",
                        placeholder="Filter Receivers",
                        multi=True,
                        clearable=True,
                        options=pathways["receiver"].unique(),
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
                        options=pathways["ligand"].unique(),
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="receptor-select",
                        placeholder="Filter Receptors",
                        multi=True,
                        clearable=True,
                        options=pathways["receptor"].unique(),
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
                        options=pathways["em"].unique(),
                        className="filter",
                    ),
                    dcc.Dropdown(
                        id="target-select",
                        placeholder="Filter Target Genes",
                        multi=True,
                        clearable=True,
                        options=pathways["target"].unique(),
                        className="filter",
                    ),
                ],
                className="filterColumn",
            ),
            html.Div(
                [
                    dcc.Dropdown(
                        id="any-role-select",
                        placeholder="Filter All",
                        multi=True,
                        clearable=True,
                        options=all_molecules,
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
        id="filter-container",
        className="filterContainer sidebarElement",
    )


def slider_container(
    has_rna_score,
    has_final_score,
    has_p_value,
):

    def _slider(id: str, minval: int, maxval: int, step: int, value: int, label: str):

        tooltip_format = {
            "placement": "left",
            "always_visible": True,
        }

        return html.Div(
            [
                dcc.Slider(
                    min=minval,
                    max=maxval,
                    step=step,
                    value=value,
                    marks=None,
                    tooltip=tooltip_format,
                    id=id,
                    className="slider",
                ),
                html.Span(label),
            ],
            className="sliderContainer",
        )

    def _range_slider(
        id: str, minval: int, maxval: int, step: int, value: list, label: str
    ):

        tooltip_format = {
            "placement": "left",
            "always_visible": True,
        }

        return html.Div(
            [
                dcc.RangeSlider(
                    min=minval,
                    max=maxval,
                    step=step,
                    value=value,
                    marks=None,
                    tooltip=tooltip_format,
                    id=id,
                    className="slider",
                ),
                html.Span(label),
            ],
            className="sliderContainer",
        )

    sliders = [
        _slider(
            {"type": "numerical-filter", "index": "sigweight"},
            0,
            1,
            0.01,
            0.7,
            "SigWeight",
        )
    ]
    if has_p_value:
        sliders.append(
            _slider(
                {"type": "numerical-filter", "index": "p-value"},
                minval=0,
                maxval=1,
                step=0.01,
                value=1,
                label="P-Value",
            )
        )
    if has_rna_score:
        sliders.append(
            _range_slider(
                {"type": "numerical-filter", "index": "rna-score"},
                minval=-2,
                maxval=2,
                step=0.01,
                value=[-2, 2],
                label="RNA Score",
            )
        )
    if has_final_score:
        sliders.append(
            _range_slider(
                {"type": "numerical-filter", "index": "final-score"},
                minval=-2,
                maxval=2,
                step=0.01,
                value=[-2, 2],
                label="Final Score",
            )
        )
    return html.Div(
        sliders,
        className="sidebarElement allSlidersContainer",
        id="allSlidersContainer",
    )


# import dash_bootstrap_components as dbc
# from dash import Input, Output, State, html

# collapse = html.Div(
#     [
#         ,
#         dbc.Collapse(
#             dbc.Card(dbc.CardBody("This content is hidden in the collapse")),
#             id="collapse",
#             is_open=False,
#         ),
#     ]
# )


# @app.callback(
#     Output("collapse", "is_open"),
#     [Input("collapse-button", "n_clicks")],
#     [State("collapse", "is_open")],
# )
# def toggle_collapse(n, is_open):
#     if n:
#         return not is_open
#     return is_open
