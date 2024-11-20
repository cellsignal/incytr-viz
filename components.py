import pandas as pd
from dash import html, dcc
import dash_cytoscape as cyto
import plotly.graph_objects as go

from util import *
from stylesheet import cytoscape_styles


def hist(df, col, title, num_bins=20):
    try:
        return html.Div(
            dcc.Graph(
                figure=px.histogram(
                    df,
                    x=col,
                    title=title,
                    nbins=num_bins,
                    histfunc="count",
                    width=300,
                    height=300,
                ),
            ),
            className="hist",
        )
    except:
        return html.Div(
            dcc.Graph(
                figure=px.histogram(
                    pd.DataFrame({title: []}), title=title, width=300, height=300
                )
            ),
            className="hist",
        )


def hist_container(group_id, *histograms):
    return html.Div(
        histograms,
        id=f"hist-container-{group_id}",
        className="histContainer",
    )


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
            ),
        ],
        className="cytoscapeContainer",
    )


def sankey_container(ids, labels, source, target, value, title, group_id):

    return html.Div(
        [
            html.H2(title),
            dcc.Graph(
                style={
                    "height": "400px",
                    "width": "100%",
                },
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
                        link=dict(source=source, target=target, value=value),
                    ),
                ),
                id=f"sankey-{group_id}",
            ),
        ],
        className="sankeyContainer",
    )


def radio_container() -> html.Div:
    return html.Div(
        dcc.RadioItems(
            [
                {
                    "label": html.Div(
                        ["Network View"],
                    ),
                    "value": "network",
                },
                {
                    "label": html.Div(
                        ["Pathways View"],
                    ),
                    "value": "sankey",
                },
            ],
            value="network",
            id="view-radio",
            labelClassName="radioLabel",
            className="radioContainer sidebarElement",
        ),
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
            dcc.Dropdown(
                id="sender-select",
                placeholder="Filter Senders",
                multi=True,
                clearable=True,
                options=pathways["sender"].unique(),
            ),
            dcc.Dropdown(
                id="receiver-select",
                placeholder="Filter Receivers",
                multi=True,
                clearable=True,
                options=pathways["receiver"].unique(),
            ),
            dcc.Dropdown(
                id="ligand-select",
                placeholder="Filter Ligands",
                multi=True,
                clearable=True,
                options=pathways["ligand"].unique(),
            ),
            dcc.Dropdown(
                id="receptor-select",
                placeholder="Filter Receptors",
                multi=True,
                clearable=True,
                options=pathways["receptor"].unique(),
            ),
            dcc.Dropdown(
                id="em-select",
                placeholder="Filter Effectors",
                multi=True,
                clearable=True,
                options=pathways["em"].unique(),
            ),
            dcc.Dropdown(
                id="target-select",
                placeholder="Filter Target Genes",
                multi=True,
                clearable=True,
                options=pathways["target"].unique(),
            ),
            dcc.Dropdown(
                id="any-role-select",
                placeholder="Filter All",
                multi=True,
                clearable=True,
                options=all_molecules,
            ),
            dcc.Dropdown(
                id="umap-select",
                disabled=False,
                options=[],
            ),
        ],
        id="filter-container",
    )


def slider_container(
    has_rna_score,
    has_final_score,
    has_p_value,
):

    def _slider(id: str, minval: int, maxval: int, step: int, value: int, label: str):

        tooltip_format = {
            "placement": "bottom",
            "always_visible": True,
        }

        return html.Div(
            [
                html.Div(
                    dcc.Slider(
                        min=minval,
                        max=maxval,
                        step=step,
                        value=value,
                        marks=None,
                        tooltip=tooltip_format,
                        id=id,
                    ),
                    className="slider",
                ),
                html.H4(label),
            ],
            className="sliderContainer",
        )

    def _range_slider(
        id: str, minval: int, maxval: int, step: int, value: list, label: str
    ):

        tooltip_format = {
            "placement": "bottom",
            "always_visible": True,
        }

        return html.Div(
            [
                html.Div(
                    dcc.RangeSlider(
                        min=minval,
                        max=maxval,
                        step=step,
                        value=value,
                        marks=None,
                        tooltip=tooltip_format,
                        id=id,
                    ),
                    className="slider",
                ),
                html.H4(label),
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
    return html.Div(sliders, className="sidebarElement", id="allSlidersContainer")


def sidebar(
    pathways: pd.DataFrame,
    has_rna_score: bool,
    has_final_score: bool,
    has_p_value: bool,
):

    return html.Div(
        [
            html.Div(
                children=[
                    html.Div(
                        [
                            filter_container(pathways),
                            radio_container(),
                        ],
                    ),
                    slider_container(
                        has_rna_score=has_rna_score,
                        has_final_score=has_final_score,
                        has_p_value=has_p_value,
                    ),
                ],
            ),
        ],
        className="sidebar",
    )


def group_component(group):
    return html.Div(
        [
            hist_container(group),
            html.Div([], id=f"figure-{group}-container"),
        ]
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
