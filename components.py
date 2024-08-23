import pandas as pd
import pdb
from dash import html, dcc
import dash_cytoscape as cyto
import plotly.graph_objects as go
import numpy as np

from util import *
import incytr_stylesheets


def hist_container(container_style={}):

    return html.Div(
        [
            html.Div(dcc.Graph(id="sw-hist")),
            html.Div(dcc.Graph(id="rnas-hist")),
            html.Div(dcc.Graph(id="fs-hist")),
        ],
        style=container_style,
    )


def get_cytoscape_component(
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
                stylesheet=incytr_stylesheets.cytoscape,
            ),
        ],
        style={
            "width": "50%",
        },
    )


def pathways_df_to_sankey(
    sankey_df: pd.DataFrame,
    always_include_target_genes: bool = False,
) -> tuple:

    min_score, max_score = (
        sankey_df["final_score"].min(),
        sankey_df["final_score"].max(),
    )
    if np.isnan(min_score) and np.isnan(max_score):
        min_score, max_score = 0, 0

    def _get_values(
        df: pd.DataFrame, source_colname: str, target_colname: str
    ) -> pd.DataFrame:
        out = (
            df.groupby(source_colname)[target_colname]
            .value_counts()
            .reset_index(name="value")
        )
        out.rename(
            columns={source_colname: "Source", target_colname: "Target"}, inplace=True
        )
        out["source_id"] = out["Source"] + "_" + source_colname
        out["target_id"] = out["Target"] + "_" + target_colname

        return out

    l_r = _get_values(sankey_df, "Ligand", "Receptor")
    r_em = _get_values(sankey_df, "Receptor", "EM")
    em_t = _get_values(sankey_df, "EM", "Target")

    included_links = [l_r, r_em]

    ## auto-determine if target genes should be included
    def _should_display_targets() -> bool:
        num_targets = len(em_t["Target"].unique())

        return True if always_include_target_genes else num_targets <= 75

    if _should_display_targets():
        included_links.append(em_t)

    links = pd.concat(included_links, axis=0).reset_index(drop=True)
    # ids allow for repeating labels in ligand, receptor, etc. without pointing to same node
    ids = list(set(pd.concat([links["source_id"], links["target_id"]])))
    labels = [x.split("_")[0] for x in ids]
    source = [next(i for i, e in enumerate(ids) if e == x) for x in links["source_id"]]
    target = [next(i for i, e in enumerate(ids) if e == x) for x in links["target_id"]]
    value = links["value"]

    # direct_targets = links.apply(
    #     lambda x: list(set(links[links["source_id"] == x["source_id"]]["Target"])),
    #     axis=1,
    # )
    return (ids, labels, source, target, value)


def get_sankey_component(pathways, id, title):

    ids, labels, source, target, value = pathways_df_to_sankey(
        sankey_df=pathways, always_include_target_genes=False
    )

    # customdata = [ids, [", ".join(x) for x in direct_targets.values]]

    return html.Div(
        [
            html.H3(title),
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
                        link=dict(source=source, target=target, value=value),
                    ),
                ),
                id=id,
            ),
        ]
    )


def radio_container(container_style={}) -> html.Div:
    return html.Div(
        dcc.RadioItems(
            [
                {
                    "label": html.Div(
                        ["Network View"],
                        style={"fontSize": 20},
                    ),
                    "value": "network",
                },
                {
                    "label": html.Div(
                        ["Pathways View"],
                        style={"fontSize": 20},
                    ),
                    "value": "pathways",
                },
            ],
            value="network",
            id="view-radio",
            style=container_style,
        ),
    )


def filter_container(pathways, container_style={}):

    all_molecules = pd.concat(
        [
            pathways["Ligand"],
            pathways["Receptor"],
            pathways["EM"],
            pathways["Target"],
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
                options=pathways["Sender"].unique(),
            ),
            dcc.Dropdown(
                id="receiver-select",
                placeholder="Filter Receivers",
                multi=True,
                clearable=True,
                options=pathways["Receiver"].unique(),
            ),
            dcc.Dropdown(
                id="ligand-select",
                placeholder="Filter Ligands",
                multi=True,
                clearable=True,
                options=pathways["Ligand"].unique(),
            ),
            dcc.Dropdown(
                id="receptor-select",
                placeholder="Filter Receptors",
                multi=True,
                clearable=True,
                options=pathways["Receptor"].unique(),
            ),
            dcc.Dropdown(
                id="em-select",
                placeholder="Filter Effectors",
                multi=True,
                clearable=True,
                options=pathways["EM"].unique(),
            ),
            dcc.Dropdown(
                id="target-select",
                placeholder="Filter Target Genes",
                multi=True,
                clearable=True,
                options=pathways["Target"].unique(),
            ),
            dcc.Dropdown(
                id="all-molecules-select",
                placeholder="Filter All",
                multi=True,
                clearable=True,
                options=all_molecules,
            ),
        ],
        style=container_style,
    )


def slider_container(slider_container_style={}):

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
                    style={
                        "width": "70%",
                    },
                ),
                html.H4(label),
            ],
            style={
                "display": "flex",
                "alignItems": "center",
                "justifyContent": "space-between",
            },
        )

    return html.Div(
        children=[
            html.Div(
                [
                    _slider("sw-slider", 0, 1, 0.01, 0.7, "SigWeight"),
                    _slider("rnas-slider", 0, 1, 0.01, 0, "RNA Score (abs)"),
                    _slider("fs-slider", 0, 2, 0.01, 0, "Final Score (abs)"),
                ],
            ),
        ],
        style=slider_container_style,
    )


def pathway_filter_components(pathways: pd.DataFrame):

    return html.Div(
        [
            html.Div(
                children=[
                    html.Div(
                        [
                            filter_container(
                                pathways, container_style={"width": "200px"}
                            ),
                            radio_container(),
                        ],
                        style={
                            "display": "flex",
                            "flexDirection": "column",
                            "margin": "20px",
                        },
                    ),
                    slider_container(
                        {"width": "200px", "display": "flex", "flexDirection": "column"}
                    ),
                ],
            ),
        ],
        style={
            "display": "flex",
            "width": "25%",
        },
    )
