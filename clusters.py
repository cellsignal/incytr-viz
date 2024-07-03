import random
import logging

from typing import Optional

import pandas as pd
import pdb

from dash import html, dcc
import dash_cytoscape as cyto
from dash.dependencies import Input, Output, State
from util import *

logger = logging.getLogger(__name__)


def apply_cytoscape_callbacks(app, full_pathways: pd.DataFrame, clusters: pd.DataFrame):
    @app.callback(
        Output("cytoscape-tapNodeData-output", "children"),
        Input("cytoscape-figure", "tapNodeData"),
    )
    def displayTapNodeData(data):
        if data:
            return "Cluster: " + data["label"] + "\nSize: " + str(data["cluster_size"])

    @app.callback(
        Output("cytoscape-figure", "elements"),
        Input("sender-select", "value"),
        Input("receiver-select", "value"),
        Input("ligand-select", "value"),
        Input("receptor-select", "value"),
        Input("em-select", "value"),
        Input("target-select", "value"),
        Input("threshold-slider", "value"),
        Input("direction-select", "value"),
    )
    def update_cytoscape(
        sender_select,
        receiver_select,
        ligand_select,
        receptor_select,
        em_select,
        target_select,
        threshold,
        direction_select,
    ):

        filtered_pathways = filter_pathways(
            full_pathways,
            filter_senders=sender_select,
            filter_receivers=receiver_select,
            filter_ligands=ligand_select,
            filter_receptors=receptor_select,
            filter_em=em_select,
            filter_target_genes=target_select,
            threshold=threshold,
            direction=direction_select,
        )
        nodes, edges = load_nodes(clusters), load_edges(filtered_pathways)
        print([e for e in edges if e["data"]["id"] == "MicrogliaMicroglia"])
        return nodes + edges


def random_color():
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    return f"rgb({r}, {g}, {b})"


def clean_clusters(df) -> pd.DataFrame:

    # TODO handle duplicates
    # duplicates = df[df.duplicated(subset="Type", keep=False)]
    return df


def node_size_map(cluster_count: int, total_count: int):
    """map cell count to node diameter. diameters will add to 50% of graph width"""
    proportion = (cluster_count / total_count) / 0.5
    return str(proportion * 100) + "%"


def edge_width_map(pathways: int, max_paths: int, max_width_px: int = 5):
    return str(pathways / max_paths * max_width_px) + "px"


def load_nodes(clusters: pd.DataFrame) -> dict:
    """{'cluster_name': count}

    Args:
        df (_type_): _description_

    Returns:
        dict: _description_
    """

    nodes = []

    # TODO clean clusters
    # clusters = clean_clusters(clusters)

    total_cells = clusters["Population"].sum()

    for _, s in clusters.iterrows():

        data, style = {}, {}
        data["id"] = s.Type
        data["label"] = s.Type
        data["cluster_size"] = s.Population
        style["background-color"] = random_color()
        style["width"] = node_size_map(s.Population, total_cells)
        style["height"] = node_size_map(s.Population, total_cells)
        nodes.append({"data": data, "style": style})

    return nodes


def load_edges(
    pathways: str,
    differential: bool = False,
):
    """add pathways from source to target"""
    edges = []

    pathways = pathways.copy()

    sr_pairs = pathways.groupby(["Sender.group", "Receiver.group"]).size().to_dict()
    for sr, num_pathways in sr_pairs.items():
        source_id, target_id = sr
        data = dict()
        data["id"] = source_id + target_id
        data["source"] = source_id
        data["target"] = target_id
        data["weight"] = num_pathways
        data["label"] = str(num_pathways)

        edges.append({"data": data})
    if len(edges) == 0:
        return edges
    # else:
    #     max_paths = max([e["data"]["weight"] for e in edges])
    #     for e in edges:
    #         e["style"]["width"] = edge_width_map(e["data"]["weight"], max_paths)
    return edges


def get_cytoscape_component(
    pathways: pd.DataFrame,
    clusters: pd.DataFrame,
    layout_name="circle",
    differential: bool = False,
):

    stylesheet = [
        {
            "selector": "node",
            "style": {
                "label": "data(id)",
                "font-size": "10vh",
                "text-wrap": "ellipsis",
                "text-valign": "center",
                "text-halign": "center",
                "font-size": "10px",
            },
        },
        {
            "selector": "edge",
            "style": {
                "curve-style": "bezier",
                "target-arrow-shape": "triangle",
                "arrow-scale": ".5",
                "label": "data(label)",
            },
        },
    ]

    # nodes, edges = load_nodes(clusters), load_edges(pathways, differential=differential)

    return (
        cyto.Cytoscape(
            id="cytoscape-figure",
            elements=[],
            style={"width": "100%", "height": "90vh"},
            layout={"name": layout_name},
            stylesheet=stylesheet,
        ),
    )
    # html.P(id="cytoscape-tapNodeData-output"),
    # html.P(id="cytoscape-tapEdgeData-output"),
