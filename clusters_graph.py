import os
import pdb
import math
import random
import logging

from typing import Optional

import numpy as np
import pandas as pd

from dash import Dash, html
import dash_cytoscape as cyto
from dash.dependencies import Input, Output

import dtypes

logger = logging.getLogger(__name__)


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


# def get_pathway_files(paths_dir: str):
#     """get all csvs in a directory"""
#     return [
#         os.path.join(paths_dir, obj)
#         for obj in os.listdir(paths_dir)
#         if obj.endswith(
#             (
#                 ".csv",
#                 ".tsv",
#             )
#         )
#         and os.path.isfile(os.path.join(paths_dir, obj))
#     ]


def load_nodes(clusters_filepath) -> dict:
    """{'cluster_name': count}

    Args:
        df (_type_): _description_

    Returns:
        dict: _description_
    """

    nodes = []

    clusters = pd.read_csv(
        clusters_filepath, dtype={"Type": str, "Population": int}
    ).reset_index(drop=True)

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
    pathways_file: str,
    direction: Optional[str] = None,
    threshold: Optional[float] = None,
):
    """add pathways from source to target"""

    edges = []

    df = pd.read_csv(pathways_file)

    if direction == "up":
        df = df[df["final_score"] > 0]
    elif direction == "down":
        df = df[df["final_score"] < 0]

    if threshold:
        df = df[df["final_score"].abs() > threshold]

    sr_pairs = df.groupby(["Sender.group", "Receiver.group"]).size().to_dict()
    for sr, num_pathways in sr_pairs.items():
        source_id, target_id = sr
        try:
            edge_index = next(
                i
                for i, e in enumerate(edges)
                if e["data"]["source"] == source_id and e["data"]["target"] == target_id
            )
            edges[edge_index]["data"]["weight"] += num_pathways
        except StopIteration:
            data = dict()
            data["id"] = source_id + target_id
            data["source"] = source_id
            data["target"] = target_id
            data["weight"] = num_pathways

            edges.append({"data": data, "style": dict()})

    max_paths = max([e["data"]["weight"] for e in edges])
    for e in edges:
        e["style"]["width"] = edge_width_map(e["data"]["weight"], max_paths)
    return edges


def incytr_cytoscape_app(nodes, edges, layout_name="circle"):

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
            },
        },
    ]

    app = Dash(__name__)
    app.layout = html.Div(
        [
            cyto.Cytoscape(
                id="cytoscape-clusters-graph",
                elements=nodes + edges,
                style={"width": "100%", "height": "90vh"},
                layout={"name": layout_name},
                stylesheet=stylesheet,
            ),
            html.P(id="cytoscape-tapNodeData-output"),
            html.P(id="cytoscape-tapEdgeData-output"),
        ]
    )

    @app.callback(
        Output("cytoscape-tapNodeData-output", "children"),
        Input("cytoscape-clusters-graph", "tapNodeData"),
    )
    def displayTapNodeData(data):
        if data:
            return "Cluster: " + data["label"] + "\nSize: " + str(data["cluster_size"])

    @app.callback(
        Output("cytoscape-tapEdgeData-output", "children"),
        Input("cytoscape-clusters-graph", "tapEdgeData"),
    )
    def displayTapEdgeData(data):
        if data:
            return (
                "Source: "
                + data["source"]
                + "\nTarget: "
                + data["target"]
                + "\nPathways: "
                + str(data["weight"])
            )

    return app


if __name__ == "__main__":

    CLUSTERS_PATH = "data/cluster_pop.csv"
    PATHWAYS_FILE = "data/Allpaths_061524.csv"

    nodes, edges = load_nodes(CLUSTERS_PATH), load_edges(PATHWAYS_FILE)

    incytr_cytoscape_app(nodes, edges).run(debug=True)
