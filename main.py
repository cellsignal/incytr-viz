import pandas as pd
from dash import Dash, html
import logging

from dtypes import pathway_dtypes
from callbacks import (
    apply_filter_callback,
    apply_sankey_callbacks,
    apply_cluster_edge_callback,
)
from components import (
    hist_container,
    pathway_filter_components,
)

from util import *

logger = logging.getLogger(__name__)


def load_clusters_files(clusters_a_path, clusters_b_path):
    def _load_clusters_file(clusters_path):

        clusters_dtypes = {"type": str, "population": int}

        clusters = pd.read_csv(clusters_path, dtype=clusters_dtypes).reset_index(
            drop=True
        )

        clusters.columns = clusters.columns.str.lower().str.strip()

        if not all(c in clusters.columns for c in clusters_dtypes.keys()):
            raise ValueError(
                f"Invalid clusters file: missing one of {clusters_dtypes.keys()}"
            )

        clusters = clusters[list(clusters_dtypes.keys())].reset_index(drop=True)
        return clusters.set_index("type")

    clusters_a = _load_clusters_file(clusters_a_path)
    clusters_b = _load_clusters_file(clusters_b_path)

    merged_clusters = clusters_a.join(
        clusters_b, on="type", how="outer", lsuffix="_a", rsuffix="_b"
    )

    merged_clusters["population_a"] = (
        merged_clusters["population_a"].fillna(0).astype(int)
    )
    merged_clusters["population_b"] = (
        merged_clusters["population_b"].fillna(0).astype(int)
    )

    return merged_clusters


def apply_callbacks(app, full_pathways, full_clusters):

    apply_filter_callback(app, full_pathways, full_clusters)
    apply_sankey_callbacks(app)
    apply_cluster_edge_callback(app)

    return app


def format_full_pathways(full_pathways: pd.DataFrame) -> pd.DataFrame:

    full_pathways.columns = full_pathways.columns.str.strip()
    full_pathways["Ligand"] = full_pathways["Path"].str.split("*").str[0]
    full_pathways["Receptor"] = full_pathways["Path"].str.split("*").str[1]
    full_pathways["EM"] = full_pathways["Path"].str.split("*").str[2]
    full_pathways["Target"] = full_pathways["Path"].str.split("*").str[3]

    full_pathways["SigWeight"] = full_pathways.apply(
        lambda row: (
            row["SigWeight_X"] if row["final_score"] > 0 else row["SigWeight_Y"]
        ),
        axis=1,
    )

    TO_KEEP = [
        "Path",
        "Ligand",
        "Receptor",
        "EM",
        "Target",
        "SigWeight",
        "final_score",
        "SigWeight_X",
        "SigWeight_Y",
        "RNA_score",
        "Sender",
        "Receiver",
        "adjlog2FC",
    ]

    return full_pathways[TO_KEEP]


def incytr_app(pathways_file, clusters_a_filepath, clusters_b_filepath):

    app = Dash(__name__, suppress_callback_exceptions=True)
    logger.info("loading pathways....")
    full_pathways: pd.DataFrame = format_full_pathways(
        pd.read_csv(pathways_file, dtype=pathway_dtypes)
    )

    clusters = load_clusters_files(clusters_a_filepath, clusters_b_filepath)

    app.layout = html.Div(
        [
            pathway_filter_components(full_pathways),
            html.Div(
                id="data-container",
                children=[
                    hist_container({"display": "flex", "flexDirection": "row"}),
                    html.Div(id="figures-container"),
                ],
                style={"display": "flex", "flexDirection": "column"},
            ),
        ],
        id="app-container",
        style={"display": "flex", "width": "100vw"},
    )

    return apply_callbacks(app, full_pathways, clusters)


if __name__ == "__main__":

    CLUSTERS_A_FILE = "data/mc38_fake/population_a.csv"
    CLUSTERS_B_FILE = "data/mc38_fake/population_b.csv"
    PATHWAYS_FILE = "data/mc38_fake/pathways_small_short.csv"

    app = incytr_app(PATHWAYS_FILE, CLUSTERS_A_FILE, CLUSTERS_B_FILE)

    app.run(debug=True)
