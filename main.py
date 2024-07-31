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


def apply_callbacks(app, full_pathways, full_clusters):

    apply_filter_callback(app, full_pathways, full_clusters)
    apply_sankey_callbacks(app)
    apply_cluster_edge_callback(app)

    return app


def format_full_pathways(full_pathways: pd.DataFrame) -> pd.DataFrame:

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


def incytr_app(pathways_file, clusters_file):

    app = Dash(__name__, suppress_callback_exceptions=True)
    logger.info("loading pathways....")
    full_pathways: pd.DataFrame = format_full_pathways(
        pd.read_csv(pathways_file, dtype=pathway_dtypes)
    )

    full_clusters: pd.DataFrame = pd.read_csv(clusters_file)

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

    return apply_callbacks(app, full_pathways, full_clusters)


if __name__ == "__main__":

    CLUSTERS_FILE = "data/mc38/population.csv"
    PATHWAYS_FILE = "data/mc38/mc38_incytr_out_kinase_proteomicsinput.csv"

    app = incytr_app(PATHWAYS_FILE, CLUSTERS_FILE)

    app.run(debug=True)
