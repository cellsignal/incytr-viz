import logging
import pandas as pd
from dash import Dash, html

from util import *
from dtypes import pathway_dtypes
import argparse

from callbacks import (
    apply_filter_callback,
    apply_sankey_callbacks,
    apply_cluster_edge_callback,
)
from components import (
    sidebar,
)

logger = logging.getLogger(__name__)

from data import load_pathways_input, load_cell_populations


def apply_callbacks(
    app, full_pathways, full_clusters, has_rna_score, has_final_score, has_p_value
):

    apply_filter_callback(
        app,
        full_pathways,
        full_clusters,
        has_rna_score=has_rna_score,
        has_final_score=has_final_score,
        has_p_value=has_p_value,
    )
    apply_sankey_callbacks(app)
    apply_cluster_edge_callback(app)

    return app


def incytr_app(pathways_file, clusters_a_filepath, clusters_b_filepath):

    logger.info("loading pathways....")
    full_pathways: pd.DataFrame = load_pathways_input(
        pd.read_csv(pathways_file, dtype=pathway_dtypes)
    )

    has_rna_score = CN.rna_score_available(full_pathways)
    has_final_score = CN.final_score_available(full_pathways)
    has_p_value = CN.p_value_available(full_pathways)

    clusters = load_cell_populations(clusters_a_filepath, clusters_b_filepath)

    app = Dash(
        __name__,
        serve_locally=True,
        suppress_callback_exceptions=True,
    )
    app.layout = html.Div(
        [
            sidebar(
                full_pathways,
                has_final_score=has_final_score,
                has_p_value=has_p_value,
                has_rna_score=has_rna_score,
            ),
            html.Div(
                [
                    html.Div([], id="group-a-container", className="groupContainer"),
                    html.Div([], id="group-b-container", className="groupContainer"),
                ],
                style={"display": "flex", "flexDirection": "row"},
                id="main-container",
                className="mainContainer",
            ),
        ],
        id="app-container",
        className="app",
    )

    return apply_callbacks(
        app, full_pathways, clusters, has_rna_score, has_final_score, has_p_value
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run the InCytr visualization app.")
    parser.add_argument(
        "--group_a_popluations",
        type=str,
        required=True,
        help="Path to clusters A CSV file",
    )
    parser.add_argument(
        "--group_b_populations",
        type=str,
        required=True,
        help="Path to clusters B CSV file",
    )
    parser.add_argument(
        "--pathways", type=str, required=True, help="Path to pathways CSV file"
    )

    args = parser.parse_args()

    CLUSTERS_A_FILE = args.group_a_popluations
    CLUSTERS_B_FILE = args.group_b_populations
    PATHWAYS_FILE = args.pathways

    app = incytr_app(PATHWAYS_FILE, CLUSTERS_A_FILE, CLUSTERS_B_FILE)

    app.run(debug=True)
