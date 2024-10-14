import logging
import pandas as pd
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash import Dash, html

from util import *
import i_o
import argparse

from callbacks import (
    apply_filter_callback,
    apply_modal_callbacks,
    apply_sankey_callbacks,
    apply_cluster_edge_callback,
    apply_node_tap_callback,
)
from components import (
    sidebar,
)

logger = logging.getLogger(__name__)


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
    apply_modal_callbacks(app)
    apply_cluster_edge_callback(app)
    apply_node_tap_callback(
        app,
        outputs=[Output("metrics-a", "children")],
        inputs=[Input("cytoscape-a", "tapNode")],
    )
    apply_node_tap_callback(
        app,
        outputs=[Output("metrics-b", "children")],
        inputs=[Input("cytoscape-b", "tapNode")],
    )

    return app


def incytr_app(pathways_path, clusters_a_filepath, clusters_b_filepath):

    paths, has_rna_score, has_final_score, has_p_value = i_o.load_pathways(
        pathways_path
    )

    clusters_a = i_o.load_cell_clusters(clusters_a_filepath)
    clusters_b = i_o.load_cell_clusters(clusters_b_filepath)

    clusters_merged = clusters_a.merge(
        clusters_b, on="type", how="outer", suffixes=("_a", "_b")
    )

    app = Dash(
        __name__,
        serve_locally=True,
        suppress_callback_exceptions=True,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
    )
    app.layout = html.Div(
        [
            sidebar(
                paths,
                has_final_score=has_final_score,
                has_p_value=has_p_value,
                has_rna_score=has_rna_score,
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [], id="group-a-container", className="groupContainer"
                            ),
                            html.Div([], id="metrics-a", className="metricsContainer"),
                        ]
                    ),
                    html.Div(
                        [
                            html.Div(
                                [], id="group-b-container", className="groupContainer"
                            ),
                            html.Div([], id="metrics-b", className="metricsContainer"),
                        ]
                    ),
                    html.Div(
                        [
                            dbc.Button("Help", id="open", n_clicks=0),
                            dbc.Modal(
                                [
                                    dbc.ModalHeader(dbc.ModalTitle("Header")),
                                    dbc.ModalBody("This is the content of the modal"),
                                    dbc.ModalFooter(
                                        dbc.Button(
                                            "Close",
                                            id="close",
                                            className="ms-auto",
                                            n_clicks=0,
                                        )
                                    ),
                                ],
                                id="modal",
                                size="lg",
                                is_open=False,
                            ),
                        ],
                        className="modalContainer",
                    ),
                ],
                id="main-container",
                className="mainContainer",
            ),
        ],
        id="app-container",
        className="app",
    )

    return apply_callbacks(
        app, paths, clusters_merged, has_rna_score, has_final_score, has_p_value
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run the InCytr visualization app.")
    parser.add_argument(
        "--group_a_populations",
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

    CLUSTERS_A_FILE = args.group_a_populations
    CLUSTERS_B_FILE = args.group_b_populations
    PATHWAYS_FILE = args.pathways

    app = incytr_app(PATHWAYS_FILE, CLUSTERS_A_FILE, CLUSTERS_B_FILE)

    app.run(debug=True)
