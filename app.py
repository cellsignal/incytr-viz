import logging
import pandas as pd
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash import Dash, html

from util import *
from dtypes import pathway_dtypes
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


def incytr_app(pathways_file, clusters_a_filepath, clusters_b_filepath):

    logger.info("loading pathways....")

    if ".csv" in pathways_file:
        sep = ","
    elif ".tsv" in pathways_file:
        sep = "\t"
    else:
        raise ValueError("Pathways file must be a CSV or TSV -- check filename")

    full_pathways: pd.DataFrame = load_pathways_input(
        pd.read_csv(pathways_file, dtype=pathway_dtypes, sep=sep)
    )

    has_rna_score = CN.rna_score_available(full_pathways)
    has_final_score = CN.final_score_available(full_pathways)
    has_p_value = CN.p_value_available(full_pathways)

    TO_KEEP = [
        get_cn("path"),
        get_cn("ligand"),
        get_cn("receptor"),
        get_cn("em"),
        get_cn("target"),
        get_cn("sender"),
        get_cn("receiver"),
        CN.SIGWEIGHT(cols=full_pathways.columns, group="a"),
        CN.SIGWEIGHT(cols=full_pathways.columns, group="b"),
    ]

    if has_final_score:
        TO_KEEP.append(get_cn("final_score"))
    if has_rna_score:
        TO_KEEP.append(get_cn("rna_score"))
    if has_p_value:
        TO_KEEP += [
            CN.PVAL(full_pathways.columns, "a"),
            CN.PVAL(full_pathways.columns, "b"),
        ]

    full_pathways = full_pathways[TO_KEEP]

    clusters = load_cell_populations(clusters_a_filepath, clusters_b_filepath)

    app = Dash(
        __name__,
        serve_locally=True,
        suppress_callback_exceptions=True,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
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
                            dbc.Button("Open modal", id="open", n_clicks=0),
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
        app, full_pathways, clusters, has_rna_score, has_final_score, has_p_value
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
