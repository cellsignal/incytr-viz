import logging
import dash_bootstrap_components as dbc
from dash import Dash, html, dcc

from util import *
import i_o
import argparse


logger = logging.getLogger(__name__)


def incytr_app(pathways_path, clusters_a_filepath, clusters_b_filepath):

    paths, has_rna_score, has_final_score, has_p_value = i_o.load_pathways(
        pathways_path
    )

    clusters_a = i_o.load_cell_clusters(clusters_a_filepath)
    clusters_b = i_o.load_cell_clusters(clusters_b_filepath)

    all_clusters = clusters_a.merge(
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
            html.Div(
                [
                    dcc.Store(id="has-rna", data=has_rna_score),
                    dcc.Store(id="has-final", data=has_final_score),
                    dcc.Store(id="has-p-value", data=has_p_value),
                    html.Div(
                        children=[
                            html.Div(
                                [
                                    html.Div([], id="filter-container"),
                                    html.Div(
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
                                    ),
                                ],
                            ),
                            html.Div(
                                [],
                                id="slider-container",
                                className="sliderContainer sidebarElement",
                            ),
                        ],
                    ),
                ],
                className="sidebar",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div([], id="hist-a-container"),
                                    html.Div([], id="figure-a-container"),
                                ],
                                id="group-a-container",
                                className="groupContainer",
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div([], id="hist-b-container"),
                                    html.Div([], id="figure-b-container"),
                                ],
                                id="group-b-container",
                                className="groupContainer",
                            ),
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

    return app
    # return apply_callbacks(
    #     app, paths, all_clusters, has_rna_score, has_final_score, has_p_value
    # )


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
