import logging
import dash_bootstrap_components as dbc
from dash import Dash, html, dcc
import pdb

from callbacks import apply_callbacks
from components import umap_container

from util import *
from components import *
import i_o
import argparse


logger = logging.getLogger(__name__)


def incytr_app(pathways_path, clusters_a_filepath, clusters_b_filepath):

    (
        paths,
        has_rna_score,
        has_final_score,
        has_p_value,
        has_umap,
        group_a_name,
        group_b_name,
    ) = i_o.load_pathways(pathways_path)

    clusters = i_o.load_cell_clusters(clusters_a_filepath, clusters_b_filepath)

    pf = PathwaysFilter(
        all_paths=paths,
        group_a_name=group_a_name,
        group_b_name=group_b_name,
    )

    app = Dash(
        __name__,
        serve_locally=True,
        suppress_callback_exceptions=True,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
    )

    app.layout = html.Div(
        [
            dcc.Store(id="has-rna", data=has_rna_score),
            dcc.Store(id="has-final", data=has_final_score),
            dcc.Store(id="has-p-value", data=has_p_value),
            dcc.Store(id="has-umap", data=has_umap),
            dcc.Store(id="group-a-name", data=group_a_name),
            dcc.Store(id="group-b-name", data=group_b_name),
            html.Div(
                children=[
                    slider_container(
                        has_rna_score=has_rna_score,
                        has_final_score=has_final_score,
                        has_p_value=has_p_value,
                    ),
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
                            value="sankey",
                            id="view-radio",
                            labelClassName="radioLabel",
                            className="radioContainer sidebarElement",
                        ),
                    ),
                    html.Div(
                        [
                            dbc.Checkbox(id="show-umap", label="Show UMAP"),
                            dcc.Dropdown(
                                id="sankey-color-flow-dropdown",
                                placeholder="Color Sankey Flow By",
                                multi=False,
                                clearable=True,
                                options=["sender", "receiver"],
                                className="filter",
                            ),
                        ],
                        className="sidebarElement",
                    ),
                    filter_container(paths),
                    html.Div(
                        [
                            html.Button("Download Current Paths", id="btn_csv"),
                            dcc.Download(id="download-dataframe-a-csv"),
                            dcc.Download(id="download-dataframe-b-csv"),
                        ]
                    ),
                ],
                className="sidebar",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            umap_container(
                                group_id="a",
                                group_name=group_a_name,
                                has_umap=has_umap,
                                all_pathways=pf.a_data,
                            ),
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
                            umap_container(
                                group_id="b",
                                group_name=group_b_name,
                                has_umap=has_umap,
                                all_pathways=pf.b_data,
                            ),
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
                ],
                className="mainContainer",
                id="main-container",
            ),
        ],
        id="app-container",
        className="app",
    )

    return apply_callbacks(app, paths, clusters)


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

    # html.Div(
    #     [
    #         dbc.Button("Help", id="open", n_clicks=0),
    #         dbc.Modal(
    #             [
    #                 dbc.ModalHeader(dbc.ModalTitle("Header")),
    #                 dbc.ModalBody("This is the content of the modal"),
    #                 dbc.ModalFooter(
    #                     dbc.Button(
    #                         "Close",
    #                         id="close",
    #                         className="ms-auto",
    #                         n_clicks=0,
    #                     )
    #                 ),
    #             ],
    #             id="modal",
    #             size="lg",
    #             is_open=False,
    #         ),
    #     ],
    #     className="modalContainer",
    # ),
