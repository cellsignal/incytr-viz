import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import incytr_viz.i_o as i_o
from incytr_viz.callbacks import apply_callbacks
from incytr_viz.components import *
from incytr_viz.components import slider, umap_container
from dash import Dash, dcc, html
from incytr_viz.modal_content import content
from incytr_viz.util import *

# import dash_bootstrap_components as dbc
# from dash import Input, Output, State, html


# @app.callback(
#     Output("collapse", "is_open"),
#     [Input("collapse-button", "n_clicks")],
#     [State("collapse", "is_open")],
# )
# def toggle_collapse(n, is_open):
#     if n:
#         return not is_open
#     return is_open


def create_app(pathways, clusters_a, clusters_b):

    clusters = i_o.load_cell_clusters(clusters_a, clusters_b)

    pi = i_o.process_input_data(pathways)

    pf = PathwaysFilter(
        all_paths=pi.paths,
        group_a_name=pi.group_a,
        group_b_name=pi.group_b,
    )

    cyto.load_extra_layouts()

    app = Dash(
        __name__,
        suppress_callback_exceptions=True,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
    )

    app.layout = html.Div(
        [
            dcc.Store(id="has-rna", data=pi.has_rna),
            dcc.Store(id="has-final", data=pi.has_final),
            dcc.Store(id="has-p-value", data=pi.has_p_value),
            dcc.Store(id="has-umap", data=pi.has_umap),
            dcc.Store(id="group-a-name", data=pi.group_a),
            dcc.Store(id="group-b-name", data=pi.group_b),
            html.Div(
                children=[
                    slider_container(
                        has_rna_score=pi.has_rna,
                        has_final_score=pi.has_final,
                        has_p_value=pi.has_p_value,
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
                            value="network",
                            id="view-radio",
                            labelClassName="radioLabel",
                            className="radioContainer sidebarElement",
                        ),
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H6("Network Options"),
                                    dbc.Checkbox(
                                        id="show-network-weights",
                                        label="Show Network Weights",
                                    ),
                                    html.Div(
                                        [
                                            dbc.Button(
                                                "Open collapse",
                                                id="collapse-button",
                                                className="mb-3",
                                                color="primary",
                                                n_clicks=0,
                                            ),
                                            dbc.Collapse(
                                                dbc.Card(
                                                    dbc.CardBody(
                                                        [
                                                            dcc.Slider(
                                                                id="node-scale-factor",
                                                                min=1.1,
                                                                max=10,
                                                                step=0.01,
                                                                value=2,
                                                                marks=None,
                                                                # label="Node Scale Factor",
                                                                vertical=False,
                                                            ),
                                                            dcc.Slider(
                                                                id="edge-scale-factor",
                                                                min=0.1,
                                                                max=3,
                                                                step=0.1,
                                                                value=1,
                                                                marks=None,
                                                                # label="Edge Scale Factor",
                                                                vertical=False,
                                                            ),
                                                        ],
                                                        style={
                                                            "height": "175px",
                                                            "width": "200px",
                                                            "border": "1px solid #000",
                                                        },
                                                    ),
                                                ),
                                                id="collapse",
                                                is_open=True,
                                            ),
                                        ]
                                    ),
                                ]
                            ),
                            html.Div(
                                [
                                    html.H6("Sankey Options"),
                                    dcc.Dropdown(
                                        id="sankey-color-flow-dropdown",
                                        placeholder="Color Sankey Flow By",
                                        multi=False,
                                        clearable=True,
                                        options=["sender", "receiver", "kinase"],
                                        className="filter",
                                    ),
                                ]
                            ),
                        ],
                        className="sidebarElement figureSpecificOptions",
                    ),
                    filter_container(pi.paths),
                    html.Div(
                        [
                            dbc.Button("Help", id="open", n_clicks=0),
                            dbc.Modal(
                                [
                                    dbc.ModalHeader(
                                        dbc.ModalTitle("Incytr Data Visualization")
                                    ),
                                    dbc.ModalBody(content),
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
                                size="xl",
                                is_open=False,
                            ),
                            html.Div(
                                [
                                    html.Button("Download Current Paths", id="btn_csv"),
                                    dcc.Download(id="download-dataframe-a-csv"),
                                    dcc.Download(id="download-dataframe-b-csv"),
                                ]
                            ),
                        ],
                        className="sidebarElement",
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
                                    html.H3(
                                        pi.group_a,
                                        style={"textTransform": "uppercase"},
                                    ),
                                    html.Div(
                                        [
                                            html.Span("Pathways Displayed: "),
                                            html.Span(0, id="pathways-count-a"),
                                        ],
                                        style={
                                            "width": "20%",
                                            "display": "flex",
                                            "justify-content": "space-between",
                                        },
                                    ),
                                ],
                                className="groupTitle",
                            ),
                            umap_container(
                                group_id="a",
                                group_name=pi.group_a,
                                has_umap=pi.has_umap,
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
                        ],
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H3(
                                        pi.group_b,
                                        style={"textTransform": "uppercase"},
                                    ),
                                    html.Div(
                                        [
                                            html.Span("Pathways Displayed: "),
                                            html.Span(0, id="pathways-count-b"),
                                        ],
                                        style={
                                            "width": "20%",
                                            "display": "flex",
                                            "justify-content": "space-between",
                                        },
                                    ),
                                ],
                                className="groupTitle",
                            ),
                            umap_container(
                                group_id="b",
                                group_name=pi.group_b,
                                has_umap=pi.has_umap,
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

    return apply_callbacks(app, pi.paths, clusters)


def get_server(pathways, clusters_a, clusters_b):
    return create_app(pathways, clusters_a, clusters_b).server
