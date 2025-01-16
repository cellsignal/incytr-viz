import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import incytr_viz.i_o as i_o
from incytr_viz.callbacks import apply_callbacks
from incytr_viz.components import *
from incytr_viz.components import slider, umap_graph
from dash import Dash, dcc, html
from incytr_viz.modal_content import content
from incytr_viz.util import *


def create_app(raw_pathways, raw_clusters):
    print(ascii())
    # cyto.load_extra_layouts()

    clusters, groups = i_o.load_clusters(raw_clusters)

    pi = i_o.load_pathways(raw_pathways, groups)

    pf = PathwaysFilter(
        all_paths=pi.paths,
        group_a_name=pi.group_a,
        group_b_name=pi.group_b,
    )

    app = Dash(
        __name__,
        suppress_callback_exceptions=True,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
    )

    app.layout = html.Div(
        [
            dcc.Store(id="has-tprs", data=pi.has_tprs),
            dcc.Store(id="has-prs", data=pi.has_prs),
            dcc.Store(id="has-p-value", data=pi.has_p_value),
            dcc.Store(id="has-umap", data=pi.has_umap),
            dcc.Store(id="group-a-name", data=pi.group_a),
            dcc.Store(id="group-b-name", data=pi.group_b),
            html.Div(
                children=[
                    html.Div(
                        [
                            dbc.RadioItems(
                                options=[
                                    {
                                        "label": html.Div(
                                            ["Network View"],
                                        ),
                                        "value": "network",
                                    },
                                    {
                                        "label": html.Div(
                                            ["River View"],
                                        ),
                                        "value": "sankey",
                                    },
                                ],
                                value="network",
                                id="view-radio",
                                style={
                                    "display": "flex",
                                    "justifyContent": "center",
                                    "alignItems": "center",
                                    "border": "1px solid #d3d3d3",
                                },
                                inputClassName="btn-check",
                                labelClassName="btn btn-outline-primary",
                                labelCheckedClassName="active",
                            ),
                            dbc.DropdownMenu(
                                label="Options",
                                children=[
                                    html.Div(
                                        [
                                            dbc.Checkbox(
                                                id="show-network-weights",
                                                label="Show Network Weights",
                                            ),
                                            dbc.Checkbox(
                                                id="show-umap",
                                                label="Show UMAP",
                                                value=False,
                                                disabled=not pi.has_umap,
                                            ),
                                            dcc.Slider(
                                                id="node-scale-factor",
                                                min=1.1,
                                                max=10,
                                                step=0.01,
                                                value=2,
                                                marks=None,
                                                # label="Node Scale Factor",
                                            ),
                                            dcc.Slider(
                                                id="edge-scale-factor",
                                                min=0.1,
                                                max=3,
                                                step=0.1,
                                                value=1,
                                                marks=None,
                                                # label="Edge Scale Factor",
                                            ),
                                            dcc.Slider(
                                                id="label-scale-factor",
                                                min=8,
                                                max=24,
                                                step=1,
                                                value=12,
                                                marks=None,
                                                # label="Edge Scale Factor",
                                            ),
                                        ]
                                    ),
                                    html.Div(
                                        [
                                            dcc.Dropdown(
                                                id="sankey-color-flow-dropdown",
                                                placeholder="Color Sankey Flow By",
                                                multi=False,
                                                clearable=True,
                                                options=[
                                                    "sender",
                                                    "receiver",
                                                    "kinase",
                                                ],
                                                className="filter",
                                            ),
                                        ]
                                    ),
                                ],
                            ),
                        ]
                    ),
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
                                    html.Button(
                                        "Download Current Paths",
                                        id="btn_csv",
                                        className="btn btn-primary",
                                    ),
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
            slider_container(
                has_tprs=pi.has_tprs,
                has_prs=pi.has_prs,
                has_p_value=pi.has_p_value,
            ),
            html.Div(filter_container(pi.paths), className="sidebar"),
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
                            html.Div(
                                umap_graph(
                                    group_id="a",
                                    has_umap=pi.has_umap,
                                    show_umap=pi.has_umap,
                                    all_pathways=pf.a_data,
                                ),
                                className="umapContainer",
                                id="umap-a-container",
                                style={"display": "none"},
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
                            html.Div(
                                umap_graph(
                                    group_id="b",
                                    has_umap=pi.has_umap,
                                    show_umap=pi.has_umap,
                                    all_pathways=pf.b_data,
                                ),
                                className="umapContainer",
                                id="umap-b-container",
                                style={"display": "none"},
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


def get_server(raw_pathways, raw_clusters):
    return create_app(raw_pathways, raw_clusters).server
