import json
from typing import Optional

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
from dash import ALL, Dash, callback, ctx, dcc, html
from dash.dependencies import Input, Output, State
from flask import current_app

from incytr_viz.components import (
    create_hist_figure,
    cytoscape_container,
    filter_container,
    sankey_container,
    slider_container,
    umap_graph,
)
from incytr_viz.util import *

logger = create_logger(__name__)


def create_app(pathways_file, clusters_file):
    app = Dash(
        __name__,
        suppress_callback_exceptions=True,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
    )

    incytr_input = IncytrInput(clusters_path=clusters_file, pathways_path=pathways_file)

    app.server.config["INCYTR_INPUT"] = incytr_input

    app.layout = html.Div(
        [
            dbc.NavbarSimple(
                children=[
                    dbc.NavItem(
                        html.Div(
                            dbc.RadioItems(
                                options=[
                                    {
                                        "label": "Network View",
                                        "value": "network",
                                    },
                                    {
                                        "label": "River View",
                                        "value": "sankey",
                                    },
                                ],
                                value="network",
                                id="view-radio",
                                className="btn-group",
                                inputClassName="btn btn-check",
                                labelClassName="btn btn-outline-primary",
                                labelCheckedClassName="active",
                            ),
                        )
                    ),
                    dbc.NavItem(
                        children=[
                            dbc.Button(
                                "Reset Filters",
                                id="reset",
                                className="btn btn-primary",
                            ),
                        ],
                    ),
                    dbc.NavItem(
                        dbc.DropdownMenu(
                            label="Options",
                            children=html.Div(
                                [
                                    dbc.Checkbox(
                                        id="show-network-weights",
                                        label="Show Network Weights",
                                    ),
                                    dbc.Checkbox(
                                        id="show-population-fractions",
                                        label="Show Cluster Population Sizes",
                                    ),
                                    dbc.Checkbox(
                                        id="restrict-afc",
                                        label="Restrict on aFC direction (recommended)",
                                        value=True,
                                    ),
                                    dbc.Checkbox(
                                        id="show-umap",
                                        label="Show UMAP",
                                        value=False,
                                        disabled=not incytr_input.has_umap,
                                    ),
                                    html.Div(
                                        [
                                            dcc.Slider(
                                                id="node-scale-factor",
                                                min=1.1,
                                                max=10,
                                                step=0.01,
                                                value=2,
                                                marks=None,
                                                className="scaleFactor",
                                            ),
                                            html.Div("Scale Network Nodes"),
                                        ],
                                        className="optionSlider",
                                    ),
                                    html.Div(
                                        [
                                            dcc.Slider(
                                                id="edge-scale-factor",
                                                min=0.1,
                                                max=3,
                                                step=0.1,
                                                value=1,
                                                marks=None,
                                                className="scaleFactor",
                                            ),
                                            html.Div("Scale Network Edges"),
                                        ],
                                        className="optionSlider",
                                    ),
                                    # html.Div(
                                    #     [
                                    #         dcc.Slider(
                                    #             id="label-scale-factor",
                                    #             min=8,
                                    #             max=24,
                                    #             step=1,
                                    #             value=12,
                                    #             marks=None,
                                    #             className="scaleFactor",
                                    #         ),
                                    #         html.Div("Scale Label Size"),
                                    #     ],
                                    #     className="optionSlider",
                                    # ),
                                    dcc.Dropdown(
                                        id="sankey-color-flow-dropdown",
                                        placeholder="Color River Flow",
                                        multi=False,
                                        clearable=True,
                                        options=(
                                            [
                                                "sender",
                                                "receiver",
                                                "kinase",
                                            ]
                                            if incytr_input.has_kinase
                                            else ["sender", "receiver"]
                                        ),
                                    ),
                                ],
                                style={"padding": "5px 5px", "width": "300px"},
                            ),
                        )
                    ),
                    dbc.NavItem(dbc.Button("Help", id="open", n_clicks=0)),
                    dbc.Modal(
                        [
                            dbc.ModalHeader(
                                dbc.ModalTitle("Incytr Data Visualization")
                            ),
                            dbc.ModalBody(dcc.Markdown(children=get_help_file())),
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
                    dbc.NavItem(
                        html.Button(
                            "Download Current Paths",
                            id="btn_csv",
                            className="btn btn-primary",
                        ),
                    ),
                    dcc.Download(id="download-dataframe-a-csv"),
                    dcc.Download(id="download-dataframe-b-csv"),
                ],
                brand="Incytr Pathway Visualization",
                brand_href="#",
                color="primary",
                dark=True,
            ),
            html.Div(
                slider_container(
                    has_tpds=incytr_input.has_tpds,
                    has_ppds=incytr_input.has_ppds,
                    has_p_value=incytr_input.has_p_value,
                ),
                id="slider-container",
            ),
            html.Div(
                filter_container(
                    sender=list(incytr_input.unique_senders),
                    receiver=list(incytr_input.unique_receivers),
                    ligand=list(incytr_input.unique_ligands),
                    receptor=list(incytr_input.unique_receptors),
                    em=list(incytr_input.unique_em),
                    target=list(incytr_input.unique_targets),
                ),
                className="sidebar",
                id="filter-container",
            ),
            dcc.Loading(
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Div(
                                            [
                                                html.Span(
                                                    "⊕",
                                                ),
                                                html.Span(incytr_input.group_a.title()),
                                            ],
                                            className="groupTitle",
                                        ),
                                        html.Div(
                                            [
                                                html.Span(
                                                    "Pathways Displayed:",
                                                    style={"paddingRight": "5px"},
                                                ),
                                                html.Span(0, id="pathways-count-a"),
                                            ],
                                            className="pathwaysCount",
                                        ),
                                    ],
                                    className="groupHeader",
                                ),
                                html.Div(
                                    umap_graph(
                                        "a", incytr_input.has_umap, incytr_input.paths
                                    ),
                                    className="umapContainer",
                                    id="umap-a-container",
                                    style={"display": "none"},
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            [dcc.Graph()], id="figure-a-container"
                                        ),
                                        html.Div(
                                            [dcc.Graph(id="hist-a-graph")],
                                            id="hist-a-container",
                                            className="histContainer",
                                        ),
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
                                        html.Div(
                                            [
                                                html.Span(
                                                    "⊖",
                                                ),
                                                html.Span(incytr_input.group_b.title()),
                                            ],
                                            className="groupTitle",
                                        ),
                                        html.Div(
                                            [
                                                html.Span(
                                                    "Pathways Displayed:",
                                                    style={"paddingRight": "5px"},
                                                ),
                                                html.Span(0, id="pathways-count-b"),
                                            ],
                                            className="pathwaysCount",
                                        ),
                                    ],
                                    className="groupHeader",
                                ),
                                html.Div(
                                    umap_graph(
                                        "b", incytr_input.has_umap, incytr_input.paths
                                    ),
                                    className="umapContainer",
                                    id="umap-b-container",
                                    style={"display": "none"},
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            [dcc.Graph()], id="figure-b-container"
                                        ),
                                        html.Div(
                                            [dcc.Graph(id="hist-b-graph")],
                                            id="hist-b-container",
                                            className="histContainer",
                                        ),
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
                id="loading",
                delay_show=500,
            ),
        ],
        id="app-container",
        className="app",
    )

    return app.server


def load_nodes(clusters: pd.DataFrame, node_scale_factor) -> list[dict]:
    """
    Generate cytoscape nodes from clusters file

    clusters: clusters df with expected column names:

    type
    population_a
    population_b
    rgb_colors

    Output:

    [nodes_a, nodes_b]

    Each member of the list is a list of dicts, each dict is a node

    nodes_a ~ [{"data": {...node_data}}, .....]

    """
    # TODO clean clusters
    # clusters = clean_clusters(clusters)

    def calculate_node_diameters(clusters, node_scale_factor):

        clusters = clusters.copy()
        if (clusters["population"] <= 0).all():

            clusters.loc[:, "node_diameter"] = 0
            return clusters

        min_pop = clusters[clusters["population"] > 0]["population"].min()
        clusters.loc[:, "pop_min_ratio"] = clusters.loc[:, "population"] / min_pop
        clusters.loc[:, "normalized_ratio"] = (
            clusters.loc[:, "pop_min_ratio"] * node_scale_factor
        )
        clusters.loc[:, "node_area"] = np.round(
            400 * (log_base(clusters["normalized_ratio"], node_scale_factor)),
            4,
        )
        clusters.loc[:, "node_diameter"] = np.round(
            np.sqrt(4 * clusters["node_area"] / np.pi), 4
        )
        return clusters

    clusters = calculate_node_diameters(clusters, node_scale_factor)

    def _add_node(row: pd.Series) -> dict:

        node_type = row.name
        node_label = row.type_userlabel
        node_population = row["population"]

        if (not node_population) or (np.isnan(node_population)):
            return np.nan

        data = dict()
        data["id"] = node_type
        data["label"] = node_label
        data["cluster_size"] = node_population
        data["label_with_size"] = f"{node_label} ({node_population:.4f})"
        data["width"] = row["node_diameter"]
        data["height"] = row["node_diameter"]
        data["background_color"] = row["color"]
        return {"data": data}

    return list(
        clusters.apply(
            lambda row: _add_node(row),
            axis=1,
        ).dropna()
    )


def load_edges(
    nodes: list[dict],
    pathways: pd.DataFrame,
    global_max_paths: int,
    edge_scale_factor: float,
):
    """add pathways from source to target"""
    edges = []

    ## filter pathways if sender/receiver not in nodes
    node_ids = pd.Series([x["data"]["id"] for x in nodes])
    pathways = pathways[
        (pathways["sender"].isin(node_ids)) & (pathways["receiver"].isin(node_ids))
    ]

    ## filter pathways that are below sigprob threshold

    if len(pathways) == 0:
        return edges

    s: pd.Series = pathways.groupby(["sender", "receiver"]).size()

    sr_pairs = s.to_dict()
    for sr, weight in sr_pairs.items():
        source_id, target_id = sr
        data = dict()
        data["id"] = source_id + target_id
        data["source"] = source_id
        data["target"] = target_id
        data["weight"] = weight
        data["label"] = str(weight)
        data["line_color"] = next(
            x["data"]["background_color"] for x in nodes if x["data"]["id"] == source_id
        )

        edges.append({"data": data})

    if edges:
        for e in edges:
            e["data"]["width"] = edge_width_map(
                abs(e["data"]["weight"]),
                edge_scale_factor=edge_scale_factor,
                global_max_paths=global_max_paths,
            )

    return edges


def pathways_df_to_sankey(
    sankey_df: pd.DataFrame,
    all_clusters: pd.DataFrame,
    sankey_color_flow: Optional[str] = None,  # sender or receiver
) -> tuple:

    def _get_values(
        df: pd.DataFrame, source_colname: str, target_colname: str
    ) -> pd.DataFrame:

        def _kinase_color_map(row, source, target):

            mapper = kinase_color_map()
            if (source == "receptor") and (target == "em"):
                if row["sik_r_of_em"] and row["sik_em_of_r"]:
                    return mapper["Bidirectional"]
                elif row["sik_r_of_em"]:
                    return mapper["Receptor --> Effector"]
                elif row["sik_em_of_r"]:
                    return mapper["Effector --> (Receptor/Target Gene)"]
            elif (source == "em") and (target == "target"):
                if row["sik_em_of_t"] and row["sik_t_of_em"]:
                    return mapper["Bidirectional"]
                elif row["sik_em_of_t"]:
                    return mapper["Effector --> (Receptor/Target Gene)"]
                elif row["sik_t_of_em"]:
                    return mapper["Target --> Effector"]
            return "lightgrey"

        if sankey_color_flow in ["sender", "receiver"]:
            color_grouping_column = sankey_color_flow
            out = (
                df.groupby([source_colname, color_grouping_column])[target_colname]
                .value_counts()
                .reset_index(name="value")
            )
            out[color_grouping_column] = (
                out[color_grouping_column].astype(str).str.lower()
            )
            out["color"] = out[color_grouping_column].map(
                dict(zip(all_clusters.index, all_clusters["color"]))
            )
        elif sankey_color_flow == "kinase":

            # only coloring adjacent
            adjacent_kinase_cols = [
                "sik_r_of_em",
                "sik_em_of_r",
                "sik_em_of_t",
                "sik_t_of_em",
            ]

            kinase_mask = df[adjacent_kinase_cols].any(axis=1)
            df.loc[kinase_mask, "kinase_color_map"] = df.loc[kinase_mask].apply(
                lambda row: _kinase_color_map(row, source_colname, target_colname),
                axis=1,
            )
            df.loc[~kinase_mask, "kinase_color_map"] = "lightgrey"
            out = (
                df.groupby([source_colname, "kinase_color_map"])[target_colname]
                .value_counts()
                .reset_index(name="value")
            )
            out["color"] = out["kinase_color_map"]

        else:
            out = (
                df.groupby([source_colname])[target_colname]
                .value_counts()
                .reset_index(name="value")
            )
            out["color"] = "lightgrey"

        out.rename(
            columns={
                source_colname: "source",
                target_colname: "target",
            },
            inplace=True,
        )
        out["source_id"] = out["source"] + "_" + source_colname
        out["target_id"] = out["target"] + "_" + target_colname

        return out

    l_r = _get_values(sankey_df, "ligand", "receptor")
    r_em = _get_values(sankey_df, "receptor", "em")
    em_t = _get_values(sankey_df, "em", "target")

    included_links = [l_r, r_em]

    ## auto-determine if target genes should be included
    def _should_display_targets() -> bool:
        num_targets = len(em_t["target"].unique())

        return num_targets <= 400
        # return True

    if _should_display_targets():
        included_links.append(em_t)

    links = pd.concat(included_links, axis=0).reset_index(drop=True)
    # ids allow for repeating labels in ligand, receptor, etc. without pointing to same node
    ids = list(set(pd.concat([links["source_id"], links["target_id"]])))
    labels = [x.split("_")[0] for x in ids]

    source = [next(i for i, e in enumerate(ids) if e == x) for x in links["source_id"]]
    target = [next(i for i, e in enumerate(ids) if e == x) for x in links["target_id"]]
    value = links["value"]

    color = links["color"]
    return (ids, labels, source, target, value, color)


def pathway_component_filter_inputs(state=False):
    klass = State if state else Input
    return dict(
        sender_select=klass("sender-select", "value"),
        receiver_select=klass("receiver-select", "value"),
        ligand_select=klass("ligand-select", "value"),
        receptor_select=klass("receptor-select", "value"),
        em_select=klass("em-select", "value"),
        target_select=klass("target-select", "value"),
        any_role_select=klass("any-role-select", "value"),
        sankey_color_flow=klass("sankey-color-flow-dropdown", "value"),
        umap_select_a=klass("umap-select-a", "value"),
        umap_select_b=klass("umap-select-b", "value"),
        kinase_select=klass("kinase-select", "value"),
        restrict_afc=klass("restrict-afc", "value"),
    )


def network_style_inputs(state=False):
    klass = State if state else Input
    return dict(
        node_scale_factor=klass("node-scale-factor", "value"),
        edge_scale_factor=klass("edge-scale-factor", "value"),
        # label_scale_factor=klass("label-scale-factor", "value"),
    )


@callback(
    output=dict(
        hist_a=Output("hist-a-graph", "figure"),
        hist_b=Output("hist-b-graph", "figure"),
        figure_a=Output("figure-a-container", "children"),
        figure_b=Output("figure-b-container", "children"),
        num_paths_a=Output("pathways-count-a", "children"),
        num_paths_b=Output("pathways-count-b", "children"),
    ),
    inputs=dict(
        pcf=pathway_component_filter_inputs(),
        nsi=network_style_inputs(),
        slider_changed=Input({"type": "numerical-filter", "index": ALL}, "value"),
        sliders_container_children=State("allSlidersContainer", "children"),
        view_radio=Input("view-radio", "value"),
    ),
    state=dict(
        show_network_weights=State("show-network-weights", "value"),
    ),
    # prevent_initial_call=True,
)
def update_figure_and_histogram(
    pcf,
    nsi,
    slider_changed,
    sliders_container_children,
    view_radio,
    show_network_weights,
):

    incytr_input = current_app.config["INCYTR_INPUT"]
    clusters = incytr_input.clusters

    filter_umap_a = parse_umap_filter_data(pcf.get("umap_select_a"))
    filter_umap_b = parse_umap_filter_data(pcf.get("umap_select_b"))

    slider_values = parse_slider_values_from_tree(sliders_container_children)

    pf = PathwaysFilter(
        all_paths=incytr_input.paths,
        group_a_name=incytr_input.group_a,
        group_b_name=incytr_input.group_b,
        filter_afc_direction=pcf.get("restrict_afc"),
        filter_umap_a=filter_umap_a,
        filter_umap_b=filter_umap_b,
        filter_senders=pcf.get("sender_select"),
        filter_receivers=pcf.get("receiver_select"),
        filter_ligands=pcf.get("ligand_select"),
        filter_receptors=pcf.get("receptor_select"),
        filter_kinase=pcf.get("kinase_select"),
        filter_em=pcf.get("em_select"),
        filter_target_genes=pcf.get("target_select"),
        filter_all_molecules=pcf.get("any_role_select"),
        ppds_bounds=incytr_input.has_ppds and slider_values.get("ppds"),
        sp_threshold=slider_values.get("sigprob"),
        tppds_bounds=incytr_input.has_tpds and slider_values.get("tpds"),
        pval_threshold=incytr_input.has_p_value and slider_values.get("p-value"),
    )

    a_pathways = pf.filter("a", should_filter_umap=bool(filter_umap_a))
    b_pathways = pf.filter("b", should_filter_umap=bool(filter_umap_b))

    def _get_group_figures(
        filtered_group_paths: pd.DataFrame,
        clusters: pd.DataFrame,
        group_name: str,
        group_id: str,
        global_max_paths: int,
    ):

        if view_radio == "network":
            nodes = load_nodes(
                clusters.loc[clusters["group"] == group_name],
                node_scale_factor=nsi.get("node_scale_factor", 2),
            )
            edges = load_edges(
                nodes,
                filtered_group_paths,
                global_max_paths,
                edge_scale_factor=nsi.get(
                    "edge_scale_factor",
                ),
            )

            cytoscape = cytoscape_container(
                f"cytoscape-{group_id}",
                group_name,
                nodes + edges,
                show_network_weights=show_network_weights,
            )

            graph_container = cytoscape

        elif view_radio == "sankey":

            ids, labels, source, target, value, color = pathways_df_to_sankey(
                sankey_df=filtered_group_paths,
                sankey_color_flow=pcf.get("sankey_color_flow"),
                all_clusters=clusters,
            )

            sankey = sankey_container(
                clusters,
                ids,
                labels,
                source,
                target,
                value,
                color,
                group_id,
                color_flow=pcf.get("sankey_color_flow"),
            )

            graph_container = sankey
        return [
            graph_container,
            create_hist_figure(
                paths=filtered_group_paths,
                has_tpds=incytr_input.has_tpds,
                has_ppds=incytr_input.has_ppds,
                has_p_value=incytr_input.has_p_value,
            ),
        ]

    a_max_paths = np.max(a_pathways.groupby(["sender", "receiver"]).size())
    b_max_paths = np.max(b_pathways.groupby(["sender", "receiver"]).size())

    if np.isnan(a_max_paths):
        a_max_paths = 0
    if np.isnan(b_max_paths):
        b_max_paths = 0

    global_max_paths = max(a_max_paths, b_max_paths)

    group_a_figs = _get_group_figures(
        filtered_group_paths=a_pathways,
        clusters=clusters,
        global_max_paths=global_max_paths,
        group_name=incytr_input.group_a,
        group_id="a",
    )
    group_b_figs = _get_group_figures(
        filtered_group_paths=b_pathways,
        clusters=clusters,
        global_max_paths=global_max_paths,
        group_name=incytr_input.group_b,
        group_id="b",
    )
    num_paths_a = len(a_pathways)
    num_paths_b = len(b_pathways)

    return dict(
        hist_a=group_a_figs[1],
        hist_b=group_b_figs[1],
        figure_a=group_a_figs[0],
        figure_b=group_b_figs[0],
        num_paths_a=num_paths_a,
        num_paths_b=num_paths_b,
    )


def _relayout_umap(relayoutData):
    """
    {
        'xaxis.range[0]': -0.6369249007630504,
        'xaxis.range[1]': 6.965720316453904,
        'yaxis.range[0]': 3.7282259393124537,
        'yaxis.range[1]': 9.59742380103187
    }
    """
    if relayoutData and "xaxis.range[0]" in relayoutData:
        return json.dumps(relayoutData)
    else:
        return None


@callback(
    Output("umap-select-a", "value"),
    inputs=Input("umap-graph-a", "relayoutData"),
    prevent_initial_call=True,
)
def relayout_umap_a(
    relayoutData,
):
    return _relayout_umap(relayoutData)


@callback(
    Output("umap-select-b", "value"),
    inputs=Input("umap-graph-b", "relayoutData"),
    prevent_initial_call=True,
)
def relayout_umap_b(
    relayoutData,
):
    return _relayout_umap(relayoutData)


@callback(
    Output("umap-a-container", "style"),
    Output("umap-b-container", "style"),
    inputs=Input("show-umap", "value"),
    prevent_initial_call=True,
)
def show_umap(
    show_umap,
):

    style = {} if show_umap else {"display": "none"}

    return style, style


@callback(
    Output("cytoscape-a", "stylesheet"),
    Output("cytoscape-b", "stylesheet"),
    inputs=[
        Input("show-network-weights", "value"),
        Input("show-population-fractions", "value"),
    ],
    state=State("cytoscape-a", "stylesheet"),
    prevent_initial_call=True,
)
def show_network_weights_callback(
    show_network_weights, show_population_fractions, stylesheet
):

    if ctx.triggered_id == "show-network-weights":
        label_value = "data(label)" if show_network_weights else ""

        for i, el in enumerate(stylesheet):
            if el["selector"] == "edge":
                stylesheet[i] = {**el, "style": {**el["style"], "label": label_value}}

    elif ctx.triggered_id == "show-population-fractions":
        label_value = (
            "data(label_with_size)" if show_population_fractions else "data(label)"
        )
        for i, el in enumerate(stylesheet):
            if el["selector"] == "node":
                stylesheet[i] = {**el, "style": {**el["style"], "label": label_value}}

    return (stylesheet, stylesheet)


@callback(
    Output("sender-select", "value"),
    Output("receiver-select", "value"),
    Output("view-radio", "value"),
    Input("cytoscape-a", "tapEdgeData"),
    Input("cytoscape-b", "tapEdgeData"),
    State("sender-select", "value"),
    State("receiver-select", "value"),
    State("view-radio", "value"),
    prevent_initial_call=True,
)
def cluster_edge_callback(
    cs_down_data, cs_up_data, sender_select, receiver_select, view_radio
):
    data = cs_down_data or cs_up_data
    if data:
        return (
            update_filter_value([], data["source"]),
            update_filter_value([], data["target"]),
            "sankey",
        )
    else:
        return sender_select, receiver_select, view_radio


@callback(
    Output(
        "ligand-select",
        "value",
    ),
    Output(
        "receptor-select",
        "value",
    ),
    Output(
        "em-select",
        "value",
    ),
    Output(
        "target-select",
        "value",
    ),
    Input("sankey-a", "clickData"),
    Input("sankey-b", "clickData"),
    State("ligand-select", "value"),
    State("receptor-select", "value"),
    State("em-select", "value"),
    State("target-select", "value"),
    prevent_initial_call=True,
    allow_duplicate=True,
)
def update_filters_click_node(
    click_data_a,
    click_data_b,
    ligand_select,
    receptor_select,
    em_select,
    target_select,
):

    def _update(current, new):
        return list(set(current + [new]) if isinstance(current, list) else set([new]))

    click_data = click_data_a or click_data_b
    if click_data:
        try:
            customdata = click_data["points"][0]["customdata"]
            node_label = customdata.split("_")[0]
            node_type = customdata.split("_")[1]
            if node_type == "ligand":
                ligand_select = _update(ligand_select, node_label)
            elif node_type == "receptor":
                receptor_select = _update(receptor_select, node_label)
            elif node_type == "em":
                em_select = _update(em_select, node_label)
            elif node_type == "target":
                target_select = _update(target_select, node_label)
        except Exception as e:
            pass

    return (
        ligand_select,
        receptor_select,
        em_select,
        target_select,
    )


@callback(
    output=dict(
        ligand_select=Output("ligand-select", "value", allow_duplicate=True),
        receptor_select=Output("receptor-select", "value", allow_duplicate=True),
        em_select=Output("em-select", "value", allow_duplicate=True),
        target_select=Output("target-select", "value", allow_duplicate=True),
        sender_select=Output("sender-select", "value", allow_duplicate=True),
        receiver_select=Output("receiver-select", "value", allow_duplicate=True),
        kinase_select=Output("kinase-select", "value"),
        any_role_select=Output("any-role-select", "value"),
        sigprob=Output({"type": "numerical-filter", "index": "sigprob"}, "value"),
        p_value=Output({"type": "numerical-filter", "index": "p-value"}, "value"),
        tpds=Output({"type": "numerical-filter", "index": "tpds"}, "value"),
        ppds=Output({"type": "numerical-filter", "index": "ppds"}, "value"),
    ),
    inputs=[Input("reset", "n_clicks")],
    prevent_initial_call=True,
)
def update_filters_click_node(nclicks):

    return filter_defaults()


@callback(
    Output("download-dataframe-a-csv", "data"),
    Output("download-dataframe-b-csv", "data"),
    inputs=dict(
        n_clicks=Input("btn_csv", "n_clicks"),
    ),
    state=dict(
        pcf=pathway_component_filter_inputs(state=True),
        sliders_container_children=State("allSlidersContainer", "children"),
    ),
    prevent_initial_call=True,
)
def download(
    n_clicks: int,
    pcf: dict,
    sliders_container_children,
):

    incytr_input = current_app.config["INCYTR_INPUT"]

    if n_clicks and n_clicks > 0:

        slider_values = parse_slider_values_from_tree(sliders_container_children)
        pf = PathwaysFilter(
            all_paths=incytr_input.paths,
            group_a_name=incytr_input.group_a,
            group_b_name=incytr_input.group_b,
            filter_afc_direction=pcf.get("restrict_afc"),
            filter_umap_a=parse_umap_filter_data(pcf.get("umap_select_a")),
            filter_umap_b=parse_umap_filter_data(pcf.get("umap_select_b")),
            filter_senders=pcf.get("sender_select"),
            filter_receivers=pcf.get("receiver_select"),
            filter_ligands=pcf.get("ligand_select"),
            filter_receptors=pcf.get("receptor_select"),
            filter_em=pcf.get("em_select"),
            filter_target_genes=pcf.get("target_select"),
            filter_all_molecules=pcf.get("any_role_select"),
            ppds_bounds=incytr_input.has_ppds and slider_values.get("ppds"),
            sp_threshold=slider_values.get("sigprob"),
            tppds_bounds=incytr_input.has_tpds and slider_values.get("tpds"),
            pval_threshold=incytr_input.has_p_value and slider_values.get("p-value"),
        )

        a_pathways = pf.filter("a", should_filter_umap=incytr_input.has_umap)
        b_pathways = pf.filter("b", should_filter_umap=incytr_input.has_umap)

        return (
            dcc.send_data_frame(a_pathways.to_csv, f"{incytr_input.group_a}.csv"),
            dcc.send_data_frame(b_pathways.to_csv, f"{incytr_input.group_b}.csv"),
        )


@callback(
    Output("modal", "is_open"),
    [Input("open", "n_clicks"), Input("close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open
