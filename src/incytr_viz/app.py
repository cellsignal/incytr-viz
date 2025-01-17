import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, ALL, ctx, callback
from incytr_viz.modal_content import content
from incytr_viz.components import create_hist_figure, filter_container, slider_container
from flask_caching import Cache
import pandas as pd
from incytr_viz.util import create_logger
from incytr_viz.dtypes import clusters_dtypes, pathways_dtypes
from tabulate import tabulate
import os
import pdb

import numpy as np
from typing import Optional
import matplotlib.pyplot as plt

from dash.dependencies import Input, Output, State

from incytr_viz.components import (
    cytoscape_container,
    sankey_container,
)


from incytr_viz.util import (
    edge_width_map,
    update_filter_value,
    parse_slider_values_from_tree,
    parse_umap_filter_data,
    log_base,
    PathwaysFilter,
)


logger = create_logger(__name__)


app = Dash(
    __name__,
    suppress_callback_exceptions=True,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
)


cache = Cache(app.server, config={"CACHE_TYPE": "SimpleCache"})


def format_headers(headers):
    return (
        headers.str.strip()
        .str.lower()
        .str.replace(" ", "_")
        .str.replace("sender.group", "sender")
        .str.replace("receiver.group", "receiver")
    )


# @cache.cached(timeout=None, key_prefix="clusters")
@cache.memoize(timeout=None)
def get_clusters(fpath):
    if ".csv" in fpath:
        sep = ","
    elif ".tsv" in fpath:
        sep = "\t"
    else:
        raise ValueError(
            f"Pathways file suffix must be in [.csv,.tsv] -- check filename: {fpath}"
        )

    logger.info(
        "Loading cluster populations from {} as {}".format(
            fpath, {"\t": "TSV", ",": "CSV"}[sep]
        )
    )

    df = pd.read_csv(fpath, dtype=clusters_dtypes, sep=sep, compression="infer")

    df.columns = df.columns.str.lower().str.strip()
    if not all(c in df.columns for c in clusters_dtypes.keys()):
        raise ValueError(
            f"Invalid cell populations file: ensure the following columns are present: {clusters_dtypes.keys()}"
        )

    df = df[list(clusters_dtypes.keys())].reset_index(drop=True)

    df["type"] = df["type"].str.strip().str.lower()
    df["group"] = df["condition"].str.strip().str.lower()
    df = df.set_index("type")

    df["population"] = df["population"].fillna(0)
    df["pop_min_ratio"] = df["population"] / (
        df[df["population"] > 0]["population"].min()
    )

    df.drop(columns=["condition"], inplace=True)
    # assign colors to each cell type
    cmap = plt.get_cmap("tab20")

    cell_types = df.index.unique()

    plt_colors = cmap(np.linspace(0, 1, len(cell_types)))

    # 256 would not be websafe value
    rgb_colors = [[int(x * 255) for x in c[0:3]] for c in plt_colors]

    colors = {t: rgb_colors[i] for i, t in enumerate(cell_types)}
    df["color"] = df.index.map(colors)
    df["color"] = df["color"].apply(lambda x: f"rgb({x[0]},{x[1]},{x[2]})")

    if len(df["group"].unique()) != 2:
        raise ValueError(
            f"Expected exactly 2 groups in cluster populations file, found {len(df['group'].unique())}"
        )
    return df, df["group"].unique()


def parse_pathway_headers(headers, group_a, group_b):

    formatted = format_headers(headers)
    mapper = list(zip(headers, formatted))

    required = [
        "path",
        "sender",
        "receiver",
        "afc",
        "sigprob_" + group_a,
        "sigprob_" + group_b,
    ]

    optional = [
        "p_value_" + group_a,
        "p_value_" + group_b,
        "tprs",
        "prs",
        "kinase_r_of_em",
        "kinase_r_of_t",
        "kinase_em_of_t",
        "umap1",
        "umap2",
    ]

    logger.info("scanning pathways file for required and optional columns")

    required_df = pd.DataFrame.from_dict(
        {"colname": required, "required": True, "found": False}
    )
    optional_df = pd.DataFrame.from_dict(
        {"colname": optional, "required": False, "found": False}
    )
    columns_df = pd.concat([required_df, optional_df], axis=0)

    for row in columns_df.iterrows():
        col = row[1]["colname"]
        columns_df.loc[columns_df["colname"] == col, "found"] = col in formatted

    logger.info(
        "Pathways file column summary\n"
        + tabulate(columns_df, headers="keys", tablefmt="fancy_grid", showindex=False)
    )

    if (columns_df["required"] & ~columns_df["found"]).any():
        raise ValueError(
            f"Required columns not found in pathways file: {columns_df[columns_df['required'] & ~columns_df['found']]['colname'].values}"
        )

    if (~columns_df["found"] & ~columns_df["required"]).any():
        logger.warning(
            f"Optional columns missing in pathways file: {columns_df[~columns_df['found'] & ~columns_df['required']]['colname'].values}"
        )

    return [
        x[0]
        for x in mapper
        if x[1] in columns_df[columns_df["found"]]["colname"].values
    ]


class PathwayInput:

    def __init__(self, group_a, group_b, paths):
        self.group_a = group_a
        self.group_b = group_b
        self.paths = paths
        self.has_tprs = "tprs" in self.paths.columns
        self.has_prs = "prs" in self.paths.columns
        self.has_p_value = all(
            x in self.paths.columns
            for x in ["p_value_" + self.group_a, "p_value_" + self.group_b]
        )
        self.has_umap = all(x in self.paths.columns for x in ["umap1", "umap2"])
        self.unique_senders = self.paths["sender"].unique()
        self.unique_receivers = self.paths["receiver"].unique()
        self.unique_ligands = self.paths["ligand"].unique()
        self.unique_receptors = self.paths["receptor"].unique()
        self.unique_em = self.paths["em"].unique()
        self.unique_targets = self.paths["target"].unique()


# @cache.cached(timeout=None, key_prefix="pathways")
@cache.memoize(timeout=None)
def get_pathways(fpath, group_a, group_b):
    if ".csv" in fpath:
        sep = ","
    elif ".tsv" in fpath:
        sep = "\t"
    else:
        raise ValueError(
            f"Pathways file suffix must be in [.csv,.tsv] -- check filename {fpath}"
        )

    logger.info(
        "Loading pathways from {} as {}".format(fpath, {"\t": "TSV", ",": "CSV"}[sep])
    )

    headers = pd.read_csv(fpath, nrows=0, sep=sep).columns
    to_keep = parse_pathway_headers(headers, group_a, group_b)

    paths = pd.read_csv(fpath, dtype=pathways_dtypes, usecols=to_keep, sep=sep)
    paths.columns = format_headers(paths.columns)

    paths = paths.astype(
        {k: v for k, v in pathways_dtypes.items() if k in paths.columns}
    )

    num_invalid = 0

    incomplete_paths = paths["path"].str.strip().str.split("*").str.len() != 4

    if incomplete_paths.sum() > 0:
        logger.warning(
            f"{incomplete_paths.sum()} rows with invalid pathway format found. Expecting form L*R*EM*T"
        )
        logger.warning("First 10 invalid paths:")
        logger.warning(paths[incomplete_paths]["path"].head().values)

    num_invalid += len(incomplete_paths)

    paths = paths.loc[~incomplete_paths]

    paths["ligand"] = paths["path"].str.split("*").str[0].str.strip()
    paths["receptor"] = paths["path"].str.split("*").str[1].str.strip()
    paths["em"] = paths["path"].str.split("*").str[2].str.strip()
    paths["target"] = paths["path"].str.split("*").str[3].str.strip()
    paths["sender"] = paths["sender"].str.strip().str.lower()
    paths["receiver"] = paths["receiver"].str.strip().str.lower()
    paths["path"] = (
        paths["path"]
        .str.cat(paths["sender"], sep="*")
        .str.cat(paths["receiver"], sep="*")
    )

    duplicates_mask = paths.duplicated()
    if duplicates_mask.sum() > 0:
        logger.warning(f"{duplicates_mask.sum()} duplicate rows found")

    is_na_mask = (
        paths[["afc", "sigprob_" + group_a, "sigprob_" + group_b]].isna().any(axis=1)
    )
    if is_na_mask.sum() > 0:
        logger.info(
            f"{is_na_mask.sum()} rows with invalid values found in required columns"
        )

    invalid = duplicates_mask | is_na_mask

    if invalid.sum() > 0:
        logger.info(f"Removing {invalid.sum()} duplicate or invalid rows")

    paths = paths[~invalid].reset_index(drop=True)

    return PathwayInput(
        group_a=group_a,
        group_b=group_b,
        paths=paths,
    )


clusters, groups = get_clusters(os.environ["INCYTR_CLUSTERS"])
pi = get_pathways(os.environ["INCYTR_PATHWAYS"], groups[0], groups[1])

app.layout = html.Div(
    [
        html.Div(
            children=[
                dcc.Store(id="dummy-store", data=None),
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
                                            disabled=False,
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
        html.Div(
            slider_container(
                has_tprs=pi.has_tprs, has_prs=pi.has_prs, has_p_value=pi.has_p_value
            ),
            id="slider-container",
        ),
        html.Div(
            filter_container(
                sender=list(pi.unique_senders),
                receiver=list(pi.unique_receivers),
                ligand=list(pi.unique_ligands),
                receptor=list(pi.unique_receptors),
                em=list(pi.unique_em),
                target=list(pi.unique_targets),
            ),
            className="sidebar",
            id="filter-container",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.H3(
                                    pi.group_a.title(),
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
                            [],
                            className="umapContainer",
                            id="umap-a-container",
                            style={"display": "none"},
                        ),
                        html.Div(
                            [
                                html.Div([dcc.Graph()], id="figure-a-container"),
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
                                html.H3(
                                    pi.group_b.title(),
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
                            [],
                            className="umapContainer",
                            id="umap-b-container",
                            style={"display": "none"},
                        ),
                        html.Div(
                            [
                                html.Div([dcc.Graph()], id="figure-b-container"),
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
    ],
    id="app-container",
    className="app",
)


server = app.server
# def get_server(pathways_fpath, clusters_fpath):
#     os.environ["INCYTR_PATHWAYS"] = pathways_fpath
#     os.environ["INCYTR_CLUSTERS"] = clusters_fpath
#     return app.server


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
        node_population = row["population"]

        if (not node_population) or (np.isnan(node_population)):
            return np.nan

        data = dict()
        data["id"] = node_type
        data["label"] = node_type
        data["cluster_size"] = node_population
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
    node_labels = pd.Series([x["data"]["label"] for x in nodes])
    pathways = pathways[
        (pathways["sender"].isin(node_labels))
        & (pathways["receiver"].isin(node_labels))
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
            x["data"]["background_color"]
            for x in nodes
            if x["data"]["label"] == source_id
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

        return num_targets <= 200

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


# def store_data_inputs(state=False):

#     klass = State if state else Input
#     return dict(
#         has_tprs=klass("has-tprs", "data"),
#         has_prs=klass("has-prs", "data"),
#         has_p_value=klass("has-p-value", "data"),
#         has_umap=klass("has-umap", "data"),
#         group_a_name=klass("group-a-name", "data"),
#         group_b_name=klass("group-b-name", "data"),
#     )


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
    )


def network_style_inputs(state=False):
    klass = State if state else Input
    return dict(
        node_scale_factor=klass("node-scale-factor", "value"),
        edge_scale_factor=klass("edge-scale-factor", "value"),
        label_scale_factor=klass("label-scale-factor", "value"),
    )


@app.callback(
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
    state=dict(show_network_weights=State("show-network-weights", "value")),
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

    clusters, groups = get_clusters(os.environ["INCYTR_CLUSTERS"])
    pi = get_pathways(os.environ["INCYTR_PATHWAYS"], groups[0], groups[1])

    filter_umap_a = parse_umap_filter_data(pcf.get("umap_select_a"))
    filter_umap_b = parse_umap_filter_data(pcf.get("umap_select_b"))

    slider_values = parse_slider_values_from_tree(sliders_container_children)

    pf = PathwaysFilter(
        all_paths=pi.paths,
        group_a_name=pi.group_a,
        group_b_name=pi.group_b,
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
        prs_bounds=pi.has_prs and slider_values.get("prs"),
        sp_threshold=slider_values.get("sigprob"),
        tprs_bounds=pi.has_tprs and slider_values.get("tprs"),
        pval_threshold=pi.has_p_value and slider_values.get("p-value"),
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
            warn = ids and not any(x.endswith("target") for x in ids)

            sankey = sankey_container(
                clusters,
                ids,
                labels,
                source,
                target,
                value,
                color,
                "Sankey " + group_name,
                group_id,
                warn,
                color_flow=pcf.get("sankey_color_flow"),
            )

            graph_container = sankey
        return [
            graph_container,
            create_hist_figure(
                paths=filtered_group_paths,
                has_tprs=pi.has_tprs,
                has_prs=pi.has_prs,
                has_p_value=pi.has_p_value,
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
        group_name=pi.group_a,
        group_id="a",
    )
    group_b_figs = _get_group_figures(
        filtered_group_paths=b_pathways,
        clusters=clusters,
        global_max_paths=global_max_paths,
        group_name=pi.group_b,
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


# @app.callback(
#     Output("filter-container", "children"),
#     Input("dummy-store", "data"),
# )
# def update_filter_container(store_data):

#     _, groups = get_clusters(os.environ["INCYTR_CLUSTERS"])
#     pi = get_pathways(os.environ["INCYTR_PATHWAYS"], groups[0], groups[1])

#     return filter_container(
#         sender=pi.unique_senders,
#         receiver=pi.unique_receivers,
#         ligand=pi.unique_ligands,
#         receptor=pi.unique_receptors,
#         em=pi.unique_em,
#         target=pi.unique_targets,
#     )


# @app.callback(
#     Output("slider-container", "children"),
#     Input("filter-container", "id"),
# )
# def update_slider_container(id):

#     _, groups = get_clusters(os.environ["INCYTR_CLUSTERS"])
#     pi = get_pathways(os.environ["INCYTR_PATHWAYS"], groups[0], groups[1])

#     return slider_container(
#         has_tprs=pi.has_tprs, has_prs=pi.has_prs, has_p_value=pi.has_p_value
#     )


# def _relayout_umap(relayoutData):
#     """
#     {
#         'xaxis.range[0]': -0.6369249007630504,
#         'xaxis.range[1]': 6.965720316453904,
#         'yaxis.range[0]': 3.7282259393124537,
#         'yaxis.range[1]': 9.59742380103187
#     }
#     """
#     if relayoutData and "xaxis.range[0]" in relayoutData:
#         return json.dumps(relayoutData)
#     else:
#         return None


# @app.callback(
#     Output("umap-select-a", "value"),
#     inputs=Input("scatter-plot-a", "relayoutData"),
#     prevent_initial_call=True,
# )
# def relayout_umap_a(
#     relayoutData,
# ):
#     return _relayout_umap(relayoutData)


# @app.callback(
#     Output("umap-select-b", "value"),
#     inputs=Input("scatter-plot-b", "relayoutData"),
#     prevent_initial_call=True,
# )
# def relayout_umap_b(
#     relayoutData,
# ):
#     return _relayout_umap(relayoutData)


@app.callback(
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


@app.callback(
    Output("cytoscape-a", "stylesheet"),
    Output("cytoscape-b", "stylesheet"),
    inputs=Input("show-network-weights", "value"),
    state=State("cytoscape-a", "stylesheet"),
    prevent_initial_call=True,
)
def show_network_weights_callback(show_network_weights, stylesheet):

    label_value = "data(label)" if show_network_weights else ""

    for i, el in enumerate(stylesheet):
        if el["selector"] == "edge":
            stylesheet[i] = {**el, "style": {**el["style"], "label": label_value}}
    return (stylesheet, stylesheet)


@app.callback(
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


@app.callback(
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


@app.callback(
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

    clusters, groups = get_clusters(os.environ["INCYTR_CLUSTERS"])
    pi = get_pathways(os.environ["INCYTR_PATHWAYS"], groups[0], groups[1])

    if n_clicks and n_clicks > 0:

        slider_values = parse_slider_values_from_tree(sliders_container_children)
        pf = PathwaysFilter(
            all_paths=pi.paths,
            group_a_name=pi.group_a,
            group_b_name=pi.group_b,
            filter_umap_a=parse_umap_filter_data(pcf.get("umap_select_a")),
            filter_umap_b=parse_umap_filter_data(pcf.get("umap_select_b")),
            filter_senders=pcf.get("sender_select"),
            filter_receivers=pcf.get("receiver_select"),
            filter_ligands=pcf.get("ligand_select"),
            filter_receptors=pcf.get("receptor_select"),
            filter_em=pcf.get("em_select"),
            filter_target_genes=pcf.get("target_select"),
            filter_all_molecules=pcf.get("any_role_select"),
            prs_bounds=pi.has_prs and slider_values.get("prs"),
            sp_threshold=slider_values.get("sigprob"),
            tprs_bounds=pi.has_tprs and slider_values.get("tprs"),
            pval_threshold=pi.has_p_value and slider_values.get("p-value"),
        )

        a_pathways = pf.filter("a", should_filter_umap=pi.has_umap)
        b_pathways = pf.filter("b", should_filter_umap=pi.has_umap)

        return (
            dcc.send_data_frame(a_pathways.to_csv, f"{pi.group_a}.csv"),
            dcc.send_data_frame(b_pathways.to_csv, f"{pi.group_b}.csv"),
        )


# @app.callback(
#     Output("modal", "is_open"),
#     [Input("open", "n_clicks"), Input("close", "n_clicks")],
#     [State("modal", "is_open")],
# )
# def toggle_modal(n1, n2, is_open):
#     if n1 or n2:
#         return not is_open
#     return is_open


# @app.callback(
#     Output("umap", "is_open"),
#     Input("show-umap", "value"),
# )
# def toggle_umap(n1, n2, is_open):
#     if n1 or n2:
#         return not is_open
#     return is_open
