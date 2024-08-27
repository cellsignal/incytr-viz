import numpy as np
from dash import Dash, html
from dash.dependencies import Input, Output, State
from util import *
from components import (
    get_cytoscape_component,
    get_sankey_component,
    get_hist,
)
import matplotlib.pyplot as plt


def has_rna_score(full_pathways):
    return


def load_nodes(clusters: pd.DataFrame) -> list[dict]:
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

    cmap = plt.get_cmap("viridis")

    # rgba arrays, values 0-1
    plt_colors = cmap(np.linspace(0, 1, len(clusters)))
    clusters["rgb_colors"] = [[int(x * 256) for x in c[0:3]] for c in plt_colors]

    def _nodes_from_clusters_file(
        row: pd.Series, pop_colname: str, total_cells: int
    ) -> dict:

        node_type = row.name
        node_population = row[pop_colname]
        node_rgb_color = row["rgb_colors"]

        if (not node_population) or (np.isnan(node_population)):
            return np.nan

        data = dict()
        data["id"] = node_type
        data["label"] = node_type
        data["cluster_size"] = node_population
        data["width"] = node_size_map(node_population, total_cells)
        data["height"] = node_size_map(node_population, total_cells)
        data["background_color"] = "rgb({}, {}, {})".format(*node_rgb_color)
        return {"data": data}

    clusters["nodes_a"] = clusters.apply(
        lambda row: _nodes_from_clusters_file(
            row, "population_a", clusters["population_a"].sum()
        ),
        axis=1,
    )
    clusters["nodes_b"] = clusters.apply(
        lambda row: _nodes_from_clusters_file(
            row, "population_b", clusters["population_b"].sum()
        ),
        axis=1,
    )

    # omit nodes that are na -- these nodes do not exist or have pop 0
    return list(clusters["nodes_a"].dropna()), list(clusters["nodes_b"].dropna())


def load_edges(
    nodes: list[dict],
    pathways: pd.DataFrame,
    sigweight_filter_column,
    sighweight_threshold: float,
    global_max_paths: int,
):
    """add pathways from source to target"""
    edges = []

    ## filter pathways if sender/receiver not in nodes
    node_labels = pd.Series([x["data"]["label"] for x in nodes])
    pathways = pathways[
        (pathways["Sender"].isin(node_labels))
        & (pathways["Receiver"].isin(node_labels))
    ]

    ## filter pathways that are below sigweight threshold
    pathways = pathways[pathways[sigweight_filter_column] >= sighweight_threshold]

    if len(pathways) == 0:
        return edges

    s: pd.Series = pathways.groupby(["Sender", "Receiver"]).size()

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
                abs(e["data"]["weight"]), global_max_paths
            )

    return edges


def apply_filter_callback(
    app, full_pathways, full_clusters, has_rna_score, has_final_score
):
    @app.callback(
        Output("figures-container", "children"),
        Output("sw-a-hist", "figure"),
        Output("sw-b-hist", "figure"),
        Output("rnas-hist", "figure"),
        Output("fs-hist", "figure"),
        Input("sender-select", "value"),
        Input("receiver-select", "value"),
        Input("ligand-select", "value"),
        Input("receptor-select", "value"),
        Input("em-select", "value"),
        Input("target-select", "value"),
        Input("all-molecules-select", "value"),
        Input("fs-slider", "value"),
        Input("sw-slider", "value"),
        Input("rnas-slider", "value"),
        Input("view-radio", "value"),
    )
    def update_figures(
        sender_select,
        receiver_select,
        ligand_select,
        receptor_select,
        em_select,
        target_select,
        all_mols,
        fs_threshold,
        sw_threshold,
        rnas_threshold,
        view_radio,
    ):

        filtered_pathways = filter_pathways(
            full_pathways,
            filter_senders=sender_select,
            filter_receivers=receiver_select,
            filter_ligands=ligand_select,
            filter_receptors=receptor_select,
            filter_em=em_select,
            filter_target_genes=target_select,
            filter_all_molecules=all_mols,
            fs_threshold=fs_threshold,
            sw_threshold=sw_threshold,
            rnas_threshold=rnas_threshold,
        )

        if view_radio == "network":
            global_max_paths = np.max(
                filtered_pathways.groupby(["Sender", "Receiver"]).size()
            )

            nodes_a, nodes_b = load_nodes(full_clusters)
            edges_a = load_edges(
                nodes_a,
                filtered_pathways,
                CN.SIGWEIGHT_A(filtered_pathways),
                sw_threshold,
                global_max_paths,
            )
            edges_b = load_edges(
                nodes_b,
                filtered_pathways,
                CN.SIGWEIGHT_B(filtered_pathways),
                sw_threshold,
                global_max_paths,
            )

            cytoscape_a = get_cytoscape_component(
                "cytoscape-a", "Exp. Condition", nodes_a + edges_a
            )
            cytoscape_b = get_cytoscape_component(
                "cytoscape-b", "WT Condition", nodes_b + edges_b
            )

            graphs = html.Div(
                children=[cytoscape_a, cytoscape_b],
                id="cytoscape-container",
                style={"display": "flex", "flexDirection": "row", "width": "100%"},
            )

        elif view_radio == "pathways":
            sankey_a = get_sankey_component(
                filtered_pathways,
                "sankey-a",
                CN.SIGWEIGHT_A(filtered_pathways),
                sw_threshold,
                "Exp. Condition",
            )
            sankey_b = get_sankey_component(
                filtered_pathways,
                "sankey-b",
                CN.SIGWEIGHT_B(filtered_pathways),
                sw_threshold,
                "WT Condition",
            )

            graphs = html.Div(
                children=[sankey_a, sankey_b],
                id="sankey-container",
                style={"width": "100%"},
            )

        sw_a_hist = get_hist(filtered_pathways, CN.SIGWEIGHT_A(filtered_pathways), 20)
        sw_b_hist = get_hist(filtered_pathways, CN.SIGWEIGHT_B(filtered_pathways), 20)
        rnas_hist = (
            get_hist(filtered_pathways, get_cn("rna_score"), 20)
            if has_rna_score
            else empty_hist(get_cn("rna_score"))
        )
        fs_hist = (
            get_hist(filtered_pathways, get_cn("final_score"), 20)
            if has_final_score
            else empty_hist(get_cn("final_score"))
        )

        return [graphs, sw_a_hist, sw_b_hist, rnas_hist, fs_hist]

    return app


def apply_cluster_edge_callback(app):
    @app.callback(
        Output("sender-select", "value"),
        Output("receiver-select", "value"),
        Output("view-radio", "value"),
        Input("cytoscape-b", "tapEdgeData"),
        Input("cytoscape-a", "tapEdgeData"),
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
                "pathways",
            )
        else:
            return sender_select, receiver_select, view_radio

    return app


def apply_sankey_callbacks(
    app: Dash,
):
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
        click_data_up,
        click_data_down,
        ligand_select,
        receptor_select,
        em_select,
        target_select,
    ):

        def _update(current, new):
            return list(
                set(current + [new]) if isinstance(current, list) else set([new])
            )

        click_data = click_data_up or click_data_down

        if click_data:

            try:
                customdata = click_data["points"][0]["customdata"]
                node_label = customdata.split("_")[0]
                node_type = customdata.split("_")[1]
                if node_type == "Ligand":
                    ligand_select = _update(ligand_select, node_label)
                elif node_type == "Receptor":
                    receptor_select = _update(receptor_select, node_label)
                elif node_type == "EM":
                    em_select = _update(em_select, node_label)
                elif node_type == "Target":
                    target_select = _update(target_select, node_label)
            except Exception as e:
                print(e)

        return (
            ligand_select,
            receptor_select,
            em_select,
            target_select,
        )
