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


def load_nodes(clusters: pd.DataFrame, group) -> list[dict]:
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

    if group == "a":
        return list(
            clusters.apply(
                lambda row: _nodes_from_clusters_file(
                    row, "population_a", clusters["population_a"].sum()
                ),
                axis=1,
            ).dropna()
        )
    elif group == "b":
        return list(
            clusters.apply(
                lambda row: _nodes_from_clusters_file(
                    row, "population_b", clusters["population_b"].sum()
                ),
                axis=1,
            ).dropna()
        )


def load_edges(
    nodes: list[dict],
    pathways: pd.DataFrame,
    global_max_paths: int,
):
    """add pathways from source to target"""
    edges = []

    ## filter pathways if sender/receiver not in nodes
    node_labels = pd.Series([x["data"]["label"] for x in nodes])
    pathways = pathways[
        (pathways[get_cn("sender")].isin(node_labels))
        & (pathways[get_cn("receiver")].isin(node_labels))
    ]

    ## filter pathways that are below sigweight threshold

    if len(pathways) == 0:
        return edges

    s: pd.Series = pathways.groupby([get_cn("sender"), get_cn("receiver")]).size()

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
    app, full_pathways, clusters, has_rna_score, has_final_score, has_p_value
):

    @app.callback(
        output=[
            Output("figures-container", "children"),
            Output("sw-a-hist", "figure"),
            Output("sw-b-hist", "figure"),
            Output("pval-a-hist", "figure"),
            Output("pval-b-hist", "figure"),
            Output("rnas-hist", "figure"),
            Output("fs-hist", "figure"),
        ],
        inputs=dict(
            pcf=pathway_component_filter_inputs(),
            pvf=pathway_value_filter_inputs(
                has_rna_score, has_final_score, has_p_value
            ),
            view_radio=view_radio_input(),
        ),
    )
    def update_filtered_elements(pcf, pvf, view_radio):

        sender_select = pcf.get("sender_select", None)
        receiver_select = pcf.get("receiver_select", None)
        ligand_select = pcf.get("ligand_select", None)
        receptor_select = pcf.get("receptor_select", None)
        em_select = pcf.get("em_select", None)
        target_select = pcf.get("target_select", None)
        all_mols_select = pcf.get("all_mols_select", None)
        sw_threshold = pvf.get("sw_threshold", None)
        pval_threshold = pvf.get("pval_threshold", None)
        rnas_bounds = pvf.get("rnas_bounds", None)
        fs_bounds = pvf.get("fs_bounds", None)

        filtered_pathways_a = filter_pathways(
            full_pathways=full_pathways,
            group="a",
            filter_senders=sender_select,
            filter_receivers=receiver_select,
            filter_ligands=ligand_select,
            filter_receptors=receptor_select,
            filter_em=em_select,
            filter_target_genes=target_select,
            filter_all_molecules=all_mols_select,
            fs_bounds=fs_bounds,
            sw_threshold=sw_threshold,
            rnas_bounds=rnas_bounds,
            pval_threshold=pval_threshold,
        )

        filtered_pathways_b = filter_pathways(
            full_pathways,
            "b",
            filter_senders=sender_select,
            filter_receivers=receiver_select,
            filter_ligands=ligand_select,
            filter_receptors=receptor_select,
            filter_em=em_select,
            filter_target_genes=target_select,
            filter_all_molecules=all_mols_select,
            fs_bounds=fs_bounds,
            sw_threshold=sw_threshold,
            rnas_bounds=rnas_bounds,
            pval_threshold=pval_threshold,
        )

        filtered_pathways_total = pd.concat(
            [filtered_pathways_a, filtered_pathways_b]
        ).drop_duplicates(subset=[get_cn("path")])

        if view_radio == "network":
            a_max_paths = np.max(
                filtered_pathways_a.groupby(
                    [get_cn("sender"), get_cn("receiver")]
                ).size()
            )
            b_max_paths = np.max(
                filtered_pathways_b.groupby(
                    [get_cn("sender"), get_cn("receiver")]
                ).size()
            )
            global_max_paths = max(a_max_paths, b_max_paths)

            a_nodes = load_nodes(clusters, "a")
            b_nodes = load_nodes(clusters, "b")

            a_edges = load_edges(a_nodes, filtered_pathways_a, global_max_paths)
            b_edges = load_edges(b_nodes, filtered_pathways_b, global_max_paths)

            cytoscape_a = get_cytoscape_component(
                "cytoscape-a", get_group_name(full_pathways, "a"), a_nodes + a_edges
            )
            cytoscape_b = get_cytoscape_component(
                "cytoscape-b", get_group_name(full_pathways, "b"), b_nodes + b_edges
            )

            graph_a, graph_b = cytoscape_a, cytoscape_b
        elif view_radio == "sankey":
            sankey_a = get_sankey_component(
                filtered_pathways_a, "sankey-a", get_group_name(full_pathways, "a")
            )
            sankey_b = get_sankey_component(
                filtered_pathways_b, "sankey-b", get_group_name(full_pathways, "b")
            )

            graph_a, graph_b = sankey_a, sankey_b

        sw_a_hist = get_hist(filtered_pathways_a, CN.SIGWEIGHT(full_pathways, "a"), 20)

        sw_b_hist = get_hist(filtered_pathways_b, CN.SIGWEIGHT(full_pathways, "b"), 20)

        if has_rna_score:
            rnas_hist = (
                get_hist(filtered_pathways_total, get_cn("rna_score"), 20)
                if has_rna_score
                else empty_hist(get_cn("rna_score"))
            )
        else:
            rnas_hist = empty_hist(get_cn("rna_score"))
        if has_final_score:
            fs_hist = (
                get_hist(filtered_pathways_total, get_cn("final_score"), 20)
                if has_final_score
                else empty_hist(get_cn("final_score"))
            )
        else:
            fs_hist = empty_hist(get_cn("final_score"))
        if has_p_value:
            pval_a_hist = (
                get_hist(filtered_pathways_total, CN.PVAL(full_pathways, "a"), 20)
                if has_p_value
                else empty_hist(CN.PVAL(full_pathways, "a"))
            )
            pval_b_hist = (
                get_hist(filtered_pathways_total, CN.PVAL(full_pathways, "b"), 20)
                if has_p_value
                else empty_hist(CN.PVAL(full_pathways, "b"))
            )
        else:
            pval_a_hist = empty_hist(f"pval_{get_group_name(full_pathways, 'a')}")
            pval_b_hist = empty_hist(f"pval_{get_group_name(full_pathways, 'b')}")

        return (
            [graph_a, graph_b],
            sw_a_hist,
            sw_b_hist,
            pval_a_hist,
            pval_b_hist,
            rnas_hist,
            fs_hist,
        )


def apply_cluster_edge_callback(app):
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
                if node_type == get_cn("ligand"):
                    ligand_select = _update(ligand_select, node_label)
                elif node_type == get_cn("receptor"):
                    receptor_select = _update(receptor_select, node_label)
                elif node_type == get_cn("em"):
                    em_select = _update(em_select, node_label)
                elif node_type == get_cn("target"):
                    target_select = _update(target_select, node_label)
            except Exception as e:
                print(e)

        return (
            ligand_select,
            receptor_select,
            em_select,
            target_select,
        )
