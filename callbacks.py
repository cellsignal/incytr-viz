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

    cmap = plt.get_cmap("tab20")

    # rgba arrays, values 0-1
    plt_colors = cmap(np.linspace(0, 1, len(clusters)))
    clusters["rgb_colors"] = [[int(x * 256) for x in c[0:3]] for c in plt_colors]

    def _add_node(row: pd.Series, pop_colname: str, total_cells: int) -> dict:

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

    total_cells = clusters["population_a"].sum() + clusters["population_b"].sum()

    if group == "a":
        return list(
            clusters.apply(
                lambda row: _add_node(row, "population_a", total_cells),
                axis=1,
            ).dropna()
        )
    elif group == "b":
        return list(
            clusters.apply(
                lambda row: _add_node(row, "population_b", total_cells),
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
            Output("rnas-a-hist", "figure"),
            Output("rnas-b-hist", "figure"),
            Output("fs-a-hist", "figure"),
            Output("fs-b-hist", "figure"),
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

        a_max_paths = np.max(
            filtered_pathways_a.groupby([get_cn("sender"), get_cn("receiver")]).size()
        )
        b_max_paths = np.max(
            filtered_pathways_b.groupby([get_cn("sender"), get_cn("receiver")]).size()
        )
        if np.isnan(a_max_paths):
            a_max_paths = 0
        if np.isnan(b_max_paths):
            b_max_paths = 0

        global_max_paths = max(a_max_paths, b_max_paths)

        def _get_group_figures(
            filtered_pathways: pd.DataFrame, global_max_paths: int, group
        ):

            suffix = get_group_name(filtered_pathways, group)

            if view_radio == "network":

                nodes = load_nodes(clusters, group)

                edges = load_edges(nodes, filtered_pathways, global_max_paths)

                cytoscape = get_cytoscape_component(
                    f"cytoscape-{group}", suffix, nodes + edges
                )

                graph = cytoscape

            elif view_radio == "sankey":
                sankey = get_sankey_component(
                    filtered_pathways, "sankey-{group}", suffix
                )

                graph = sankey

            sw_hist = get_hist(
                filtered_pathways, CN.SIGWEIGHT(full_pathways, group), "sigweight"
            )
            rnas_hist = get_hist(filtered_pathways, get_cn("rna_score"), "rna_score")
            fs_hist = get_hist(filtered_pathways, get_cn("final_score"), "final_score")
            pval_hist = get_hist(
                filtered_pathways, CN.PVAL(full_pathways, group), "p_val"
            )

            return (
                graph,
                sw_hist,
                pval_hist,
                rnas_hist,
                fs_hist,
            )

        group_a = _get_group_figures(filtered_pathways_a, global_max_paths, "a")

        group_b = _get_group_figures(filtered_pathways_b, global_max_paths, "b")

        return (
            [group_a[0], group_b[0]],
            group_a[1],
            group_b[1],
            group_a[2],
            group_b[2],
            group_a[3],
            group_b[3],
            group_a[4],
            group_b[4],
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
