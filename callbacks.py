import numpy as np
from dash import Dash, html
from dash.dependencies import Input, Output, State
from util import *
from components import (
    load_edges,
    load_nodes,
    get_cytoscape_component,
    get_sankey_component,
    get_hist,
)


def apply_filter_callback(app, full_pathways, full_clusters):
    @app.callback(
        Output("figures-container", "children"),
        Output("sw-hist", "figure"),
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
    def filter_figures(
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

        up_pathways = filtered_pathways[filtered_pathways["final_score"] > 0]
        down_pathways = filtered_pathways[filtered_pathways["final_score"] < 0]

        if view_radio == "network":
            global_max_paths = np.max(
                filtered_pathways.groupby(["Sender.group", "Receiver.group"]).size()
            )

            nodes = load_nodes(full_clusters)
            up_edges = load_edges(nodes, up_pathways, global_max_paths)
            down_edges = load_edges(nodes, down_pathways, global_max_paths)

            up_cytoscape = get_cytoscape_component("cytoscape-up", nodes + up_edges)
            down_cytoscape = get_cytoscape_component(
                "cytoscape-down", nodes + down_edges
            )

            graphs = html.Div(
                children=[up_cytoscape, down_cytoscape],
                id="cytoscape-container",
                style={"display": "flex", "flexDirection": "row", "width": "100%"},
            )

        elif view_radio == "pathways":
            up_sankey = get_sankey_component(up_pathways, "sankey-up")
            down_sankey = get_sankey_component(down_pathways, "sankey-down")

            graphs = html.Div(
                children=[up_sankey, down_sankey],
                id="sankey-container",
                style={"width": "100%"},
            )

        sw_hist = get_hist(filtered_pathways, "SigWeight_WT", 20)
        rnas_hist = get_hist(filtered_pathways, "adjlog2FC", 20)
        fs_hist = get_hist(filtered_pathways, "final_score", 20)

        return [graphs, sw_hist, rnas_hist, fs_hist]

    return app


def apply_cluster_edge_callback(app):
    @app.callback(
        Output("sender-select", "value"),
        Output("receiver-select", "value"),
        Output("view-radio", "value"),
        Input("cytoscape-down", "tapEdgeData"),
        Input("cytoscape-up", "tapEdgeData"),
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
        Input("sankey-up", "clickData"),
        Input("sankey-down", "clickData"),
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
