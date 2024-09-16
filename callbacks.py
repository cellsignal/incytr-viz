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
from data import filter_pathways, load_nodes, load_edges


def apply_filter_callback(
    app, full_pathways, clusters, has_rna_score, has_final_score, has_p_value
):

    @app.callback(
        Output("group-a-container", "children"),
        Output("group-b-container", "children"),
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
                    filtered_pathways, f"sankey-{group}", suffix
                )

                graph = sankey

            sw_hist = get_hist(
                filtered_pathways,
                CN.SIGWEIGHT(full_pathways.columns, group),
                "sigweight",
            )
            rnas_hist = get_hist(filtered_pathways, get_cn("rna_score"), "rna_score")
            fs_hist = get_hist(filtered_pathways, get_cn("final_score"), "final_score")
            pval_hist = get_hist(
                filtered_pathways, CN.PVAL(full_pathways.columns, group), "p_val"
            )

            return [
                html.Div(
                    [
                        sw_hist,
                        pval_hist,
                        rnas_hist,
                        fs_hist,
                    ],
                    className="histContainer",
                ),
                graph,
            ]

        return (
            _get_group_figures(filtered_pathways_a, global_max_paths, "a"),
            _get_group_figures(filtered_pathways_b, global_max_paths, "b"),
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
