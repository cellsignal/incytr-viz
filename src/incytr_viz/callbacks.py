from typing import Optional

import numpy as np
import pdb
import json
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt

from dash import Dash, ALL, dcc, ctx
from dash.dependencies import Input, Output, State

from incytr_viz.components import (
    cytoscape_container,
    sankey_container,
    hist_container,
)


from incytr_viz.util import (
    edge_width_map,
    update_filter_value,
    parse_slider_values_from_tree,
    parse_umap_filter_data,
    log_base,
    PathwaysFilter,
)


def load_nodes(
    clusters: pd.DataFrame, node_scale_factor, edge_size_factor
) -> list[dict]:
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

    ## filter pathways that are below sigweight threshold

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


def store_data_inputs(state=False):

    klass = State if state else Input
    return dict(
        has_rna=klass("has-rna", "data"),
        has_final=klass("has-final", "data"),
        has_p_value=klass("has-p-value", "data"),
        has_umap=klass("has-umap", "data"),
        group_a_name=klass("group-a-name", "data"),
        group_b_name=klass("group-b-name", "data"),
    )


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
    )


def apply_callbacks(app: Dash, all_pathways, clusters):

    @app.callback(
        Output("group-a-container", "children"),
        Output("group-b-container", "children"),
        Output("pathways-count-a", "children"),
        Output("pathways-count-b", "children"),
        inputs=dict(
            sdi=store_data_inputs(),
            pcf=pathway_component_filter_inputs(),
            nsi=network_style_inputs(),
            slider_changed=Input({"type": "numerical-filter", "index": ALL}, "value"),
            sliders_container_children=State("allSlidersContainer", "children"),
            view_radio=Input("view-radio", "value"),
        ),
        state=dict(show_network_weights=State("show-network-weights", "value")),
    )
    def update_figure_and_histogram(
        sdi,
        pcf,
        nsi,
        slider_changed,
        sliders_container_children,
        view_radio,
        show_network_weights,
    ):

        filter_umap_a = parse_umap_filter_data(pcf.get("umap_select_a"))
        filter_umap_b = parse_umap_filter_data(pcf.get("umap_select_b"))

        slider_values = parse_slider_values_from_tree(sliders_container_children)

        pf = PathwaysFilter(
            all_paths=all_pathways,
            group_a_name=sdi.get("group_a_name"),
            group_b_name=sdi.get("group_b_name"),
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
            fs_bounds=slider_values.get("final-score"),
            sw_threshold=slider_values.get("sigweight"),
            rnas_bounds=slider_values.get("rna-score"),
            pval_threshold=slider_values.get("p-value"),
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
                    edge_size_factor=nsi.get("edge_scale_factor", 0.1),
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
                hist_container(group_id, filtered_group_paths),
            ]

        a_max_paths = np.max(a_pathways.groupby(["sender", "receiver"]).size())
        b_max_paths = np.max(b_pathways.groupby(["sender", "receiver"]).size())

        if np.isnan(a_max_paths):
            a_max_paths = 0
        if np.isnan(b_max_paths):
            b_max_paths = 0

        global_max_paths = max(a_max_paths, b_max_paths)

        return (
            _get_group_figures(
                filtered_group_paths=a_pathways,
                clusters=clusters,
                global_max_paths=global_max_paths,
                group_name=sdi.get("group_a_name"),
                group_id="a",
            ),
            _get_group_figures(
                filtered_group_paths=b_pathways,
                clusters=clusters,
                global_max_paths=global_max_paths,
                group_name=sdi.get("group_b_name"),
                group_id="b",
            ),
            len(a_pathways),
            len(b_pathways),
        )

    def _umap_callback(relayoutData):
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

    @app.callback(
        Output("umap-select-a", "value"),
        inputs=Input("scatter-plot-a", "relayoutData"),
        prevent_initial_call=True,
    )
    def umap_callback_a(
        relayoutData,
    ):
        return _umap_callback(relayoutData)

    @app.callback(
        Output("umap-select-b", "value"),
        inputs=Input("scatter-plot-b", "relayoutData"),
        prevent_initial_call=True,
    )
    def umap_callback_b(
        relayoutData,
    ):
        return _umap_callback(relayoutData)

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
            return list(
                set(current + [new]) if isinstance(current, list) else set([new])
            )

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
            sdi=store_data_inputs(state=True),
            pcf=pathway_component_filter_inputs(state=True),
            sliders_container_children=State("allSlidersContainer", "children"),
        ),
        prevent_initial_call=True,
    )
    def download(
        n_clicks: int,
        sdi: dict,
        pcf: dict,
        sliders_container_children,
    ):

        if n_clicks and n_clicks > 0:

            slider_values = parse_slider_values_from_tree(sliders_container_children)
            pf = PathwaysFilter(
                all_paths=all_pathways,
                group_a_name=sdi.get("group_a_name"),
                group_b_name=sdi.get("group_b_name"),
                filter_umap_a=parse_umap_filter_data(pcf.get("umap_select_a")),
                filter_umap_b=parse_umap_filter_data(pcf.get("umap_select_b")),
                filter_senders=pcf.get("sender_select"),
                filter_receivers=pcf.get("receiver_select"),
                filter_ligands=pcf.get("ligand_select"),
                filter_receptors=pcf.get("receptor_select"),
                filter_em=pcf.get("em_select"),
                filter_target_genes=pcf.get("target_select"),
                filter_all_molecules=pcf.get("any_role_select"),
                fs_bounds=slider_values.get("final-score"),
                sw_threshold=slider_values.get("sigweight"),
                rnas_bounds=slider_values.get("rna-score"),
                pval_threshold=slider_values.get("p-value"),
            )

            a_pathways = pf.filter("a", should_filter_umap=sdi.get("has_umap"))
            b_pathways = pf.filter("b", should_filter_umap=sdi.get("has_umap"))

            return (
                dcc.send_data_frame(
                    a_pathways.to_csv, f"{sdi.get('group_a_name')}.csv"
                ),
                dcc.send_data_frame(
                    b_pathways.to_csv, f"{sdi.get('group_b_name')}.csv"
                ),
            )

    @app.callback(
        Output("modal", "is_open"),
        [Input("open", "n_clicks"), Input("close", "n_clicks")],
        [State("modal", "is_open")],
    )
    def toggle_modal(n1, n2, is_open):
        if n1 or n2:
            return not is_open
        return is_open

    return app
