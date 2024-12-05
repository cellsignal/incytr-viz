from typing import Optional

import numpy as np
import pdb
import json
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt

from dash import Dash, ALL
from dash.dependencies import Input, Output, State

from components import (
    cytoscape_container,
    sankey_container,
    hist,
    hist_container,
    umap_container,
)


from util import edge_width_map, update_filter_value


def filter_pathways(
    all_pathways: pd.DataFrame,
    group_name: str,
    fs_bounds: list[float] = None,
    sw_threshold: list[float] = None,
    pval_threshold: float = None,
    rnas_bounds: list[float] = None,
    filter_senders: list[Optional[str]] = [],
    filter_receivers: list[Optional[str]] = [],
    filter_ligands: list[Optional[str]] = [],
    filter_receptors: list[Optional[str]] = [],
    filter_em: list[Optional[str]] = [],
    filter_target_genes: list[Optional[str]] = [],
    filter_all_molecules: list[Optional[str]] = [],
    filter_umap: dict = None,
) -> pd.DataFrame:

    df: pd.DataFrame = all_pathways.copy()

    if not filter_senders:
        filter_senders = all_pathways["sender"].unique()
    if not filter_receivers:
        filter_receivers = all_pathways["receiver"].unique()
    if not filter_ligands:
        filter_ligands = all_pathways["ligand"].unique()
    if not filter_receptors:
        filter_receptors = all_pathways["receptor"].unique()
    if not filter_em:
        filter_em = all_pathways["em"].unique()
    if not filter_target_genes:
        filter_target_genes = all_pathways["target"].unique()

    if filter_umap:
        df = df[
            (df["umap1"] >= filter_umap["xaxis.range[0]"])
            & (df["umap1"] <= filter_umap["xaxis.range[1]"])
            & (df["umap2"] >= filter_umap["yaxis.range[0]"])
            & (df["umap2"] <= filter_umap["yaxis.range[1]"])
        ]

    df = df[df[f"sigweight_{group_name}"] >= sw_threshold]

    if pval_threshold:
        df = df[df[f"p_value_{group_name}"] <= pval_threshold]

    if fs_bounds:
        df = df[
            (df["final_score"] >= fs_bounds[0]) & (df["final_score"] <= fs_bounds[1])
        ]

    if rnas_bounds:
        df = df[
            (df["rna_score"] >= rnas_bounds[0]) & (df["rna_score"] <= rnas_bounds[1])
        ]

    df = df[
        df["ligand"].isin(filter_ligands)
        & df["receptor"].isin(filter_receptors)
        & df["em"].isin(filter_em)
        & df["target"].isin(filter_target_genes)
        & df["sender"].isin(filter_senders)
        & df["receiver"].isin(filter_receivers)
    ]

    if filter_all_molecules:
        df = df[
            df["ligand"].isin(filter_all_molecules)
            | df["receptor"].isin(filter_all_molecules)
            | df["em"].isin(filter_all_molecules)
            | df["target"].isin(filter_all_molecules)
        ]

    return df


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

    global_min_population = clusters["population"].min() or 1

    def _node_size_mapping(population: int, min_pop, scaling_factor: int = 2000) -> str:
        log_pop, log_min_pop = (
            np.log2(population + 1),
            np.log2(min_pop + 1),
        )

        normalized = 1 + (log_pop - log_min_pop)
        area_px = np.round(normalized * scaling_factor, 4)
        diameter_px = np.round(np.sqrt(4 * area_px / np.pi), 4)

        return str(diameter_px) + "px"

    # ignore clusters with population zero
    clusters.loc[:, "population"] = clusters["population"].replace(0, np.nan)

    clusters.loc[:, "diameter"] = clusters.loc[:, "population"].apply(
        lambda x: _node_size_mapping(x, global_min_population)
    )

    def _add_node(row: pd.Series) -> dict:

        node_type = row.name
        node_population = row["population"]
        node_rgb_color = row["color"]

        if (not node_population) or (np.isnan(node_population)):
            return np.nan

        data = dict()
        data["id"] = node_type
        data["label"] = node_type
        data["cluster_size"] = node_population
        data["width"] = row["diameter"]
        data["height"] = row["diameter"]
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
                abs(e["data"]["weight"]), global_max_paths
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

        return num_targets <= 75

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


def store_data_inputs():
    return dict(
        has_rna=Input("has-rna", "data"),
        has_final=Input("has-final", "data"),
        has_p_value=Input("has-p-value", "data"),
        has_umap=Input("has-umap", "data"),
        group_a_name=Input("group-a-name", "data"),
        group_b_name=Input("group-b-name", "data"),
    )


def pathway_component_filter_inputs():
    return dict(
        sender_select=Input("sender-select", "value"),
        receiver_select=Input("receiver-select", "value"),
        ligand_select=Input("ligand-select", "value"),
        receptor_select=Input("receptor-select", "value"),
        em_select=Input("em-select", "value"),
        target_select=Input("target-select", "value"),
        any_role_select=Input("any-role-select", "value"),
        umap_select_a=Input("umap-select-a", "value"),
        umap_select_b=Input("umap-select-b", "value"),
        sankey_color_flow=Input("sankey-color-flow-dropdown", "value"),
    )


def apply_callbacks(app: Dash, all_pathways, clusters):

    # @app.callback(
    #     Output("modal", "is_open"),
    #     [Input("open", "n_clicks"), Input("close", "n_clicks")],
    #     [State("modal", "is_open")],
    # )
    # def toggle_modal(n1, n2, is_open):
    #     if n1 or n2:
    #         return not is_open
    #     return is_open

    # @app.callback()
    # def update_histograms():
    #     pass

    @app.callback(
        Output("group-a-container", "children"),
        Output("group-b-container", "children"),
        inputs=dict(
            sdi=store_data_inputs(),
            pcf=pathway_component_filter_inputs(),
            slider_changed=Input({"type": "numerical-filter", "index": ALL}, "value"),
            sliders_container_children=State("allSlidersContainer", "children"),
            view_radio=Input("view-radio", "value"),
        ),
    )
    def update_figure_and_histogram(
        sdi,
        pcf,
        slider_changed,
        sliders_container_children,
        view_radio,
    ):

        umap_select_a = (
            json.loads(pcf.get("umap_select_a"))
            if (pcf.get("umap_select_a")) and sdi.get("has_umap")
            else None
        )

        umap_select_b = (
            json.loads(pcf.get("umap_select_b"))
            if (pcf.get("umap_select_b")) and sdi.get("has_umap")
            else None
        )

        sliders = []

        # TODO validate this matches up with store parameters
        for el in sliders_container_children:
            try:
                s = el["props"]["children"][0]["props"]["children"]

                if s["type"] in ("Slider", "RangeSlider"):
                    sliders.append(s)

            except:
                pass

        sw_threshold = next(
            s for s in sliders if s["props"]["id"]["index"] == "sigweight"
        )["props"]["value"]

        if sdi.get("has_rna") == True:
            rnas_bounds = next(
                s for s in sliders if s["props"]["id"]["index"] == "rna-score"
            )["props"]["value"]
        else:
            rnas_bounds = None

        if sdi.get("has_final") == True:
            fs_bounds = next(
                s for s in sliders if s["props"]["id"]["index"] == "final-score"
            )["props"]["value"]
        else:
            fs_bounds = None

        if sdi.get("has_p_value") == True:
            pval_threshold = next(
                s for s in sliders if s["props"]["id"]["index"] == "p-value"
            )["props"]["value"]
        else:
            pval_threshold = None

        a_pathways = filter_pathways(
            all_pathways=all_pathways,
            group_name=sdi.get("group_a_name"),
            filter_senders=pcf.get("sender_select", None),
            filter_receivers=pcf.get("receiver_select", None),
            filter_ligands=pcf.get("ligand_select", None),
            filter_receptors=pcf.get("receptor_select", None),
            filter_em=pcf.get("em_select", None),
            filter_target_genes=pcf.get("target_select", None),
            filter_all_molecules=pcf.get("any_role_select", None),
            filter_umap=umap_select_a,
            fs_bounds=fs_bounds,
            sw_threshold=sw_threshold,
            rnas_bounds=rnas_bounds,
            pval_threshold=pval_threshold,
        )

        b_pathways = filter_pathways(
            all_pathways=all_pathways,
            group_name=sdi.get("group_b_name"),
            filter_senders=pcf.get("sender_select", None),
            filter_receivers=pcf.get("receiver_select", None),
            filter_ligands=pcf.get("ligand_select", None),
            filter_receptors=pcf.get("receptor_select", None),
            filter_em=pcf.get("em_select", None),
            filter_target_genes=pcf.get("target_select", None),
            filter_all_molecules=pcf.get("any_role_select", None),
            filter_umap=umap_select_b,
            fs_bounds=fs_bounds,
            sw_threshold=sw_threshold,
            rnas_bounds=rnas_bounds,
            pval_threshold=pval_threshold,
        )

        def _get_group_figures(
            filtered_group_paths: pd.DataFrame,
            clusters: pd.DataFrame,
            group_name: str,
            group_id: str,
            global_max_paths: int,
        ):

            if view_radio == "network":

                nodes = load_nodes(clusters[clusters["group"] == group_name])

                edges = load_edges(nodes, filtered_group_paths, global_max_paths)

                cytoscape = cytoscape_container(
                    f"cytoscape-{group_id}", group_name, nodes + edges
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
                    ids,
                    labels,
                    source,
                    target,
                    value,
                    color,
                    "Sankey " + group_name,
                    group_id,
                    warn,
                )

                graph_container = sankey

            sw_hist = hist(
                filtered_group_paths,
                f"sigweight_{group_name}",
                "sigweight",
            )
            rnas_hist = hist(filtered_group_paths, "rna_score", "rna_score")
            fs_hist = hist(filtered_group_paths, "final_score", "final_score")
            pval_hist = hist(filtered_group_paths, f"p_value_{group_id}", "p_val")

            return [
                graph_container,
                hist_container(
                    group_id,
                    sw_hist,
                    pval_hist,
                    rnas_hist,
                    fs_hist,
                ),
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
        )

    # @app.callback(
    #     Output("sender-select", "value"),
    #     Output("receiver-select", "value"),
    #     Output("view-radio", "value"),
    #     Input("cytoscape-a", "tapEdgeData"),
    #     Input("cytoscape-b", "tapEdgeData"),
    #     State("sender-select", "value"),
    #     State("receiver-select", "value"),
    #     State("view-radio", "value"),
    #     prevent_initial_call=True,
    # )
    # def cluster_edge_callback(
    #     cs_down_data, cs_up_data, sender_select, receiver_select, view_radio
    # ):
    #     data = cs_down_data or cs_up_data
    #     if data:
    #         return (
    #             update_filter_value([], data["source"]),
    #             update_filter_value([], data["target"]),
    #             "sankey",
    #         )
    #     else:
    #         return sender_select, receiver_select, view_radio

    def _umap_callback(relayoutData):
        # {'xaxis.range[0]': -0.6369249007630504, 'xaxis.range[1]': 6.965720316453904, 'yaxis.range[0]': 3.7282259393124537, 'yaxis.range[1]': 9.59742380103187}
        if relayoutData and "xaxis.range[0]" in relayoutData:
            return json.dumps(relayoutData)
        else:
            return None

    @app.callback(
        Output("umap-select-a", "value"),
        inputs=Input("scatter-plot-a", "relayoutData"),
        prevent_initial_call=True,
    )
    def umap_callback(
        relayoutData,
    ):
        return _umap_callback(relayoutData)

    @app.callback(
        Output("umap-select-b", "value"),
        inputs=Input("scatter-plot-b", "relayoutData"),
        prevent_initial_call=True,
    )
    def umap_callback(
        relayoutData,
    ):
        return _umap_callback(relayoutData)

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

    return app

    # # @app.callback(*outputs, *inputs)
    # # def display_tooltip(node):
    # #     if node:

    # #         id = node["data"]["id"]
    # #         title = node["data"]["label"]
    # #         size = node["data"]["cluster_size"]
    # #         sum_outward = sum(
    # #             [e["weight"] for e in node["edgesData"] if e["source"] == id]
    # #         )
    # #         sum_inward = sum(
    # #             [e["weight"] for e in node["edgesData"] if e["target"] == id]
    # #         )

    # #         return [
    # #             html.H4(f"{title}"),
    # #             html.P(f"Node Size: {size}"),
    # #             html.P(f"Sum of Outward Edges: {sum_outward}"),
    # #             html.P(f"Sum of Inward Edges: {sum_inward}"),
    # #         ]
