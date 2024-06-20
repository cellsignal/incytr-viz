import io
import base64
from typing import Optional

import numpy as np
import pandas as pd
import pdb
import plotly.graph_objects as go
from dash import Dash, html, dcc, ctx
from dash.dependencies import Input, Output, State
from dtypes import pathway_dtypes
from util import *

from clusters_graph import cytoscape_container, load_edges, load_nodes
from sankey import prep_sankey_data, sankey_container


CLUSTERS_FILE = "data/cluster_pop.csv"
PATHWAYS_FILE = "data/Allpaths_061524.csv"

full_pathways_df = pd.read_csv(PATHWAYS_FILE, dtype=pathway_dtypes)
clusters_df = pd.read_csv(CLUSTERS_FILE)

app = Dash(__name__, suppress_callback_exceptions=True)

# pdb.set_trace()

app.layout = html.Div(
    [cytoscape_container(full_pathways_df, clusters_df)], id="app-container"
)


@app.callback(
    Output("cytoscape-tapNodeData-output", "children"),
    Input("cytoscape-clusters-graph", "tapNodeData"),
)
def displayTapNodeData(data):
    if data:
        return "Cluster: " + data["label"] + "\nSize: " + str(data["cluster_size"])


@app.callback(
    Output("app-container", "children"),
    Input("cytoscape-clusters-graph", "tapEdgeData"),
    prevent_initial_call=True,
)
def showPathways(data):
    if data:
        return sankey_container(
            full_pathways_df,
            default_senders=[data["source"]],
            default_receivers=[data["target"]],
        )


@app.callback(
    Output("sankey-graph", "figure"),
    Output("ligand-select", "value"),
    Output("receptor-select", "value"),
    Output("em-select", "value"),
    Output("target-select", "value"),
    Output("num-pathways-displayed", "children"),
    # Input("upload-data", "contents"),
    # Input("upload-data", "filename"),
    Input("sankey-graph", "clickData"),
    # Input("sankey-graph", "restyleData"),
    Input("sender-select", "value"),
    Input("receiver-select", "value"),
    Input("ligand-select", "value"),
    Input("receptor-select", "value"),
    Input("direction-select", "value"),
    Input("em-select", "value"),
    Input("target-select", "value"),
    Input("threshold-slider", "value"),
    State("num-pathways-displayed", "children"),
)
def update_sankey(
    # upload_contents,
    # upload_filename,
    click_data,
    sender_select,
    receiver_select,
    ligand_select,
    receptor_select,
    direction_select,
    em_select,
    target_select,
    threshold,
    pathways_displayed,
):

    # if upload_contents:
    #     print(upload_contents)
    #     content_type, content_string = upload_contents.split(",")
    #     decoded = base64.b64decode(content_string)
    #     try:
    #         raw = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep=None)
    #     except Exception as e:
    #         print(e)

    direction = None
    filter_senders = None
    filter_receivers = None
    filter_ligands = None
    filter_receptors = None
    filter_em = None
    filter_target_genes = None

    if direction_select:
        direction = direction_select
    if sender_select:
        filter_senders = sender_select
    if receiver_select:
        filter_receivers = receiver_select
    if ligand_select:
        filter_ligands = ligand_select
    if receptor_select:
        filter_receptors = receptor_select
    if em_select:
        filter_em = em_select
    if target_select:
        filter_target_genes = target_select

    def _update(current, new):
        return list(set(current + [new]) if isinstance(current, list) else set([new]))

    if click_data and ctx.triggered_id == "sankey-graph":

        try:
            customdata = click_data["points"][0]["customdata"]
            node_label = customdata.split("_")[0]
            node_type = customdata.split("_")[1]
            if node_type == "Ligand":
                filter_ligands = _update(filter_ligands, node_label)
            elif node_type == "Receptor":
                filter_receptors = _update(filter_receptors, node_label)
            elif node_type == "EM":
                filter_em = _update(filter_em, node_label)
            elif node_type == "Target":
                filter_target_genes = _update(filter_target_genes, node_label)
        except Exception as e:
            print(e)

    df = filter_pathways_df(
        full_pathways_df,
        filter_senders=filter_senders,
        filter_receivers=filter_receivers,
        filter_ligands=filter_ligands,
        filter_receptors=filter_receptors,
        filter_em=filter_em,
        filter_target_genes=filter_target_genes,
        threshold=threshold,
        direction=direction,
    )

    ids, labels, source, target, value, min_score, max_score = prep_sankey_data(
        sankey_df=df, always_include_target_genes=False
    )

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="fixed",
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=labels,
                    customdata=ids,
                    hovertemplate="Node %{customdata} has total value %{value}<extra></extra>",
                    color=get_node_colors(ids),
                ),
                link=dict(source=source, target=target, value=value),
            )
        ]
    )

    fig.update_layout(title_text="my graph", font_size=10)

    return (
        fig,
        filter_ligands,
        filter_receptors,
        filter_em,
        filter_target_genes,
        len(df),
    )


@app.callback(
    Output("cytoscape-clusters-graph", "elements"),
    # Input("upload-data", "contents"),
    # Input("upload-data", "filename"),
    # Input("sankey-graph", "restyleData"),
    Input("cyto-sender-select", "value"),
    Input("cyto-receiver-select", "value"),
    prevent_initial_call=True,
)
def update_cytoscape(
    cyto_sender_select,
    cyto_receiver_select,
):

    filter_senders = None
    filter_receivers = None

    if cyto_sender_select:
        filter_senders = cyto_sender_select
    if cyto_receiver_select:
        filter_receivers = cyto_receiver_select

    df = filter_pathways_df(
        full_pathways_df,
        filter_senders=filter_senders,
        filter_receivers=filter_receivers,
    )

    nodes, edges = load_nodes(clusters_df), load_edges(df)

    return nodes + edges


if __name__ == "__main__":
    app.run(debug=True)
