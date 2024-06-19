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

#### RIVER INPUTS #####

RIVER_INPUT_FILE = "/home/icossentino/code/incytr-viz/data/Allpaths_061524.csv"
raw = pd.read_csv(RIVER_INPUT_FILE, dtype=pathway_dtypes).dropna(subset=["final_score"])


def get_node_colors(ids):

    colors = {
        "Ligand": "red",
        "Receptor": "blue",
        "EM": "green",
        "Target": "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]


def filter_pathways_df(
    raw: pd.DataFrame,
    threshold: Optional[float] = 0.95,
    direction: Optional[str] = None,
    filter_senders: Optional[list[str]] = None,
    filter_receivers: Optional[list[str]] = None,
    filter_ligands: Optional[list[str]] = None,
    filter_receptors: Optional[list[str]] = None,
    filter_em: Optional[list[str]] = None,
    filter_target_genes: Optional[list[str]] = None,
    always_include_target_genes: bool = False,
) -> pd.DataFrame:

    df = raw.copy()
    if direction == "down":
        df = df[df["final_score"] <= 0]
    if direction == "up":
        df = df[df["final_score"] >= 0]
    if threshold:
        df = df[df["final_score"].abs() >= threshold]
    if filter_senders:
        df = df[df["Sender.group"].isin(filter_senders)]
    if filter_receivers:
        df = df[df["Receiver.group"].isin(filter_receivers)]
    if filter_ligands:
        df = df[df["Ligand"].isin(filter_ligands)]
    if filter_receptors:
        df = df[df["Receptor"].isin(filter_receptors)]
    if filter_em:
        df = df[df["EM"].isin(filter_em)]
    if filter_target_genes:
        df = df[df["Target"].isin(filter_target_genes)]
    if direction:
        if direction == "positive":
            df = df[df["final_score"]] > 0
        elif direction == "negative":
            df = df[df["final_score"]] < 0

    return df


def prep_sankey_data(
    sankey_df: pd.DataFrame,
    always_include_target_genes: bool = False,
) -> tuple:

    min_score, max_score = (
        sankey_df["final_score"].min(),
        sankey_df["final_score"].max(),
    )
    if np.isnan(min_score) and np.isnan(max_score):
        min_score, max_score = 0, 0

    def _get_values(
        df: pd.DataFrame, source_colname: str, target_colname: str
    ) -> pd.DataFrame:
        out = (
            df.groupby(source_colname)[target_colname]
            .value_counts()
            .reset_index(name="value")
        )
        out.rename(
            columns={source_colname: "Source", target_colname: "Target"}, inplace=True
        )
        out["source_id"] = out["Source"] + "_" + source_colname
        out["target_id"] = out["Target"] + "_" + target_colname

        return out

    l_r = _get_values(sankey_df, "Ligand", "Receptor")
    r_em = _get_values(sankey_df, "Receptor", "EM")
    em_t = _get_values(sankey_df, "EM", "Target")

    included_links = [l_r, r_em]

    ## auto-determine if target genes should be included
    def _should_display_targets() -> bool:
        num_targets = len(em_t["Target"].unique())

        return True if always_include_target_genes else num_targets <= 50

    if _should_display_targets():
        included_links.append(em_t)

    links = pd.concat(included_links, axis=0).reset_index(drop=True)
    # ids allow for repeating labels in ligand, receptor, etc. without pointing to same node
    ids = list(set(pd.concat([links["source_id"], links["target_id"]])))
    labels = [x.split("_")[0] for x in ids]
    source = [next(i for i, e in enumerate(ids) if e == x) for x in links["source_id"]]
    target = [next(i for i, e in enumerate(ids) if e == x) for x in links["target_id"]]
    value = links["value"]

    return (ids, labels, source, target, value, min_score, max_score)


app = Dash(__name__)


app.layout = html.Div(
    [
        # html.Div(
        #     dcc.Upload(
        #         id="upload-data",
        #         children=html.Div(["Drag and Drop or ", html.A("Select Files")]),
        #         multiple=False,
        #     ),
        # ),
        html.Div(
            [
                html.Div(children=0, id="num-pathways-displayed"),
                dcc.Dropdown(
                    id="sender-select",
                    multi=True,
                    clearable=True,
                    placeholder="filter senders",
                    options=raw["Sender.group"].unique(),
                ),
                dcc.Dropdown(
                    id="receiver-select",
                    multi=True,
                    clearable=True,
                    placeholder="filter receivers",
                    options=raw["Receiver.group"].unique(),
                ),
                dcc.Dropdown(
                    id="ligand-select",
                    multi=True,
                    clearable=True,
                    placeholder="filter ligands",
                    options=raw["Ligand"].unique(),
                ),
                dcc.Dropdown(
                    id="receptor-select",
                    multi=True,
                    clearable=True,
                    placeholder="filter receptors",
                    options=raw["Receptor"].unique(),
                ),
                dcc.Dropdown(
                    id="em-select",
                    multi=True,
                    clearable=True,
                    placeholder="filter effectors",
                    options=raw["EM"].unique(),
                ),
                dcc.Dropdown(
                    id="target-select",
                    multi=True,
                    clearable=True,
                    placeholder="filter target genes",
                    options=raw["Target"].unique(),
                ),
                html.Div(
                    [
                        html.H3("Up/Down Regulated"),
                        dcc.RadioItems(
                            id="direction-select",
                            options=[
                                {"label": "Up", "value": "up"},
                                {"label": "Down", "value": "down"},
                                {"label": "All", "value": "all"},
                            ],
                            value="all",  # Set default value to 'all'
                        ),
                    ]
                ),
                html.Div(
                    [
                        dcc.Slider(
                            min=0,
                            max=1,
                            step=0.01,
                            value=0.95,
                            marks=None,
                            tooltip={"placement": "bottom", "always_visible": True},
                            id="threshold-slider",
                        )
                    ]
                ),
            ],
            style={
                "justify-content": "space-between",
            },
        ),
        dcc.Graph(id="sankey-graph"),
        html.Div(id="file-path"),
    ]
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
        raw,
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

    fig.update_layout(title_text=RIVER_INPUT_FILE, font_size=10)

    return (
        fig,
        filter_ligands,
        filter_receptors,
        filter_em,
        filter_target_genes,
        len(df),
    )


if __name__ == "__main__":

    app.run(debug=True)
