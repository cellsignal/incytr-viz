import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, html, dcc, ctx
from dash.dependencies import Input, Output, State
from util import *


def pathways_df_to_sankey(
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


def apply_sankey_callbacks(app: Dash, full_pathways: pd.DataFrame):
    @app.callback(
        Output("ligand-select", "value", allow_duplicate=True),
        Output("receptor-select", "value", allow_duplicate=True),
        Output("em-select", "value", allow_duplicate=True),
        Output("target-select", "value", allow_duplicate=True),
        Input("pathways-figure", "clickData"),
        State("ligand-select", "value"),
        State("receptor-select", "value"),
        State("em-select", "value"),
        State("target-select", "value"),
        prevent_initial_call=True,
    )
    def update_filters_click_node(
        click_data,
        ligand_select,
        receptor_select,
        em_select,
        target_select,
    ):

        def _update(current, new):
            return list(
                set(current + [new]) if isinstance(current, list) else set([new])
            )

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

    @app.callback(
        Output("pathways-figure", "figure"),
        Output("num-pathways-displayed", "children"),
        Input("sender-select", "value"),
        Input("receiver-select", "value"),
        Input("ligand-select", "value"),
        Input("receptor-select", "value"),
        Input("direction-select", "value"),
        Input("em-select", "value"),
        Input("target-select", "value"),
        Input("threshold-slider", "value"),
        prevent_initial_call=True,
    )
    def update_sankey(
        sender_select,
        receiver_select,
        ligand_select,
        receptor_select,
        direction_select,
        em_select,
        target_select,
        threshold,
    ):

        filtered_pathways = filter_pathways(
            full_pathways,
            filter_senders=sender_select,
            filter_receivers=receiver_select,
            filter_ligands=ligand_select,
            filter_receptors=receptor_select,
            filter_em=em_select,
            filter_target_genes=target_select,
            threshold=threshold,
            direction=direction_select,
        )

        ids, labels, source, target, value, min_score, max_score = (
            pathways_df_to_sankey(
                sankey_df=filtered_pathways, always_include_target_genes=False
            )
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
            len(filtered_pathways),
        )


# html.Div(
#     dcc.Upload(
#         id="upload-data",
#         children=html.Div(["Drag and Drop or ", html.A("Select Files")]),
#         multiple=False,
#     ),
# ),


# if upload_contents:
#     print(upload_contents)
#     content_type, content_string = upload_contents.split(",")
#     decoded = base64.b64decode(content_string)
#     try:
#         raw = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep=None)
#     except Exception as e:
#         print(e)
