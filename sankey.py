import numpy as np
import pandas as pd
from dash import html, dcc


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


def sankey_container(full_pathways_df, default_senders=None, default_receivers=None):

    return html.Div(
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
                        value=default_senders,
                        multi=True,
                        clearable=True,
                        placeholder="filter senders",
                        options=full_pathways_df["Sender.group"].unique(),
                    ),
                    dcc.Dropdown(
                        id="receiver-select",
                        value=default_receivers,
                        multi=True,
                        clearable=True,
                        placeholder="filter receivers",
                        options=full_pathways_df["Receiver.group"].unique(),
                    ),
                    dcc.Dropdown(
                        id="ligand-select",
                        multi=True,
                        clearable=True,
                        placeholder="filter ligands",
                        options=full_pathways_df["Ligand"].unique(),
                    ),
                    dcc.Dropdown(
                        id="receptor-select",
                        multi=True,
                        clearable=True,
                        placeholder="filter receptors",
                        options=full_pathways_df["Receptor"].unique(),
                    ),
                    dcc.Dropdown(
                        id="em-select",
                        multi=True,
                        clearable=True,
                        placeholder="filter effectors",
                        options=full_pathways_df["EM"].unique(),
                    ),
                    dcc.Dropdown(
                        id="target-select",
                        multi=True,
                        clearable=True,
                        placeholder="filter target genes",
                        options=full_pathways_df["Target"].unique(),
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
