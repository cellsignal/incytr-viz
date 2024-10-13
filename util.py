import pdb
import enum
from math import pi

import pandas as pd
import plotly.express as px
from dash import html, dcc
from dash.dependencies import Input


def pathway_component_filter_inputs():
    return dict(
        sender_select=Input("sender-select", "value"),
        receiver_select=Input("receiver-select", "value"),
        ligand_select=Input("ligand-select", "value"),
        receptor_select=Input("receptor-select", "value"),
        em_select=Input("em-select", "value"),
        target_select=Input("target-select", "value"),
        all_mols_select=Input("all-molecules-select", "value"),
    )


def pathway_value_filter_inputs(
    has_rna_score: bool, has_final_score: bool, has_p_value: bool
):
    pvf = dict(
        sw_threshold=Input("sw-slider", "value"),
    )

    if has_p_value:
        pvf["pval_threshold"] = Input("pval-slider", "value")
    if has_final_score:
        pvf["fs_bounds"] = Input("fs-slider", "value")
    if has_rna_score:
        pvf["rnas_bounds"] = Input("rnas-slider", "value")

    return pvf


def view_radio_input():
    return Input("view-radio", "value")


class CN(enum.Enum):
    PATH = "path"
    LIGAND = "ligand"
    RECEPTOR = "receptor"
    EM = "em"
    TARGET = "target"
    FINAL_SCORE = "final_score"
    RNA_SCORE = "rna_score"
    SENDER = "sender"
    RECEIVER = "receiver"
    ADJLOG2FC = "adjlog2fc"

    @classmethod
    def group_suffix(cls, cols, group):
        return cls.SIGWEIGHT(cols, group).split("_")[1]

    @classmethod
    def SIGWEIGHT(cls, cols, group):
        try:
            if group == "a":
                return next(c for c in cols if "sigweight" in c)
            elif group == "b":
                return next(c for c in cols[::-1] if "sigweight" in c)
        except StopIteration:
            raise ValueError("No sigweight column found")

    @classmethod
    def PVAL(cls, cols, group):
        return "pval_" + cls.group_suffix(cols, group)

    @classmethod
    def rna_score_available(cls, full_pathways):
        return get_cn("rna_score") in full_pathways.columns

    @classmethod
    def final_score_available(cls, full_pathways):
        return get_cn("final_score") in full_pathways.columns

    @classmethod
    def p_value_available(cls, full_pathways):
        return any("pval" in c for c in full_pathways.columns)


def get_cn(enum_choice: str, **kwargs):
    if not kwargs:
        return CN[enum_choice.upper()].value
    return eval(f"CN.{enum_choice.upper()}")(**kwargs)


def update_filter_value(current, new):
    return list(set(current + [new]) if isinstance(current, list) else set([new]))


def clean_clusters(df) -> pd.DataFrame:

    # TODO handle duplicates
    # duplicates = df[df.duplicated(subset="Type", keep=False)]
    return df


def edge_width_map(pathways: int, global_max_paths: int, max_width_px: int = 10):
    floor = 0.5
    pixels = max((pathways / global_max_paths * max_width_px), floor)
    return str(pixels) + "px"


def get_hist(df, col, title, num_bins=20):
    try:
        return html.Div(
            dcc.Graph(
                figure=px.histogram(
                    df,
                    x=col,
                    title=title,
                    nbins=num_bins,
                    histfunc="count",
                    width=300,
                    height=300,
                ),
            ),
            className="hist",
        )
    except:
        return html.Div(
            dcc.Graph(
                figure=px.histogram(
                    pd.DataFrame({title: []}), title=title, width=300, height=300
                )
            ),
            className="hist",
        )


def get_group_name(full_pathways, group):
    return CN.SIGWEIGHT(full_pathways.columns, group).split("_")[1]


def get_node_colors(ids):

    colors = {
        get_cn("ligand"): "red",
        get_cn("receptor"): "blue",
        get_cn("em"): "green",
        get_cn("target"): "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]
