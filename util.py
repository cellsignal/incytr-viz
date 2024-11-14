import pdb
import enum
from math import pi

import pandas as pd
import plotly.express as px
from dash import html, dcc
from dash.dependencies import Input


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
        if group == "a":
            return next(c for c in cols if "sigweight" in c).split("_")[1]

        elif group == "b":
            return next(c for c in cols[::-1] if "sigweight" in c).split("_")[1]

    @classmethod
    def SIGWEIGHT(cls, cols, group):
        suffix = cls.group_suffix(cols, group)
        try:
            return "sigweight_" + suffix
        except StopIteration:
            raise ValueError("No sigweight column found")

    @classmethod
    def PVAL(cls, cols, group):
        return "p_value_" + cls.group_suffix(cols, group)

    @classmethod
    def rna_score_available(cls, full_pathways):
        return "rna_score" in full_pathways.columns

    @classmethod
    def final_score_available(cls, full_pathways):
        return "final_score" in full_pathways.columns

    @classmethod
    def p_value_available(cls, full_pathways):
        return any("p_value_" in c for c in full_pathways.columns)


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
    return CN.group_suffix(full_pathways.columns, group)


def get_node_colors(ids):

    colors = {
        "ligand": "red",
        "receptor": "blue",
        "em": "green",
        "target": "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]
