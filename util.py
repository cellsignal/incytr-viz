import pandas as pd
import pdb
from typing import Optional, Literal
from dash import html
import plotly.express as px
import enum
import typing
from dash.dependencies import Input, Output, State


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
    def SIGWEIGHT(self, full_pathways, group):
        try:
            if group == "a":
                return next(c for c in full_pathways.columns if "sigweight" in c)
            elif group == "b":
                return next(c for c in full_pathways.columns[::-1] if "sigweight" in c)
            return ValueError("Invalid group provided", group)
        except StopIteration:
            return None

    @classmethod
    def PVAL(self, full_pathways, group):
        try:
            if group == "a":
                return next(c for c in full_pathways.columns if "p_value" in c)
            elif group == "b":
                return next(c for c in full_pathways.columns[::-1] if "p_value" in c)
            return ValueError("Invalid group provided", group)
        except StopIteration:
            return None

    @classmethod
    def group_suffix(self, full_pathways_group):
        return CN.SIGWEIGHT(full_pathways_group, "a").split("_")[1]


def get_cn(enum_choice: str):
    return CN[enum_choice.upper()].value


def rna_score_available(full_pathways):
    return get_cn("rna_score") in full_pathways.columns


def final_score_available(full_pathways):
    return get_cn("final_score") in full_pathways.columns


def p_value_available(full_pathways):
    try:
        has_p_value = CN.PVAL(full_pathways, "a") and CN.PVAL(full_pathways, "b")
    except StopIteration:
        has_p_value = False
    return has_p_value


def update_filter_value(current, new):
    return list(set(current + [new]) if isinstance(current, list) else set([new]))


def clean_clusters(df) -> pd.DataFrame:

    # TODO handle duplicates
    # duplicates = df[df.duplicated(subset="Type", keep=False)]
    return df


def node_size_map(cluster_count: int, total_count: int):
    """map cell count to node diameter. diameters will add to 50% of graph width"""
    proportion = (cluster_count / total_count) / 0.4
    return str(proportion * 100) + "%"


def edge_width_map(pathways: int, global_max_paths: int, max_width_px: int = 10):
    floor = 0.5
    pixels = max((pathways / global_max_paths * max_width_px), floor)
    return str(pixels) + "px"


def get_hist(df, col, title, num_bins=20):
    try:
        return px.histogram(
            df,
            x=col,
            title=title,
            nbins=num_bins,
            histfunc="count",
            width=300,
            height=300,
        )
    except:
        return px.histogram(
            pd.DataFrame({title: []}), title=title, width=300, height=300
        )


def get_group_name(full_pathways, group):
    return CN.SIGWEIGHT(full_pathways, group).split("_")[1]


def get_node_colors(ids):

    colors = {
        get_cn("ligand"): "red",
        get_cn("receptor"): "blue",
        get_cn("em"): "green",
        get_cn("target"): "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]


def filter_pathways(
    full_pathways: pd.DataFrame,
    group: typing.Literal["a", "b"],
    fs_bounds: list[float],
    sw_threshold: list[float],
    pval_threshold: float,
    rnas_bounds: float,
    filter_senders: list[Optional[str]] = [],
    filter_receivers: list[Optional[str]] = [],
    filter_ligands: list[Optional[str]] = [],
    filter_receptors: list[Optional[str]] = [],
    filter_em: list[Optional[str]] = [],
    filter_target_genes: list[Optional[str]] = [],
    filter_all_molecules: list[Optional[str]] = [],
    always_include_target_genes: bool = False,
) -> pd.DataFrame:

    df: pd.DataFrame = full_pathways.copy()

    if not filter_senders:
        filter_senders = full_pathways[get_cn("sender")].unique()
    if not filter_receivers:
        filter_receivers = full_pathways[get_cn("receiver")].unique()
    if not filter_ligands:
        filter_ligands = full_pathways[get_cn("ligand")].unique()
    if not filter_receptors:
        filter_receptors = full_pathways[get_cn("receptor")].unique()
    if not filter_em:
        filter_em = full_pathways[get_cn("em")].unique()
    if not filter_target_genes:
        filter_target_genes = full_pathways[get_cn("target")].unique()

    df = df[df[CN.SIGWEIGHT(full_pathways, group)] >= sw_threshold]

    if pval_threshold:
        df = df[df[CN.PVAL(full_pathways, group)] <= pval_threshold]

    if fs_bounds:
        df = df[
            (df[get_cn("final_score")] >= fs_bounds[0])
            & (df[get_cn("final_score")] <= fs_bounds[1])
        ]

    if rnas_bounds:
        df = df[
            (df[get_cn("rna_score")] >= rnas_bounds[0])
            & (df[get_cn("rna_score")] <= rnas_bounds[1])
        ]

    df = df[
        df[get_cn("ligand")].isin(filter_ligands)
        & df[get_cn("receptor")].isin(filter_receptors)
        & df[get_cn("em")].isin(filter_em)
        & df[get_cn("target")].isin(filter_target_genes)
        & df[get_cn("sender")].isin(filter_senders)
        & df[get_cn("receiver")].isin(filter_receivers)
    ]

    if not filter_all_molecules:
        return df
    else:
        crosstalk_df = df[
            df[get_cn("ligand")].isin(filter_all_molecules)
            | df[get_cn("receptor")].isin(filter_all_molecules)
            | df[get_cn("em")].isin(filter_all_molecules)
            | df[get_cn("target")].isin(filter_all_molecules)
        ]

        return crosstalk_df
