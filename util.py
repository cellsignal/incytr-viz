import pandas as pd
import pdb
from typing import Optional, Literal
import plotly.express as px


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
    floor = 1
    pixels = max((pathways / global_max_paths * max_width_px), floor)
    return str(pixels) + "px"


def get_hist(df, col, num_bins):
    return px.histogram(
        df,
        x=col,
        nbins=num_bins,
        histfunc="count",
        width=300,
        height=300,
    )


def get_node_colors(ids):

    colors = {
        "Ligand": "red",
        "Receptor": "blue",
        "EM": "green",
        "Target": "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]


def filter_pathways(
    full_pathways: pd.DataFrame,
    fs_threshold: float = 0,
    sw_threshold: float = 0,
    rnas_threshold: float = 0,
    filter_senders: list[Optional[str]] = [],
    filter_receivers: list[Optional[str]] = [],
    filter_ligands: list[Optional[str]] = [],
    filter_receptors: list[Optional[str]] = [],
    filter_em: list[Optional[str]] = [],
    filter_target_genes: list[Optional[str]] = [],
    filter_all_molecules: list[Optional[str]] = [],
    always_include_target_genes: bool = False,
) -> pd.DataFrame:

    if not filter_senders:
        filter_senders = full_pathways["Sender"].unique()
    if not filter_receivers:
        filter_receivers = full_pathways["Receiver"].unique()
    if not filter_ligands:
        filter_ligands = full_pathways["Ligand"].unique()
    if not filter_receptors:
        filter_receptors = full_pathways["Receptor"].unique()
    if not filter_em:
        filter_em = full_pathways["EM"].unique()
    if not filter_target_genes:
        filter_target_genes = full_pathways["Target"].unique()

    df = full_pathways.copy()

    df = df[df["SigWeight"] >= sw_threshold]
    df = df[df["final_score"].abs() >= fs_threshold]
    df = df[df["adjlog2FC"].abs() >= rnas_threshold]

    if filter_all_molecules:
        crosstalk_df = df[
            (
                df["Ligand"].isin(filter_all_molecules)
                | df["Receptor"].isin(filter_all_molecules)
                | df["EM"].isin(filter_all_molecules)
                | df["Target"].isin(filter_all_molecules)
            )
            & df["Sender"].isin(filter_senders)
            & df["Receiver"].isin(filter_receivers)
        ]

        return crosstalk_df

    return df[
        df["Ligand"].isin(filter_ligands)
        & df["Receptor"].isin(filter_receptors)
        & df["EM"].isin(filter_em)
        & df["Target"].isin(filter_target_genes)
        & df["Sender"].isin(filter_senders)
        & df["Receiver"].isin(filter_receivers)
    ]


# def scores_to_hist_data(scores: pd.Series):
# hist_data = pd.cut(
#     full_pathways["final_score"], bins=np.linspace(-1.0, 1.0, num=21)
# ).value_counts()
