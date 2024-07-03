import pandas as pd
from typing import Optional


def get_node_colors(ids):

    colors = {
        "Ligand": "red",
        "Receptor": "blue",
        "EM": "green",
        "Target": "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]


def filter_pathways(
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
