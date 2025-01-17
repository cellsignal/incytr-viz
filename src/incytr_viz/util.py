import pandas as pd
import numpy as np
import json
import pdb
import re
import logging
from dataclasses import dataclass, field


def filter_defaults():

    return {
        "sender_select": [],
        "receiver_select": [],
        "ligand_select": [],
        "receptor_select": [],
        "em_select": [],
        "target_select": [],
        "any_role_select": [],
        "kinase_select": None,
        "sigprob": 0.9,
        "p_value": 0.05,
        "prs": [-2, 2],
        "tprs": [-2, 2],
    }


def create_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    ch.setFormatter(logging.Formatter("%(asctime)s - %(message)s"))

    logger.addHandler(ch)

    return logger


def parse_umap_filter_data(umap_json_str):
    if umap_json_str:
        out = json.loads(umap_json_str)
        if out.get("xaxis.range[0]"):
            return out
    return {}


def parse_slider_values_from_tree(children):

    sliders = []
    for c in children:
        sliders.extend(c["props"]["children"][0]["props"]["children"])
        sliders.extend(c["props"]["children"][1]["props"]["children"])

    def _get_slider_value(sliders, index):
        return next(
            s for s in sliders if s["props"].get("id", {}).get("index", "") == index
        )["props"]["value"]

    slider_ids = ["sigprob", "tprs", "prs", "p-value"]

    out = {id: None for id in slider_ids}
    for id in slider_ids:
        try:
            out[id] = _get_slider_value(sliders, id)
        except StopIteration:
            continue

    return out


@dataclass
class PathwaysFilter:

    NAMESPACED_COLUMNS = ["sigprob", "p_value", "siks_score"]

    all_paths: pd.DataFrame
    group_a_name: str
    group_b_name: str
    sp_threshold: float = 0
    pval_threshold: float = None
    prs_bounds: list[float] = field(default_factory=list)
    tprs_bounds: list[float] = field(default_factory=list)
    filter_kinase: str = ""
    filter_senders: list[str] = field(default_factory=list)
    filter_receivers: list[str] = field(default_factory=list)
    filter_ligands: list[str] = field(default_factory=list)
    filter_receptors: list[str] = field(default_factory=list)
    filter_em: list[str] = field(default_factory=list)
    filter_target_genes: list[str] = field(default_factory=list)
    filter_all_molecules: list[str] = field(default_factory=list)
    filter_umap_a: dict = field(default_factory=dict)
    filter_umap_b: dict = field(default_factory=dict)

    def __post_init__(self):

        self.a_suffix = f"_{self.group_a_name}"
        self.b_suffix = f"_{self.group_b_name}"

        if not self.filter_senders:
            self.filter_senders = self.all_paths["sender"].unique()

        if not self.filter_receivers:
            self.filter_receivers = self.all_paths["receiver"].unique()

        if not self.filter_ligands:
            self.filter_ligands = self.all_paths["ligand"].unique()

        if not self.filter_receptors:
            self.filter_receptors = self.all_paths["receptor"].unique()

        if not self.filter_em:
            self.filter_em = self.all_paths["em"].unique()

        if not self.filter_target_genes:
            self.filter_target_genes = self.all_paths["target"].unique()

    def get_namespaced_columns(self):
        return [
            c
            for c in self.all_paths.columns
            if any(c.startswith(ns) for ns in self.NAMESPACED_COLUMNS)
        ]

    @property
    def a_data(self):
        df = self.all_paths.loc[:, ~self.all_paths.columns.str.endswith(self.b_suffix)]

        pattern = re.compile(f"_{self.group_a_name}$")

        return df.rename(
            columns=lambda x: (
                re.sub(pattern, "", x) if x in self.get_namespaced_columns() else x
            )
        )

    @property
    def b_data(self):
        df = self.all_paths.loc[
            :, ~self.all_paths.columns.str.endswith(self.a_suffix)
        ].copy()

        pattern = re.compile(f"_{self.group_b_name}$")

        return df.rename(
            columns=lambda x: (
                re.sub(pattern, "", x) if x in self.get_namespaced_columns() else x
            )
        )

    def filter(self, group_id, should_filter_umap=False):
        if group_id == "a":
            df = self.a_data
            filter_umap = self.filter_umap_a
        elif group_id == "b":
            df = self.b_data
            filter_umap = self.filter_umap_b

        if should_filter_umap:
            if filter_umap.get("xaxis.range[0]"):
                df = df.loc[
                    (
                        (df["umap1"] >= filter_umap["xaxis.range[0]"])
                        & (df["umap1"] <= filter_umap["xaxis.range[1]"])
                    ),
                    :,
                ]
            if filter_umap.get("yaxis.range[0]"):
                df = df.loc[
                    (
                        (df["umap2"] >= filter_umap["yaxis.range[0]"])
                        & (df["umap2"] <= filter_umap["yaxis.range[1]"])
                    ),
                    :,
                ]

        # pdb.set_trace()
        df = df[df["sigprob"] >= self.sp_threshold]
        if self.pval_threshold:
            df = df[df["p_value"] <= self.pval_threshold]
        if self.prs_bounds:
            df = df[
                (df["prs"] >= self.prs_bounds[0]) & (df["prs"] <= self.prs_bounds[1])
            ]
        if self.tprs_bounds:
            df = df[
                (df["tprs"] >= self.tprs_bounds[0])
                & (df["tprs"] <= self.tprs_bounds[1])
            ]
        df = df[
            df["ligand"].isin(self.filter_ligands)
            & df["receptor"].isin(self.filter_receptors)
            & df["em"].isin(self.filter_em)
            & df["target"].isin(self.filter_target_genes)
            & df["sender"].isin(self.filter_senders)
            & df["receiver"].isin(self.filter_receivers)
        ]
        if self.filter_all_molecules:
            df = df[
                df["ligand"].isin(self.filter_all_molecules)
                | df["receptor"].isin(self.filter_all_molecules)
                | df["em"].isin(self.filter_all_molecules)
                | df["target"].isin(self.filter_all_molecules)
            ]
        if self.filter_kinase:
            val = self.filter_kinase
            if not all(
                x in df.columns
                for x in ["kinase_r_of_em", "kinase_r_of_t", "kinase_em_of_t"]
            ):
                df = df.iloc[0:0]
            elif val == "r_em":
                df = df[~df["kinase_r_of_em"].isna()]
            elif val == "r_t":
                df = df[~df["kinase_r_of_t"].isna()]
            elif val == "em_t":
                df = df[~df["kinase_em_of_t"].isna()]

        return df


def update_filter_value(current, new):
    return list(set(current + [new]) if isinstance(current, list) else set([new]))


def clean_clusters(df) -> pd.DataFrame:

    # TODO handle duplicates
    # duplicates = df[df.duplicated(subset="Type", keep=False)]
    return df


def edge_width_map(
    pathways: int, global_max_paths: int, edge_scale_factor, max_width_px: int = 10
):
    floor = 2
    pixels = (
        max((pathways / global_max_paths * max_width_px), floor) ** edge_scale_factor
    )
    return str(pixels) + "px"


def get_node_colors(ids):

    colors = {
        "ligand": "red",
        "receptor": "blue",
        "em": "green",
        "target": "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]


def log_base(x, base):
    """
    Calculates the logarithm of x to the given base.

    Args:
      x: The number for which to calculate the logarithm.
      base: The base of the logarithm.

    Returns:
      The logarithm of x to the base 'base'.
    """
    return np.log(x) / np.log(base)


def ascii():
    return """
  _____                      _          __      __ _      
 |_   _|                    | |         \ \    / /(_)     
   | |   _ __    ___  _   _ | |_  _ __   \ \  / /  _  ____
   | |  | '_ \  / __|| | | || __|| '__|   \ \/ /  | ||_  /
  _| |_ | | | || (__ | |_| || |_ | |       \  /   | | / / 
 |_____||_| |_| \___| \__, | \__||_|        \/    |_|/___|
                       __/ |                              
                      |___/                               
"""
