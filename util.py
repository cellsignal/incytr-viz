import pandas as pd
import numpy as np
from typing import Optional
from dash.dependencies import Input, Output, State
import json
import pdb
import re


def parse_umap_filter_data(umap_json_str):
    if umap_json_str:
        out = json.loads(umap_json_str)
        if out.get("xaxis.range[0]"):
            return out
    return None


def parse_slider_values_from_tree(children):

    sliders = []
    for c in children:
        try:
            s = c["props"]["children"][0]

            if s["type"] in ("Slider", "RangeSlider"):
                sliders.append(s)

        except KeyError:
            continue

    def _get_slider_value(sliders, index):
        return next(s for s in sliders if s["props"]["id"]["index"] == index)["props"][
            "value"
        ]

    slider_ids = ["sigweight", "rna-score", "final-score", "p-value"]

    out = {id: None for id in slider_ids}

    for id in slider_ids:
        try:
            out[id] = _get_slider_value(sliders, id)
        except StopIteration:
            continue

    return out


class PathwaysFilter:

    NAMESPACED_COLUMNS = ["sigweight", "p_value", "siks_score"]

    def __init__(
        self,
        all_paths,
        group_a_name,
        group_b_name,
        sw_threshold: float = None,
        pval_threshold: float = None,
        fs_bounds: list[float] = None,
        rnas_bounds: list[float] = None,
        filter_senders=None,
        filter_receivers=None,
        filter_ligands=None,
        filter_receptors=None,
        filter_em=None,
        filter_target_genes=None,
        filter_all_molecules=None,
        filter_umap_a: dict = None,
        filter_umap_b: dict = None,
    ):
        self.all_paths = all_paths
        self.group_a_name = group_a_name
        self.group_b_name = group_b_name
        self.a_suffix = f"_{group_a_name}"
        self.b_suffix = f"_{group_b_name}"
        self.sw_threshold = sw_threshold
        self.pval_threshold = pval_threshold
        self.fs_bounds = fs_bounds
        self.rnas_bounds = rnas_bounds
        self.filter_umap_a = filter_umap_a
        self.filter_umap_b = filter_umap_b
        self.filter_all_molecules = filter_all_molecules
        self.filter_senders = filter_senders
        self.filter_receivers = filter_receivers
        self.filter_ligands = filter_ligands
        self.filter_receptors = filter_receptors
        self.filter_em = filter_em
        self.filter_target_genes = filter_target_genes

        if not self.filter_senders:
            self.filter_senders = all_paths["sender"].unique()

        if not self.filter_receivers:
            self.filter_receivers = all_paths["receiver"].unique()

        if not self.filter_ligands:
            self.filter_ligands = all_paths["ligand"].unique()

        if not self.filter_receptors:
            self.filter_receptors = all_paths["receptor"].unique()

        if not self.filter_em:
            self.filter_em = all_paths["em"].unique()

        if not self.filter_target_genes:
            self.filter_target_genes = all_paths["target"].unique()

    def get_namespaced_columns(self):
        return [
            c
            for c in self.all_paths.columns
            if any(c.startswith(ns) for ns in self.NAMESPACED_COLUMNS)
        ]

    @property
    def a_data(self):

        df = self.all_paths.loc[
            :, ~self.all_paths.columns.str.endswith(self.b_suffix)
        ].copy()

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
            filter_umap = self.filter_umap_a

        if should_filter_umap:
            df = df[
                (df["umap1"] >= filter_umap["xaxis.range[0]"])
                & (df["umap1"] <= filter_umap["xaxis.range[1]"])
                & (df["umap2"] >= filter_umap["yaxis.range[0]"])
                & (df["umap2"] <= filter_umap["yaxis.range[1]"])
            ]

        df = df[df["sigweight"] >= self.sw_threshold]

        if self.pval_threshold:
            df = df[df["p_value"] <= self.pval_threshold]

        if self.fs_bounds:
            df = df[
                (df["final_score"] >= self.fs_bounds[0])
                & (df["final_score"] <= self.fs_bounds[1])
            ]

        if self.rnas_bounds:
            df = df[
                (df["rna_score"] >= self.rnas_bounds[0])
                & (df["rna_score"] <= self.rnas_bounds[1])
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

        return df


def update_filter_value(current, new):
    return list(set(current + [new]) if isinstance(current, list) else set([new]))


def clean_clusters(df) -> pd.DataFrame:

    # TODO handle duplicates
    # duplicates = df[df.duplicated(subset="Type", keep=False)]
    return df


def edge_width_map(pathways: int, global_max_paths: int, max_width_px: int = 10):
    floor = 2
    pixels = max((pathways / global_max_paths * max_width_px), floor)
    return str(pixels) + "px"


def get_node_colors(ids):

    colors = {
        "ligand": "red",
        "receptor": "blue",
        "em": "green",
        "target": "purple",
    }
    return [colors[x.split("_")[1]] for x in ids]
