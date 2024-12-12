import pandas as pd
import json
import pdb
import re
from dataclasses import dataclass, field


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


@dataclass
class PathwaysFilter:

    NAMESPACED_COLUMNS = ["sigweight", "p_value", "siks_score"]

    all_paths: pd.DataFrame
    group_a_name: str
    group_b_name: str
    sw_threshold: float = 0
    pval_threshold: float = 1
    fs_bounds: list[float] = field(default_factory=list)
    rnas_bounds: list[float] = field(default_factory=list)
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
