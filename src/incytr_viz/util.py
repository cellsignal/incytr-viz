import json
import logging
from importlib import resources as impresources
import pdb
import re
from flask_caching import Cache
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tabulate import tabulate
from typing import Literal, Callable

from incytr_viz import assets
from incytr_viz.dtypes import clusters_dtypes, pathways_dtypes


def get_help_file():
    helpfile = impresources.files(assets) / "help.md"

    with helpfile.open("rt") as f:
        return f.read()


cache = Cache()


def create_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    ch.setFormatter(logging.Formatter("%(asctime)s - %(message)s"))

    logger.addHandler(ch)

    return logger


logger = create_logger(__name__)


# def cache_func(
#     func: Callable, cache_method: Literal["cached", "memoize"], **cache_kwargs
# ):

#     if cache_method == "cached":
#         method = cache.cached
#     elif cache_method == "memoize":
#         method = cache.memoize

#     return method(**cache_kwargs)(func)


def format_headers(headers):
    return (
        headers.str.strip()
        .str.lower()
        .str.replace(" ", "_")
        .str.replace("sender.group", "sender")
        .str.replace("receiver.group", "receiver")
    )


def get_clusters(fpath):
    if ".csv" in fpath:
        sep = ","
    elif ".tsv" in fpath:
        sep = "\t"
    else:
        raise ValueError(
            f"Pathways file suffix must be in [.csv,.tsv] -- check filename: {fpath}"
        )

    logger.info(
        "Loading cluster populations from {} as {}".format(
            fpath, {"\t": "TSV", ",": "CSV"}[sep]
        )
    )

    df = pd.read_csv(fpath, dtype=clusters_dtypes, sep=sep, compression="infer")

    df.columns = df.columns.str.lower().str.strip()
    if not all(c in df.columns for c in clusters_dtypes.keys()):
        raise ValueError(
            f"Invalid cell populations file: ensure the following columns are present: {clusters_dtypes.keys()}"
        )

    df = df[list(clusters_dtypes.keys())].reset_index(drop=True)

    df["type_userlabel"] = df["type"]
    df["condition_userlabel"] = df["condition"]

    df["type"] = df["type"].str.strip().str.lower()
    df["group"] = df["condition"].str.strip().str.lower()
    df = df.set_index("type")

    df["population"] = df["population"].fillna(0)
    df["pop_min_ratio"] = df["population"] / (
        df[df["population"] > 0]["population"].min()
    )

    df.drop(columns=["condition"], inplace=True)
    # assign colors to each cell type
    cmap = plt.get_cmap("tab20")

    cell_types = df.index.unique()

    plt_colors = cmap(np.linspace(0, 1, len(cell_types)))

    # 256 would not be websafe value
    rgb_colors = [[int(x * 255) for x in c[0:3]] for c in plt_colors]

    colors = {t: rgb_colors[i] for i, t in enumerate(cell_types)}
    df["color"] = df.index.map(colors)
    df["color"] = df["color"].apply(lambda x: f"rgb({x[0]},{x[1]},{x[2]})")

    if len(df["group"].unique()) != 2:
        raise ValueError(
            f"Expected exactly 2 groups in cluster populations file, found {len(df['group'].unique())}"
        )
    return df, df["group"].unique()


def kinase_color_map():
    return {
        "sik_r_of_em": "red",
        "sik_em_of_t": "blue",
        "sik_em_of_r": "green",
        "sik_t_of_em": "yellow",
        "bidirectional": "black",
    }


def parse_pathway_headers(headers, group_a, group_b):

    formatted = format_headers(headers)

    # Remove duplicate columns
    mapper = list(zip(headers, formatted))

    required = [
        "path",
        "sender",
        "receiver",
        "afc",
        "sigprob_" + group_a,
        "sigprob_" + group_b,
    ]

    optional = [
        "p_value_" + group_a,
        "p_value_" + group_b,
        "tprs",
        "prs",
        "sik_r_of_em",
        "sik_r_of_t",
        "sik_em_of_t",
        "sik_em_of_r",
        "sik_t_of_r",
        "sik_t_of_em",
        "umap1",
        "umap2",
    ]

    logger.info("scanning pathways file for required and optional columns")

    required_df = pd.DataFrame.from_dict(
        {"colname": required, "required": True, "found": False}
    )
    optional_df = pd.DataFrame.from_dict(
        {"colname": optional, "required": False, "found": False}
    )
    columns_df = pd.concat([required_df, optional_df], axis=0)

    for row in columns_df.iterrows():
        col = row[1]["colname"]
        columns_df.loc[columns_df["colname"] == col, "found"] = col in formatted

    logger.info(
        "Pathways file column summary\n"
        + tabulate(columns_df, headers="keys", tablefmt="fancy_grid", showindex=False)
    )

    if (columns_df["required"] & ~columns_df["found"]).any():
        raise ValueError(
            f"Required columns not found in pathways file: {columns_df[columns_df['required'] & ~columns_df['found']]['colname'].values}"
        )

    if (~columns_df["found"] & ~columns_df["required"]).any():
        logger.warning(
            f"Optional columns missing in pathways file: {columns_df[~columns_df['found'] & ~columns_df['required']]['colname'].values}"
        )

    return [
        x[0]
        for x in mapper
        if x[1] in columns_df[columns_df["found"]]["colname"].values
    ]


class PathwayInput:

    def __init__(self, group_a, group_b, paths):
        self.group_a = group_a
        self.group_b = group_b
        self.paths = paths
        self.has_tprs = "tprs" in self.paths.columns
        self.has_prs = "prs" in self.paths.columns
        self.has_p_value = all(
            x in self.paths.columns
            for x in ["p_value_" + self.group_a, "p_value_" + self.group_b]
        )
        self.has_umap = all(x in self.paths.columns for x in ["umap1", "umap2"])
        self.unique_senders = self.paths["sender"].unique()
        self.unique_receivers = self.paths["receiver"].unique()
        self.unique_ligands = self.paths["ligand"].unique()
        self.unique_receptors = self.paths["receptor"].unique()
        self.unique_em = self.paths["em"].unique()
        self.unique_targets = self.paths["target"].unique()


def get_pathways(fpath, group_a, group_b):
    if ".csv" in fpath:
        sep = ","
    elif ".tsv" in fpath:
        sep = "\t"
    else:
        raise ValueError(
            f"Pathways file suffix must be in [.csv,.tsv] -- check filename {fpath}"
        )

    logger.info(
        "Detected pathways at path {} as {}".format(
            fpath, {"\t": "TSV", ",": "CSV"}[sep]
        )
    )

    headers = pd.read_csv(fpath, nrows=0, sep=sep).columns
    to_keep = parse_pathway_headers(headers, group_a, group_b)

    logger.info("Loading pathways............")
    paths = pd.read_csv(fpath, dtype=pathways_dtypes, usecols=to_keep, sep=sep)
    paths.columns = format_headers(paths.columns)

    paths = paths.astype(
        {k: v for k, v in pathways_dtypes.items() if k in paths.columns}
    )

    num_invalid = 0

    incomplete_paths = paths["path"].str.strip().str.split("*").str.len() != 4

    if incomplete_paths.sum() > 0:
        logger.warning(
            f"{incomplete_paths.sum()} rows with invalid pathway format found. Expecting form L*R*EM*T"
        )
        logger.warning("First 10 invalid paths:")
        logger.warning(paths[incomplete_paths]["path"].head().values)

    num_invalid += len(incomplete_paths)

    paths = paths.loc[~incomplete_paths]

    paths["ligand"] = paths["path"].str.split("*").str[0].str.strip()
    paths["receptor"] = paths["path"].str.split("*").str[1].str.strip()
    paths["em"] = paths["path"].str.split("*").str[2].str.strip()
    paths["target"] = paths["path"].str.split("*").str[3].str.strip()
    paths["sender"] = paths["sender"].str.strip().str.lower()
    paths["receiver"] = paths["receiver"].str.strip().str.lower()
    paths["path"] = (
        paths["path"]
        .str.cat(paths["sender"], sep="*")
        .str.cat(paths["receiver"], sep="*")
    )

    duplicates_mask = paths.duplicated()
    if duplicates_mask.sum() > 0:
        logger.warning(f"{duplicates_mask.sum()} duplicate rows found")

    is_na_mask = (
        paths[["afc", "sigprob_" + group_a, "sigprob_" + group_b]].isna().any(axis=1)
    )
    if is_na_mask.sum() > 0:
        logger.info(
            f"{is_na_mask.sum()} rows with invalid values found in required columns"
        )

    invalid = duplicates_mask | is_na_mask

    if invalid.sum() > 0:
        logger.info(f"Removing {invalid.sum()} duplicate or invalid rows")

    paths = paths[~invalid].reset_index(drop=True)

    kinase_cols = [
        "sik_r_of_em",
        "sik_r_of_t",
        "sik_em_of_t",
        "sik_em_of_r",
        "sik_t_of_r",
        "sik_t_of_em",
    ]
    for col in kinase_cols:
        if col in paths.columns:
            paths[col] = paths[col].replace([0, "NA", "nan", False], "")

    return PathwayInput(
        group_a=group_a,
        group_b=group_b,
        paths=paths,
    )


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
        "sigprob": 0.7,
        "p_value": 0.05,
        "prs": [-0.5, 0.5],
        "tprs": [-0.5, 0.5],
    }


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
    filter_afc_direction: bool
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
        self.cache = cache

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
        if self.filter_afc_direction:
            df = df.loc[df["afc"] > 0]
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

        if self.filter_afc_direction:
            df = df.loc[df["afc"] < 0]

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

        df = df[df["sigprob"] >= self.sp_threshold]
        if self.pval_threshold:
            df = df[df["p_value"] <= self.pval_threshold]

        if self.prs_bounds:
            df = df[
                (df["prs"] <= self.prs_bounds[0]) | (df["prs"] >= self.prs_bounds[1])
            ]
        if self.tprs_bounds:
            df = df[
                (df["tprs"] <= self.tprs_bounds[0])
                | (df["tprs"] >= self.tprs_bounds[1])
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

            try:
                if val == "r_em":
                    df = df[~(df["sik_r_of_em"] == "")]
                elif val == "r_t":
                    df = df[~(df["sik_r_of_t"] == "")]
                elif val == "em_t":
                    df = df[~(df["sik_em_of_t"] == "")]
                elif val == "em_r":
                    df = df[~(df["sik_em_of_r"] == "")]
                elif val == "t_r":
                    df = df[~(df["sik_r_of_r"] == "")]
                elif val == "t_em":
                    df = df[~(df["sik_t_of_em"] == "")]
            except KeyError:
                logger.warning(
                    f"kinase column not detected for {val} -- please check input"
                )
                df = df.iloc[0:0]

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
