import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import logging
from tabulate import tabulate
import pdb
from incytr_viz.util import create_logger


logger = create_logger(__name__)

pathway_dtypes = {
    "ligand": str,
    "receptor": str,
    "em": str,
    "target": str,
    "sender": str,
    "receiver": str,
    "ligand_sclog2fc": np.float64,
    "receptor_sclog2fc": np.float64,
    "em_sclog2fc": np.float64,
    "target_sclog2fc": np.float64,
    "path": str,
    "log2fc": np.float64,
    "afc": np.float64,
    "ligand_pr_log2fc": np.float64,
    "ligand_pr_afc": np.float64,
    "receptor_pr_log2fc": np.float64,
    "receptor_pr_afc": np.float64,
    "em_pr_log2fc": np.float64,
    "em_pr_afc": np.float64,
    "target_pr_log2fc": np.float64,
    "target_pr_afc": np.float64,
    "ligand_ps_log2fc": np.float64,
    "ligand_ps_afc": np.float64,
    "receptor_ps_log2fc": np.float64,
    "receptor_ps_afc": np.float64,
    "em_ps_log2fc": np.float64,
    "em_ps_afc": np.float64,
    "target_ps_log2fc": np.float64,
    "target_ps_afc": np.float64,
    "ligand_py_log2fc": np.float64,
    "ligand_py_afc": np.float64,
    "receptor_py_log2fc": np.float64,
    "receptor_py_afc": np.float64,
    "em_py_log2fc": np.float64,
    "em_py_afc": np.float64,
    "target_py_log2fc": np.float64,
    "target_py_afc": np.float64,
    "sc_up": np.int8,
    "sc_down": np.int8,
    "pr_up": np.int8,
    "pr_down": np.int8,
    "ps_up": np.int8,
    "ps_down": np.int8,
    "py_up": np.int8,
    "py_down": np.int8,
    "tprs": np.float64,
    "pr_score": np.float64,
    "ps_score": np.float64,
    "py_score": np.float64,
    "prs": np.float64,
    "final_score_wkinase": np.float64,
    "kinase_r_of_em": str,
    "kinase_r_of_t": str,
    "kinase_em_of_t": str,
    "kinase_r_of_em_eicondition1": str,
    "kinase_r_of_em_eicondition2": str,
    "kinase_r_of_t_eicondition1": str,
    "kinase_r_of_t_eicondition2": str,
    "kinase_em_of_t_eicondition1": str,
    "kinase_em_of_t_eicondition2": str,
    "kinase_score_5x": np.float64,
    "kinase_score_wt": np.float64,
}

clusters_dtypes = {
    "type": str,
    "condition": str,
    "population": np.float64,  # fraction of total
}


class PathwayInput:
    def __init__(self, raw, group_a, group_b):
        self.raw = raw
        self.group_a = group_a
        self.group_b = group_b
        self.paths = validate_pathways(self.raw, self.group_a, self.group_b)
        self.has_tprs = "tprs" in self.paths.columns
        self.has_prs = "prs" in self.paths.columns
        self.has_p_value = all(
            x in self.paths.columns
            for x in ["p_value_" + self.group_a, "p_value_" + self.group_b]
        )
        self.has_umap = all(x in self.paths.columns for x in ["umap1", "umap2"])


def validate_pathways(raw, group_a, group_b):

    raw.columns = raw.columns.str.strip().str.lower()
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
        "kinase_r_of_em",
        "kinase_r_of_t",
        "kinase_em_of_t",
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
        columns_df.loc[columns_df["colname"] == col, "found"] = col in raw.columns

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

    paths = raw[[c for c in columns_df["colname"] if c in raw.columns]]

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

    duplicates = paths.duplicated()
    logger.warning(f"{duplicates.sum()} duplicate rows found")

    required_cols = columns_df[columns_df["required"]]["colname"].values

    is_na = paths[required_cols].isna().any(axis=1)

    logger.info(f"{is_na.sum()} rows with invalid values found in required columns")

    invalid = duplicates | is_na

    logger.info(f"Removing {invalid.sum()} duplicate or invalid rows")
    paths = paths[~invalid].reset_index(drop=True)

    return paths


def load_clusters(clusters_path) -> pd.DataFrame:

    if ".csv" in clusters_path:
        sep = ","
    elif ".tsv" in clusters_path:
        sep = "\t"
    else:
        raise ValueError(
            f"Pathways file suffix must be in [.csv,.tsv] -- check filename: {clusters_path}"
        )

    logger.info(
        "Loading cluster populations from {} as {}".format(
            clusters_path, {"\t": "TSV", ",": "CSV"}[sep]
        )
    )

    df = pd.read_csv(clusters_path, dtype=clusters_dtypes, sep=sep, compression="infer")

    df.columns = df.columns.str.lower().str.strip()
    if not all(c in df.columns for c in clusters_dtypes.keys()):
        raise ValueError(
            f"Invalid cell populations file: ensure the following columns are present: {clusters_dtypes.keys()}"
        )

    df = df[list(clusters_dtypes.keys())].reset_index(drop=True)

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


def load_pathways(pathways_path, groups):

    if ".csv" in pathways_path:
        sep = ","
    elif ".tsv" in pathways_path:
        sep = "\t"
    else:
        raise ValueError(
            f"Pathways file suffix must be in [.csv,.tsv] -- check filename {pathways_path}"
        )

    logger.info(
        "Loading pathways from {} as {}".format(
            pathways_path, {"\t": "TSV", ",": "CSV"}[sep]
        )
    )

    raw = pd.read_csv(
        pathways_path,
        dtype=pathway_dtypes,
        sep=sep,
        compression="infer",
        low_memory=False,
    )

    return PathwayInput(raw, groups[0].lower(), groups[1].lower())
