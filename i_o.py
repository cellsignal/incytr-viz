import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
import re


class PathwayInput:
    def __init__(self, raw, group_a, group_b):
        self.raw = raw
        self.group_a = group_a
        self.group_b = group_b
        self.paths = filter_pathways(self.raw, self.group_a, self.group_b)
        self.has_rna = "rna_score" in self.paths.columns
        self.has_final = "final_score" in self.paths.columns
        self.has_p_value = all(
            x in self.paths.columns
            for x in ["p_value_" + group_a, "p_value_" + group_b]
        )
        self.has_umap = all(x in self.paths.columns for x in ["umap1", "umap2"])


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
    "adjlog2fc": np.float64,
    "ligand_pr_log2fc": np.float64,
    "ligand_pr_adjlog2fc": np.float64,
    "receptor_pr_log2fc": np.float64,
    "receptor_pr_adjlog2fc": np.float64,
    "em_pr_log2fc": np.float64,
    "em_pr_adjlog2fc": np.float64,
    "target_pr_log2fc": np.float64,
    "target_pr_adjlog2fc": np.float64,
    "ligand_ps_log2fc": np.float64,
    "ligand_ps_adjlog2fc": np.float64,
    "receptor_ps_log2fc": np.float64,
    "receptor_ps_adjlog2fc": np.float64,
    "em_ps_log2fc": np.float64,
    "em_ps_adjlog2fc": np.float64,
    "target_ps_log2fc": np.float64,
    "target_ps_adjlog2fc": np.float64,
    "ligand_py_log2fc": np.float64,
    "ligand_py_adjlog2fc": np.float64,
    "receptor_py_log2fc": np.float64,
    "receptor_py_adjlog2fc": np.float64,
    "em_py_log2fc": np.float64,
    "em_py_adjlog2fc": np.float64,
    "target_py_log2fc": np.float64,
    "target_py_adjlog2fc": np.float64,
    "sc_up": np.int8,
    "sc_down": np.int8,
    "pr_up": np.int8,
    "pr_down": np.int8,
    "ps_up": np.int8,
    "ps_down": np.int8,
    "py_up": np.int8,
    "py_down": np.int8,
    "rna_score": np.float64,
    "pr_score": np.float64,
    "ps_score": np.float64,
    "py_score": np.float64,
    "final_score": np.float64,
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
    "population": np.float64,  # fraction of total
}


def validate_pathways(raw):

    raw.columns = raw.columns.str.strip().str.lower()

    required_fixed = [
        "ligand",
        "receptor",
        "em",
        "target",
        "path",
        "sender",
        "receiver",
        "adjlog2fc",
    ]

    print("scanning input for required columns: ")
    print(" ".join(required_fixed))

    missing = [c for c in required_fixed if c not in raw.columns]

    if missing:
        raise ValueError(
            f"Could not detect one or more required columns (check input): {missing}"
        )

    print("scanning input for condition-specific columns: ")

    sigweight_cols = [c for c in raw.columns if "sigweight" in c]

    if not len(sigweight_cols) == 2:
        raise ValueError(
            "Expected exactly 2 sigweight columns in input data (check input)"
        )

    print("sigweight columns found: ")
    print(" ".join(sigweight_cols))

    print("parsing sigweight columns for condition names")
    matches = [re.match(r"sigweight_(.+)", x) for x in sigweight_cols]

    if not all(matches):
        raise ValueError(
            "Sigweight columns not in the expected format (sigweight_<groupname>)"
        )

    group_a, group_b = [m.group(1) for m in matches]

    print("condition 1: ", group_a)
    print("condition 2: ", group_b)

    return group_a, group_b
    # assert


def load_cell_clusters(*clusters_filepaths) -> pd.DataFrame:

    out = pd.DataFrame()

    for p in clusters_filepaths:

        group_name = p.split("/")[-1].split("_")[0].lower()

        df = pd.read_csv(p, dtype=clusters_dtypes)
        df.columns = df.columns.str.lower().str.strip()

        if not all(c in df.columns for c in clusters_dtypes.keys()):
            raise ValueError(
                f"Invalid cell populations file: ensure the following columns are present: {clusters_dtypes.keys()}"
            )

        df = df[list(clusters_dtypes.keys())].reset_index(drop=True)

        df["type"] = df["type"].str.strip().str.lower()
        df["group"] = group_name
        df = df.set_index("type")

        df["population"] = df["population"].fillna(0)
        df["pop_min_ratio"] = df["population"] / (
            df[df["population"] > 0]["population"].min()
        )

        out = pd.concat([out, df], axis=0)

    # assign colors to each cell type
    cmap = plt.get_cmap("tab20")

    cell_types = out.index.unique()

    plt_colors = cmap(np.linspace(0, 1, len(cell_types)))

    # 256 would not be websafe value
    rgb_colors = [[int(x * 255) for x in c[0:3]] for c in plt_colors]

    colors = {t: rgb_colors[i] for i, t in enumerate(cell_types)}
    out["color"] = out.index.map(colors)
    out["color"] = out["color"].apply(lambda x: f"rgb({x[0]},{x[1]},{x[2]})")

    return out


def filter_pathways(paths: pd.DataFrame, group_a, group_b):

    paths.columns = paths.columns.str.strip().str.lower()

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

    sigweight_cols = ["sigweight_" + group_a, "sigweight_" + group_b]
    pval_cols = ["p_value_" + group_a, "p_value_" + group_b]

    TO_KEEP = [
        "path",
        "ligand",
        "receptor",
        "em",
        "target",
        "sender",
        "receiver",
        "adjlog2fc",
        "rna_score",
        "final_score",
        "kinase_r_of_em",
        "kinase_r_of_t",
        "kinase_em_of_t",
        "umap1",
        "umap2",
        *sigweight_cols,
        *pval_cols,
    ]

    paths = paths[[c for c in TO_KEEP if c in paths.columns]]

    duplicates = paths.duplicated()
    print(f"{duplicates.sum()} duplicate rows found")

    is_na = paths[
        paths[
            [
                "path",
                "ligand",
                "receptor",
                "em",
                "target",
                "sender",
                "receiver",
                "adjlog2fc",
                *sigweight_cols,
            ]
        ].isna()
    ].any(axis=1)

    print(f"{is_na.sum()} rows with invalid values found in required columns")

    invalid = duplicates | is_na

    print(f"Removing {invalid.sum()} duplicate or invalid rows")
    paths = paths[~invalid].reset_index(drop=True)

    return paths


def process_input_data(pathways_path):

    if ".csv" in pathways_path:
        sep = ","
    elif ".tsv" in pathways_path:
        sep = "\t"
    else:
        raise ValueError("Pathways file must be a CSV or TSV -- check filename")

    raw = pd.read_csv(
        pathways_path,
        dtype=pathway_dtypes,
        sep=sep,
        compression="infer",
        low_memory=False,
    )

    group_a, group_b = validate_pathways(raw)

    pi = PathwayInput(raw, group_a, group_b)
    return pi
