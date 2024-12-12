import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb


# pathway_dtypes = {
#     "ligand": str,
#     "receptor": str,
#     "em": str,
#     "target": str,
#     "path": str,
#     "sender.group": str,
#     "receiver.group": str,
#     "sigprob_5x": np.float64,
#     "sigprob_wt": np.float64,
#     "log2fc": np.float64,
#     "afc": np.float64,
#     "p_value_5x": np.float64,
#     "p_value_wt": np.float64,
#     "ligand_pr_log2fc": np.float64,
#     "ligand_pr_afc": np.float64,
#     "receptor_pr_log2fc": np.float64,
#     "receptor_pr_afc": np.float64,
#     "em_pr_log2fc": np.float64,
#     "em_pr_afc": np.float64,
#     "target_pr_log2fc": np.float64,
#     "target_pr_afc": np.float64,
#     "ligand_ps_log2fc": np.float64,
#     "ligand_ps_afc": np.float64,
#     "receptor_ps_log2fc": np.float64,
#     "receptor_ps_afc": np.float64,
#     "em_ps_log2fc": np.float64,
#     "em_ps_afc": np.float64,
#     "target_ps_log2fc": np.float64,
#     "target_ps_afc": np.float64,
#     "ligand_py_log2fc": np.float64,
#     "ligand_py_afc": np.float64,
#     "receptor_py_log2fc": np.float64,
#     "receptor_py_afc": np.float64,
#     "em_py_log2fc": np.float64,
#     "em_py_afc": np.float64,
#     "target_py_log2fc": np.float64,
#     "target_py_afc": np.float64,
#     "pr_up": np.int8,
#     "pr_down": np.int8,
#     "ps_up": np.int8,
#     "ps_down": np.int8,
#     "py_up": np.int8,
#     "py_down": np.int8,
#     "siks_r_of_em": str,
#     "siks_r_of_t": str,
#     "siks_em_of_t": str,
#     "siks_r_of_em_eicondition1": str,
#     "siks_r_of_em_eicondition2": str,
#     "siks_r_of_t_eicondition1": str,
#     "siks_r_of_t_eicondition2": str,
#     "siks_em_of_t_eicondition1": str,
#     "siks_em_of_t_eicondition2": str,
#     "siks_score_5x": np.float64,
#     "siks_score_wt": np.float64,
#     "tprs": np.float64,
#     "pprs": np.float64,
#     "phprs_ps": np.float64,
#     "phprs_py": np.float64,
#     "multimodel_score": np.float64,
#     "prs": np.float64,
#     "id_1": str,
#     "id_2": str,
# }

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
    "population": np.float64,
}


def rna_score_available(full_pathways):
    return "rna_score" in full_pathways.columns


def final_score_available(full_pathways):
    return "final_score" in full_pathways.columns


def p_value_available(full_pathways):
    return any("p_value_" in c for c in full_pathways.columns)


def umap_available(full_pathways):
    return ("umap1" in full_pathways.columns) and ("umap2" in full_pathways.columns)


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

        df["population"] = df["population"].fillna(0)
        df["type"] = df["type"].str.strip().str.lower()
        df["group"] = group_name
        df = df.set_index("type")

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


def format_pathways(paths) -> pd.DataFrame:

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

    return paths


def filter_duplicates(paths):
    duplicates = paths.duplicated()
    print(f"{duplicates.sum()} duplicate rows found")

    is_na = paths.isna().any(axis=1)
    print(f"{is_na.sum()} rows with invalid values found in relevant columns")

    invalid = duplicates | is_na

    print(f"Removing {invalid.sum()} duplicate or invalid rows")

    paths = paths[~invalid].reset_index(drop=True)

    return paths


def filter_pathway_data(paths):

    sigweight_cols = [c for c in paths.columns if "sigweight" in c]
    pval_cols = [c for c in paths.columns if "p_value" in c]
    rna_score_cols = [c for c in paths.columns if "rna_score" in c]
    final_score_cols = [c for c in paths.columns if "final_score" in c]

    if (not len(sigweight_cols) == 2) or (not len(pval_cols) in [0, 2]):
        raise ValueError(
            "Ambiguous input. Are there exactly 2 SigWeight/P-Value columns?"
        )

    group_a, group_b = [c.split("_")[1] for c in sigweight_cols]

    has_rna = rna_score_available(paths)
    has_final = final_score_available(paths)
    has_p_value = p_value_available(paths)
    has_umap = umap_available(paths)

    TO_KEEP = [
        "path",
        "ligand",
        "receptor",
        "em",
        "target",
        "sender",
        "receiver",
        "umap1",
        "umap2",
        "adjlog2fc",
        *sigweight_cols,
        *pval_cols,
        *rna_score_cols,
        *final_score_cols,
    ]

    paths = paths[[c for c in TO_KEEP if c in paths.columns]]

    paths = filter_duplicates(paths)

    return [paths, has_rna, has_final, has_p_value, has_umap, group_a, group_b]


def process_input_data(pathways_path):

    if ".csv" in pathways_path:
        sep = ","
    elif ".tsv" in pathways_path:
        sep = "\t"
    else:
        raise ValueError("Pathways file must be a CSV or TSV -- check filename")

    paths = pd.read_csv(
        pathways_path, dtype=pathway_dtypes, sep=sep, compression="infer"
    )

    return filter_pathway_data(format_pathways(paths))
