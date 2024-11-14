import numpy as np
import pandas as pd
from util import CN
import pdb

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


def load_cell_clusters(clusters_path):

    clusters = pd.read_csv(clusters_path, dtype=clusters_dtypes)

    clusters.columns = clusters.columns.str.lower().str.strip()

    if not all(c in clusters.columns for c in clusters_dtypes.keys()):
        raise ValueError(
            f"Invalid cell populations file: ensure the following columns are present: {clusters_dtypes.keys()}"
        )

    clusters = clusters[list(clusters_dtypes.keys())].reset_index(drop=True)

    clusters["population"] = clusters["population"].fillna(0)
    clusters["type"] = clusters["type"].str.strip()

    return clusters.set_index("type")


def load_pathways(pathways_path) -> list:

    if "csv" in pathways_path:
        sep = ","
    elif "tsv" in pathways_path:
        sep = "\t"
    else:
        raise ValueError("Pathways file must be a CSV or TSV -- check filename")

    paths = pd.read_csv(
        pathways_path, dtype=pathway_dtypes, sep=sep, compression="infer"
    )

    paths.columns = paths.columns.str.strip().str.lower()
    paths["ligand"] = paths["path"].str.split("*").str[0]
    paths["receptor"] = paths["path"].str.split("*").str[1]
    paths["em"] = paths["path"].str.split("*").str[2]
    paths["target"] = paths["path"].str.split("*").str[3]

    has_rna = CN.rna_score_available(paths)
    has_final = CN.final_score_available(paths)
    has_p_value = CN.p_value_available(paths)

    TO_KEEP = [
        "path",
        "ligand",
        "receptor",
        "em",
        "target",
        "sender",
        "receiver",
        CN.SIGWEIGHT(cols=paths.columns, group="a"),
        CN.SIGWEIGHT(cols=paths.columns, group="b"),
    ]

    if has_rna:
        TO_KEEP.append("rna_score")
    if has_final:
        TO_KEEP.append("final_score")
    if has_p_value:
        TO_KEEP += [
            CN.PVAL(paths.columns, "a"),
            CN.PVAL(paths.columns, "b"),
        ]

    out = paths[TO_KEEP]

    return [out, has_rna, has_final, has_p_value]
