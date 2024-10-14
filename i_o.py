import pandas as pd
from dtypes import pathway_dtypes
from util import CN, get_cn
import pdb


def load_cell_clusters(clusters_path):

    fields = {"type": str, "population": int}

    clusters = pd.read_csv(clusters_path, dtype=fields)

    clusters.columns = clusters.columns.str.lower().str.strip()

    if not all(c in clusters.columns for c in fields.keys()):
        raise ValueError(
            f"Invalid cell populations file: ensure the following columns are present: {fields.keys()}"
        )

    clusters = clusters[list(fields.keys())].reset_index(drop=True)
    clusters["population"] = clusters["population"].fillna(0).astype(float)
    clusters["type"] = clusters["type"].str.strip()

    return clusters.set_index("type")


def load_pathways(pathways_path):

    if ".csv" in pathways_path:
        sep = ","
    elif ".tsv" in pathways_path:
        sep = "\t"
    else:
        raise ValueError("Pathways file must be a CSV or TSV -- check filename")

    paths = pd.read_csv(pathways_path, dtype=pathway_dtypes, sep=sep)

    paths.columns = paths.columns.str.strip().str.lower()
    paths["ligand"] = paths[get_cn("path")].str.split("*").str[0]
    paths["receptor"] = paths[get_cn("path")].str.split("*").str[1]
    paths["em"] = paths[get_cn("path")].str.split("*").str[2]
    paths["target"] = paths[get_cn("path")].str.split("*").str[3]

    has_rna = CN.rna_score_available(paths)
    has_final = CN.final_score_available(paths)
    has_p_value = CN.p_value_available(paths)

    TO_KEEP = [
        get_cn("path"),
        get_cn("ligand"),
        get_cn("receptor"),
        get_cn("em"),
        get_cn("target"),
        get_cn("sender"),
        get_cn("receiver"),
        CN.SIGWEIGHT(cols=paths.columns, group="a"),
        CN.SIGWEIGHT(cols=paths.columns, group="b"),
    ]

    if has_rna:
        TO_KEEP.append(get_cn("rna_score"))
    if has_final:
        TO_KEEP.append(get_cn("final_score"))
    if has_p_value:
        TO_KEEP += [
            CN.PVAL(paths.columns, "a"),
            CN.PVAL(paths.columns, "b"),
        ]

    out = paths[TO_KEEP]

    return [out, has_rna, has_final, has_p_value]
