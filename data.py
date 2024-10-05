import pandas as pd
import numpy as np
import pdb
import matplotlib.pyplot as plt

from util import get_cn, CN, edge_width_map
from typing import Optional, Literal


def load_cell_populations(clusters_a_path, clusters_b_path):
    def _load(clusters_path):

        fields = {"type": str, "population": int}

        clusters = pd.read_csv(clusters_path, dtype=fields)

        clusters.columns = clusters.columns.str.lower().str.strip()

        if not all(c in clusters.columns for c in fields.keys()):
            raise ValueError(
                f"Invalid cell populations file: ensure the following columns are present: {fields.keys()}"
            )

        clusters = clusters[list(fields.keys())].reset_index(drop=True)
        clusters["population"] = clusters["population"].fillna(0).astype(int)

        return clusters.set_index("type")

    clusters_a = _load(clusters_a_path)
    clusters_b = _load(clusters_b_path)

    return clusters_a.join(
        clusters_b, on="type", how="outer", lsuffix="_a", rsuffix="_b"
    )


def load_pathways_input(full_pathways: pd.DataFrame) -> pd.DataFrame:

    full_pathways.columns = full_pathways.columns.str.strip().str.lower()
    full_pathways["ligand"] = full_pathways[get_cn("path")].str.split("*").str[0]
    full_pathways["receptor"] = full_pathways[get_cn("path")].str.split("*").str[1]
    full_pathways["em"] = full_pathways[get_cn("path")].str.split("*").str[2]
    full_pathways["target"] = full_pathways[get_cn("path")].str.split("*").str[3]

    TO_KEEP = [
        get_cn("path"),
        get_cn("ligand"),
        get_cn("receptor"),
        get_cn("em"),
        get_cn("target"),
        get_cn("final_score"),
        get_cn("rna_score"),
        get_cn("sender"),
        get_cn("receiver"),
        get_cn("adjlog2fc"),
        CN.SIGWEIGHT(cols=full_pathways.columns, group="a"),
        CN.SIGWEIGHT(cols=full_pathways.columns, group="b"),
    ]

    if CN.PVAL(full_pathways, "a") and CN.PVAL(full_pathways, "b"):
        TO_KEEP += [
            CN.PVAL(full_pathways, "a"),
            CN.PVAL(full_pathways, "b"),
        ]

    TO_KEEP = [c for c in TO_KEEP if c in full_pathways.columns]
    return full_pathways[TO_KEEP]


def filter_pathways(
    full_pathways: pd.DataFrame,
    group: Literal["a", "b"],
    fs_bounds: list[float],
    sw_threshold: list[float],
    pval_threshold: float,
    rnas_bounds: float,
    filter_senders: list[Optional[str]] = [],
    filter_receivers: list[Optional[str]] = [],
    filter_ligands: list[Optional[str]] = [],
    filter_receptors: list[Optional[str]] = [],
    filter_em: list[Optional[str]] = [],
    filter_target_genes: list[Optional[str]] = [],
    filter_all_molecules: list[Optional[str]] = [],
) -> pd.DataFrame:

    df: pd.DataFrame = full_pathways.copy()

    if not filter_senders:
        filter_senders = full_pathways[get_cn("sender")].unique()
    if not filter_receivers:
        filter_receivers = full_pathways[get_cn("receiver")].unique()
    if not filter_ligands:
        filter_ligands = full_pathways[get_cn("ligand")].unique()
    if not filter_receptors:
        filter_receptors = full_pathways[get_cn("receptor")].unique()
    if not filter_em:
        filter_em = full_pathways[get_cn("em")].unique()
    if not filter_target_genes:
        filter_target_genes = full_pathways[get_cn("target")].unique()

    df = df[df[CN.SIGWEIGHT(full_pathways.columns, group)] >= sw_threshold]

    if pval_threshold:
        df = df[df[CN.PVAL(full_pathways.columns, group)] <= pval_threshold]

    if fs_bounds:
        df = df[
            (df[get_cn("final_score")] >= fs_bounds[0])
            & (df[get_cn("final_score")] <= fs_bounds[1])
        ]

    if rnas_bounds:
        df = df[
            (df[get_cn("rna_score")] >= rnas_bounds[0])
            & (df[get_cn("rna_score")] <= rnas_bounds[1])
        ]

    df = df[
        df[get_cn("ligand")].isin(filter_ligands)
        & df[get_cn("receptor")].isin(filter_receptors)
        & df[get_cn("em")].isin(filter_em)
        & df[get_cn("target")].isin(filter_target_genes)
        & df[get_cn("sender")].isin(filter_senders)
        & df[get_cn("receiver")].isin(filter_receivers)
    ]

    if not filter_all_molecules:
        return df
    else:
        return df[
            df[get_cn("ligand")].isin(filter_all_molecules)
            | df[get_cn("receptor")].isin(filter_all_molecules)
            | df[get_cn("em")].isin(filter_all_molecules)
            | df[get_cn("target")].isin(filter_all_molecules)
        ]


def load_nodes(clusters: pd.DataFrame, group) -> list[dict]:
    """
    Generate cytoscape nodes from clusters file

    clusters: clusters df with expected column names:

    type
    population_a
    population_b
    rgb_colors

    Output:

    [nodes_a, nodes_b]

    Each member of the list is a list of dicts, each dict is a node

    nodes_a ~ [{"data": {...node_data}}, .....]

    """
    # TODO clean clusters
    # clusters = clean_clusters(clusters)

    def _node_size_mapping(population: int, min_pop, scaling_factor: int = 300) -> str:
        log_pop, log_min_pop = (
            np.log2(population + 1),
            np.log2(min_pop + 1),
        )

        normalized = 1 + (log_pop - log_min_pop)
        area_px = np.round(normalized * scaling_factor, 4)
        diameter_px = np.round(np.sqrt(4 * area_px / np.pi), 4)

        return str(diameter_px) + "px"

    cmap = plt.get_cmap("tab20")

    # rgba arrays, values 0-1
    plt_colors = cmap(np.linspace(0, 1, len(clusters)))
    clusters["rgb_colors"] = [[int(x * 256) for x in c[0:3]] for c in plt_colors]

    # Ignore zeros for minimum population calculation
    min_pop = np.min(pd.concat([clusters["population_a"], clusters["population_b"]]))

    # ignore clusters with population zero
    clusters["population_a"] = clusters["population_a"].replace(0, np.nan)
    clusters["population_b"] = clusters["population_b"].replace(0, np.nan)

    min_pop = np.min(pd.concat([clusters["population_a"], clusters["population_b"]]))
    clusters["diameter_a"] = clusters["population_a"].apply(
        lambda x: _node_size_mapping(x, min_pop)
    )
    clusters["diameter_b"] = clusters["population_b"].apply(
        lambda x: _node_size_mapping(x, min_pop)
    )

    def _add_node(row: pd.Series, group: str) -> dict:

        node_type = row.name
        node_population = row["population_" + group]
        node_rgb_color = row["rgb_colors"]

        if (not node_population) or (np.isnan(node_population)):
            return np.nan

        data = dict()
        data["id"] = node_type
        data["label"] = node_type
        data["cluster_size"] = node_population
        data["width"] = row["diameter_" + group]
        data["height"] = row["diameter_" + group]
        data["background_color"] = "rgb({}, {}, {})".format(*node_rgb_color)
        return {"data": data}

    if group == "a":
        return list(
            clusters.apply(
                lambda row: _add_node(row, "a"),
                axis=1,
            ).dropna()
        )
    elif group == "b":
        return list(
            clusters.apply(
                lambda row: _add_node(row, "b"),
                axis=1,
            ).dropna()
        )


def load_edges(
    nodes: list[dict],
    pathways: pd.DataFrame,
    global_max_paths: int,
):
    """add pathways from source to target"""
    edges = []

    ## filter pathways if sender/receiver not in nodes
    node_labels = pd.Series([x["data"]["label"] for x in nodes])
    pathways = pathways[
        (pathways[get_cn("sender")].isin(node_labels))
        & (pathways[get_cn("receiver")].isin(node_labels))
    ]

    ## filter pathways that are below sigweight threshold

    if len(pathways) == 0:
        return edges

    s: pd.Series = pathways.groupby([get_cn("sender"), get_cn("receiver")]).size()

    sr_pairs = s.to_dict()
    for sr, weight in sr_pairs.items():
        source_id, target_id = sr
        data = dict()
        data["id"] = source_id + target_id
        data["source"] = source_id
        data["target"] = target_id
        data["weight"] = weight
        data["label"] = str(weight)
        data["line_color"] = next(
            x["data"]["background_color"]
            for x in nodes
            if x["data"]["label"] == source_id
        )

        edges.append({"data": data})

    if edges:
        for e in edges:
            e["data"]["width"] = edge_width_map(
                abs(e["data"]["weight"]), global_max_paths
            )

    return edges


def pathways_df_to_sankey(
    sankey_df: pd.DataFrame,
    always_include_target_genes: bool = False,
) -> tuple:

    def _get_values(
        df: pd.DataFrame, source_colname: str, target_colname: str
    ) -> pd.DataFrame:
        out = (
            df.groupby(source_colname)[target_colname]
            .value_counts()
            .reset_index(name="value")
        )
        out.rename(
            columns={source_colname: "Source", target_colname: get_cn("target")},
            inplace=True,
        )
        out["source_id"] = out["Source"] + "_" + source_colname
        out["target_id"] = out[get_cn("target")] + "_" + target_colname

        return out

    l_r = _get_values(sankey_df, get_cn("ligand"), get_cn("receptor"))
    r_em = _get_values(sankey_df, get_cn("receptor"), get_cn("em"))
    em_t = _get_values(sankey_df, get_cn("em"), get_cn("target"))

    included_links = [l_r, r_em]

    ## auto-determine if target genes should be included
    def _should_display_targets() -> bool:
        num_targets = len(em_t[get_cn("target")].unique())

        return True if always_include_target_genes else num_targets <= 75

    if _should_display_targets():
        included_links.append(em_t)

    links = pd.concat(included_links, axis=0).reset_index(drop=True)
    # ids allow for repeating labels in ligand, receptor, etc. without pointing to same node
    ids = list(set(pd.concat([links["source_id"], links["target_id"]])))
    labels = [x.split("_")[0] for x in ids]
    source = [next(i for i, e in enumerate(ids) if e == x) for x in links["source_id"]]
    target = [next(i for i, e in enumerate(ids) if e == x) for x in links["target_id"]]
    value = links["value"]

    return (ids, labels, source, target, value)
