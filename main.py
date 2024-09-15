import logging
import pandas as pd
from dash import Dash, html

from util import *
from dtypes import pathway_dtypes

from callbacks import (
    apply_filter_callback,
    apply_sankey_callbacks,
    apply_cluster_edge_callback,
)
from components import (
    hist_container,
    pathway_filter_components,
)

logger = logging.getLogger(__name__)


def load_clusters_files(clusters_a_path, clusters_b_path):
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


def apply_callbacks(
    app, full_pathways, full_clusters, has_rna_score, has_final_score, has_p_value
):

    apply_filter_callback(
        app,
        full_pathways,
        full_clusters,
        has_rna_score=has_rna_score,
        has_final_score=has_final_score,
        has_p_value=has_p_value,
    )
    apply_sankey_callbacks(app)
    apply_cluster_edge_callback(app)

    return app


def format_full_pathways(full_pathways: pd.DataFrame) -> pd.DataFrame:

    full_pathways.columns = full_pathways.columns.str.strip().str.lower()
    full_pathways["ligand"] = full_pathways[get_cn("path")].str.split("*").str.str[0]
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
        CN.SIGWEIGHT(full_pathways, "a"),
        CN.SIGWEIGHT(full_pathways, "b"),
    ]

    if CN.PVAL(full_pathways, "a") and CN.PVAL(full_pathways, "b"):
        TO_KEEP += [
            CN.PVAL(full_pathways, "a"),
            CN.PVAL(full_pathways, "b"),
        ]

    TO_KEEP = [c for c in TO_KEEP if c in full_pathways.columns]
    return full_pathways[TO_KEEP]


def incytr_app(pathways_file, clusters_a_filepath, clusters_b_filepath):

    app = Dash(__name__, suppress_callback_exceptions=True)
    logger.info("loading pathways....")
    full_pathways: pd.DataFrame = format_full_pathways(
        pd.read_csv(pathways_file, dtype=pathway_dtypes)
    )

    has_rna_score = rna_score_available(full_pathways)
    has_final_score = final_score_available(full_pathways)
    has_p_value = p_value_available(full_pathways)

    clusters = load_clusters_files(clusters_a_filepath, clusters_b_filepath)

    app.layout = html.Div(
        [
            pathway_filter_components(
                full_pathways,
                has_final_score=has_final_score,
                has_p_value=has_p_value,
                has_rna_score=has_rna_score,
                component_style={
                    "display": "flex",
                    "flexDirection": "column",
                    "flex": "0 0 1",
                    "width": "20vw",
                },
            ),
            html.Div(
                id="data-container",
                children=[
                    hist_container(
                        {
                            "display": "flex",
                            "flexDirection": "row",
                            "width": "50%",
                            "justifyContent": "space-between",
                            "margin": "1rem",
                            "padding": "1rem",
                        },
                    ),
                    html.Div(id="figures-container", style={"display": "flex"}),
                ],
                style={"display": "flex", "flexDirection": "column", "flex": "0 0 1"},
            ),
        ],
        id="app-container",
        style={"display": "flex", "width": "95vw", "flexDirection": "row"},
    )

    return apply_callbacks(
        app, full_pathways, clusters, has_rna_score, has_final_score, has_p_value
    )


if __name__ == "__main__":

    CLUSTERS_A_FILE = "data/mc38_fake/population_a.csv"
    CLUSTERS_B_FILE = "data/mc38_fake/population_b.csv"
    PATHWAYS_FILE = "data/mc38_fake/pathways_small_short.csv"

    app = incytr_app(PATHWAYS_FILE, CLUSTERS_A_FILE, CLUSTERS_B_FILE)

    app.run(debug=True)
