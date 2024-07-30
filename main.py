import pandas as pd
from dash import Dash, html


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

from util import *


def apply_callbacks(app, full_pathways, full_clusters):

    apply_filter_callback(app, full_pathways, full_clusters)
    apply_sankey_callbacks(app)
    apply_cluster_edge_callback(app)

    return app


def incytr_app(pathways_file, clusters_file):

    app = Dash(__name__, suppress_callback_exceptions=True)

    full_pathways: pd.DataFrame = pd.read_csv(pathways_file, dtype=pathway_dtypes)
    full_clusters: pd.DataFrame = pd.read_csv(clusters_file)

    app.layout = html.Div(
        [
            pathway_filter_components(full_pathways),
            html.Div(
                id="data-container",
                children=[
                    hist_container({"display": "flex", "flexDirection": "row"}),
                    html.Div(id="figures-container"),
                ],
                style={"display": "flex", "flexDirection": "column"},
            ),
        ],
        id="app-container",
        style={"display": "flex", "width": "100vw"},
    )

    return apply_callbacks(app, full_pathways, full_clusters)


if __name__ == "__main__":

    CLUSTERS_FILE = "data/cluster_pop.csv"
    PATHWAYS_FILE = "data/Allpaths_061524.csv"

    app = incytr_app(PATHWAYS_FILE, CLUSTERS_FILE)

    app.run(debug=True)
