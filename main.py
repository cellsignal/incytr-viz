import io
import base64
from typing import Optional

import numpy as np
import pandas as pd
import pdb
import plotly.graph_objects as go
from dash import Dash, html, dcc, ctx
from dash.dependencies import Input, Output, State
from dtypes import pathway_dtypes
from util import *

from clusters import (
    get_cytoscape_component,
    apply_cytoscape_callbacks,
)
from sankey import pathways_df_to_sankey, apply_sankey_callbacks
from filters import pathway_filter_components


def incytr_app(pathways_file, clusters_file):

    app = Dash(__name__, suppress_callback_exceptions=True)

    full_pathways: pd.DataFrame = pd.read_csv(pathways_file, dtype=pathway_dtypes)
    clusters: pd.DataFrame = pd.read_csv(clusters_file)

    app.layout = html.Div(
        [
            html.Div(
                [
                    html.Div(children=0, id="num-pathways-displayed"),
                ],
                style={
                    "justify-content": "space-between",
                },
            ),
            pathway_filter_components(full_pathways),
            html.Div(
                id="figures-container",
                children=[],
            ),
        ],
        id="app-container",
    )

    apply_sankey_callbacks(app, full_pathways)
    apply_cytoscape_callbacks(app, full_pathways, clusters)

    @app.callback(
        Output("figures-container", "children"),
        Input("view-radio", "value"),
    )
    def update_view(view):
        if view == "network":
            return get_cytoscape_component(full_pathways, clusters)

        elif view == "pathways":
            return (dcc.Graph(id="pathways-figure"),)

    return app


if __name__ == "__main__":

    CLUSTERS_FILE = "data/cluster_pop.csv"
    PATHWAYS_FILE = "data/Allpaths_061524.csv"

    app = incytr_app(PATHWAYS_FILE, CLUSTERS_FILE)

    app.run(debug=True)
