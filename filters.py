import pandas as pd
from dash import html, dcc


def pathway_filter_components(pathways: pd.DataFrame):

    return html.Div(
        [
            dcc.RadioItems(
                [
                    {
                        "label": html.Div(
                            ["Network View"], style={"color": "Gold", "font-size": 20}
                        ),
                        "value": "network",
                    },
                    {
                        "label": html.Div(
                            ["Pathways View"],
                            style={"color": "MediumTurqoise", "font-size": 20},
                        ),
                        "value": "pathways",
                    },
                ],
                value="network",
                id="view-radio",
            ),
            dcc.RadioItems(
                [
                    {
                        "label": html.Div(
                            ["Differential View"],
                            style={"color": "Gold", "font-size": 20},
                        ),
                        "value": True,
                    },
                    {
                        "label": html.Div(
                            ["Aggregate View"],
                            style={"color": "MediumTurqoise", "font-size": 20},
                        ),
                        "value": False,
                    },
                ],
                value=False,
                id="differential-radio",
            ),
            dcc.Dropdown(
                id="sender-select",
                multi=True,
                clearable=True,
                placeholder="filter senders",
                options=pathways["Sender.group"].unique(),
            ),
            dcc.Dropdown(
                id="receiver-select",
                multi=True,
                clearable=True,
                placeholder="filter receivers",
                options=pathways["Receiver.group"].unique(),
            ),
            dcc.Dropdown(
                id="ligand-select",
                multi=True,
                clearable=True,
                placeholder="filter ligands",
                options=pathways["Ligand"].unique(),
            ),
            dcc.Dropdown(
                id="receptor-select",
                multi=True,
                clearable=True,
                placeholder="filter receptors",
                options=pathways["Receptor"].unique(),
            ),
            dcc.Dropdown(
                id="em-select",
                multi=True,
                clearable=True,
                placeholder="filter effectors",
                options=pathways["EM"].unique(),
            ),
            dcc.Dropdown(
                id="target-select",
                multi=True,
                clearable=True,
                placeholder="filter target genes",
                options=pathways["Target"].unique(),
            ),
            html.Div(
                [
                    html.H3("Up/Down Regulated"),
                    dcc.RadioItems(
                        id="direction-select",
                        options=[
                            {"label": "Up", "value": "up"},
                            {"label": "Down", "value": "down"},
                            {"label": "All", "value": "all"},
                        ],
                        value="all",  # Set default value to 'all'
                    ),
                ]
            ),
            html.Div(
                [
                    dcc.Slider(
                        min=0,
                        max=1,
                        step=0.01,
                        value=0.95,
                        marks=None,
                        tooltip={"placement": "bottom", "always_visible": True},
                        id="threshold-slider",
                    )
                ]
            ),
        ],
        id="pathway-filters-container",
    )
