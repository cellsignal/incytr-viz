import os
import pytest
from incytr_viz.util import PathwaysFilter, IncytrInput
from incytr_viz.app import load_edges, load_nodes, create_app

from incytr_viz.__main__ import run_wsgi


@pytest.fixture
def clusters():
    return "./tests/clusters.csv"


@pytest.fixture
def pathways():
    return "./tests/pathways.csv"


@pytest.fixture
def incytr_input(clusters, pathways):
    return IncytrInput(clusters_path=clusters, pathways_path=pathways)


@pytest.fixture
def pcf():
    return {
        "sender_select": ["cancer cells"],
        "receiver_select": None,
        "ligand_select": [],
        "receptor_select": None,
        "em_select": None,
        "target_select": None,
        "any_role_select": None,
        "sankey_color_flow": None,
        "umap_select_a": None,
        "umap_select_b": None,
        "kinase_select": None,
    }


@pytest.fixture
def slider_values():
    return {"sigprob": 0.7, "tpds": [-2, 2], "ppds": [-2, 2], "p-value": 1}


def test_get_pathways(incytr_input):

    paths = incytr_input.paths

    assert all(
        [
            x in paths.columns
            for x in [
                "path",
                "sender",
                "receiver",
                "sigprob_5x",
                "sigprob_wt",
                "afc",
                "tpds",
                "ppds",
            ]
        ]
    )

    assert incytr_input.has_p_value == True
    assert incytr_input.has_ppds == True
    assert incytr_input.has_tpds == True
    assert incytr_input.has_umap == False
    assert incytr_input.group_a == "5x"
    assert incytr_input.group_b == "wt"


def test_filter_pathways(incytr_input):

    paths = incytr_input.paths

    pf = PathwaysFilter(
        all_paths=paths,
        group_a_name=incytr_input.group_a,
        group_b_name=incytr_input.group_b,
        sp_threshold=0.7,
        filter_afc_direction=True,
    )

    filtered_a = pf.filter("a", should_filter_umap=False)
    filtered_b = pf.filter("b", should_filter_umap=False)

    assert (
        filtered_a.shape[0]
        == paths[(paths.sigprob_5x >= 0.7) & (paths.afc > 0)].shape[0]
    )
    assert (
        filtered_b.shape[0]
        == paths[(paths.sigprob_wt >= 0.7) & (paths.afc < 0)].shape[0]
    )

    assert filtered_a.shape[0] != filtered_b.shape[0] != len(paths)


def test_serve_app(clusters, pathways):
    run_wsgi(pathways=pathways, clusters=clusters)


# def test_filter_umap():
#     pass


# def test_remove_afc_filter():
#     pass


# def test_callback_outputs():
#     pass


# def test_cytoscape_nodes_edges(incytr_input):
#     clusters = incytr_input.clusters
#     paths = incytr_input.paths

#     pf = PathwaysFilter(
#         all_paths=paths,
#         group_a_name=incytr_input.group_a,
#         group_b_name=incytr_input.group_b,
#         sp_threshold=0.7,
#         filter_afc_direction=True,
#     )

#     filtered_a = pf.filter("a", should_filter_umap=False)
#     filtered_b = pf.filter("b", should_filter_umap=False)

#     nodes_a = load_nodes(clusters, node_scale_factor=2)

#     assert (
#         filtered_a.shape[0]
#         == paths[(paths.sigprob_5x >= 0.7) & (paths.afc > 0)].shape[0]
#     )
#     assert (
#         filtered_b.shape[0]
#         == paths[(paths.sigprob_wt >= 0.7) & (paths.afc < 0)].shape[0]
#     )
