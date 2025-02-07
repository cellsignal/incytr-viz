import os

import pandas as pd
import pytest

import incytr_viz.dtypes
from incytr_viz.app import create_app, load_edges, load_nodes
from incytr_viz.util import IncytrInput, PathwaysFilter


@pytest.fixture
def datadir():
    return "./local/tests/demo"


@pytest.fixture
def clusters(datadir):
    return os.path.join(datadir, "clusters.csv")


@pytest.fixture
def pathways(datadir):
    return os.path.join(datadir, "pathways.csv")


@pytest.fixture
def no_p_value(clusters):
    return IncytrInput(
        clusters_path=clusters, pathways_path="./local/tests/demo/no_p_value.csv"
    )


@pytest.fixture
def formatted_pathways():
    return pd.read_csv(
        "./local/tests/demo/formatted_paths.csv",
        dtype=incytr_viz.dtypes.pathways_dtypes,
    )


@pytest.fixture
def incytr_input(clusters, pathways):
    return IncytrInput(clusters_path=clusters, pathways_path=pathways)


@pytest.fixture
def base_pathway_filter(incytr_input):
    return PathwaysFilter(
        all_paths=incytr_input.paths,
        group_a_name=incytr_input.group_a,
        group_b_name=incytr_input.group_b,
        filter_afc_direction=True,
    )


@pytest.fixture
def a_paths(base_pathway_filter):
    return base_pathway_filter.a_data


@pytest.fixture
def b_paths(base_pathway_filter):
    return base_pathway_filter.b_data


@pytest.fixture
def slider_values():
    return {"sigprob": 0.7, "tpds": [-2, 2], "ppds": [-2, 2], "p-value": 1}


def test_create_app(clusters, pathways):
    create_app(clusters_file=clusters, pathways_file=pathways)


def test_incytr_input(incytr_input, formatted_pathways):
    assert len(incytr_input.unique_senders) == len(
        formatted_pathways["sender"].unique()
    )
    assert len(incytr_input.unique_ligands) == len(
        formatted_pathways["ligand"].unique()
    )


def test_get_pathways(incytr_input):

    expected_columns = [
        "path",
        "sender",
        "receiver",
        "sigprob_5x",
        "sigprob_wt",
        "afc",
        "p_value_5x",
        "p_value_wt",
        "sik_r_of_em",
        "sik_r_of_t",
        "sik_em_of_t",
        "sik_em_of_r",
        "sik_t_of_r",
        "sik_t_of_em",
        "tpds",
        "ppds",
        "ligand",
        "receptor",
        "em",
        "target",
    ]

    assert all(x in incytr_input.paths.columns for x in expected_columns)

    assert incytr_input.has_p_value == True
    assert incytr_input.has_ppds == True
    assert incytr_input.has_tpds == True
    assert incytr_input.has_umap == False
    assert incytr_input.group_a == "5x"
    assert incytr_input.group_b == "wt"


def test_no_p_value(no_p_value):

    assert no_p_value.has_p_value == False
    assert no_p_value.has_ppds == True


def test_filter_pathways_sigprob(base_pathway_filter: PathwaysFilter):

    pf = base_pathway_filter

    pf.sp_threshold = 0.7

    filtered_a = pf.filter("a", should_filter_umap=False)
    filtered_b = pf.filter("b", should_filter_umap=False)

    paths = pf.all_paths

    assert (
        filtered_a.shape[0]
        == paths[(paths.sigprob_5x >= 0.7) & (paths.afc > 0)].shape[0]
    )
    assert (
        filtered_b.shape[0]
        == paths[(paths.sigprob_wt >= 0.7) & (paths.afc < 0)].shape[0]
    )


def test_nodes_edges(base_pathway_filter: PathwaysFilter, incytr_input):
    a_clusters = incytr_input.clusters[
        incytr_input.clusters["group"] == incytr_input.group_a
    ]

    a_paths = base_pathway_filter.a_data
    a_nodes = load_nodes(a_clusters, node_scale_factor=2)
    a_edges = load_edges(a_nodes, a_paths, 1000, 3)

    astro_astro = [x for x in a_edges if x["data"]["id"] == "astrocytesastrocytes"]

    assert len(astro_astro) == 1

    edge = astro_astro[0]
    assert edge["data"]["weight"] == len(
        a_paths[
            (a_paths["sender"] == "astrocytes") & (a_paths["receiver"] == "astrocytes")
        ]
    )
    assert len(a_nodes) == len(a_clusters.index.unique())


# def test_filter_umap():
#     pass


# def test_remove_afc_filter():
#     pass


# def test_callback_outputs():
#     pass


@pytest.fixture
def no_kinase(datadir):
    return os.path.join(datadir, "no_kinase.csv")


@pytest.fixture
def no_ppds(datadir):
    return os.path.join(datadir, "no_ppds.csv")


@pytest.fixture
def no_sigprob(datadir):
    return os.path.join(datadir, "no_sigprob.csv")


@pytest.fixture
def no_tpds(datadir):
    return os.path.join(datadir, "no_tpds.csv")


@pytest.fixture
def one_p_value(datadir):
    return os.path.join(datadir, "one_p_value.csv")


def test_create_apps(no_kinase, no_ppds, no_sigprob, no_tpds, one_p_value, clusters):
    create_app(clusters_file=clusters, pathways_file=no_kinase)
    create_app(clusters_file=clusters, pathways_file=no_ppds)
    with pytest.raises(ValueError):
        create_app(clusters_file=clusters, pathways_file=no_sigprob)
    create_app(clusters_file=clusters, pathways_file=no_tpds)
    create_app(clusters_file=clusters, pathways_file=one_p_value)
