import os
import pytest
from incytr_viz.util import PathwaysFilter

os.environ["INCYTR_PATHWAYS"] = (
    "/home/icossentino/code/incytr-viz/tests/data/mc38/mc38_hegs_degs_proteomics_jan2025_ligand-target_sample.tsv"
)
os.environ["INCYTR_CLUSTERS"] = (
    "/home/icossentino/code/incytr-viz/tests/data/mc38/pop.csv"
)

from incytr_viz.app import get_pathways, cache


@pytest.fixture
def pi():
    group_a = "10days"
    group_b = "14days"
    fpath = os.environ["INCYTR_PATHWAYS"]
    yield get_pathways(fpath=fpath, group_a=group_a, group_b=group_b)
    cache.clear()


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


def test_get_pathways(pi):

    paths = pi.paths

    assert all(
        [
            x in paths.columns
            for x in [
                "path",
                "sender",
                "receiver",
                "sigprob_14days",
                "sigprob_10days",
                "afc",
                "tpds",
                "ppds",
                "ligand",
                "receptor",
                "em",
                "target",
            ]
        ]
    )

    assert pi.has_p_value == True
    assert pi.has_ppds == True
    assert pi.has_tpds == True
    assert pi.has_umap == False
    assert pi.group_a == "10days"
    assert pi.group_b == "14days"
    assert (paths == cache.get("pathways").paths).all().all()


def test_filter_pathways(pi):

    pf = PathwaysFilter(
        all_paths=pi.paths,
        group_a_name=pi.group_a,
        group_b_name=pi.group_b,
        sp_threshold=0.7,
    )

    filtered_a = pf.filter("a", should_filter_umap=False)
    filtered_b = pf.filter("b", should_filter_umap=False)

    assert filtered_a.shape[0] == pi.paths[pi.paths.sigprob_10days >= 0.7].shape[0]
    assert filtered_b.shape[0] == pi.paths[pi.paths.sigprob_14days >= 0.7].shape[0]
    assert filtered_a.shape[0] != filtered_b.shape[0] != len(pi.paths)
