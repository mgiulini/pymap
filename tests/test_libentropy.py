"""Test the libentropy library."""
import numpy as np
import pandas as pd
import pytest

from pymap.libclust import get_clust

from pymap.libentropy import (
    calculate_pbar,
    calculate_smap,
    calculate_smap_inf,
    calculate_entropies,
)


def test_calculate_pbar():
    """Test correct calculation of pbar."""
    full_d = {
        'a': [0, 1, 1, 1],
        'b': [4, -1, -1, 4]
        }
    df = pd.DataFrame(full_d)
    at_mapping = np.array([0, 1])
    at_df = get_clust(df, at_mapping)
    nobs = df.shape[0]
    mapping = np.array([0])
    cg_df = get_clust(df, mapping)
    # keeping only "a" gives rise to confs 0, with 1 at conf (1 original), and
    # 1 with 2 at confs (4 original). pbar([0,4]) = 0.25,  pbar([1,x]) = 0.75/2
    expected_pbar = {
        "a": [0, 1],
        "omega_1": [1, 2],
        "p_bar": [0.25, 0.375]
        }
    expected_pbar_df = pd.DataFrame(expected_pbar)
    pbar = calculate_pbar(at_df, cg_df, nobs, mapping)
    assert pbar.equals(expected_pbar_df)


def test_smap_zero():
    """Test correct calculation of smap."""
    df = pd.DataFrame(
        {
            'a': [0, 1, 1, 0],
            'b': [4, -1, -1, 4]
        }
    )
    at_mapping = np.array([0, 1])
    at_df = get_clust(df, at_mapping)
    pr = at_df["records"]/df.shape[0]
    mapping = np.array([0])
    p_bar = pd.DataFrame(
        {
            "a": [0, 1],
            "omega_1": [2, 2],
            "p_bar": [0.5, 0.5]
        }
    )
    smap = calculate_smap(at_df, mapping, pr, p_bar)
    expected_smap = 0.0
    assert expected_smap == smap


def test_smap_inf_error():
    """Test correct exceptions during the calculation of smap_inf."""
    nat, ncg = 9, 10
    with pytest.raises(ValueError, match="n (9) < N (10)"):
        calculate_smap_inf(nat, ncg, 0.0, 0.0, 2)
    nat, ncg = 10, 9
    with pytest.raises(ValueError, match="hs_at (0.0) < hs_cg (1.0)"):
        calculate_smap_inf(nat, ncg, 0.0, 1.0, 2)


def test_smap_inf_zero():
    """Test correct calculation of smap_inf."""
    nat, ncg = 2, 1
    V = 3
    hs_at = np.log(9)
    hs_cg = np.log(3)
    smap_inf = calculate_smap_inf(nat, ncg, hs_at, hs_cg, V)
    expected_smap_inf = 0.0
    assert expected_smap_inf == smap_inf


def test_hs_hk_error():
    """Test error in calculate_entropies."""
    cg_clust = pd.DataFrame({'a': [0, 1], 'b': [4, -1]})
    exp_str = "cg_clust does not have a 'records' column"
    with pytest.raises(ValueError, match=exp_str):
        calculate_entropies(cg_clust)
    cg_clust = pd.DataFrame({'a': [0, 1], 'b': [4, -1], 'records': [1, -1]})
    exp_str = "'records' column in cg_clust contains negative values"
    with pytest.raises(ValueError, match=exp_str):
        calculate_entropies(cg_clust)


def test_hs_hk():
    """Test correct calculation of hs and hk."""
    cg_clust = pd.DataFrame({'a': [0, 1], 'b': [4, -1], 'records': [1, 1]})
    hs, hk = calculate_entropies(cg_clust)
    exp_hs = np.log(2)
    exp_hk = 0.0
    assert exp_hs == hs
    assert exp_hk == hk
