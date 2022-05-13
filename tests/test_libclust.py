import numpy as np
import pandas as pd
import pytest

from pymap.libclust import (
    check_volume,
    get_clust,
    validate_clust,
)

def test_volume():
    """test correct calculation of the volume in a trivial case"""
    d = {'a' : [0,0,1], "b" : [1,0,1]}
    df = pd.DataFrame(d)
    assert check_volume(df,df.shape[1]) == 2

def test_at_clust():
    """test correct at clustering in a trivial case"""
    d = {'a': [1, 0, 1, 1, 1], 'b': [-1, 4, -1, 4, 0]}
    full_df = pd.DataFrame(d)
    mapping = np.array([0,1])
    at_df = get_clust(full_df,mapping)
    expected_data = {'a': [0,1,1,1], 'b': [4,-1,0,4], 'records': [1,2,1,1]}
    expected_at_df = pd.DataFrame(expected_data)
    assert at_df.equals(expected_at_df)

def test_cg_clust():
    """test correct cg clustering in a trivial case"""
    at_d = {'a': [0,1,1,1], 'b': [4,-1,0,4]}
    at_df = pd.DataFrame(at_d)
    mapping = np.array([0])
    cg_df = get_clust(at_df,mapping)
    expected_data = {'a': [0,1], 'records': [1,3]}
    expected_cg_df = pd.DataFrame(expected_data)
    assert cg_df.equals(expected_cg_df)

def test_cg_clust_error():
    """cg clustering should raise an error if equal rows are present"""
    # the row (1,-1) is a duplicate
    wrong_cg_clust = {'a': [1, 0, 1, 1, 1],
                      'b': [-1, 4, -1, 4, 0],
                      "records": [1,2,3,4,5]}
    df = pd.DataFrame(wrong_cg_clust)
    with pytest.raises(Exception, match="Duplicate row detected in dataframe"):
        validate_clust(df)