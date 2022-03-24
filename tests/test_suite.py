""" test_suite

This file contains a few python tests to check the correct installation of pymap
"""

import os
import numpy as np
import pandas as pd
# import utils modules
from utils import check_volume, get_clust

def test_volume():
    """test correct calculation of the volume in a trivial case"""
    np_arr = np.array([[0,1],[0,0],[1,1]])
    df = pd.DataFrame(np_arr)
    print(df)
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

#def test_cg_clust_error():
    """cg clustering should raise an error if equal rows are present"""
