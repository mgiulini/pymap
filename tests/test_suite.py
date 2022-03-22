""" test_suite

This file contains a few python tests to check the correct installation of pymap
"""

import os
import numpy as np
import pandas as pd
# import utils modules
from utils import check_volume

def test_volume():
    np_arr = np.array([[0,1],[0,0],[1,1]])
    df = pd.DataFrame(np_arr)
    print(df)
    assert check_volume(df,df.shape[1]) == 2