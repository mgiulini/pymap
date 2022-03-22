import numpy as np 
import pandas as pd

def check_volume(dataframe,ncols):
    """
    calculate the volume associated to each degree of freedom
    """
    V = np.unique([len(np.unique(dataframe.iloc[:,n])) for n in range(ncols-1)], return_counts=True)
    if len(V[1]) != 1:
        print("non-constant volume detected, V = ", V, "choosing the most likely one")
        agmax = np.argmax(V[1])
        V = V[0][agmax]
    else:
        print("constant volume detected, V = ", V)
        V = V[0][0]
    print("V = ", V)
    return V