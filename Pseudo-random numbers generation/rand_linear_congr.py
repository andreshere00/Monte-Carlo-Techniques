import time
import pandas as pd

def rand_linear_congr(seed=int(time.time()), a=5, b=3, m=16, N=10):
    if type(seed)!=int or int(seed)%1 != 0:
        raise Exception('The seed is not a natural number.')

    X = int(seed)
    Y = a * X + b

    data = pd.DataFrame({"i": [0], "X": [X], "Y": [Y]})

    for i in range(1, N + 1):
        X = Y % m
        Y = a * X + b
        data.loc[len(data.index)] = [i, X, Y] 

    data = pd.DataFrame(data)
    data = data.set_index('i')

    return data