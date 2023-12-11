import time
import pandas as pd

def rand_lehmer(seed=int(time.time()), n=4, mu=76, k=2, N=10):
    """
    Generates a sequence of pseudo-random numbers using the Lehmer method.

    Parameters:
    - seed (int): Initial seed for pseudo-random number generation. By default, 
    the current time is used.
    - n (int): Number of digits in the random number.
    - mu (int): Parameter for the Lehmer method.
    - k (int): Number of digits of mu.
    - N (int): Total number of iterations to generate the sequence. Default 
    is 10.

    Returns:
    - pandas.DataFrame: A DataFrame containing the generated sequence with 
    columns 'i', 'X' (current number), 'X * mu' (product of current number and 
    mu), 'Y', and 'Z'.

    Raises:
    - Exception if the seed is not a natural number.
    """
    if type(seed)!=int or int(seed)%1 != 0:
        raise Exception('The seed is not a natural number.')

    aux = 10 ** n
    seed = float(seed) % 10 ** n
    mu = mu % 10 ** k

    X = int(seed)
    Xmu = X * mu
    Y = Xmu // aux
    Z = Xmu % aux

    data = pd.DataFrame({"i": [0], "X": [X], "X * mu": [Xmu],
                         "Y": [Y], "Z": [Z]})

    for i in range(1, N + 1):
        X = Z - Y
        Xmu = X * mu
        Y = Xmu // aux
        Z = Xmu % aux
        data.loc[len(data.index)] = [i, X, Xmu, Y, Z] 

    data = data.set_index('i')
    return data