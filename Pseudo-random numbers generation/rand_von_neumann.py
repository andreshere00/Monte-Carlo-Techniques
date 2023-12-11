import time
import pandas as pd

def rand_von_neumann(seed=int(time.time()), n2=4, N=10):
    """
    Generates a sequence of pseudo-random numbers using the Von Neumann method.

    Parameters:
    - seed (int): Initial seed for pseudo-random number generation. By default, 
    the current time is used.
    - n2 (int): Number of central digits to consider in each iteration. Default 
    is 4.
    - N (int): Total number of iterations to generate the sequence. Default is 
    10.

    Returns:
    - pandas.DataFrame: A DataFrame containing the generated sequence with 
    columns 'i', 'X' (current number), and 'X^2' (pseudo-random number).

    Raises:
    - Exception: If the seed is not a natural number.
    """
    
    if int(seed)%1 != 0:
        raise Exception("The seed value must be a natural number.")
        seed = int(seed) % (10 ** n2)

    aux1 = 10 ** (2 * n2 - n2 / 2)
    aux2 = 10 ** (n2 / 2)

    X = seed
    X2 = X ** 2

    data = pd.DataFrame({"i": [0], "X": [X], "X^2": [X2]})
    for i in range(1,N):
        X = int((X2 - (X2 // aux1) * aux1) // aux2)
        X2 = X ** 2
        data.loc[len(data.index)] = [i, X, X2] 

    data = data.set_index('i')
    return data