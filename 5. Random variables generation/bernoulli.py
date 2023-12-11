import numpy as np

def bernoulli(n=10,p=0.5):
    """
    Generates a Bernoulli distribution of the following parameters:
    Parameters:
    - n (list): number of samples to generate. Default is 10.
    - p (float): probability of success. Default is 0.5.
    Returns:
    - b (list): with the result of the experiment
    """
    b = []
    for i in range(n):
        u = np.random.uniform(0,1)
        b.append(1) if u <= p else b.append(0)
    return b

bernoulli(10,0.5)