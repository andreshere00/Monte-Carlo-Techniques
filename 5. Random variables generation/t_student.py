import numpy as np
import matplotlib.pyplot as plt

def t_student(n = 1000):
    """
    Description:
    - Simulate a T-Student distribution.
    Parameters:
    - n (int): number of trials.
    Returns:
    - x (list): a list with the sample size (n)
    """
    x = np.array([])
    for i in range(1000):
        z = np.random.normal(0, 1, 1)
        y = np.random.chisquare(n, 1)
        x = np.concatenate((x, z / np.sqrt(y / n)))

    plt.hist(x, bins=30, edgecolor='black')
    plt.show()

    return x