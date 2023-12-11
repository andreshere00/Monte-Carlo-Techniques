import numpy as np
import matplotlib.pyplot as plt

def F_distribution(n = 10, m = 9):
    """
    Simulates samples from an F-distribution and visualizes the resulting 
    histogram.

    Parameters:
    - n (int): Numerator degrees of freedom for the F-distribution. Default=10.
    - m (int): Denominator degrees of freedom for the F-distribution. Default=9.

    Returns:
    - x (numpy.ndarray): Array of simulated samples from the F-distribution.
    """
    x = np.array([])

    for i in range(1000):
        z = np.random.chisquare(n, 1)
        y = np.random.chisquare(m, 1)
        x = np.concatenate((x, (y/m) / (z/n)))

    plt.hist(x, bins=30, edgecolor='black')
    plt.title('Histogram of F-distribution')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.show()

    return x