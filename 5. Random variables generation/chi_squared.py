import numpy as np
import matplotlib.pyplot as plt

def chi_squared(n=20, N=1000):
    """
    Description:
    simulates a chi-squared distribution with a specified number of degrees of
    freedom and a specified number of iterations (trials).

    Parameters:
    - n (int): Degrees of freedom for the chi-squared distribution. 
    Default = 20.
    - N (int): Number of repetitions, i.e., the number of random samples to 
    generate. Default = 1000

    Returns:
    - chi2 (list): a list of generated samples.
    """
    chi2 = np.array([])

    for i in range(1, N + 1):
        if n % 2 == 0:
            if np.floor(n / 2) > 1:
                u = np.random.uniform(0, 1, int(n / 2))
                chi2 = np.concatenate((chi, [-2 * np.sum(np.log(u))]))
        else:
            if np.floor(n / 2) >= 1:
                u = np.random.uniform(0, 1, int((n - 1) / 2))
                z = np.random.normal(0, 1, 1)
                chi2 = np.concatenate((chi, [-2 * np.sum(np.log(u)) + z**2]))

    plt.hist(chi2, bins=30, edgecolor='black')
    plt.show()

    return chi2