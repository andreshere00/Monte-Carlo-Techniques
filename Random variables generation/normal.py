import numpy as np
import matplotlib.pyplot as plt

def normal(n=1000):
    """
    Normal distribution with parameters (mu = 0, sigma = 1).

    Parameters:
    - n: number of samples to generate.

    Returns:
    - Xs: generated samples.
    """
    Xs = np.zeros(n)

    for i in range(len(Xs)):
        u_1 = np.random.uniform(0, 1)
        u_2 = np.random.uniform(0, 1)
        x_i = -np.log(1/u_1 - 1)
        y_i = u_2 * 4/np.sqrt(2*np.pi) * np.exp(-x_i) / ((1+np.exp(-x_i))**2)

        while y_i > np.exp(-(x_i**2)/2) / np.sqrt(2*np.pi):  
          # dnorm(x_i, mean=0, sd=1) in R is equivalent 
          # to exp(-(x_i^2)/2) / sqrt(2*pi) in Python
            u_1 = np.random.uniform(0, 1)
            u_2 = np.random.uniform(0, 1)
            x_i = -np.log(1/u_1 - 1)
            y_i = u_2 * 4/np.sqrt(2*np.pi) * np.exp(-x_i) /((1+np.exp(-x_i))**2)

        Xs[i] = x_i

    # Plot the histogram
    plt.hist(Xs, bins = int(n/50), density=True, alpha=0.5)
    plt.title(r'Normal distribution $N(0,1)$')
    plt.xlabel(r'x')
    plt.ylabel('Probability')
    plt.show()

    return Xs