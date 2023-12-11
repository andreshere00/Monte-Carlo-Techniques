import numpy as np
import matplotlib.pyplot as plt

def geometric(p=0.1, n=1000):
    """
    Generate a discrete geometric random variable.

    Parameters:
    - p (float, optional): Probability of success. Default is 0.1.
    - n (int, optional): Number of iterations. Default is 1000.

    Returns:
    - x (list): A list of generated random variables.
    """
    x = []

    for j in range(1, n + 1):
        y = 0
        while True:
            u = np.random.rand(1)
            if u <= p:
                x.append(y)
                break
            else:
                y = y + 1

    plt.hist(x, bins='auto')
    plt.title('Geometric Random Variables distribution')
    plt.show()

    return x