import numpy as np
import matplotlib.pyplot as plt

def uniform_cont(a = 10, b = 25, n = 100):
    """
    Simulate a continuous uniform distribution based on "a" and "b" parameters.

    Parameters:
    - a (int): minimum value of the uniform distribution.
    - b (int): maximum value of the uniform distribution.
    - n (int): trials of the experiment, size of the generated sample.

    Returns:
    - x (list): a list which contains the generated random variables.
    """
    x = []
    for i in range(n):
        u = np.random.uniform(0,1)
        x.append(a + u*(b-a))

    plt.plot(x, alpha=0.7)
    plt.title(fr'Uniform distribution $U({a},{b})$')
    plt.ylim(a*0.8, b*1.2)
    plt.xlabel('Frequency')
    plt.ylabel(r'x')
    plt.grid(True)
    plt.show()
    return x