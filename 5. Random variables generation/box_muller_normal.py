import numpy as np
import matplotlib.pyplot as plt

def box_muller_normal(n):
    """
    Generate n random numbers from a standard normal distribution using Box-Muller method.

    Parameters:
    - n: Number of random numbers to generate.

    Returns:
    - random_numbers: Array of generated random numbers.
    """
    u1 = np.random.rand(n)
    u2 = np.random.rand(n)

    z1 = np.sqrt(-2 * np.log(u1)) * np.cos(2 * np.pi * u2)
    z2 = np.sqrt(-2 * np.log(u1)) * np.sin(2 * np.pi * u2)

    random_numbers = np.concatenate([z1, z2])

    plt.hist(random_numbers, bins=30, density=True, alpha=0.5)
    plt.title('Generated Normal Random Numbers')
    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.show()

    return random_numbers