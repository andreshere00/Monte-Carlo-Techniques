import numpy as np
import matplotlib.pyplot as plt

def marsaglia_normal(n):
    """
    Generate n random numbers from a standard normal distribution using Marsaglia's polar method.

    Parameters:
    - n: Number of random numbers to generate.

    Returns:
    - random_numbers: Array of generated random numbers.
    """
    random_numbers = np.zeros(n)
    i = 0

    while i < n:
        u = 2 * np.random.rand() - 1
        v = 2 * np.random.rand() - 1
        s = u**2 + v**2

        if 0 < s < 1:
            x = u * np.sqrt(-2 * np.log(s) / s)
            y = v * np.sqrt(-2 * np.log(s) / s)

            random_numbers[i] = x
            i += 1

            if i < n:
                random_numbers[i] = y
                i += 1

    plt.hist(random_numbers, bins=30, density=True, alpha=0.5)
    plt.title("Generated Normal Random Numbers (Marsaglia's Method)")
    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.show()

    return random_numbers