import numpy as np
import matplotlib.pyplot as plt

def rejection_generalized_method(f, g, a = 2, n = 1000):
    """
    Generalized rejection sampling algorithm.

    Parameters:
    - f: Function representing the target probability density function.
    - g: Function representing the envelope (safety) probability density function.
    - a: Upper bound for the envelope function (default=2).
    - n: Number of samples to generate (default=1000).

    Returns:
    - samples: List of n samples generated using generalized rejection sampling.
    """

    samples = []

    while len(samples) < n:
        x = np.random.uniform()  # Generate x from the proposal distribution g

        y = np.random.uniform(0, a * g(x))  # Generate y from U(0, a * g(x))

        if y <= f(x):
            samples.append(x)

    plt.hist(samples, bins=30, density=True, alpha=0.5, label='Generated Samples')
    x_values = np.linspace(0, 3, 100)
    plt.plot(x_values, target_distribution(x_values), label='Target Distribution')
    plt.title('Generalized Rejection Sampling')
    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()

    return samples

def target_distribution(x):
    return np.exp(-x**2)
def envelope_distribution(x):
    return 1