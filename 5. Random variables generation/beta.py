import numpy as np
import matplotlib.pyplot as plt

def target_beta_distribution(x, alpha=2, beta=2):
    """ 
    Mathematical expression of the fdp for the beta distribution.
    """
    return (x**(alpha-1) * (1-x)**(beta-1)) / ((2**(alpha+beta-1)) \
    / (np.math.factorial(alpha-1) * np.math.factorial(beta-1)))

def beta(n=100, alpha=2, beta=2):
    """
    Generate samples from the beta distribution using the rejection sampling method.

    Parameters:
    - n: Number of samples to generate (default=100).
    - alpha: Shape parameter (default=2).
    - beta: Shape parameter (default=2).

    Returns:
    - samples: List of n samples generated using rejection sampling.
    """
    samples = []

    # Constant for the envelope function
    c = target_beta_distribution(1, alpha, beta) 

    while len(samples) < n:
        u1, u2 = np.random.uniform(0, 1, 2)

        x = u1
        y = c * u2

        if y <= target_beta_distribution(x, alpha, beta):
            samples.append(x)

    plt.hist(samples, bins=30, density=False, alpha=0.5, label='Generated Samples')
    x_values = np.linspace(0, 1, 100)
    plt.title('Rejection Sampling for Beta Distribution')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.show()

    return samples