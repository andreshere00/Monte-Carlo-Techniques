import numpy as np
import matplotlib.pyplot as plt

def rejection_sampling(target_distribution, proposal_distribution = 1, 
                       a = 0, b = 5, c = 1, n = 1000):
    """
    Rejection sampling algorithm.

    Parameters:
    - target_distribution: Function representing the target probability density 
    function.
    - proposal_distribution: Function representing the proposal probability 
    density function.
    - a: Lower bound of the range for x.
    - b: Upper bound of the range for x.
    - c: Upper bound of the range for y (considered as the maximum value of the 
    target distribution).
    - n: Number of samples to generate.

    Returns:
    - samples: List of n samples generated using rejection sampling.
    """

    samples = []

    while len(samples) < n:
        u1 = np.random.uniform(0, 1)
        u2 = np.random.uniform(0, 1)

        x = a + (b - a) * u1
        y = c * u2

        if y <= target_distribution(x) / proposal_distribution(x):
            samples.append(x)

    return samples

def target_distribution(x):
    return np.exp(-x)
def proposal_distribution(x):
    return 1

a = 0; b = 5; c = 1; n = 1000
samples = rejection_sampling(target_distribution, proposal_distribution, a, b, c, n)

plt.hist(samples, bins=30, density=True, alpha=0.5, label='Generated Samples')
x_values = np.linspace(a, b, 100)
plt.plot(x_values, target_distribution(x_values), label='Target Distribution')
plt.title('Rejection Sampling')
plt.xlabel('Value')
plt.ylabel('Probability Density')
plt.legend()
plt.show()