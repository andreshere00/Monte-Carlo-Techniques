import random
import matplotlib.pyplot as plt

def negative_binomial(m = 10, p = 0.5, size = 10000):
    """
    Simulate a negative binomial distribution and plot its histogram.

    Parameters:
    - m (int, optional): Number of successes. Defaults to 10.
    - p (float, optional): Probability of success on each trial. Defaults to 0.5.
    - size (int, optional): Number of simulations. Defaults to 10000.

    Returns:
    - random_values (list): A list of random variables generated from the 
    negative binomial distribution.
    """

    random_values = []
    for _ in range(size):
        trials = 0  # Count of Bernoulli trials
        successes = 0  # Count of successes

        while successes < m:
            trial_result = random.random()  # Generate a random number between 0 and 1
            if trial_result < p:
                successes += 1
            trials += 1

        random_values.append(trials)

    plt.hist(random_values, bins=range(min(random_values), max(random_values) + 1, 1), alpha=0.7, rwidth=0.85)
    plt.title(f'Negative Binomial Distribution (m={m}, p={p})')
    plt.xlabel('Number of Trials')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

    return random_values