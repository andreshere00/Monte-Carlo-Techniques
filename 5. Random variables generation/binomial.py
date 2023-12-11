import matplotlib.pyplot as plt

def binomial(n=10, p=0.5, trials=10):
    """
    Generates a Binomial distribution using the Bernoulli distribution.

    Parameters:
    - n (int): Number of samples to generate. Default is 10.
    - p (float): Probability of success for each Bernoulli trial. Default is 0.5.
    - trials (int): Number of Bernoulli trials per sample. Default is 1.

    Returns:
    - binomial_dist (list): List with the result of the Binomial distribution.
    """
    binomial_dist = []

    for _ in range(n):
        # Summing up the results of 'trials' number of Bernoulli trials
        result = sum(bernoulli(trials, p))
        binomial_dist.append(result)

    # Plot histogram
    plt.hist(binomial_dist, bins='auto')
    plt.title('Binomial Random Distribution')
    plt.show()

    return binomial_dist