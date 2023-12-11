def generate_rv(n, X, p):
    """
    Generate random variables based on a discrete probability distribution.
    If X is binary and p has only two possible values, it is a Bernouilli 
    distribution.

    Parameters:
    - n (int): Number of random variables to generate.
    - X (list): List of unique values of the random variable.
    - p (list): List of probabilities associated with each value in X.

    Returns:
    list: A list of n random variables sampled from the given distribution.
    """

    FD = list((zip(*sorted(zip(X, p), key=lambda x: x[1], reverse=True))))
    FD[1] = np.cumsum(FD[1])

    x = []

    u_values = np.random.uniform(0, 1, n)
    x = [FD[0][np.searchsorted(FD[1], u)] for u in u_values]
    return x

generate_rv(10, [0, 1, 2, 3], [0.1, 0.2, 0.5, 0.2])