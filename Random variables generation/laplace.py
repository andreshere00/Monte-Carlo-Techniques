def laplace(L=1, b=1, n=100):
    """
    Generate random samples from a Laplace distribution.

    Parameters:
    - L: Location parameter (default=1).
    - b: Scale parameter (default=1).
    - n: Number of samples to generate (default=100).

    Returns:
    - x: Array of generated samples.
    """
    u1 = np.random.uniform(0, 1, n)
    u2 = np.random.uniform(0, 1, n)
    x = np.where(u1 <= 0.5, mu - b * np.log(u2), mu + b * np.log(u2))

    plt.hist(x, bins = int(n/10), density=True, alpha=0.5)
    plt.title(fr'Laplace distribution $L({L},{b})$')
    plt.xlabel(r'x')
    plt.ylabel('Probability')
    plt.show()
    return x