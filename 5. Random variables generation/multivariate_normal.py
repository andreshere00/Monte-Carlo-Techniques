import numpy as np

def multivariate_normal(mu = [0,0], sigma = [[2,1],[1,1]], n = 2, m = 10):
    """
    Generates samples from a multivariate normal distribution.

    Parameters:
    - mu (numpy.ndarray): Mean vector of the multivariate normal distribution.
    Default = 0
    - sigma (numpy.ndarray): Covariance matrix of the multivariate normal 
    distribution. Default = 0.5
    - n (int): Number of dimensions in the multivariate normal distribution. 
    Default = 10.
    - size (int): Number of samples to generate. Default = 1000.

    Returns:
    - samples (list): List of generated samples.
    """

    # Cholesky decomposition of the covariance matrix
    L = np.linalg.cholesky(sigma)

    # Initialize list
    samples = []

    # Generate samples
    for i in range(m):
        z = np.random.normal(0, 1, n)
        x = mu + np.dot(L, z)
        samples.append(x.tolist())

    sns.histplot(samples, bins=30, kde=True).\
                 set(title=rf'Multivariate Normal distribution $N({mu},{sigma})$')

    return samples