def weibull(alpha = 1, beta = 1, n = 100):
    """
    Description:
    The weibull function generates a sample from the Weibull distribution with 
    shape parameter alpha and scale parameter beta.

    Parameters:
    - alpha (optional, default=1): The shape parameter of the Weibull 
    distribution. It determines the shape of the distribution curve.
    - beta (optional, default=1): The scale parameter of the Weibull 
    distribution. It influences the spread or width of the distribution.
    - n (optional, default=100): The number of samples to generate from 
    the Weibull distribution.

    Returns:
    - x: A list containing n samples from the Weibull distribution.
    """
    x = [];
    for i in range(n):
        u = np.random.uniform(0,1)
        x.append((-(np.log(u))**(1/alpha))/beta)

    plt.plot(x, alpha=0.7, label = 'exponential rv')
    plt.title(fr'Weibull distribution $W({alpha},{beta})$')
    plt.xlabel('Frequency')
    plt.ylabel(r'x')
    plt.grid(True)
    plt.show()
    return x