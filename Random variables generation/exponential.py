def exponential(L = 1, n = 100):
    """
    Description:
    The exponential function generates a sample from the exponential 
    distribution with rate parameter L.

    Parameters:
    - L (optional, default=1): The rate parameter of the exponential 
    distribution. It represents the average rate of events occurring per 
    unit of time.
    - n (optional, default=100): The number of samples to generate from the 
    exponential distribution.

    Returns:
    - x: A list containing n samples from the exponential distribution.
    """
    x = [];
    for i in range(n):
        u = np.random.uniform(0,1)
        x.append(-np.log(u)/L)

    plt.plot(x, alpha=0.7, label = 'exponential rv')
    plt.title(fr'Exponential distribution $EXP({L})$')
    plt.xlabel('Frequency')
    plt.ylabel(r'x')
    plt.grid(True)
    plt.show()
    return x