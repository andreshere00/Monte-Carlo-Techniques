def erlang(p = 100, L = 1, n = 100):
    """
    Description:
    The erlang function generates a sample from the Erlang distribution, also 
    known as the gamma distribution with an integer shape parameter p.

    Parameters:
    - p (optional, default=100): The shape parameter of the Erlang distribution.
     It represents the number of events for which the waiting time is modeled.
    - L (optional, default=1): The rate parameter of the Poisson process. 
    It indicates the average rate of events per unit of time.
    - n (optional, default=100): The number of samples to generate from the 
    Erlang distribution.

    Returns:
    - X: A list containing n samples from the Erlang distribution.
    """
    X =[]
    for i in range(n):
        x = 0
        for i in range(p):
            u = np.random.uniform(0,1)
            x -= np.log(u)/L
        X.append(x)

    plt.plot(X)
    plt.title(fr'Erlang distribution $\gamma({p},{L})$')
    plt.xlabel('Frequency')
    plt.ylabel(r'x')
    plt.grid()
    plt.show()

    return X