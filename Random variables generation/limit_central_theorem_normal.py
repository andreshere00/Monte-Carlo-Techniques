def limit_central_theorem_normal(n):
    """
    Generate n random numbers from a standard normal distribution using the Central Limit Theorem.

    Parameters:
    - n: Number of random numbers to generate.

    Returns:
    - random_numbers: Array of generated random numbers.
    """
    random_numbers = np.zeros(n)

    for i in range(n):
        # Suma de 12 variables aleatorias uniformes [-0.5, 0.5]
        sum_uniforms = np.sum(np.random.uniform(-0.5, 0.5, 12))
        random_numbers[i] = sum_uniforms

    # Plot the histogram of generated random numbers
    plt.hist(random_numbers, bins=30, density=True, alpha=0.5)
    plt.title("Generated Normal Random Numbers (Central Limit Theorem)")
    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.show()

    return random_numbers