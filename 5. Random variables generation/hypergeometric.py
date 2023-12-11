import matplotlib.pyplot as plt
import random

def hypergeometric(N=500, D=50, n=50, size=1000):
    """
    Simulate a hypergeometric distribution and plot its histogram.

    Parameters:
    - N (int, optional): Population size. Defaults to 500.
    - D (int, optional): Number of success states in the population. Defaults to 50.
    - n (int, optional): Number of draws. Defaults to 50.
    - size (int, optional): Number of simulations. Defaults to 1000.

    Returns:
    - random_values (list): A list of random variables generated from the 
    hypergeometric distribution.
    """
    random_values = []
    for _ in range(size):
        x = 0
        N_prime = N
        d = D
        c = N - D

        for _ in range(n):
            u = random.uniform(0, 1)  # Generar un número aleatorio en el intervalo [0, 1)

            if u <= d / N_prime:
                x += 1
                N_prime -= 1
                d -= 1
            else:
                N_prime -= 1
                c -= 1
        random_values.append(x)
    # Crear un histograma para visualizar los datos
    plt.hist(random_values, bins=range(min(random_values), max(random_values) + 1, 1), alpha=0.7, rwidth=0.85)
    plt.title(f'Distribución Hipergeométrica (N={N}, D={D}, n={n})')
    plt.xlabel('Número de Éxitos')
    plt.ylabel('Frecuencia')
    plt.grid(True)
    plt.show()

    return random_values