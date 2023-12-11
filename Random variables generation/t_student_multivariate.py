import numpy as np
import matplotlib.pyplot as plt




def t_student_multivariante(m=100, n=3, sigma=1, plot=True):
  def generate_x(m, sigma):
    return np.random.normal(0, sigma, m)

  def generate_y(n):
    return np.random.chisquare(n-1)

  def calculate_t(x, y, n):
    return x / y / n

  def plot_distribution(values):
    plt.hist(values, bins=20, density=True, alpha=0.7, color='g')
    plt.title('t-Student Multivariante')
    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.show()

  x_values = generate_x(m, sigma)
  y_value = generate_y(n)

  t_values = [calculate_t(x, y_value, n) for x in x_values]

  x_output = tuple(t_values)

  if plot == True:
    plot_distribution(x_output)


  return x_output

resultados= t_student_multivariante(plot=True)1