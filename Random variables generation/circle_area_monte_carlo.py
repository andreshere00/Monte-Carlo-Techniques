import random
import numpy as np

def circle_area_monte_carlo(n=1e6):
    """
    Estimate the area of a unit circle using the Monte Carlo method.

    Parameters:
    - n (int): The number of random points to generate for the estimation. 
    Default is 1e6.

    Returns:
    float: Estimated area of the unit circle.

    Note: Increasing the value of 'n' generally improves the accuracy of the 
    estimation.
    """
    inside_circle = 0
    for _ in range(n):
        x = random.uniform(-1, 1)
        y = random.uniform(-1, 1)
        # Verify if the point (x,y) is inside the circle.
        if x**2 + y**2 <= 1:
            inside_circle += 1

    est_circle_area = inside_circle / n
    est_square_area = 4
    est_area = est_square_area * est_circle_area

    return est_area