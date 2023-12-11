import numpy as np
from scipy.stats import ksone

def kolmogorov_smirnov_test(data, alpha = 0.05):
    """
    Perform the Kolmogorov-Smirnov (KS) test on a given dataset.

    Parameters:
    - data: List or array of the sample data.
    - alpha: Significance level (default is 0.05).

    Returns:
    - D: Kolmogorov-Smirnov statistic.
    - D_critical: Critical value for the given significance level.
    - critical_region: Critical region for the test.
    - p_value: p-value for the KS test.
    """
    n = len(data)
    data =  sorted(data)
    D_plus = 0; D_minus = 0
    dist = []

    for i in range(n):
        D_plus = abs((i+1)/n - data[i])
        D_minus = abs(i/n - data[i])
        dist.append(max(D_plus, D_minus))

    # outputs
    D = max(dist)
    critical_region = [ksone.ppf(1-alpha/2,n),np.infty]
    D_critical = ksone.ppf(1-alpha/2,n) # tabulate this function is complicated
    p_value = 1 - ksone.sf(D * np.sqrt(n), n)

    return D, D_critical, critical_region, p_value