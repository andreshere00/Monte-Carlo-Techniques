def runs_wald_wolfowitz_test(data,alpha=0.05):
    """
    Perform a Wald-Wolfowitz runs test on a given dataset.

    Parameters:
    - data: List or array of the sample values.
    - alpha: Significance level for the test (default is 0.05).

    Returns:
    - Z: Wald-Wolfowitz Z statistic.
    - p_value: p-value for the Wald-Wolfowitz test.
    - Z_critical: Critical value for the given significance level.
    - critical_region: Critical region for the test.

    Notes:
    - The function standardizes the data to have mean 0 and standard deviation 1.
    - Calculates the number of runs (R) in the sequence.
    - Assesses if the number of runs significantly deviates from the expected number of runs.
    """

    n = len(data)
    mean = np.mean(data)
    sigma = np.std(data)
    data = [(x - mean) / sigma for x in data]
    R = 1; N_neg=0; N_pos=0

    if data[0] < 0: N_neg += 1
    else: N_pos += 1

    for i in range(1,n-1):
        if data[i+1] < data[i] and data[i] >= data[i-1] or \
            data[i+1] >= data[i] and data[i] < data[i-1]:
            R += 1

    N_neg = sum([1 for x in data if x < 0])
    N_pos = sum([1 for x in data if x>= 0])

    mu = (2*N_pos*N_neg)/(N_pos + N_neg) + 1
    desv2= ((mu - 1)*(mu - 2))/(N_pos + N_neg - 1)

    Z = (R-mu)/np.sqrt(desv2)
    p_value = norm.cdf(Z)
    CR = [[np.infty, -norm.ppf(alpha/2)], [norm.ppf(alpha/2), np.infty]]
    Z_critical = CR[1][0]

    return Z, p_value, Z_critical, CRimport numpy as np
from scipy.stats import norm

def runs_random_test(data,alpha=0.05):
    n = len(data)
    R = 1
    for i in range(1,n-1):
        if data[i+1] < data[i] and data[i] >= data[i-1] or \
        data[i+1] >= data[i] and data[i] < data[i-1]:
            R += 1

    mu = (2*n-1)/3
    desv2= (16*n - 29)/90
    Z = (R-mu)/np.sqrt(desv2)

    p_value = norm.cdf(Z)
    CR = [[np.infty, -norm.ppf(alpha/2)], [norm.ppf(alpha/2), np.infty]]
    Z_critical = CR[1][0]

    return Z, Z_critical, p_value, CR