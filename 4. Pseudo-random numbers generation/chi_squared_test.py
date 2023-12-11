import numpy as np
from scipy.stats import chi2

def chi_squared_test(data, n_intervals = 5, alpha = 0.05):
    """
    Perform a chi-squared goodness-of-fit test on a given dataset.

    Parameters:
    - data: List or array of the observed values.
    - n_intervals: Number of intervals for expected frequencies (default is 5).
    - alpha: Significance level for the test (default is 0.05).

    Returns:
    - chi_sq: Chi-squared statistic.
    - critical_region: Critical region for the test.
    - chi_2_critical: Critical value for the given significance level.
    - p_value: p-value for the chi-squared test.

    Notes:
    - The function divides the range of the data into 'n_intervals' intervals 
    and calculates observed frequencies.
    - Conducts a chi-squared test to assess if the observed frequencies differ 
    significantly from expected frequencies.
    """

    n = len(data)
    data = sorted(data)
    exp_freq = np.linspace(0, 1, num=n_intervals+1)
    obs_freq = [[i for i in data if exp_freq[k] <= i < exp_freq[k+1]] \
                for k in range(n_intervals)]
    chi_sq = 0
    for j in range(n_intervals):
        chi_sq += (len((obs_freq[j]))-n*(exp_freq[1]))**2/(n*(exp_freq[1]))

    # outputs
    critical_region = [chi2.ppf(1 - alpha, n_intervals-1),np.infty]
    chi_2_critical = critical_region[0]
    p_value = chi2.sf(chi_sq * np.sqrt(n_intervals), n_intervals)

    return chi_sq, critical_region, chi_2_critical, p_value