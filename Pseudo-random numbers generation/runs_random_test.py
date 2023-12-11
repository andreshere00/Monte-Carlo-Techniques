import numpy as np
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