# -*- coding: utf-8 -*-
"""
Created on Sun May 28 15:38:54 2023

@author: jeremy
"""

import math
from scipy.stats import norm

def margrabeEuroExchageOptions(S1, S2, sigma1, sigma2, D1, D2, T, r, rho):
    # Input parameter validation
    if S1 < 0:
        raise ValueError("S1 - must be positive.")
    elif S2 < 0:
        raise ValueError("S2 - must be positive.")
    elif sigma1 <= 0:
        raise ValueError("sigma1 must be greater than 0.")
    elif sigma2 <= 0:
        raise ValueError("sigma2 must be greater than 0.")
    elif T <= 0:
        raise ValueError("T must be a positive number.")
    elif r <= 0:
        raise ValueError("r must be a positive number.")
    elif D1 <= 0:
        raise ValueError("D1 must be a positive number.")
    elif D2 <= 0:
        raise ValueError("D2 must be a positive number.")
    elif rho < 0:
        raise ValueError("rho must be at least zero.")

    sigmaTot = math.sqrt(sigma1 ** 2 + sigma2 ** 2 - 2 * rho * sigma1 * sigma2)
    d1 = (math.log(S1 / S2) + (D2 - D1 + 0.5 * sigmaTot ** 2) * T) / (sigmaTot * math.sqrt(T))
    d2 = d1 - sigmaTot * math.sqrt(T)

    price = S1 * math.exp(-D1 * T) * norm.cdf(d1) - S2 * math.exp(-D2 * T) * norm.cdf(d2)
    return price

# price = margrabeEuroExchageOptions(130, 105, 0.2, 0.2, 0.06, 0.04, 1, 0.1, 0.5)
# print(price)  # Output: 23.50829459339569
