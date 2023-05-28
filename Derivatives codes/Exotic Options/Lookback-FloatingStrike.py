# -*- coding: utf-8 -*-
"""
Created on Sun May 28 15:41:25 2023

@author: jeremy
"""
# You can use the fltStrikeLookback function in Python by passing the required input parameters:
# callPut, S, S_max_min, T, D, r, and sigma. 
# The function will return the calculated call or put option price based on the provided inputs.

import math
from scipy.stats import norm

def fltStrikeLookback(callPut, S, S_max_min, T, D, r, sigma):
    # Input parameter validation
    if callPut not in (0, 1):
        raise ValueError("callPut - must be a logical value of either 1 (call option) or 0 (put option).")
    elif S_max_min < 0:
        raise ValueError("S_max_min - must be positive.")
    elif S < 0:
        raise ValueError("S - must be positive.")
    elif T <= 0:
        raise ValueError("T must be a positive number.")
    elif r <= 0:
        raise ValueError("r must be a positive number.")
    elif D <= 0:
        raise ValueError("D must be a positive number.")
    elif sigma <= 0:
        raise ValueError("sigma must be a positive number.")
    elif r == D:
        raise ValueError("D must not be equal to r.")

    price = None
    term1 = None
    term2 = None
    brackCoeff = S * math.exp(-r * T) * (sigma * sigma) / (2 * (r - D))
    brack = None
    brack2 = None
    temp = None
    dt = sigma * math.sqrt(T)

    if callPut:
        a1 = (math.log(S / S_max_min) + (r - D + 0.5 * sigma * sigma) * T) / dt
        a2 = a1 - dt
        term1 = S * math.exp(-D * T) * norm.cdf(a1)
        term2 = S_max_min * math.exp(-r * T) * norm.cdf(a2)
        temp = -a1 + 2 * (r - D) * math.sqrt(T) / sigma
        brack = (S / S_max_min) ** (-2 * (r - D) / (sigma * sigma)) * norm.cdf(temp)
        brack2 = math.exp((r - D) * T) * norm.cdf(-a1)
        price = term1 - term2 + brackCoeff * (brack - brack2)
    else:
        b1 = (math.log(S_max_min / S) + (D + r + 0.5 * sigma * sigma) * T) / dt
        b2 = b1 - dt
        term1 = S_max_min * math.exp(-r * T) * norm.cdf(-b2)
        term2 = S * math.exp(-D * T) * norm.cdf(-b1)
        temp = b1 - 2 * (r - D) * math.sqrt(T) / sigma
        brack = (S / S_max_min) ** (2 * (r - D) / (sigma * sigma)) * norm.cdf(temp)
        brack2 = math.exp((r - D) * T) * norm.cdf(b1)
        price = term1 - term2 + brackCoeff * (-brack + brack2)

    return price

# Examples
# call = fltStrikeLookback(1, 100, 98, 0.501, 0.03, 0.1, 0.3)
# print(call)  # Output: 17.18954547977881

# put = fltStrikeLookback(0, 100, 98, 0.501, 0.03, 0.1, 0.3)
# print(put)  # Output: 14.04526346904223