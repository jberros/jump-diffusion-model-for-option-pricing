# -*- coding: utf-8 -*-
"""
Created on Sun May 28 15:53:39 2023

@author: jeremy
"""

# You can use the BS function in Python by passing the required input parameters: 
# S0, t, K, T, Rgrow, Rdisc, and sigma. 
# The function will return the calculated call option price based on the Black-Scholes formula.

import math
from scipy.stats import norm

def BS(S0, t, K, T, Rgrow, Rdisc, sigma):
    F = S0 * math.exp(Rgrow * T)
    d1 = (math.log(F / K) / (sigma * math.sqrt(T - t))) + (sigma * math.sqrt(T) / 2)
    d2 = (math.log(F / K) / (sigma * math.sqrt(T - t))) - (sigma * math.sqrt(T) / 2)
    Call = math.exp(-Rdisc * T) * (F * norm.cdf(d1) - K * norm.cdf(d2))
    return Call


# Example usage: 
# call_price = BS(100, 0, 90, 1, 0.05, 0.03, 0.2)
# print(call_price)

# Please note that the Python code uses the norm.cdf function from the scipy.stats module 
# to calculate the cumulative distribution function of the standard normal distribution.