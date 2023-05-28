# -*- coding: utf-8 -*-
"""
Created on Sun May 28 13:58:21 2023

@author: jeremy
"""

import math
from scipy.special import gamma, hypergeom
from scipy.integrate import quad

def temp(n, x):
    return x + math.sqrt(x**2 + 4*n) / 2

def Hh(n, x):
    def integrand(t, n, x):
        return (t - x)**n * math.exp(-t**2 / 2)

    if -6 <= x < 10:
        result, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), x, math.inf)
        Hh = 1 / math.factorial(n) * result
    else:
        temp_val = temp(n, x)
        result1, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), x, temp_val - 3)
        result2, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val - 3, temp_val - 1)
        result3, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val - 1, temp_val)
        result4, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val, temp_val + 1)
        result5, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val + 1, temp_val + 3)
        result6, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val + 3, math.inf)
        Hh = (1 / math.factorial(n)) * (result1 + result2 + result3 + result4 + result5 + result6)

    return Hh
