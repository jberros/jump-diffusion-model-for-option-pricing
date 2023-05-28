# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:02:29 2023

@author: jeremy
"""

from scipy.special import erf
import math

def phi(x):
    if -10 < x < 10:
        phi = (1 + erf(x / math.sqrt(2))) / 2
    else:
        phi = 0
    return phi
