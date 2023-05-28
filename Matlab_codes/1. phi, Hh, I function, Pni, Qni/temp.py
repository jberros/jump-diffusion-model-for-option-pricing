# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:08:11 2023

@author: jeremy
"""

import math

def temp(n, x):
    temp = (x + math.sqrt(x**2 + 4*n) / 2)
    return temp