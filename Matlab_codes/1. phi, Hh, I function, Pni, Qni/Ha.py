# -*- coding: utf-8 -*-
"""
Created on Sun May 28 13:57:08 2023

@author: jeremy
"""

import math

def Ha(m, y):
    if m == -1:
        return math.exp(-y**2/2)
    elif m == 0:
        return math.sqrt(2*math.pi)*phi(-y)
    else:
        return (-y/m)*Ha(m - 1, y) + (1/m)*Ha(m - 2, y)
