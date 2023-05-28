# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:52:23 2023

@author: jeremy
"""

import math

def PiN(n, la, bigT):
    PiN = math.exp(-la*bigT) * ((la*bigT)**n) / math.factorial(n)
    return PiN