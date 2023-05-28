# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:56:06 2023

@author: jeremy
"""

def sec(nStep):
    sec = 0
    for n in range(1, nStep+1):
        for k in range(1, n+1):
            sec += Sum(PiNPni(n, k) * IITwo(k))
    return sec