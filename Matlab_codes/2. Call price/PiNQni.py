# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:54:54 2023

@author: jeremy
"""

def PiNQni(nStep):
    PiNQni = 0
    for n in range(1, nStep+1):
        for k in range(1, n+1):
            PiNQni += Table(PiN(n) * Qni(n, k, p, eta1, eta2) * ((sig*math.sqrt(bigT)*eta2)**k))
    return PiNQni
