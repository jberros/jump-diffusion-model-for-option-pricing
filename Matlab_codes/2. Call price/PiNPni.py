# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:53:23 2023

@author: jeremy
"""

def PiNPni(nStep):
    PiNPni = 0
    for n in range(1, nStep+1):
        for k in range(1, n+1):
            PiNPni += Table(PiN(n) * Pni(n, k, p, eta1, eta2) * ((sig*math.sqrt(bigT)*eta1)**k))
    return PiNPni