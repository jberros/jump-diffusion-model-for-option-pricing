# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:49:37 2023

@author: jeremy
"""

def IIFour(jj, ll, aa, bb, dd, nStep):
    IIFour = 0
    for k in range(1, nStep+1):
        IIFour += Table(II(k-1, aa - mu * bigT, eta2, 1/(sig*math.sqrt(bigT)), -(sig*math.sqrt(bigT))*eta2))
    return IIFour