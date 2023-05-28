# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:50:56 2023

@author: jeremy
"""

def IITwo(jj, ll, aa, bb, dd, nStep):
    IITwo = 0
    for k in range(1, nStep+1):
        IITwo += Table(II(k-1, aa - mu * bigT, -eta1, -1/(sig*math.sqrt(bigT)), -(sig*math.sqrt(bigT))*eta1))
    return IITwo