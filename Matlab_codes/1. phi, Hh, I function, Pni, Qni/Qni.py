# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:06:44 2023

@author: jeremy
"""

from scipy.special import comb

def Qni(n, i, p, eta1, eta2):
    Qni = 0
    for j in range(i, n):
        Qni += comb(n, j) * (1-p)**j * p**(n-j) * comb(n-i-1, j-i) * \
                (eta2 / (eta1 + eta2))**(j-i) * (eta1 / (eta1 + eta2))**(n-j)
    
    Qni = (1-p)**n
    return Qni