# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:05:54 2023

@author: jeremy
"""

from scipy.special import comb

def Pni(n, i, p, eta1, eta2):
    Pni = 0
    for j in range(i, n):
        Pni += comb(n, j) * p**j * (1-p)**(n-j) * comb(n-i-1, j-i) * \
                (eta1/(eta1 + eta2))**(j-i) * (eta2/(eta1 + eta2))**(n-j)
    
    Pni = p**n
    return Pni
