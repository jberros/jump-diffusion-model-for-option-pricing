# -*- coding: utf-8 -*-
"""
Created on Sun May 28 15:05:36 2023

@author: jeremy
"""

def zetaaOR(p, eta1, eta2):
    zetaaOR = p * eta1 / (eta1 - 1) + (1 - p) * eta2 / (eta2 + 1) - 1
    return zetaaOR