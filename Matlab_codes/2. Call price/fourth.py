# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:45:24 2023

@author: jeremy
"""

def fourth(nStep):
    fourth = 0
    for n in range(1, nStep+1):
        for k in range(1, n+1):
            fourth += PiNQni((n, k)) * IIFour((k))
    return fourth