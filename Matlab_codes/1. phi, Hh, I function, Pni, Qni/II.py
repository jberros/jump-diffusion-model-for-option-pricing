# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:00:33 2023

@author: jeremy
"""

import math
from scipy.special import gamma

def Table(func, values):
    return [func(i) for i in values]

def Hh(m, y):
    if m == -1:
        return math.exp(-y**2/2)
    elif m == 0:
        return math.sqrt(2*math.pi)*phi(-y)
    else:
        return (-y/m)*Hh(m - 1, y) + (1/m)*Hh(m - 2, y)

def phi(x):
    # Implement the phi function according to your requirements
    # This is a placeholder, replace it with the actual implementation
    return 0

def II(jj, ll, aa, bb, dd):
    if bb > 0 and aa != 0:
        II = -(math.exp(aa * ll) / aa) * sum([Table(lambda i: (bb / aa) ** (jj - i), range(jj + 1))] *
                                             [Table(lambda i: Hh(i, bb * ll - dd), range(jj + 1))]) + \
             ((bb / aa) ** (jj + 1)) * (math.sqrt(2 * math.pi) / bb) * math.exp(aa * dd / bb + (1 / 2) * (aa / bb) ** 2) * \
             phi(-bb * ll + dd + aa / bb)

    elif bb < 0 and aa < 0:
        II = -(math.exp(aa * ll) / aa) * sum([Table(lambda i: (bb / aa) ** (jj - i), range(jj + 1))] *
                                             [Table(lambda i: Hh(i, bb * ll - dd), range(jj + 1))]) - \
             ((bb / aa) ** (jj + 1)) * (math.sqrt(2 * math.pi) / bb) * math.exp(aa * dd / bb + (1 / 2) * (aa / bb) ** 2) * \
             phi(bb * ll - dd - aa / bb)

    elif bb > 0 and aa == 0:
        II = Hh(jj + 1, bb * ll - dd) / bb

    return II
