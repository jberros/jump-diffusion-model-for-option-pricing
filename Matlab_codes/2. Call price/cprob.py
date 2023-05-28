# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:13:26 2023

@author: jeremy
"""

from scipy.special import factorial
from scipy.stats import norm
import math

def cprob(mu, eta1, eta2, la, p, sig, aa, bigT, nStep):
    def IITwo(jj, ll, aa, bb, dd):
        IITwo = 0
        for k in range(1, nStep+1):
            IITwo += II(k-1, aa - mu * bigT, -eta1, -1/(sig*math.sqrt(bigT)), -(sig*math.sqrt(bigT))*eta1)
        return IITwo

    def IIFour(jj, ll, aa, bb, dd):
        IIFour = 0
        for k in range(1, nStep+1):
            IIFour += II(k-1, aa - mu * bigT, eta2, 1/(sig*math.sqrt(bigT)), -(sig*math.sqrt(bigT))*eta2)
        return IIFour

    def PiN(n):
        return math.exp(-la*bigT) * ((la*bigT)**n) / factorial(n)

    def PiNPni():
        PiNPni = 0
        for n in range(1, nStep+1):
            for k in range(1, n+1):
                PiNPni += PiN(n) * Pni(n, k, p, eta1, eta2) * ((sig*math.sqrt(bigT)*eta1)**k)
        return PiNPni

    def PiNQni():
        PiNQni = 0
        for n in range(1, nStep+1):
            for k in range(1, n+1):
                PiNQni += PiN(n) * Qni(n, k, p, eta1, eta2) * ((sig*math.sqrt(bigT)*eta2)**k)
        return PiNQni

    def sec():
        sec = 0
        for n in range(1, nStep+1):
            for k in range(1, n+1):
                sec += PiNPni((n,k)) * IITwo((k))
        return sec

    def fourth():
        fourth = 0
        for n in range(1, nStep+1):
            for k in range(1, n+1):
                fourth += PiNQni((n, k)) * IIFour((k))
        return fourth

    cprob = (sec() * math.exp(((sig*eta1)**2)*bigT/2) + fourth() * math.exp(((sig*eta2)**2)*bigT/2)) / \
            (math.sqrt(2*math.pi) * sig * math.sqrt(bigT)) + math.exp(-la*bigT) * norm.cdf(-(aa- mu*bigT)/(sig*math.sqrt(bigT)))
    return cprob
