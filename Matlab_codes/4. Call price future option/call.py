# -*- coding: utf-8 -*-
"""
Created on Sun May 28 14:57:14 2023

@author: jeremy
"""

def call(eta1, eta2, la, p, sig, bond, bigF, bigK, bigT, nStep):
    def zetaa():
        return p * eta1 / (eta1 - 1) + (1 - p) * eta2 / (eta2 + 1) - 1

    def tempaa1():
        return sig**2 / 2 - la * zetaa()

    def tempaa2():
        return tempaa1() - sig**2

    call = bond * (bigF * cprob(tempaa1(), eta1 - 1, eta2 + 1, la * (1 + zetaa()), p * eta1 / ((1 + zetaa()) * (eta1 - 1)), sig, math.log(bigK / bigF), bigT, nStep) -
                   bigK * cprob(tempaa2(), eta1, eta2, la, p, sig, math.log(bigK / bigF), bigT, nStep))
    return call
