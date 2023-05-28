# -*- coding: utf-8 -*-
"""
Created on Sun May 28 15:57:48 2023

@author: jeremy
"""

# You can use the JDimpv function in Python by passing the required input parameters: 
# S0, X, r, T, a, b, lambda, and value. 

# The function will return the implied volatility of the underlying asset 
# based on the provided European option prices.

# Please note that the brentq function from the scipy.optimize module is used 
# to solve the nonlinear equation in the objective function. 

# The function JDprice is not provided, so you will need to define it separately in your code 
# to calculate the call value based on the input parameters.

import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson
from scipy.stats import norm
import math
from mibian import BS

def blsprice(S0, X, r, T, vol):
    option = BS([S0, X, r, T], volatility=vol)
    return option.callPrice

def BS1(S0, t, K, T, Rgrow, Rdisc, sigma):
    F = S0 * math.exp(Rgrow * T)
    d1 = (math.log(F / K) / (sigma * math.sqrt(T - t))) + (sigma * math.sqrt(T) / 2)
    d2 = (math.log(F / K) / (sigma * math.sqrt(T - t))) - (sigma * math.sqrt(T) / 2)
    Call = math.exp(-Rdisc * T) * (F * norm.cdf(d1) - K * norm.cdf(d2))
    return Call

def JDprice(S0, X, r, T, vol, a, b, _lambda):
    m = 100  # Number of Monte-Carlo simulations
    if _lambda == 0:
        JDCallPrice = blsprice(S0, X, r, T, vol)
        std_err = 0
    else:
        alpha = 1e-6
        Jav = (np.exp(b) - np.exp(a)) / (b - a) - 1
        K = poisson.ppf(1 - alpha, _lambda * T)
        p0 = poisson.pmf(0, _lambda * T)
        scaling = poisson.cdf(K, _lambda * T)
        BS_k0 = BS1(S0 * np.exp(-_lambda * Jav * T), 0, X, T, r, r, vol)  # k=0, no jumps
        p = poisson.pmf(np.arange(1, K + 1), _lambda * T)
        U = np.random.rand(K, m)
        Sk = np.zeros((K, m))
        Sk_a = np.zeros((K, m))
        S0_poiss = np.zeros((K, m))
        S0_poiss_a = np.zeros((K, m))
        BS_k = np.zeros((K, m))
        BS_k_a = np.zeros((K, m))
        for k in range(1, K + 1):
            Sk[k - 1, :] = k * a + (b - a) * np.sum(U[0:k, :], axis=0)
            Sk_a[k - 1, :] = (a + b) * k - Sk[k - 1, :]
            S0_poiss[k - 1, :] = S0 * np.exp(Sk[k - 1, :] - _lambda * Jav * T)
            S0_poiss_a[k - 1, :] = S0 * np.exp(Sk_a[k - 1, :] - _lambda * Jav * T)
            BS_k[k - 1, :] = BS1(S0_poiss[k - 1, :], 0, X, T, r, r, vol)
            BS_k_a[k - 1, :] = BS1(S0_poiss_a[k - 1, :], 0, X, T, r, r, vol)
        sample = (p0 * BS_k0 + np.dot(p, BS_k)) / scaling
        sample_a = (p0 * BS_k0 + np.dot(p, BS_k_a)) / scaling
        x = 0.5 * (sample + sample_a)
        y = 0.5 * (p * np.exp(Sk) + p * np.exp(Sk_a))
        VARy = 0.5 * (np.exp(_lambda * T * Jav) - 2 * np.exp(2 * _lambda * T * Jav) + np.exp(_lambda * T * (np.exp(a + b) - 1)))
        beta = m / (m - 1) * (np.mean(x * y) - np.mean(x) * np.mean(y)) / VARy
        Z = x - beta * (y - np.exp(_lambda * T * Jav))
        y2 = y * (2 * np.mean(y) - y)
        bias = 1 / (m - 1) * (np.mean(x * y2) - np.mean(x) * np.mean(y2)) / VARy
        JDCallPrice = np.mean(Z) - bias
        std_err = np.std(Z) / np.sqrt(m)
    
    return JDCallPrice, std_err

def JDimpv(S0, X, r, T, a, b, _lambda, value):
    def objfcn(volatility, S0, X, r, T, a, b, _lambda, value):
        callValue = JDprice(S0, X, r, T, volatility, a, b, _lambda)
        delta = value - callValue
        return delta

    volatility = np.empty_like(value)
    
    for i in range(value.size):
        if S0 > 0 and X > 0 and T > 0:
            try:
                result = brentq(objfcn, 0, 9, args=(S0, X, r, T, a, b, _lambda, value[i]))
                volatility[i] = result
            except:
                volatility[i] = np.nan
    
    volatility = volatility.reshape(value.shape)
    return volatility
