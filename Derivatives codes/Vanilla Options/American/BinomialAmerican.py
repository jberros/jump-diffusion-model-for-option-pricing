# -*- coding: utf-8 -*-
"""
Created on Sun May 28 16:18:18 2023

@author: jeremy
"""

# binomial_american function calculates the price of the American option 
# using the Cox, Ross & Rubinstein binomial tree method. 
# The function returns the price of the option. 

# Computes the Cox, Ross & Rubinstein (1979) Binomial Tree for American Call/Put Option Values based
# on the following inputs:
# CallPut           =       Call = 1, Put = 0
# AssetP            =       Underlying Asset Price
# Strike            =       Strike Price of Option
# RiskFree          =       Risk Free rate of interest
# Div               =       Dividend Yield of Underlying
# Time              =       Time to Maturity
# Vol               =       Volatility of the Underlying
# nSteps            =       Number of Time Steps for Binomial Tree to take
# Please note that the use of this code is not restricted in anyway.
# However, referencing the author of the code would be appreciated.
# To run this program, simply use the function defined in the 1st line.
# http://www.global-derivatives.com 
# info@global-derivatives.com
# Kevin Cheng (Nov 2003)

import math

def binomial_american(CallPut, AssetP, Strike, RiskFree, Div, Time, Vol, nSteps):
    dt = Time / nSteps

    if CallPut:
        b = 1
    else:
        b = -1

    RR = pow(math.e, RiskFree * dt)
    Up = pow(math.e, Vol * math.sqrt(dt))
    Down = 1 / Up
    P_up = (pow(math.e, (RiskFree - Div) * dt) - Down) / (Up - Down)
    P_down = 1 - P_up
    Df = pow(math.e, -RiskFree * dt)

    # Sets up the asset movements on the binomial tree
    Value = [0] * (nSteps + 1)
    for i in range(nSteps + 1):
        State = i + 1
        St = AssetP * pow(Up, i) * pow(Down, nSteps - i)
        Value[State] = max(0, b * (St - Strike))

    # Works backwards recursively to determine the price of the option
    for TT in range(nSteps - 1, -1, -1):
        for i in range(TT + 1):
            State = i + 1
            Value[State] = max((b * (AssetP * pow(Up, i) * pow(Down, abs(i - TT)) - Strike)),
                              (P_up * Value[State + 1] + P_down * Value[State]) * Df)

    Binomial = Value[1]
    return Binomial

# Example usage
# CallPut = 1  # 1 for Call, 0 for Put
# AssetP = 100
# Strike = 105
# RiskFree = 0.05
# Div = 0.02
# Time = 1
# Vol = 0.2
# nSteps = 100

# option_price = binomial_american(CallPut, AssetP, Strike, RiskFree, Div, Time, Vol, nSteps)
# print(option_price)
