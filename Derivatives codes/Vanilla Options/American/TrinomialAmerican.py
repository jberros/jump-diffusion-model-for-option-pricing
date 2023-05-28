# -*- coding: utf-8 -*-
"""
Created on Sun May 28 16:26:18 2023

@author: jeremy
"""

# trinomial_american function calculates the price of the American option 
# using the Boyle trinomial tree method. 
# The function returns the price of the option. 

# % Computes the Boyle (1986) Trinomial Tree for American Call/Put Option Values based
# % on the following inputs:
# % CallPut           =       Call = 1, Put = 0
# % AssetP            =       Underlying Asset Price
# % Strike            =       Strike Price of Option
# % RiskFree          =       Risk Free rate of interest
# % Div               =       Dividend Yield of Underlying
# % Time              =       Time to Maturity
# % Vol               =       Volatility of the Underlying
# % nSteps            =       Number of Time Steps for Trinomial Tree to take
# % Please note that the use of this code is not restricted in anyway.
# % However, referencing the author of the code would be appreciated.
# % To run this program, simply use the function defined in the 1st line.
# % http://www.global-derivatives.com 
# % info@global-derivatives.com
# % Kevin Cheng (Nov 2003)

import math

def trinomial_american(CallPut, AssetP, Strike, RiskFree, Div, Time, Vol, nSteps):
    dt = Time / nSteps
    cc = RiskFree - Div

    if CallPut:
        b = 1
    else:
        b = -1

    RR = math.exp(RiskFree * dt)
    Up = math.exp(Vol * math.sqrt(2 * dt))
    Down = 1 / Up
    P_up = ((math.exp(cc * dt / 2) - math.exp(-Vol * math.sqrt(dt / 2))) / (math.exp(Vol * math.sqrt(dt / 2)) - math.exp(-Vol * math.sqrt(dt / 2)))) ** 2
    P_down = ((math.exp(Vol * math.sqrt(dt / 2)) - math.exp(cc * dt / 2)) / (math.exp(Vol * math.sqrt(dt / 2)) - math.exp(-Vol * math.sqrt(dt / 2)))) ** 2
    P_mid = 1 - P_up - P_down
    Df = math.exp(-RiskFree * dt)

    # Sets up the asset movements on the trinomial tree
    Value = [0] * (2 * nSteps + 1)
    for i in range(2 * nSteps + 1):
        State = i + 1
        Value[State] = max(0, b * (AssetP * pow(Up, max(i - nSteps, 0)) * pow(Down, max(nSteps * 2 - nSteps - i, 0)) - Strike))

    # Works backwards recursively to determine the price of the option
    for TT in range(nSteps - 1, -1, -1):
        for i in range(TT * 2 + 1):
            State = i + 1
            Value[State] = max((P_up * Value[State + 2] + P_mid * Value[State + 1] + P_down * Value[State]) * Df,
                              b * (AssetP * pow(Up, max(i - TT, 0)) * pow(Down, max(TT * 2 - TT - i, 0)) - Strike))

    Trinomial = Value[1]
    return Trinomial

# Example usage
CallPut = 1  # 1 for Call, 0 for Put
AssetP = 100
Strike = 105
RiskFree = 0.05
Div = 0.02
Time = 1
Vol = 0.2
nSteps = 100

option_price = trinomial_american(CallPut, AssetP, Strike, RiskFree, Div, Time, Vol, nSteps)
print(option_price)
