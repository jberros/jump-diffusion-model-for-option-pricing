#!/usr/bin/env python
# coding: utf-8

# In[50]:


"""
Created on Sun May 28 15:13:56 2023

@author: jeremy
"""


# In[51]:


import math
import pandas as pd
# import scipy.special as sp
# from scipy.special import gamma, hypergeom
from scipy.integrate import quad


# In[52]:


def phi(x):
    if -10 < x < 10:
        return (1 + math.erf(x / math.sqrt(2))) / 2
    else:
        return 0


# In[53]:


def temp(n, x):
    return x + math.sqrt(x**2 + 4*n) / 2


# In[54]:


def Hh(n, x):
    def integrand(t, n, x):
        return (t - x)**n * math.exp(-t**2 / 2)

    if -6 <= x < 10:
        result, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), x, math.inf)
        Hh = 1 / math.factorial(n) * result
    else:
        temp_val = temp(n, x)
        result1, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), x, temp_val - 3)
        result2, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val - 3, temp_val - 1)
        result3, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val - 1, temp_val)
        result4, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val, temp_val + 1)
        result5, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val + 1, temp_val + 3)
        result6, _ = quad(lambda t: (t - x)**n * math.exp(-t**2 / 2), temp_val + 3, math.inf)
        Hh = (1 / math.factorial(n)) * (result1 + result2 + result3 + result4 + result5 + result6)

    return Hh


# In[55]:


def Ha(m, y):
    if m == -1:
        return math.exp(-y**2/2)
    elif m == 0:
        return math.sqrt(2*math.pi)*phi(-y)
    else:
        return (-y/m)*Ha(m - 1, y) + (1/m)*Ha(m - 2, y)


# In[56]:


def II(jj, ll, aa, bb, dd):
    if bb > 0 and aa != 0:
        table1 = {i: (bb/aa)**(jj-i) for i in range(jj+1)}
        table2 = {i: Hh(i, bb*ll-dd) for i in range(jj+1)}
        II = -(math.exp(aa*ll)/aa) * sum(table1[i] * table2[i] for i in range(jj+1)) +             ((bb/aa)**(jj+1)) * (math.sqrt(2*math.pi)/bb) * math.exp(aa*dd/bb + (1/2) * (aa/bb)**2) * phi(-bb*ll + dd + aa/bb)
    elif bb < 0 and aa < 0:
        table1 = {i: (bb/aa)**(jj-i) for i in range(jj+1)}
        table2 = {i: Hh(i, bb*ll-dd) for i in range(jj+1)}
        II = -(math.exp(aa*ll)/aa) * sum(table1[i] * table2[i] for i in range(jj+1)) -             ((bb/aa)**(jj+1)) * (math.sqrt(2*math.pi)/bb) * math.exp(aa*dd/bb + (1/2) * (aa/bb)**2) * phi(bb*ll - dd - aa/bb)
    elif bb > 0 and aa == 0:
        II = Hh(jj + 1, bb*ll - dd)/bb
    return II


# In[57]:


def Pni(n, i, p, eta1, eta2):
    Pni = 0
    for j in range(i, n):
        Pni += math.comb(n, j) * (p ** j) * ((1 - p) ** (n - j)) * math.comb(n - i - 1, j - i) *                ((eta1 / (eta1 + eta2)) ** (j - i)) * ((eta2 / (eta1 + eta2)) ** (n - j))
    return Pni


# In[58]:


def Qni(n, i, p, eta1, eta2):
    Qni = 0
    for j in range(i, n):
        Qni += math.comb(n, j) * ((1 - p) ** j) * (p ** (n - j)) * math.comb(n - i - 1, j - i) *                ((eta2 / (eta1 + eta2)) ** (j - i)) * ((eta1 / (eta1 + eta2)) ** (n - j))
    return (1 - p) ** n


# In[79]:


def cprob(mu, eta1, eta2, la, p, sig, aa, bigT, nStep):
    def IITwo(k):
        IITwo = 0
        for k in range(1, nStep + 1):
            IITwo += II(k-1, aa - mu * bigT, -eta1, -1/(sig*math.sqrt(bigT)), -(sig*math.sqrt(bigT))*eta1)
        return IITwo

    def IIFour(k):
        IIFour = 0
        for k in range(1, nStep+1):
            IIFour += II(k-1, aa - mu * bigT, eta2, 1/(sig*math.sqrt(bigT)), -(sig*math.sqrt(bigT))*eta2)
        return IIFour

    def PiN(n):
        return math.exp(-la*bigT) * ((la*bigT)**n) / math.factorial(n)

    def PiNPni():
        result = 0
        for n in range(1, nStep+1):
            for k in range(1, n+1):
                result += PiN(n) * Pni(n, k, p, eta1, eta2) * ((sig*math.sqrt(bigT)*eta1)**k)
        return result

    def PiNQni():
        result = 0
        for n in range(1, nStep+1):
            for k in range(1, n+1):
                result += PiN(n) * Qni(n, k, p, eta1, eta2) * ((sig*math.sqrt(bigT)*eta2)**k)
        return result

    sec = 0
    for n in range(1, nStep+1):
        for k in range(1, n+1):
            sec += PiNPni() * IITwo(k)

    fourth = 0
    for n in range(1, nStep+1):
        for k in range(1, n+1):
            fourth += PiNQni() * IIFour(k)


    cprob = (sec * math.exp(((sig*eta1)**2)*bigT/2) + fourth * math.exp(((sig*eta2)**2)*bigT/2)) / (math.sqrt(2*math.pi) * sig * math.sqrt(bigT)) +             math.exp(-la*bigT) * phi(-(aa- mu*bigT)/(sig*math.sqrt(bigT)))

    return cprob


# In[80]:


def callOR(eta1, eta2, la, p, sig, rr, bigS, bigK, bigT, nStep):
    def zetaaOR():
        return p * eta1 / (eta1 - 1) + (1 - p) * eta2 / (eta2 + 1) - 1
    
    def tempaa1OR():
        return rr + sig * sig / 2 - la * zetaaOR()
    
    def tempaa2OR():
        return tempaa1OR() - sig * sig
    
    callOR = bigS * cprob(tempaa1OR(), 
                          eta1 - 1, 
                          eta2 + 1, 
                          la * (1 + zetaaOR()),
                          p * eta1 / ((1 + zetaaOR()) * (eta1 - 1)), 
                          sig, 
                          math.log(bigK / bigS), 
                          bigT, 
                          nStep) 
    - bigK * math.exp(-rr * bigT) * cprob(tempaa2OR(), 
                                          eta1, 
                                          eta2, 
                                          la, 
                                          p, 
                                          sig, 
                                          math.log(bigK / bigS), 
                                          bigT, 
                                          nStep)
    
    return callOR


# In[81]:


def call(eta1, eta2, la, p, sig, bond, bigF, bigK, bigT, nStep):
    def zetaa():
        return p * eta1 / (eta1 - 1) + (1 - p) * eta2 / (eta2 + 1) - 1
    
    def tempaa1():
        return sig * sig / 2 - la * zetaa()
    
    def tempaa2():
        return tempaa1() - sig * sig
    
    call = bond * (bigF * cprob(tempaa1(), 
                                eta1 - 1, 
                                eta2 + 1, 
                                la * (1 + zetaa()),
                                p * eta1 / ((1 + zetaa()) * (eta1 - 1)), 
                                sig, 
                                math.log(bigK / bigF), 
                                bigT, 
                                nStep) 
                   - bigK * cprob(tempaa2(), 
                                  eta1, 
                                  eta2, 
                                  la, 
                                  p, 
                                  sig, 
                                  math.log(bigK / bigF), 
                                  bigT, 
                                  nStep))
    
    return call


# In[83]:


call(10, 5, 1, 0.4, 0.16, 0.05, 100, 98, 0.5, 20 )


# In[ ]:




