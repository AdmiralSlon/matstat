import random
import matplotlib.pyplot as plt
from math import log, e, sqrt, ceil
import numpy as np
from collections import Counter
from data import *

h = 2
ls, res = [], []
mx = 12


def Uniform():
    return random.randint(1, mx)

def Expon():
    x = random.uniform(0, 1)
    return (-log(1 - x)) / h

def distr_func_Exp(x):
    return 1 - (e) ** (-h * x)

def distr_func_Uni(x):
    return x / mx

def distr_func_Uni_Tetta(x):
    return x / max(UniSample)

def makeVariationVb(sz):
    vb = []
    for j in range(sz):
        vb.append(Uniform())
    return vb


UniSample.sort()
ExpSample.sort()

#хи-квадрат пирсона

def Chi_Square_Uni(k, l, r):
    n = 0
    chi = 0
    for j in np.arange(l, r, (r - l) / k):
        for i in UniSample:
            if(j < i <= j + (r - l) / k):
                n += 1
            if(i > r):
                break
        p = distr_func_Uni_Tetta(j + (r - l) / k) - distr_func_Uni_Tetta(j)
        chi += (n - len(UniSample) * p) ** 2 / (len(UniSample) * p)
        n = 0
    return chi
print(Chi_Square_Uni(11, 1, 12))


def Chi_Square_Exp(k, l, r):
    n = 0
    m = 0
    chi = 0
    for j in np.arange(l, r, (r - l) / k):
        for i in range(len(ExpSample)):
            if(j < ExpSample1[i] <= j + (r - l) / k):
                n += 1
            if(ExpSample1[i] > r):
                break
            if (j < ExpSample2[i] <= j + (r - l) / k):
                n += 1
            if (ExpSample2[i] > r):
                break
        chi += (n/ (len(ExpSample1)) - m/ (len(ExpSample1))) ** 2 / (n + m)
        n = 0
        m = 0
    return chi

print(Chi_Square_Exp(10, 0, 3))



#хи-квадрат для сложной гипотезы



def Chi_Square_Exp_Tetta(k, l, r):
    n = 0
    chi = 0
    for j in np.arange(l, r, (r - l) / k):
        for i in ExpSample_Tetta:
            if(j < i <= j + (r - l) / k):
                n += 1
            if(i > r):
                break
        p = distr_func_Exp_Tetta(j + (r - l) / k) - distr_func_Exp_Tetta(j)
        chi += (n - len(ExpSample_Tetta) * p) ** 2 / (len(ExpSample_Tetta) * p)
        n = 0
    return chi
print(Chi_Square_Exp_Tetta(10, 0, 3))



#колмогоров

Tetta = 1/(sum(ExpSample)/len(ExpSample))

def distr_func_Exp_Tetta(x):
    return 1 - (e) ** (-Tetta * x)

def Expon_Tetta():
    x = random.uniform(0, 1)
    return (-log(1 - x)) / Tetta

ExpSample_Tetta = []
for i in range(1000):
     ExpSample_Tetta.append(Expon_Tetta())
ExpSample_Tetta.sort()
m1 = 0
m2 = 0
for i in range(len(ExpSample_Tetta)):
    x = ExpSample_Tetta[i]
    if abs((i + 1) / len(ExpSample_Tetta) - distr_func_Exp_Tetta(x)) > m1:
        m1 = abs((i + 1) / len(ExpSample_Tetta) - distr_func_Exp_Tetta(x))
    if abs(distr_func_Exp_Tetta(x) - (i) / len(ExpSample_Tetta)) > m2:
        m2 = abs(distr_func_Exp_Tetta(x) - (i) / len(ExpSample_Tetta))
d = max(m1, m2)
print(d * sqrt(1000))


#колмогоров для сложной гипотезы




#смирнов однородность


def near_value(it, value):
    return min(it, key = lambda x: abs(x - value))

def func_Smir(sample,x):
    m = 0
    temp = near_value(sample, x)
    if temp >= x:
         m = (sample.index(temp) + 1)/len(sample)
    else:
        m = (sample.index(temp))/len(sample)
    return m

sample1 = []
sample2 = []
for i in range(1000):
    sample1.append(Expon())
    sample2.append(Expon())
sample1.sort()
sample2.sort()

D1 = 0
D2 = 0
for r in range(1000):
    if D1 < ((r + 1) / 1000) - func_Smir(sample1, sample2[r]):
        D1 = ((r + 1) / 1000) - func_Smir(sample1, sample2[r])
    if D2 < (func_Smir(sample1,sample2[r]) - r / 1000):
        D2 = (func_Smir(sample1, sample2[r]) - r / 1000)

print(max(D1, D2))

print(sqrt(1.3581 * (1 / 1000 + 1 / 1000)))
хи-квадрат

sample1 = []
sample2 = []
for i in range(1000):
    sample1.append(Expon())
    sample2.append(Expon())
sample1.sort()
sample2.sort()

def Chi_Square_ODN(k, l, r):
    n = 0
    m = 0
    chi = 0
    for j in np.arange(l, r, (r - l) / k):
        for i in range(len(ExpSample)):
            if(j < sample1[i] <= j + (r - l) / k):
                n += 1
            if(sample1[i] > r):
                break
            if (j < sample2[i] <= j + (r - l) / k):
                m += 1
            if (sample2[i] > r):
                break
        chi += ((n/ (len(sample1)) - m / (len(sample1))) ** 2) / (n + m)
        n = 0
        m = 0
    return (len(sample1)) * (len(sample2)) * chi

print(Chi_Square_ODN(23, 0, 2.5))


