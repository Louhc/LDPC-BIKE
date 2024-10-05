from math import *
from scipy import integrate # type: ignore

alpha = 1.58
k = 3
X = [22,16,10,5,2]
T = 10000

def f( dfr, x, T ):
    return (T + 1) * comb(T, x) * pow(dfr, x) * pow(1 - dfr, T - x)

def calc( x, T ):
    lower_bound = x / T * (1 / alpha)
    upper_bound = x / T * alpha
    res = integrate.quad(f, lower_bound, upper_bound, (x, T))
    return res

for i in range(len(X)):
    print(calc(X[i] * k, T * k))