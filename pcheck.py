from math import *
from scipy import integrate # type: ignore

alpha = 1.58
k = 1
T = 10000
X = 3

def f( dfr, x, T ):
    return (T + 1) * comb(T, x) * pow(dfr, x) * pow(1 - dfr, T - x)

def calc( x, T, alpha ):
    lower_bound = x / T * (1 / alpha)
    upper_bound = min(x / T * alpha, 1)
    res = integrate.quad(f, lower_bound, upper_bound, (x, T))
    return res

l = 1
r = 20
for i in range(100):
    mid = (l + r) / 2
    s = calc(X, T, mid)[0]
    if  s < 0.99:
        l = mid
    else:
        r = mid
print(l, r)