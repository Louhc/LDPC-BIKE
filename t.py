from math import *

X = 0.07376633381293914
r = 12323
t = 134
w = 142
d = w // 2
n = r * 2


def VAR_TH_FCT( x ):
    pi1 = (x + X) / t / d
    pi0 = (d * 2 * x - X) / (n - t) / d
    T = ceil((log((n - t) / t) + d * log((1 - pi0) / (1 - pi1))) / (log(pi1 / pi0) + log((1 - pi0) / (1 - pi1))))
    return max(T, (d + 1) / 2.0)

def f( x ):
    return max(13.530 + 0.0069722 * (x), 36)

for i in range(4500, 5000):
    if abs(VAR_TH_FCT(i) - f(i)) >=1:
        print("{}  {}".format(VAR_TH_FCT(i), f(i)))