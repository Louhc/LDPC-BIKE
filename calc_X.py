from math import *
from functools import lru_cache

w = 142
t = 134
rr = [8689,8707,8731,8747,8779,8803,8821,8837,8861,8887,8923,8951,9091,9109]
X = []

@lru_cache(None)  # Cache the results of binomial coefficient calculations
def cached_comb(n, k):
    return comb(n, k)

def compute_sum_odd_l_optimized(n, t, w):
    X = 0
    total_comb_nt = cached_comb(n, t)  # Compute this once since it's used in every term
    # Iterate over odd l from 1 to t
    for l in range(1, min(t, w) + 1, 2):
        if l > w: break  # If l is greater than w, comb(w, l) will be 0
        numerator = (l - 1) * (n // 2) * cached_comb(w, l) * cached_comb(n - w, t - l)
        term = numerator / total_comb_nt
        X += term
    return X

def VAR_TH_FCT( x, THR_X, T1, DV, N_BITS ):
    pi1 = (x + THR_X) / T1 / DV
    pi0 = (DV * 2 * x - THR_X) / (N_BITS - T1) / DV
    # assert(pi1 >= pi0)
    T = ceil((log((N_BITS - T1) / T1) + DV * log((1 - pi0) / (1 - pi1))) / (log(pi1 / pi0) + log((1 - pi0) / (1 - pi1))))
    return max(T, (DV + 1) / 2.0)

def f( x ):
    return max(13.530 + 0.0069722 * (x), 36)


for r in rr:
    X.append(compute_sum_odd_l_optimized(r * 2, t, w))

print(X)

# r = 9547
# X = compute_sum_odd_l_optimized(r * 2, t, w)

# for i in range(4500, 5000):
#     print(f(i), VAR_TH_FCT(i, X, t, w // 2, r * 2))