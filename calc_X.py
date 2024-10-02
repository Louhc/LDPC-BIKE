from math import comb
from functools import lru_cache

w = 142
t = 134
rr = [11239, 11273, 11317, 11369, 11423, 11471, 11503, 11579, 11621, 11689, 11731, 11789, 11827, 11867, 11923, 11953, 11987, 12043, 12101, 12143, 12197, 12241, 12277]
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
        numerator = (l - 1) * cached_comb(w, l) * cached_comb(n - w, t - l)
        term = numerator / total_comb_nt
        X += term
    return X

for r in rr:
    X.append(compute_sum_odd_l_optimized(r * 2, t, w))

print(X)
