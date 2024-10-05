import os
import time
import threading
from Crypto.Random import get_random_bytes
from concurrent.futures import ThreadPoolExecutor

dir = "bgf_data_withrandominitial_10k"
os.system(f"mkdir {dir}")

def get_seed( len ):
    seed = get_random_bytes(len)
    s = "{"
    for i in seed:
        s += f"(unsigned char){i},"
    return s[:-1] + "}"

# W_DECODER = 0  :  SP
# W_DECODER = 1  :  MS
# W_DECODER = 2  :  BGF

s = "g++ -m64 -O3 ldpc_tests/test_singlethread.c *.c ntl.cpp FromNIST/rng.c -I. -I/include -L/lib -std=c++11 -lcrypto -lssl -lm -ldl -lntl -lgmp -lgf2x -lpthread -DVERBOSE=0 -DNIST_RAND=1"

def run_test( w, t, r, X, name ):
    os.system(f"rm -f LDPC_test_{name}")
    para = f"-DR_BITS={r}ULL -DDV={w//2}ULL -DT1={t}ULL -DW_DECODER={2} -DTHR_X={X} -DINIT_SEED='{get_seed(48)}'"
    os.system(f"{s} {para} -o LDPC_test_{name}")
    os.system(f"./LDPC_test_{name} > {dir}/output_{name}")
    os.system(f"rm -f LDPC_test_{name}")

w = 142  # the hamming weight of one row in G
t = 134  # the hamming weight of the error vector
r = [9349,9547,9749,9803,9859,9883,9901,9907,9923,9941,9949,10037,10067,10069,10091,10093,10099,10133,10139,10141,10181, 10253, 10259]
X = [0.13861183391136014, 0.13234455765248448, 0.12633178399815584, 0.12478595003428355, 0.12320922738734197, 0.12254157018682327, 0.12204397059262725, 0.12187869922501204, 0.1214394232200582, 0.12094774219501649, 0.12073006407433258, 0.11836951706590952, 0.11757877677183516, 0.11752630999438585, 0.11695121811627009, 0.11689912214356364, 0.11674301865361224, 0.11586363093179719, 0.1157093571428282, 0.11565799305380617, 0.11463702751289369, 0.11282919085792652, 0.11268024844471754]

assert len(r) == len(X)

Ntests = 1000
T = 10000
th = []

# for i in range(0, len(r)):
#     for j in range(0, T // Ntests):
#         th.append(threading.Thread(target = run_test, args=(w, t, r[i], X[i], f"bgf_g{i}_p{j}")))

# for i in range(len(th)):
#     print(i)
#     th[i].start()

# for i in range(len(th)):
#     th[i].join()

with ThreadPoolExecutor(max_workers = 128) as executor:
    for i in range(0, len(r)):
        for j in range(0, T // Ntests):
            executor.submit(run_test, w, t, r[i], X[i], f"bgf_g{i}_p{j}")