import os
from Crypto.Random import get_random_bytes
from concurrent.futures import ThreadPoolExecutor

dir = "sp_data1_withrandominitial_10k"
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

def run_test( w, t, r, name ):
    os.system(f"rm -f LDPC_test_{name}")
    para = f"-DR_BITS={r}ULL -DDV={w//2}ULL -DT1={t}ULL -DW_DECODER={0} -DINIT_SEED='{get_seed(48)}'"
    os.system(f"{s} {para} -o LDPC_test_{name}")
    os.system(f"./LDPC_test_{name} > {dir}/output_{name}")
    os.system(f"rm -f LDPC_test_{name}")

w = 142  # the hamming weight of one row in G
t = 134  # the hamming weight of the error vector
r = [9349,9547,9749,9803,9859,9883,9901,9907,9923,9941,9949,10037,10067,10069,10091,10093,10099,10133,10139,10141,10181, 10253, 10259]

Ntests = 1000
T = 10000
th = []

with ThreadPoolExecutor(max_workers = 128) as executor:
    for i in range(0, len(r)):
        for j in range(0, T // Ntests):
            executor.submit(run_test, w, t, r[i], f"sp_g{i}_p{j}")