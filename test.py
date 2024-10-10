import os
import time
import threading
from Crypto.Random import get_random_bytes
from concurrent.futures import ThreadPoolExecutor

alg_map = {'hyb': -1, 'sp': 0, 'ms': 1, 'bgf': 2, 'bf': 3}
num_map = {1: '1', 10: '10', 100: '100', 1000: '1k', 10000: '10k', 100000: '100k', 1000000: '1m', 10000000: '10m', 100000000: '100m', 48000000: '48m'}
dataset = {
    1: (142, 
        134, 
        [9349,9547,9749,9803,9859,9883,9901,9907,9923,9941,9949,10037,10067,10069,10091,10093,10099,10133,10139,10141,10181, 10253, 10259],
        [1295.8820352373061, 1263.4934919082689, 1231.6085621980212, 1223.2766681860812, 1214.7197728118047, 1211.0783381563747, 1208.3573528376025, 1207.452273222194, 1205.0433966126377, 1202.341505160659, 1201.1434074755352, 1188.0748427905337, 1183.665545762064, 1183.372415333471, 1180.1547420112815, 1179.8628397949878, 1178.9877453828299, 1174.0461722319012, 1173.177172071135, 1172.8877075586481, 1167.1195771087705, 1156.837693866321, 1155.9866687943572]),
    3: (206,
        199,
        [19139, 19141, 19157, 19163, 19181, 19219, 19237, 19259, 19301, 19333, 19373, 19379, 19387, 19403, 19427, 19469, 19483, 19501, 19507,19541],
        [2978.865075785719, 2978.5030368395746, 2975.6089189212425, 2974.524629336766, 2971.275044183069, 2964.4309508822867, 2961.196639692408, 2957.250234854245, 2949.736423729387, 2944.0293849634227, 2936.9171178328756, 2935.8523359193123, 2934.433460419369, 2931.5985650755138, 2927.353349989053, 2919.9447502974745, 2917.481007407138, 2914.317581801897, 2913.2641662363017, 2907.3048000408808]),
    0: (142,
        134,
        [12323],
        [909.0225315768487]),
    2: (142,
        134,
        [9349],
        [1295.8820352373061])
}


# -------P-A-N-E-L----------
ALG         = "bf"
T_TEST      = 10000        # number of tests in a thread
T           = 48000000     # number of total tests
W_DATASET   = 1
MAX_THREAD  = 200
NbIter      = 100
# --------------------------

s = "g++ -m64 -O3 ldpc_tests/test.c *.c ntl.cpp FromNIST/rng.c -I. -I/include -L/lib -std=c++11 -lcrypto -lssl -lm -ldl -lntl -lgmp -lgf2x -lpthread -DVERBOSE=0 -DNIST_RAND=1"
dir = f"{ALG}_data{W_DATASET}_{num_map[T]}"

def get_seed(len):
    seed = get_random_bytes(len)
    s = "{"
    for i in seed:
        s += f"(unsigned char){i},"
    return s[:-1] + "}"

def run_test(w, t, r, X, name):
    os.system(f"rm -f LDPC_test_{name}")
    para = f"-DR_BITS={r}ULL -DDV={w//2}ULL -DT1={t}ULL -DW_DECODER={alg_map[ALG]} -DNbIter={NbIter} -DTHR_X={X} -DT_TEST={T_TEST} -DINIT_SEED='{get_seed(48)}'"
    os.system(f"{s} {para} -o LDPC_test_{name}")
    os.system(f"./LDPC_test_{name} > {dir}/output_{name}")
    os.system(f"rm -f LDPC_test_{name}")

os.system(f"mkdir {dir}")

w, t, r, X = dataset[W_DATASET]

begin_t = time.time()
with ThreadPoolExecutor(max_workers = MAX_THREAD) as executor:
    for i in range(2, len(r)):
        for j in range(0, T // T_TEST):
            executor.submit(run_test, w, t, r[i], X[i], f"{ALG}_g{i}_p{j}")
end_t = time.time()

print(f"Running time: {end_t - begin_t}")
print(f"Number of tests per second: {T * len(r) / (end_t - begin_t)}")
