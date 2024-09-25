import os

w = 142
t = 134
r = [9349,9547,9749,9803,9859,9883,9901,9907,9923,9941,9949,10037,10067,10069,10091,10093,10099,10133,10139,10141,10181, 10253, 10259]

s = "g++ -m64 -O3 ldpc_tests/test.c *.c ntl.cpp FromNIST/rng.c -I. -I/include -L/lib -std=c++11 -lcrypto -lssl -lm -ldl -lntl -lgmp -lgf2x -lpthread -DVERBOSE=0 -DNIST_RAND=1"

# W_DECODER = 0  :  SP
# W_DECODER = 1  :  MS
# W_DECODER = 2  :  BGF

for i in range(0, len(r)):
    os.system("rm -f test_{}".format(i))
    os.system("{} -DR_BITS={}ULL -DDV={}ULL -DT1={}ULL -DW_DECODER={} -o test_{}".format(s, r[i], w//2, t, 2, i))
    os.system("./test_{}".format(i))