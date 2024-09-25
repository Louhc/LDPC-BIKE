import os

# W_DECODER = 0  :  SP
# W_DECODER = 1  :  MS
# W_DECODER = 2  :  BGF

w = 206
t = 199
r = [19139, 19141, 19157, 19163, 19181, 19219, 19237, 19259, 19301, 19333, 19373, 19379, 19387, 19403, 19427, 19469, 19483, 19501, 19507,19541]

s = "g++ -m64 -O3 ldpc_tests/test.c *.c ntl.cpp FromNIST/rng.c -I. -I/include -L/lib -std=c++11 -lcrypto -lssl -lm -ldl -lntl -lgmp -lgf2x -lpthread -DVERBOSE=0 -DNIST_RAND=1"

for i in range(0, len(r)):
    os.system("rm -f LDPC_test_ms3_{}".format(i))
    os.system("{} -DR_BITS={}ULL -DDV={}ULL -DT1={}ULL -DW_DECODER={} -o LDPC_test_ms3_{}".format(s, r[i], w//2, t, 1, i))
    os.system("./LDPC_test_ms3_{} > output_ms3_{}".format(i, i))