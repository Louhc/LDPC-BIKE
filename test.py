import os
import time

w = 206
t = 199
r = 19139

s = "g++ -m64 -O3 ldpc_tests/test.c *.c ntl.cpp FromNIST/rng.c -I. -I/include -L/lib -std=c++11 -lcrypto -lssl -lm -ldl -lntl -lgmp -lgf2x -lpthread -DVERBOSE=0 -DNIST_RAND=1"

# W_DECODER = 0  :  SP
# W_DECODER = 1  :  MS
# W_DECODER = 2  :  BGF

start_t = time.time()
os.system("rm -f test_time_1")
os.system("{} -DR_BITS={}ULL -DDV={}ULL -DT1={}ULL -DW_DECODER={} -o test_time_1".format(s, r, w//2, t, 0))
os.system("./test_time_1 > output1")
end_t = time.time()

print("SP : {}".format(end_t - start_t))

# start_t = time.time()
# os.system("rm -f test_time_2")
# os.system("{} -DR_BITS={}ULL -DDV={}ULL -DT1={}ULL -DW_DECODER={} -o test_time_2".format(s, r, w//2, t, 2))
# os.system("./test_time_2 > output2")
# end_t = time.time()

# print("BGF: {}".format(end_t - start_t))