import os
import time

# W_DECODER = 0  :  SP
# W_DECODER = 1  :  MS
# W_DECODER = 2  :  BGF

w = 142
t = 134
r = [11239, 11273, 11317, 11369, 11423, 11471, 11503, 11579, 11621, 11689, 11731, 11789, 11827, 11867, 11923, 11953, 11987, 12043, 12101, 12143, 12197, 12241, 12277]
X = [0.0914770912009593, 0.09084042223745625, 0.09002519223018442, 0.0890742014270221, 0.08810068573265407, 0.0872471520527872, 0.08668421490460705, 0.08536641484570438, 0.08464954431537437, 0.08350573395613703, 0.08280948974331248, 0.08186061409165839, 0.08124674597749425, 0.0806071494528153, 0.07972288905502715, 0.07925446822929251, 0.07872799546790363, 0.07787093365082587, 0.07699626464667421, 0.07637100569826508, 0.07557696263379413, 0.0749380456944261, 0.07442061654818351]

s = "g++ -m64 -O3 ldpc_tests/test.c *.c ntl.cpp FromNIST/rng.c -I. -I/include -L/lib -std=c++11 -lcrypto -lssl -lm -ldl -lntl -lgmp -lgf2x -lpthread -DVERBOSE=0 -DNIST_RAND=1"

begin_t = time.time()
for i in range(0, len(r)):
    os.system("rm -f LDPC_test_spt1_{}".format(i))
    para = "-DR_BITS={}ULL -DDV={}ULL -DT1={}ULL -DW_DECODER={} -DTHR_X={}".format(r[i], w//2, t, 2, X[i])
    os.system("{} {} -o LDPC_test_spt1_{}".format(s, para, i))
    os.system("./LDPC_test_spt1_{} > output_sp1_t{}".format(i, i))
    end_t = time.time()
    print(end_t - begin_t)
