#include "stdio.h"
#include "kem.h"
#include "utilities.h"
#include "measurements.h"
#include "pthread.h"

#define T_TEST 48000
#define NUM_THREADS 100

void test_once(){
    sk_t sk    = {0};
    pk_t pk    = {0};
    ct_t ct    = {0};
    ss_t k_enc = {0};
    ss_t k_dec = {0};

    status_t res = SUCCESS;

    res = static_cast<status_t>(crypto_kem_keypair(pk.raw, sk.raw));
    if(res != SUCCESS)
    {
        MSG("Keypair failed with error: %d\n", res);
        return;
    }

    uint32_t dec_rc = 0;
    
    res = static_cast<status_t>(crypto_kem_enc(ct.raw, k_enc.raw, pk.raw));
    if(res != SUCCESS)
    {
        MSG("encapsulate failed with error: %d\n", res);
        return;
    }

    dec_rc = crypto_kem_dec(k_dec.raw, ct.raw, sk.raw);
    if (dec_rc != 0)
    {
        MSG("Decoding failed after %d code tests!\n", i);
    }
    else
    {
        if (safe_cmp(k_enc.raw, k_dec.raw, sizeof(k_dec)/sizeof(uint64_t)))
        {
            MSG("Success! decapsulated key is the same as encapsulated key!\n");
        } else {
            MSG("Failure! decapsulated key is NOT the same as encapsulated key!\n");
        }
    }
}

void* test(void *threadid) {
    int tid = (long)threadid;
    for ( int i = 0; i * NUM_THREADS + tid < T_TEST; ++i ){
        test_once();
    }
    pthread_exit(NULL);
}

int main(){
    MSG("BIKE LDPC Test:\n");

    pthread_t threads[NUM_THREADS];

    for ( int i = 0; i < NUM_THREADS; ++i ){
        int rc = pthread_create(threads + i, NULL, test, (void *)((long)i));
        if ( rc ){
            printf( "ERROR: return code from pthread_create is %d\n", rc);
            exit(-1);
        }
    }

    pthread_exit(NULL);

    return 0;
}