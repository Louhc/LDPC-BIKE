#include "stdio.h"
#include "kem.h"
#include "utilities.h"
#include "measurements.h"
#include "pthread.h"
#include "FromNIST/rng.h"
#include <time.h>

#define T_TEST 10000

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
}

int main(){
    // unsigned char entropy_input[48] = INIT_SEED;
    // unsigned char personalization_string[48];
    // memset(personalization_string, 0x00, 48);
    // randombytes_init(entropy_input, personalization_string, 0);

    MSG("BIKE LDPC Test:\n");

    for ( int i = 0; i < T_TEST; ++i ){
        test_once();
    }

    return 0;
}
