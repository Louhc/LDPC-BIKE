#include "stdio.h"
#include "kem.h"
#include "utilities.h"
#include "measurements.h"

#define T_TEST 48000000

int main(){
    sk_t sk    = {0};
    pk_t pk    = {0};
    ct_t ct    = {0};
    ss_t k_enc = {0};
    ss_t k_dec = {0};

    MSG("BIKE LDPC Test:\n");

    for (uint32_t i=1; i <= T_TEST; ++i)
    {
        status_t res = SUCCESS;

        MSG("Test %d:\n", i);

        res = static_cast<status_t>(crypto_kem_keypair(pk.raw, sk.raw));
        if(res != SUCCESS)
        {
            MSG("Keypair failed with error: %d\n", res);
            continue;
        }

        uint32_t dec_rc = 0;
        
        res = static_cast<status_t>(crypto_kem_enc(ct.raw, k_enc.raw, pk.raw));
        if(res != SUCCESS)
        {
            MSG("encapsulate failed with error: %d\n", res);
            continue;
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

    return 0;
}