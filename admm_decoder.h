#ifndef __ADMM_DECODER__
#define __ADMM_DECODER__

#include "types.h"
#include "conversions.h"

int ADMM_decoder(uint8_t e[R_BITS*2],
        uint8_t s[R_BITS],
        uint32_t h0_compact[DV],
        uint32_t h1_compact[DV]);

#endif