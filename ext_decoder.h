#ifndef __EXT_DECODER__
#define __EXT_DECODER__

#include "types.h"
#include "conversions.h"

int SP_decoder(uint8_t e[R_BITS*2],
        uint8_t s[R_BITS],
        uint32_t h0_compact[DV],
        uint32_t h1_compact[DV]);

int MS_decoder(uint8_t e[R_BITS*2],
        uint8_t s[R_BITS],
        uint32_t h0_compact[DV],
        uint32_t h1_compact[DV]);

int H_decoder(uint8_t e[R_BITS*2],
        uint8_t s[R_BITS],
        uint32_t h0_compact[DV],
        uint32_t h1_compact[DV]);

#endif