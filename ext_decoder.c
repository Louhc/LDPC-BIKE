#include "ext_decoder.h"+
#include "decode.h"
#include "utilities.h"

#include "kem.h"
#include "sampling.h"

#include "ring_buffer.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

const int MAXITER = 50;
double error_prob = 1. * T1 / N_BITS;
double error_prob_llr = log((1 - error_prob) / error_prob);

int pos[R_BITS][DV*2], S[R_BITS][DV*2];
double M[R_BITS][DV*2], E[R_BITS][DV*2], M_tmp[R_BITS][DV*2];
double L[R_BITS*2];
uint8_t ds[R_BITS];

// Algorithm SP - Sum-Product Decoder
int SP_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{

    memset(e, 0, R_BITS*2);

    for ( int i = 0; i < R_BITS; ++i ){
        for ( int j = 0; j < DV; ++j ){
            pos[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
            pos[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
            
            M[i][j] = M[i][j + DV] = error_prob_llr;
        }
    }

    for ( int T = 0; T < MAXITER; ++T ){

        for ( int i = 0; i < R_BITS; ++i ){
            
            double temp = 1.0;
            
            for ( int j = 0; j < DV * 2; ++j ){
                E[i][j] = temp;
                temp *= tanh(M[i][j] / 2);
            }

            temp = 1.0;
            for ( int j = DV * 2 - 1; j >= 0; --j ){
                E[i][j] *= temp;
                E[i][j] = (s[i] == 0 ? 1 : -1) * log((1 + E[i][j]) / (1 - E[i][j]));
                temp *= tanh(M[i][j] / 2);
            }

        }
        
        for ( int i = 0; i < R_BITS * 2; ++i )
            L[i] = error_prob_llr;
        for ( int i = 0; i < R_BITS; ++i ){
            for ( int j = 0; j < DV * 2; ++j ){
                M[i][j] = L[pos[i][j]];
                L[pos[i][j]] += E[i][j];
            }
        }
        for ( int j = 0; j < R_BITS * 2; ++j ){
            e[j] = L[j] > 0 ? 0 : 1;
            L[j] = 0;
        }
        for ( int i = R_BITS - 1; i >= 0; --i ){
            for ( int j = 0; j < DV * 2; ++j ){
                M[i][j] += L[pos[i][j]];
                L[pos[i][j]] += E[i][j];
            }
        }

        for ( int i = 0; i < R_BITS; ++i ){
            ds[i] = s[i];
            for ( int j = 0; j < DV * 2; ++j ){
                ds[i] ^= e[pos[i][j]];
            }
        }

        if ( getHammingWeight(ds, R_BITS) == 0 )
            return 0; // SUCCESS
    }

    return 1; // FAILURE
}


// Algorithm MS - Min-Sum Decoder
int MS_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{

    memset(e, 0, R_BITS*2);

    for ( int i = 0; i < R_BITS; ++i ){
        for ( int j = 0; j < DV; ++j ){
            pos[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
            pos[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
            
            M[i][j] = M[i][j + DV] = error_prob_llr;
        }
    }

    for ( int T = 1; T <= MAXITER; ++T ){

        double alpha = 1.0 - pow(2.0, -1.0 * T);

        for ( int i = 0; i < R_BITS; ++i ){
            
            double temp = 1e308;
            int sgn = s[i];
            
            for ( int j = 0; j < DV * 2; ++j ){
                E[i][j] = temp;
                S[i][j] = sgn;
                if ( abs(M[i][j]) < temp )
                    temp = abs(M[i][j]);
                if ( M[i][j] <= 0 )
                    ++sgn;
            }

            temp = 1e308;
            sgn = 0;
            for ( int j = DV * 2 - 1; j >= 0; --j ){
                if ( temp < E[i][j] )
                    E[i][j] = temp;
                S[i][j] += sgn;

                E[i][j] *= ((S[i][j] & 1) ? -1 : 1) * alpha;

                if ( abs(M[i][j]) < temp )
                    temp = abs(M[i][j]);
                if ( M[i][j] <= 0 )
                    ++sgn;
            }
        }
        
        for ( int i = 0; i < R_BITS * 2; ++i )
            L[i] = error_prob_llr;
        for ( int i = 0; i < R_BITS; ++i ){
            for ( int j = 0; j < DV * 2; ++j ){
                M_tmp[i][j] = L[pos[i][j]];
                L[pos[i][j]] += E[i][j];
            }
        }
        for ( int j = 0; j < R_BITS * 2; ++j ){
            e[j] = L[j] > 0 ? 0 : 1;
            L[j] = 0;
        }
        for ( int i = R_BITS - 1; i >= 0; --i ){
            for ( int j = 0; j < DV * 2; ++j ){
                M_tmp[i][j] += L[pos[i][j]];
                L[pos[i][j]] += E[i][j];
            }
        }
        for ( int i = 0; i < R_BITS; ++i ){
            for ( int j = 0; j < DV * 2; ++j ){
                // M[i][j] = M_tmp[i][j];
                if ( M[i][j] * M_tmp[i][j] >= -1e-15 ) M[i][j] = M_tmp[i][j];
                else M[i][j] = 0;
            }
        }

        for ( int i = 0; i < R_BITS; ++i ){
            ds[i] = s[i];
            for ( int j = 0; j < DV * 2; ++j ){
                ds[i] ^= e[pos[i][j]];
            }
        }

        if ( getHammingWeight(ds, R_BITS) == 0 )
            return 0; // SUCCESS
    }

    return 1; // FAILURE
}