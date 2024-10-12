#include "ext_decoder.h"
#include "decode.h"
#include "utilities.h"

#include "kem.h"
#include "sampling.h"

#include "ring_buffer.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

const int MAXITER = NbIter;
double error_prob = 1. * T1 / N_BITS;
double error_prob_llr = log((1 - error_prob) / error_prob);

// Algorithm SP - Sum-Product Decoder
int SP_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    int pos[R_BITS][DV*2], S[R_BITS][DV*2];
    double M[R_BITS][DV*2], E[R_BITS][DV*2], M_tmp[R_BITS][DV*2];
    double L[R_BITS*2];
    uint8_t ds[R_BITS];

    memset(e, 0, R_BITS*2);

    for ( int i = 0; i < R_BITS; ++i ){
        for ( int j = 0; j < DV; ++j ){
            pos[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
            pos[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
            
            M[i][j] = M[i][j + DV] = error_prob_llr;
        }
    }

    int T;
    for ( T = 1; T <= MAXITER; ++T ){

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

        if ( getHammingWeight(ds, R_BITS) == 0 ){
            printf("%d\n", T);
            return 0; // SUCCESS
        }
    }
    printf("101\n");

    return 1; // FAILURE
}


// Algorithm MS - Min-Sum Decoder
int MS_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    int pos[R_BITS][DV*2], S[R_BITS][DV*2];
    double M[R_BITS][DV*2], E[R_BITS][DV*2], M_tmp[R_BITS][DV*2];
    double L[R_BITS*2];
    uint8_t ds[R_BITS];

    memset(e, 0, R_BITS*2);

    for ( int i = 0; i < R_BITS; ++i ){
        for ( int j = 0; j < DV; ++j ){
            pos[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
            pos[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
            
            M[i][j] = M[i][j + DV] = error_prob_llr;
        }
    }

    int T;

    for ( T = 1; T <= MAXITER; ++T ){

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

        if ( getHammingWeight(ds, R_BITS) == 0 ){
            printf("%d\n", T);
            return 0; // SUCCESS
        }
    }

    printf("101\n");
    return 1; // FAILURE
}

// Algorithm Hybrid
int H_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    memset(e, 0, R_BITS*2);
    uint8_t ss[R_BITS];
    memcpy(ss, s, R_BITS);

    // computing the first column of each parity-check block:
    uint32_t h0_compact_col[DV] = {0};
    uint32_t h1_compact_col[DV] = {0};
    getCol(h0_compact_col, h0_compact);
    getCol(h1_compact_col, h1_compact);

    uint8_t black[R_BITS*2] = {0};
    uint8_t gray[R_BITS*2] = {0};

    for (int i = 1; i <= NbIter; i++)
    {
        memset(black, 0, R_BITS*2);
        memset(gray, 0, R_BITS*2);

        uint32_t T = floor(VAR_TH_FCT(getHammingWeight(s, R_BITS)));

        BFIter(e, black, gray, s, T, h0_compact, h1_compact, h0_compact_col, h1_compact_col);

        if (i == 1)
        {
            BFMaskedIter(e, s, black, (DV+1)/2 + 1, h0_compact, h1_compact, h0_compact_col, h1_compact_col);
            BFMaskedIter(e, s, gray, (DV+1)/2 + 1, h0_compact, h1_compact, h0_compact_col, h1_compact_col);
        }
        if (getHammingWeight(s, R_BITS) == 0){
            printf( "%d\n", i + 5 );
            return 0; // SUCCESS
        }
    }
    uint8_t e1[R_BITS*2] = {0};
    int t = SP_decoder(e1, ss, h0_compact, h1_compact);
    for ( int i = 0; i < R_BITS * 2; ++i )
        e[i] = e1[i];
    return t;
}

// Backflip

const double alpha = 0.45, beta = 1.1;
const int MAX_TTL = 5;
double ttl(double delta){
    return fmax(1, fmin(MAX_TTL, ceil(alpha * delta + beta)));
}

void BackflipIter(uint8_t e[R_BITS*2],
    int D[R_BITS*2],
    int time_stamp,
    uint8_t s[R_BITS],
    uint32_t T,
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV],
    uint32_t h0_compact_col[DV],
    uint32_t h1_compact_col[DV])
{

    uint8_t pos[R_BITS*2] = {0};

    for (uint32_t j = 0; j < R_BITS; j++)
    {
        uint32_t counter = ctr(h0_compact_col, j, s);
        if (counter >= T)
        {
            flipAdjustedErrorPosition(e, j);
            pos[j] = 1; D[j] = time_stamp + ttl(counter - T);
        }
    }
    for (uint32_t j = 0; j < R_BITS; j++)
    {
        uint32_t counter = ctr(h1_compact_col, j, s);

        if (counter >= T)
        {
            flipAdjustedErrorPosition(e, R_BITS+j);
            pos[R_BITS+j] = 1; D[R_BITS+j] = time_stamp + ttl(counter - T);
        }
    }

    for(uint32_t j=0; j < 2*R_BITS; j++){
        if(pos[j] == 1){
            recompute_syndrome(s, j, h0_compact, h1_compact);
        }
    }
}

int Backflip_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    memset(e, 0, R_BITS*2);

    // computing the first column of each parity-check block:
    uint32_t h0_compact_col[DV] = {0};
    uint32_t h1_compact_col[DV] = {0};
    getCol(h0_compact_col, h0_compact);
    getCol(h1_compact_col, h1_compact);

    int D[R_BITS*2] = {0};

    for (int i = 1; i <= NbIter; i++)
    {
        for ( int j = 0; j < 2 * R_BITS; ++j ){
            if ( D[j] == i ){
                flipAdjustedErrorPosition(e, j);
                recompute_syndrome(s, j, h0_compact, h1_compact);
            }
        }

        uint32_t T = floor(VAR_TH_FCT(getHammingWeight(s, R_BITS)));
        BackflipIter(e, D, i, s, T, h0_compact, h1_compact, h0_compact_col, h1_compact_col);
        
        if (getHammingWeight(s, R_BITS) == 0){
            printf( "%d\n", i );
            return 0; // SUCCESS
        }
    }
    printf( "101\n" );
    return 1; // FAILURE
}
