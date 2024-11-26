#include "ext_decoder.h"
#include "decode.h"
// #include "admm_decode.h"
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
    double M[R_BITS][DV*2], E[R_BITS][DV*2];
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
            for ( int i = 0; i < 2 * R_BITS; ++i )
                s[i] = 0;
            printf("%d\n", T);
            return 0; // SUCCESS
        }
    }
    for ( int i = 0; i < R_BITS; ++i )
        s[i] = ds[i];
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
    int t = 1;
    // printf("SP "); t = SP_decoder(e, s, h0_compact, h1_compact);
    // if ( t == 1 ){
    //     uint8_t e0[R_BITS * 2];
    //     printf("BGF "); t = BGF_decoder(e0, s, h0_compact, h1_compact);
    //     for ( int i = 0; i < R_BITS * 2; ++i ) e[i] ^= e0[i];
    // }
    printf("BGF "); t = BGF_decoder(e, s, h0_compact, h1_compact);
    if ( t == 1 ){
        uint8_t e0[R_BITS * 2];
        printf("SP "); t = SP_decoder(e0, s, h0_compact, h1_compact);
        for ( int i = 0; i < R_BITS * 2; ++i ) e[i] ^= e0[i];
    }
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

// ----------------

void quicksort_dec(int first, int last, int *idx, int *dat){
    int i, j, tmp; int x;
    x = dat[idx[(first + last) / 2]];
    i = first, j = last;
    for (;;) {
        while (dat[idx[i]] > x) i++;
        while (x > dat[idx[j]]) j--;
        if (i >= j) break;
        tmp = idx[i], idx[i] = idx[j], idx[j] = tmp;
        i++, j--;
    }
    if (first < i - 1) quicksort_dec(first, i - 1, idx, dat);
    if (j + 1 < last) quicksort_dec(j + 1, last, idx, dat);
}

void SSUP(int pos[R_BITS][DV * 2], uint8_t s[R_BITS], int np, int nt, int *S, int &n){
    int sigma[2 * R_BITS] = {0}, idx[2 * R_BITS] = {0};
    uint8_t a[2 * R_BITS] = {0}, b[2 * R_BITS] = {0};

    for ( uint32_t i = 0; i < R_BITS; ++i ){
        if ( s[i] == 0 ) continue;
        for (uint32_t j = 0; j < 2 * DV; ++j )
            ++sigma[pos[i][j]];
    }
    for ( uint32_t i = 0; i < 2 * R_BITS; ++i )
        idx[i] = i;
    quicksort_dec(0, 2 * R_BITS, idx, sigma);
    n = 0;
    for ( uint32_t i = 0; i < np; ++i ){
        S[n++] = idx[i];
        a[idx[i]] = 1;
    }
    for (uint32_t i = 0; i < R_BITS; ++i ){
        uint8_t flag = 0;
        for ( uint32_t j = 0; j < 2 * DV; ++j ){
            if ( a[pos[i][j]] ){
                flag = 1; break;
            }
        }
        if ( flag ){
            for ( uint32_t j = 0; j < 2 * DV; ++j ){
                if ( !a[pos[i][j]] ) b[pos[i][j]] = 1;
            }
        }
    }

    for ( uint32_t i = 0; i < 2 * R_BITS && n < nt; ++i ){
        if ( b[i] ) S[n++] = i;
    }
}

int AAVN_decoder(uint8_t ee[R_BITS*2],
    uint8_t ss[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    uint8_t e0[R_BITS*2] = {0};
    printf("BGF ");
    int status = BGF_decoder(e0, ss, h0_compact, h1_compact);
    if ( status == 0 ){
        for ( int i = 0; i < R_BITS * 2; ++i ) ee[i] = e0[i];
        return 0;
    }

    // ----------------------------------------

    int n = R_BITS, m = R_BITS * 2, n0 = 0, nmax = R_BITS + T1;
    int ns = 0;
    int pos0[R_BITS][DV*2], SUP[3 * T1];

    for ( int i = 0; i < R_BITS; ++i ){
        for ( int j = 0; j < DV; ++j ){
            pos0[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
            pos0[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
        }
    }
    SSUP(pos0, ss, T1, T1 * 3, SUP, ns);
    
    int a[T1 * 3][DV * 2] = {0};
    int cw[T1 * 3] = {0}, inSUP[R_BITS * 2];
    int pos[R_BITS + T1][DV * 2 + T1] = {0}, w[R_BITS + T1] = {0}; uint8_t s[R_BITS + T1] = {0};
    uint8_t f[R_BITS] = {0}, fa[T1 * 3][R_BITS] = {0};
    int S[R_BITS + T1][DV*2 + T1];
    double M[R_BITS + T1][DV*2 + T1], E[R_BITS + T1][DV*2 + T1];
    double L[R_BITS*2 + T1];
    uint8_t ds[R_BITS + T1], e[R_BITS * 2 + T1];

    for ( uint32_t i = 0; i < R_BITS * 2; ++i )
        inSUP[i] = -1;
    for ( uint32_t i = 0; i < ns; ++i ) inSUP[SUP[i]] = i;
    for ( uint32_t i = 0; i < R_BITS; ++i ){
        for ( uint32_t j = 0; j < 2 * DV; ++j ){
            if ( inSUP[pos0[i][j]] >= 0 ) a[inSUP[pos0[i][j]]][cw[inSUP[pos0[i][j]]]++] = i;
        }
    }
    for ( uint32_t i = 0; i < ns && n + n0 < nmax; ++i ){
        for ( uint32_t j = i + 1; j < ns && n + n0 < nmax; ++j ){
            int b[DV * 2], cnt = 0;
            for ( uint32_t k = 0; k < cw[j]; ++k )
                if ( a[j][k] >= 0 ) f[a[j][k]] = 0;
            for ( uint32_t k = 0; k < cw[i]; ++k )
                if ( a[i][k] >= 0 ) f[a[i][k]] = 1;
            for ( uint32_t k = 0; k < cw[j]; ++k )
                if ( a[j][k] >= 0 && f[a[j][k]] ) b[cnt++] = a[j][k], f[a[j][k]] = 2;
            if ( cnt >= 2 ){
                for ( uint32_t k = 0; k < cnt; ++k )
                    fa[i][b[k]] = fa[j][b[k]] = 1;
                for ( uint32_t k = 0; k < cw[i]; ++k )
                    if ( a[i][k] >= 0 && f[a[i][k]] == 2 ) pos[a[i][k]][w[a[i][k]]++] = 2 * R_BITS + n0, a[i][k] = -1;
                for ( uint32_t k = 0; k < cw[j]; ++k )
                    if ( a[j][k] >= 0 && f[a[j][k]] == 2 ) a[j][k] = -1;
                pos[R_BITS + n0][w[R_BITS + n0]++] = SUP[i];
                pos[R_BITS + n0][w[R_BITS + n0]++] = SUP[j];
                pos[R_BITS + n0][w[R_BITS + n0]++] = 2 * R_BITS + n0;
                s[R_BITS + n0] = 0;
                ++n0;
            }
        }
    }
    for ( uint32_t i = 0; i < R_BITS; ++i ){
        s[i] = ss[i];
        for ( uint32_t j = 0; j < 2 * DV; ++j ){
            if ( inSUP[pos0[i][j]] < 0 || !fa[inSUP[pos0[i][j]]][i] ){
                pos[i][w[i]++] = pos0[i][j];
            }
        }
    } n += n0; m += n0;

    for ( uint32_t i = 0; i < n; ++i )
        for ( uint32_t j = 0; j < w[i]; ++j )
            M[i][j] = pos[i][j] < 2 * R_BITS ? error_prob_llr : 0;
    return 1;

    // ----------------------------------------

    int T;
    for ( T = 1; T <= MAXITER; ++T ){

        for ( int i = 0; i < n; ++i ){
            
            double temp = 1.0;
            for ( int j = 0; j < w[i]; ++j ){
                E[i][j] = temp;
                temp *= tanh(M[i][j] / 2);
            }

            temp = 1.0;
            for ( int j = w[i] - 1; j >= 0; --j ){
                E[i][j] *= temp;
                E[i][j] = (s[i] == 0 ? 1 : -1) * log((1 + E[i][j]) / (1 - E[i][j]));
                temp *= tanh(M[i][j] / 2);
            }

        }
        
        for ( int i = 0; i < R_BITS * 2; ++i )
            L[i] = error_prob_llr;
        for ( int i = R_BITS * 2; i < m; ++i )
            L[i] = 0;
        for ( int i = 0; i < n; ++i ){
            for ( int j = 0; j < w[i]; ++j ){
                M[i][j] = L[pos[i][j]];
                L[pos[i][j]] += E[i][j];
            }
        }
        for ( int j = 0; j < m; ++j ){
            e[j] = L[j] > 0 ? 0 : 1;
            L[j] = 0;
        }
        for ( int j = 0; j < R_BITS * 2; ++j ) ee[j] = e[j] ^ e0[j];

        for ( int i = n - 1; i >= 0; --i ){
            for ( int j = 0; j < w[i]; ++j ){
                M[i][j] += L[pos[i][j]];
                L[pos[i][j]] += E[i][j];
            }
        }
        
        for ( int i = 0; i < n; ++i ){
            ds[i] = s[i];
            for ( int j = 0; j < w[i]; ++j ){
                ds[i] ^= e[pos[i][j]];
            }
        }

        if ( getHammingWeight(ds, n) == 0 ){
            printf("SPA %d\n", T);
            return 0; // SUCCESS
        }
    }
    printf("SPA 101\n");
    return 1; // FAILURE
}

// ----------------

void SPn_Iter(int n, int pos[R_BITS][DV * 2], uint8_t s[R_BITS], uint8_t e[R_BITS * 2]){
    double M[R_BITS][DV*2], E[R_BITS][DV*2];
    double L[R_BITS*2];
    uint8_t ds[R_BITS];
    int T;
    for ( int i = 0; i < R_BITS; ++i )
        for ( int j = 0; j < 2 * DV; ++j )
            M[i][j] = error_prob_llr;

    for ( T = 1; T <= n; ++T ){
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
            break;
    }
    for ( int i = 0; i < R_BITS; ++i )
        s[i] = ds[i];
}

int RSP_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    int pos[R_BITS][DV*2];
    memset(e, 0, R_BITS*2);

    for ( int i = 0; i < R_BITS; ++i ){
        for ( int j = 0; j < DV; ++j ){
            pos[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
            pos[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
        }
    }
    for ( int T = 1; T <= NbIter; ++T ){
        uint8_t e0[R_BITS * 2] = {0};
        SPn_Iter(50, pos, s, e0);
        for ( int i = 0; i < R_BITS * 2; ++i )
            e[i] ^= e0[i];
        if ( getHammingWeight(s, R_BITS) == 0 ){
            printf("%d\n", T);
            return 0;
        }
    }
    printf("%d\n", NbIter + 1);

    return 1; // FAILURE
}

// -------------------

// double ext_threshold(int i, const uint8_t s[R_BITS + T1], const uint8_t s_tilte[R_BITS]){
//     uint32_t T_prime = floor(NEW_VAR_TH_FCT(getHammingWeight(s, R_BITS)));
//     uint32_t M = (DV + 1)/2;
//     uint32_t T;
//     if ( i==1 ){
//         T = floor(T_prime + delta);
//     }
//     if ( i==2 ){
//         T = floor(2*T_prime + M)/3 + delta;
//     }
//     if ( i==3 ){
//         T = floor(T_prime + 2*M)/3 + delta;
//     }
//     if ( i>=4 ){
//         T = M + delta;
//     }
//     return MAX(VAR_TH_FCT(getHammingWeight(s_tilte, R_BITS)), T);
// }

// int extBGF_decoder(uint8_t ee[R_BITS*2],
//     uint8_t ss[R_BITS],
//     uint32_t h0_compact[DV],
//     uint32_t h1_compact[DV])
// {
//     uint8_t e0[R_BITS*2] = {0};
//     printf("BGF ");
//     int status = BGF_decoder(e0, ss, h0_compact, h1_compact);
//     if ( status == 0 ){
//         for ( int i = 0; i < R_BITS * 2; ++i ) ee[i] = e0[i];
//         return 0;
//     }

//     // ----------------------------------------

//     int n = R_BITS, m = R_BITS * 2, n0 = 0, nmax = R_BITS + T1;
//     int ns = 0;
//     int pos0[R_BITS][DV*2], SUP[3 * T1];

//     for ( int i = 0; i < R_BITS; ++i ){
//         for ( int j = 0; j < DV; ++j ){
//             pos0[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
//             pos0[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
//         }
//     }
//     SSUP(pos0, ss, T1, T1 * 3, SUP, ns);
    
//     int a[T1 * 3][DV * 2] = {0};
//     int cw[T1 * 3] = {0}, inSUP[R_BITS * 2];
//     int pos[R_BITS + T1][DV * 2 + T1] = {0}, w[R_BITS + T1] = {0}; uint8_t s[R_BITS + T1] = {0};
//     uint8_t f[R_BITS] = {0}, fa[T1 * 3][R_BITS] = {0};
//     int S[R_BITS + T1][DV*2 + T1];
//     double M[R_BITS + T1][DV*2 + T1], E[R_BITS + T1][DV*2 + T1];
//     double L[R_BITS*2 + T1];
//     uint8_t ds[R_BITS + T1], e[R_BITS * 2 + T1];

//     for ( uint32_t i = 0; i < R_BITS * 2; ++i )
//         inSUP[i] = -1;
//     for ( uint32_t i = 0; i < ns; ++i ) inSUP[SUP[i]] = i;
//     for ( uint32_t i = 0; i < R_BITS; ++i ){
//         for ( uint32_t j = 0; j < 2 * DV; ++j ){
//             if ( inSUP[pos0[i][j]] >= 0 ) a[inSUP[pos0[i][j]]][cw[inSUP[pos0[i][j]]]++] = i;
//         }
//     }
//     for ( uint32_t i = 0; i < ns && n + n0 < nmax; ++i ){
//         for ( uint32_t j = i + 1; j < ns && n + n0 < nmax; ++j ){
//             int b[DV * 2], cnt = 0;
//             for ( uint32_t k = 0; k < cw[j]; ++k )
//                 if ( a[j][k] >= 0 ) f[a[j][k]] = 0;
//             for ( uint32_t k = 0; k < cw[i]; ++k )
//                 if ( a[i][k] >= 0 ) f[a[i][k]] = 1;
//             for ( uint32_t k = 0; k < cw[j]; ++k )
//                 if ( a[j][k] >= 0 && f[a[j][k]] ) b[cnt++] = a[j][k], f[a[j][k]] = 2;
//             if ( cnt >= 2 ){
//                 for ( uint32_t k = 0; k < cnt; ++k )
//                     fa[i][b[k]] = fa[j][b[k]] = 1;
//                 for ( uint32_t k = 0; k < cw[i]; ++k )
//                     if ( a[i][k] >= 0 && f[a[i][k]] == 2 ) pos[a[i][k]][w[a[i][k]]++] = 2 * R_BITS + n0, a[i][k] = -1;
//                 for ( uint32_t k = 0; k < cw[j]; ++k )
//                     if ( a[j][k] >= 0 && f[a[j][k]] == 2 ) a[j][k] = -1;
//                 pos[R_BITS + n0][w[R_BITS + n0]++] = SUP[i];
//                 pos[R_BITS + n0][w[R_BITS + n0]++] = SUP[j];
//                 pos[R_BITS + n0][w[R_BITS + n0]++] = 2 * R_BITS + n0;
//                 s[R_BITS + n0] = 0;
//                 ++n0;
//             }
//         }
//     }
//     for ( uint32_t i = 0; i < R_BITS; ++i ){
//         s[i] = ss[i];
//         for ( uint32_t j = 0; j < 2 * DV; ++j ){
//             if ( inSUP[pos0[i][j]] < 0 || !fa[inSUP[pos0[i][j]]][i] ){
//                 pos[i][w[i]++] = pos0[i][j];
//             }
//         }
//     } n += n0; m += n0;

//     for ( uint32_t i = 0; i < n; ++i )
//         for ( uint32_t j = 0; j < w[i]; ++j )
//             M[i][j] = pos[i][j] < 2 * R_BITS ? error_prob_llr : 0;
//     return 1;

//     // ----------------------------------------

//     int T;
//     for ( T = 1; T <= MAXITER; ++T ){

        

//         if ( getHammingWeight(ds, n) == 0 ){
//             printf("SPA %d\n", T);
//             return 0; // SUCCESS
//         }
//     }
//     printf("SPA 101\n");
//     return 1; // FAILURE
// }