#include "ext_decoder.h"
#include "decode.h"
#include "utilities.h"
#include "admm_decoder.h"
#include "ntl.h"

#include "kem.h"
#include "sampling.h"

#include "ring_buffer.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

const int MU = 6, ALPHA2 = 15;
const int dc = DV * 2, dv = DV, ERRN = T1;
const int n = R_BITS * 2, m = R_BITS;
const int ITRN = NbIter;
bool flag = 0;

inline double cutoff(double x) { return x > 1 ? 1 : (x < 0 ? 0 : x); }
inline double sign(double x) { return x > 0 ? 1 : -1; }
inline int hardd(double x) { return x >= 0 ? 0 : 1; }
inline int abs_(int x) { return x >= 0 ? x : -x; }

void quicksort(int first, int last, int *idx, double *dat)
{
    int i, j, tmp;
    double x;

    x = dat[idx[(first + last) / 2]];
    i = first;
    j = last;
    for (;;)
    {
        while (dat[idx[i]] < x)
            i++;
        while (x < dat[idx[j]])
            j--;
        if (i >= j)
            break;
        tmp = idx[i];
        idx[i] = idx[j];
        idx[j] = tmp;
        i++;
        j--;
    }
    if (first < i - 1)
        quicksort(first, i - 1, idx, dat);
    if (j + 1 < last)
        quicksort(j + 1, last, idx, dat);
}
void quicksort_dec(int first, int last, int *idx, double *dat)
{
    int i, j, tmp;
    double x;

    x = dat[idx[(first + last) / 2]];
    i = first;
    j = last;
    for (;;)
    {
        while (dat[idx[i]] > x)
            i++;
        while (x > dat[idx[j]])
            j--;
        if (i >= j)
            break;
        tmp = idx[i];
        idx[i] = idx[j];
        idx[j] = tmp;
        i++;
        j--;
    }
    if (first < i - 1)
        quicksort_dec(first, i - 1, idx, dat);
    if (j + 1 < last)
        quicksort_dec(j + 1, last, idx, dat);
}

int proj(double v[m][dc], double *z, int row) {
    z[0] = 1;
    int i, j, k, beta_max, r, L, R, M, plus;
    double tmp, frz, frz_pre, frzero, beta_opt, z_zero_sum;

    double z_hat[dc], fr[dc], Beta_Set[dc * 4];
    int idx_p[dc * 4];

    /* r */
    tmp = 0;
    for (k = 0; k < dc; k++) tmp += (z_hat[k] = cutoff(v[row][k]));

    r = (int)tmp;
    if (r & 1) r ^= 1;

    if (r == dc || r == dc - 1)
    {
        for (k = 0; k < dc; k++) z[k] = z_hat[k];
        return (1);
    }

    /* fr */
    for ( k = 0; k <= r; ++k ) fr[k] = 1;
    for ( k = r + 1; k < dc; ++k) fr[k] = -1;

    tmp = 0;
    for (k = 0; k < dc; k++) tmp += fr[k] * z_hat[k];
    if (tmp <= r + 0.00001)
    { // 13:
        for (k = 0; k < dc; k++) z[k] = z_hat[k];
        return (1);
    }

    //****************************************************************

    tmp = 0.5 * (v[row][r] - v[row][r + 1]) + 0.00000000001;
    Beta_Set[0] = 0;

    plus = 1;
    for (k = 0; k < dc; k++)
    { // 16:
        if (k > 0 && v[row][k] > v[row][k - 1] - 0.00000000001 && v[row][k] < v[row][k - 1] + 0.00000000001)
            continue;
        if (k <= r)
        {
            if ((v[row][k] - 1 > 0) && ((v[row][k] - 1) < tmp))
            {
                Beta_Set[plus] = v[row][k] - 1;
                plus++;
            }
        }
        else
        {
            if ((-v[row][k] > 0) && ((-v[row][k]) < tmp))
            {
                Beta_Set[plus] = -v[row][k];
                plus++;
            }
        }
    }

    for (k = 0; k < dc; k++)
    { // 16:
        if (k <= r)
        {
            if ((v[row][k] > 0) && (v[row][k] < tmp))
            {
                Beta_Set[plus] = v[row][k];
                plus++;
            }
        }
        else
        {
            if ((1 - v[row][k] > 0) && ((1 - v[row][k]) < tmp))
            {
                Beta_Set[plus] = 1 - v[row][k];
                plus++;
            }
        }
    }

    Beta_Set[plus] = tmp;
    for (k = 0; k < 2 * dc; k++) idx_p[k] = k;
    quicksort(0, plus, idx_p, Beta_Set);

    L = 0;
    R = plus;
    M = (L + R) / 2;
    while (R - L != 1)
    { // Bi-section Method
        for (j = 0; j <= r; j++)
        {
            z_hat[j] = v[row][j] - Beta_Set[idx_p[M]];
        }
        for (j = r + 1; j < dc; j++)
        {
            z_hat[j] = v[row][j] + Beta_Set[idx_p[M]];
        }
        tmp = 0;
        for (j = 0; j < dc; j++)
        {
            tmp += fr[j] * cutoff(z_hat[j]);
        }
        if (tmp < r) R = M;
        else L = M;
        M = (L + R) >> 1;
    }

    for (j = 0; j <= r; j++)
    {
        z_hat[j] = v[row][j] - Beta_Set[idx_p[R]];
    }
    for (j = r + 1; j < dc; j++)
    {
        z_hat[j] = v[row][j] + Beta_Set[idx_p[R]];
    }
    tmp = 0;
    for (j = 0; j < dc; j++)
    {
        tmp += fr[j] * cutoff(z_hat[j]);
    }
    frz = tmp;
    for (j = 0; j <= r; j++)
    {
        z_hat[j] = v[row][j] - Beta_Set[idx_p[L]];
    }
    for (j = r + 1; j < dc; j++)
    {
        z_hat[j] = v[row][j] + Beta_Set[idx_p[L]];
    }
    tmp = 0;
    for (j = 0; j < dc; j++)
    {
        tmp += fr[j] * cutoff(z_hat[j]);
    }
    frz_pre = tmp;
    beta_opt = (r - frz) * (Beta_Set[idx_p[R]] - Beta_Set[idx_p[L]]) / (frz - frz_pre) + Beta_Set[idx_p[R]];

    for (j = 0; j <= r; j++) z[j] = cutoff(v[row][j] - beta_opt);
    for (j = r + 1; j < dc; j++) z[j] = cutoff(v[row][j] + beta_opt);
    return (1);
}

void convert( uint32_t *h ){
    for ( int i = 0; i < DV; ++i )
        h[i] = (R_BITS - h[i]) % R_BITS;
}

double getETime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}


int ADMM_decoder(uint8_t e[R_BITS * 2],
                 uint8_t s[R_BITS],
                 uint32_t h0_compact[DV],
                 uint32_t h1_compact[DV])
{
    uint8_t ts[R_BITS];
    int i, j, k, itr_num = 0, err, posit;
    double v[m][dc], z[m][dc], lambda[m][dc];
    double tmp, tmp2, prev = 0;
    int idx[dc], iidx[dc];
    double v_tmp[dc], z_tmp[dc];
    double admmsol[n];
    int admmbin_sol[n];
    double Obj[n], receivedword[n];
    double sigma = 0;
    int Pj[m][dc], POS[n][dc];
    uint8_t h0[R_SIZE] = {0}, inv_h0[R_SIZE] = {0}, y0[R_SIZE] = {0}, ss[R_SIZE] = {0}, y0b[R_BITS] = {0};

    for ( i = 0; i < R_BITS; ++i ) ts[i] = s[i];
    for ( i = 0; i < R_BITS; ++i ) s[i] = ts[(R_BITS - i) % R_BITS];

    for (i = 0; i < DV; ++i)
        h0[h0_compact[i] >> 3] |= ((uint8_t)1) << (h0_compact[i] & 7);
    convertBinaryToByte(ss, s, R_BITS);
    ntl_mod_inv(inv_h0, h0);
    ntl_mod_mul(y0, ss, inv_h0);
    convertByteToBinary(y0b, y0, R_BITS);
    for (i = 0; i < R_BITS; ++i)
        receivedword[i] = y0b[i];
    for (i = R_BITS; i < n; ++i)
        receivedword[i] = 0;
    

    for ( j = 0; j < n; ++j ){
        if ( (int)receivedword[j] < 0.5 )
            Obj[j] = log((double)(n - ERRN) / (double)ERRN);
        else
            Obj[j] = log((double)ERRN / (double)(n - ERRN));
    }

    for (i = 0; i < R_BITS; ++i)
    {
        for ( j = 0; j < DV; ++j ) Pj[i][j] = (R_BITS - h0_compact[DV - j - 1] + i) % R_BITS;
        for ( j = 0; j < DV; ++j ) Pj[i][DV + j] = (R_BITS - h1_compact[DV - j - 1] + i) % R_BITS + R_BITS;
    }

    int cnt[R_BITS * 2] = {0};
    for ( i = 0; i < R_BITS; ++i ){
        for ( j = 0; j < DV * 2; ++j ){
            POS[Pj[i][j]][cnt[Pj[i][j]]] = i;
            POS[Pj[i][j]][cnt[Pj[i][j]] + 1] = j;
            cnt[Pj[i][j]] += 2;
        }
    }

    for (i = 0; i < m; ++i)
    {
        for (k = 0; k < dc; ++k)
        {
            lambda[i][k] = 0;
            z[i][k] = receivedword[Pj[i][k]];
        }
    }

    while (itr_num < ITRN)
    {

        for (j = 0; j < n; ++j)
        {
            tmp = 0;
            for (i = 0; i < 2 * dv; i += 2)
            {
                tmp += z[POS[j][i]][POS[j][i + 1]] - (1 / ((double)MU)) * lambda[POS[j][i]][POS[j][i + 1]];
            }

            tmp = (tmp - (1 / ((double)MU)) * (Obj[j] + ALPHA2)) / (double)(dv - 2 * ALPHA2 / (double)MU);
            admmsol[j] = cutoff(tmp);
        }

        for (i = 0; i < m; ++i)
        {
            for (k = 0; k < dc; ++k)
            {
                v[i][k] = admmsol[Pj[i][k]] + lambda[i][k] / (double)MU;
            }
            for (k = 0; k < dc; ++k)
            {
                idx[k] = k;
                v_tmp[k] = v[i][k];
            }

            quicksort_dec(0, dc - 1, idx, v_tmp);

            for (k = 0; k < dc; ++k)
            {
                iidx[idx[k]] = k;
                v[i][k] = v_tmp[idx[k]];
            }
            err = proj(v, z_tmp, i);

            for (k = 0; k < dc; ++k)
            {
                z[i][k] = z_tmp[iidx[k]];
            }

            for (k = 0; k < dc; ++k)
            {
                lambda[i][k] += ((double)MU) * (admmsol[Pj[i][k]] - z[i][k]);
            }
        }

        ++itr_num;

        for (j = 0; j < n; ++j)
        {
            if (admmsol[j] > 0.5)
                admmbin_sol[j] = 1;
            else
                admmbin_sol[j] = 0;
        }

        int er[R_BITS * 2] = {0};
        // if ( itr_num == 50 ){
        //     FILE *pt = fopen("error", "r");
        //     for ( int i = 0; i < 2 * R_BITS; ++i ){
        //         // int t; if ( i < R_BITS ) t = (R_BITS - i) % R_BITS; else t = (2 * R_BITS - i) % R_BITS + R_BITS;
        //         fscanf( pt, "%d", er + i); 
        //     }
        //     for ( int i = 0; i < 2 * R_BITS; ++i ){
        //         admmbin_sol[i] = er[i] ^ (int)receivedword[i];
        //     }
        // }

        int total = 0;
        for (i = 0; i < R_BITS; ++i)
        {
            int ds = 0;
            for (j = 0; j < DV * 2; ++j)
            {
                ds ^= admmbin_sol[Pj[i][j]];
            }
            total += ds;
        }

        if (total == 0)
        {
            printf("%d\n", itr_num);
            for (i = 0; i < 2 * R_BITS; ++i)
                e[i] = admmbin_sol[i] ^ (int)receivedword[i];
            // for (i = 0; i < 2 * R_BITS; ++i ) printf( "%d ", e[i]); printf("\n");
            return 0;
        }
    }
    printf("%d\n", ITRN + 1);
    return 1; // FAILURE
}