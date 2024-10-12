#include "ext_decoder.h"
#include "decode.h"
#include "utilities.h"
#include "admm_decoder.h"

#include "kem.h"
#include "sampling.h"

#include "ring_buffer.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

const int MU = 6, ALPHA2 = 15;
const int dc = DV * 2, dv = DV;
const int ITRN = 10000

void makeObj(double *Obj, double *receivedword, double sigma, int n){
   int j;
   double Pr0, Pr1, s;
   
   for(j = 0; j < n; j++){
     if((int)(receivedword[j])<0.5){
       Obj[j] = log((double)(n-ERRN)/(double)ERRN);
     }
     else{
       Obj[j] = log(((double)ERRN)/((double)(n-ERRN)));
     }
   }
}

double cutoff(double x){
    if(x>1){ return(1);}
    else if(x<0){ return(0);}
    else{ return(x);}
}

double sign(double x){
  if(x>0){
    return(1);
  }else{
    return(-1);
  }
}

int hardd(double x){
  if(x>=0){
    return(0);
  }else{
    return(1);
  }
}

int abs_(int x){
  if(x>=0){return(x);
  }
  else{
    return(-x);
  }     
}

int MemberQ(int j, int i){
  int k;
  for(k=0;k<dc;k++){
    //printf("%d ",Pj[i][k]);
    if(j==Pj[i][k]){
      return(k);
    }
  }
  return(-1);
}

void quicksort(int first, int last, int *idx, double *dat)
{
  int i, j, tmp;
  double x;
  
  x = dat[idx[(first+last)/2]];
  i = first;  j = last;
  for ( ; ; ) {
    while (dat[idx[i]] < x) i++;
    while (x < dat[idx[j]]) j--;
    if (i >= j) break;
    tmp = idx[i];  idx[i] = idx[j];  idx[j] = tmp;
    i++;  j--;
  }
  if (first  < i - 1) quicksort(first , i - 1, idx, dat);
  if (j + 1 < last) quicksort(j + 1, last, idx, dat);
}

void quicksort_dec(int first, int last, int *idx, double *dat)
{
  int i, j, tmp;
  double x;
  
  x = dat[idx[(first+last)/2]];
  i = first;  j = last;
  for ( ; ; ) {
    while (dat[idx[i]] > x) i++;
    while (x > dat[idx[j]]) j--;
    if (i >= j) break;
    tmp = idx[i];  idx[i] = idx[j];  idx[j] = tmp;
    i++;  j--;
  }
  if (first  < i - 1) quicksort_dec(first , i - 1, idx, dat);
  if (j + 1 < last) quicksort_dec(j + 1, last, idx, dat);
}

int proj(double **v, double *z, int row){
  int i,j,k,beta_max,r,posi,L,R,M,plus;
  double tmp, frz, frz_pre, frzero, beta_opt, z_zero_sum;
  
  double z_hat[dc], fr[dc], Beta_Set[dc * 4], idx_p[dc * 4];

  /* r */
  tmp = 0;
  for(k=0;k<dc;k++){
    if(v[row][k]>1){
      tmp += 1;
      z_hat[k] = 1;
    }
    else if(v[row][k]<0){
      tmp += 0;
      z_hat[k] = 0;
    }
    else{
      tmp += v[row][k];
      z_hat[k] = v[row][k];
    }
  }
  
  r = (int)tmp;
  if(r%2==1) r -= 1; 
  
  if(r==dc||r==dc-1){
    for(k=0;k<dc;k++){z[k] = z_hat[k];}
    return(1);
  }

   /* fr */
  for(k=0;k<dc;k++){
    if(k<=r){
      fr[k] = 1;
    }else{
      fr[k] = -1;
    }
  }
 
  tmp=0;
  for(k=0;k<dc;k++){tmp += fr[k]*z_hat[k];}
  if(tmp<=r+0.00001){                                    //13:
    for(k=0;k<dc;k++){
      z[k] = z_hat[k];
    }
    return(1);
  }
 
  //****************************************************************
  
  tmp = 0.5*(v[row][r]-v[row][r+1])+0.00000000001;
  Beta_Set[0] = 0;

  plus = 1;
  for(k=0;k<dc;k++){                                       //16:
    if(k>0&&v[row][k]>v[row][k-1]-0.00000000001&&v[row][k]<v[row][k-1]+0.00000000001) continue;
    if(k<=r){
      if((v[row][k]-1>0)&&((v[row][k]-1)<tmp)){
	Beta_Set[plus]= v[row][k] - 1;
	plus++;
      }
    }else{
      if((-v[row][k]>0)&&((-v[row][k])<tmp)){
	Beta_Set[plus] = -v[row][k];
	plus++;
      }
    }
  }
  
  for(k=0;k<dc;k++){                                       //16:
    //if(k>0&&v[row][k]>v[row][k-1]-0.00000000001&&v[row][k]<v[row][k-1]+0.00000000001) continue;
    if(k<=r){
      if((v[row][k]>0)&&(v[row][k]<tmp)){
	Beta_Set[plus] = v[row][k];
	plus++;
      }	
    }else{
      if((1-v[row][k]>0)&&((1-v[row][k])<tmp)){
	Beta_Set[plus] = 1-v[row][k];
	plus++;
      }
    }
  }

  Beta_Set[plus] = tmp;

  for(k=0;k<2*dc;k++){ idx_p[k] = k;}
  quicksort(0,plus,idx_p,Beta_Set); 

  posi = 0;  

  L = posi;
  //R = 2*dc-1;
  R = plus;
  M = (L + R)/2;  
  while(R-L!=1){// Bi-section Method
    for(j=0;j<=r;j++){ z_hat[j] = v[row][j] - Beta_Set[idx_p[M]]; }
    for(j=r+1;j<dc;j++){ z_hat[j] = v[row][j] + Beta_Set[idx_p[M]]; }
    tmp = 0;
    for(j=0;j<dc;j++){ tmp += fr[j]*cutoff(z_hat[j]); }
    if(tmp<r){
      R = M;
    }else{
      L = M;
    }
    M = (L + R)/2;
  }

  for(j=0;j<=r;j++){ z_hat[j] = v[row][j] - Beta_Set[idx_p[R]]; }
  for(j=r+1;j<dc;j++){ z_hat[j] = v[row][j] + Beta_Set[idx_p[R]]; }
  tmp = 0;
  for(j=0;j<dc;j++){ tmp += fr[j]*cutoff(z_hat[j]); }
  frz = tmp;
  for(j=0;j<=r;j++){ z_hat[j] = v[row][j] - Beta_Set[idx_p[L]]; }
  for(j=r+1;j<dc;j++){ z_hat[j] = v[row][j] + Beta_Set[idx_p[L]]; }
  tmp = 0;
  for(j=0;j<dc;j++){ tmp += fr[j]*cutoff(z_hat[j]); }
  frz_pre = tmp;
  beta_opt = (r-frz)*(Beta_Set[idx_p[R]]-Beta_Set[idx_p[L]])/(frz-frz_pre) + Beta_Set[idx_p[R]];
  
  for(j=0;j<=r;j++){ z[j] = cutoff(v[row][j] - beta_opt); } //optimal sol
  for(j=r+1;j<dc;j++){ z[j] = cutoff(v[row][j] + beta_opt);}
  return(1);
  
  printf("error: something is wrong?\n");
  exit(1);
 
}

int ADMM_decoder(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    const int n = R_BITS * 2, m = R_BITS;

    int i,j,k, itr_num=1, err, posit;
    double **v, **z, **lambda;
    double tmp, tmp2, prev=0;
    int *idx, *iidx;
    double *v_tmp, *z_tmp;
    z = (double **)malloc(m*sizeof(double *));
    for(i = 0; i < m; i++) z[i] = (double *)calloc(dc,sizeof(double));
    lambda = (double **)malloc(m*sizeof(double *));
    for(i = 0; i < m; i++) lambda[i] = (double *)calloc(dc,sizeof(double));
    v = (double **)malloc(m*sizeof(double *));
    for(i = 0; i < m; i++) v[i] = (double *)calloc(dc,sizeof(double));
    idx = (int *)calloc(dc,sizeof(int));
    iidx = (int *)calloc(dc,sizeof(int));
    v_tmp = (double *)calloc(dc,sizeof(double));
    z_tmp = (double *)calloc(dc,sizeof(double));
    double *admmsol; int *admmbin_sol;
    admmsol = (double *)malloc(n*sizeof(double));
    admmbin_sol = (int *)malloc(n*sizeof(int));
    double *Obj, *receivedword; double sigma = 0;
    Obj = (double *)malloc(n*sizeof(double));
    receivedword = (double *)malloc(n*sizeof(double));
    for ( i = 0; i < n; ++i ) receivedword[i] = 0;
    makeObj(Obj, receivedword, sigma, n);
    int **Pj, **POS;
    Pj = (int **)malloc(m * sizeof(int*));
    for ( i = 0; i < m; ++i ) Pj[i] = (int *)calloc(dc, sizeof(int));
    POS = (int **)malloc(m * sizeof(int*));
    for ( i = 0; i < m; ++i ) POS[i] = (int *)calloc(dc, sizeof(int));
    
    for ( int i = 0; i < R_BITS; ++i ){
        for ( int j = 0; j < DV; ++j ){
            Pj[i][j] = (R_BITS - (h0_compact[j] + i) % R_BITS) % R_BITS;
            Pj[i][j + DV] = (R_BITS - (h1_compact[j] + i) % R_BITS) % R_BITS + R_BITS;
        }
    }
    for ( int i1 = 0; i1 < n; i1++ ){
        int i2 = 0;
        for ( i = 0; i < m; i++ ){
            if(MemberQ(i1,i)>=0){
                POS[i1][i2] = i;
                POS[i1][i2+1] = MemberQ(i1,i);
                i2 = i2 + 2;
            }
        }
    }

    for ( i = 0; i < m; ++i ){
        for ( k = 0; k < dc; ++k ){
            lambda[i][k] = 0;
            z[i][k] = 0; // !
        }
    }

    while (itr_num <= ITRN){
        for ( j = 0; j < n; ++j ){
            tmp = 0;
            for ( i = 0; i < 2 * dv; i += 2 ){
                tmp += z[POS[j][i]][POS[j][i+1]]-(1/((double)MU))*lambda[POS[j][i]][POS[j][i+1]];
            }

            tmp = (tmp - (1/((double)MU))*(Obj[j]+ALPHA2))/(double)(Nv[j]-2*ALPHA2/(double)MU);
            admmsol[j]= cutoff(tmp);
        }

        for ( i = 0; i < m; ++i ){
            for ( k = 0; k < dc; ++k ){
                v[i][k] = admmsol[Pj[i][k]] + lambda[i][k]/(double)MU;
            }
            for ( k = 0; k < dc; ++k ){
                idx[k] = k;
                v_tmp[k] = v[i][k];
            }

            quicksort(0, dc - 1, idx, v_tmp);

            for ( k = 0; k < dc; ++k ){
                iidx[idx[k]] = k;
                v[i][k] = v_tmp[idx[k]];
            }

            err = proj(v, z_tmp, i);
            
            for ( k = 0; k < dc; ++k ){
                z[i][k] = z_tmp[iidx[k]];
            }

            for ( k = 0; k < dc; ++k ){
                lambda[i][k] += ((double)MU)*(admmsol[Pj[i][k]] - z[i][k]);
            }
        }

        ++itr_num;

        for ( j = 0; j < n; ++j ){
            if ( admmsol[j] > 0.5 ) admmbin_sol[j] = 1;
            else admmbin_sol[j] = 0;
        }

        int total = 0;
        for ( i = 0; i < R_BITS; ++i ){
            int ds = s[i];
            for ( j = 0; j < DV * 2; ++j ){
                ds ^= admmbin_sol[Pj[i][j]];
            }
            total += ds;
        }

        if ( total == 0 ){
            for(i = 0; i < m; i++){free(z[i]), free(lambda[i]), free(v[i]);}
            free(z),free(lambda), free(v);
            free(idx), free(iidx), free(v_tmp), free(z_tmp), free(admmsol), free(admmbin_sol);
            printf("%d\n", itr_num);
            return 0;
        }

    }

    for(i = 0; i < m; i++){free(z[i]), free(lambda[i]), free(v[i]);}
    free(z),free(lambda), free(v);
    free(idx), free(iidx), free(v_tmp), free(z_tmp), free(admmsol), free(admmbin_sol);
    printf("%d\n", ITRN + 1);
    return 1; // FAILURE
}