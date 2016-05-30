#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"

#include "get_zk_partial.h"
#include "fill_zk_partial__order2.h"
#include "fill_zk_partial__order1.h"
/*******************************************************************/
/* mexFunction interface avec Matlab 

%%%%% BEWARE %%%%%
MATLAB code to compile it : mex stats_Psigma_PK__Pkalpha.c fill_zk_partial__order1.c fill_zk_partial__order2.c get_zk_partial.c -lm
(not only mex MCMC_nsteps.c)

Here Hess_m is the Hessian only for hi_l, since ak_l can be fitted analytically at each step
*/

void get_Hess_m_onlyh(double *Hess_m, double **expfki1_m, int N, int Kmax, int dimH, double Z, double **cZi_k, double *Pi_l){
		/* here only the hessian for parameters h is computed */
		/* deals with one k at a time : uses less memory */
		int k, i, j;
				
		/* mexPrintf("Beginning Hess_m 3 \n"); */
		for (i=0; i<N; i++){
			Hess_m[i+ i*dimH] = Pi_l[i] - Pi_l[i]*Pi_l[i];  /* Var P(i) */
		}	

		/* ---------- */
		
		for (i=0; i<N; i++){ 
			for (j=0; j<i; j++){
				Hess_m[i+ j*dimH] = 0;/* cov P(j)  P(i) */
			}
		}
		/* ----------- */

		double* cZij_kk;
		cZij_kk = mxMalloc(N*N*sizeof(double));
		for (k=1; k<Kmax; k++){ /* mexPrintf("Beginning Zij_k for k = %d \n", k+1); */ /* under 2 = 1+1, it is useless */  
			fill_zk_partial__order2(expfki1_m[k], N, k+1, cZij_kk); /*  cov P(k)  P(k)   */ /* k+1: index for C */
			for (i=0; i<N; i++){
				for (j=0; j<i; j++){
					cZij_kk[i+j*N] = cZij_kk[i+j*N]*expfki1_m[k][i]*expfki1_m[k][j];
					cZij_kk[j+i*N] = cZij_kk[i+j*N];
				}
				cZij_kk[i+i*N] = cZi_k[k][i]; /* used for correlations involving si*si */
			}

			for (i=0; i<N; i++){ 
				for (j=0; j<i; j++){
					Hess_m[i+ j*dimH] = Hess_m[i+ j*dimH]  + cZij_kk[i+j*N];   	
				}
			}
		}

		/* ---------- */
		for (i=0; i<N; i++){ 
			for (j=0; j<i; j++){
				Hess_m[i+ j*dimH] = Hess_m[i+ j*dimH]/Z;
				Hess_m[i+ j*dimH] = Hess_m[i+ j*dimH] - Pi_l[j]*Pi_l[i];
			}
		}		

		/* ----------- */
		
		/* mexPrintf("Beginning Hess_m 5 \n"); */
		for (i=0; i<dimH; i++){ /* completes the second half of the symmetric matric Hess_m */
			for (j=(i+1); j<dimH; j++){
				Hess_m[i+j*dimH] = Hess_m[j+i*dimH];
			}
		}
	
	mxFree(cZij_kk);
	return;
}



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	/* conversion of input arguments */
	double *Pk_l_targ;
	Pk_l_targ = mxGetData(prhs[0]);
	
	double *h_l;
	h_l = mxGetData(prhs[1]);

	int N;
	N = mxGetNumberOfElements(prhs[1]); /* number of neurons */
	if (N==0){ mexErrMsgTxt("h_l has length 0"); }

	int Kmax;
	Kmax = mxGetNumberOfElements(prhs[0]); /* maximum number of K (population rate) */
	if (Kmax==0){ mexErrMsgTxt("Pk_l_targ has length 0"); }

	/* output arguments */

	plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL); 
	double *Pi_l;
	Pi_l = mxGetData(plhs[0]);

	plhs[1]=mxCreateDoubleMatrix(Kmax,1,mxREAL);
	double *a_l_sol;
	a_l_sol = mxGetData(plhs[1]);
	
	int dimH = N;
	double *Hess_m;
	if (nlhs>2){
		plhs[2]=mxCreateDoubleMatrix(dimH,dimH,mxREAL); 
		Hess_m = mxGetData(plhs[2]);
	}

	/*mexPrintf("k = %d \n",k); */
	/* mexPrintf("here 1 \n"); */

	int i;
	int k;
	/* single exponential weights */
	double mean_h;
	mean_h = 0;
	for (i=0; i<N; i++){ mean_h = mean_h + h_l[i]; }
	mean_h = mean_h/N;

	double *exph_cent; 
	exph_cent = mxMalloc(N*sizeof(double)); /* mexPrintf("here 111 1\n"); */
	for (i=0; i<N; i++){
		exph_cent[i] = exp(h_l[i] - mean_h);
	}

	/*mexPrintf("here 2 \n");*/
	/* zk_l */
	double *zk_l_part;
	zk_l_part = mxMalloc(Kmax*sizeof(double));
	for (k=0; k<Kmax; k++){
		zk_l_part[k] = get_zk_part(exph_cent, N, k+1); /* k = (k+1) : index for C */
	}

	/* Z */
	/*double Z_part;
	Z_part = 1; /* 1 : set silence to enery 0 */
	/*for (k=0; k<Kmax; k++){
		Z_part = Z_part + zk_l_part[k];
	} */

	double P0_targ;
	P0_targ = 1;
	for (k=0; k<Kmax; k++){ P0_targ = P0_targ - Pk_l_targ[k]; }
	
	double Z;	
	Z = 1/P0_targ;

	/* double expa_sol; */
	double expa_k;
	double **expfki1_m; 
	expfki1_m = mxMalloc(Kmax*sizeof(double*)); /* mexPrintf("here 111 1\n"); */
	for (k=0; k<Kmax; k++){ /* here we note f_ki = a_k/k + h_i + k*g_k */
		
		expa_k = Pk_l_targ[k]/(zk_l_part[k]*P0_targ);

		a_l_sol[k] = log(expa_k) - (k+1)*mean_h;
		/*expa_sol = exp(a_l_sol[k]/(k+1));*/

		expfki1_m[k] = mxMalloc(N*sizeof(double)); /* mexPrintf("here 112 \n"); */
		for (i=0; i<N; i++){
			expfki1_m[k][i] = exp( a_l_sol[k]/(k+1) + h_l[i] );
		}	
	}

	/*mexPrintf("here 3 \n");*/
	/* zk partials */
	double** cZi_k;
	cZi_k = mxMalloc(Kmax*sizeof(double*));
	for (k=0; k<Kmax; k++){
		cZi_k[k] = mxMalloc(N*sizeof(double*));
		fill_zk_partial__order1(expfki1_m[k], N, k+1, cZi_k[k]); /* k + 1  : index for k in C */
		for (i=0; i<N;i++){
			cZi_k[k][i] = cZi_k[k][i]*expfki1_m[k][i];		
		}
	}
	
	for (i=0; i<N; i++){
		Pi_l[i] = 0;
		for (k=0; k<Kmax; k++){      /*mexPrintf("here 31, i = %d, k = %d \n",i,k);*/
			Pi_l[i] = Pi_l[i] + cZi_k[k][i];
			/*mexPrintf("i = %d , k = %d ,  coef = %f , Z = %f \n",i+1,k+1, coef, Z); */
		}
		Pi_l[i] = Pi_l[i]/Z;
	}

	/*mexPrintf("here 4 \n");*/

	if (nlhs>2){/* Hessian */
		/* get_Hess_m_fast(Hess_m, expfki1_m, Kmax, N, dimH, Z, cZi_k, Pk_l, Pi_l); */
		get_Hess_m_onlyh(Hess_m, expfki1_m, N, Kmax, dimH, Z, cZi_k, Pi_l);
	}
	 /* Free array adresses */
	for (k=0; k<Kmax; k++){
		mxFree(expfki1_m[k]);
	}
	mxFree(expfki1_m);
	for (k=0; k<Kmax; k++){
		mxFree(cZi_k[k]);
	}
	mxFree(cZi_k);
	mxFree(zk_l_part);
	mxFree(exph_cent);
	return;			
}







