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
MATLAB code to compile it : mex stats_Psigma_PK_meansigmaK__alphabetagamma.c get_zk_partial.c fill_zk_partial__order1.c fill_zk_partial__order2.c  -lm 
(not only mex MCMC_nsteps.c)
*/


void get_Hess_m(double *Hess_m, int Kmax, int N, int dimH, double **expfki_m, double *expak_l, double Z, double **cZi_k, double *Pk_l, double *Pi_l, double *mean_Ki_l){
		int k, i, j;
		double** cZij_k;
		cZij_k = mxMalloc(Kmax*sizeof(double*)); 
		for (k=1; k<Kmax; k++){ /* mexPrintf("Beginning Zij_k for k = %d \n", k+1); */ /* under 2 = 1+1, it is useless */  
			cZij_k[k] = mxMalloc(N*N*sizeof(double));
			fill_zk_partial__order2(expfki_m[k], N, k+1, cZij_k[k]); /*  cov P(k)  P(k)   */ /* k+1: index for C */
			for (i=0; i<N; i++){
				for (j=0; j<i; j++){
					cZij_k[k][i+j*N] = cZij_k[k][i+j*N]*expak_l[k]*expfki_m[k][i]*expfki_m[k][j];
					cZij_k[k][j+i*N] = cZij_k[k][i+j*N];
				}
				cZij_k[k][i+i*N] = cZi_k[k][i]; /* used for correlations involving si*si */
			}
		}

		for (k=0; k<Kmax; k++){ 
			Hess_m[k+k*dimH] = Pk_l[k] - Pk_l[k]*Pk_l[k];  /*  var P(k)   */
		}
		for (k=0; k<Kmax; k++){ /* for clarity, k should be called l here */
			for (j=0; j<k; j++){
				Hess_m[k+j*dimH] = - Pk_l[k]*Pk_l[j];	 /*  cov P(k)  P(k')   */
			}
		}
		
		/* mexPrintf("Beginning Hess_m 1 \n"); */
		for (k=0; k<Kmax; k++){
			for (j=0; j<N; j++){
				Hess_m[(Kmax+j)+ k*dimH] = cZi_k[k][j];  /* cov P(k) P(j) */	
				Hess_m[(Kmax+j)+ k*dimH] = Hess_m[(Kmax+j)+ k*dimH]/Z;
				Hess_m[(Kmax+j)+ k*dimH] = Hess_m[(Kmax+j)+ k*dimH] - Pk_l[k]*Pi_l[j];
	
				Hess_m[(Kmax+N+j) + k*dimH] = (k+1)*cZi_k[k][j];  /* cov P(k)  < K sj > */	
				Hess_m[(Kmax+N+j) + k*dimH] = Hess_m[(Kmax+N+j) + k*dimH]/Z;
				Hess_m[(Kmax+N+j) + k*dimH] = Hess_m[(Kmax+N+j) + k*dimH] - Pk_l[k]*mean_Ki_l[j];
			}
		}	

		/* mexPrintf("Beginning Hess_m 2 \n"); */
		for (i=0; i<N; i++){
			for (j=0; j<N; j++){
				if ( i==j){  Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = cZi_k[0][i];  }   /* cov P(j) < Ksi > */
				else{	Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = 0;  } 
				for (k=1; k<Kmax; k++){  	Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH]   + (k+1)*cZij_k[k][i+j*N];	}
				Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH]/Z;
				Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] - Pi_l[j]*mean_Ki_l[i];
			}
		}

		/* mexPrintf("Beginning Hess_m 3 \n"); */
		for (i=0; i<N; i++){
			Hess_m[(Kmax+i)+ (Kmax+i)*dimH] = Pi_l[i] - Pi_l[i]*Pi_l[i];  /* Var P(i) */
		
			Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = 0; /* Var < K si >  */
			for (k=0; k<Kmax; k++){  	Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH]   +  (k+1)*(k+1)*cZi_k[k][i]; 	}
			Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH]/Z;
			Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] - mean_Ki_l[i]*mean_Ki_l[i];
		}

		/* mexPrintf("Beginning Hess_m 4 \n"); */
		for (i=0; i<N; i++){ 
			for (j=0; j<i; j++){
				Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = 0;/* cov P(j)  P(i) */
				for (k=1; k<Kmax; k++){  	Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(Kmax+i)+(Kmax+j)*dimH]  + cZij_k[k][i+j*N];    }
				Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(Kmax+i)+ (Kmax+j)*dimH]/Z;
				Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(Kmax+i)+ (Kmax+j)*dimH] - Pi_l[j]*Pi_l[i];

				Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = 0;/* cov < K sj >  < K si > */
				for (k=1; k<Kmax; k++){  	Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH]  +  (k+1)*(k+1)*cZij_k[k][i+j*N];	}
				Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH]/Z;
				Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] - mean_Ki_l[i]*mean_Ki_l[j];
			}
		}
		
		/* mexPrintf("Beginning Hess_m 5 \n"); */
		for (i=0; i<dimH; i++){ /* completes the second half of the symmetric matric Hess_m */
			for (j=(i+1); j<dimH; j++){
				Hess_m[i+j*dimH] = Hess_m[j+i*dimH];
			}
		}
	
		/* mexPrintf("Free 1 \n"); */
		for (k=1; k<Kmax; k++){
			mxFree(cZij_k[k]);
		}
		mxFree(cZij_k);
	return;
}

void get_Hess_m_fast(double *Hess_m, double **expfki1_m, int Kmax, int N, int dimH, double Z, double **cZi_k, double *Pk_l, double *Pi_l, double *mean_Ki_l){
		/* deals with one k at a time : uses less memory */
		int k, i, j;
		for (k=0; k<Kmax; k++){ 
			Hess_m[k+k*dimH] = Pk_l[k] - Pk_l[k]*Pk_l[k];  /*  var P(k)   */
		}
		for (k=0; k<Kmax; k++){ /* for clarity, k should be called l here */
			for (j=0; j<k; j++){
				Hess_m[k+j*dimH] = - Pk_l[k]*Pk_l[j];	 /*  cov P(k)  P(k')   */
			}
		}
		
		/* mexPrintf("Beginning Hess_m 1 \n"); */
		for (k=0; k<Kmax; k++){
			for (j=0; j<N; j++){
				Hess_m[(Kmax+j)+ k*dimH] = cZi_k[k][j];  /* cov P(k) P(j) */	
				Hess_m[(Kmax+j)+ k*dimH] = Hess_m[(Kmax+j)+ k*dimH]/Z;
				Hess_m[(Kmax+j)+ k*dimH] = Hess_m[(Kmax+j)+ k*dimH] - Pk_l[k]*Pi_l[j];
	
				Hess_m[(Kmax+N+j) + k*dimH] = (k+1)*cZi_k[k][j];  /* cov P(k)  < K sj > */	
				Hess_m[(Kmax+N+j) + k*dimH] = Hess_m[(Kmax+N+j) + k*dimH]/Z;
				Hess_m[(Kmax+N+j) + k*dimH] = Hess_m[(Kmax+N+j) + k*dimH] - Pk_l[k]*mean_Ki_l[j];
			}
		}

		/* mexPrintf("Beginning Hess_m 3 \n"); */
		for (i=0; i<N; i++){
			Hess_m[(Kmax+i)+ (Kmax+i)*dimH] = Pi_l[i] - Pi_l[i]*Pi_l[i];  /* Var P(i) */
		
			Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = 0; /* Var < K si >  */
			for (k=0; k<Kmax; k++){  	Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH]   +  (k+1)*(k+1)*cZi_k[k][i]; 	}
			Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH]/Z;
			Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+i)*dimH] - mean_Ki_l[i]*mean_Ki_l[i];
		}	

		/* ---------- */
		for (i=0; i<N; i++){
			for (j=0; j<N; j++){
				if ( i==j){  Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = cZi_k[0][i];  }   /* cov P(j) < Ksi > */
				else{	Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = 0;  } 
			}
		}

		for (i=0; i<N; i++){ 
			for (j=0; j<i; j++){
				Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = 0;/* cov P(j)  P(i) */
				Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = 0;/* cov < K sj >  < K si > */
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
				for (j=0; j<N; j++){ 
					Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH]   + (k+1)*cZij_kk[i+j*N];
				}
			}

			for (i=0; i<N; i++){ 
				for (j=0; j<i; j++){
					Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(Kmax+i)+(Kmax+j)*dimH]  + cZij_kk[i+j*N];   
					Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH]  +  (k+1)*(k+1)*cZij_kk[i+j*N];	
				}
			}
		}

		/* ---------- */
		for (i=0; i<N; i++){
			for (j=0; j<N; j++){ 
				Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH]/Z;
				Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (Kmax+j)*dimH] - Pi_l[j]*mean_Ki_l[i];
			}
		}

		for (i=0; i<N; i++){ 
			for (j=0; j<i; j++){
				Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(Kmax+i)+ (Kmax+j)*dimH]/Z;
				Hess_m[(Kmax+i)+ (Kmax+j)*dimH] = Hess_m[(Kmax+i)+ (Kmax+j)*dimH] - Pi_l[j]*Pi_l[i];

				Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH]/Z;
				Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] = Hess_m[(N+Kmax+i)+ (N+Kmax+j)*dimH] - mean_Ki_l[i]*mean_Ki_l[j];
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
	double *a_l;
	a_l = mxGetData(prhs[0]);
	
	double *h_l;
	h_l = mxGetData(prhs[1]);

	double *g_l;
	g_l = mxGetData(prhs[2]);
	
	int N;
	N = mxGetNumberOfElements(prhs[1]); /* number of neurons */
	if (N==0){ mexErrMsgTxt("h_l has length 0"); }

	int Kmax;
	Kmax = mxGetNumberOfElements(prhs[0]); /* maximum number of K (population rate) */
	if (Kmax==0){ mexErrMsgTxt("a_l has length 0"); }

	/* output arguments */
	plhs[0]=mxCreateDoubleMatrix(1,Kmax,mxREAL); 
	double *Pk_l;
	Pk_l = mxGetData(plhs[0]);

	plhs[1]=mxCreateDoubleMatrix(1,N,mxREAL); 
	double *Pi_l;
	Pi_l = mxGetData(plhs[1]);

	plhs[2]=mxCreateDoubleMatrix(1,N,mxREAL); 
	double *mean_Ki_l;
	mean_Ki_l = mxGetData(plhs[2]);

	int dimH = Kmax+2*N;
	double *Z_out;
	double *Hess_m;
	if (nlhs>3){
		plhs[3]=mxCreateDoubleMatrix(dimH,dimH,mxREAL); 
		Hess_m = mxGetData(plhs[3]);

		plhs[4]=mxCreateDoubleMatrix(1,1,mxREAL); 
		Z_out = mxGetData(plhs[4]);
	}

	/*mexPrintf("k = %d \n",k); */
	/* mexPrintf("here 1 \n"); */

	int i;
	int k;
	/* single exponential weights */
	double **expfki1_m; 
	expfki1_m = mxMalloc(Kmax*sizeof(double*)); /* mexPrintf("here 111 1\n"); */
	for (k=0; k<Kmax; k++){ /* here we note f_ki = a_k/k + h_i + k*g_k */
		expfki1_m[k] = mxMalloc(N*sizeof(double)); /* mexPrintf("here 112 \n"); */
		for (i=0; i<N; i++){
			expfki1_m[k][i] = exp(a_l[k]/(k+1) + h_l[i] + (k+1)*g_l[i]);
		}
	}

	/*mexPrintf("here 2 \n");*/
	/* zk_l */
	double *zk_l;
	zk_l = mxMalloc(Kmax*sizeof(double));
	for (k=0; k<Kmax; k++){
		zk_l[k] = get_zk_part(expfki1_m[k], N, k+1); /* k = (k+1) : index for C */

	}

	/* Z */
	double Z;
	Z = 1; /* 1 : set silence to enery 0 */
	for (k=0; k<Kmax; k++){
		Z = Z + zk_l[k];
	}
	if (nlhs>3){
		Z_out[0] = Z; /* for output */
	}

	for (k=0; k<Kmax; k++){
		Pk_l[k] = zk_l[k]/Z;
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
		mean_Ki_l[i] = 0;
		Pi_l[i] = 0;
		for (k=0; k<Kmax; k++){      /*mexPrintf("here 31, i = %d, k = %d \n",i,k);*/
			Pi_l[i] = Pi_l[i] + cZi_k[k][i];
			mean_Ki_l[i] = mean_Ki_l[i] + (k+1)*cZi_k[k][i]; /* k+1: index for k in C */
			/*mexPrintf("i = %d , k = %d ,  coef = %f , Z = %f \n",i+1,k+1, coef, Z); */
		}
		Pi_l[i] = Pi_l[i]/Z;
		mean_Ki_l[i] = mean_Ki_l[i]/Z;
	}

	/*mexPrintf("here 4 \n");*/

	if (nlhs>3){/* Hessian */
		get_Hess_m_fast(Hess_m, expfki1_m, Kmax, N, dimH, Z, cZi_k, Pk_l, Pi_l, mean_Ki_l);
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
	mxFree(zk_l);
	return;			
}







