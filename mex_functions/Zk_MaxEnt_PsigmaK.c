#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"
/*******************************************************************/
/* 
Computes Zk coefficients

%%%%% BEWARE %%%%%
MATLAB code to compile it : mex Zk_MaxEnt_PsigmaK.c -lm 
(not sufficient to write "mex MCMC_nsteps.c" )

Polynoms are described by their coefficient in decreasing order (3X^2 + 8X is [3 8 0])
*/

int min(int x, int y)
{
  return (x < y) ? x : y;
}

double* temp_conv(double *P1, int P1_len, double a){
    /* 
    Convolves polynom P_in with polynom [1 a]
	P1 is convolved with {1, a} 
    */
	double *P2;
	P2 = mxMalloc((P1_len+1)*sizeof(double));

	P2[0] = P1[0];

	int i;
	for (i=1;i<P1_len;i++){
		P2[i]=P1[i] + P1[i-1]*a;	
	}

	P2[P1_len] = P1[P1_len-1]*a;
	return P2;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	/* conversion of input arguments */
	double *exph_l;
	exph_l = mxGetData(prhs[0]);

	double *k_d;
	k_d = mxGetData(prhs[1]);
	int k;
	k = (int) k_d[0];

	int N;
	N = mxGetNumberOfElements(prhs[0]);

	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
	double *zk;
	zk = mxGetData(plhs[0]);

	if (k==0){
		zk[0] = 1;
	}
	else{	
		if (N==0){ mexErrMsgTxt("exph_l has length 0 but k>0"); }

		double *P_in;
		P_in = mxMalloc(2*sizeof(double));
		P_in[0] = 1;
		P_in[1] = exph_l[0];

		double *P_out;
		int i;
		int j;
		int L_max;
		int P_in_len = 2;
		/* 1 because we have already used exph_l[0] in P_in */
		for (i=1;i<N;i++){
			P_out = temp_conv(P_in, P_in_len, exph_l[i]);
			L_max = min(P_in_len+1,k+1);
			mxFree(P_in);
			P_in = mxMalloc(L_max*sizeof(double));
			for (j=0;j<L_max;j++){
				P_in[j] = P_out[j];
			}
			P_in_len = L_max;
			mxFree(P_out);
		}
		zk[0] = P_in[k];
		mxFree(P_in);
	}
	
	return;		
}

