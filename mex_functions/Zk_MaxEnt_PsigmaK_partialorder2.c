#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"

#include "fill_zk_partial__order2.h"
/*******************************************************************/
/* mexFunction interface avec Matlab 

%%%%% BEWARE %%%%%
MATLAB code to compile it : mex Zk_MaxEnt_PsigmaK_partialorder2.c fill_zk_partial__order2.c -lm  
(not only mex MCMC_nsteps.c)

Polynoms are described by their coefficient in decreasing order (3X^2 + 8X is [3 8 0])
*/

int min(int x, int y)
{
  return (x < y) ? x : y;
}

double* temp_conv(double *P1, int P1_len, double a){
	/* P1 is convolved with {1, a} */
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

double get_zk_part(double *exph_part_l, int N_part, int k_part){

	double *P_in;
	P_in = mxMalloc(2*sizeof(double));
	P_in[0] = 1;
	P_in[1] = exph_part_l[0];

	double *P_out;
	int i;
	int j;
	int L_max;
	int P_in_len = 2;

	for (i=1;i<N_part;i++){
		/*mexPrintf("Here : %f\n",P_in[P_in_len-1]);*/
		P_out = temp_conv(P_in, P_in_len, exph_part_l[i]);
		L_max = min(P_in_len+1,k_part+1);
		/*mexPrintf("P_in_len = %d \n", P_in_len);
		mexPrintf("L_max = %d \n",L_max);*/
		mxFree(P_in);
		P_in = mxMalloc(L_max*sizeof(double));
		for (j=0;j<L_max;j++){
			P_in[j] = P_out[j];
		}
		P_in_len = L_max;
		mxFree(P_out);
	}
	
	double zk_part;
	zk_part = P_in[k_part];
	mxFree(P_in);
	return zk_part;

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

	/* plhs[0]=mxCreateNumericMatrix(N,1,mxUINT32_CLASS,mxREAL); */
	plhs[0]=mxCreateDoubleMatrix(N,N,mxREAL); 
	double *zk_m;
	zk_m = mxGetData(plhs[0]);

	/*mexPrintf("k = %d \n",k); 
	mexPrintf("N = %d \n",N); */
	
	fill_zk_partial__order2(exph_l, N, k,zk_m); /* does not fill the diagonal */
	return;		
}

