#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"

#include "fill_zk_partial__order1.h"
/*******************************************************************/
/* mexFunction interface avec Matlab 

%%%%% BEWARE %%%%%
MATLAB code to compile it : mex Zk_MaxEnt_PsigmaK_partialorder1.c  get_zk_partial.c fill_zk_partial__order1.c -lm 
(not only mex MCMC_nsteps.c)

Polynoms are described by their coefficient in decreasing order (3X^2 + 8X is [3 8 0])
*/


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

	plhs[0]=mxCreateDoubleMatrix(1,N,mxREAL); 
	double *zk_l;
	zk_l = mxGetData(plhs[0]);
	
	fill_zk_partial__order1(exph_l, N, k, zk_l);
	return;		
}

