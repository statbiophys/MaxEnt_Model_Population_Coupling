#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"

#include "get_zk_partial.h"
#include "fill_zk_partial__order1.h"
/*******************************************************************/


void fill_zk_partial__order1(double *exph_l, int N, int k, double *zk_l){  /* mexPrintf("Naive \n"); */
	int i;
	double *exph_part_l;
	exph_part_l = mxMalloc((N-1)*sizeof(double));

	int r;

	if (k==0){
		for (i=0; i<N; i++){
			zk_l[i] = 1;
		}
	}
	else{	
		if (N==0){ mexErrMsgTxt("exph_l has length 0 but k>0"); }
		
		for (i=0; i<N; i++){	
			for (r=0; r<i; r++){
				exph_part_l[r] = exph_l[r];
			}
			for (r=i+1; r<N; r++){
				exph_part_l[r-1] = exph_l[r];
			}
			zk_l[i] = get_zk_part(exph_part_l, N-1, k-1);
			/*mexPrintf("i = %d , j = %d ,  zk_part = %f \n",i,j, zk_l[i+j*N]); */
		}
	}
	mxFree(exph_part_l);
	return;
}


/* supposedly 2 times less computation than naive solution, but actually longer */
void fill_zk_partial__order1_naive1(double *exph_l, int N, int k, double *zk_l){ 
	int i;
	double *exph_part_l;
	exph_part_l = mxMalloc((N-1)*sizeof(double));

	int r;

	if (k==0){
		for (i=0; i<N; i++){
			zk_l[i] = 1;
		}
	}
	else{	
		if (N==0){ mexErrMsgTxt("exph_l has length 0 but k>0"); }
		
		double **P_in;
		double **P_out;
		P_in = mxMalloc(N*sizeof(double*));
		P_out = mxMalloc(N*sizeof(double*));
		
		P_in[0] = mxMalloc(1*sizeof(double));
		P_in[0][0] = 1;
		int P_in_len;
		P_in_len = 1;

		int step_i;
		int i;
		int j;	
		int L_max;
		for (step_i=0; step_i<N; step_i++){
			/* convolve polynoms for i with i<step_i */
			P_in_len = step_i;
			L_max = min(P_in_len+1,k); /* here zk has k-1 terms */
			for (i=0;i<step_i;i++){ 
				/* mexPrintf("Here : step_i = %d , i = %d \n", step_i, i); */
				P_out[i] = temp_conv(P_in[i], P_in_len, exph_l[step_i]);
				mxFree(P_in[i]);
				P_in[i] = mxMalloc(L_max*sizeof(double));
				for (j=0;j<L_max;j++){
					P_in[i][j] = P_out[i][j];
				}
				mxFree(P_out[i]);
			}
			/* don't touch i for i=step_i*/
			if (step_i<(N-1)){ /* if it is not yet the last step, do i = step_i+1 */
				P_in_len = step_i+1;
				L_max = min(P_in_len+1,k); /* here zk has k-1 terms */
				P_out[step_i+1] = temp_conv(P_in[step_i], P_in_len, exph_l[step_i]); 
				P_in[step_i+1] = mxMalloc(L_max*sizeof(double));
				for (j=0;j<L_max;j++){ 
					P_in[step_i+1][j] = P_out[step_i+1][j];
				}
				mxFree(P_out[step_i+1]);
			}
		}	
		
		for (i=0; i<N; i++){
			zk_l[i] = P_in[i][k-1];
			mxFree(P_in[i]);	
		}
		mxFree(P_in);	
		mxFree(P_out);
	}
	return;
}

void fill_zk_partial__order1_naive2(double *exph_l, int N, int k, double *zk_l){ 
	int i;
	double *exph_part_l;
	exph_part_l = mxMalloc((N-1)*sizeof(double));

	int r;

	if (k==0){
		for (i=0; i<N; i++){
			zk_l[i] = 1;
		}
	}
	else{	
		if (N==0){ mexErrMsgTxt("exph_l has length 0 but k>0"); }
		
		double **P_in;
		double **P_out;
		P_in = mxMalloc(N*sizeof(double*));
		P_out = mxMalloc(N*sizeof(double*));
		
		int i;
		for (i=0; i<N; i++){ /* could be made smaller, depending on k : didn't work */
			P_in[i] = mxMalloc(N*sizeof(double));
			P_out[i] = mxMalloc(N*sizeof(double));			
		}
		P_in[0][0] = 1;
		int P_in_len;
		P_in_len = 1;

		int step_i;
		int j;	
		int L_max;
		for (step_i=0; step_i<N; step_i++){
			/* convolve polynoms for i with i<step_i */
			P_in_len = step_i;
			L_max = min(P_in_len+1,k); /* here zk has k-1 terms */
			for (i=0;i<step_i;i++){  /*   mexPrintf("Here : step_i = %d , i = %d \n", step_i, i); */
				fill_conv(P_in[i], P_in_len, exph_l[step_i], P_out[i]);
				
				
				for (j=0;j<L_max;j++){
					P_in[i][j] = P_out[i][j];
				}
				
			}
			/* don't touch i for i=step_i*/
			if (step_i<(N-1)){  /*   mexPrintf("Here last step for step_i = %d ", step_i);  /* if it is not yet the last step, do i = step_i+1 */
				P_in_len = step_i+1;
				L_max = min(P_in_len+1,k); /* here zk has k-1 terms */
				fill_conv(P_in[step_i], P_in_len, exph_l[step_i], P_out[step_i+1]); 
			
				for (j=0;j<L_max;j++){ 
					P_in[step_i+1][j] = P_out[step_i+1][j];
				}
			
			}
		}	
		
		for (i=0; i<N; i++){
			zk_l[i] = P_in[i][k-1]; 
		}
		for (i=0; i<N; i++){ /* mexPrintf("Free i = %d \n", i); */
			mxFree(P_in[i]);
			mxFree(P_out[i]);	
		}
		mxFree(P_in);	
		mxFree(P_out);
		
	}
	return;
}
