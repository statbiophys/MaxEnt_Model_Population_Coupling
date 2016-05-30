#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"

#include "get_zk_partial.h"
#include "fill_zk_partial__order2.h"
/*******************************************************************/


void fill_zk_partial__order2(double *exph_l, int N, int k, double *zk_m){
	int i;
	int j;	
	double *exph_part_l;
	exph_part_l = mxMalloc((N-2)*sizeof(double));

	int r;

	if (k==0){
		for (i=0; i<N; i++){
			for (j=i+1;j<N;j++){
				zk_m[i,j] = 1;
			}
		}
	}
	else{	
		if (N==0){ mexErrMsgTxt("exph_l has length 0 but k>0"); }
		
		for (i=0; i<N; i++){	
			for (j=i+1;j<N;j++){
				for (r=0; r<i; r++){
					exph_part_l[r] = exph_l[r];
				}
				for (r=i+1; r<j; r++){
					exph_part_l[r-1] = exph_l[r];
				}
				for (r=j+1; r<N; r++){
					exph_part_l[r-2] = exph_l[r];
				}
				zk_m[i+j*N] = get_zk_part(exph_part_l, N-2, k-2);
				zk_m[j+i*N] = zk_m[i+j*N];
				/*mexPrintf("i = %d , j = %d ,  zk_part = %f \n",i,j, zk_m[i+j*N]); */
			}
		}
	}
	mxFree(exph_part_l);

	return;
}
