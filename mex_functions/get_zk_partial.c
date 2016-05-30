#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"

#include "get_zk_partial.h" /* header file */

int min(int x, int y)
{
  return (x < y) ? x : y;
}

void fill_conv(double *P1, int P1_len, double a, double* P2){
	/* P1 is convolved with {1, a} */
	/* P2 must have at least length length(P1) + 1 */
	P2[0] = P1[0];

	int i;
	for (i=1;i<P1_len;i++){
		P2[i]=P1[i] + P1[i-1]*a;	
	}

	P2[P1_len] = P1[P1_len-1]*a;
	return;
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

double get_zk_part_naive(double *exph_part_l, int N_part, int k_part){
	/* N_part : number of elements in exph_part_l
	 k_part : numbe of term in each product */
	if (k_part == 0){
		return 1;
	}
	else{
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
			/* mexPrintf("Here : %f\n",P_in[P_in_len-1]); */
			P_out = temp_conv(P_in, P_in_len, exph_part_l[i]);
			L_max = min(P_in_len+1,k_part+1);
			/* mexPrintf("P_in_len = %d \n", P_in_len); */
			/* mexPrintf("L_max = %d \n",L_max); */
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
}

double get_zk_part(double *exph_part_l, int N_part, int k_part){
	/* N_part : number of elements in exph_part_l
	 k_part : numbe of term in each product */
	if (k_part == 0){
		return 1;
	}
	else{
		double *P_in;
		double *P_out;
		P_in = mxMalloc((k_part+1)*sizeof(double));
		P_out = mxMalloc((k_part+2)*sizeof(double));
		P_in[0] = 1;
		P_in[1] = exph_part_l[0];

		int i;
		int j;
		int L_max;
		int P_in_len = 2;

		for (i=1;i<N_part;i++){
			/* mexPrintf("Here : %f\n",P_in[P_in_len-1]); */
			fill_conv(P_in, P_in_len, exph_part_l[i], P_out);
			L_max = min(P_in_len+1,k_part+1);
			/* mexPrintf("P_in_len = %d \n", P_in_len); */
			/* mexPrintf("L_max = %d \n",L_max); */
			for (j=0;j<L_max;j++){
				P_in[j] = P_out[j];
			}
			P_in_len = L_max;
		}
	
		double zk_part;
		zk_part = P_in[k_part];
		mxFree(P_in);
		mxFree(P_out);
		return zk_part;
	}
}
