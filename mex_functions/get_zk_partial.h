#ifndef get_zk_partial_HEADER_FILE
#define get_zk_partial_HEADER_FILE

double get_zk_part(double *exph_part_l, int N_part, int k_part);
double get_zk_part_naive(double *exph_part_l, int N_part, int k_part);
double* temp_conv(double *P1, int P1_len, double a);
void fill_conv(double *P1, int P1_len, double a, double *P2);

#endif
