function [ Pk_l, Pi_l, mKi_l, Pki_m, P0, Kmax_new, Nki_emp_m ] = regularize( spikewords , lambda_reg, Kmax_add)
%REGULARIZE
% the regularized Kmax is set to max(K) + Kmax_add
[Nt, Nneu] = size(spikewords);

%% Prior
[ Pki_indep_m, ~ ] = probas_indep_Ising( spikewords );

%% count empirical event
K_l = sum(spikewords,2);
Kmax = max(K_l);

Nki_emp_m = zeros(Kmax, Nneu);
for k = 1:Kmax
    keep = (K_l == k);
    Nki_emp_m(k,:) = sum(bsxfun(@times,double(spikewords),keep(:)),1);
end

%% regularization
Kmax_new = min(Kmax + Kmax_add, Nneu);
Pki_m = regularization_pseudocount_Pki( Nki_emp_m, Nt, Pki_indep_m, lambda_reg, Kmax_new );

%% computation of statistics
Pk_l = bsxfun(@times, sum(Pki_m,2),1./((1:Kmax_new)'));
P0 = 1-sum(Pk_l);
Pi_l = sum(Pki_m,1);
mKi_l = (1:Kmax_new)*Pki_m;
end

