function [ Pki_m, Pk_l, P0, r_l, exph_l ] = probas_indep_Ising( spikewords )
%INDEPENDENT_MODEL_PROBAS
% probabilities for a model of independant neurons reproducing only firing rates
% the model has form   P(sigma) = 1/Z exp( sum_i h_i sigma_i  )
%
% INPUT
% spikewords : (Ntime, Nneuron) matrix

% OUTPUT
% Pki_m : [(Nneuron, Nneuron) matrix]  Pki_m(i,k) = P(sigma_i = 1, K=k)   for independent model
% Pk_l : [(Nneuron) list]  Pk_l(k) = P(K=k)   for independent model
% P0 : [number] = P(K=0)   for independent model
% r_l : [number] = P(sigma_i = 1);   for independent model
% exph_l : [(Nneuron) list] exph_l(i) = exp(h_i)

Nneu = size(spikewords,2);

r_l = mean(spikewords,1); % empirical firing proba
exph_l = r_l./(1-r_l);

Z_indep = Z_indep_Ising(exph_l);

P0 = 1/Z_indep;
Pki_m = zeros(Nneu, Nneu);

for k = 1:Nneu
    for n = 1:Nneu
        Pki_m(k,n) = exph_l(n)*Zk_MaxEnt_PsigmaK(exph_l([1:(n-1), (n+1):Nneu]),k-1);
    end
end
Pki_m = Pki_m/Z_indep; % Pki_indep_m(k,i) is the proba of sigma_i = 1 and K=k in the independent model with empirical firing rates
Pk_l = sum(Pki_m,2)./(1:Nneu)';

end

