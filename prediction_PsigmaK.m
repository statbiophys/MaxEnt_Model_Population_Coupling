function [ P0, Pk_l, Pki_m, Z ] = prediction_PsigmaK( h_m)
%DIRECT_PROBLEM_MAXENT_PK_PSIGMAK
% INPUT
% h_ki_m : (Nneu, Neu) mat of the h_ki. 1 raw : 1 k; 1 column : 1 neuron

% OUTPUT
% Pki_m(k,i) : proba of (sigma_i =1, K=k)
%
% Pk_l(k) : proba of (K = k)
%
% P0 : proba of (K = 0)

[Kmax, Nneu] = size(h_m);
Zk_l = zeros(Kmax,1);

exph_ki_m = exp(h_m);

for k = 1:Kmax
    Zk_l(k) = Zk_MaxEnt_PsigmaK( exph_ki_m(k,:), k );
end

Z = 1 + sum(Zk_l); % partition function.  silence has energy 1

P0 = 1/Z;
Pk_l = Zk_l/Z;

Pki_m = ones(Kmax,Nneu);
for i_c = 1:Nneu
    for k = 2:Kmax
        Pki_m(k, i_c) = Zk_MaxEnt_PsigmaK( exph_ki_m(k,[1:(i_c-1), (i_c+1):end]), k-1);
    end
end
Pki_m = bsxfun(@times, Pki_m, exph_ki_m);
Pki_m = Pki_m/Z;

end

