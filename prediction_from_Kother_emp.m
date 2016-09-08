function [ Pi_condKo_m, errstd_Pi_condK, Nko_m ] = prediction_from_Kother_emp( spikewords, pseudocnt_lambda, Kmax_add )
%PREDICTION_FROM_KOTHER_EMP
% Pi_condKo_l(Koi+1, neu_i) is the proba that neuron neu_i spikes if
% the other neurons have Koi spikes

[ Pk_l,~,~, Pki_m, P0 ] = regularize( spikewords, pseudocnt_lambda, Kmax_add);

[Kmax,Nneu] = size(Pki_m);
Pk_i_plus = zeros(Kmax, Nneu);
Pk_i_minus = zeros(Kmax, Nneu);

for k=0:(Kmax-1)
    Pk_i_plus(k+1,:) = Pki_m(k+1,:);
end

Pk_i_minus(1,:) = P0;
for k=1:(Kmax-1)
    Pk_i_minus(k+1,:) = Pk_l(k)-Pki_m(k,:);
end

Pi_condKo_m = Pk_i_plus./(Pk_i_plus+Pk_i_minus);


k_l = sum(spikewords, 2);
Kmax = max(k_l);
%
% for k=1:Kmax
%     Pi_condKo_m(k,:) = mean(spikewords(k_l == k,:),1)*mean(k_l == k);
%     Pi_condKo_m(k,:) = Pi_condKo_m(k,:)./mean( bsxfun(@plus,k_l, - double(spikewords)) ==(k-1));
% end


% regularized version
if nargin>1
    Nko_m = zeros(Kmax, Nneu);
end

for neu_i = 1:Nneu
    if nargin>1
        n_l = histc(k_l,0:(Kmax-1));
        Nko_m(:,neu_i) = n_l(:);
    end
end

if nargin>1
    errstd_Pi_condK = Pi_condKo_m - Pi_condKo_m.^2;
    errstd_Pi_condK = sqrt( errstd_Pi_condK(1:Kmax,:)./Nko_m );
end
end

