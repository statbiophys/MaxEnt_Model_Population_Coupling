function Pi_condKo_m = prediction_from_Kothers( exph_ki_m )
%PREDICTION__FROM_KOTHERS
% Pi_condKo_l(Koi+1, neu_i) is the proba that neuron neu_i spikes if 
% the other neurons have Koi spikes 

if any(exph_ki_m(:)<0)
    error('input must be exph_ki_m, but negative number found');
end

[Kmax,Nneu] = size(exph_ki_m);
Zk_i_plus = zeros(Kmax, Nneu);
Zk_i_minus = zeros(Kmax, Nneu);

for k=0:(Kmax-1)
    Zk_i_plus(k+1,:) = Zk_MaxEnt_PsigmaK_partialorder1(exph_ki_m(k+1,:), k+1);
end
Zk_i_plus = Zk_i_plus.*exph_ki_m;

Zk_i_minus(1,:) = 1;
for k=1:(Kmax-1)
    Zk_i_minus(k+1,:) = Zk_MaxEnt_PsigmaK_partialorder1(exph_ki_m(k,:), k+1);
end

Pi_condKo_m = Zk_i_plus./(Zk_i_plus+Zk_i_minus);
end

