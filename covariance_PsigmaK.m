function [Cpred_m, r_pred] = covariance_PsigmaK( h_m, Z )
%GET_COV_MAXENT_PSIGK covariance and firing rate predicted from model
Nneu = size(h_m,2);
Cpred_cond = zeros(Nneu, Nneu, Nneu);

exph_ki_m = exp(h_m);

for k = 2:size(h_m,1) % 2 because for k = 0,1 the neurons cannot fire at the same time
    for ii = 1:Nneu
        for jj = (ii+1):Nneu
            Cpred_cond(k,ii,jj) = exp(h_m(k,ii)+h_m(k,jj))*Zk_MaxEnt_PsigmaK(exph_ki_m(k,[1:(ii-1), (ii+1:jj-1), (jj+1:Nneu)]), k-2);
            
        end
    end
    fprintf([' k : ' int2str(k) ' done \n']);
end

r_pred = firingrate_PsigmaK( h_m, Z );

Cpred_m = squeeze(sum(Cpred_cond,1))/Z;
Cpred_m = Cpred_m + Cpred_m';

for ii = 1:Nneu % because it was only computed for ii <> jj
    Cpred_m(ii,ii) = r_pred(ii);
end

Cpred_m = Cpred_m - r_pred(:)*r_pred(:)';

end

