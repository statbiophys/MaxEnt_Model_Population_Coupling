function [r_pred, r_iK_pred ] = firingrate_PsigmaK( h_m, Z )
%GET_R_CONDK_MAXENT_PSIGK
Nneu = size(h_m,2);
r_iK_pred = zeros(Nneu, Nneu);

exph_ki_m = exp(h_m);
for k = 1:size(h_m,1)
%     fprintf([' beginning k : ' int2str(k) ' \n']);
    for ii = 1:Nneu
        r_iK_pred(k,ii) = exp(h_m(k,ii))*Zk_MaxEnt_PsigmaK(exph_ki_m(k,[1:(ii-1), (ii+1:Nneu)]), k-1);
        
    end
    
    %end
%     fprintf([' k : ' int2str(k) ' done \n']);
end

r_pred = sum(r_iK_pred,1)/Z;

end

