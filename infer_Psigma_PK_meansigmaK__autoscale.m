function [ h_m, likel_l] = infer_Psigma_PK_meansigmaK__autoscale( Pk_l, Pi_l, mKi_l, a_start, error_max, params)

a = a_start;
h_m =Inf;
while any(isnan(h_m(:))) || any(isinf(h_m(:)))
    if nargout> 1
        [ h_m, likel_l ] = infer_Psigma_PK_meansigmaK(Pk_l, Pi_l, mKi_l, a, error_max, params);
    else
        [ h_m ] =  infer_Psigma_PK_meansigmaK(Pk_l, Pi_l, mKi_l, a, error_max, params);
    end
    a = a/1.3;
    
    if any(isnan(h_m(:))) || any(isinf(h_m(:)))
        fprintf(['---> Trying new set of parameters : lower a set to a = '  num2str(a) '\n']);
    end
end

end

