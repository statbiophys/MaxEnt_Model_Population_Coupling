function [ h_m, likel_l ] = infer_Psigma_PK__autoscale( Pk_l, Pi_l, a_start, error_max, params )
%INFER_ALL_MAXENT_PSIGMA_PK
a = a_start;
h_m =Inf;
while any(isnan(h_m(:))) || any(isinf(h_m(:)))
    if nargout> 1
        [ h_m, likel_l ] = infer_Psigma_PK(Pk_l, Pi_l, a, error_max, params);
    else
        [ h_m ] = infer_Psigma_PK(Pk_l, Pi_l, a, error_max, params);
    end
    a = a/1.3;
    
    if any(isnan(h_m(:))) || any(isinf(h_m(:)))
        fprintf(['---> Trying new set of parameters : lower a set to a = '  num2str(a) '\n']);
    end
end

end

