function h_m = get_h_PkPimKi( beta_l, alpha_l, gamma_l )

Kmax =length(beta_l);
Nneu = length(alpha_l);

h_m = repmat(beta_l(:)./((1:Kmax)'),[1 Nneu]);
h_m = h_m + repmat(alpha_l(:)',[Kmax, 1]);
h_m = h_m + (1:Kmax)'*gamma_l(:)';

end

