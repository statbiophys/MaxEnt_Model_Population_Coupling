function h_m = get_h_PkPi(beta_l, alpha_l )

Nneu = length(alpha_l);
Kmax_PkPi = length(beta_l);

h_m = repmat(beta_l(:)./((1:Kmax_PkPi)'),[1 Nneu]);
h_m = h_m + repmat(alpha_l(:)',[Kmax_PkPi, 1]);

end

