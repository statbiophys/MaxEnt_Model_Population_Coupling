function [  beta_l, alpha_l, gamma_l ] = resize_params_Pk_Pi_meanKi(  beta_l0 , alpha_l0, gamma_l0 )
%REFORMAT_PARAMS__MAXENT_PK_PI_MEANKI
% sets the parameters in " standard " format :
% we suppose ak_l(0) = 0, so Z = 1/P0
%        sum( hi_l = 0)
%        sum( gi_l = 0)

beta_l = beta_l0;
alpha_l = alpha_l0;
gamma_l = gamma_l0;

mh = mean(alpha_l);
alpha_l = alpha_l - mh;

mg = mean(gamma_l);
gamma_l = gamma_l - mg;

k_l = zeros(size(beta_l)); % make k_l in the same shape than a_l
k_l(1:end) = 1:length(k_l); 

beta_l = beta_l + mh*k_l + mg*(k_l.^2); 

end

