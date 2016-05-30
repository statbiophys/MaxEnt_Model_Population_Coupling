function [  beta_l, alpha_l] = resize_params_Pk_Pi( beta_l0 , alpha_l0)
% sets the parameters in " standard " format :
% we suppose betak_l(0) = 0, so Z = 1/P0
%        sum( alphai_l = 0)

beta_l = beta_l0;
alpha_l = alpha_l0;

mh = mean(alpha_l);
alpha_l = alpha_l - mh;

k_l = zeros(size(beta_l)); % make k_l in the same shape than beta_l
k_l(1:end) = 1:length(k_l); 

beta_l = beta_l + mh*k_l; 

end

