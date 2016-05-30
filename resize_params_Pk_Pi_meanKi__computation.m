function [  beta_l, alpha_l, gamma_l ] = resize_params_Pk_Pi_meanKi__computation( beta_l0 , alpha_l0, gamma_l0 )
% sets the parameters in " computable" format (convenient for computation) :

beta_l = beta_l0;
alpha_l = alpha_l0;
gamma_l = gamma_l0;

ma = mean(beta_l);
mh = mean(alpha_l);
mg = mean(gamma_l);

Kmax=length(beta_l0);

Ru = (Kmax+1)/2; % K(K+1)/2   /K  (mean of the squares)
Rv = Ru*(2*Kmax+1)/3;

Ca = 1/(Ru - Kmax - Rv/Kmax);
Ch = Ca*(Rv/Kmax + Kmax);
Cg = Ca*Rv;

u = Ca*ma + Ch*mh + Cg*mg;

v = -((mh+u)/Kmax + mg);

k_l = zeros(size(beta_l)); % make k_l in the same shape than a_l
k_l(1:end) = 1:length(k_l); 
beta_l = beta_l -u*k_l -v*(k_l.^2); 

alpha_l = alpha_l+u;
gamma_l = gamma_l+v;

end

