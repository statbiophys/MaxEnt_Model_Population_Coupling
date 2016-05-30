function Pki_reg_m = regularization_pseudocount_Pki( Nki_emp_m, Nt, Pki_indep_m_full, lambda_reg, Kmax )
% we set to 0 all probabilities P(K>Kmax). By default Kmax=Nneu
Nneu = size(Nki_emp_m,2);
if size(Pki_indep_m_full,2) ~= Nneu
    error('Should have the same number of neurons');
end

if (size(Nki_emp_m,1)>Kmax) && any(any(Nki_emp_m((Kmax+1):end,:)))
    error('Nki_emp_m is larger than Kmax : data would be lost');
end

if nargin == 4
    Kmax = Nneu; % by default
    Pki_indep = Pki_indep_m_full;
    
elseif nargin == 5
    Pki_indep = Pki_indep_m_full(1:Kmax,1:Nneu);
    for neu_i = 1:Nneu % cut the probability for k>Kmax, but keep same probability for silence
        Pki_indep(:,neu_i) = Pki_indep(:,neu_i)*norm(Pki_indep_m_full(:,neu_i),1)/norm(Pki_indep(:,neu_i),1);
    end
    
end

Nki_emp_new = Nki_emp_m;
Nki_emp_new((end+1):Kmax,:) = 0;

%% regularize P(K)
Pk_l_indep = sum(Pki_indep,2)./((1:Kmax)');
Nk_emp_l = sum(Nki_emp_new,2)./((1:Kmax)');
Pk_sm_l = (Nk_emp_l+lambda_reg*Pk_l_indep)/(Nt+lambda_reg);

%% regularize P(sigma_i | K)
Pi_condK_indep = bsxfun(@times, Pki_indep, 1./Pk_l_indep);
Pi_condK_sm = bsxfun(@times,(Nki_emp_new(1:Kmax,:) + lambda_reg*Pi_condK_indep),1./(Nk_emp_l + lambda_reg));

%% joint probability P(sigma, K)
Pki_reg_m = bsxfun(@times, Pk_sm_l, Pi_condK_sm);

end

