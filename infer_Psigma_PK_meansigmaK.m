function [ h_m, likel_l] = infer_Psigma_PK_meansigmaK(Pk_l, Pi_l, mean_Ki, a, error_max, params)
%INFER_MAXENT_PSIGMA_PK_MEANSIGMAK
% fit MaxEnt reproducing P(sigma_i), P(K) and <K sigma_i> with K the population rate

% we set ak_l(0) = 0, so Z = 1/P0
%        sum( hi_l = 0)
%        sum( gi_l = 0)

% here K>length(Pk_l) are not considered (makes it faster to compute) (not used in practice)

%% additional parameters

% params.update_type : type of algorithm used. Must be one of the following
%      'gradient' : simple gradient ascent     x_(n+1) =  x_n + a*(stats_target - stats_current)
%
%      'momentum' : gradient ascent with momentum    x_(n+1) =  x_n + a*(stats_target - stats_current) + p*(x_n - x_(n-1))
%                         (see "On the momentum term in gradient descent learning algorithms", Qian 99)
%
%
%      'Hess' : gradient ascent with Hessian   x_(n+1) =  x_n + a*(Hessian^-1)*(stats_target - stats_current)
%
%      'Hess_update_a' :   gradient ascent with Hessian   x_(n+1) =  x_n + a*(Hessian^-1)*(stats_target - stats_current)
%                  with  'a' scaled automatically (see parameter params.coef_a)
%
%       (typ. 'Hess')
%
% params.init_mode : mode of parameters initialization
%     'all0'  : all parameters at 0
%     'nonzeroBeta'  :  alpha and gamma at 0, and beta in order to fit Pk_l in this case
%
% params.Nstep_max  % maximum number of steps before optimization ends (typ. 2000000)
% params.disp_info_rate  % display ongoing results every disp_info_rate step (typ. 50)
%
% parameters not common to all methods:
% params.update_params_and_Hessian_rate : parameters (e.g. Hessian) updated every update_params_and_Hessian_rate steps
%       (typ. 200 for update_type = 'Hess' and 160 neurons)
% params.p :  only if update_type='momentum', momentum of the gradient (typ 0.1)
% params.coef_a : only if update_type='Hess_update_a', a is updated to coef_a*2/( min(V)+ max(V))
%       with V the eigenvalues of the hessian    (typically <1, eg 0.8 for update_params_and_Hessian_rate=200 )
% params.lambda : only if 'Hess' : regularization factor of the hessian. 0 : no regularization (typ. 0, or 0.01)
%
% a : typ 0.15 for N=160 neurons, update_type='Hess' update_params_and_Hessian_rate = 100
%     typ 0.8 for N=20 neurons, update_type='Hess', lambda = 0.001, update_params_and_Hessian_rate = 2000

%%
if any(Pk_l == 0)
    error('Pk_l = 0 for certain values : must be avoided');
end
if any(Pi_l == 0)
    error('Pi_l == 0 for certain values : must be avoided');
end
if any(mean_Ki == 0)
    error('mean_Ki == 0 for certain values : must be avoided');
end
if nargout>1
    fprintf('WARNING: list of likelihoods at each step requested as output ! Makes algorithm was slower.  Only ouput h_m to increase speed.');
end

Nneu = length(Pi_l);
Kmax = length(Pk_l);

%% orthogonal projection from the space of parameters without redundancy

l1 = zeros(1,Kmax + 2*Nneu);
l2 = zeros(1,Kmax + 2*Nneu);

l1(Kmax+1) = 1; % parameter i = 0 are fixed
l2(Kmax+Nneu+1) = 1;

P = null([l1; l2]);
PP = P*P';

%%
vals_targ = [Pk_l(:); Pi_l(:); mean_Ki(:)]; %target statistics

params_curr = zeros(Kmax + 2*Nneu,1);
switch params.init_mode
    case 'all0'
        fprintf('Initialization with all parameters 0 \n');
    case 'nonzeroBeta'
        warning('off','MATLAB:nchoosek:LargeCoefficient'); % we don't need a precise estimate here
        for k = 1:Kmax % first estimation of beta_k(k) if all alpha and gamma are 0
            params_curr(k) = log(Pk_l(k)/nchoosek(Nneu, k));
        end
        warning('on','MATLAB:nchoosek:LargeCoefficient');
        fprintf('Initialization with non-zero beta_k, alpha_i=0 and gamma_i = 0 \n');
end

switch params.update_type
    case {'gradient','momentum'}
        [Pk_l_curr, Pi_l_curr, mean_Ki_curr] = stats_Psigma_PK_meansigmaK__alphabetagamma(params_curr(1:Kmax), params_curr((Kmax+1):(Kmax+Nneu)), params_curr((Kmax+Nneu+1):end));
    case {'Hess_update_a', 'Hess'}
        [Pk_l_curr, Pi_l_curr, mean_Ki_curr, Hess_m] = stats_Psigma_PK_meansigmaK__alphabetagamma(params_curr(1:Kmax), params_curr((Kmax+1):(Kmax+Nneu)), params_curr((Kmax+Nneu+1):end));
        switch params.update_type
            case 'Hess'
                if params.lambda == 0; %1e-5; % for regularization %0 was fine for 160 neurons
                    U = P'*Hess_m*P;
                else
                    [~,V] = eig(Hess_m,'vector');
                    U = P'*(Hess_m+params.lambda*max(V)*eye(length(Hess_m)))*P;
                end
                
            case 'Hess_update_a'
                U = P'*Hess_m*P;
                [~,V] = eig(U,'vector');
                V = sort(V,'ascend');
                a = params.coef_a*2/(V(1)+V(end));
        end
end
vals_curr = [Pk_l_curr(:); Pi_l_curr(:); mean_Ki_curr(:)];
err = norm(vals_curr - vals_targ, Inf);

fprintf('Beginning inference of MaxEnt for P(sigma_i), P(K) and < sigma_i K > \n');
if nargout>1
    likel_l = [];
end

step_i = 1; % iteration number
switch params.update_type
    case 'momentum'
        last_Delta =  a*(P*P')*(vals_targ - vals_curr);
end
warning('off','MATLAB:nearlySingularMatrix');
while (err>error_max) && (step_i<params.Nstep_max)
    switch params.update_type
        case {'gradient','Hess_update_a'}
            d = PP*(vals_targ - vals_curr);
            last_Delta = a*d;
            
        case 'momentum'% add momentum term
            d = PP*(vals_targ - vals_curr);
            last_Delta = a*d +params.p*last_Delta;
            
        case 'Hess'            
            d = U\(P'*(vals_targ - vals_curr));
            last_Delta = a*(P*d);
            
    end
    params_curr = params_curr + last_Delta;
    
    [ beta_kl , alpha_il, gamma_il ] = resize_params_Pk_Pi_meanKi__computation(params_curr(1:Kmax), params_curr((Kmax+1):(Kmax+Nneu)), params_curr((Kmax+Nneu+1):end) );
    if ~mod(step_i, params.update_params_and_Hessian_rate) % avoids going out of the range of double precision when using exponential
        params_curr = [beta_kl(:); alpha_il(:); gamma_il(:)];
    end
    
    switch params.update_type
        case {'gradient','momentum'}
            if nargout>1 % need Z for the likelihood
                [Pk_l_curr, Pi_l_curr, mean_Ki_curr, ~, Z ] = ...
                    stats_Psigma_PK_meansigmaK__alphabetagamma(beta_kl, alpha_il, gamma_il);
            else
                [Pk_l_curr, Pi_l_curr, mean_Ki_curr ] = ...
                    stats_Psigma_PK_meansigmaK__alphabetagamma(beta_kl, alpha_il, gamma_il);
            end
            
        case {'Hess','Hess_update_a'}
            if mod(step_i-1, params.update_params_and_Hessian_rate) 
                if nargout>1 % need Z for the likelihood
                    [Pk_l_curr, Pi_l_curr, mean_Ki_curr, ~, Z ] = ...
                        stats_Psigma_PK_meansigmaK__alphabetagamma(beta_kl, alpha_il, gamma_il);
                else
                    [Pk_l_curr, Pi_l_curr, mean_Ki_curr ] = ...
                        stats_Psigma_PK_meansigmaK__alphabetagamma(beta_kl, alpha_il, gamma_il);
                end
                
            else
                [Pk_l_curr, Pi_l_curr, mean_Ki_curr, Hess_m, Z ] = ...
                    stats_Psigma_PK_meansigmaK__alphabetagamma(params_curr(1:Kmax), params_curr((Kmax+1):(Kmax+Nneu)), params_curr((Kmax+Nneu+1):end));
                
                switch params.update_type
                    case 'Hess'
                        if params.lambda == 0; 
                            U = P'*Hess_m*P;
                        else
                            [~,V] = eig(Hess_m,'vector');
                            U = P'*(Hess_m+params.lambda*max(V)*eye(length(Hess_m)))*P;
                        end
                        
                    case 'Hess_update_a'
                        U = P'*Hess_m*P;
                        [~,V] = eig(U,'vector');
                        V = sort(V,'ascend');
                        a = params.coef_a*2/(V(1)+V(end));
                end
            end
    end
    
    if nargout>1 % list of the log likelihood tested
        likel_l(end+1) = - log(Z) + sum(Pk_l(:).*beta_kl(:)) + sum(Pi_l(:).*alpha_il(:)) + sum(mean_Ki(:).*gamma_il(:));
    end
    
    vals_curr = [Pk_l_curr(:); Pi_l_curr(:); mean_Ki_curr(:)];
    err = norm(vals_curr - vals_targ, Inf);
    
    if ~mod(step_i, params.disp_info_rate)
        fprintf([' Step ' int2str(step_i) '    error : ' num2str(err) '   >< ' num2str(error_max) '    - param norm : ' num2str(norm(params_curr,1))]);
        fprintf('\n');
    end
    step_i = step_i +1;
end
warning('on','MATLAB:nearlySingularMatrix');

if isnan(err)
    fprintf('Algorithm failed : NaN error found\n');
    fprintf('All parameters set to Inf or NaN to avoid further usage \n');
    fprintf('Try with more Hessian updates or smaller step coefficient (a)\n');
    
    params_curr = Inf*ones(size(params_curr));
end

if step_i == params.Nstep_max
    fprintf(['Max number of iterations: ' int2str(params.Nstep_max) '   reached : algorithm stopped.\n']);
    fprintf('All parameters set to NaN to avoid further usage \n');
    
    params_curr = Inf*ones(size(params_curr));
end

beta_kl = params_curr(1:Kmax);
alpha_il = params_curr((Kmax+1):(Kmax+Nneu))';
gamma_il = params_curr((Kmax+ Nneu+1):end)';
[ beta_kl , alpha_il, gamma_il ] = resize_params_Pk_Pi_meanKi( beta_kl , alpha_il, gamma_il );
h_m = get_h_PkPimKi( beta_kl , alpha_il, gamma_il );
end


