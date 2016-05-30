function [h_m, likel_l] = infer_Psigma_PK(Pk_l, Pi_l, a, error_max, params)
%INFER_PSIGMA_PK_MEANSIGMAK
% fit MaxEnt reproducing P(sigma_i) and P(K)  with K the population rate

% we set betak_l(0) = 0, so Z = 1/P0
%        sum( alphai_l = 0)

% here K> length(Pk_l) = Kmax are not considered (faster to compute)

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
% params.Nstep_max  % maximum number of steps before optimization ends (typ. 2000000)
% params.disp_info_rate  % display ongoing results every disp_info_rate step (typ. 20)
%
% parameters not common to all methods:
% params.update_params_and_Hessian_rate : parameters (e.g. Hessian) updated every update_params_and_Hessian_rate steps
%       (typ. 200 for update_type = 'Hess' and 160 neurons)
% params.p :  only if update_type='momentum', momentum of the gradient (typ 0.1)
% params.coef_a : only if update_type='Hess_update_a', a is updated to coef_a*2/( min(V)+ max(V))
%       with V the eigenvalues of the hessian    (typically <1, eg 0.8 for update_params_and_Hessian_rate=200 )
% params.lambda : only if 'Hess' : regularization factor of the hessian. 0 : no regularization (typ. works well without regularization, or 0.01)

% OUPUT:
% likel_l : list of likelihood at each step
%%
if any(Pk_l == 0)
    error('Pk_l = 0 for certain values : must be avoided');
end
if any(Pi_l == 0)
    error('Pi_l == 0 for certain values : must be avoided');
end
if nargout>1
    fprintf('WARNING: list of likelihoods at each step requested as output ! Makes algorithm was slower.  Only ouput h_m to increase speed.');
end

Nneu = length(Pi_l);

%% orthogonal projection to the space of solution without redundancy
l1 = ones(1,Nneu); % must be horizontal
P = null(l1);
PP = P*P';

%%
Pk_l_targ = Pk_l(:);
Pi_l_targ = Pi_l(:);

alpha_il_curr = zeros(Nneu,1);

switch params.update_type % initialization
    case {'gradient','momentum'}
        [ Pi_l_curr, beta_kl_curr] = stats_Psigma_PK__Pkalpha(Pk_l_targ, alpha_il_curr);
    case {'Hess_update_a', 'Hess'}
        [ Pi_l_curr, beta_kl_curr, Hess_m] = stats_Psigma_PK__Pkalpha(Pk_l_targ, alpha_il_curr);
        switch params.update_type
            case 'Hess'
                %                 U = P'*eye(length(Hess_m))*P;
                if params.lambda == 0; % for regularization %0 was fine for 160 neurons
                    U = P'*Hess_m*P;
                else
                    %                     U = P'*(Hess_m+params.lambda*eye(length(Hess_m)))*P;
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

err = norm(Pi_l_curr(:) - Pi_l_targ(:), Inf); % Pk_l has exact fit

fprintf('Beginning inference of MaxEnt for P(sigma_i) and P(K) \n');
if nargout>1 % list of the log likelihood tested
    likel_l = [];
end

step_i = 1; % iteration number
switch params.update_type
    case 'momentum'
        last_Delta =  a*(P*P')*(Pi_l_targ - Pi_l_curr);
end
warning('off','MATLAB:nearlySingularMatrix');
while (err>error_max) && (step_i<params.Nstep_max) && ~any(isinf(Pi_l_curr)) && ~any(isinf(alpha_il_curr))
    switch params.update_type
        case {'gradient','Hess_update_a'}
            d = PP*(Pi_l_targ - Pi_l_curr);
            last_Delta = a*d;
            
        case 'momentum'
            d = PP*(Pi_l_targ - Pi_l_curr);
            last_Delta = a*d +params.p*last_Delta; % add momentum term
            
        case 'Hess'
            d = U\(P'*(Pi_l_targ - Pi_l_curr));
            last_Delta = a*(P*d);
            
    end
    alpha_il_curr = alpha_il_curr + last_Delta;
    
    switch params.update_type
        case {'gradient','momentum'}
            [ Pi_l_curr, beta_kl_curr] = stats_Psigma_PK__Pkalpha(Pk_l_targ, alpha_il_curr);
            
        case {'Hess','Hess_update_a'}
            if mod(step_i-1, params.update_params_and_Hessian_rate)
                [ Pi_l_curr, beta_kl_curr] = stats_Psigma_PK__Pkalpha(Pk_l_targ, alpha_il_curr);
                
            else
                [ Pi_l_curr, beta_kl_curr, Hess_m] = stats_Psigma_PK__Pkalpha(Pk_l_targ, alpha_il_curr);
                
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
    
    if nargout>1
        likel_l(end+1) = sum(Pk_l_targ(:).*beta_kl_curr(:)) +  sum(Pi_l_targ(:).*alpha_il_curr(:) );
    end
    
    err = norm(Pi_l_curr - Pi_l_targ, Inf);
    if any(isnan(err)) ||  any(isinf(beta_kl_curr)) || ~mod(step_i+1,460)
    end
    if ~mod(step_i, params.disp_info_rate)
        fprintf([' Step ' int2str(step_i) '    error : ' num2str(err) '   >< ' num2str(error_max) ...
            '    - param norm a ' num2str(norm(beta_kl_curr,1)) '  h ' num2str(norm(alpha_il_curr,1))]);
        fprintf('\n');
    end
    step_i = step_i +1;
end
warning('on','MATLAB:nearlySingularMatrix');

if isnan(err) % check errors
    fprintf('Algorithm failed : NaN error found\n');
    fprintf('All parameters set to Inf to avoid further usage \n');
    fprintf('Try with more Hessian updates (lower update_params_and_Hessian_rate) or smaller step coefficient (lower a)\n');
    
    beta_kl_curr = Inf*ones(size(beta_kl_curr));
    alpha_il_curr = Inf*ones(size(alpha_il_curr));
end

if step_i == params.Nstep_max % check stopped because of max number of iteration reached
    fprintf(['Max number of iterations: ' int2str(params.Nstep_max) '   reached : algorithm stopped.\n']);
    fprintf('All parameters set to NaN to avoid further usage \n');
    beta_kl_curr = NaN*ones(size(beta_kl_curr));
    alpha_il_curr = NaN*ones(size(alpha_il_curr));
end

[ beta_kl , alpha_il] = resize_params_Pk_Pi( beta_kl_curr , alpha_il_curr);
h_m = get_h_PkPi(beta_kl , alpha_il);

if nargout>1
    likel_l = likel_l -log(1-sum(Pk_l));
end
end


