function [ h_m ] = infer_PsigmaK__autoscale( P0, Pki_m, a_start, error_max, params)
% for each k independently, (h_ik)_i is optimized by likelihood gradient assent,
% using the Hessian.
% x_(n+1) = x_n + a*(Hessian^-1)*(stats_target - stats_current)
%
% parameters a is ajusted automatically if the algorithm fails

%% Parameters
% a_start  (typ. 0.5 for N=160 neurons)
% params.use_parfor : if 1, function will use parfor (increases algorithm). if 0 it will not
%                     if parallel computing toolbox not installed, set to 0
%           (typ. 1)
%
% params.Nstep_max : if this number of step is reached, a is adjusted
%           (typ. 50000 for N=160 neurons)
%
% params.update_params_and_Hessian_rate : every update_params_and_Hessian_rate steps the Hessian is updated
%           (typ. 70 for N=160 neurons)

%%
fprintf('Beginning inference of MaxEnt for P(sigma_i,K) \n');
warning('off','MATLAB:nearlySingularMatrix');

if any(Pki_m == 0)
    error('Pki_m == 0 for certain values : must be avoided');
end

[Kmax, Nneu] = size(Pki_m);
h_m = zeros(Kmax, Nneu);
zk_l = ones(Kmax, 1);

if params.use_parfor; parforArg = Inf; else parforArg = 0; end
update_Hessian_rate = params.update_params_and_Hessian_rate; % avoid broadast of params variable in parfor loop
Nstep_max = params.Nstep_max;

parfor (k = 1:min(Kmax, Nneu-1), parforArg) % (if parforArg=0, equivalent to for)  (Neu-1) because for Neu, all neuron fire
    targ = k*Pki_m(k,:)./norm(Pki_m(k,:),1);
    
    h_l_start = zeros(1,Nneu);
    h_l_opt = NaN;
    a = a_start;
    
    while any(isnan(h_l_opt)) || any(isinf(h_l_opt))
        [h_l_opt,  Nstep] = solve_fixed_k(Nneu, h_l_start, targ,k, a, error_max, update_Hessian_rate, Nstep_max);
        if Nstep<Nstep_max %failed because a too large
            a = a/2;
        else
            a =a*1.3; % failed because a too small (took too long)
            h_l_opt = NaN;
        end
    end
    
    h_m(k,:)= h_l_opt;
    zk_l(k) = Zk_MaxEnt_PsigmaK( exp(h_l_opt), k );
    
    fprintf(['Done : k = ' int2str(k) ' with ' int2str(Nstep-1) ' steps \n']);
end

if Kmax == Nneu
    zk_l(Nneu) = Zk_MaxEnt_PsigmaK( exp(h_m(Nneu,:)), Nneu );
    
end

%% optimize the Vk_l, and add Vk_l/k to each h_ki to compensate
Pk_l = sum(Pki_m,2)./(1:Kmax)';

Z = 1/P0; % set to 0 the energy of silence ( sigma_i = 0 for all i )
Vk_l = log(Z*Pk_l./zk_l);
h_m = bsxfun(@plus, h_m, Vk_l(:)./(1:Kmax)');

warning('on','MATLAB:nearlySingularMatrix');
fprintf('End of computation \n');
end

function [h_l, step_i] = solve_fixed_k(Nneu, h_l_init, mean_si_l_targ, k, a, error_max, update_Hessian_rate, Nstep_max)

if k == 1
    h_l = log(mean_si_l_targ);
    h_l = h_l - mean(h_l);
    step_i = 1;
    
else
    h_l = h_l_init;
    [~, mean_si_l_curr, cov_sij] = get_stats_condk(Nneu, h_l, k);
    
    step_i = 1;
    error_val = max(abs(mean_si_l_targ - mean_si_l_curr));
    warning('off','MATLAB:nearlySingularMatrix');
    while error_max<error_val
        
        % here, h_l(1) is kept fixed, in order to remove the degeneration of the solution
        grad = (mean_si_l_targ(2:end) - mean_si_l_curr(2:end)); % gradient for loglikelihood
        
        h_l(2:end) = h_l(2:end) + (a*(cov_sij(2:end,2:end)\(grad')))';
        h_l = h_l - mean(h_l);
        
        if ~mod(step_i - floor(update_Hessian_rate), update_Hessian_rate)
            [~, mean_si_l_curr, cov_sij] = get_stats_condk(Nneu, h_l, k); % update curvature (covariance) matrix
        else
            [~, mean_si_l_curr] = get_stats_condk(Nneu, h_l, k);
        end
        
        if any(isnan(mean_si_l_curr))
            error_val = NaN;
        else
            error_val = max(abs(mean_si_l_targ - mean_si_l_curr));
        end
        
        if ~mod(step_i,20000) % not too many information, because parfoor loops sends it to promp asynchronously (thus not understandable)
            fprintf(['Processing                          k = ' int2str(k) ' reached ' int2str(step_i) ...
                ' steps      a : ' num2str(a) '     error ' num2str(error_val) '\n']);
        end
        
        step_i = step_i+1;
        if step_i == Nstep_max
            break;
        end
    end
    warning('on','MATLAB:nearlySingularMatrix');
    
    if isnan(error_val) || isinf(error_val)
        h_l = NaN;
    end
end
end


function [Zk, mean_si_l, cov_sij] = get_stats_condk(Nneu, h_l, k)
exph_l = exp(h_l);

Zk = Zk_MaxEnt_PsigmaK( exph_l, k );
mean_si_l = Zk_MaxEnt_PsigmaK_partialorder1( exph_l, k);
mean_si_l = mean_si_l.*exph_l/Zk;

if nargout>2 % cov sigma_i sigma_j
    if any(diff(h_l))
        mean_sisj = Zk_MaxEnt_PsigmaK_partialorder2(exph_l,k);
        mean_sisj = mean_sisj.*(exph_l(:)*exph_l(:)');
        mean_sisj = mean_sisj./Zk;
        mean_sisj = mean_sisj + diag(mean_si_l(:)); % for i==j
        
    else
        x = Zk_MaxEnt_PsigmaK(exph_l(3:end), k-2)*(exph_l(1)^2)/Zk; % for i >< j
        mean_sisj = x*ones(Nneu, Nneu);
        mean_sisj = mean_sisj + diag(mean_si_l(:) - x);  % for i==j  (the others are 0)
        
    end
    
    cov_sij = mean_sisj - (mean_si_l(:)*mean_si_l(:)'); % covariance
end
end
