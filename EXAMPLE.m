% .mex files must be compiled. Use COMPILE_mex_files.m for this.

% The algorithms target the maximum of models likelihood, but the error if the maximum error
% of models predictions. Thus while the likelihood is suppose to increase during the optimization
% the error we measure is not always necessarily decreasing (and can increase, for example at the beginning
% of the algorithm). The error is communicated in the prompt during the computation. 

% One can also get the likelihood at each step of the experiments for the minimal and linear-coupling models
% by taking output [h, likelihood_list] = ... instead of [ h ] ...     eg for the minimal model:
%     [ h_minim, likelihood_list ] = infer_Psigma_PK__autoscale(Pk_l, Pi_l, a_start, error_max,params_minim);
% instead of 
%     [ h_minim ] = infer_Psigma_PK__autoscale(Pk_l, Pi_l, a_start, error_max,params_minim);
% But this slows down the computation a lot, so it shouldn't be use unless for parameters adjustment

% TEST: in order to test the code, you can generate random data using:
% spikewords = rand(20000,50)>(1-0.2);   % 20000 responses from 50 neurons, with firing probability of 0.2

[Nt, Nneu] = size(spikewords); % spikewords is a matrix of size (Ntimes, Nneurons) with 1 if a neuron spikes, 0 else

%% regularization

lambda_reg = 1;
Kmax_add = 3;
[ Pk_l, Pi_l, mKi_l, Pki_m, P0, Kmax ] = regularize( spikewords,  lambda_reg, Kmax_add);

%% learn minimal model      [ see parameters explanations in function code ]

error_max = 10^(-6);
a_start = 1; % this parameter is adjusted automatically: if the algorithm fails, it starts with a smaller a

params_minim.update_type = 'Hess';
params_minim.disp_info_rate = 20; 

params_minim.Nstep_max = 2000000; 

params_minim.update_params_and_Hessian_rate = 200;
params_minim.p = 0.1; 
params_minim.coef_a = 0.8; 
params_minim.lambda = 0;

tic % inference
[ h_minim ] = infer_Psigma_PK__autoscale(Pk_l, Pi_l, a_start, error_max,params_minim);
toc

% prediction and Z
[ P0_pred_minim, Pk_l_pred_minim, Pki_m_pred_minim, Z_minim ] = prediction_PsigmaK( h_minim );

%% learn linear-coupling model      [ see parameters explanations in function code ]

error_max = 10^(-6);
a_start = 0.15; % this parameter is adjusted automatically : if the algorithm fails, it starts with a smaller a

params_lin.Nstep_max = 2000000; 

params_lin.update_type = 'Hess';
params_lin.disp_info_rate = 50; 
params_lin.init_mode = 'nonzeroBeta';

params_lin.update_params_and_Hessian_rate = 100;

params_lin.p = 0.8; 
params_lin.coef_a = 0.8; 
params_lin.lambda = 0;

tic % inference
[h_lin] = infer_Psigma_PK_meansigmaK__autoscale(Pk_l, Pi_l, mKi_l,  a_start, error_max, params_lin);
toc

% prediction and Z
[ P0_pred_lin, Pk_l_pred_lin, Pki_m_pred_lin, Z_lin ] = prediction_PsigmaK( h_lin );

%% learn complete coupling model      [ see parameters explanations in function code ]

error_max = 10^(-6);
a_start = 0.5;

params_comp.use_parfor = 1;
params_comp.update_params_and_Hessian_rate = 70;
params_comp.Nstep_max = 50000;

tic % inference
[ h_comp ] =  infer_PsigmaK__autoscale( P0, Pki_m, a_start, error_max, params_comp );
toc

% prediction and Z
[ P0_pred_comp, Pk_l_pred_comp, Pki_m_pred_comp, Z_comp ] = prediction_PsigmaK( h_comp );

%% pairwise covariance prediction 

[Cov_prediction_minim] = covariance_PsigmaK(h_minim, Z_minim);
[Cov_prediction_lin] = covariance_PsigmaK(h_lin, Z_lin);
[Cov_prediction_comp] = covariance_PsigmaK(h_comp, Z_comp);
