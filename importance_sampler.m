function [] = importance_sampler(raw_data, analysis_settings)

% [] = IMPORTANCE_SAMPLER(RAW_DATA, ANALYSIS_SETTINGS)
% 
% Purpose
% 
% Recovers a curve that best explains the relationship between the predictor and dependent variables
% 
% Input
%
% --raw_data: The data matrix (total number of trials x 6 columns). Refer to RUN_IMPORTANCE_SAMPLER()
% --analysis_settings: A struct that holds algorithm relevant settings. Refer to RUN_IMPORTANCE_SAMPLER()
%
% Output
%
% --Saves a .mat file in current_path/analysis_id/analysis_id_importance_sampler.mat
%
% Example usage:
%
% importance_sampler(raw_data_matrix, analysis_settings_struct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the P-CIT toolbox released under the BSD license.
% Copyright (c) 2012, Princeton University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in %the documentation and/or other materials provided with the distribution.
% Neither the name of the Princeton University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT %NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = clock;
disp(sprintf('Start time %d/%d %d:%d', time(2), time(3), time(4), time(5)));

% Checks if the correct number of arguments are passed in
if nargin < 2, error('Missing input parameters'); end

% Resetting the random number seed based on the clock
rand('twister', sum(100 * clock));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing the data matrix and updating the analysis_settings struct with additional/missing information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[preprocessed_data, ana_opt] = preprocessing_setup(raw_data, analysis_settings);
% Removing away the old information
clear raw_data; clear analysis_settings;

%%%%%%%%%%%%%%%
% House keeping
%%%%%%%%%%%%%%%
importance_sampler = struct(); % Creating the output struct
hold_betas_per_iter = NaN(ana_opt.em_iterations+1, 2); % Matrix to hold the betas over em iterations
exp_max_fval = NaN(ana_opt.em_iterations, 1); % Matrix to hold the fvals over em iterations
normalized_w = NaN(ana_opt.em_iterations+1, ana_opt.particles);  % Vector to hold the normalized weights
global tau
global bounds
global w
global net_effects
global dependent_var

%%%%%%%%%%%%%%%%%%
% Fetch parameters
%%%%%%%%%%%%%%%%%%
tau = ana_opt.tau; % Store the tau for convenience
bounds = family_of_curves(ana_opt.curve_type, 'get_bounds'); % Get the curve parameter absolute bounds
nParam = family_of_curves(ana_opt.curve_type, 'get_nParams'); % Get the number of curve parameters
hold_betas = [ana_opt.beta_0, ana_opt.beta_1]; % Store the betas into a vector

for em = 1:ana_opt.em_iterations % For every em iteration

	hold_betas_per_iter(em, :) = hold_betas; % Store the logreg betas over em iterations
	disp(sprintf('Betas: %0.4f, %0.4f', hold_betas(1), hold_betas(2)));
	disp(sprintf('EM Iteration: %d', em));

	% Initialize the previous iteration curve parameters, weight vector, net_effects and dependent_var matrices
	prev_iter_curve_param = NaN(ana_opt.particles, family_of_curves(ana_opt.curve_type, 'get_nParams')); % Matrix to hold the previous iteration curve parameters
	w = NaN(1, ana_opt.particles); % Vector to hold normalized weights
	net_effects = NaN(length(ana_opt.net_effect_clusters), ana_opt.particles); % Matrix to hold the predictor variables (taking net effects if relevant) over all particles
	dependent_var = []; % This vector cannot be initialized in advance since we don't know the length of this vector (dropping outliers)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Sampling curve parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if em == 1 % only for the first em iteration
		param = common_to_all_curves(ana_opt.curve_type, 'initial_sampling', ana_opt.particles, ana_opt.resolution); % Good old uniform sampling
	else % for em iterations 2, 3, etc
		% Sample curve parameters from previous iteration's curve parameters based on normalized weights
		prev_iter_curve_param = param; % storing the previous iteration's curve parameters since we need them to compute likelihood
		% Here we sample curves (with repetitions) based on the weights
		param = prev_iter_curve_param(randsample([1:ana_opt.particles], ana_opt.particles, true, normalized_w(em-1, :)), :);
		
		% Add Gaussian noise since some curves are going to be identical due to the repetitions
		% NOISE: Sample from truncated normal distribution using individual curve parameter bounds, mean = sampled curve parameters and sigma = tau
		for npm = 1:nParam
			param(:, npm) = truncated_normal(bounds(npm, 1), bounds(npm, 2), param(:, npm), tau, ana_opt.particles);
		end
	end
	param = common_to_all_curves(ana_opt.curve_type, 'check_if_exceed_bounds', param); % Check whether curve parameters lie within the upper and lower bounds
	if strcmp(ana_opt.curve_type, 'horz_indpnt')
		param = common_to_all_curves(ana_opt.curve_type, 'sort_horizontal_params', param); % Check if the horizontal curve parameters are following the right trend i.e. x1 < x2
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Compute the likelihood over all subjects (i.e. log probability mass function if logistic regression)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% This is where we use the chunking trick II
	for ptl_idx = 1:size(ana_opt.ptl_chunk_idx, 1)
		output_struct = family_of_curves(ana_opt.curve_type, 'compute_likelihood', ana_opt.net_effect_clusters, ana_opt.ptl_chunk_idx(ptl_idx, 3),...
				param(ana_opt.ptl_chunk_idx(ptl_idx, 1):ana_opt.ptl_chunk_idx(ptl_idx, 2), :), hold_betas, preprocessed_data,...
				ana_opt.distribution, ana_opt.dist_specific_params, ana_opt.data_matrix_columns);

                w(ana_opt.ptl_chunk_idx(ptl_idx, 1):ana_opt.ptl_chunk_idx(ptl_idx, 2)) = output_struct.w; % Gather weights

		net_effects(:, ana_opt.ptl_chunk_idx(ptl_idx, 1):ana_opt.ptl_chunk_idx(ptl_idx, 2)) = output_struct.net_effects; % Gather predictor variable
		if ptl_idx == 1, dependent_var = output_struct.dependent_var; end % Gather dependent variable only once, since it is the same across all ptl_idx
	end

	clear output_struct;
	if any(isnan(w(:))), error('NaNs in normalized weight vector w!'); end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Compute the p(theta) and q(theta) weights
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if em > 1
		[p_theta_minus_q_theta] = compute_weights(ana_opt.curve_type, ana_opt.particles, normalized_w(em-1, :), prev_iter_curve_param, param, ana_opt.wgt_chunks, ana_opt.resolution);
		w = w + p_theta_minus_q_theta;
	end

	w = exp(w - logsumexp(w, 2)); % Normalize the weights using logsumexp to avoid numerical underflow
	normalized_w(em, :) = w; % Store the normalized weights

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Optimize betas using fminunc
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	options = optimset('Display', 'iter', 'largescale', 'off', 'gradobj', 'on');
	optimizing_function = family_of_distributions(ana_opt.distribution, 'fminunc_both_betas', w, net_effects, dependent_var, ana_opt.dist_specific_params);

	[hold_betas, fval] = fminunc(optimizing_function, hold_betas, options);

	exp_max_fval(em, 1) = fval; % gather the fval over em iterations
end

hold_betas_per_iter(em+1, :) = hold_betas; % Store away the last em iteration betas
disp(sprintf('>>>>>>>>> Final Betas: %0.4f, %0.4f <<<<<<<<<', hold_betas(1), hold_betas(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flipping the vertical curve parameters if beta_1 is negative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
importance_sampler.flip = false;
neg_beta_idx = hold_betas(2) < 0;
if neg_beta_idx
        disp(sprintf('!!!!!!!!!!!!!!!!!!!! Beta 1 is flipped !!!!!!!!!!!!!!!!!!!!'));
        hold_betas(2) = hold_betas(2) .* -1;
	param = common_to_all_curves(ana_opt.curve_type, 'flip_vertical_params', param);
        importance_sampler.flip = true;
end

w = NaN(1, ana_opt.particles); % Clearing the weight vector
w_null_hypothesis = NaN(1, ana_opt.particles); % Used for a likelihoods ratio test to see if our beta1 value is degenerate
null_hypothesis_beta = [hold_betas(1) 0]; % The null hypothesis for the likelihoods ratio test states that our model y_hat = beta_0 + beta_1 * predictor variable is no different than the simpler model y_hat = beta_0 + beta_1 * predictor variable WHERE BETA_1 = ZERO i.e. our model is y_hat = beta_0

for ptl_idx = 1:size(ana_opt.ptl_chunk_idx, 1)
	output_struct = family_of_curves(ana_opt.curve_type, 'compute_likelihood', ana_opt.net_effect_clusters, ana_opt.ptl_chunk_idx(ptl_idx, 3),...
			param(ana_opt.ptl_chunk_idx(ptl_idx, 1):ana_opt.ptl_chunk_idx(ptl_idx, 2), :), hold_betas, preprocessed_data,...
			ana_opt.distribution, ana_opt.dist_specific_params, ana_opt.data_matrix_columns);
	w(ana_opt.ptl_chunk_idx(ptl_idx, 1):ana_opt.ptl_chunk_idx(ptl_idx, 2)) = output_struct.w;

end

% this code computes the log likelihood of the data under the null hypothesis i.e. using null_hypothesis_beta instead of hold_betas -- it's "lazy" because, unlike the alternative hypothesis, we don't have to compute the data likelihood for each particle because it's exactly the same for each particle (b/c compute_likelihood uses z = beta_1 * x + beta_0, but (recall that our particles control the value of x in this equation) beta_1 is zero for the null hypothesis) that's why we pass in the zero vector representing a single particle with irrelevant weights so we don't have to do it for each particle unnecessarily
output_struct_null_hypothesis_lazy = family_of_curves(ana_opt.curve_type, 'compute_likelihood', ana_opt.net_effect_clusters, 1, [0 0 0 0 0 0], null_hypothesis_beta, preprocessed_data, ana_opt.distribution, ana_opt.dist_specific_params, ana_opt.data_matrix_columns);
	
data_likelihood_null_hypothesis = output_struct_null_hypothesis_lazy.w;

data_likelihood_alternative_hypothesis = w;

w = w + p_theta_minus_q_theta;
if any(isnan(w(:))), error('NaNs in normalized weight vector w!'); end

w = exp(w - logsumexp(w, 2)); % Normalize the weights using logsumexp to avoid numerical underflow
normalized_w(em+1, :) = w; % Store the normalized weights

% Added for debugging chi-sq, might remove eventually
importance_sampler.data_likelihood_alternative_hypothesis = data_likelihood_alternative_hypothesis;
importance_sampler.data_likelihood_null_hypothesis = data_likelihood_null_hypothesis;

% we calculate the data_likelihood over ALL particles by multiplying the data_likelihood for each particle by that particle's importance weight
[dummy_var, importance_sampler.likratiotest] = likratiotest( w * data_likelihood_alternative_hypothesis', data_likelihood_null_hypothesis, 2, 1);

if any(isnan(normalized_w(:))), error('NaNs in normalized weights vector!'); end
if any(isnan(exp_max_fval(:))), error('NaNs in Expectation maximilzation fval matrix!'); end
if any(isnan(hold_betas_per_iter(:))), error('NaNs in hold betas matrix!'); end

importance_sampler.normalized_weights = normalized_w;
importance_sampler.exp_max_fval = exp_max_fval;
importance_sampler.hold_betas_per_iter = hold_betas_per_iter;
importance_sampler.curve_params = param;
importance_sampler.analysis_settings = ana_opt;

if ana_opt.bootstrap
	save(sprintf('%s/%s_b%d_importance_sampler.mat', ana_opt.target_dir, ana_opt.analysis_id, ana_opt.bootstrap_run), '-struct', 'importance_sampler');
elseif ana_opt.scramble
	save(sprintf('%s/%s_s%d_importance_sampler.mat', ana_opt.target_dir, ana_opt.analysis_id, ana_opt.scramble_run), '-struct', 'importance_sampler');
else
	save(sprintf('%s/%s_importance_sampler.mat', ana_opt.target_dir, ana_opt.analysis_id), '-struct', 'importance_sampler');
end
disp(sprintf('Results are stored in be stored in %s', ana_opt.target_dir));

time = clock;
disp(sprintf('Finish time %d/%d %d:%d', time(2), time(3), time(4), time(5)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[p_theta_minus_q_theta] = compute_weights(curve_name, nParticles, normalized_w, prev_iter_curve_param, param, wgt_chunks, resolution)

% [P_THETA_MINUS_Q_THETA] = COMPUTE_WEIGHTS(CURVE_NAME, NPARTICLES, NORMALIZED_W, PREV_ITER_CURVE_PARAM, PARAM, WGT_CHUNKS)
% 
% Purpose
% 
% To compute (P_theta - Q_theta)
%  
% Input
%
% --curve_name: Name of the family of curves (explicitly passed in)
% --nParticles: Number of particles to be used (explicitly passed in)
% --normalized_w: Previous iteration's normalized weights
% --prev_iter_curve_param: Curve parameters held for the previous iteration
% --param: Curve parameters held for the current iteration
% --wgt_chunks: Size of chunk. Purely for limited RAM purposes we break up the giant matrix into smaller matrices to compute the weights (explicitly passed in)
% --resolution: Resolution to which the activations are rounded of
%
% Output
%
% --p_theta_minus_q_theta: Vector of length P (particles)
% 

if nargin < 7, error('Missing input aruguments!'); end

global which_param

total_vol = common_to_all_curves(curve_name, 'curve_volumes', resolution); % Get the curve volumes (Lesbesgue measure)
nParam = family_of_curves(curve_name, 'get_nParams'); % Get the number of curve parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing q(theta), i.e. what is the probability of a curve given all curves from the previous iteration
% P(theta|old_theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_theta = zeros(nParticles, 1);
reduced_nParticles = nParticles / wgt_chunks;
reduced_nParticles_idx = [1:reduced_nParticles:nParticles; (1:reduced_nParticles:nParticles)+reduced_nParticles-1];

for idx = 1:size(reduced_nParticles_idx, 2)
	prob_grp_lvl_curve = zeros(nParticles, reduced_nParticles);
	target_indices = reduced_nParticles_idx(1, idx):reduced_nParticles_idx(2, idx);
	for npm = 1:nParam
		which_param = npm;
		nth_grp_lvl_param = repmat(param(:, npm), 1, reduced_nParticles);
		nth_prev_iter_curve_param = prev_iter_curve_param(target_indices, npm)';
		prob_grp_lvl_curve = bsxfun(@plus, prob_grp_lvl_curve, bsxfun(@compute_trunc_likes, nth_grp_lvl_param, nth_prev_iter_curve_param));
	end
	if any(isnan(prob_grp_lvl_curve(:))), error('NaNs in probability of group level curves matrix!'); end
	q_theta = bsxfun(@plus, q_theta, (exp(prob_grp_lvl_curve) * normalized_w(target_indices)'));
	clear prob_grp_lvl_curve; clear target_indices; clear nth_grp_lvl_param; clear nth_prev_iter_curve_param;
end
if any(isnan(q_theta(:))), error('NaNs in q_theta vector!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing p(theta) prior i.e. what is the probability of a curve in the curve space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_theta = ones(nParticles, 1);
p_theta = bsxfun(@times, p_theta, (1 / total_vol));
if length(unique(p_theta)) ~= 1, error('p_theta is NOT unique!'); end
if any(isnan(p_theta(:))), error('NaNs in p_theta vector!'); end

p_theta_minus_q_theta = log(p_theta)' - log(q_theta)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[log_likelihood] = compute_trunc_likes(x, mu)

% [LOG_LIKELIHOOD] = COMPUTE_TRUNC_LIKES(X, MU)
% 
% Purpose
% 
% To compute the log likelihood of truncated normal distribution
% 
% Input
%
% --x: Kth curve parameter in the current iteration
% --mu: All curve parameters from the previous iteration
% --tau: Similar to sigma in Gaussian distribution (via global)
% --bounds and which_param: Are used to fetch the respective curve parameter's absolute abounds (via global). This is required since we are computing likelihoods for truncated normal
% NOTE: tau and bounds are NOT modified throughout this code once initially set 
% 
% Output
%
% --log_likelihood of that curve given all the P curves from the previous iteration
%

if nargin < 2, error('Missing input aruguments!'); end

global tau; global bounds; global which_param;

if tau <= 0, error('Tau is <= 0!'); end

% This ugly thing below is a manifestation of log(1 ./ (tau .* (normcdf((bounds(which_param, 2) - mu) ./ tau) - normcdf((bounds(which_param, 1) - mu) ./ tau))) .* normpdf((x - mu) ./ tau))
% Refer to http://en.wikipedia.org/wiki/Truncated_normal_distribution for the truncated normal distribution
log_likelihood = -(log(tau) + log((0.5 .* erfc(-(((bounds(which_param, 2)-mu) ./ (tau.*sqrt(2)))))) - (0.5 .* erfc(-(((bounds(which_param, 1)-mu) ./ (tau .* sqrt(2)))))))) + ...
						(-0.5 .* (log(2) + log(pi)) - (0.5 .* ((x - mu) ./ tau) .^ 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
