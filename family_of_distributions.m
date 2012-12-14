function[output] = family_of_distributions(distribution_name, get_info, varargin)

% [OUTPUT] = FAMILY_OF_DISTRIBUTIONS(DISTRIBUTION_NAME, GET_INFO, VARARGIN)
%
% Purpose
%
% For each of the family of distributions this script performs specific computations like number of pdf/pmf, etc
%
% Input
%
% --distribution_name: distribution name, string, e.g. 'bernoulli', 'normal'
% --get_info: Cues for specific information / computation, string, e.g. 'get_nParams'
% --varargin: Is either empty or has arguments depending on the computation
%
% Output
%
% --output: Holds the output of all computations
%
% Note: If you need to create a new family of distributions, make sure to have all functions with '--> (*)' replicated for the new family of distributions.
% The case id's will need to match as well
%
% Example usage:
%
% family_of_distributions('normal', 'compute_densities', varargin)
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

% Checks if the correct number of arguments are passed in
if nargin < 2, error('Missing input parameters'); end

switch distribution_name
	case 'bernoulli', 	output = bernoulli_distribution(get_info, varargin);
	case 'normal', 	        output = normal_distribution(get_info, varargin);
	otherwise, 		error('Invalid distribution!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = bernoulli_distribution(get_info, input_params)

	switch get_info
		case 'compute_densities' % --> (1), Compute the log densities. NOTE: We compute the log(probability function)
			if length(input_params) <= 1, error('Missing input parameters!'); end

			z = input_params{1};
			y = input_params{2};
			clear input_params;
	
			% Compute fz = 1 / (1 + exp(-z) - Logistic function
			fz = 1 ./ (1 + exp(-z));
			clear z;
			fz = max(fz, eps);
			fz = min(fz, 1-eps);

			% Compute bern_log_pmf = p ^ k + (1 - p) ^ (1 - k). http://en.wikipedia.org/wiki/Bernoulli_distribution
			% Here p = fz and k = y. Taking the log results in y x log(fz) + (1 - y) x log(1 - fz). This is written below in bsxfun syntax
			out = sum(bsxfun(@times, log(fz), y) + bsxfun(@times, log(1 - fz), bsxfun(@minus, 1, y)));

		case 'fminunc_both_betas' % --> (2), This fetches the right function handle for the fminunc
			if length(input_params) <= 2, error('Missing input parameters!'); end
			out = @(betas)fminunc_bernoulli_both(betas, input_params{1}, input_params{2}, input_params{3});

		otherwise
			error('Invalid operation!');
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Note: this function will need to live in this space as a nested function since I am settimg up a handle for it
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function [f, g] = fminunc_bernoulli_both(betas, w, net_effects, dependent_var)

		% [F, G] = FMINUNC_BERNOULLI_BOTH(BETAS)
		% 
		% Purpose
		% 
		% To optimize logistic regression betas using cost function F
		%  
		% Input
		%
		% --betas: The current betas that were used to compute likelihoods
		% --w: Weight vector that holds the normalized weights for P particles
		% --net_effects: Predictor variable Matrix (number of trials x particles)
		% --dependent_var: Dependent variable Matrix (number of trials x 1)
		% 
		% Output
		%
		% --f: Scalar, Objective function
		% --g: Vector of length 2 i.e. gradients with respect to beta_0 and beta_1
		%

		if nargin < 4, error('Missing input aruguments!'); end

		beta_0 = betas(1);
		beta_1 = betas(2);

		z = beta_1 .* net_effects + beta_0;
		fz = 1 ./ (1 + exp(-z));
		if any(isinf(fz(:))), error('Inf in fz matrix!'); end
		fz = max(fz, eps);
		fz = min(fz, 1-eps);

		% Cost function
		% We will need to maximize the betas but fminunc minimizes hence a -ve.
		% Here we compute the log pmf over all trials and then component multiply by the weights and then sum them up over all particles
		f = -sum(w .* sum(bsxfun(@times, log(fz), dependent_var) + bsxfun(@times, log(1 - fz), bsxfun(@minus, 1, dependent_var))));

		% Here we take the partial derivative of log pmf over beta_0 and beta_1 respectively, component multiply by the weights and sum them up over all paricles
		g(1) = -sum(w .* sum(bsxfun(@minus, dependent_var, (exp(z) ./ (1 + exp(z))))));
		g(2) = -sum(w .* sum(bsxfun(@times, net_effects, dependent_var) - ((net_effects .* exp(z)) ./ (1 + exp(z)))));
		if any(isinf(g(:))), error('Inf in partial derivative!'); end
		if any(isnan(g(:))), error('NaN in partial derivative!'); end
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = normal_distribution(get_info, input_params)

	switch get_info
		case 'compute_densities' % --> (1), Compute the log densities. NOTE: We compute the log(probability function)
			if length(input_params) <= 2, error('Missing input parameters!'); end

			mu = input_params{1};
			y = input_params{2};
			dist_specific_params = input_params{3};
			clear input_params;
			sigma = dist_specific_params.sigma;

			% Compute log_pdf http://en.wikipedia.org/wiki/Normal_distribution
			out = sum(bsxfun(@minus, ((1 ./ sigma .^ 2) .* (bsxfun(@minus, bsxfun(@times, y, mu), bsxfun(@plus, bsxfun(@times, 0.5, bsxfun(@power, mu, 2)),...
					                 bsxfun(@times, 0.5, bsxfun(@power, y, 2)))))), (0.5 .* log(2 .* pi .* sigma .^ 2))));

		case 'fminunc_both_betas' % --> (2), This fetches the right function handle for the fminunc
			if length(input_params) <= 3, error('Missing input parameters!'); end
			out = @(betas)fminunc_normal_both(betas, input_params{1}, input_params{2}, input_params{3}, input_params{4});

		otherwise
			error('Invalid operation!');
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Note: this function will need to live in this space as a nested function since I am settimg up a handle for it
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function [f, g] = fminunc_normal_both(betas, w, net_effects, dependent_var, dist_specific_params)

		% [F, G] = FMINUNC_NORMAL_BOTH(BETAS)
		% 
		% Purpose
		% 
		% To optimize logistic regression betas using cost function F
		%  
		% Input
		%
		% --betas: The current betas that were used to compute likelihoods
		% --w: Weight vector that holds the normalized weights for P particles
		% --net_effects: Predictor variable Matrix (number of trials x particles)
		% --dependent_var: Dependent variable Matrix (number of trials x 1)
		% --sigma: Used to specify variance in the Normal distribution
		% 
		% Output
		%
		% --f: Scalar, Objective function
		% --g: Vector of length 2 i.e. gradients with respect to beta_0 and beta_1
		%

		if nargin < 5, error('Missing input aruguments!'); end

		beta_0 = betas(1);
		beta_1 = betas(2);
		sigma = dist_specific_params.sigma;

		mu = beta_1 .* net_effects + beta_0;

		% Cost function
		% We will need to maximize the betas but fminunc minimizes hence a -ve.
		% Here we compute the log pdf over all trials and then component multiply by the weights and then sum them up over all particles
		f = -sum(w .* sum(bsxfun(@minus, bsxfun(@times, (1 ./ sigma .^ 2), bsxfun(@minus, bsxfun(@times, dependent_var, mu), bsxfun(@plus,...
			bsxfun(@times, 0.5, bsxfun(@power, mu, 2)), bsxfun(@times, 0.5, bsxfun(@power, dependent_var, 2))))), (0.5 .* log(2 .* pi .* sigma .^ 2)))));

		% Here we take the partial derivative of log pdf over beta_0 and beta_1 respectively, component multiply by the weights and sum them up over all paricles
		g(1) = -sum(w .* sum((1 ./ sigma .^ 2) .* bsxfun(@minus, dependent_var, bsxfun(@plus, beta_0, bsxfun(@times, beta_1, net_effects)))));
		g(2) = -sum(w .* sum(bsxfun(@times, bsxfun(@rdivide, net_effects, (sigma .^ 2)), bsxfun(@minus, dependent_var,...
									bsxfun(@plus, beta_0, bsxfun(@times, beta_1, net_effects))))));

		if any(isinf(g(:))), error('Inf in partial derivative!'); end
		if any(isnan(g(:))), error('NaN in partial derivative!'); end
	end
end
end
