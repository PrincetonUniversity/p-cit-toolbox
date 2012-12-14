function[out] = common_to_all_curves(curve_type, get_info, varargin)

% [OUT] = COMMON_TO_ALL_CURVES(CURVE_TYPE, GET_INFO, VARARGIN)
%
% Purpose
%
% This function fetches information that is common to all curves
%
% Input
%
% --curve_type: Family of curves, string, e.g. 'free_rmpw'
% --get_info: Cues for specific information / computation, string, e.g. 'initial_sampling'
% --varargin: Has arguments depending on the computation
%
% Output
%
% --The output of all computations sit in out
%
% Example usage:
%
% common_to_all_curves('horz_indpnt', 'initial_sampling', 1000, 4)
% common_to_all_curves('horz_indpnt', 'check_if_exceed_bounds', some_matrix)
% common_to_all_curves('horz_indpnt', 'curve_volumes', 5)
% common_to_all_curves('horz_indpnt', 'flip_vertical_params', some_matrix)
% common_to_all_curves('horz_indpnt', 'sort_horizontal_params', some_matrix)
% common_to_all_curves('horz_indpnt', 'draw_bcm_curves', [0.2, 0.5, 0.1, 0.1, 0.1, 0.1], 4)
% common_to_all_curves('horz_indpnt', 'auto_generate', 'con', 4)
% common_to_all_curves('horz_indpnt', 'weighted_curve', importance_sampler_mat, 0.9, 4)
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

% Checks if input arguments are passed in
if isempty(varargin)
	error('No input arguments!');
end

switch get_info
	case 'initial_sampling'
		if length(varargin) < 2
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		nParticles = varargin{1};
		if nParticles <= 0
			error(sprintf('Number of particles will need to > 0!'));
		end
		resolution = varargin{2};
		if resolution <= 0
			error(sprintf('Resolution will need to > 0!'));
		end

		bounds = family_of_curves(curve_type, 'get_bounds');
		nParams = family_of_curves(curve_type, 'get_nParams');
		out = NaN(nParticles, nParams);
	
		% Uniform sampling each curve parameter bounded by its respective bounds
		for i = 1:nParams
			out(:, i) = unifrnd(bounds(i, 1), bounds(i, 2), nParticles, 1);
		end

		out = round_to(out, resolution);
		if any(isnan(out(:))), error('NaNs in initial sampling output matrix!'); end

	case 'check_if_exceed_bounds'
		if length(varargin) < 1
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		out = varargin{1};
		nParams = family_of_curves(curve_type, 'get_nParams');
		if isempty(out) | size(out, 2) ~= nParams
			error(sprintf('Not a valid input matrix!'));
		end

		bounds = family_of_curves(curve_type, 'get_bounds');
		nParams = family_of_curves(curve_type, 'get_nParams');

		% If a curve parameter is found to exceeding bounds then it is set to the bounds
		% For instance if a vertical parameter is -1.02 then it is set to -1 since -1 is the lower bound for vertical parameters
		for i = 1:nParams
			out(:, i) = max(out(:, i), bounds(i, 1));
			out(:, i) = min(out(:, i), bounds(i, 2));
		end

	case 'curve_volumes'
		if length(varargin) < 1
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		resolution = varargin{1};
		if resolution <= 0
			error(sprintf('Resolution will need to > 0!'));
		end

		bounds = family_of_curves(curve_type, 'get_bounds');
		nParams = family_of_curves(curve_type, 'get_nParams');
		total_vol = 1;

		% Lebesgue measure http://en.wikipedia.org/wiki/Lebesgue_measure
		for i = 1:nParams
			total_vol = total_vol * length(bounds(i, 1):(1/10^resolution):bounds(i, 2));
		end
		out(1) = total_vol;

	case 'flip_vertical_params'
		if length(varargin) < 1
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		input_params = varargin{1};
		nParams = family_of_curves(curve_type, 'get_nParams');
		if isempty(input_params) | size(input_params, 2) ~= nParams
			error(sprintf('Not a valid input matrix!'));
		end
		out = input_params;
		vertical_params = family_of_curves(curve_type, 'get_vertical_params_only');

		% Flipping vertical parameters of the curve. If a y1 = -0.4, flipping it will result in 0.4
		for i = 1:length(vertical_params)
			out(:, vertical_params(i)) = bsxfun(@times, input_params(:, vertical_params(i)), -1);
		end

	case 'sort_horizontal_params'
		if length(varargin) < 1
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		input_params = varargin{1};
		nParams = family_of_curves(curve_type, 'get_nParams');
		if isempty(input_params) | size(input_params, 2) ~= nParams
			error(sprintf('Not a valid input matrix!'));
		end
		out = input_params;
		horizontal_params = family_of_curves(curve_type, 'get_horizontal_params_only');
		if length(horizontal_params) ~= 2
			error('Incorrect horizontal parameters count for %s family of curves', curve_type);
		end

		% This piece of code ensures that x1 <= x2 especially for the horz_indpnt family of curves
		idx = input_params(:, horizontal_params(1)) > input_params(:, horizontal_params(2));
		out(idx, horizontal_params(1)) = input_params(idx, horizontal_params(2));
		out(idx, horizontal_params(2)) = input_params(idx, horizontal_params(1));

		if ~all(out(:, horizontal_params(1)) <= out(:, horizontal_params(2)))
			error('Horizontal parameter 1 is NOT <= Horizontal parameter 2 in %s family of curves', curve_type);
		end

	case 'draw_bcm_curve'
		if length(varargin) < 2
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		input_params = varargin{1};
		resolution = varargin{2};
		if resolution <= 0
			error(sprintf('Resolution will need to > 0!'));
		end
		% This draws a BCM curve for you. If you passed in the 'input_params' as 'con' then it randomly draws a theory consistent curve; 'inc' - theory inconsistent curve
		if strcmp(input_params, 'con') | strcmp(input_params, 'inc')
			input_params = common_to_all_curves(curve_type, 'auto_generate', input_params, resolution);
		end
		nParams = family_of_curves(curve_type, 'get_nParams');
		if isempty(input_params) | size(input_params, 2) ~= nParams
			error(sprintf('Not a valid input matrix!'));
		end

		% If instead you passed in [y1, x1, x2, y2, y3 and y4] into 'input_params' then it draws a curve directly rather then randomly generating one for you
		out = family_of_curves(curve_type, 'get_curve_xy_vals', input_params);

		figure(); set(gcf, 'Position', [800, 900, 600, 600]);
		plot(out.xval, out.yval, 'LineWidth', 3); grid on;
		xlabel('Activation', 'FontSize', 18); ylabel('Change in Memory Strength', 'FontSize', 18);
		ylim([-1.2, 1.2]); title(out.title_string);

	case 'auto_generate'
		if length(varargin) < 2
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		input_params = varargin{1};
		resolution = varargin{2};
		if resolution <= 0
			error(sprintf('Resolution will need to > 0!'));
		end
		nSamples = 100;
		nParam = family_of_curves(curve_type, 'get_nParams');
		params = NaN(nSamples, nParam);
		out = NaN(1, nParam);

		% Generate 100 curves and randomly pick a theory consistent or inconsistent curve depending on the request
		params = common_to_all_curves(curve_type, 'initial_sampling', nSamples, resolution);
		if strcmp(curve_type, 'horz_indpnt') % Enforce the right ordering for the horizontal curve parameters i.e. x1 < x2
			params = common_to_all_curves(curve_type, 'sort_horizontal_params', params);
		end

		if any(isnan(params(:))), error('NaNs in curve parameter matrix!'); end
		params_indices = family_of_curves(curve_type, 'count_particles', params);

		if strcmp(input_params, 'con')
			th_con_params_indices = [find(params_indices)]; % Finding the theory consistent trial indices
			if length(th_con_params_indices) <= 0, error('Did not generate any theory consistent indices!'); end
			th_con_params_indices = th_con_params_indices(randperm(size(th_con_params_indices, 1))); % Randomly permuting the th_con trial indices
			out = params(th_con_params_indices(1), :); % picking one consistent particle
		elseif strcmp(input_params, 'inc')
			th_inc_params_indices = [find(~params_indices)]; % Finding the theory inconsistent trial indices
			if length(th_inc_params_indices) <= 0, error('Did not generate any theory inconsistent indices!'); end
			th_inc_params_indices = th_inc_params_indices(randperm(size(th_inc_params_indices, 1))); % Randomly permuting the th_inc trial indices
			out = params(th_inc_params_indices(1), :); % picking one inconsistent particle
		else
			error(sprintf('Invalid string! valid ones include ''con'' or ''inc'' only'));
		end
		if any(isnan(out(:))), error('NaNs in curve parameters!'); end

	case 'weighted_curve'
		% Get the weighted curve and the associated credible interval
		if length(varargin) < 3
			error(sprintf('Missing input parameters in %s computation!', get_info));
		end
		% importance_sampler_mat: .mat file
		importance_sampler_mat = varargin{1};
		% credible_interval: the credible interval http://en.wikipedia.org/wiki/Credible_interval
		credible_interval = varargin{2};
		if credible_interval <= 0
			error(sprintf('Resolution will need to > 0!'));
		end
		resolution = varargin{3};
		if resolution <= 0
			error(sprintf('Resolution will need to > 0!'));
		end

		% We the fetch the x and y vals for each of the particles. For instance if there are 100 particles we the x and y vals for 100 different curves
		out = family_of_curves(importance_sampler_mat.analysis_settings.curve_type, 'get_curve_xy_vals', importance_sampler_mat.curve_params, resolution);

		% Getting the y-axis to plot - final em iteration weights(1 x particles) * out.yval (particles x Y vals).
		% Note this is matrix multiplication so we should get a 1 x y vals matrix, which is our final y vals
		y_final = importance_sampler_mat.normalized_weights(end, :) * out.yval;

		% Getting the x-axis to plot
		x_final = out.xval(1, :);

		% Computing the credible interval envelope
		% The way this works is:
		% 1. For each x value take all the associated y values
		% 2. Sort the y values and get the associated sorted weights. For instance the y values look like 0.8,  0.7, -0.1,  -0.5, etc and the associated weights
		% 											look like 0.04, 0.2,  0.3, 0.003, etc
		% 3. Take the cumulative sum of weights. For the example above this will look like 0.04, 0.24, 0.54, 0.543, etc
		% 4. Identify the 0.05 weight location and the 0.95 weight location. This ensures that 90% the weights are concentrated within those bounds

		interval = NaN(2, length(x_final));
		for c = 1:size(out.yval, 2)
			[sorted_yval, sorted_yval_idx] = sort(out.yval(:, c));
			sorted_wgts = importance_sampler_mat.normalized_weights(end, sorted_yval_idx);
			cumsum_wgts = cumsum(sorted_wgts);
			lower_bnd = find(cumsum_wgts >= (1 - credible_interval) / 2); % i.e. 0.0500
			upper_bnd = find(cumsum_wgts <= (1 - (1 - credible_interval) / 2)); % i.e. 0.9500
			interval(1, c) = sorted_yval(lower_bnd(1));
			interval(2, c) = sorted_yval(upper_bnd(end));
		end
		assert(~any(isnan(interval(:))));

		% The x and y vals for the blue line and the coordinates of the gray envelope
		out = struct();
		out.x_final = x_final;
		out.y_final = y_final;
		out.interval = interval;

	otherwise
		error('Invalid operation!');
end

