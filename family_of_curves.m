function[output] = family_of_curves(curve_type, get_info, varargin)

% [OUTPUT] = FAMILY_OF_CURVES(CURVE_TYPE, GET_INFO, VARARGIN)
%
% Purpose
%
% For each of the family of curves this script has curve specific computations like number of curve parameters, boundaries, log likelihood, etc
%
% Input
%
% --curve_type: Family of curves, string, e.g. 'horz_indpnt'
% --get_info: Cues for specific information / computation, string, e.g. 'get_nParams'
% --varargin: Is either empty or has arguments depending on the computation
%
% Output
%
% --output: Holds the output of all computations
%
% Note: If you need to create a new family of curves, make sure to have all functions with '--> (*)' replicated for the new family of curves. The case id's will need to match as well
%
% Example usage:
%
% family_of_curves('horz_indpnt', 'get_bounds')
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

switch curve_type
	case 'horz_indpnt', 	output = horz_indpnt_curve(get_info, varargin);
	otherwise, 		error('Invalid curve!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = horz_indpnt_curve(get_info, input_params)

% Order of curve parameters y1, x1, x2, y2, y3 and y4
% Note always x1 will need to precede x2 (both when passing in curve parameters as well as when plotting)

nParam = 6;
curve_type = 'horz_indpnt';

switch get_info
	case 'get_nParams' % --> (1), number of curve parameters
		out = nParam;

	case 'get_bounds' % --> (2), Absolute boundaries of the curve parameters
		out = [-1 1; 0 1; 0 1; -1 1; -1 1; -1 1]; % These bounds will always need to be absolute bounds; it does not vary with consistency or inconsistency

	case 'get_vertical_params_only' % --> (3), Indices of vertical parameters i.e. the ones corresponding to v* in the order of curve parameters
		out = [1, 4, 5, 6];

	case 'get_horizontal_params_only' % --> (4), Indices of horizontal parameters i.e. the ones corresponding to h* in the order of curve parameters
		out = [2, 3];

	case 'compute_likelihood' % --> (5), Get the curve y-vals (with or without net effects) for each of the P particle curves and then compute the log probability mass function (pmf)
		if length(input_params) <= 5, error('Missing input parameters!'); end

		net_effect_clusters = input_params{1};
		particles = input_params{2};
		y1 = input_params{3}(:, 1);
		x1 = input_params{3}(:, 2);
		x2 = input_params{3}(:, 3);
		y2 = input_params{3}(:, 4);
		y3 = input_params{3}(:, 5);
		y4 = input_params{3}(:, 6);
		b0 = input_params{4}(1);
		b1 = input_params{4}(2);
		data = input_params{5};
		distribution = input_params{6};
		dist_specific_params = input_params{7};

		data_matrix_columns = input_params{8};
		predictor_var_column = data_matrix_columns.predictor_var;
		dependent_var_column = data_matrix_columns.dependent_var;
		net_effect_clusters_column = data_matrix_columns.net_effect_clusters;

		clear input_params;
		if ~all(all(x1 <= x2)), error('Horizontal parameter 1 is NOT <= Horizontal parameter 2 in %s family of curves', curve_type); end

		x = NaN(length(net_effect_clusters), particles);
		y = []; % This vector cannot be initialized in advance since we don't know the length of this vector

		% In this loop we map the predictor variables to the associated y vals for all curves / particles simulataneously
		for i = 1:length(net_effect_clusters)
			cluster_idx = find(data(:, net_effect_clusters_column) == net_effect_clusters(i));
			X = zeros(length(cluster_idx), particles);
			for j = 1:length(cluster_idx)
				if isnan(data(cluster_idx(j), predictor_var_column))
					x(i, :) = 0;
				else
					% If an activation is falling in the third segment of the curve then get the associated y val
					ix3 = data(cluster_idx(j), predictor_var_column) > x2;
					X(j, ix3) = (((y4(ix3) - y3(ix3)) ./ (1 - x2(ix3))) .* (data(cluster_idx(j), predictor_var_column) - 1)) + y4(ix3);

					% If an activation is falling in the second segment of the curve then get the associated y val
					ix2 = ~ix3 & (data(cluster_idx(j), predictor_var_column) > x1); % segment #2
					X(j, ix2) = (((y3(ix2) - y2(ix2)) ./ (x2(ix2) - x1(ix2))) .* (data(cluster_idx(j), predictor_var_column) - x1(ix2))) + y2(ix2);

					% If an activation is falling in the first segment of the curve then get the associated y val
					ix1 = ~ix3 & ~ix2 & data(cluster_idx(j), predictor_var_column) > 0; % segment #1
					X(j, ix1) = (((y2(ix1) - y1(ix1)) ./ x1(ix1)) .* data(cluster_idx(j), predictor_var_column)) + y1(ix1);

					% If an activation is at the intercept of the curve then get the associated y val
					ix0 = ~ix3 & ~ix2 & ~ix1 & data(cluster_idx(j), predictor_var_column) == 0; % Intercept (Boundary condition)
					X(j, ix0) = y1(ix0);

					% If an item has net effects then taking the sum below will compute the net effects.
					% If an item has no net effect then this loop will be executed only once and the sum has no effect
					x(i, :) = sum(X, 1);
				end
			end

			% Our model enforces that the dependent variable will need to be unique for items within a net effect cluster i.e. all 1's or all 0's
			if length(unique(data(cluster_idx, dependent_var_column))) ~= 1
				error('Dependent var is NOT unique for net effect cluster %d', i);
			end
			% We accumulate the dependent variables for each net effect cluster
			y = [y; unique(data(cluster_idx, dependent_var_column))];
		end

		clear X; clear ix0; clear ix1; clear ix2; clear ix3; clear x1; clear x2; clear y1; clear y2; clear y3; clear y4;
		clear data; clear data_matrix_columns;
		if any(isnan(x(:))), error('NaNs in trials x particles matrix!'); end
		if any(isinf(x(:))), error('Inf in trials x particles matrix!'); end

		% Compute z = beta_0 + beta_1 x predictor variable
		z = b1 .* x + b0;

		out = struct();
		out.w = family_of_distributions(distribution, 'compute_densities', z, y, dist_specific_params);
		out.net_effects = x;
		out.dependent_var = y;

	case 'count_particles' % --> (6), Use some criterion to carve out the curve space into theory consistent and theory inconsistent
		if length(input_params) <= 0, error('Missing input parameters!'); end
		if ~all(all(input_params{1}(:, 2) <= input_params{1}(:, 3))), error('Horizontal parameter 1 is NOT <= Horizontal parameter 2 in %s family of curves', curve_type); end

		%{
		The different branches of theory-consistent curves can be defined in terms of the location of the dip in the curve (the point that anchors the
		weakening part of the curve) and the rise in the curve (the point that anchors the strengthening part of the curve). More formally,

		-- The dip in a theory-consistent curve is a point that is located horizontally between the left edge of the curve and the rise.
		Within this horizontal range, the dip is the lowest point on the curve; it also has to fall below zero on the y-axis.

		-- The rise in a theory-consistent curve is a point that is located horizontally to the right of the dip. Within this horizontal
		range, the rise is the highest point on the curve; it also has to fall above zero on the y-axis.

		Branch I: y2 defines the dip and y3 defines the rise
		-1 <= y2 < 0, y2 is the dip so it must fall below zero
		0 < y3 <= 1, y3 is the rise so it must fall above zero
		-1 <= y4 <= y3, y4 can hold any value that is below the rise (y3)
		y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)

		Branch II: y2 defines the dip and y4 defines the rise
		-1 <= y2 < 0, y2 is the dip so it must fall below zero
		0 < y4 <= 1, y4 is the rise so it must fall above zero
		y2 <= y3 <= y4, y3 can hold any value between the dip and the rise
		y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)

		Branch III: y3 defines the dip and y4 defines the rise
		-1 <= y3 < 0, y3 is the dip so it must fall below zero
		0 < y4 <= 1, y4 is the rise so it must fall above zero
		y3 < y1 <= 1, y1 can hold any value that is above the dip (y3)
		y3 <= y2 <= 1, y2 can hold any value that is above the dip (y3)
		%}

		% Fetch the indices of the theory consistent particles
		% First two lines ensure that the horizontal parameters cannot be 0 OR 1, since that would eliminate a line segment altogether
		out = input_params{1}(:, 2) ~= 0  			& input_params{1}(:, 2) ~= 1 & ...
		      input_params{1}(:, 3) ~= 0  			& input_params{1}(:, 3) ~= 1 & ...
		    ((input_params{1}(:, 4) >= -1 			& input_params{1}(:, 4) < 0 & ... % Branch I
		      input_params{1}(:, 5) > 0 			& input_params{1}(:, 5) <= 1 & ...
		      input_params{1}(:, 6) >= -1			& input_params{1}(:, 6) <= input_params{1}(:, 5) & ...
		      input_params{1}(:, 1) > input_params{1}(:, 4) 	& input_params{1}(:, 1) <= 1) | ...
		     (input_params{1}(:, 4) >= -1 			& input_params{1}(:, 4) < 0 & ... % Branch II
		      input_params{1}(:, 6) > 0 			& input_params{1}(:, 6) <= 1 & ...
		      input_params{1}(:, 5) >= input_params{1}(:, 4) 	& input_params{1}(:, 5) <= input_params{1}(:, 6) & ...
		      input_params{1}(:, 1) > input_params{1}(:, 4) 	& input_params{1}(:, 1) <= 1) | ...
		     (input_params{1}(:, 5) >= -1 			& input_params{1}(:, 5) < 0 & ... % Branch III
		      input_params{1}(:, 6) > 0 			& input_params{1}(:, 6) <= 1 & ...
		      input_params{1}(:, 1) > input_params{1}(:, 5) 	& input_params{1}(:, 1) <= 1 & ...
		      input_params{1}(:, 4) >= input_params{1}(:, 5) 	& input_params{1}(:, 4) <= 1));

	case 'get_curve_xy_vals' % --> (7), This is the same as compute_likelihood in the sense map the predictor variable to the curve y val but there are couple of differences ...
				 % 1. We only know the curve parameters and we have to map all the xvals (0-to-1) to the curve yvals where as in compute_likelihood we had
				 %    specific xvals (predictor variables) and curve parameters
				 % 2. There is NO net effect cluster concept here
				 % 3. We DO NOT compute the pmf as well
				 % Hence parts of the code will look similar but we felt these two chunks of code will need to be separate

		if length(input_params) <= 0, error('Missing input parameters!'); end

		if length(input_params) > 1, resolution = input_params{2};
		else, resolution = 4; end

		particles = size(input_params{1}(:, 1), 1);

		if any(input_params{1}(:, [1, 4, 5, 6]) < -1) | any(input_params{1}(:, [1, 4, 5, 6]) > 1)
			error('Vertical curve parameters exceed bounds [-1, 1]!');
		end

		if any(input_params{1}(:, [2, 3]) < 0) | any(input_params{1}(:, [2, 3]) > 1)
			error('Horizontal curve parameters exceed bounds [0, 1]!');
		end

		xval = 0:(1/10^resolution):1; % This will need to start from zero since we are scaling activations between 0 and 1. So there will definitely be a zero in every analyzes
		xval = repmat(xval, particles, 1);
		yval = NaN(size(xval));
		out = struct();

		y1 = repmat(input_params{1}(:, 1), 1, size(xval, 2));
		x1 = repmat(input_params{1}(:, 2), 1, size(xval, 2));
		x2 = repmat(input_params{1}(:, 3), 1, size(xval, 2));
		y2 = repmat(input_params{1}(:, 4), 1, size(xval, 2));
		y3 = repmat(input_params{1}(:, 5), 1, size(xval, 2));
		y4 = repmat(input_params{1}(:, 6), 1, size(xval, 2));
		if ~all(all(x1 <= x2)), error('Horizontal parameter 1 is NOT <= Horizontal parameter 2 in %s family of curves', curve_type); end

		ix3 = xval > x2; % segment #3
		yval(ix3) = (((y4(ix3) - y3(ix3)) ./ (1 - x2(ix3))) .* (xval(ix3) - 1)) + y4(ix3);

		ix2 = ~ix3 & (xval > x1); % segment #2
		yval(ix2) = (((y3(ix2) - y2(ix2)) ./ (x2(ix2) - x1(ix2))) .* (xval(ix2) - x1(ix2))) + y2(ix2);

		ix1 = ~ix3 & ~ix2 & xval > 0; % segment #1
		yval(ix1) = (((y2(ix1) - y1(ix1)) ./ x1(ix1)) .* xval(ix1)) + y1(ix1);

		ix0 = ~ix3 & ~ix2 & ~ix1 & xval == 0; % Intercept (Boundary condition)
		yval(ix0) = y1(ix0);

		if any(isnan(yval(:))), error('NaNs in trials x particles matrix!'); end
		if any(isinf(yval(:))), error('Inf in trials x particles matrix!'); end
		out.xval = xval;
		out.yval = yval;
		if particles == 1
			out.curve_params = input_params{1};
			out.title_string = sprintf('y1=%0.2f, x1=%0.2f, x2=%0.2f y2=%0.2f, y3=%0.2f, y4=%0.2f', y1(1), x1(1), x2(1), y2(1), y3(1), y4(1));
		end

	otherwise
		error('Invalid operation!');
end

