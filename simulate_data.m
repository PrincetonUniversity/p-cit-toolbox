function [] = simulate_data(analysis_id, noise_knob, curve_type, yval_distribution, net_effects, varargin)

% [] = SIMULATE_DATA(ANALYSIS_ID, NOISE_KNOB, CURVE_TYPE, YVAL_DISTRIBUTION, NET_EFFECTS, VARARGIN)
%
% Purpose:
%
% Generate simulated data from a ground truth curve
% 
% Inputs:
%
% analysis_id: Valid analysis Id. NOTE will need to end in sim
% noise_knob: variance in noise added to activations
% curve_type: valid curve type listed in show_bcm_curve_beta.m
% yval_distribution: distribution of yvals. Right now the code supports 'bernoulli' and 'normal'
% net effects: If > 1, those many repetitions will be sampled per item. If <= 0 then only one repetition per item
% varargin = vector of inputs depending on the curve / 'con' or 'inc'
%
% Outputs:
%
% simulated_data.mat
%
% Example usage:
%
% Note the order of the curve parameters [y1, x1, x2, y2, y3, y4]
% simulate_data('my_analysis_id', 0.001, 'horz_indpnt', 'bernoulli', 10, [0.6, 0.2, 0.5, -0.4, 0.5, 0.9])
% simulate_data('my_analysis_id', 0.001, 'horz_indpnt', 'bernoulli', 10, 'con') - for a consistent curve
% simulate_data('my_analysis_id', 0.001, 'horz_indpnt', 'bernoulli', 10, 'inc') - for an inconsistent curve
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Checks if the correct number of arguments are passed in
if nargin < 6, error('Missing input parameters'); end

% Resetting the random number seed based on the clock
rand('twister', sum(100 * clock));

% Setting the target directory
results_dir = fullfile(pwd, 'results');
if ~exist(results_dir), mkdir(results_dir); end
target_dir = fullfile(results_dir, analysis_id);
if ~exist(target_dir), mkdir(target_dir); end

% Setting the number of subjects
nSubjects = 35;
% Setting the number of samples per subject
nItems = 8;
% Setting the item_spread = 0.15. This ensure that each samples' repetitions i.e. items have atleast 0.15 spread
item_spread = 0.15;
% Setting the resolution i.e. 10 ^ -(resolution)
resolution = 4;

% Draw a curve and get the x, y values and curve parameters
if strcmp(varargin{1}, 'con') | strcmp(varargin{1}, 'inc') % If 'con' OR 'inc' then it randomly generates a consistent or inconsistent curve
	curve_params = common_to_all_curves(curve_type, 'auto_generate', varargin{1}, resolution);
elseif isnumeric(varargin{1}) % If a vector is passed in order [y1, x1, x2, y2, y3, y4] then it generates a curve using these curve parameters. Note the order of parameters is important
	curve_params = varargin{1};
else
	error('Invalid varargin! the valid arguments are ''con'', ''inc'' or [y1, x1, x2, y2, y3, y4]');
end
out = family_of_curves(curve_type, 'get_curve_xy_vals', curve_params);

% House keeping
subj_id_list = [];
item_list = [];
predictor_var_list = [];
dependent_var_list = [];
net_effects_list = [];
net_eff_counter = 1;

for s = 1:nSubjects % For each subject
	subj_obs = [];
	yvals = [];
	if net_effects > 0
		tmp_net_effects_list = [];
		% Uniformly sample 'nItems' items per subject between item_spread and (1- item_spread).
		% Now the [item_spread, (1-item_spread)] boundary case comes from the fact you are going to sample net effects for each item as 'each item +/- item_spread' bin boundary
		subj_act_means = unifrnd(item_spread, (1-item_spread), 1, nItems);
		% For each item sample 'net_effects' number of net effects
		for t = 1:length(subj_act_means)
			subj_obs = [subj_obs, unifrnd(subj_act_means(t)-item_spread, subj_act_means(t)+item_spread, 1, net_effects)];

			% We need to come up with unique net effect values for each net effect cluster. The below two lines are doing that
			tmp_net_effects_list = [tmp_net_effects_list, repmat(net_eff_counter, 1, net_effects)];
			net_eff_counter = net_eff_counter + 1;
		end

		% At this point we have all the activations for a subject i.e. samples x net effects. We scale the data within each subject to be between 0 and 1
		subj_obs = round_to(scale_data(subj_obs, 0, 1), 4);

		switch yval_distribution
		case 'bernoulli'
			% We use the scaled data to get the y vals from the ground truth curve
			% Why do we need y values? We need the y vals to assign a binary variable (dependent variable) to each item
			for u = unique(tmp_net_effects_list)
				idx = find(u == tmp_net_effects_list);
				yvals = [yvals, sum(out.yval(round_to((subj_obs(idx) + 1e-4) * 1e+4, 0)))];
			end
			% We assign 1's to y vals > median(y vals) and 0 otherwise. This ensures that 50% are assigned 1's and 50% assigned 0's
			% Note the reshape and repmat makes sure that the same dependent variable is assigned to all net effect vals
			dependent_var_list = [dependent_var_list, reshape(repmat((yvals > median(yvals)), net_effects, 1), 1, net_effects * nItems)];
		case 'normal'
			for u = unique(tmp_net_effects_list)
				idx = find(u == tmp_net_effects_list);
				%get y vals of each of our net effects and average that
				mean_yval = mean(out.yval(round_to((subj_obs(idx) + 1e-4) * 1e+4, 0)));
				yvals = [yvals, mean_yval];
			end
			
			dependent_var_list = [dependent_var_list, reshape(repmat(yvals,net_effects,1),1,net_effects * nItems)];

		otherwise
			error('Invalid yval_distribution! valid distributions are ''bernoulli'' and ''normal''');
		end
				% Gather the net effects in a vector
		net_effects_list = [net_effects_list, tmp_net_effects_list];

		% Create a raw_data matrix. The dimension of this matrix is different if we don't include net effects
		if s == nSubjects, raw_data = NaN(nSubjects * nItems * net_effects, 6); end
	else
		% Uniformly sample 'nItems' items per subject between 0 and 1
		subj_obs = round_to(unifrnd(0, 1, 1, nItems), 4);

		% We need to come up with unique net effect values for each net effect cluster. The below two lines are doing that
		net_effects_list = [net_effects_list, net_eff_counter:net_eff_counter + nItems - 1];
		net_eff_counter = net_eff_counter + nItems;

		% At this point we have all the activations for a subject i.e. samples x net effects. We scale the data within each subject to be between 0 and 1
		subj_obs = round_to(scale_data(subj_obs, 0, 1), 4);

		% We use the scaled data to get the y vals from the ground truth curve
		% Why do we need y values? We need the y vals to assign a binary variable (dependent variable) to each item
		yvals = out.yval(round_to((subj_obs + 1e-4) * 1e+4, 0));

		switch yval_distribution
		case 'bernoulli'
			% We assign 1's to y vals > median(y vals) and 0 otherwise. This ensures that 50% are assigned 1's and 50% assigned 0's
			dependent_var_list = [dependent_var_list, yvals > median(yvals)];
		case 'normal'
			% accumulate generated yvals, no need to add noise here
			dependent_var_list = [dependent_var_list, yvals];
		otherwise
			error('Invalid yval_distribution! valid distributions are ''bernoulli'' and ''normal''');
		end

		% Create a raw_data matrix. The dimension of this matrix is different if we have net effects
		if s == nSubjects, raw_data = NaN(nSubjects * nItems, 6); end
	end

	% Now that we have all the information we need we add truncated Gaussian noise to the predictor variable
	% Reason:- since we used the median to assign labels, a perfect (straight) logistic function that clssifies this dataset will have slope infinity i.e. slope = rise / run, here run = 0
	predictor_var_list = [predictor_var_list; truncated_normal(0, 1, subj_obs', noise_knob, length(subj_obs))];

	% We populate the subject Id and item list
	subj_id_list = [subj_id_list, repmat(s, 1, length(subj_obs))];
	item_list = [item_list, 1:length(subj_obs)];
end

% Save all the information into a struct
simulated_data = struct();
raw_data(:, 1) = subj_id_list;
raw_data(:, 2) = item_list;
raw_data(:, 3) = -1;
raw_data(:, 4) = predictor_var_list;
raw_data(:, 5) = dependent_var_list;
raw_data(:, 6) = net_effects_list;
simulated_data.raw_data = raw_data;
simulated_data.curve_xval = out.xval;
simulated_data.curve_yval = out.yval;
simulated_data.curve_params = out.curve_params;
simulated_data.nSubjects = nSubjects;
simulated_data.nItems = nItems;
simulated_data.net_effects = net_effects;
simulated_data.analysis_id = analysis_id;
simulated_data.noise_knob = noise_knob;
simulated_data.curve_type = curve_type;
simulated_data.yval_distribution = yval_distribution;

save(sprintf('%s/%s_simulated_data.mat', target_dir, analysis_id), 'simulated_data');

% Visualize the simulated data
plot_sim_data(analysis_id);

disp(sprintf('Simulated data and figures are stored in %s', target_dir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sim_data(analysis_id)

% [] = PLOT_SIM_DATA(ANALYSIS_ID)
%
% Purpose:
%
% Plot the simulated data
%
% Inputs:
%
% analysis_id: Valid analysis id
%
% Outputs:
%
% bunch of plots
%
% Example usage:
%
% plot_sim_data('my_analysis_id')
%

close all;

% Setting the target directory
results_dir = fullfile(pwd, 'results');
target_dir = fullfile(results_dir, analysis_id);
load(sprintf('%s/%s_simulated_data.mat', target_dir, analysis_id));
visible_off = true;

% I have hard coded each figure to have only 15 subjects i.e. 3 x 5 subplots
new_fig_prompt = 1:15:simulated_data.nSubjects;
save_fig_prompt = new_fig_prompt(2:end) - 1;

% For each subject we are going to create a sub plot and draw the ground truth curve, below the ground truth curve we are going to plot the predictor var for that subject
for s = 1:simulated_data.nSubjects
	if any(s == new_fig_prompt)
		if visible_off
			figure('visible', 'off');
		else
			figure();
		end
		set(gcf, 'Position', [50, 900, 1600, 1000]);
		fig_idx = 1;
	end

	subplot(3, 5, fig_idx);
	% Draw the ground truth curve
	plot(simulated_data.curve_xval, simulated_data.curve_yval); hold on;
	xlabel('Activation'); ylabel('Change in memory strength'); title(sprintf('subj %d', s));
	fig_idx = fig_idx + 1;

	% Recall the raw_data = [subject id, item, category, predictor var, dependent var, net effect cluster]
	% Here we are fetching the target subject idx
	subject_idx = find(simulated_data.raw_data(:, 1) == s);
	% Here we are fetching the associated net effect index for that subject
	net_effect_idx = unique(simulated_data.raw_data(subject_idx, 6));

	for i = 1:length(net_effect_idx)
		% Using the net effect indices we fish out each item, check if the associated dependent variable is 0, 1 OR normally distributed and then plot accordingly
		target_idx = intersect(find(simulated_data.raw_data(:, 6) == net_effect_idx(i)), subject_idx);

		switch simulated_data.yval_distribution
		case 'bernoulli'
			% Setting the y limits
			y_lim_lower = -3; y_lim_upper = 1;
			% This line of code is loop invariant but needs to be in here since the y limits are assigned only in the previous line
			y_coordinate = linspace(min(simulated_data.curve_yval)-0.1, y_lim_lower+0.1, length(net_effect_idx));

			if unique(simulated_data.raw_data(target_idx, 5)) > 0
				% If dependent variable is 1 plot that item in green
				plot(simulated_data.raw_data(target_idx, 4), repmat(y_coordinate(i), 1, length(target_idx)), 'go'); hold on;
			else
				% If dependent variable is 0 plot that item in red
				plot(simulated_data.raw_data(target_idx, 4), repmat(y_coordinate(i), 1, length(target_idx)), 'ro'); hold on;
			end
			if length(simulated_data.raw_data(target_idx, 4)) > 1, plot(mean(simulated_data.raw_data(target_idx, 4)), y_coordinate(i), 'k*'); end
		case 'normal'
			% Setting the y limits
			y_lim_lower = -1; y_lim_upper = 1;
			% If dependent variable is normally distributed then plot in magenta on top of the curve
			plot(simulated_data.raw_data(target_idx, 4), simulated_data.raw_data(target_idx, 5), 'm.');
		otherwise
			error('Invlaid yval_distribution! valid distributions are ''bernoulli'' and ''normal''');
		end
	end
	ylim([y_lim_lower, y_lim_upper]);

	if any(s == save_fig_prompt) | (s == simulated_data.nSubjects)
		savesamesize(gcf, 'file', sprintf('%s/%s_%d_simulated_data', target_dir, simulated_data.analysis_id, s), 'format', 'png');
	end
end

