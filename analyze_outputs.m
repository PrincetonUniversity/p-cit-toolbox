function[] = analyze_outputs(analysis_id, varargin)

% ANALYZE_OUTPUTS(...)
%
% Purpose:
%
% Analyse the expectation maximization plus importance sampling results
%
% Inputs:
%
% analysis_id: Valid analysis Id
% varargin: If 'bootstrap' then it also expects the total number of bootstrap runs; optional - original run analysis id
% 	    If 'scramble' then it also expects the total number of scramble runs; optional - original run analysis id
% 	-- Example 1: analyze_outputs('my_bootstrap_id', 'bootstrap', 100, 'my_original_id')
% 	-- Example 2: analyze_outputs('my_bootstrap_id', 'bootstrap', 100)
% 	-- Example 3: analyze_outputs('my_scramble_id', 'scramble', 100, 'my_original_id')
% 	-- Example 4: analyze_outputs('my_scramble_id', 'scramble', 100)
%
% Outputs:
%
% Bunch of plots
%
% Example usage:
%
% analyze_outputs('my_analysis_id')
% analyze_outputs('my_bootstrap_id', 'bootstrap', 100, 'my_original_id')
% analyze_outputs('my_bootstrap_id', 'bootstrap', 100)
% analyze_outputs('my_scramble_id', 'scramble', 100, 'my_original_id')
% analyze_outputs('my_scramble_id', 'scramble', 100)
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

close all;

% Checks if the correct number of arguments are passed in
if nargin < 1, error('Missing input parameters'); end

% Set paths
% Setting the target directory
root_dir = pwd;
results_dir = fullfile(root_dir, 'results');
read_dir = fullfile(results_dir, analysis_id);
write_dir = fullfile(results_dir, analysis_id);

resolution = 2;
credible_interval = 0.9; % credible interval a 0.9 implies that 90% of samples should lie in that interval
n_bins = 10; % number of bins to plot histgrams
image_format = 'png';
alpha = 0.05;
visible_off = true; % control whether or not the generated figures pop up on screen or not before being saved to disk.

% Check length of varargin to determine the type of analysis we're plotting, either one of the following: { regular, bootstrap, scrambled}
if length(varargin) > 1 & length(varargin) <= 3
	nRuns = varargin{2};
	x = 0:(1/10^resolution):1;
	y_all_runs = NaN(nRuns, length(x));
	weight_all_runs = NaN(nRuns, 1);
	switch varargin{1}
	case 'bootstrap'
		% Individually generate figures for each run, accumulating particle weights and y values for each run
		for b = 1:nRuns
			disp(sprintf('Bootstrap run %d', b));
			[weight_all_runs(b, :), y_all_runs(b, :)] = plot_figures(sprintf('%s_b%d', analysis_id, b), read_dir, write_dir, resolution,...
												credible_interval, n_bins, image_format, visible_off);
		end

		% Generate a bootstrap specific plot
		if (any(isnan(y_all_runs(:))) | any(isnan(weight_all_runs(:)))), error('Nan''s in vector!'); end
		plot_boot_results(analysis_id, write_dir, x, y_all_runs, weight_all_runs, alpha, resolution, credible_interval,...
												image_format, n_bins, varargin{:});
	case 'scramble'
		% Individually generate figures for each run, accumulating particle weights and y values for each run
		for s = 1:nRuns
			disp(sprintf('Scramble run %d', s));
			weight_all_runs(s, :) = plot_figures(sprintf('%s_s%d', analysis_id, s), read_dir, write_dir, resolution,...
												credible_interval, n_bins, image_format, visible_off);
		end
		% Generate a scramble specific plot
		plot_scram_results(analysis_id, write_dir, weight_all_runs, resolution, image_format, n_bins, varargin{:});

	otherwise, error('Invalid analysis type! Valid ones include ''bootstrap'' and ''scramble'' only');
	end
elseif length(varargin) == 0
	plot_figures(analysis_id, read_dir, write_dir, resolution, credible_interval, n_bins, image_format, visible_off);
else
	error('Invalid number of input arguments failed to trigger any of the analyses'' pipelines!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[p_wgts, y_final] = plot_figures(analysis_id, read_dir, write_dir, resolution, credible_interval, n_bins, image_format, visible_off)

% Load the target .mat file
importance_sampler_mat = load(sprintf('%s/%s_importance_sampler.mat', read_dir, analysis_id));

% Find the particles / curves indices that are theory consistent
th_con_idx = family_of_curves(importance_sampler_mat.analysis_settings.curve_type, 'count_particles', importance_sampler_mat.curve_params);

% Using the above indices fetch the sum of the particle weights. This is what we call / refer throughout as p(theory consistent)
% Here I am fetch the weights from the last em eiteration hence the 'end'
p_wgts = sum(importance_sampler_mat.normalized_weights(end, th_con_idx));


% Calculate the p-value that our beta1 is significantly different from 0
likelihood_ratio_test_for_beta1 = importance_sampler_mat.likratiotest;

% Calculate relative changes in fminunc over time
fminunc_fvals = importance_sampler_mat.exp_max_fval;
[max_fval, max_fval_iteration] = max(fminunc_fvals);
[min_fval, min_fval_iteration] = min(fminunc_fvals);
fminunc_val_change_dir = 'DEC';
fminunc_relative_value_name = 'max';
if max_fval_iteration > min_fval_iteration
	fminunc_val_change_dir = 'INC';
    fminunc_relative_value_name = 'min';
	fminunc_relative_change = ((max_fval - min_fval) / min_fval);
else
	fminunc_relative_change = ((max_fval - min_fval) / max_fval);
end



% Get the weighted sum of the curves over all particles and the associated credible interval. This is the blue line and the grey envelope around the blue line
weighted_curve_struct = common_to_all_curves(importance_sampler_mat.analysis_settings.curve_type, 'weighted_curve', importance_sampler_mat, credible_interval, resolution);
x_final = weighted_curve_struct.x_final;
y_final = weighted_curve_struct.y_final;
interval = weighted_curve_struct.interval;

% Plot the curve
if visible_off
	figure('visible', 'off');
else
	figure();
end
set(gcf, 'Position', [50, 900, 500, 500]);

color = [107, 107, 107] ./ 255; transparency = 0.4;
hhh = jbfill(x_final, interval(1, :), interval(2, :), color, color, 0, transparency); hold on;
plot(x_final, y_final, 'b-', 'LineWidth', 2);

ylabel('Change in Memory Strength', 'FontSize', 15, 'FontName', 'Helvetica');
ylim([-1, 1]);
xlabel('Activation', 'FontSize', 15, 'FontName', 'Helvetica');
xlim([0, 1]);

grid on; set(gca, 'Layer', 'top');
title(sprintf('P(theory consistent) = %0.4f', p_wgts));
file_name = sprintf('%s/%s_weighted_curve', write_dir, analysis_id);
print(gcf, sprintf('-d%s', image_format), '-painters', file_name);
disp(sprintf('1. Recovered curve plot is saved as %s.%s', file_name, image_format));

if visible_off
	figure('visible', 'off');
else
	figure();
end
set(gcf, 'Position', [50, 900, 1000, 1000]);

subplot(2, 2, 1);
beta_1 = importance_sampler_mat.hold_betas_per_iter(:, 2);
plot(0:importance_sampler_mat.analysis_settings.em_iterations, beta_1, 'bo-', 'MarkerFaceColor', 'b');

ylabel('beta 1');
ylim([-0.2, 2.2]);
xlabel('EM iterations');

title({'Beta 1 over em iterations' sprintf('P(Beta 1 = 0) = %0.6f',likelihood_ratio_test_for_beta1 )});
grid on; set(gca, 'Layer', 'top');

subplot(2, 2, 2);
hist(importance_sampler_mat.normalized_weights(end, :), n_bins); hold on;
h1 = plot(max(importance_sampler_mat.normalized_weights(end, :)), 1, 'ro', 'MarkerFaceColor', 'r');

ylabel('count', 'FontSize', 12, 'FontName', 'Helvetica');
xlabel('posterior weights' , 'FontSize', 12, 'FontName', 'Helvetica');

title({ sprintf('Distribution of posterior weights (%d iteration)', importance_sampler_mat.analysis_settings.em_iterations)  sprintf('sample max weight=%0.4f', max(importance_sampler_mat.normalized_weights(end, :))) });

grid on; set(gca, 'Layer', 'top');
legend([h1], 'Max weight'); 

subplot(2, 2, [3 4]);
fminunc_fvals = importance_sampler_mat.exp_max_fval;
plot(1:importance_sampler_mat.analysis_settings.em_iterations, fminunc_fvals, 'bo-', 'MarkerFaceColor', 'b');

ylabel('fminunc fval', 'FontSize', 12, 'FontName', 'Helvetica');
xlabel('EM iterations', 'FontSize', 12, 'FontName', 'Helvetica');

title({'fminunc fval over em iterations' sprintf('fminunc %s over iterations by %0.4f of %s fval',fminunc_val_change_dir,fminunc_relative_change,fminunc_relative_value_name)});
grid on; set(gca, 'Layer', 'top');

file_name = sprintf('%s/%s_report_plot', write_dir, analysis_id);
savesamesize(gcf, 'file', file_name, 'format', image_format);
disp(sprintf('2. Toolbox report plot is saved as %s.%s', file_name, image_format));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_boot_results(analysis_id, write_dir, x, y_all_runs, weight_all_runs, alpha, resolution, credible_interval, image_format, n_bins, varargin)

root_dir = pwd;
total_bootstrap_runs = varargin{2};

sample_idx(1) = floor(total_bootstrap_runs * (alpha / 2));
if sample_idx(1) == 0, sample_idx(1) = ceil(total_bootstrap_runs * (alpha / 2)); end
sample_idx(2) = ceil(total_bootstrap_runs * (alpha / 2));
sample_idx(3) = floor(total_bootstrap_runs * (1 - (alpha / 2)));
if sample_idx(3) == 0, ceil(total_bootstrap_runs * (1 - (alpha / 2))); end
sample_idx(4) = ceil(total_bootstrap_runs * (1 - (alpha / 2)));

sorted_y_all_runs = sort(y_all_runs);
envelope_bounds = NaN(2, length(x));
envelope_bounds(1, :) = (sorted_y_all_runs(sample_idx(1), :) + sorted_y_all_runs(sample_idx(2), :)) / 2;
envelope_bounds(2, :) = (sorted_y_all_runs(sample_idx(3), :) + sorted_y_all_runs(sample_idx(4), :)) / 2;

legend_str = '';
if length(varargin) == 3
	results_dir = fullfile(root_dir, 'results');
	read_dir = fullfile(results_dir, varargin{3});
	% Load the original .mat file
	importance_sampler_mat = load(sprintf('%s/%s_importance_sampler.mat', read_dir, varargin{3}));
	% Find the particles / curves indices that are theory consistent
	th_con_idx = family_of_curves(importance_sampler_mat.analysis_settings.curve_type, 'count_particles', importance_sampler_mat.curve_params);
	% Using the above indices fetch the sum of the particle weights. This is what we call / refer throughout as p(theory consistent)
	% Here I am fetch the weights from the last em eiteration hence the 'end'
	weight_original = sum(importance_sampler_mat.normalized_weights(end, th_con_idx));

	% Get the weighted sum of the curves over all particles and the associated credible interval. This is the blue line and the grey envelope around the blue line
	weighted_curve_struct = common_to_all_curves(importance_sampler_mat.analysis_settings.curve_type, 'weighted_curve', importance_sampler_mat, credible_interval, resolution);
	y_original = weighted_curve_struct.y_final;
	legend_str = 'Original recovered curve';
else
	weight_original = mean(weight_all_runs);
	y_original = mean(y_all_runs);
	legend_str = 'Mean(bootstrap runs)';
end
proportion = sum(weight_all_runs > 0.5) / length(weight_all_runs);
proportion_str = sprintf('proportion(runs)>0.5=%0.2f', proportion);

% Plot the curve
figure(); set(gcf, 'Position', [50, 900, 1200, 500]);

subplot(1, 2, 1);
color = [107, 107, 107] ./ 255; transparency = 0.4;
hhh = jbfill(x, envelope_bounds(1, :), envelope_bounds(2, :), color, color, 0, transparency); hold on;
% hAnnotation = get(hhh, 'Annotation');
% hLegendEntry = get(hAnnotation', 'LegendInformation');
% set(hLegendEntry, 'IconDisplayStyle', 'off');
h2 = plot(x, y_original, 'b-', 'LineWidth', 2);

ylabel('Change in Memory Strength', 'FontSize', 12, 'FontName', 'Helvetica');
ylim([-1, 1]);
xlabel('Activation', 'FontSize', 12, 'FontName', 'Helvetica');
xlim([0, 1]);

grid on; set(gca, 'Layer', 'top');
title(sprintf('%d bootstrap runs, recovered curves', total_bootstrap_runs));
legend([h2(1), hhh(1)], legend_str, sprintf('%%95 confidence interval'), 'Orientation', 'Horizontal');

sorted_weights = sort(weight_all_runs);
wgt_lb = (sorted_weights(sample_idx(1)) + sorted_weights(sample_idx(2))) / 2;
wgt_ub = (sorted_weights(sample_idx(3)) + sorted_weights(sample_idx(4))) / 2;

subplot(1, 2, 2);
numbers = hist(weight_all_runs, n_bins);
hist(weight_all_runs, n_bins); hold on;
h = findobj(gca, 'Type', 'patch'); set(h, 'FaceColor', [255,193,193] ./ 255, 'EdgeColor', 'w');
hAnnotation = get(h, 'Annotation');
hLegendEntry = get(hAnnotation', 'LegendInformation');
set(hLegendEntry, 'IconDisplayStyle', 'off');

h1 = plot(weight_original, 0:3, 'b*');
h2 = plot([repmat(wgt_lb, 1, 4), repmat(wgt_ub, 1, 4)], [0:3, 0:3], 'b*', 'Color', color, 'Marker', '*');

ylabel('Count', 'FontSize', 12, 'FontName', 'Helvetica');
ylim([0, max(numbers)+5]);
xlabel('Distribution of posterior weights', 'FontSize', 12, 'FontName', 'Helvetica');
xlim([-0.02, 1.02]);
title(sprintf('%d bootstrap runs, posterior weights\n%s', total_bootstrap_runs, proportion_str));
legend([h1(1), h2(1)], legend_str, sprintf('%%95 confidence interval'), 'Location', 'NorthWest', 'Orientation', 'Vertical');

file_name = sprintf('%s/%s_bootstrap_results', write_dir, analysis_id);
savesamesize(gcf, 'file', file_name, 'format', image_format);
disp(sprintf('Bootstrap results plot is saved as %s.%s', file_name, image_format));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_scram_results(analysis_id, write_dir, weight_all_runs, resolution, image_format, n_bins, varargin)

root_dir = pwd;
total_scramble_runs = varargin{2};

legend_str = '';
if length(varargin) == 3
	results_dir = fullfile(root_dir, 'results');
	read_dir = fullfile(results_dir, varargin{3});
	% Load the original .mat file
	importance_sampler_mat = load(sprintf('%s/%s_importance_sampler.mat', read_dir, varargin{3}));
	% Find the particles / curves indices that are theory consistent
	th_con_idx = family_of_curves(importance_sampler_mat.analysis_settings.curve_type, 'count_particles', importance_sampler_mat.curve_params);
	% Using the above indices fetch the sum of the particle weights. This is what we call / refer throughout as p(theory consistent)
	% Here I am fetch the weights from the last em eiteration hence the 'end'
	weight_original = sum(importance_sampler_mat.normalized_weights(end, th_con_idx));
	legend_str = 'Original recovered curve';
else
	weight_original = mean(weight_all_runs);
	legend_str = 'Mean(scramble runs)';
end
pval = sum(weight_all_runs >= weight_original) / (total_scramble_runs + 1);

% Plot the curve
figure(); set(gcf, 'Position', [50, 900, 600, 600]);

numbers = hist(weight_all_runs, n_bins);
hist(weight_all_runs, n_bins); hold on;
h = findobj(gca, 'Type', 'patch'); set(h, 'FaceColor', [255,193,193] ./ 255, 'EdgeColor', 'w');
hAnnotation = get(h, 'Annotation');
hLegendEntry = get(hAnnotation', 'LegendInformation');
set(hLegendEntry, 'IconDisplayStyle', 'off');

plot(weight_original, 0:3, 'b*');

ylabel('Count', 'FontSize', 15, 'FontName', 'Helvetica');
ylim([0, max(numbers)+5]);
xlabel('Distribution of posterior weights', 'FontSize', 15, 'FontName', 'Helvetica');
xlim([-0.02, 1.02]);
title(sprintf('%d scramble samples posterior weights\npval=%0.4f', total_scramble_runs, pval));
legend(legend_str, 'Orientation', 'Horizontal');

file_name = sprintf('%s/%s_scramble_results', write_dir, analysis_id);
savesamesize(gcf, 'file', file_name, 'format', image_format);
disp(sprintf('Scramble results plot is saved as %s.%s', file_name, image_format));
