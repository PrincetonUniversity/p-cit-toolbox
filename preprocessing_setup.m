function[data, analysis_settings] = preprocessing_setup(data, analysis_settings)

% [PREPROCESSED_DATA, UPDATED_ANALYSIS_SETTINGS] = PREPROCESSING_SETUP(RAW_DATA, ANALYSIS_SETTINGS)
% 
% Purpose
% 
% To peform sanity checks on the input data and the algorithm parameter struct. Massage the data (i.e. drop outliers, zscore data, etc)
% 
% Input
%
% --data: Input data matrix (total number of trials x 6 columns)
% --analysis_settings: Struct with algorithm parameters
%
% Output
%
% --data: Input data matrix (if applicable, outlier free, zscored, category specific data only, etc)
% --analysis_settings: Struct with algorithm parameters; some additional parameters are added to this struct as well
% 
% Example usage:
%
% preprocessing_setup(data_matrix, analysis_settings_struct)
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

%%%%%%%%%%%%%%%
% Sanity checks
%%%%%%%%%%%%%%%

disp(sprintf('********** START OF MESSAGES **********'));

% Checks if the correct number of arguments are passed in
if nargin < 2, error('Missing input parameters'); end

% Checks if the data matrix has 6 columns
nCols = size(data, 2);
if nCols ~= 6
	error('Incorrect number of columns (%d) in the input matrix!', nCols);
end

% Registering which column in the data matrix is carrying which piece of information
if ~isfield(analysis_settings, 'data_matrix_columns') | isempty(analysis_settings.data_matrix_columns)
	% Setting it to the default
	analysis_settings.data_matrix_columns = struct();
	analysis_settings.data_matrix_columns.subject_id = 1;
	analysis_settings.data_matrix_columns.trials = 2;
	analysis_settings.data_matrix_columns.category = 3;
	analysis_settings.data_matrix_columns.predictor_var = 4;
	analysis_settings.data_matrix_columns.dependent_var = 5;
	analysis_settings.data_matrix_columns.net_effect_clusters = 6;
end
subject_id_column = analysis_settings.data_matrix_columns.subject_id;
trials_column = analysis_settings.data_matrix_columns.trials;
category_column = analysis_settings.data_matrix_columns.category;
predictor_var_column = analysis_settings.data_matrix_columns.predictor_var;
dependent_var_column = analysis_settings.data_matrix_columns.dependent_var;
net_effect_clusters_column = analysis_settings.data_matrix_columns.net_effect_clusters;

% Checks if the em iterations is specified; if not specified then it is set to a default of 20
if ~isfield(analysis_settings, 'em_iterations') | analysis_settings.em_iterations <= 0
	analysis_settings.em_iterations = 20;
	disp(sprintf('Missing number of iterations! It is set to a default of %d', analysis_settings.em_iterations));
end

% Checks if the no. of particles is specified; if not specified then it is set to a default of 1000
if ~isfield(analysis_settings, 'particles') | analysis_settings.particles <= 0
	analysis_settings.particles = 100000;
	disp(sprintf('Missing number of particles! It is set to a default of %d', analysis_settings.particles));
end

% Checks if the family of curves is specified; if not specified then it is set to a default of 'horz_indpnt' (Refer to family of curves)
if ~isfield(analysis_settings, 'curve_type') | isempty(analysis_settings.curve_type)
	analysis_settings.curve_type = 'horz_indpnt';
	disp(sprintf('Missing family of curves! It is set to a default of %s', analysis_settings.curve_type));
end

% Checks if the family of curves exist by fetching the number of curve parameters. This is just a sanity check
if ~isnumeric(family_of_curves(analysis_settings.curve_type, 'get_nParams'))
	error('%s - Does not exist! Check family_of_curves.m script', analysis_settings.curve_type);
end

% Checks if the distribution is specified;
% If not specified and if the dependent variable is binary then it is set to 'bernoulli'; otherwise it is set to to 'normal'
if ~isfield(analysis_settings, 'distribution') | isempty(analysis_settings.distribution)
	if length(unique(data(:, dependent_var_column))) == 2
		analysis_settings.distribution = 'bernoulli';
	else
		analysis_settings.distribution = 'normal';
	end
	disp(sprintf('Missing distribution! based on the dependent variable it is set to %s', analysis_settings.distribution));
end

% Checks if the distribution specific parameters exist
if ~isfield(analysis_settings, 'dist_specific_params') | isempty(analysis_settings.dist_specific_params)
	switch analysis_settings.distribution
	case 'bernoulli'
		% For a Bernoulli dist there are no parameters so it is empty. We still need the struct to exist
		analysis_settings.dist_specific_params = struct();
	case 'normal'
		% For normal distribution the additional parameter is sigma. We pass in sigma here.
		analysis_settings.dist_specific_params = struct();
		analysis_settings.dist_specific_params.sigma = 1; % Default is 1
		disp(sprintf('Missing sigma for normal distribution! It is set to %d', analysis_settings.dist_specific_params.sigma));
	otherwise
		disp(sprintf('Do you need additional parameters for your %s distribution? You can add them here. If none then create an empty struct like case ''bernoulli''', analysis_settings.distribution));
		keyboard;
	end
end

% Checks if normal distribution specific parameter is valid i.e. sigma > 0
if strcmp(analysis_settings.distribution, 'normal') & analysis_settings.dist_specific_params.sigma <= 0
	error('Normal distribution sigma will need to > 0! sigma = %0.4f', analysis_settings.dist_specific_params.sigma);
end

% Checks if beta_0 is specified; if not specified then it is set to a default of 0
if ~isfield(analysis_settings, 'beta_0')
	analysis_settings.beta_0 = 0;
	disp(sprintf('Missing initial setting for beta_0! It is set to a default of %d', analysis_settings.beta_0));
end

% Checks if beta_1 is specified; if not specified then it is set to a default of 1
if ~isfield(analysis_settings, 'beta_1')
	analysis_settings.beta_1 = 1;
	disp(sprintf('Missing initial setting for beta_1! It is set to a default of %d', analysis_settings.beta_1));
end

% Checks if tau is specified; if not specified then it is set to a default of 0.05
if ~isfield(analysis_settings, 'tau') | analysis_settings.tau <= 0
	analysis_settings.tau = 0.05;
	disp(sprintf('Missing tau! It is set to a default of %0.4f', analysis_settings.tau));
end

% Checks if this is a bootstrap run; if field does not exist then it is set to false
if ~isfield(analysis_settings, 'bootstrap')
	analysis_settings.bootstrap = false;
end

% Checks if bootstrap flag is boolean
if ~islogical(analysis_settings.bootstrap)
	error('analysis_settings.bootstrap field will need to be boolean!');
end

% Checks if this is a scramble run; if field does not exist then it is set to false
if ~isfield(analysis_settings, 'scramble')
	analysis_settings.scramble = false;
end

% Checks if scramble flag is boolean
if ~islogical(analysis_settings.scramble)
	error('analysis_settings.scramble field will need to be boolean!');
end

% Errors if both bootstrap and scramble flags exist
if analysis_settings.bootstrap & analysis_settings.scramble
	error('Cannot run both scramble AND bootstrap analyses at the same time! Set any one flag to be false');
end

% Builds a bootstrap data matrix from the original data matrix
if analysis_settings.bootstrap & ~analysis_settings.scramble
	% We need a bootstrap sample number
	if ~isfield(analysis_settings, 'bootstrap_run') | isempty(analysis_settings.bootstrap_run)
		error('Missing bootstrap sample number! set analysis_settings.bootstrap_run to a valid sample number');
	end

	bootstrap_data = [];
	new_cluster_count = 1;
	new_subject_count = 1;

	% Get the number of subjects from the data matrix
	nSubjects = length(unique(data(:, subject_id_column)));
	% Randomly sample with replacement the number of subjects thus generating our bootstrap sample
	subj_num_with_replacement = randsample([1:nSubjects], nSubjects, true);
	% For each subject in our bootstrap sample gather all relevant information
	for i = 1:length(subj_num_with_replacement)
		subj_idx = find(data(:, subject_id_column) == subj_num_with_replacement(i));

		% Recreate a new net effect cluster since this will need to be unique in the data matrix
		% (by repeatedly sampling subjects we could be repeating the net effect clusters)
		cluster_vector = data(subj_idx, net_effect_clusters_column);
		cluster_numbers = unique(cluster_vector);
		for j = 1:length(cluster_numbers)
			target_idx = find(data(subj_idx, net_effect_clusters_column) == cluster_numbers(j));
			cluster_vector(target_idx) = new_cluster_count;
			new_cluster_count = new_cluster_count + 1;
		end

		% Recreate a new subject id
		% (by repeatedly sampling subjects we could be repeating the subject id's)
		% Gather all information into a bootstrap_data matrix
		bootstrap_data = [bootstrap_data; [repmat(new_subject_count, length(subj_idx), 1), data(subj_idx, trials_column:dependent_var_column), cluster_vector]];
		new_subject_count = new_subject_count + 1;
	end

	% Perform some sanity checks to ensure that the bootstrap_data matrix is similar to the actual data matrix
	if ~isequal(size(bootstrap_data), size(data))
		error('Size of bootstrap dataset NOT the same as original data!');
	end
	if ~isequal(length(unique(data(:, net_effect_clusters_column))), length(unique(bootstrap_data(:, net_effect_clusters_column))))
		error('The number of clusters are not the same in the original and bootstrap sample!');
	end
	if ~isequal(data(:, subject_id_column), bootstrap_data(:, subject_id_column))
		error('The ordering of subjects are not the same in the original and bootstrap sample!');
	end

	% Store away the bootstrap sample subject information for future reference
	analysis_settings.bootstrap_run_subj_id = subj_num_with_replacement;
	data = bootstrap_data;
end

% Checks if this analysis will be need to performed for a specific category; if not specified then it is set to a default of [] i.e. NOT category specific
if ~isfield(analysis_settings, 'category')
	analysis_settings.category = [];
	disp(sprintf('Missing category specific analyses information! We are going to ignore the category dimension i.e. all trials from all categories will be analysed'));
end

% If this analysis is to be performed for a specific category then it filters out data from other irrelavant categories
if length(analysis_settings.category) > 0
	target_cat_idx = [];
	data_cat = unique(data(:, category_column));
	for c = 1:length(analysis_settings.category)
		cat_exist = find(data_cat == analysis_settings.category(c));
		if isempty(cat_exist)
			error('Category does not exist! You have set analysis_settings.category[%d]=%d', c, analysis_settings.category(c));
		end
		target_cat_idx = [target_cat_idx, find(data(:, category_column) == analysis_settings.category(c))];
	end
	data = data(target_cat_idx, :);
end

% Checks if outliers (i.e. data trials) will need to dropped; if not specified then it is set to a default of 'DO NOT DROP OUTLIERS'
if ~isfield(analysis_settings, 'drop_outliers')
	analysis_settings.drop_outliers = 3;
	disp(sprintf('Missing drop_outliers specific information! We are dropping outliers that are %d standard deviations away from the group mean', analysis_settings.drop_outliers));
end

% If this analysis requires the outliers to be dropped, then the code below drops the data trials within XXXX standard deviations from the GROUP MEAN
if analysis_settings.drop_outliers > 0
	% NaN's do not qualify as outliers so we filter them out and add them at the end of this step
	nan_free_idx = ~isnan(data(:, predictor_var_column));
	% NaN free data
	nan_free_data = data(nan_free_idx, :);
	std_dev_predictor_var = std(nan_free_data(:, predictor_var_column)) * analysis_settings.drop_outliers;
	mean_predictor_var = mean(nan_free_data(:, predictor_var_column));
	predictor_var_idx = (nan_free_data(:, predictor_var_column) > (mean_predictor_var - std_dev_predictor_var)) &...
			      (nan_free_data(:, predictor_var_column) < (mean_predictor_var + std_dev_predictor_var));
	disp(sprintf('%d trials are dropped since they are regarded as outliers', size(nan_free_data, subject_id_column) - sum(predictor_var_idx)));
	nan_free_data_outlier_dropped = nan_free_data(predictor_var_idx, :);
	% NaN's trials
	nan_data = data(~nan_free_idx, :);
	% Combine the NaN data with the outlier free data
	data = [nan_free_data_outlier_dropped; nan_data];
end

% Following the 'filter by category' and 'drop outliers', if applicable, we check if the data matrix is empty
nTrials = size(data, subject_id_column);
if nTrials <= 0, error('No input data!'); end

% Checks if we need to zscore predictor var within subjects, if not specified then it is set to default of FALSE
if ~isfield(analysis_settings, 'zscore_within_subjects')
	analysis_settings.zscore_within_subjects = 0;
	disp(sprintf('Missing zscore_within_subjects information! We are NOT zscoring within subjects'));
end

% Verifies if zscore within subjects is boolean
if ~islogical(analysis_settings.zscore_within_subjects)
	error('zscore_within_subjects field will need to be boolean!');
end

% Zscore the predictor variable within each subject
if analysis_settings.zscore_within_subjects
	% NaN's do not qualify to be zscored
	nan_free_idx = ~isnan(data(:, predictor_var_column));
	% NaN free data
	nan_free_data = data(nan_free_idx, :);
	% We get the list of subject id's (we use this cell array in zscoring the data within each subject, if applicable)
	subject_id_list = unique(nan_free_data(:, subject_id_column));
	% We get the number of subjects
	nSubjects = length(subject_id_list);
	if nSubjects <= 0, error('Not valid number of subjects!'); end
	for s = 1:nSubjects
		subject_idx = find(nan_free_data(:, subject_id_column) == subject_id_list(s));
		nan_free_data(subject_idx, predictor_var_column) = zscore(nan_free_data(subject_idx, predictor_var_column));
	end
	disp(sprintf('Predictor variables within each subject are zscored!'));
	% NaN's trials
	nan_data = data(~nan_free_idx, :);
	% Combine the NaN data with the outlier free data
	data = [nan_free_data; nan_data];
end

% Checks if the resolution is specified, if not specified then it is set to default of 4. This translates to 1e-4 = 0.0001
if ~isfield(analysis_settings, 'resolution') | analysis_settings.resolution <= 0
	analysis_settings.resolution = 4;
	disp(sprintf('Missing resolution! It is set to a default of %d', analysis_settings.resolution));
end

% if we have normally distributed data, we want to z-score the dependent variable
if strcmp(analysis_settings.distribution, 'normal')
	data(:,dependent_var_column) = zscore(data(:, dependent_var_column));
end

% We scale the predictor var to be between 0 and 1 and round it to 4 digits
nan_free_idx = ~isnan(data(:, predictor_var_column));
nan_free_data = data(nan_free_idx, :);
nan_free_data(:, predictor_var_column) = round_to(scale_data(nan_free_data(:, predictor_var_column), 0, 1), analysis_settings.resolution);
nan_data = data(~nan_free_idx, :);
data = [nan_free_data; nan_data];

% Scrambling the data matrix
if analysis_settings.scramble & ~analysis_settings.bootstrap
	if ~isfield(analysis_settings, 'scramble_run') | isempty(analysis_settings.scramble_run)
		error('Missing scramble sample number! set analysis_settings.scramble_run to a valid sample number');
	end
	if ~isfield(analysis_settings, 'scramble_style') | isempty(analysis_settings.scramble_style)
		analysis_settings.scramble_style = 'within_subjects_within_categories'; % The most conservative among all scramble techniques
		disp(sprintf('Missing scramble style! It is set a default of %s', analysis_settings.scramble_style));
	end

	% We get the list of subject id's
	subject_id_list = unique(data(:, subject_id_column));
	% We get the number of subjects in this analysis
	nSubjects = length(subject_id_list);
	if nSubjects <= 0, error('Not valid number of subjects!'); end

	switch analysis_settings.scramble_style
	case 'within_subjects_within_categories'
		% Here we scramble all dependent variables WHILE respecting the net effect boundaries, subject groupings and category groupings
		categories = unique(data(:, category_column));
		for s = 1:nSubjects
			for c = 1:length(categories)
				subject_category_idx = find(data(:, subject_id_column) == subject_id_list(s) & data(:, category_column) == categories(c));
				if length(subject_category_idx) > 1
					data(subject_category_idx, dependent_var_column) = scramble_dependent_variable(...
							data(subject_category_idx, dependent_var_column), data(subject_category_idx, net_effect_clusters_column));
				end
			end
		end

	case 'within_subjects_across_categories'
		% Here we scramble all dependent variables WHILE respecting the net effect boundaries and subject groupings
		for s = 1:nSubjects
			subject_idx = find(data(:, subject_id_column) == subject_id_list(s));
			if length(subject_idx) > 1
				data(subject_idx, dependent_var_column) = scramble_dependent_variable(...
							data(subject_idx, dependent_var_column), data(subject_idx, net_effect_clusters_column));
			end
		end

	case 'across_subjects_across_categories'
		% Here we scramble all dependent variables WHILE respecting the net effect boundaries
		all_idx = 1:size(data, 1);
		if length(all_idx) > 1
			data(all_idx, dependent_var_column) = scramble_dependent_variable(...
							data(all_idx, dependent_var_column), data(all_idx, net_effect_clusters_column));
		end

	otherwise error('Invalid analysis_settings.scramble_style=%s', analysis_settings.scramble_style)
	end
end

% Our data matrix looks like data = [subject id, item, category, predictor var, dependent var, net effect cluster]
% We verify if the subject id and dependent var columns are unique for the net effect clusters
% Below is a example of a valid data matrix (note dependent variable is unique within net effect cluster 111)
% data(1, :) = [24, 1, 1, 0.3333, 0, 111]
% data(2, :) = [24, 2, 2, 0.2222, 0, 111]
% data(3, :) = [24, 3, 1, 0.4444, 0, 111]
% Below is a example of an invalid data matrix (note dependent variable is not unique within net effect cluster 111)
% data(1, :) = [24, 1, 1, 0.3333, 0, 111]
% data(2, :) = [24, 2, 2, 0.2222, 1, 111]
% data(3, :) = [24, 3, 1, 0.4444, 0, 111]

% Fetching the net effect clusters
net_effect_clusters = unique(data(:, net_effect_clusters_column));
analysis_settings.net_effect_clusters = net_effect_clusters;

% If net effect clusters exist verify if the Subject Id and dependent variable are unique for those clusters
if length(net_effect_clusters) ~= size(data, 1)
	for i = 1:length(net_effect_clusters)
		cluster_idx = find(data(:, net_effect_clusters_column) == net_effect_clusters(i));
		if size(unique(data(cluster_idx, [subject_id_column, dependent_var_column]), 'rows'), 1) ~= 1
			error('Subject Id and/or dependent variable not unique for net effect cluster %d! Check the data matrix', net_effect_clusters(i));
		end
	end
else
	% If net effect clusters DO NOT exist then we treat each row as a net effect cluster by itself
	disp(sprintf('Each row will be treated separately. We will NOT be computing the net effect of any rows'));
end

% We create an analysis id unique to this analysis
if ~isfield(analysis_settings, 'analysis_id') | isempty(analysis_settings.analysis_id)
	time = clock();
	analysis_settings.analysis_id = sprintf('%d-%d-%d-%d-%d', time(1), time(2), time(3), time(4), time(5));
end

% We create a results directory if no specific target directory is mentioned
if ~isfield(analysis_settings, 'target_dir') | isempty(analysis_settings.target_dir)
	results_dir = fullfile(pwd, 'results');
	if ~exist(results_dir), mkdir(results_dir); end
	analysis_settings.target_dir = results_dir;
end
% target_directory = 'results/analysis_id'
analysis_settings.target_dir = fullfile(analysis_settings.target_dir, analysis_settings.analysis_id);
if ~exist(analysis_settings.target_dir), mkdir(analysis_settings.target_dir); end

% Due to memory constraints we perform two chunking tricks

% Chunking trick I
% In the curve fitting algorithm we need to compute the p(current iteration curves | previous iteration curves). This matrix is huge when the number of particles (curves) is 
% large, say 100,000. Even with a 8 Gb RAM, dedicated to Matlab, we still get a out of memory errors. To avoid this problem we chunk the matrix into smaller, more manageable matrices.
% Setting the chunk size to be particles x 0.05 -> 100,000 x 0.05 = 5000, translstes to p(current iteration curves(5000 curves at a time) | previous iteration curves).
analysis_settings.wgt_chunks = analysis_settings.particles * 0.05;
% If the chunk size is less then 5000 we set it be the number of particles itself
if analysis_settings.wgt_chunks < 5000
	analysis_settings.wgt_chunks = analysis_settings.particles;
end

% Chunking trick II
if ~isfield(analysis_settings, 'particle_chunks')
	analysis_settings.particle_chunks = 2;
	disp(sprintf('Missing particle chunks! It is set to a default of %d', analysis_settings.particle_chunks));
end
% Depending on the number of particle chunks we get start, end points and the number of particles within each chunk. For instance 1000 particles divided into 4 chunks will look like,
% 1	250	250
% 251	500	250
% 501	750	250
% 751	1000	250
analysis_settings.ptl_chunk_idx(:, 1) = 1:(analysis_settings.particles / analysis_settings.particle_chunks):analysis_settings.particles;
analysis_settings.ptl_chunk_idx(:, 2) = [analysis_settings.ptl_chunk_idx(2:end, 1) - 1; analysis_settings.particles];
analysis_settings.ptl_chunk_idx(:, 3) = analysis_settings.ptl_chunk_idx(:, 2) - analysis_settings.ptl_chunk_idx(:, 1) + 1;

% Storing analysis relevant information into the analysis_settings struct
% We get the list of subject id's
subject_id_list = unique(data(:, subject_id_column));
% We get the number of subjects in this analysis
analysis_settings.nSubjs = length(subject_id_list);
if analysis_settings.nSubjs <= 0, error('Not valid number of subjects!'); end

disp(sprintf('********** END OF MESSAGES **********'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[scrambled_vector] = scramble_dependent_variable(target_dependent_variables, target_net_effect_clusters)

% [SCRAMBLED_VECTOR] = SCRAMBLE_DEPENDENT_VARIABLE(TARGET_DEPENDENT_VARIABLES, TARGET_NET_EFFECT_CLUSTERS)
% 
% Purpose
% 
% To take a dependent variable vector and scramble it such that the net effect cluster groupings are NOT violated
% 
% Input
%
% --target_dependent_variables: The vector you would like scrambled
% --target_net_effect_clusters: The groupings that you would like to NOT violate. Follow the example below
%
% Output
%
% --scrambled_vector: A scrambled vector
% 
% Example usage:
%
% scramble_dependent_variable(target_dependent_variables, target_net_effect_clusters)
%

if ~isequal(size(target_dependent_variables), size(target_net_effect_clusters))
	error('Size of input vectors must be the same!');
end

% Detailed example
% example data matrix: target_dependent_variables = [1, 0, 1, 0, 0, 0, 1] and target_net_effect_clusters = [3, 5, 3, 7, 7, 5, 8]

% Fetch the sorted list of net effect clusters and their respective locations
% e.g. for [3, 5, 3, 7, 7, 5, 8] will return [3, 3, 5, 5, 7, 7, 8] and [1, 3, 2, 6, 4, 5, 7]
[sorted_neteff_clusters, sorted_neteff_clusters_indices] = sort(target_net_effect_clusters);
just_ones = ones(size(sorted_neteff_clusters)); % Populate a vector full of ones
% compute the length of each net effect cluster
% e.g. for [3, 5, 3, 7, 7, 5, 8] will return [2, 2, 2, 1] i.e. 3 is repeated twice and so on
length_each_neteff_cluster = accumarray(sorted_neteff_clusters(:), just_ones(:), [], @sum)';
length_each_neteff_cluster = length_each_neteff_cluster(logical(length_each_neteff_cluster));

% Get the unique list of clusters (i.e. excluding repetitions if any) e.g. [3, 5, 7, 8]
[unique_neteff_clusters, unique_indices] = unique(target_net_effect_clusters);
% Get the associated dependent variables (one per cluster; recall it is unique within a cluster) e.g. [1, 0, 0, 1]
associated_dependent_variables = target_dependent_variables(unique_indices);
% scramble the dependent variables e.g. [0, 0, 1, 1]
scrambled_indices = randperm(length(associated_dependent_variables));
scrambled_dependent_variables = associated_dependent_variables(scrambled_indices);

% Now we will need to repeat each scrambled dependent variable for the length of that net effect cluster. The next three lines will result in
% [0, 0, 0, 0, 1, 1, 1] corresponding to [3, 3, 5, 5, 7, 7, 8] since the scrambled dependent variable looks like [0, 0, 1, 1] for [3, 5, 7, 8]
cumsum_clutsers = cumsum(length_each_neteff_cluster);
indicator_vector = zeros(1, cumsum_clutsers(end));
indicator_vector([1, cumsum_clutsers(1:end-1)+1]) = 1;

% Store the scrambled dependent variable in the respective cluster locations
% The original vector looked like [3, 5, 3, 7, 7, 5, 8] so the scrambled vector will look like [0, 0, 0, 1, 1, 0, 1]
scrambled_vector = NaN(size(sorted_neteff_clusters_indices));
scrambled_vector(sorted_neteff_clusters_indices) = scrambled_dependent_variables(cumsum(indicator_vector));

if any(isnan(scrambled_vector(:))), error('Nan''s in scrambled dependent variable vector!'); end

