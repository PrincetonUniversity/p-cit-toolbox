function [] = run_importance_sampler()

% [] = RUN_IMPORTANCE_SAMPLER()
% 
% Purpose
% 
% This scripts sets up the data matrix (number of samples x 6 columns) and the 'analysis_settings' struct with algorithm parameters
% 
% Input
%
% --None
% 
% Output
%
% --None
% 
% Example usage:
%
% run_importance_sampler()
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populating the analysis_settings struct with algorithm settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analysis_settings = struct(); % Creating a struct

analysis_settings.analysis_id = 'my_analysis_id'; % analysis_id: specifies the target directory into which the output .mat will be located; if empty then the target directory is the timestamp of the form - YYYY-MM-DD-HH-MM

analysis_settings.em_iterations = 20; % Number of expectation maximization iterations
analysis_settings.particles = 100000; % Number of particles to be used in the importance sampling algorithm
analysis_settings.curve_type = 'horz_indpnt'; % Name of the family of curves to be used. Refer to the family_of_curves.m file for more info

analysis_settings.distribution = 'bernoulli'; % Name of the distribution (and the default canonical link function which maps the predictor variable to the dependent variable)
analysis_settings.dist_specific_params = struct(); % For normal distribution the additional parameter is sigma. We pass in sigma here.
analysis_settings.dist_specific_params.sigma = 1;

analysis_settings.beta_0 = 0; % Initializing beta_0 for linear predictor
analysis_settings.beta_1 = 1; % Initializing beta_1 for linear predictor
analysis_settings.tau = 0.05; % Specifies the radius to sample curves in the curve space
analysis_settings.category = []; % Specifies if the analyses will need to run on a specific category. Vector length Should be greater than 0. For instance [2] will cause the analyses to be run only on the second category; [] will run the analyses on all categories

analysis_settings.drop_outliers = 3; % specifies how many std dev away from group mean will the predictor variable outliers need to be dropped
analysis_settings.zscore_within_subjects = false; % if TRUE, the independednt variables will be zscored within each suibject
% Registering which column in the data matrix is carrying which piece of information
analysis_settings.data_matrix_columns = struct();
analysis_settings.data_matrix_columns.subject_id = 1;
analysis_settings.data_matrix_columns.trials = 2;
analysis_settings.data_matrix_columns.category = 3;
analysis_settings.data_matrix_columns.predictor_var = 4;
analysis_settings.data_matrix_columns.dependent_var = 5;
analysis_settings.data_matrix_columns.net_effect_clusters = 6;

analysis_settings.resolution = 4; % Denotes the resolution in which the data will be processed
analysis_settings.particle_chunks = 2; % Denotes the number of chunks you plan to partition the trials x particles matrix. An example chunk size will be 2 for a 3000 x 50,000 matrix

analysis_settings.bootstrap = false; % indicates that this run is a bootstrap run
analysis_settings.bootstrap_run = -1; % will need to specify a bootstrap sample number. This will need to be unique for each sample

analysis_settings.scramble = false; % indicates that this run is a scramble run
analysis_settings.scramble_run = -1; % will need to specify a scramble sample number. This will need to be unique for each sample
analysis_settings.scramble_style = -1; % choosing the appropriate scramble option from three options below
if analysis_settings.scramble_style > 0
	switch analysis_settings.scramble_style
	case 1, analysis_settings.scramble_style = 'within_subjects_within_categories';
	case 2, analysis_settings.scramble_style = 'within_subjects_across_categories';
	case 3, analysis_settings.scramble_style = 'across_subjects_across_categories';
	otherwise, error('Invalid scramble style!');
	end
end

%%%%%%%%%%%%%%%%%%%%%
% Reading in the data
%%%%%%%%%%%%%%%%%%%%%
% The three lines below load the simulated data into the raw_data matrix. Replace these two lines of the code with code to load your actual data

results_dir = fullfile(pwd, 'results');
load(sprintf('%s/%s/%s_simulated_data.mat', results_dir, analysis_settings.analysis_id, analysis_settings.analysis_id));
raw_data = simulated_data.raw_data;

importance_sampler(raw_data, analysis_settings);

