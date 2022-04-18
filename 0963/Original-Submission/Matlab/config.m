%CONFIG Global configuration script
%
%    Creates the cgmm_config variable which is used to influence a lot of
%    different aspects of the program. Modify this file to change the
%    program's default behavior.
%
%
% created by Benedikt Rudolph
% DATE: 03-Dec-2012

% directories
cgmm_config.directories.project = pwd;
cgmm_config.directories.data = [ cgmm_config.directories.project '/data' ];
cgmm_config.directories.main = [ cgmm_config.directories.project '/main' ];
cgmm_config.directories.mat = [ cgmm_config.directories.project '/mat' ];
cgmm_config.directories.latex = [ cgmm_config.directories.project '/output/latex' ];
cgmm_config.directories.plot = [ cgmm_config.directories.project '/plot' ];
cgmm_config.directories.png = [ cgmm_config.directories.project '/output/png' ];
cgmm_config.directories.print = [ cgmm_config.directories.project '/print' ];

% time series to be used
cgmm_config.time_series.name = '_example'; % for storing the results
cgmm_config.time_series.file = [cgmm_config.directories.data '/example_series.csv']; % for loading the data

% Heuristic estimation configuration
cgmm_config.heuristic_estimation.simulation_runs = 1000;

% CGMM configuration
cgmm_config.cgmm.grid_min = -1000;
cgmm_config.cgmm.grid_max = 1000;
cgmm_config.cgmm.grid_res = 50;

% Monte Carlo study configuration
cgmm_config.monte_carlo.time_steps = [0.2 8 40]; % for fixed dt
cgmm_config.monte_carlo.simulation_runs = 100;

% Feature test configuration
cgmm_config.feature_test.M0 = 80;
cgmm_config.feature_test.N = 3000;

% extra options (the less you use, the better...)
cgmm_config.extra_options.wasc = '';
cgmm_config.extra_options.pcsv = '';

% files containing estimates
cgmm_config.estimates.wasc = [ cgmm_config.directories.mat '/estimates_wasc' ...
                               cgmm_config.time_series.name ...
                               cgmm_config.extra_options.wasc '.mat'];

cgmm_config.estimates.pcsv = [ cgmm_config.directories.mat '/estimates_pcsv' ...
                               cgmm_config.time_series.name ...
                               cgmm_config.extra_options.pcsv '.mat'];
                               
cgmm_config.estimates.pcsv1d = [ cgmm_config.directories.mat '/estimates_pcsv1d' ...
                               cgmm_config.time_series.name ...
                               cgmm_config.extra_options.pcsv '.mat'];

cgmm_config.estimates.pcsv_partial = [ cgmm_config.directories.mat '/estimates_pcsv_partial' ...
                               cgmm_config.time_series.name ...
                               cgmm_config.extra_options.pcsv '.mat'];

% files containing monte carlo study results
cgmm_config.monte_carlo.wasc = [ cgmm_config.directories.mat '/mc_wasc' ...
                                 cgmm_config.time_series.name ...
                                 cgmm_config.extra_options.wasc '.mat'];

cgmm_config.monte_carlo.pcsv = [ cgmm_config.directories.mat '/mc_pcsv' ...
                                 cgmm_config.time_series.name ...
                                 cgmm_config.extra_options.pcsv '.mat'];

cgmm_config.monte_carlo.pcsv1d = [ cgmm_config.directories.mat '/mc_pcsv1d' ...
                                 cgmm_config.time_series.name ...
                                 cgmm_config.extra_options.pcsv '.mat'];

cgmm_config.monte_carlo.pcsv_partial = [ cgmm_config.directories.mat '/mc_pcsv_partial' ...
                                 cgmm_config.time_series.name ...
                                 cgmm_config.extra_options.pcsv '.mat'];

cgmm_config.monte_carlo.gbm  = [ cgmm_config.directories.mat '/mc_gbm' ...
                                 cgmm_config.time_series.name ...
                                 cgmm_config.extra_options.pcsv '.mat'];

% files containing feature test results
                             
cgmm_config.feature_test_cmm = 1;
estimator = 'cmm';
if ~cgmm_config.feature_test_cmm
  estimator = 'cgmm';
end

% helper that contains model independent configuration
feature_test_file_config = [ cgmm_config.time_series.name ...
                             '_' estimator ...
                             '_M0-' num2str(cgmm_config.feature_test.M0) ...
                             '_N-' num2str(cgmm_config.feature_test.N) ];

cgmm_config.feature_test.wasc = [ cgmm_config.directories.mat '/featuretest_wasc' ...
                                  feature_test_file_config ...
                                  cgmm_config.extra_options.wasc '.mat'];
                                  
cgmm_config.feature_test.pcsv = [ cgmm_config.directories.mat '/featuretest_pcsv' ...
                                  feature_test_file_config ...
                                  cgmm_config.extra_options.pcsv '.mat'];

cgmm_config.feature_test.pcsv_partial = [ cgmm_config.directories.mat '/featuretest_pcsv_partial' ...
                                  feature_test_file_config ...
                                  cgmm_config.extra_options.pcsv '.mat'];

cgmm_config.feature_test.gbm = [ cgmm_config.directories.mat '/featuretest_gbm' ...
                                 feature_test_file_config '.mat'];

% files containing heuristic estimation monte carlo results
cgmm_config.heur_mc.wasc = [ cgmm_config.directories.mat '/heur_mc_wasc' ...
                                 cgmm_config.time_series.name ...
                                 cgmm_config.extra_options.wasc '.mat'];
                                 
cgmm_config.heur_mc.pcsv = [ cgmm_config.directories.mat '/heur_mc_pcsv' ...
                                 cgmm_config.time_series.name ...
                                 cgmm_config.extra_options.pcsv '.mat'];

% plot files
cgmm_config.plots.device = '-dpng';
cgmm_config.plots.mre.all = @(model) [ cgmm_config.directories.png '/mre_all_' model '.png' ];
cgmm_config.plots.mre.individual = @(model) [ cgmm_config.directories.png '/mre_individual_' model '.png' ];
cgmm_config.plots.feature_test = @(model) [ cgmm_config.directories.png '/feature_test_' estimator '_' model '.png' ];
cgmm_config.plots.simulation = @(model) [ cgmm_config.directories.png '/sim_' model '.png' ];
cgmm_config.plots.cf = @(model) [ cgmm_config.directories.png '/cf_' model '.png' ];
