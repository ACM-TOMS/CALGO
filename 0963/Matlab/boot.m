%BOOT Startup script for all applications of CGMM estimation
%
%    - Prints a copyright notice.
%    - Adds all subfolders to the current path so that function files
%      can be called from scripts and functions in different folders.
%
% See also PATH.
%
% created by Benedikt Rudolph
% DATE: 16-Aug-2012

disp(repmat('-',1,80));
disp('Starting up CGMM estimation...');
disp(strcat('This program has been developed in a research project' ...
    , ' at Risklab, Toronto.'));
disp('(c) Benedikt Rudolph, 2013');
disp(repmat('-',1,80));
disp('');

% Add all subfolders to the current path
addpath(genpath(pwd));

% OCTAVE specific requirements
if exist('OCTAVE_VERSION')
  % load package financial for cov2corr
  pkg load financial;
end

% create configuration
config;
