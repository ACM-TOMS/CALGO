function gpopsIntlabSetup(intlab_home_directory);

% -----------------------------------------------
% This function sets up INTLAB for Use with GPOPS
% -----------------------------------------------
% 
% Input: intlab_home_directory=Base Directory of INTLAB
% (e.g. Windows: c:\Intlab;  Unix: /home/user/Intlab)
					 dirs = {'intval','gradient','hessian','slope','utility','long','AccSum'};
if ~isequal(intlab_home_directory(end),'/'),
  intlab_home_directory = strcat(intlab_home_directory,'/');
end;
for i=1:length(dirs)
  currdir = strcat(intlab_home_directory,dirs{i});
  addpath(currdir);
end;
savepath;
