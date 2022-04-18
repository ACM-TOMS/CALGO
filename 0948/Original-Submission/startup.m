function startup
% function DAESAstartup
% DAESA startup file
% Set the directories below to Matlab's search path

global DAESA_HOME;
DAESA_HOME=pwd;
addpath(DAESA_HOME,'-begin');
addpath([DAESA_HOME,'/src'],'-begin');
addpath([DAESA_HOME,'/src/implementation'],'-begin');
addpath([DAESA_HOME,'/src/implementation/common'],'-begin');
addpath([DAESA_HOME,'/src/implementation/lapdm'],'-begin');
addpath([DAESA_HOME,'/src/implementation/visualization'],'-begin');
addpath([DAESA_HOME,'/src/implementation/dataoutput'],'-begin');
addpath([DAESA_HOME,'/src/interface'],'-begin');
addpath([DAESA_HOME,'/examples'],'-begin');
addpath([DAESA_HOME,'/examples/DAEs'],'-begin');


end
