% render_into_file(basename, plot_params)
%
% a utility for saving plots to a file from matlab.
% intended for use with the paramotopy suite of functions.
%
% input:
%  plotmode -- a string which forms the basis of a filename

%  plot_params -- a structure with data members for setting up the save.
%      members:  format, a string for the format to save.  see the print
%                         help.
%                format_flag, again a string for the save driver to use.
%
%  outputs a file of the appropriate format to the current working dir.
%
%
%

% daniel brake
% colorado state, north carolina state, notre dame
% mathematics
% 2013-15
% danielthebrake@gmail.com


function plot_params = render_into_file(varargin)
%

[~, default_filename, ~] = fileparts(pwd);
default_filename = default_filename(~isspace(default_filename));
			
display(varargin)
if isempty(varargin)
	
	plot_params.fontsize = 16;
	plot_params.format = 'eps';
	plot_params.format_flag = 'psc2';
	plot_params.basename = default_filename;
elseif and(ischar(varargin{1} ), length(varargin)==1)
	
	plot_params.format = 'eps';
	plot_params.format_flag = 'psc2';
	
	plot_params.basename = varargin{1};
	
elseif isstruct(varargin{1})
	plot_params = varargin{1};
	
	if ~isfield(plot_params,'basename')
		error('incomplete plot_params.  add field ''basename''');
	end
else
	plot_params.window = get(gcf,'Position');
	plot_params.format = 'eps';
	plot_params.format_flag = 'psc2';
	plot_params.basename = default_filename;
end

fig1 = gcf;

set(fig1,'PaperPositionMode','auto');

if isfield(fig1,'PaperSize')
	fig1.PaperSize = fig1.PaperPosition(3:4);
end



currname = increment_name(plot_params.basename);
nameforfile = sprintf('%s.%s',currname,plot_params.format);

varargin
if strcmp(plot_params.format,'eps')
	varargin = [varargin, {'-fillpage'}]
end

evalme = sprintf('print(fig1,''%s'',''-d%s''',nameforfile,plot_params.format_flag);
for ii = 2:length(varargin)
    evalme = sprintf('%s,''%s''',evalme,varargin{ii});
end
evalme = sprintf('%s);',evalme);
display(evalme);
eval(evalme);
end

