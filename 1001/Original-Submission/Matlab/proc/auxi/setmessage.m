%% setmessage
% Display a message of parameter.
%
%% Syntax
%
%   setmessage(struct,fieldname)
%   setmessage(struct,fieldname,dispDepth)
%
%% Description
%
% Display the message of the set value of a parameter.
% This is used when the parameter was automatically set to a default value. 
%
%% Example
%
%   setmessage(seti,'testdirname');
%
% Then Output is:
%
% |Parameter "testdirname" was not set (correctly). Set testdirname = output/20160816T102134.|
%
%% Input Arguments
%
% * |struct|      : structure array with data
% * |fieldname| : name of the field, i.e. call is seti.(fieldname)
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayes messages (0: no, 1 or greater: yes).
%
%% Output Arguments
%
% A message, e.g.:
%
% |Parameter "dirOutput" was not set. Set dirOutput = output.|
%
%% Code
%
function setmessage(struct,fieldname,varargin)
if nargin == 3
    dispDepth = varargin{1};
else
    dispDepth = 0;
end
if dispDepth >= 1
    fprintf('   Set parameter %s = %s.\n',fieldname,num2str(struct.(fieldname)));
end

end
