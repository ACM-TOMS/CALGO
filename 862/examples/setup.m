echo on

% Print our answers in compact format so that more fits on the
% screen. Use set(o,' FormatSpacing', 'loose') to change to back to
% normal. 
set(0, 'FormatSpacing', 'compact');

% Add the directory that contains all the classes to the path. In
% this case, it should be the directory one up.
p = [pwd '/../'];
addpath(p);

echo off
