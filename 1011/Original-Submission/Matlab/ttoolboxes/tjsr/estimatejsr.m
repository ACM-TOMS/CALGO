function JSR = estimatejsr(M);  
% JSR = estimatejsr(M);  
% Rough estimate of the joint spectral radius.
% Input
%   M       cell array of matrices
%
% Output
%   JSR     interval containing the JSR
%
% E.g.: estimatejsr({2, 3})
%
% See also: tjsr, findsmp
%
% Written by tommsch, 2018

d = max( 3, ceil(log(1000)/log(size(M,2))) );
[~,~,info] = findsmp( M, 'gripenberg', 'maxsmpdepth',d, 'delta',.8, 'verbose',0, 'maxtime',0 );
%[~,~,info] = findsmp(M,'gripenberg','maxsmpdepth',inf,'delta',.95,'verbose',0,'maxtime',1);
JSR = info.jsrbound;

end 

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   