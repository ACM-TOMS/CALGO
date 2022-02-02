function [ fout , sglout ] = transf( fin , sglin , ival )

%TRANSF Transformation from a real interval onto the interval [-1,1]
%   [FOUT,SGLOUT]=TRANSF(FIN,SGL,IVAL) transforms a (matrix-valued) 
%   function FIN with singularities in SGLIN, defined on the interval(s) 
%   IVAL of the form [a_i,b_i], with -inf =< a_i < b_i =< inf, into the 
%   function(s) FOUT_i with singularities SGLOUT_i defined on the interval 
%   [-1,1] such that int(FIN(x),x=a_i..b_i) == int(FOUT_i(x),x=-1..1). IVAL 
%   should be a matrix of size N x 2, where N denotes the number of 
%   intervals, with IVAL(i,1)<IVAL(i,2) for i=1,...,N.

sglin = sglin(:);
sizeI = size(ival);
if ~(sizeI(2)==2),
    error('IVAL needs to be a vector/matrix of size N x 2.');
elseif max(~(ival(:,1)<ival(:,2))),
    error('IVAL should be of the form [a,b] with a<b.');
else
    if ischar(fin),
        F = str2func(fin);
    else
        F = fcnchk(fin);
    end
    sizeF = size(F((ival(1,1)+ival(1,2))/2));
end

ival1 = ival(1,:);
[ fout , sglout ] = tfsingle( fin , sglin , ival1 );
for k=2:sizeI(1),
    ival1 = ival(k,:);
    [ f , s ] = tfsingle( fin , sglin , ival1 );
    if sizeF(1)>sizeF(2),
        fout = @(x) [fout(x), f(x)];
    else 
        fout = @(x) [fout(x); f(x)];
    end
    sglout = [sglout;s];
end

end

% -------------------------------------------------------------------------
function [ fout , sglout ] = tfsingle( fin , sglin , ival )
% Transformation for one single interval

k = isinf(ival(1))+2*isinf(ival(2));
switch k,
    case 0,
        su = sum(ival); di = diff(ival);
        tau = @(x) (di*x+su)/2; 
        taup = @(x) di/2;
        sglout = (2*sglin-su)/di; 
    case 1,
        su = ival(2)-1; di = ival(2)+1;
        tau = @(x) (di*x+su)./(x+1); 
        taup = @(x) 2./((x+1).^2);
        sglout = (sglin-su)./(di-sglin); 
    case 2,
        su = -(ival(1)+1); di = ival(1)-1;
        tau = @(x) (di*x+su)./(x-1); 
        taup = @(x) 2./((x-1).^2);
        sglout = (sglin+su)./(-di+sglin);
    otherwise
        tau = @(x) x./(1-x.^2);
        taup = @(x) (1+x.^2)./((1-x.^2).^2);
        di = length(sglin);
        sglout = (-1+sqrt(1+4*sglin.^2))./(2*sglin);
        sglout = [sglout;(-1-sqrt(1+4*sglin.^2))./(2*sglin)];
end
fout = @(x) fin(tau(x)).*taup(x);

end