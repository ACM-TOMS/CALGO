function M = tbuildproduct(varargin)
% M = tbuildproduct( A, oo, [options] )
% Constructs the product, of matrices in the cell A corresponding to the sequence in mo
% This function has a big overhead. For many computations use tbuildproduct_fast
%
% Input:
%   A       cell array of matrices
%   oo      The ordering. Either 
%               ordering: {[1xN],[1xM]}, if A is a cell array of matrices, N or M can be zero, or 
%               row or column vector, if A is a cell array of matrices
%               (experimental) ordering: {[2xN],[2xM]}, if A is a cell-array of cell-arrays of matrices, N or M can be zero, or
%
% Options:
%         'sym'              the matrices A are converted to symbolic matrices prior to computation
%         'l',val            approximate length of the expanded product sequence (default 1000)
%         'reverse'          reverses the product (e.g. may be used in the computation with transition matrices)
% Info:
%   Constructs the product, of matrices in the cell A corresponding to the sequence in mo
%   i.e.: tbuildproduct({M1 M2},[1 2 2]) == M2*M2*M1
%         tbuildproduct({M1,M2}, {[],[1 2]}) == (M2*M1)^inf
%   if oo(i)<=0, A{oo(i)} is replaced by the identity
%
% E.g.: tbuildproduct({[1 1/2; 0 1/2],[1/2 0; 1/2 1]}, [1 2 2 1])
%       tbuildproduct({[0.8 0.2; 0.2 0.8],[0.4 0.6;0.4 0.6]}, {[],[1 2]})
%     
%
% See also: tbuildproduct_fast, ordering2num
%
% Written by tommsch, 2018

% XX: tbuildproduct(M,ooo,'reverse') berechnet symbolisch, obwohl es das nicht sollte
% XX Infinite products not implemented correctly

[reverseflag,varargin]=parsem('reverse',varargin);
[symflag,varargin]=parsem('sym',varargin);
[l,varargin]=parsem('l',varargin,1000);

A=varargin{1};
oo=varargin{2};
varargin(1:2)=[];
parsem(varargin,'test');

if(~iscell(oo) && iscolumn(oo)); 
    oo=oo.'; end; %make oo to row vector
if(~iscell(oo)); 
    val=oo; oo={val,[]}; end;  %make oo to {oo,[]}
if(size(oo,2)==1); 
    oo={oo{1},[]}; end;      %make {oo1} to {oo,[]}
LEN{1} = size(oo{1},2); LEN{2} = size(oo{2},2); %Length of sequences
l=ceil(l/LEN{2}); %sic!, if LEN{2}==0, then we do not need this value

if(iscell(A{1}))
    dim = size(A{1}{1},1);
    if(size(oo{1},1)==1); 
        error('tbuildproduct: ''oo'' needs a 2xN array.'); end;    
else
    dim = size(A{1},1);
    val=A; A={}; A{1}=val;
    if(~isempty(oo{1})); 
        oo{1}(2,:)=1; end; 
    oo{1}=flip(oo{1},1);
    if(~isempty(oo{2})); 
        oo{2}(2,:)=1; end; 
    oo{2}=flip(oo{2},1);
end
if(symflag); 
    for i=1:size(A,2); 
        for j=1:numel(A{i});A{i}{j}=sym(A{i}{j}); end; end; 
    M1 = sym(eye(dim));
    M2 = sym(eye(dim));
else
    M1 = eye(dim);
    M2 = eye(dim);
end;



if(reverseflag); 
    oo{1}=flip(oo{1},2); 
    oo{2}=flip(oo{2},2); 
end;

for t=1:LEN{1};
    j=oo{1}(1,t); %number of set 
    i=oo{1}(2,t); %number of transition in set j    
    if(i>0 && j>0);
        M1 = A{j}{i}*M1; end;
end;   
for t=1:LEN{2};
    j=oo{2}(1,t); %number of set 
    i=oo{2}(2,t); %number of transition in set j    
    if(i>0 && j>0);
        M2 = A{j}{i}*M2; end;
end;

if(~isempty(oo{2}))
    if(reverseflag);
        M=M1*M2^l;
    else;
        M=M2^l*M1;
    end;
else
    M=M1;
end

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   