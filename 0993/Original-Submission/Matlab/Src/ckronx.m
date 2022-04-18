% CKRONX The product of repeated Kronecker products and a matrix. 
% USAGE      
%   [C,order]=ckronx(A,B,options);
% INPUTS
%   A       : a d-element cell array with element i an m(i) x n(i) matrix
%   B       : a compatible matrix prod(n) x p (prod(m) x p if transpose=1)
%   options : an options structure (see below)
%
% Fields in options structure:
%   ind       : a selection vector of numbers on 1..d that reorders 
%                 the elements of A [default: 1:d]
%   transpose : d-vectors of logicals: 1 to use the transpose of A(i) 
%                 a scalar entry will be expanded to a d-vector [default: 0]
%   optorder  : use optimal ordering of operations
%   forward   : 1 forces use of the forward algorithm, 0 forces the backward
%   print     : print information about the operation
% OUTPUTS  
%   C         :  prod(m) x p matrix (prod(n) x p if transpose=1)
%   order     :  1 x d vector with the optimal order of variables
% Solves (A1 x A2 x...x Ad)*B
% where x denotes Kronecker (tensor) product.
% The Ai are passed as a cell array A. 
% A must be a vector cell array containing 2-D numerical arrays (matrices).
% If A is a matrix the function returns A*B (or A'*B if transpose=1)
%
% Example:
%   m=[50 50]; n=[30 60]; 
%   d=length(m);
%   A=cell(1,d);
%   for i=1:d, A{i}=randn(m(i),n(i)); end
%   B=randn(prod(m),5);
%   options=struct('transpose',ones(1,d));
%   C=ckronx(A,B,options);
%
% Alternative input syntax (for backward compatibility):
%   C=ckronx(A,B,ind,transpose);

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2018, Paul L. Fackler (paul_fackler@ncsu.edu)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without  
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%        this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%        this list of conditions and the following disclaimer in the 
%        documentation and/or other materials provided with the distribution.
%    * Neither the name of the North Carolina State University nor of Paul L. 
%        Fackler may be used to endorse or promote products derived from this 
%        software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% For more information, see the Open Source Initiative OSI site:
%   http://www.opensource.org/licenses/bsd-license.php

% Adapted from the ckronx function in the CompEcon Toolbox
%(c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% Modified in 2018 by Paul L. Fackler

function [C,order]=ckronx(A,B,options,transpose)

if nargin<2
  error('At least two parameters must be passed')
end
if ~exist('transpose','var')
  transpose=false; 
end
optorder=0;
useforward=[];
print=0;
order=[];

% handle case if a single numeric matix is passed as A
if ~iscell(A)                         % A is a matrix: return A*B
  if exist('options','var') && ~isempty(options) && isfield(options,'transpose')
    transpose = options.transpose;
  end
  if transpose
    if size(A,1)~=size(B,1)
      error('A and B are not conformable')
    end
    C=A'*B;
  else
    if size(A,2)~=size(B,1)
      error('A and B are not conformable')
    end
    C=A*B;
  end
  return
end

% handle cell array case
d=numel(A);
% get default order
if ~exist('options','var')  || isempty(options)
  ind=1:d;
else
  if isstruct(options)
    ind=1:d;       
    if isfield(options,'ind'),       ind=options.ind;              end  
    if isfield(options,'transpose'), transpose=options.transpose;  end   
    if isfield(options,'optorder'),  optorder=options.optorder;    end   
    if isfield(options,'forward'),   useforward=options.forward;   end   
    if isfield(options,'print'),     print=options.print;          end           
  else
    ind=options;
  end
end
if length(transpose)==1
  transpose=repmat(transpose,1,d); 
end
A=A(ind);
transpose=transpose(ind);
m=zeros(1,d);  % # of rows (cols if transpose(i)=1)
n=zeros(1,d);  % # of cols (rows if transpose(i)=1)
q=zeros(1,d);  % # of non-zeros
for i=1:d, 
  if transpose(i)
    [n(i),m(i)]=size(A{i}); 
  else
    [m(i),n(i)]=size(A{i});
  end
  if issparse(A{i})
    q(i)=nnz(A{i}); 
  else
    q(i)=numel(A{i});
  end
end
if prod(n)~=size(B,1)
  error('A and B are not conformable')
end
% handle single matix case
if d==1, 
  if transpose(1)
    C=A{1}'*B; 
  else
    C=A{1}*B;
  end
  return; 
end

% check if use of either forward or backward is forced
if ~isempty(useforward)
  if useforward~=0 % use forward approach
    C=forward(A,B,d,n,transpose);
  else  % use backward approach
    C=backward(A,B,d,n,transpose);
  end
  return
end
  
cost=@(m,n,q)sum(q.*cumprod([1 m(1:end-1)]).*fliplr(cumprod([1 fliplr(n(2:end))])));
fcost=cost(m,n,q);
bcost=cost(n,m,q);

if optorder
  [~,order]=sort((m-n)./q);
  gcost=cost(m(order),n(order),q(order));
  if gcost+prod(m)+prod(n) < min(fcost,bcost) 
    p=size(B,2);
    C=permute(reshape(B,[fliplr(n) p]),[d+1 d+1-fliplr(order)]);
    n=n(order);
    transpose=transpose(order);
    for i=1:d
      C=reshape(C,numel(C)/n(i),n(i));
      if transpose(i), C=A{order(i)}'*C';
      else             C=A{order(i)}*C';
      end
    end
    C=ipermute(reshape(C,[fliplr(m(order)) p]),[d+1-fliplr(order) d+1]);
    C=reshape(C,numel(C)/size(B,2),size(B,2));
    if print
      fprintf('order: '); fprintf('%1.0f ',order); fprintf('\n');
      disp('reordered cost')
      fprintf('%25.0f\n',gcost);
      disp('reorder costs')
      fprintf('%25.0f\n',[prod(n)*p; prod(m)*p]);
    end
    return
  end
end

% if reordering is not done
if print
  disp('forward & backward costs')
  fprintf('%25.0f\n',[fcost;bcost])
  disp('full Kronecker cost')
  fprintf('%25.0f\n',prod(q));
  if ~optorder
    [~,order]=sort((m-n)./q);
    gcost=cost(m(order),n(order),q(order));
    fprintf('order: '); fprintf('%1.0f ',order); fprintf('\n');
    disp('reordered cost')
    fprintf('%25.0f\n',gcost);
    disp('reorder costs')
    p=size(B,2);
    fprintf('%25.0f\n',[prod(n)*p; prod(m)*p]);
  end
end
if fcost<bcost || (fcost==bcost && all(~transpose)) % use forward approach
  C=forward(A,B,d,n,transpose);
else  % use backward approach
  C=backward(A,B,d,n,transpose);
end

function C=forward(A,B,d,n,transpose)
  C=B';
  for i=1:d
    C=reshape(C,[],n(i));
    if transpose(i), C=A{i}'*C';
    else             C=A{i}*C';
    end
  end
  C=reshape(C,[],size(B,2));
    
    
function C=backward(A,B,d,n,transpose)
  C=B;
  for i=d:-1:1
    C=reshape(C,n(i),[]);
    if transpose(i), C=C'*A{i};
    else             C=C'*A{i}';
    end
  end
  C=reshape(C,size(B,2),[])';