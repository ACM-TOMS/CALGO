%% gmresKelley
% GMRES linear equation solver.
%
%% Syntax
%
%   x = gmresKelley(x0, b, atv, params)
%
%% Description
%
% |x = gmresKelley(x0, b, atv, params)| computes the solution x of linear
% equation |atv(x) = b| with matrix-vector product routine |atv|, and right
% hand side |b|. The initial vector is x0. Further parameters are given in
% params.
%
% Implementation following Saad-Schultz.
%
% C. T. Kelley, July 10, 1994
%
% This code comes with no guarantee or warranty of any kind.
%
% Output Arguments "error" and "total_iters" can be added by changing the 
% code in |function [x, error, total_iters] = gmres(x0, b, atv, params)|.
%
% Requires the function |givapp|, that is given as subfunction.
%
%% Example
%
%   x0 = [1; 1];
%   A = [1 2; 3 4];
%   b = [4.5; 9.5];
%   atv = @(x) A*x;
%   tol = 1E-6;
%   x = gmresKelley(x0, b, atv , [tol ,100,1]);
%
% _Result:_
%
%   x =
%
%    0.5000
%    2.0000
%
%
%% Input Arguments
%
% * x0  :   initial iterate
% * b   :   right hand side
% * atv :   a matrix-vector product routine
%           atv must return Ax when x is input
%           the format for atv is
%           function ax = atv(x).
%           Note that for GMRES we incorporate any 
%           preconditioning into the atv routine.
% * params : three dimensional vector to control iteration
% * params(1) = relative residual reduction factor, 
% * params(2) = max number of iterations, 
% * params(3) (Optional) = reorthogonalization method:
%                   1 -- Brown/Hindmarsh condition (default), 
%                   2 -- Never reorthogonalize (not recommended), 
%                   3 -- Always reorthogonalize (not cheap!).
%
%% Output Arguments
%
% * x   :   solution
%
% _Additional outputs after code changing:_
%
% * error   :   vector of residual norms for the history of
%               the iteration.
% * total_iters : number of iterations.
%
%
%% License Information
%
% Note that the source code of this file was not written by the authors of
% the toolbox IPscatt. It contains the GMRES by Kelley that is faster than the standard GMRES.
%
% This method is describe in Ch. 3.4 of the following book:
% C. Tim Kelley. Iterative Methods for Linear and Nonlinear Equations. Society for Industrial and Applied Mathematics, 1995.
% 
% The source code is available at:
% C. Tim Kelley. Iterative methods for linear and nonlinear equations. companion software, 2002.
% <https://de.mathworks.com/matlabcentral/fileexchange/2198-iterative-methods-for-linear-and-nonlinear-equations/content/kelley/gmres.m>
% Accessed: March 2017.
% 
% The corresponding license is available at:
% <https://de.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=2198>
% 
% In the following we will cite this license (Accessed: May 2017):
% 
% *License*
% 
% Copyright (c) 2016, C.T. Kelley
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%% Code
%
function x = gmresKelley(x0, b, atv, params)
%
% initialization
%
n=length(b);
errtol=params(1);
kmax=params(2);
reorth=1;
if length(params) == 3
    reorth=params(3);
end
x=x0;
%
%
h=zeros(kmax);
v=zeros(n,kmax);
c=zeros(kmax+1,1);
s=zeros(kmax+1,1);
if norm(x) ~=0
   r = b-feval(atv,x);
else
   r = b;
end
rho=norm(r);
g=rho*eye(kmax+1,1);
errtol=errtol*norm(b);
error=[];
%
% test for termination on entry
%
error=[error,rho];
total_iters=0;
if(rho < errtol) 
    return
end
%
v(:,1)=r/rho;
beta=rho;
k=0;
%
% GMRES iteration
%
while((rho > errtol) & (k < kmax))
    k=k+1;
    v(:,k+1)=feval(atv,v(:,k));
    normav=sqrt(v(:,k+1)'*v(:,k+1)); %originally norm(v(:,k+1));
%
% Modified Gram-Schmidt
%
    for j=1:k
        h(j,k)=v(:,j)'*v(:,k+1);
        v(:,k+1)=v(:,k+1)-h(j,k)*v(:,j);
    end
    h(k+1,k)=sqrt(v(:,k+1)'*v(:,k+1)); %originally norm(v(:,k+1));
    normav2=h(k+1,k);
%
% Reorthogonalize?
%
if  (reorth == 1 & normav + .001*normav2 == normav) | reorth ==  3
    for j=1:k
        hr=v(:,j)'*v(:,k+1);
        h(j,k)=h(j,k)+hr;
        v(:,k+1)=v(:,k+1)-hr*v(:,j);
    end
    h(k+1,k)=norm(v(:,k+1));
end
%
%   watch out for happy breakdown 
%
    if(h(k+1,k) ~= 0)
         v(:,k+1)=v(:,k+1)/h(k+1,k);
    end
%
%   Form and store the information for the new Givens rotation
%
    if k > 1
        h(1:k,k)=givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1);
    end
    nu=sqrt(h(k:k+1,k)'*h(k:k+1,k)); %originally norm(h(k:k+1,k));
    if nu~=0
%        c(k)=h(k,k)/nu;
        c(k)=conj(h(k,k)/nu);
        s(k)=-h(k+1,k)/nu;
        h(k,k)=c(k)*h(k,k)-s(k)*h(k+1,k);
        h(k+1,k)=0;
        g(k:k+1)=givapp(c(k),s(k),g(k:k+1),1);
    end
%
% Update the residual norm
%
    rho=abs(g(k+1));
    error=[error,rho];
end
%
% At this point either k > kmax or rho < errtol.
% It's time to compute x and leave.
%
y=h(1:k,1:k)\g(1:k);
total_iters=k;
x = x0 + v(1:n,1:k)*y;

return

%%
% *Code: subfunction givapp*
%
function vrot=givapp(c,s,vin,k)
%  Apply a sequence of k Givens rotations, used within gmres codes
% 
%  C. T. Kelley, July 10, 1994
%
% This code comes with no guarantee or warranty of any kind.
%
%  function vrot=givapp(c, s, vin, k)
%
vrot=vin;
for i=1:k
    w1=c(i)*vrot(i)-s(i)*vrot(i+1);
%    w2=s(i)*vrot(i)+c(i)*vrot(i+1);
    w2=s(i)*vrot(i)+conj(c(i))*vrot(i+1);
    vrot(i:i+1)=[w1,w2];
end