function p = fl(k)
%
%  generalization of Farmer-Loizou example:
%
%    (x-1)^4k * (x-2)^3k * (x-3)^2K * (X-4)^k
%

p = poly([ones(1,4*k),2*ones(1,3*k),3*ones(1,2*k),4*ones(1,k)]);
