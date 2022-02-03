function [cql, fql] = getQLdata(sadata)
%getQLdata returns information about the quasilinearity of each diagonal 
%block in the coarse and fine block-triangularizations of a DAE structure.
%cql = getQLdata(sadata)
%[cql,fql] = getQLdata(sadata)
%sadata is an object returned by daeSA.
%
%cql and fql are logical row vectors associated with the coarse and fine 
%block triangularizations, respectively. Assume that the DAE is structurally
%well-posed (SWP), and let
%[pe,pv,cb,fb] = getBTF(sadata);
%
%For the coarse block triangularization, if cql(i)=true, then the block
%comprising pe(cb(i):cb(i+1)-1) in variables pv(cb(i):cb(i+1)-1) is quasilinear
%in leading derivatives of these variables; otherwise,  if cql(i)=false, 
%this block is non quasilinear.
%
%For the fine block triangularization, if fql(i)=true, then the block 
%comprising pe(fb(i):fb(i+1)-1) in variables pv(fb(i):fb(i+1)-1) is
%quasilinear in leading derivatives of these variables;
%otherwise, if fql(i)=false, it is non quasilinear.
%
%If the DAE is not SWP, then cql=[] and fql=[].
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       [pe,pv,cb,fb] = getBTF(sadata)
%       [cql,fql] = getQLdata(sadata)
%
%produces
%pe =
%     5     4     6     3     2     1
%pv =
%     5     6     4     1     2     3
%cb =
%     1     4     7
%fb =
%     1     2     3     4     7
%cql =
%     0     1
%fql =
%     0     1     0     1
%
%
%For the coarse block triangularization, this implies that the block 
%comprising equations pe(1:3) in variables pv(1:3) is non quasilinear,
%and the block comprising equations pe(4:6) in variables pv(4:6) is 
%quasilinear.
%
%For the fine block triangularization, this implies
%  pe(1:1) = 5 in pv(1:1) = 5, fql(1) = 0, non quasilinear
%  pe(2:2) = 4 in pv(6:6) = 6  fql(2) = 1, quasilinear 
%  pe(3:3) = 6 in pv(4,4) = 4  fql(3) = 0, non quasilinear
%  pe(4:6) = 3,2,1 in pv(4:6) = 1,2,3 fql(4) = 1, quasilinear
%
% See also daeSA, isSWP, getBTF.
%
% Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

[cql, fql] = DAESAgetQLdata(sadata);
end
