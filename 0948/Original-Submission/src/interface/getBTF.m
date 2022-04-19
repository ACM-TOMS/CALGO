function [pe,pv,cb,fb] = getBTF(sadata)
%getBTF returns block triangular form (BTF) of a DAE.
%[pe,pv,cb,fb] = getBTF(sadata)
%sadata is an object returned by daeSA.
%
%(1) If the DAE is structurally well-posed (SWP), getBTF produces coarse 
%and fine BTFs: irreducible BTFs based on the sparsity patterns S of the 
%signature matrix and S_0 of the Jacobian. They are unique up to possible 
%reordering. See the algorithm paper.
%
%Assume S=getSigma(sadata). Then pe and pv are row and column permutation 
%vectors, respectively, such that B=S(pe,pv) is in both coarse and fine 
%BTF. cb and fb are row vectors indicating coarse and fine block 
%boundaries, respectively. That is,
%- the (i,j)th coarse block is B(cb(i):cb(i+1)-1,cb(j):cb(j+1)-1), and
%- the (i,j)th fine block is B(fb(i):fb(i+1)-1,fb(j):fb(j+1)-1).
%
%(2) If the DAE is structurally ill-posed (SIP), getBTF produces diagnostic
%BTF, which identifies a structurally under-determined set of variables and
%a structurally over-determined set of equations.
%
%Also assume S=getSigma(sadata). B=S(pe,pv) is of the form (blanks are
%blocks of -inf's)
%
% B = [B11 B12 B13 B14
%              B23 B24
%                  B34
%                  B44],
%
%where the block boundaries are described by fb, a 2-by-5 positive integer
%matrix: the (i,j)th block is B(fb(1,i):fb(1,i+1)-1, fb(2,j):fb(2,j+1)-1).
%In this case, cb is the empty array [].
%
%The blocks B12, B23, and B34 are square with finite diagonal entries. 
%[B11 B12] has more columns (variables) than rows (equations), and 
%[B34;B44] has more rows than columns. Hence (the subsystem leading to),
%- block B11 is under-determined,
%- blocks B12, B23, and B34 is well-determined, and
%- block B44 is over-determined.
%
%Example 1:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       [pe,pv,cb,fb] = getBTF(sadata)
%       S = getSigma(sadata);
%       B = S(pe,pv)
%
%produces
%
%pe =
%     5     4     6     3     2     1
%pv =
%     5     6     4     1     2     3
%cb =
%     1     4     7
%fb =
%     1     2     3     4     7
%B =
%     3     0  -Inf  -Inf  -Inf  -Inf
%  -Inf     0     2  -Inf  -Inf  -Inf
%     0  -Inf     0  -Inf  -Inf     2
%  -Inf  -Inf  -Inf     0     0  -Inf
%  -Inf  -Inf  -Inf  -Inf     2     0
%  -Inf  -Inf  -Inf     2  -Inf     0
%
%The course blocks are 
%B(cb(1):cb(2)-1,cb(1):cb(2)-1) = B(1:3,1:3) and 
%B(cb(2):cb(3)-1,cb(2):cb(3)-1) = B(4:6,4:6).
%
%The fine blocks are 
%B(fb(1):fb(2)-1,fb(1):fb(2)-1) = B(1:1,1:1),
%B(fb(2):fb(3)-1,fb(2):fb(3)-1) = B(2:2,2:2),
%B(fb(3):fb(4)-1,fb(3):fb(4)-1) = B(3:3,3:3), and
%B(fb(4):fb(5)-1,fb(4):fb(5)-1) = B(4:6,4:6).
%
%Example 2:
%       sadata = daeSA(@illPosed3,6);
%       [pe,pv,cb,fb] = getBTF(sadata)
%       S = getSigma(sadata);
%       B = S(pe,pv)
%
%produces
%
% pe =
%      6     1     2     3     4     5
% pv =
%      2     3     1     4     5     6
% cb =
%      []
% fb =
%      1     2     4     5     7
%      1     3     4     6     7
% B =
%      0     0     0     0     0     0
%   -Inf  -Inf  -Inf     0     0     0
%   -Inf  -Inf  -Inf     0     0     0
%   -Inf  -Inf  -Inf  -Inf  -Inf     0
%   -Inf  -Inf  -Inf  -Inf  -Inf     0
%   -Inf  -Inf  -Inf  -Inf  -Inf     0
%
%That is,
%- the under-determined system comprises f6 in x2,x3,x1,
%- the well-determined system comprises f1,f2 in x4,x5,
%and
%- the over-determined system comprises f3,f4,f5 in x6.
%
%See also daeSA, isSWP, dmperm.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

if nargin>4
    error('??? Error using ==> getBTF\nToo many output arguments.')
end

[pe,pv,cb,fb] = DAESAgetBTF(sadata);
end