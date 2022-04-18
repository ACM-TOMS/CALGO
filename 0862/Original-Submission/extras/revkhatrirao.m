function P = revkhatrirao(varargin)
%REVKHATRIRAO Reverse Khatri-Rao product
%
%   REVKHATRIRAO(A1,A2,...) computes the Khatri-Rao product of
%   multiple matrices in reverse order.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR,KHATRIRAO.

%Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
%Sandia National Laboratories, 2004.  Please address questions or
%comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
%distribute, display, and use this work, under the following
%conditions. (1) You must give the original authors credit. (2) You may
%not use or redistribute this work for commercial purposes. (3) You may
%not alter, transform, or build upon this work. (4) For any reuse or
%distribution, you must make clear to others the license terms of this
%work. (5) Any of these conditions can be waived if you get permission
%from the authors.

if nargin==1
    A = varargin{1};
else
    A = varargin;
end

P = khatrirao(A{end:-1:1});
