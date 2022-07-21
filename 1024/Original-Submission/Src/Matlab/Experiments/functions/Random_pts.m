function [mat_A]=Random_pts(m,n,varargin)
%************ Randomly generate data points *************
% 
% [mat_A]=Random_pts(m,n,vargin)
%
% Inputs:
% m=dimension of data
% n=number of points 
% vargin:
% 'normal' : Use the standard Gaussian distribution
% 'unif' : Use the standard uniform distribution
% 'unit ball': Points are uniformly generated on a unitball
% Output:
% mat_A: mxn matrix wil columns as points.
%------------------------------------------------------------------


if ~isempty(varargin)
    [Type]=deal(varargin{:});
else
    Type='normal';
end


if strcmp(Type, 'normal')==1
    mat_A=random('Normal',0,1,m,n);%generate iid gaussian 
elseif strcmp(Type, 'unif')==1
    mat_A=random('unif',0,1,m,n);  %generate iid uniform
elseif strcmp(Type, 'sparse')==1
    mat_AA=random('unif',0,1,m,n);  %generate iid uniform
    mat_AA(mat_AA>0.7)=1;
    mat_AA(mat_AA<=0.7)=0;
    mat_BB=rand(m,n);
    mat_A=mat_AA.*mat_BB;
elseif strcmp(Type, 'unit ball')==1
    mat_A_raw=random('Normal',0,1,m,n); %generate iid gaussian 
    norm_mat_A_raw=sqrt(sum(mat_A_raw.^2,1));
    mat_A= mat_A_raw*diag(1./norm_mat_A_raw);%scale points to uniball
end