function [ dim, Xmuf ] = dimVVt(varargin)
% [ dim, Xmuf ] = dimVVt( Om, [k], [options] );
% Computes the dimension of the spaces V_k and \bar{V}_k, as described in Charina, Mejstrik, 2018.
%
% Input: 
%   Om          array of column vectors, the set for which Vtilde_k and V_k shall be constructed as an 
%   k           index k, optional. Either an integer or a vector of integers (greater equal zero)
%               If k is empty, then all spaces where dim V_k>0 are computed.
%
% Options:
%   '01'        give input as a 0/1 matrix
%   'verbose'   more text output
%   'V'         only compute dimensions of spaces V_k
%   'Vt'        only compute dimensions of spaces \bar{V}_k
%
% Output:
%   dim         first line: dimensions of V_k, for k=0,1,2,3,...
%               second line: dimensions of \bar{V}_k, for k=0,1,2,3,...
%   Xmuf        vector, If set, all sets X_mu are non-empty, see [Mejstrik, PhD-thesis, 2019].
%               only returned if option 'V' is not given
%
% Note:
%  If both options 'V' and 'Vt' are given, the behaviour is undefined.
%
% E.g.: dimVVt([1 1 0 0;1 0 1 1;0 1 0 0],'01','verbose',1)
%
% See also: constructV, constructVt
%
% Written by: tommsch, 2017

%XX Write function which tests wheter $X_\mu(\Omega)$ is non-empty for all $\abs\mu=k$.

[verbose,varargin]=parsem({'verbose','v'},varargin);
[Vflag,varargin]=parsem('V',varargin);
[Vtflag,varargin]=parsem('Vt',varargin);


if(~Vtflag); 
    V=constructV(varargin{:});
    dimV=cellfun(@rank,V); 
end;
if(~Vflag); 
    [Vt,~,Xmuf]=constructVt(varargin{:});
    dimVt=cellfun(@rank,Vt); 
end;


if(verbose); 
    vdisp(varargin{1}); end;
if(~Vtflag && ~Vflag);
    if(size(dimV,2)>size(dimVt,2)); 
        dimVt(size(dimV,2))=0; end;
    dim=[dimV;dimVt];
elseif(Vtflag)
    dim=dimVt;
elseif(Vflag);
    dim=dimV;
end;

end
