function [z, Om, T] = constructPhi(varargin)
% [z, Om, T] = constructPhi( [oo], [S], [T], [p], [options]) (ordering of p and T can be interchanged)
% [z, xyzv] = constructPhi( [oo], p, 'interpolate', 'grid',val, [options])
% [z, xyzv] = constructPhi([oo], S, [p], 'interpolate', [options])
% Computes values of the vectorized basic limit function Phi.
% This function is nearly untested
%
% General Options:
%         'verbose',val     Defines the verbose level. Default=0
%         'Omega',val       dim x N integer matrix, Index set used in/for the construction of the transition matrices. 
%                                   For default:       If 'Omega' is not given, then it will get constructed.
%                                   For 'interpolate': If 'Omega' is given Phi(p) is returned, otherwise phi(p) is returned.
%
% 1) Computes the value of Phi at the point defined by the digit/subdiv-op sequence oo OR at the point p with subdiv-op sequence oo
%   Input:
%       oo                  {1xN, 1xM} if p given, {2xN,2xM} otherwise. Sequence of digits to be chosen. 
%                           Can be ommited if S consists of only one subdivision operator and p is given.
%       S, T                Subdivision operators and (non-flat!)-transition matrices. One of these two must be given. It is recommended to give both.
%       p                   point (dim x 1 - vector) where the function shall be evaluated if oo={1xN,1xM}. If oo={2xN,2xM} then p is ignored.
%
%   Options:
%       'eps',val           Used to test whether the product of transition matrices converges or not. Default: 1e-12
%       'sym'               The computation is done symbolically, return value is double
%     Output: 
%       'z'                 the value of Phi
%       'Om'                the set Omega
%       'T'                 the (non-flat)-transition matrices
%
% 3) 'interpolate': Computes the value of Phi/phi at  by interpolating the grid xyzv
%   Input:
%       oo
%       'p'                 Point where Phi/phi shall be evaluated
%       'interpolate'       Obligatory option
%    Output: 
%       'z'                 the value of Phi
%       'Om'                the set Omega (if 'Omega' is given)
%
% See also: blf, transitionmatrix, num2ordering
% 
% Written by: tommsch, 2017

% XX [z, Om, T] = constructPhi(p , T | S, [options]), machen dass man auch direkt einen Punkt geben kann. Ordering wird dann in dieser Funktion ausgerechnet.
% XX              Für diesen Fall brauche ich wohl Omega
% XX Überall Argument 'sym' entfernen. Wenn ein input symbolic ist, soll symbolisch gerechnet werden und symbolisch ausgegeben werden
% XX Über 'interpolate' nachdenken

[Om,varargin]=parsem('Omega',varargin,[]);
[verbose,varargin]=parsem({'verbose','v'},varargin,0);
[eps,varargin]=parsem('eps',varargin,1e-12);

S=[]; T=[]; oo={[],[1]}; p=[];
if(isordering(varargin{1})); oo=varargin{1}; varargin(1)=[]; end;
if(isS(varargin{1})); S=varargin{1}; varargin(1)=[]; elseif(isT(varargin{1})); T=varargin{1}; varargin(1)=[]; end
if(isS(varargin{1})); S=varargin{1}; varargin(1)=[]; elseif(isT(varargin{1})); T=varargin{1}; varargin(1)=[]; end
if(ispoint(varargin{1})); p=varargin{1}; varargin(1)=[]; end;

if(parsem('interpolate',varargin));  %If interpolation shall be used
    xyzv=parsem('grid',varargin,[]);
    
    if(isempty(xyzv));
        [~,xyzv]=blf(oo,S,'plot',0);
    end;
    if(~isempty(Om));
        z=interpolate_PHI(xyzv,p,Om);
    else;
        z=interpolate_phi(xyzv,p);
    end
    
else;  %If Transition matrices shall be used
    
    if(size(oo{1},1)==1 || size(oo{2},1)==1); 
        oo=num2ordering(oo,S,p);
    elseif(~isempty(p))
        vprintf('err','constructPhi: Both a point ''p'' and a ordering ''oo'' is given. The ordering may not coincide with the point, and thus point will be ignored.',1,verbose);
    end;
    if(isempty(T)); [T,Om]=transitionmatrix(S,'Omega',Om,'noflat'); end;

    if(parsem('sym',varargin))
        z1=tbuildproduct(T,oo{1},'reverse','sym');
        z2=tbuildproduct(T,oo{2},'reverse','sym');
        [z2,d2,p2]=eig(z2); %XX unschön gemacht mit p2 hier
        d2=d2(p2,p2);
    else
        z1=tbuildproduct(T,oo{1},'reverse');
        z2=tbuildproduct(T,oo{2},'reverse');
        [z2,d2]=eig(z2);
    end
    idx=find(abs(diag(d2)-1)<eps);
    if(length(idx)>1); vprintf('err','Subdivision operator not convergent or too slowly convergent.', 1,verbose); end;
    idx=idx(1);
    z2=z2(:,idx); 
    z=z1*z2;
    z=z/sum(z);
end;

end

function [v] = interpolate_phi( xyzv,p)
%returns the function value of phi by interpolating the grid xyzv
    dim=size(xyzv,1)-1;
    if(dim==1);
        xyzv = sortrows(xyzv',1);
        v = griddedInterpolant(xyzv(:,1),xyzv(:,2));
        v=v(p);
    else
        v=griddatan(xyzv(1:dim,:)',xyzv(dim+1,:)',p');
    end
    
    v(isnan(v))=0;
end


function [v] = interpolate_PHI( xyzv,p,Om)
%returns the function value of Phi by interpolating the grid xyzv. Calls interpolate_phi.
v=interpolate_phi(xyzv,setplus(p,Om));
end

function x = ispoint(varargin)
% x = ispoint(p,[options])
% returns true if p is nx1 vector
%   p             Point to be tests
%   'dim',val    (integer greaterequal zero) ispoint returns true if p is a dim x 1 vector. 
% 
% Eg: ispoint(zeros(0,1),'dim',0)

x=0;
p=varargin{1};
if(~isnumeric(p)); return; end;
if(~iscolumn(p)); return; end;
dim=parsem('dim',varargin,-1);
if(dim>=0 && ~isequal(size(p),[dim 1])); return; end;
x=1;
    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 