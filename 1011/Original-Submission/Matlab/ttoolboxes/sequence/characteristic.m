function [ chi, MIN, idx ] = characteristic( varargin )
% [ chi, MIN, idx ] = characteristic( D, [options] )
% Returns the characteristic function of an index set.
%
% Input:
%   D               column vectors defining the coordinates
%
% Options:
%   'amin',val      Index of first entry in D
%   'amax',val      Index of last entry in D
%
% Output:
%   chi             0-1 array of dimensions size(D,1)
%   MIN             index of first entry in chi
%   idx             index to D
% 
% Info: This function is much faster than the inverse function supp(). Use this function if possible.
%       For the use of the variable "idx", see the example below.
%
% E.g.: a=characteristic([1 0; 2 0; 1 4]')
%       [a,amin]=characteristic([1 0; 2 0; 1 -1]','amin',[-1; -1],'amax',[3;3])
%       [a,amin]=characteristic([1 0; 2 0; 1 -1]','amin',[0; 0])
%
%       D1=[0 0;2 2; 0 1;1 2;0 2]'
%       [chi,amin,idx]=characteristic(D1);
%       D2=supp(chi,2,amin);
%       D3=D2(:,idx)
%       The elements of D3 and D1 are in the same order
%
% See also: supp
%
% Written by: tommsch, 2019

%#ok<*ALIGN>

    amin=parsem('amin',varargin,[]);
    amax=parsem('amax',varargin,[]);

    D=varargin{1};

    [D,idx]=sortrows(flipud(D).');
    D=flipud(D.');
    [~,idx]=sort(idx);
    %change sorting of index such it is compatibly with linear indexing

    MIN=min(D,[],2); 
    if(~isempty(amin)); 
        if(any(amin>MIN))
            error('''amin'' too big.'); end;
        MIN=amin; 
    end;
    MAX=max(D,[],2); 
    if(~isempty(amax)); 
        if(any(amax<MAX))
            error('''amax'' too small.'); end;
        MAX=amax; 
    end;


    if(~isequal(round(D),D));  
        chi=[]; 
        return; end;
    dim=size(D,1);
    if(dim>1); 
        chi=zeros(double(1+(MAX-MIN)'));
    else; 
        chi=zeros(double(1+(MAX-MIN)'),1);
    end;

    D=setplus(D,-MIN)+1;
    Dcell=num2cell(D,1);

    for i=1:size(Dcell,2);
        Dcelli=num2cell(Dcell{i});
        chi(Dcelli{:})=1;
    end

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 