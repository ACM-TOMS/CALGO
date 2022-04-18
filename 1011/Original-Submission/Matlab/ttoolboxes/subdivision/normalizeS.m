function [ S ] = normalizeS(varargin)
% [ S ] = normalizeS(S, [options])
% Normalizes the values of the masks a in S.
%
%Input: 
%   S               Cell array of subdivision operators
%
% Options: 
%   'equal'         entries of the mask corresponding to representatives of a digit set are equal and sum up to one
%   'gauss'         entries of the mask are preserved, but gaussian like scaled
%   'scale'         entris are just scaled (Default)
%   'verbose',val   sets the verbose level
%
% Output:
%   S               Cell array of subdivision operators with normalized masks
%
% E.g.: S=normalizeS({[10 4 6 4 1]',2},'equal'); vdisp(S);
%       S=normalizeS({[10 4 6 4 1]',2},'gauss'); vdisp(S);
%       S=normalizeS({[10 4 6 4 1]',2},'scale'); vdisp(S);
%
% Written by: tommsch, 2017

 %#ok<*ALIGN>

S=varargin{1}; varargin(1)=[];
S=getS(S,'nocheck');
[equal,varargin]=parsem('equal',varargin);
[gaussian,varargin]=parsem('gauss',varargin);
[verbose,varargin]=parsem({'verbose','v'},varargin,1);
[scale,varargin]=parsem('scale',varargin); %#ok<ASGLU>

parsem( varargin, 'test' );

sizeS = size(S,1);
for i = 1:sizeS
    a = shrink(S{i,1});
    M = S{i,2};
    suppa = a.supp; %kth support
    CENTER = sum(suppa,2)/size(suppa,2); %center of mass of all entries
    Cl = constructdigit(M,'ZZ',suppa,'classify','sym'); %Class
    [C, ~, ic] = unique(Cl','rows');
    C=C';
    if( size(C,2)~=round(abs(det(M))) ); 
        warning( 'normalizeS:support', 'Cannot normalize mask, since the support of the mask is to small for index %i.', i );
        S={}; 
        return; end;
    anew = zerosm( size(a) ); %new mask
    MIN = min( suppa,[], 2 ); %support of the new mask
    MAX = max( suppa,[], 2 );

    %weight matrix
    if( ~equal ); 
        L = cell( 1, size(MIN,1) );
        for j = 1:numel(L);  
            L{j} = (MIN(j):MAX(j))-CENTER(j);
            L{j} = exp(-abs(L{j})); end;
        W = convm(L{:},'outer');
    else
        W = onesm((MAX-MIN)'+1); end;

    for j=1:size(C,2) %go through all classes
        submask = suppa(:,ic==j);
        %num=size(submask,2);
        asub = characteristic(submask,'amin',MIN,'amax',MAX).*W;

        if(gaussian); 
            ak=a.c;
            asub=asub.*ak; 
        elseif(equal); 
            %do nothing
            %asub=asub/sum(asub(:));
        else %scaled
            ak = a.c;
            asub(logical(asub)) = ak(logical(asub)); end;
        asub = asub/sum(asub(:));
        
        if( anym(~isfinite(asub)) )
            warning( 'normalizeS:failure', 'normalizeS: Cannot normalize mask for index %i.\n', i );
            S={};
            return; end;

        if( any(size(anew)-size(asub)) );
            warning( 'normalizeS:dim', 'normalizeS: Wrong dimension for index %i.', i )
            S={};
            return; end;
        anew = anew+asub; end;
    S{i,1}.c = anew; end;


end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.