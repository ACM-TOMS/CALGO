function [Mret_oret, v0, v0s,  multiplicity, oorig] = leadingeigenvector(varargin);
% [ Mret, v0, v0s, multiplicity ] = leadingeigenvector ( M, [options] );    
% [ ooret, v0, v0s, multiplicity, oorig ] = leadingeigenvector ( M, oo, [options] );    
% Returns all leading eigenvectors of matrices.
% The leading eigenvectors of a matrix are those, which correspond to a eigenvalue with highest absolute value. 
% There can be more leading eigenvectors to one matrix.
%
% Input:
%   M                   matrix or row-cell of matrices
%   oo                  cell array of orderings, The leading eigenvectors of all products M{oo{1}}(end)*...*M{oo{1}}(1), M{oo{end}}(end)*...*M{oo{end}}(1) are computed
%
% Options:
%   'epsequal',val                      if |x-y|<epsequal*max{x,y} then x==y. Default=automatically determined
%   'nomultipleleadingeigenvector'      returns only one leading eigenvector per matrix
%   'complexeigenvector',val            discards all complex eigenvectors depending on val
%                                           0: complex eigenvectors are kept
%                                           1: complex eigenvectors are removed, if there is also a real eigenvector corresponding to the same product
%                                           2: (default) complex eigenvectors are removed, if there is at least one real eigenvector among all products
%                                           3: Real vectors are computed which span the real subspace of the complex leading eigenvectors
%                                           4: complex eigenvectors are removed always
%   'verbose',val                       sets the verbose level
%   'cycle',val                         default=1, Only applicable if oo is given. Also returns the eigenvalues of all cycles.
%   'repetition',val                    default=0, Only applicable if oo is given. Also returns the eigenvalues of powers of the orderings, if lcm(lengths)<val.
%                                       Sensible value (apart from 0): 50
%   
% Output:
%   Mret{i}/oret{i}     matrix or product of matrices with ordering ooret{i}
%                       ordering of product
%   v0{i}               leading eigenvector. v0{i} is normalized to norm(v0,2)==1.
%   v0s{i}              leading eigenvector of complex transposed M{i}/product (eigenplane of M{i}/product). v0s{i} is normalized to <v0{i},v0s{i}>=1, where v0 is the standard scalar product.
%   multiplicity(i)     multiplicity of the leading eigenvalue 
%   oorig(i)            tells to which original ordering the entry corresponds, in particular: reducelength(ooret(i))==oorig(i)
%
%
%  E.g.: [M,v0,v0s,mult,oo]=leadingeigenvector({[1 1 ; 0 1],[0 0; 1 0]},{[1 2]',[1]'})
%
% Written by: tommsch, 2018

% Changelog: tommsch, 2019-05-27,  Option '3' implemented
%            tommsch, 2019-06-29,  Fixed bug in option '3'
%            tommsch, 2020-04-19,  Behaviour change of option 'cycle'

% XX I dont know if Option 3 is implemented correctly
% XX Make it work for symbolic input

%#ok<*AGROW>

[verbose,varargin] = parsem( {'verbose','v'}, varargin, 1 );
[epsilon,varargin] = parsem( {'epsequal','epsilon','eps'}, varargin, [] );
[nomultipleeigenvector,varargin] = parsem( {'nomultipleeigenvector','nomult'}, varargin,0 );
[complexeigenvector,varargin] = parsem( {'complexeigenvector','cp'}, varargin, 2 );
[cycle,varargin] = parsem( {'cycle','c'}, varargin, 1 );
[repetition,varargin] = parsem( {'repetition','rep','r','power'}, varargin, 0 );

if( ~iscell(varargin{1}) ); 
    M{1} = varargin{1}; 
else; 
    M = varargin{1}; end; %if ordering is given, then M must be a cell array, thus this line does no harm


if( size(varargin,2)==2 ); %if there are only two elements left, an ordering is given
    oo = varargin{2}; 
    oo = [oo; oo]; %copy the original ordering to the second line
    orderinggiven = 1;
    if( cycle ); 
        oo = addcycle( oo ); end; 
    if( repetition ); 
        oo = addrepetition( oo, repetition ); end;    
    noo = size(oo,2);
else; 
    orderinggiven = 0; 
    noo = numel( M ); end;


[Mret_oret, v0, v0s, oorig] = deal({}); 
multiplicity = [];

for i = 1:noo
    if( orderinggiven );
        [v0new, v0snew, multiplicitynew] = leadingeigenvector_worker( tbuildproduct(M, oo{1,i}), epsilon, nomultipleeigenvector, complexeigenvector );
        Mret_oret = [Mret_oret repcell( oo{1,i}, size(v0new) )];
        oorig = [oorig repcell( oo{2,i}, size(v0new) )];
    else;
        [v0new, v0snew, multiplicitynew]=leadingeigenvector_worker( M{i}, epsilon, nomultipleeigenvector, complexeigenvector );
        % oorig is empty in this case
        Mret_oret = [Mret_oret repcell( M{i}, size(v0new) )]; end;
    if( anym(isnan([v0snew{:}])) || anym(abs([v0snew{:}])>1/(10*eps)) );
        vprintf( 'Orthogonal eigenplane and eigenvector (left and right eigenvector) occured.\n', 'imp',[1 verbose], 'once',1 );
        for j = 1:length( v0snew ); 
            v0snew{j}(:) = NaN; end; end;
    
    v0 = [v0 v0new]; 
    v0s = [v0s v0snew];
    multiplicity = [multiplicity repmat( multiplicitynew, size(v0new) )]; end;

if( complexeigenvector==2 && any(cellfun(@isreal,v0)) );
    idx = [];
    for i = 1:numel( v0 );
        if( ~isreal(v0{i}) );
            idx = [idx i]; end; end;
    v0(idx) = []; 
    v0s(idx) = []; 
    multiplicity(idx) = []; 
    Mret_oret(idx) = []; 
    if( ~isempty(oorig) ); 
        oorig(:,idx) = []; end; end;

vprintf( 'v0: \n%v\n', v0, 'imp',[2 verbose] );
vprintf( 'v0s: \n%v\n', v0s, 'imp',[2 verbose] );


end

function [ v0, v0s, multiplicity ] = leadingeigenvector_worker( M, epsequal, nomultipleleadingeigenvector, complexeigenvector );
% [M, v0, v0s] = leadingeigenvector_worker (M);      %returns all leading eigenvectors of the matrices in M
% computes v0 und v0s for one matrix M
    
    [v, d, vs] = eig(M); %each row is one vector
    %v=v+1i*randn(size(v));  vs=vs+1i*randn(size(vs)); %DEBUG
    if(~isvector(d)); 
        d = abs( diag(d) );  end; %we are only interested in the absolute value of the eigenvalues
    if(anym(isnan(v))); 
        error('NaNs occured during computation of eigenvectors.'); end;
    
    
    %normalize eigenvectors in order to compare them, and since we need them normalized anyway
    v = normalizematrix(v,'dirnorm',[1 2]); %normalize such that 2-norm of rows==1
    v = normalizematrix(v,'positive',1);
    
    %identify equal eigenvectors (which occur e.g. for the matrix [1 1; 0 1] );
    val = [real(v); imag(v)]; 
    val=exp(val); %workaround for strange behaviour in uniquetol, uniquetol([0 1-e100]) gives [0 1e-100] YY
    if(issym(val))
        [~,idx]=unique(val.','rows');
    else
        [~,idx]=uniquetol(val.','ByRows',true); end; %indices of pairwise lin.indep. eigenvectors. XX actually we want the mutually lin.ind. vectors
    
    
    v = v(:,idx); 
    vs = vs(:,idx); 
    d = d(idx);
    if( isempty(epsequal) );
        epsequal = min(cond(v)*10e-12,10e-8); end;
    idx = abs(d/max(d))>1-epsequal; %select highest eigenvalues
    
    complexidx = any(imag(v)~=0 | imag(vs)~=0,1)';
    if(complexeigenvector==0);
        %do not remove anything
    elseif(complexeigenvector==1 || complexeigenvector==2);
        %remove if there are also real eigenvectors
        if(~isequal(complexidx&idx,idx)); %do not remove everything
            idx(complexidx)=false; end;
        
    elseif(complexeigenvector==3 && ( anym(imag(v(:,idx))) || anym(imag(vs(:,idx)))) );
        v1 = real(v); vs1 = real(vs);
        v2 = imag(v); vs2 = imag(vs);
        v = [v1 v2];  vs = [vs1 vs2]; idx = [idx; idx];
        
    elseif( complexeigenvector==4 ); %remove all complex eigenvectors
        idx(complexidx) = false; end;
    
    if(nomultipleleadingeigenvector); 
        idx(find(idx,1)+1:end)=0; end; %if we only want to return one candidate, we choose the first
    
    %choose leading eigenvectors
    v0 = v(:,idx);
    v0s = vs(:,idx);    
    multiplicity = sum(idx); % determine multiplicity
    v0 = num2cell(v0,1);
    v0s = num2cell(v0s,1);
    
    if(any(idx));  
        v0s=normalizematrix(v0s,'dotprod',v0); end; %normalize *-eigenvectors by right eigenvectors
end


function oout = addcycle( oin )
    oout = {};
    for i = 1:size( oin, 2 )
        oout{1,end+1} = oin{1,i};
        oout{2,end} = oin{2,i};
        for j = 0:length( oin{1,i} )-1
            oout{1,end+1} = circshift( oin{1,i}, j );
            oout{2,end} = oin{2,i}; end; end;
end

function oout = addrepetition( oin, maxlength )
    if( isempty(oin) ); 
        oout = cell(2,0); 
        return; end;
    L = lcmm( unique(cellfun(@length,oin(1,:))) );
    if( L>maxlength ); 
        oout = oin; 
        return; end;
    oout = {};
    for i = 1:size( oin, 2 )
        o = oin{1,i};
        for j = 1:L
            if( j*length(o)<=L ); 
                if( isrow(o) ); 
                    oout{1,end+1} = repmat( o, [1,j] ).';
                else; 
                    oout{1,end+1} = repmat( o, [j,1] ); end;
                oout{2,end} = oin{2,i}; end; end; end;
end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   