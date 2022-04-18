function [type, idx_choose,  numselected, VV, VVidx, nVVbig] = tjsr_selectvertex(type);
% [type, idxchoose,  numselected, numunselected] = tjsr_selectvertex(type);
% This function belongs to tjsr!
% Selects vertices for which the norm gets computed
% Selects the vertices of the polytope which are used to compute the norms
%
% Input: 
%   type            type-struct from tjsr
%
% Output:
%   type            type-struct for tjsr
%   idxchoose       linear index (array) of selected vertices
%   numselected     number of selected vertices/number of unselected vertices
%   VV              vertices of chosen polytope which is used to compute the norm
%   nVVbig          Number of vertices in the boundary of the polytope
%
% Written by tommsch, 2018

% XX could get rewritten


% Select vertices for which the norm gets computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idxnan = isnan([type.cyclictree.norm{:}]); %only select vertices which are not computed yet

if(type.opt.naturalselectiontype>0);
    minchoose = ceil(abs(type.opt.naturalselection)/type.counter.nummatrix)+1;
else;
    minchoose = ceil(abs(type.opt.naturalselection)); end;
maxchoose = minchoose*4;

numall = nnz(idxnan);

switch abs(type.opt.naturalselectiontype)
    case inf
         switch mod(type.counter.iteration,4);
             case {1,2,3}; 
                 val = [type.cyclictree.normest{:} ];
             case 0; 
                val = [type.cyclictree.normparent{:} ]; end
    case 1; 
        val = [ type.cyclictree.normest{:} ];
    case 2; 
        val = [ type.cyclictree.normparent{:} ];
    case 3; 
        val = [ type.cyclictree.rho{:} ];
        
    %values from 100 to 999 denote selection rules which are not suitable for production
    case 100; val = -[type.cyclictree.rho{:}]; end; %DEBUG, the algorithm does not terminate with this selection for certain matrices: B=[1 1; 0 1]; tjsr({1/2*B,1/2*B',B,B'},'ordering',{[1],[2]},'naturalselectiontype','-rho','fastnorm',0,'testspectralradius',-inf,'plot','polytope')
  
[ idx_choose, ~, barrier ] = chooseval( val, [minchoose, maxchoose] , idxnan ); 

if(type.opt.verbose>4);
    figure(1148242414); clf; 
    semilogy(sort(val),'.'); hold on;
    plot(barrier*ones(1,size(val,2)),'r-'); end;

%select all children of each parent which has at least one selected child.
if(type.opt.naturalselectiontype>0 && minchoose<inf)
    LL = cellfun(@sum,type.cyclictree.L);
    valcum = [0 cumsum(LL)];
    for i = 1:type.counter.numordering
        idx = idx_choose(valcum(i)+1:valcum(i+1));
        parentidx = unique(type.cyclictree.parent{i}(idx));
        idx_choose(valcum(i)+1:valcum(i+1)) = ismember(type.cyclictree.parent{i},parentidx) & idxnan(valcum(i)+1:valcum(i+1)); end; end;

numselected(1) = nnz(idx_choose);
numselected(2) = numall-numselected;

% Select vertices which are used to compute the norms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
if( type.opt.simplepolytope && length(type.JSR)>1 );
    distval=1+min(type.opt.simplepolytope, (type.JSR(2)/type.JSR(1)-1)/1000)-type.opt.epslinprog; %vertices nearer than distval next to each other are disregarded
else
    distval=0; end;
    
if(type.opt.naturalselectiontype>0 && minchoose<inf)
    %indices of vertices whose norm will be computed or is already computed
    idx_good = ~isnan([type.cyclictree.norm{:}]) | idx_choose;
    idx_good = mat2cell(idx_good,1, cellfun(@sum,type.cyclictree.L));

    
    idx_i=cell(1,type.counter.numordering);
    for i=1:type.counter.numordering
        % only select vertices whose childrens norms is computed or gets computed
        idx_i{i}=type.cyclictree.parent{i}(idx_good{i}==0); %linear index of vertices whose childrens norm is not known and is not computed
        val=ones(1,sum(type.cyclictree.L{i}));
        val(idx_i{i})=0;
        idx_i{i}=val; %change linear index to logical index
        
        % only select vertices which are outside or at the border of the polytope.-
        idx_i{i} = idx_i{i} & ( type.cyclictree.norm{i} > 1-type.opt.epspolytope );
        
        % only select vertices which are far away from other vertices 
        idx_i{i} = idx_i{i} & (type.cyclictree.norm{i} > distval); end;
    VVidx=[idx_i{:}];
else
    VVidx = [type.cyclictree.norm{:}] > 1-type.opt.epspolytope;
    VVidx = VVidx & [type.cyclictree.norm{:}] > distval;  end; % only select vertices which are far away from other vertices 

    
% make matrix of vertices

VV = [type.cyclictree.V{:}]; %get all vertices
VV = VV(:,VVidx); % remove all unselected vertices
VV = removezero(VV,2); % remove empty columns

%number of all vertices
nVVbig=nnz([type.cyclictree.norm{:}]>1-type.opt.epspolytope);
    
end





function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 