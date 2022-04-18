function [ type ] = tjsr_computemaxnormval( type, idx_VV );
% maxnormval = tjsr_computemaxnormval( type );
% This function belongs to tjsr!
% computes the value maxnormval, respects the value of delta and epspolyte
%
% Input:
%   type        type struct from tjsr
%   VVidx       logical indices of chosen vertices for polytope used in norm computation
%
% Output:
%   maxnormval
%
% Written by tommsch, 2019

% Changelog: 2019-03-30     added test whether polytope is ok for computation of intermediate bound

    
    %get all indices where status==0 and norm is computed 
    
    idx_freshparent = tjsr_getfreshparent( type );
    normval=[type.cyclictree.norm{:}];
    maxnormval=max(normval(idx_freshparent));
    
    %check if used polytope yields a new bound for the JSR
    idx_VV = mat2cell(idx_VV,1,cellfun(@sum,type.cyclictree.L));
    for i=1:type.counter.numordering
        idx=type.cyclictree.parent{i}(isnan(type.cyclictree.norm{i})); %linear index of vertices whose childrens norm is not computed
        val=zeros(1,sum(type.cyclictree.L{i}));
        val(idx)=1; 
        idx=val; %linear index changed to logical index
        if(any(idx & idx_VV{i})); 
            maxnormval=inf; 
            if(type.opt.naturalselectiontype>0)
                type.info.infotext = vprintf('Intermediate bound could not be computed. This should not happen and indicates a bug. Please inform the author.\n','cpr','err','str',type.info.infotext);
                type.info.errortext = vprintf('Intermediate bound could not be computed. This should not happen and indicates a bug. Please inform the author.\n','str',type.info.errortext,'npr'); end
            break; end; end; %this is actually not allowed to happen with default options
        
    %consider epspolytope and delta    
    if(isempty(maxnormval)); 
        maxnormval=1; end; %compute maxnormval
    maxnormval=max(maxnormval,1);
    maxnormval=maxnormval*(max(1,1-type.opt.epspolytope))/min(1,type.opt.delta);
    if(maxnormval<type.cyclictree.normlvl(end))
        type.cyclictree.normlvl(end+1)=maxnormval;
    else
        type.cyclictree.normlvl(end+1)=type.cyclictree.normlvl(end); end;
    
    if(type.cyclictree.normlvl(end)*type.lambda<type.JSR(2));
        type.JSR(2)=type.cyclictree.normlvl(end)*type.lambda; end;

end

function idx_freshchild = tjsr_getfreshparent( type )
% This function belongs to tjsr!
% returns linear indices of vertices which are not inside the polytope and whose children are not all computed yet
% 
% Input:
%   type        the type struct from tjsr
%   cellflag    returns cell array of indices for each tree
%
% Written by tommsch, 2019


idx_freshchild=cell(1,type.counter.numordering);
for i=1:type.counter.numordering
    idx_norm = type.cyclictree.norm{i}>1-type.opt.epslinprog; %vertices which are outside or at the border
    idx_nan = isnan(type.cyclictree.norm{i}); %vertices without norm
    idx_status = type.cyclictree.status{i}==0; %vertices which are children
    val = type.cyclictree.parent{i}(idx_nan);
    idx_partlyfreshparent = false(size(idx_nan));
    idx_partlyfreshparent(val) = true; %vertices whose children norms are not fully computed
    
    idx_freshchild{i} = idx_norm & (idx_partlyfreshparent | idx_status);
end

idx_freshchild = [idx_freshchild{:}];


end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   