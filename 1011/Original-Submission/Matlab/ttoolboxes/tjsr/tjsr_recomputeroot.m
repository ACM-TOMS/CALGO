function type = tjsr_recomputeroot(type);   
% type = tjsr_recomputeroot(type); 
% This function belongs to tjsr!
% Removes vertices from the cyclic root, which are already inside of the polytope.
% Uses 10*abs(epslinprog) as epsilon.
%
% Input: 
%   type    type-struct from tjsr containing
%               type.cyclictree.----
% Output: 
%   type    type-struct for tjsr
%
% Written by tommsch, 2018

    type.info.infotext=vprintf('Remove unnecessary vertices from root. \n','imp',[2,type.opt.verbose],'str',type.info.infotext); 
    %recompute all pn-estimates
    VV=[type.cyclictree.V{:}]; %get all vertices
    
    %choose all vertices
    new_pn=computepolytopenorm(VV, VV, type.info.algorithm, type.opt.numcore, type.opt.epslinprog, min(type.opt.verbose,1), type.opt.solver );
    idx_pn=new_pn<1-10*abs(type.opt.epslinprog); % be on the safe side
    
    idx=type.cyclictree.smpflag==0;
    val=mat2cell(idx_pn,1,[type.cyclictree.L{:}]);
    if(any([val{idx}]));
        type.info.infotext=vprintf('A vertex from the cyclic-root of a smp-candidate is inside the convex hull.\n','imp',[2,type.opt.verbose],'str',type.info.infotext);
    end
    
    if(nnz(idx_pn)); 
        type.info.infotext=vprintf('Removed %i vertices. ',nnz(idx_pn),'imp',[1,type.opt.verbose],'str',type.info.infotext); 
        type.cyclictree.norm = savetocellarray(new_pn(idx_pn), idx_pn, type.cyclictree.norm);
    end;
    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 