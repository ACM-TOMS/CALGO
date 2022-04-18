function type = tjsr_setproblemtype( type, v0properties ); 
% type = tjsr_setproblemtype( type, v0properties ); 
% This function belongs to tjsr!
% sets the value of type.info.algorithm
% Input:
%   type
%   v0properties        return struct from identifymatrix() of the set of vertices from the root
%
% Written by tommsch, 2018

    if(type.info.matrixtype.nonneg); 
        problemtype='nonneg'; %if the matrices are non-negative, by Perron-Frobenius, the leading eigenvalue is real and positive.
    elseif(~type.info.matrixtype.real); 
        problemtype='complex';
    elseif(v0properties.real); 
        problemtype='real';
    else; 
        error('Unkown error. ''type.info.problemtype.matrix'' or ''algorithm'' is invalid');end;
    
    if(isempty(type.info.algorithm) && strcmp(problemtype,'nonneg') || ~isempty(type.info.algorithm) && type.info.algorithm==TJSR_CONEFUNCT); 
        type.info.infotext=vprintf('Case (P).\n','imp',[2,type.opt.verbose],'str',type.info.infotext);
        type.info.algorithm=TJSR_CONEFUNCT;
    elseif(isempty(type.info.algorithm) && strcmp(problemtype,'real') || ~isempty(type.info.algorithm) && type.info.algorithm==TJSR_MINKFUNCT); 
        type.info.infotext=vprintf('Case (R).\n','imp',[2,type.opt.verbose],'str',type.info.infotext);
        type.info.algorithm=TJSR_MINKFUNCT;
    elseif(isempty(type.info.algorithm) && strcmp(problemtype,'complex') || ~isempty(type.info.algorithm) && type.info.algorithm==TJSR_COMPLEXFUNCT); 
        type.info.infotext=vprintf('Case (C).\n','imp',[2,type.opt.verbose],'str',type.info.infotext);
        type.info.infotext=vprintf('If this set of matrices is of practical interest, please contact the author ( tommsch@gmx.at ),\nbecause case (C) never occured so far for matrices of practical interest.\n','imp',[1,type.opt.verbose],'str',type.info.infotext,'cpr','red','npr');
        type.info.errortext=vprintf('If this set of matrices is of practical interest, please contact the author ( tommsch@gmx.at ),\nbecause case (C) never occured so far for matrices of practical interest.\n','str',type.info.errortext);
        type.info.algorithm=TJSR_COMPLEXFUNCT;
    else
        error('Unkown error. ''type.info.problemtype.matrix'' or ''algorithm'' is invalid'); end;
    
end