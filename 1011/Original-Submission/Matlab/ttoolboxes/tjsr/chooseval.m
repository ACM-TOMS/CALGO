function [idxout,valout,barrier] = chooseval(valin, minmaxval, idxin )
% [ idxout, valout, barrier ] = chooseval(num, [minmaxval], [idxin] )
% Selects highest values of a vector. 
%
% Input:
%   valin         vector, values to choose from
%   minmaxval     optional, default=[1 numel(valin)] resp. [1 nnz(idxin)] if idxin is given and nonempty, 1x2 array, 
%                 Either empty, or
%                 a vector [minval, maxval], where minval=minimal number of things to choose, maxval= maximal number of things to choose, or
%                 a number minval, where minval=minimal number of things to choose.
%   idxin         optional, default=[], logical array of size num
%                 Defines from which valin we want to choose. All entries where idx==0 are not considered in the choice-process.
%
% Output:
%   idxout        logical index-vector of size valin, values chosen
%   valout        vector of length nnz(idxout), the selected values
%   barrier       scalar, the lower bound of the chosen values
%
% Note: Selects values at the end in the vector first.
%       undefined behaviour if minimal > maximal.
%       undefined behaviour for complex input.
%       +/- Infs are treated as numbers.
%       NaNs are treated as -inf.
%
%
% E.g.: [idxout, valout, barrier]=chooseval([Inf -Inf NaN 2 4 5],[2 3],[0 1 1 1 0])
%       [idxout, valout, barrier]=chooseval(randn(1,100))
%
% Written by tommsch, 2018

 %#ok<*ALIGN>
    
    WEIGHT=1.05; %The higher this value, the more high numbers are weighted more than small numbers.
        
    if(nargin<=2 || isempty(idxin)); 
        idxin=true(size(valin)); 
    end;
    
    if(nargin<=1 || isempty(minmaxval))
        minmaxval=[1 nnz(idxin)];
    end
    if(numel(minmaxval)==1); 
        minmaxval=[minmaxval nnz(idxin)]; end; 
    
    %cases when we have to choose everything
    if(~any(idxin) || nnz(idxin)<=minmaxval(1)); 
        idxout=logical(idxin); 
        valout=valin(idxout);
        barrier=-inf;
        return; end; 
    
    %cases when we have to choose nothing
    if(minmaxval(2)==0); 
        valout=[]; idxout=false(size(valin)); 
        barrier=-inf; 
        return; end;
    
    %make valin to row vector    
    if(iscolumn(valin)); 
        valin=valin.'; 
        transposeflag=1; 
    else; 
        transposeflag=0; 
    end;
    
    valin_select=valin(logical(idxin));
    idxin=find(idxin);
    
    [values,idxout]=sort(valin_select,'ascend');    
    nanidx = isnan(values); %put nans at the beginning
    values = [ values(nanidx) values(~nanidx)];
    idxout = [ idxout(nanidx) idxout(~nanidx)];
    
    val=diff(WEIGHT.^values);
    [val,idxdiffvalues]=sort(val,'descend');
    nanidx= isnan(val); %put nans at the end
    idxdiffvalues = [ idxdiffvalues(~nanidx) idxdiffvalues(nanidx)];
    
    val=length(valin_select)-idxdiffvalues; %number of values larger than the number on the resp. position +1
    idx=find(val>=minmaxval(1) & val<=minmaxval(2),1); %first index of the elements we choose
    idx=idxdiffvalues(idx(1));
    idx=idxout(idx+1:end);
    idx=idxin(idx);
    
    %convert to logical indexing
    idxout=false(size(valin));
    idxout(idx)=true;
    
    if(nargout>=2); 
        valout=valin(idxout); end;
    if(nargout>=3); 
        barrier=min(valout); end;
    
    if(transposeflag)
        idxout=idxout.';
        if(nargout>=2); 
            valout=valout.'; end;
        if(nargout>=3); 
            idxout=idxout.'; end;
    end    
end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 