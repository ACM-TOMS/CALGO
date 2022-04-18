function [cycles]=removecombination(cycles,varargin);
% [ oclass ] = removecombination( cycles , [options] );
% Constructs a minimal set of cycles.
% Longer explanation:   Assume the product of matrices AB and BAB have the same eigenvector v1 to the eigenvalue 1.
%                       Then ABv = v and ABAv=v. This implies that Bv=v, thus Av=v, thus BAv=v. 
%                       The function returns: A, B. This is a minimal set which can construct AB and ABA
% Input:
%   cycles      (cell array of vectors) vectors of numbers which define cycles for the same eigenvector
%
% Options:
%   'verbose',val     Verbose level
%   'add'             Returns long cycles instead of short cycles
%
% Output:
%   cycles      shorter list of numbers which define a cycle
%
% E.g.: removecombination({[1 1 2],[1]})
%
% See also: addcombination
%
% Written by: tommsch, 2018
    if(parsem('add',varargin));
        cycles = addcombination(cycles);
        return; end;
    
    while(true)
        repeat=0;
        ncycles=length(cycles); %originally length_m
        if(ncycles==1); 
            return; end; %trivial case
        
        for i=1:ncycles-1;
            ooi=cycles{i};
            for j=i+1:ncycles;
                ooj=cycles{j};
                [vali, valj, shorten] = comparebegin2(ooi,ooj);
                if(shorten); 
                    ooi=vali;
                    cycles{i}=vali;
                    cycles{j}=valj;
                    repeat=1; 
                end;
            end
        end
        cycles(cellfun(@isempty,cycles))=[];
        %vdisp(cycles)
        if(~repeat); 
            break; end;
    end

end

function [cycles]=addcombination(cycles);
    if(length(cycles)==1); 
        return; end; %trivial case %length_m

    testelement=2;
    while(true)
        oo=cycles{testelement};
        for j=1:testelement-1;
            [o, or] = comparebegin(oo,cycles{j});
            if(isempty(o)); 
                continue; end;
            if(searchincellarray(or,cycles,1)); 
                continue; 
            else
                cycles{end+1}=or; %#ok<AGROW>          
            end;

        end
        testelement=testelement+1;
        if(testelement>length(cycles)); 
            break; end; %length_m
    end
end

function [o, or] = comparebegin(o1,o2)
%  compares if one of the orderings is the same as the beginning of the other ordering
%   o       the part at the beginning which is the same (thus equals either o1 or o2)
%           o is empty, if o1 and o2 are not compareable.
%   or      the rest of the ordering, where o is removed.
% Eg [o, or] = comparebegin([1 2 3],[1 2]) % o = [1 2]; or = [ 3 ];
    if(length(o1)>length(o2));  %length_m
        [o1,o2]=swap(o1,o2); 
    end;
    if(isequal(o1,o2(1:length(o1)))); %length_m
        o = o1;
        or = o2(length(o1)+1:end); %length_m
    else
        o=[];
        or=[];
    end
end


function [o1, o2, shorten] = comparebegin2(o1,o2)
%  compares if one of the orderings is the same as the beginning of the other ordering
%   o       the part at the beginning which is the same (thus equals either o1 or o2)
%           o is empty, if o1 and o2 are not compareable.
%   or      the rest of the ordering, where o is removed.
%   o1 must be shorter than o2
% Eg [o, or] = comparebegin([1 2 3],[1 2]) % o = [1 2]; or = [ 3 ];
    swapflag=0;
    shorten=0;
    
    lo1=length(o1);
    lo2=length(o2);
    
    if(lo1==0 || lo2==0); 
        return; end;
    
    if(lo1>lo2); 
        [o1,o2]=swap(o1,o2); 
        swapflag=1; 
    end;
    
    if(isequal(o1,o2(1:length(o1))));        %length_m
        o2 = o2(length(o1)+1:end); %length_m
        shorten=1;
    end
    if(swapflag);
        [o1,o2]=swap(o1,o2); 
    end
end



function [y,x]=swap(x,y)
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   