function [ oo, pp, nrep ] = findperiod(varargin)
% [ oo, pp, nrep ]   = findperiod(oo, [options])
% [ ordering, nrep ] = findperiod(oo, [options])
% Searches for periodics in sequences.
%
% Input:
%   oo              the sequence to be examined
%
% Options:
%   'everywhere'    (This option may be changed) Searches for periodics, but does not consider irregularities at the end
%   'fuzzy'         (This option may be changed) Makes a fuzzy search.
%   'verbose'       Verbose level
%
% Output:
%   ordering       cell array of the form {oo,pp}
%   oo             sequence without nonperiodic stuff at the end. Equals o in general, except when 'everywhere' is set.
%   pp             period
%   nrep           number of repetitions of period. If nrep<1, pp is likely not to be a period
%
% E.g.: [oo, pp, nrep]=findperiod([1 2 1 2 1 2; 4 5 5 5 5 5])
%       [oo, pp, nrep]=findperiod([1 2 3 3 5 3 3 3 4])
%       [oo, pp, nrep]=findperiod([1 2 3 3 5 3 3 3 4],'fuzzy')
%       [oo, pp, nrep]=findperiod([1 2 3 3 5 3 3 3 4],'everywhere')
%
% Written by: tommsch, 2018

%#ok<*ALIGN>

oo=varargin{1};
dim=size(oo,1);
loo=size(oo,2);

if(loo==0);  
    pp=[]; oo=[]; nrep=0; 
    return; end;
if(loo==1); 
    pp=oo; oo=[]; nrep=0; 
    return; end;
if(dim==0); 
    pp=NaN; oo=[]; nrep=NaN; 
    return; end;

if(parsem('everywhere',varargin,0))
    verbose=parsem({'verbose','v'},varargin,1);
    if(verbose>=1); 
        fprintf([repmat(['.'],[1 loo]) '\n\n']); end; %#ok<NBRAK>
    [oo_i,pp_i,nrep_i]=deal(cell(1,loo));
    parfor i=1:loo
        if(verbose>=1); 
            fprintf('\b|\n'); end;
        [oo_i{i},pp_i{i},nrep_i{i}]=findperiod(oo(:,1:end-i+1)); %#ok<PFBNS>
        nrep_i{i}=nrep_i{i}-(i-1)/loo; %penalty for removing the tail
    end
    [~,idx]=max(cell2mat(nrep_i));
    oo=oo_i{idx};
    pp=pp_i{idx};
    nrep=nrep_i{idx};
    
elseif(parsem('fuzzy',varargin,0))
    oo_shifts=arrayfun(@(x) [nan(dim, x) oo(:,1:end-x)], 1:loo-1,'UniformOutput',false);
    what=cellfun(@(x) all(x==oo,1),oo_shifts,'UniformOutput',false);
    counts=cellfun(@nnz, what);
    [val,lpp]=max(counts); %period length
    if(val==0);
        pp=[]; nrep=0;
    else
        idx=strfind(what{lpp},ones(1,lpp)); idx=idx(1);
        pp=oo_shifts{lpp}(:,idx:idx+lpp-1);
        idx=arrayfun(@(x) strfind(oo(x,:),pp(x,:)), 1:dim,'UniformOutput',false); %strfind only support vectors. Thus we have to split the array up into vectors
        if(size(idx,2)==1)
            idx=idx{1};
        else
            idx=intersect(idx{:});
            idx=idx(1);
        end
        oo=oo(:,1:idx-1);
        nrep=val/lpp;
    end 
else
    oo_shifts=arrayfun(@(x) [nan(dim, x) oo(:,1:end-x)], 1:loo-1,'UniformOutput',false);
    counts=cellfun(@(x) find(~flip(all(x==oo,1),2),1),oo_shifts);
    [val,lpp]=max(counts); 
    if(val==1); 
        pp=[]; nrep=0; 
    else        
        pp=oo(:,end-lpp+1:end); %period
        nrep=0;

        %remove periods from oo
        while(true)
           if(size(oo,2)>=lpp && isequal(oo(:,end-lpp+1:end),pp))
               oo(:,end-lpp+1:end)=[];
               nrep=nrep+1;
           else; break;
           end
        end

        %put residues of periods from oo to period
        while(true)
            if(~isempty(oo) && isequal(oo(:,end),pp(:,end)))
                oo(:,end)=[];
                pp=circshift(pp,[0 1]);
                nrep=nrep+1/lpp;
            else
                nrep=nrep-1; break; %do not count the period itself as occurence
            end
        end
    end
end

%Post-processing
if(nargout<=2)
    oo={oo,pp};
    pp=nrep;
end
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 