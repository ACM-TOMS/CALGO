function [ P ] = reducelength(P, flag)
% P = reducelength(P, flag)
% Takes a pattern and removes repetitions and cycles it such that it has smallest lexicographic value.
%
% Input:
%   P       row/column vectors, OR
%           cell array of row/column vectors
%   flag    what to do:
%               0   (default) removes periodics and cycles
%               1   removes only periodics
%               2   cycles only
%
% Output:
%   P     (in the same format as the Input)
%
% E.g.:   vdisp(reducelength({[2 1 2 2 1 2],[3 1]}))
%         vdisp(reducelength({[2 1 2 2 1 2],[3 1]},1))
%
% See also: lexicographic 
%
% Written by: tommsch, 2018

%#ok<*ALIGN>

if(nargin==1); 
    flag=0; end;

if(iscell(P));
    for i=1:size(P,2);
        if(isempty(P{i})); 
            continue; end;
        P{i}=reducelength_worker(P{i},flag);
    end
else
    if(isempty(P));
        P=[];
    else
        P=reducelength_worker(P,flag);
    end
end

end

function [ P ] = reducelength_worker(P,flag)
P(P==0)=[]; %remove all zeros in P
P=double(P);

n=numel(P);
if(size(P,1)>1); 
    isrowvector=1; P=P'; 
else; 
    isrowvector=0; end;
if(isempty(P)); 
    P=zeros(0,1); 
    return; end;

if(flag==0 || flag==1)
    
    %find periodics
    %%%%%%%%%%%%%%%%%%%%%
    f=factor(n);
    
    %find all divisors of n
    f=nchoosek([ones(size(f)) f], length(f));
    f=cumprod(f,2);
    f=f(:,end)';
    f=unique(f);
    f=f(1:end-1);

    %find periodics
    for i=1:length(f)
        canshorten=1;
        test=P(1:f(i));
        for j=2:numel(P)/f(i)
            if sum(P((j-1)*f(i)+1:(j)*f(i))~=test)
                canshorten=0;
                break; end;
        end
        if canshorten==1; 
            P=test;
            break; end;
    end
end

if(flag==0 || flag==2);
    %cycle lexicographic
    MAX=ones(1,size(P,2))*max(P);
    MAX=MAX.^(length(P)-1:-1:0);
    Pnew=P;
    val=sum(P.*MAX);
    for i=1:length(P)-1;
        Ptest=circshift(P,[0 i]);
        valnew=sum(Ptest.*MAX);
        if(valnew<val);
            val=sum(Ptest.*MAX);
            Pnew=Ptest; end;
    end
    P=Pnew;
end    


if(isrowvector); 
    P=P'; end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 