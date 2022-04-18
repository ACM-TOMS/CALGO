function [ X ] = setplus( varargin )
% X  = setplus( A1, A2, A3, [options] )
% Element-wise addition of vectors.
% setplus(A,B,C)={x=a+b+c:a\in A, b\in B, c\in C}
%
% Input: 
%   An              arrays of compatible size
%
% Options:
%   'rows'          setplus is applied row-wise (would be Matlab standard). Default behaviour is column-wise!!
%   'stable'        Stable-sorts the output.
%   'nounique'      Does not sort/unique output
%
% Output:
%   X               the Sum
%
% E.g.: setplus([1 0; 1 1]',[10 10; 20 20; 30 30]')
%
% Written by: tommsch, 2017


[rows,varargin]=parsem('rows',varargin); 
[stable,varargin]=parsem('stable',varargin);
[nounique,varargin]=parsem('nounique',varargin);

s=size(varargin,2);
X=varargin{s};

if(rows)
    for i=s-1:-1:1
            X=setplus_worker_rows(varargin{i},X,stable,nounique);
    end
else
    for i=s-1:-1:1
            X=setplus_worker_column(varargin{i},X,stable,nounique);
    end
end

end

function [ C ] = setplus_worker_column( A,B,stable,nounique )
    C=repmat(A,1,size(B,2))+reshape(repmat(B,size(A,2),1),size(A,1),[]);
    if(nounique); 
        return; end;
    if(stable);
        C=unique(C.','rows','stable').';
    else
        C=unique(C.','rows').';
    end
end

function [ Ct ] = setplus_worker_rows( At,Bt,stable,nounique )
    Ct=kron(ones(size(Bt,1),1),At)+kron(Bt,ones(size(At,1),1));
    if(nounique); 
        return; end;
    if(stable)
        Ct=unique(Ct,'rows','stable');
    else
        Ct=unique(Ct,'rows');
    end
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   