function C = mtimes(A,B)
% mtimes(A,B)
% Scalar * Cell : Elementwise
% Cell * Scalar : Elementwise
% Cell * Cell   : Matrix - Matrix Multiplication 
% Cell * Matrix : Matrix - Matrix Multiplication
% Matrix * Cell : Matrix - Matrix Multiplication
% Vector * Cell : Not defined
%
% E.g.: vdisp(2*{[1 2], [2;3]})
%       
% See also:cell/times, cell/mtimes
%
% Written by: tommsch, 2018

% XX {2}*{2 3;4 5} XX yields wrong result

error("mtimes for cells not yet implemented.");

if(isscalar(A) && iscell(B));
    C=cellfun(@(x) A*x,B,'UniformOutput',false);
elseif(iscell(A) && isscalar(B));
    C=cellfun(@(x) x*B,A,'UniformOutput',false);
elseif(iscell(A) && iscell(B));
    %error('Matrix - Matrix multiplication is imlemented wrong.');
    %if(size(A,2)~=size(B,1)); error('A and B have not compatible size.'); end;
    %B=B';
    L=size(A,1); M=size(A,2); N=size(B,2);
    C=cell(L,N);
    for l=1:L
        for n=1:N
            C{l,n}=0;
            for m=1:M
                if(numel(A{l,m})==1 && numel(B{m,n})==1);
                    C{l,n}=C{l,n}+A{l,m}*B{m,n};
                else
                    C{l,n}=addsequence(C{l,n},A{l,m}*B{m,n});
                end
                
            end            
        end
    end
elseif(iscell(A) && ismatrix(B))
    C=A*num2cell(B);
elseif(ismatrix(A) && iscell(B));
    C=num2cell(A)*B;
else
    error('This cell multiplication is not defined.');
end

end