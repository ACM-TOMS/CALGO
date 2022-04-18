function D = diag(C,k)
% D = diag(C,k)
% Same behaviour as diag for matrices
%
% E.g.: diag(diag(num2cell(randn(4)),2),2)
%
% See also: diag, cell/triu
%
% Written by tommsch, 2018

 %#ok<*ALIGN>

if(nargin==1); 
    k=0; end;
if(k==inf || k==-inf); 
    if(~isvector(C));
        D=cell(0,1);
        return; 
    else
        error('Kth diagonal input must be finite for the vector input form of diag().'); end;
end;
sze=size(C);

if(any(sze==1))
    sze=max(sze);
    D=cell(sze+abs(k));
    if(k>=0)
        for i=1:sze;
            D{i,i+k}=C{i};
        end
    else
        for i=1:sze
            D{i-k,i}=C{i};
        end
    end
else
    if(any(sze)==0); D={}; 
        return; end;
    if(k>=0); 
        l=min(min(sze),sze(2)-k);
        if(l<=0); 
            D=cell(0,1);
        else
            D=cell(l,1);
            for i=1:l
                D{i}=C{i,i+k};
            end
        end;   
    else;     
        l=min(min(sze),sze(1)+k);
        if(l<=0); 
            D=cell(0,1);
        else
            D=cell(l,1);
            for i=1:l
                D{i}=C{i-k,i};
            end
        end;
    end;
    
end
   

end