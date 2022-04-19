function D = dec2basem( v, M ,nin)
% D = dec2basem( v, M , [nin])
% Computes M-ary representations.
% Returns the representation of v as an array in base M.  
%
% Input:
%   v       dim x 1,non-negative integers, the number to converted
%   M       dim x dim,  dilation matrix
%   nin     (Optional) produces a representation with at least nin digits.
% 
% Output:
%   D       dim x N array of indices
%           If the representation does not stop at length nin, NaN is appended at the end
% 
% E.g.: dec2basem([101243; 23],[1 2; -1 2])
%
% Written by tommsch, 2018

% XX Add option which defines digits to be used

%#ok<*ALIGN>

if(nargin<=2); 
    nin=1; end;
r=norm(M,'inf');
i=0;
dim=size(M,2); 
if(isvector(M)); 
    dim=1; end;
if(abs(det(M))<=1); 
    error('''M'' must be expanding, i.e. abs(det(M))>1.'); end;
zero=zeros(dim,1);
D=zeros(dim,nin);
MAX=norm(v,'inf')^10;
while(~isequal(v,zero));
    if(r^size(D,2)>MAX && i>nin)  
        D=[ones(dim,1)*NaN D]; %#ok<AGROW>
        break; end;
    i=i+1;
    d=v-M*floor(M\v);
    v=M\(v-d);
    D(:,i)=d;
end
D=round(D);
D=fliplr(D);
  
end
