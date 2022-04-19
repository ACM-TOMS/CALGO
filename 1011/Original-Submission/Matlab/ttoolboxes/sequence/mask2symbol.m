function ret=mask2symbol(varargin)
% [ ret ] = mask2symbol( a, [options])
% Computes the symbol of a mask.
%
% Input: 
%       a           mask
%
% Options: 
%       'dim',val   specifies the dimension, integer
%       'amin',val  index of left upper entry, column vector. 'amin' must have the same format as 'a'
%       'var',val   Controls the name of the variables. Default: z1,z2,z3,...
%                   Format string:  char (ex: 'x')                      Variables: char,char+1,char+2, (ex: x,y,z)
%                                   char+num  (ex:'b3')                 Variables: charnum,charnum+1,charnum+2 (ex: b3,b4,b5,b6,...)
%                                   cellarray (ex:{'x','a23','fa')      Variables: cellarray{1},cellarray{2},... (ex: x,a23,fa)
%                                         If one puts values in the cell array, then the symbol is evaluated at this point.
%
% Output:
%       ret         symbol of mask
%
% Note: This function is much slower than symbol2mask()
%
% E.g.: mask2symbol([1 2 3; 4 5 6])
%       mask2symbol([1 2 3],'var','x')
%
% See also: symbol2mask
%
% Written by: tommsch, 2017

a=varargin{1}; varargin(1)=[];
[dim,varargin]=parsem('dim',varargin,   max(nestedcellfun(@(x) ndimsm(x),a))  );
[amin,varargin]=parsem('amin',varargin, unflatten(zeros(dim,1),a,'assign') );
[var,varargin]=parsem('var',varargin,'z1');

parsem(varargin,'test');

%create variables    
z = sym(zeros(1, dim)); 
idx=1;
if(iscell(var))
    for k=1:dim;
        z(idx) = sym(sprintf('%c', var{k})); 
        idx=idx+1;
    end; 
elseif(sizem(var,2)==1);
    for k=double(var(1)):dim+var(1)-1; 
        z(idx) = sym(sprintf('%c', k)); 
        idx=idx+1;
    end; 
elseif(sizem(var,2)==2);
   for k=str2double(var(2)):dim+str2double(var(2))-1; 
       z(idx) = sym(sprintf([num2str(var(1)) '%d'], k)); 
       idx=idx+1;
   end; 
end

SZE=size(a);
ret=sym(0);
IDX=cell(1,dim);
for j=1:numel(a);
    [IDX{:}]=ind2sub(SZE,j);
    IDX2=cell2mat(IDX);
    ret=ret+a(j)*prod(z.^(IDX2+amin'-1));
end  

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 