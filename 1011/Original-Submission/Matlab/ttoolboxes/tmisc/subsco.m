function ret = subsco(varargin)
% data = subsco(A,co,[options])
% Anew = subsco(A,co,data,[options])
% Indexing of matrices by coordinate-vectors.
%
% 1) data = subsco(A,co, [options])
%    Reads the data of array A at coordinates co and returns the data
%       Input:
%           A       Array to be read from
%           co      coordinates where to read from
%       Output:
%           ret     the data as a row vector
%       
%       E.g.: subsco([1 2 3; 4 5 6],[1 2; 2 3]')
%
% 2) Anew = subsrefco(A,co,data, [options])
%    Sets the values of data(i) at indicex defined by co(:,i)
%       Input:
%           A       Array to be written
%           co      where to write
%       Output:
%           A       Array with new values
%
%       E.g.: subsco([1 2 3; 4 5 6],[1 2; 2 3]',[-1 -1])
%
% Options:
%   'save'      checks for out of bounds error
%               also checks for multiple indices if subsco(A,co,data) is called and issues an error
%
% See also: sub2ind, subsref, subsasgn
%
% Written by: tommsch, 2018

A=varargin{1};
SIZE=size(A);
co=varargin{2};
%dim=max(ndimsm(A),ndimsm(co));
[saveflag,varargin]=parsem('save',varargin);
setdata=0; 
if(size(varargin,2)==3); 
    setdata=1; 
    data=varargin{3}; 
end;

% if(any(co.'>sizem(A,[],dim))); %XX das geht nicht wenn mehr als ein Element in co ist
%     if(setdata==0); error('Index out of bound.'); end;
%     val=num2cell(co.');
%     A(val{:})=data;
%     ret=A;
%     return;
% end

if(saveflag);
    idx=zeros(1,size(co,2));
    for i=1:size(co,1)
        idx=idx | co(i,:)>SIZE(i) | co(i,:) <= 0;
    end
    co(:,idx)=[];
    if(setdata)
        data(idx)=[];
        if(size(unique(co','rows'),1)~=size(co,2));
            error('subsco: Multiple indices given.'); end
    end
end

    
co=num2cell(co,2);

L = sub2ind(SIZE, co{:}); % Obtain the linear indices implied by treating co=[Row, Col, ... ] as subscript matrices

if(~setdata); %read values and return them 
    ret = A(L); 
else; %set values and return the new matrix
    A(L)=data;
    ret=A;
end

end


