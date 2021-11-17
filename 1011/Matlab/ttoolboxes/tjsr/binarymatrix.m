function T = binarymatrix(dim, N, k, flag)
% T = binarymatrix(dim, N, k, [flag])
% Returns the k^th set of binary matrices, ordered linearly
% Returns the empty set, for certain values of k
%
% Input:
%   dim         dimension
%   N           number of matrices in the set
%   k           counter
%   flag        (experimental) default=0, If set and there exists some kk<k such that binarymatrix(kk) would have the same JSR, 
%               then the empty may be returned.
% 
% Output:
%   T           1xN cell array of dimxdim matrices
%               Returns {} if a smaller k leads to the same set
%
% E.g.: vdisp(binarymatrix(2, 2, 2^5+2^3+2^2+2^1))
%       vdisp(binarymatrix(2, 2, 2^5))
%
% See also: tgallery
%
% Written by: tommsch, 2018
 
if( nargin==3 ); 
    flag = 0; end;

%generate matrix
Alin = dec2bin(k)' - '0';
Alin = [zeros(dim^2*N-length(Alin),1); Alin];
if(length(Alin)>dim^2*N); 
    T={}; 
    return; end;
Abox = reshape(Alin,dim,dim,N);

%make simple test, if this matrix corresponds to a simpler matrix
if(any(diff(summ(Abox,[1,2]))<0)); 
    T={}; 
    return; end; %each row corresponds to one matrix

T=cell(1,N);
for i=N:-1:1
    T{i}=Abox(:,:,i); end;

if(flag)
    %make more simple tests
    Atest=Abox;     
    for i=1:N; 
        Atest(:,:,i)=Abox(:,:,i).'; end;

    val=reshape(Atest,[],1); 
    if(bin2dec(num2str(Alin'))>bin2dec(num2str(val'))); 
        T={}; 
        return; end; %test for transpose

    Atest=Abox; 
    P=perms(1:dim); 
    for j=1:size(P,1); 
        for i=1:N
            Atest(:,:,i)=Abox(P(j,:),P(j,:),i); end;
        val=reshape(Atest,[],1); 
        if(bin2dec(num2str(Alin'))>bin2dec(num2str(val'))); 
            T={}; 
            return; end; end; end; %test for permutation

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 