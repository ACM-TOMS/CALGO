function [d, dmin]=upsamplem(c, idx, M, fillvalue)
% d = upsamplem(c, idx, M, [fillvalue] )
% Multivariate upsampling.
% For \alpha\in\ZZ^s: upsamplen(c,M)(\alpha)=c(\beta) for \alpha=M\beta; and fillvalue otherwise
%
% Input:
%   c           (multidimensional array) the sequence to be upsampled
%   idx         (empty OR dim x 1-vector) index of first entry in c. If cmin is empty, then cmin=0.
%   M           (dim x dim-matrix) the dilation which defines to upsample
%   fillvalue   (scalar) the value for the "holes", default=0.
%
% Output:
%   d           the upsampled array
%   dmin        index of the first entry in d
% 
% E.g.: [d,dmin]=upsamplem([1 2 3; 4 5 6], [0;0], [2 1; 0 -2],Inf)
%
% Written by: tommsch, 2018

%#ok<*ALIGN> 

if(nargin==3); 
    fillvalue=0; end;
if(isempty(idx)); 
    idx=zeros(size(c,1),1); end;

    %error checks
    if(nargin<3);
        error('d = upsamplem(c, idx, M, [fillvalue] ). Argument(s) missing.'); end;
    if(~isscalar(fillvalue)); 
        error('d = upsamplem(c, idx, M, [fillvalue] ). Optional fourth entry must be a scalar, default=0'); end;
    if(~ismatrix(M) || ~issquare(M));         
        error('d = upsamplem(c, idx, M, [fillvalue] ). Third entry must be a square matrix.'); end;
    if(~iscolumn(idx));       
        error('d = upsamplem(c, idx, M, [fillvalue] ). Second entry must be a column vector or empty.'); end;
    
    dim=size(M,2); %the dimension in which we are working
    if(dim<ndimsm(c));      
        error('Dimension of sequence is larger than dimension of upsampling matrix.'); end;
    
    %very easy case: empty sequence
    if(isempty(c)); 
        d=[]; dmin=idx; return; end;
    
    try %#ok<TRYNC> %If Toolbox is not installed
        if(dim==1);%if dim==1, we use standard upsampling 
            d=upsample(c,abs(M)); %if M<0, we upsample in the standard way ... and reverse the sequence later
            d=d(1:end-abs(M)+1); %remove zeros
            if(M<0); 
                d(:)=d(end:-1:1); 
                dmin=(length(c)-1+idx)*M;
            else
                dmin=idx*M;
            end  %we reverse the sequence, since M<0
            return; end;
    end
        
    %Make all corner points %the coordinates of the corners points
    cornersidx=mixvector(1:2,dim); %possible combinations of corner-points
    cornersidx=setplus(cornersidx,((0:dim-1)*2)'); % add offset to adress correct numbers in corners, since I use reshape for val
    val=[zeros(1,dim); sizem(c,[],dim)-1]; %each column is the first and last idx for the corresponding dimension
    val=reshape(val,1,[]);  %the possible permutation of corner-points
    corners=val(cornersidx); %coordinates of the corners
    
    valreal=setplus([zeros(1,dim); sizem(c,[],dim)-1], idx','rows','nounique'); %each column is the first and last idx for the corresponding dimension
    valreal=reshape(valreal,1,[]);  %the possible permutation of corner-points
    cornersreal=valreal(cornersidx);
    
    %Compute where they are mapped to
    CORNERS=M*corners; %the coordinates of the corner points
    CORNERSREAL=M*cornersreal; %the coordinates of the corner points
    
    %Get smallest/biggest row/column(... entries
    MIN=min(CORNERS,[],2)';
    MINREAL=min(CORNERSREAL,[],2)';
    MAX=max(CORNERS,[],2)';
    N=abs(MIN)'; %Get offsets

    if fillvalue==0;
        d=zerosm(MAX-MIN+1); %Make output array of size of upsampled sequence
    else
        d=onesm(MAX-MIN+1)*fillvalue;
    end
    
    %Fill upsampled sequence with values    
    for i=0:numel(c)-1;
        SUB=cell(1,dim);
        [SUB{:}]=ind2sub(sizem(c,[],dim),i+1);
        SUB=M*(cell2mat(SUB)-1)'+1;
        SUB=num2cell(SUB+N)';
        d(SUB{:})=c(i+1);
    end
    
    dmin=MINREAL.';
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 