% the first noe transmission eigenvalues.
function k=sptarnite(A,B,lb,noe)
%--------iteration method-------------- 
rb=lb+1;
D=[];
count=0;
k=zeros(1,noe); 
while count<noe
        while test(D)   
        [V,D] = sptarn(A, B, lb, rb);
        lb = rb; rb = lb+ 1 ;
        end  
    result=(D.*(abs(imag(D))<=0)); % set the component with nonzero imaginary part in D to 0. 
    ind=find(result~=0); % find the nonzero component in result.
    if length(ind) <= (noe-count)
        k(1,count+1:count+length(ind)) = sqrt(sort(result(ind))); % get positve eigenvalues.
    else
        k(1,count+1:noe) = sqrt(sort(result(ind(1:(noe-count)))));
    end
    count=count+length(ind); % count the total number of eigenvalues.
    D=[]; % reset the D.
end

%--- Display the lowest a few transmission eigenvalue ---------------------
function mark=test(D) 
if length(D)==0 % if D is empty, do the loop.
    mark=1;
else % if there is no real number in D, do the loop.
    mark=min(abs(imag(D))>0);
end
