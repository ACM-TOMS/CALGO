function [EVals, EVecs] = eigifp(varargin) 
%
%
%   EIGIFP (version 2.1.1):  computes a few (algebraically) smallest or largest 
%       eigenvalues/eigenvectors of the symmetric matrix eigenvalue problem:
%
%             A x = lambda x   or    A x = lambda B x
%
%       where A and B are symmetric and B is positive definite.
%
%   [d, x] = eigifp(A)      => smallest eigenvalue/eigenvector of A;
%   [d, X] = eigifp(A,k)    => k smallest eigenvalues/eigenvectors of A;
%
%   [d, x] = eigifp(A,B)    => smallest eigenvalue/eigenvector of (A,B);
%   [d, X] = eigifp(A,B,k)  => k smallest eigenvalue/eigenvector of (A,B);
%
%   A (and/or B) can be either a sparse matrix or a function that 
%   returns A*x for x. In the latter case, the dimension of A has 
%   to be input through options.size (see below). It is also desirable  
%   to have an estimate of a norm of A (or B) inputed through 
%   options.NORMA (or options.NORMB). 
%
%   To compute the largest eigenvalues, one can call as above by adding 
%   the option options.MAXEIG = 1 (see below). Alternatively, if A is 
%   numeric, one can do the following.
%
%   [d, x] = eigifp(-A); d=-d;     => largest eigenvalue/eigenvector of A;
%   [d, X] = eigifp(-A,k); d=-d;   => k largest eigenvalues/eigenvectors of A;
%
%   [d, x] = eigifp(-A,B); d=-d;   => largest eigenvalue/eigenvector of (A,B);
%   [d, X] = eigifp(-A,B,k); d=-d; => k largest eigenvalues/eigenvectors of (A,B);
%
%   Other optional inputs can be provided to use any a priori information,
%   to control computational cost, and to optimize performance.
%  
%   eigifp(A,opt), eigifp(A,k opt), eigifp(A,B,opt), eigifp(A,B,k,opt)
%   take the following optional inputs defined by opt:  
%
%     opt.SIZE:     
%       the size of A. 
%     opt.NORMA:    
%       an estimate of the 1-norm of A 
%     opt.NORMB:    
%       an estimate of the 1-norm of B
%     opt.INITIALVEC or opt.V0: 
%       a matrix whose i-th column is the i-th initial eigenvector; 
%     opt.TOLERANCE or opt.TOL: 
%       termination tolerance for the 1-norm of residual; 
%       default: 10*eps*sqrt(n)*(||A||+|lambda|*||B||)
%     opt.MAXITERATION or opt.MAXIT: 
%       set the maximum number of (outer) iterations; default: 500;
%     opt.INNERITERATION or opt.INNERIT: 
%       set a fixed inner iteration to control the memory requirement; 
%       Default: between 1 and 128 as adaptively determined.   
%     opt.USEPRECON:    
%       set preconditioning off or completely on:  
%           1.  set opt.USEPRECON='NO' to switch off preconditioning; 
%           2.  set opt.USEPRECON = an approximate eigenvalue 
%               to start preconditioned iterations with a preconditioner 
%           computed using the eigenvalue as the shift;              
%     opt.ILUTHRESH or opt.ILU:     
%       threshold between 0 and 1 used in incomplete LU for computing  
%           preconditioner (i.e. used in luinc.m); default: 1e-4; 
%       Special cases: 0 = exact LU; 1 = 0-level ILU;  
%     opt.PRECONDITIONER or opt.PRECON:  
%       user supplied preconditioner. It can be a matrix or a function.
%     opt.DISP:
%       Set to 0 to disable on-screen display of output, and to any other 
%       numerical value to enable display. Default: 1
%     opt.MAXEIG:
%       Set to 1 to compute the largest eigenvalues instead of the 
%       smallest eigenvalue of (A,B). Any other value causes the code to 
%       compute the smallest eigenvalue. Default: 0
%
%   
%   EIGIFP is based on the algorithm in: G. Golub and Q. Ye, An Inverse 
%   Free Preconditioned Krylov Subspace Method for Symmetric Generalized
%       Eigenvalue Problems,  SIAM J. Sci. Comput., 24(2002):312-334. It also
%       incorporates several enhancements. 
%   
%       It is particularly suitable for problems where 
%       a)  factorization of B is difficult; or
%       b)  factorization of A-lambda_k B is difficult; or 
%       c)  a shift near eigenvalues sought is not known. 
%
%
%   EIGIFP is developed by Qiang Ye (qye@ms.uky.edu) with the programming 
%   assistance provided by James Money. This work was supported by NSF  
%   under Grant CCR-0098133. 
%   
%   This program is provided for research and educational use and  
%   is distributed through http://www.ms.uky.edu/~qye/software. 
%   Neither redistribution nor commercial use is permitted without 
%   consent of the author. 
%   
%   Copyright (c): Qiang Ye
%   Version 2.1.1 is dated Oct. 2004
%   Version 2.1 is dated July 2004
%   Version 2.0 is dated Aug. 2003
%   Version 1 is dated Sept. 2002
%

                    % Process inputs 
A=varargin{1};
n=[]; 
if ( isnumeric(A) )
    n=size(A,1);
    if (any(size(A) ~= n))
        error('Input Error: Matrix A is not square');
    end
elseif ( ~isstr(A) ) 
    error('Input Error: A must be either a matrix or a function');
end 

B=[];
if ( nargin>1 )
    if ( (~ isstruct(varargin{2})) ) 
        if( isstr(varargin{2}) | (any(size(varargin{2}) ~= 1)))
        B=varargin{2};
        if ( isnumeric(B) )
            if ( isempty(n) )
                n=size(B,1);
            end 
            if (any(size(B) ~= n))
                    error('Input Error: Matrix B must be square and of the size of A');
            end
        elseif ( ~isstr(B) ) 
            error('Input Error: B must be either a matrix or a function');
        end         
        end
    end
end

                    % set initial or default values 
m=0;
k=1;
X=[];
L=[];
eta=1e-4;
usePrecon=1;
iterMax=500;
ksetbyuser=0;
tolerance=[];
normA=[]; 
normB=[]; 
outputYes=1;
EVals=[];
EVecs=[];
findMax=0;

if ( nargin > 1 + ~isempty(B) )
    if ( (~ isstruct(varargin{2+~isempty(B)})) & all(size(varargin{2+~isempty(B)}) == 1))
        k=varargin{2+~isempty(B)};
        ksetbyuser=1;
    end
end
if ( nargin > 1 + ~isempty(B) + ksetbyuser )
    options=varargin{2 + ~isempty(B) + ksetbyuser:nargin};
    if ( ~isstruct(options) ) 
        error('Input Error: too many input parameters.'); 
    end
    names=fieldnames(options);
    
    I=strmatch('DISP',upper(names),'exact');
    if (~isempty(I))
        outputYesStr=getfield(options,names{I});
        if (isnumeric(outputYesStr))
            if (outputYesStr==0)
                outputYes=0;
            else
                outputYes=1;
            end
        end
    end
    I=strmatch('SIZE',upper(names),'exact');
    if (~isempty(I))
        n=getfield(options,names{I});
        if( isnumeric(A) & n ~= size(A, 1) )
            error('Input Error: options.size not equal to the size of A');
        end 
        if( ~isempty(B) & isnumeric(B) & n ~= size(B, 1) )
            error('Input Error: options.size not equal to the size of B');
        end 
    end
    
    I=strmatch('NORMA',upper(names),'exact');
    if (~isempty(I))
        normA=getfield(options,names{I});
    end

    I=strmatch('NORMB',upper(names),'exact');
    if (~isempty(I))
        normB=getfield(options,names{I});
    end

    I=strmatch('INITIALVEC',upper(names),'exact');
    if (isempty(I))
        I=strmatch('V0',upper(names),'exact');
    end
    if (~isempty(I))
        X=getfield(options,names{I});
        if( isempty(n) )
            n=size(X, 1); 
        elseif ( size(X, 1) ~=n )
            error('Input Error: incorrect size of initial vectors');
        end 
        if( size(X, 2) ~=k )
            fprintf(1,'Warning: Number of initial vectors not equal to k\n');
        end 
    end
    
    I=strmatch('TOLERANCE',upper(names),'exact');
    if (isempty(I))
        I=strmatch('TOL',upper(names),'exact');
    end
    if (~isempty(I))
        tolerance=getfield(options,names{I});
        if ( tolerance <= 0 )
                error('Input Error: Invalid tolerance input.');
            end 
    end
    
    if ( isempty(n) )
        error('Input Error: need to input options.size = (dimension of A)');
    end
     
    I=strmatch('INNERITERATION',upper(names),'exact');
    if (isempty(I))
        I=strmatch('INNERIT',upper(names),'exact');
    end
    if (~isempty(I))
        m=getfield(options,names{I});
        if ( m < 0 | m >= n )
                error('Input Error: Invalid inner iteration number');
            end     
    end
    
    I=strmatch('ILUTHRESH',upper(names),'exact');
    if (isempty(I))
        I=strmatch('ILU',upper(names),'exact');
    end 
    if (~isempty(I))
        eta=getfield(options,names{I});
        if ( eta < 0 )
                error('Input Error: Invalid incomplete factorization threshold.');
            end             
    end
    
    I=strmatch('PRECONDITIONER',upper(names),'exact');
    if (isempty(I))
        I=strmatch('PRECON',upper(names),'exact');
    end
    if (~isempty(I))
        L=getfield(options,names{I});
        if ( isnumeric(L) & any(size(L) ~= n)) 
                error('Input Error: Invalid preconditioner; must be square.');
            end             
    end
    
    I=strmatch('MAXEIG',upper(names),'exact');
    if (~isempty(I))
        findMax=getfield(options,names{I});
        if ( ~isnumeric(findMax) | ~((findMax==0) | (findMax==1)))
            fprintf(1,'Warning: opt.MAXEIG should be either 0 or 1. Reset to 0.\n');
        end             
    end
    
    I=strmatch('USEPRECON',upper(names),'exact');
    if (~isempty(I))
        usePreconIn=getfield(options,names{I});
        if ( ~isempty(usePreconIn) &  ~isempty(L) )
            error('Input Error: input opt.useprecon while a preconditioner has been provided.'); 
        end
        if (isstr(usePreconIn))
            if (strcmp(upper(usePreconIn),'NO'))
                usePrecon=0;                    
            end
        elseif (~isempty(usePreconIn))      
            if ( any(size(usePreconIn) > 1) ) 
                    error('Input Error: Invalid shift; must be a real number.');
                end         
                if ( isnumeric(A) & isnumeric(B) )
                if(~issparse(A)) 
                    A=sparse(A);
                    fprintf(1,'Warning: A is not in the sparse format.\n');             
                end  
                    if (outputYes ~=0) 
                    fprintf(1,'Computing threshold ILU factorization using the shift provided.\n');
                end
                if ( ~isempty(B) )
                    if(~issparse(B)) 
                        B=sparse(B);
                        fprintf(1,'Warning: B is not in the sparse format.\n');             
                    end  
                    L=ildlte(A-usePreconIn*B,eta);
                else
                    n=size(A,1); 
                    L=ildlte(A-usePreconIn*speye(n),eta);
                end
            else
                fprintf(1,'Warning: A (or B) is a function; can not compute preconditioner using opt.useprecon.\n');            
            end
        end
    end
    I=strmatch('MAXITERATION',upper(names),'exact');
    if (isempty(I))
        I=strmatch('MAXIT',upper(names),'exact');
    end
    if (~isempty(I))
        iterMax=getfield(options,names{I});     
        if (  iterMax <= 0 )
                error('Input Error: Invalid maximum iteration.');
            end             

    end
end

if ( isempty(normA) )
    if ( isnumeric(A) )
        normA=norm(A,1); 
    else
        normA=norm(feval(A,ones(n,1)), 1)/n; 
        fprintf(1,'Warning: options.Anorm (estimate of ||A||_1) is not provided.\n');   
    end
end 
if ( isempty(normB) )
    if ( isempty(B) )
        normB=1;
    elseif ( isstr(B) )
        normB=norm(feval(B,ones(n,1)), 1)/n; 
        fprintf(1,'Warning: no options.Bnorm (estimate of ||B||_1) is not provided.\n');    
    elseif ( isnumeric(B) )
        normB=norm(B,1);
    end
end 

                        % set default tolerance 
toleranceA=10*eps*normA*sqrt(n);
if(~isempty(B)) 
   toleranceB=10*eps*normB*sqrt(n);
else 
   toleranceB=0; 
end 
if ( ~isempty(tolerance) )
    toleranceVec=[tolerance, 0];
    if ( tolerance < toleranceA )  
        fprintf(1,'Warning: Tolerance may be set too low. Suggested value: %E.\n',toleranceA);
    end 
else
    toleranceVec=[toleranceA,toleranceB];
end

if ( isstr(A) | isstr(B) ) 
    usePrecon=0;
end

if ( isempty(L) & (usePrecon ~= 0) ) 
    if(~issparse(A)) 
        A=sparse(A);
        fprintf(1,'Warning: A is not in the sparse format.\n');             
    end  
    if ( ~isempty(B) & ~issparse(B)) 
        B=sparse(B);
        fprintf(1,'Warning: B is not in the sparse format.\n');             
    end
end  
                    
if (k>n)
    error('Error: # of eigenvalues sought is greater than the matrix size');
end

                        % call the main function ifree.m

if (isnumeric(A))       %if given a matrix, use (A,B) or (-A,B) if maxeig =1
    if (findMax==1)
        if (outputYes)
            fprintf(1,'Computing the smallest eigenvalues of (-A,B) first.\n');
        end
        [EVals, EVecs] = ifree(-A,B,n,m,toleranceVec,iterMax,k,X,eta,L,usePrecon,normA,normB,outputYes);
        if (outputYes)
            fprintf(1,'  Negating the eigenvalues of (-A, B) => the largest eigenvalues of (A,B):\n');
            EVals=-EVals;
            fprintf(1,'  %e\n',EVals);
        end        
    else
        [EVals, EVecs] = ifree(A,B,n,m,toleranceVec,iterMax,k,X,eta,L,usePrecon,normA,normB,outputYes);
    end
else                    % if A is a function 
    if (findMax==1)     % here we use minusA and negate the results when done
        if (outputYes)
            fprintf(1,'Computing the smallest eigenvalues of (-A,B) first.\n');
        end
        saveAfunc=0;
        if (exist('Afunc123')==1)   % we don't want to overwrite the user's global, so we save it
            saveAfunc=1;
            tempFunc=Afunc;
            clear Afunc123;
        end
        global Afunc123;
        Afunc123=A;     % this is the function it calls in minusA()
        [EVals, EVecs] = ifree('minusA',B,n,m,toleranceVec,iterMax,k,X,eta,L,usePrecon,normA,normB,outputYes);
        if (saveAfunc==1)
            Afunc123=saveAfunc; %restore the value if the user already defined it
        end
        if (outputYes)
            fprintf(1,'  Negating the eigenvalues of (-A, B) => the largest eigenvalues of (A,B):\n');
            EVals=-EVals;
            fprintf(1,'  %e\n',EVals);
        end      
    else
        [EVals, EVecs] = ifree(A,B,n,m,toleranceVec,iterMax,k,X,eta,L,usePrecon,normA,normB,outputYes);
    end
end





function [lambda, x, r, bx, diff, rDiff, bDiff, lam_max, res] = ...
     parnob(A, B, L, lambda, x, r, bx, diff, rDiff, bDiff, m, tol, EVecs, BEVecs, lambda_1, lambda_max)
%  
%   generate B-othonormal basis V of preconditioned (by L) Krylov subspace
%   by m steps of Arnoldi algorithm and then compute the Ritz value  
%   and Ritz vectors 
%   
%   BEVecs = B*EVecs  with EVecs =  converged eigenvectors
%

                    % initialization and normalization 
n=size(x,1);                     
V = zeros(n,m);
Wr=V; 
Am = zeros(m,m);
temp=bx'*x;
if (temp <= 0)
    if ( isnumeric(B) )     % double check on x'Bx
        bx=B*x;
    else
        bx=feval(B, x);
    end
    temp=bx'*x;
    if (temp <= 0)
        error('Error: B is not positive definite.');
    end
end
temp=sqrt(temp);
V(:,1) = x / temp;
r=r/temp; 
Wr(:,1) = r;
bx=bx/temp;
if ( ~ isempty(B) )
        Wb=V; 
    Wb(:,1) = bx;
end
Am(1,1)=0;

                    % Loop for Arnoldi iteration
for i = 2:m-1,                    
                                    % Apply preconditioner if given 
    if ( isnumeric(L) )         
        r=L\r;
        r=(L')\r;    
    else 
        r=feval(L,r); 
    end
                    % generate new basis vector 
    for k = 1:(i-1)
        if ( ~ isempty(B) )
            temp= Wb(:,k)'*r;
        else
            temp= V(:,k)'*r;
        end
        r = r - temp*V(:,k);
    end
                        % reorthogonalization if m > 6 
    if( m > 6)  
        for k = 1:(i-1),        
            if ( ~ isempty(B) )
                temp= Wb(:,k)'*r;
            else
                temp= V(:,k)'*r;
            end
            r = r - temp*V(:,k);
            end
        end
    if (norm(r) == 0)
        m=i;
        break 
    end
                    % normalize and save new basis vector 
                    % as well as Bx and (A-lambda B)x 
    if ( (~isempty(B)) & isnumeric(B) )
        bx=B*r;
    elseif ( (~isempty(B)) & isstr(B) )
        bx=feval(B, r);
    else 
        bx=r; 
    end
    
    temp=bx'*r;
    if (temp <= 0)
        error('Error: B is not positive definite.');
    end
    temp = sqrt(temp);
    
    V(:,i) = r;
    if (isstr(A))
        r=feval(A, V(:,i))-lambda*bx;
    else
        r=A*V(:,i)-lambda*bx;
    end
    V(:,i) = V(:,i) / temp;
    r=r/temp;
    if ( ~ isempty(B) )
        Wb(:,i)=bx/temp;
    end
    Wr(:, i) = r;   
                    % deflation for converged e-vectors
    if ( nargin == 16 )
        count=size(BEVecs,2);
        for j=1:count
            sigma=-lambda_1(j) + lambda_max(count);
            r=r  + sigma*(BEVecs(:,j)'*V(:,i))*BEVecs(:,j);
        end 
    end
    Am(1:i,i)=V(:,1:i)'*r;
end 

                    % add the diff vector to the basis
                    % and complete the projection Am 
diffNorm=sqrt(bDiff'*diff);
for k = 1:m-1,                       
    if ( ~ isempty(B) )
        temp = Wb(:,k)'*diff;
            bDiff = bDiff - temp*Wb(:,k);
    else
        temp = V(:,k)'*diff;
    end
        diff = diff - temp*V(:,k);
    rDiff = rDiff - temp*Wr(:,k);
end
                    % reorthogonalization if m > 6 
if( m > 6) 
    for k = 1:m-1,                  
        if ( ~ isempty(B) )
            temp = Wb(:,k)'*diff;
                    bDiff = bDiff - temp*Wb(:,k);
        else
            temp = V(:,k)'*diff;
        end
            diff = diff - temp*V(:,k);
        rDiff = rDiff - temp*Wr(:,k);
    end
end 
if ( ~ isempty(B) )
    temp = bDiff'*diff;
        if (temp < 0)
        if ( isnumeric(B) )
            bx=B*x;
        else
            bx=feval(B, x);
        end
        temp=bx'*x;   
        if (temp < 0)       
            error('Error: B is not positive definite.');
        end
    end
        temp=sqrt(temp); 
else
    temp = norm(diff);
end

                    % check and add diff only if it's significant
if ( temp <= 1e-8*diffNorm | temp==0 )
        m=m-1;
elseif ( temp <= 1e-2*diffNorm )    % recompute (A-lambda B)diff if necessary
    V(:,m)=diff/temp; 
    if ( (~isempty(B)) & isnumeric(B) )
        bDiff=B*V(:,m);
        Wb(:,m)=bDiff;
    elseif ( (~isempty(B)) & isstr(B) )
        bDiff=feval(B, V(:,m));
        Wb(:,m)=bDiff;
    else 
        bDiff=V(:,m); 
    end
    if (isstr(A))
        rDiff=feval(A, V(:,m))-lambda*bDiff;
    else
        rDiff=A*V(:,m)-lambda*bDiff;
    end    
    Wr(:,m)=rDiff;
    if (nargin == 16)
        count=size(BEVecs,2);
        for j=1:count
            sigma=-lambda_1(j) + lambda_max(count);
            rDiff=rDiff+sigma*(BEVecs(:,j)'*V(:,m))*BEVecs(:,j);
        end
    end     
    Am(1:m,m)=V(:,1:m)'*rDiff;    
else     
    V(:,m)=diff/temp; 
    Wr(:,m)=rDiff/temp;
    if ( ~ isempty(B) ), Wb(:,m)=bDiff/temp;   end
    r=Wr(:,m);
    if (nargin == 16)
        count=size(BEVecs,2);
        for j=1:count
            sigma=-lambda_1(j) + lambda_max(count);
            r=r+sigma*(BEVecs(:,j)'*V(:,m))*BEVecs(:,j);
        end
    end     
    Am(1:m,m)=V(:,1:m)'*r;
end
Am=triu(Am) + triu(Am,1)';

                    % compute Ritz value and vector of projection
[U, D]=eig(Am(1:m,1:m));
[delta, Ieig]=sort(diag(D));    

U(:,Ieig(1))=U(:,Ieig(1))/norm(U(:,Ieig(1))); 
x=V(:,1:m)*U(:,Ieig(1));
if ( ~ isempty(B) ),  bx=Wb(:,1:m)*U(:,Ieig(1)); end    
r=Wr(:,1:m)*U(:,Ieig(1));
                    % post processing Ritz vector
if (nargin == 16)               
    for j=1:count
            temp=BEVecs(:,j)'*x;
            x=x-EVecs(:,j)*temp;
            if ( ~ isempty(B) ),
                    bx=bx-BEVecs(:,j)*temp;
            end
            r=r-(lambda_1(j)-lambda)*BEVecs(:,j)*temp;
    end
end
if ( isempty(B) ), bx=x; end
sigma=(x'*r)/(x'*bx);
lambda=lambda+sigma; 
r=r-sigma*bx; 
                    % update new diff and related vectors
                    % for the next iteration   
U(1,Ieig(1)) = -(U(2:m,Ieig(1))'*U(2:m,Ieig(1)))/U(1,Ieig(1)); 
diff=V(:,1:m)*U(1:m,Ieig(1));
if ( ~ isempty(B) ),  
        bDiff=Wb(:,1:m)*U(1:m,Ieig(1)); 
else
        bDiff=diff;
end 
rDiff=Wr(:,1:m)*U(1:m,Ieig(1)); 
rDiff=rDiff-sigma*bDiff;

temp=norm(x,1); 
res=norm(r,1)/temp; 
lam_max=lambda+delta(size(delta,1));     
tol0=tol(1) + abs(lambda)*tol(2);
if ( res < tol0 | m>10 )        % recompute bx and r if necessary 
    if ( (~isempty(B)) & isnumeric(B) )
        bx=B*x;  
    elseif ( (~isempty(B)) & isstr(B) )
        bx=feval(B, x);
    else
        bx=x; 
    end 
    if (isstr(A))
        r=feval(A, x);
    else
        r=A*x;
    end
    lambda = (x'*r)/(x'*bx); 
    r=r-lambda*bx;
    res=norm(r,1)/temp;
    tol0=tol(1) + abs(lambda)*tol(2);   
end

                    % use a new diff if converged. 
if ( res < tol0 )
        diff = V(:,1:m)*U(:,Ieig(2));
end


function [lambda, x, r, bx, diff, rDiff, bDiff, d2, lam_max, res] = ...
     lancb(A, B, lambda, x, r, bx, diff, rDiff, bDiff, m, tol, EVecs, BEVecs, lambda_1, lambda_max)
%  
%   generate othonormal basis V of Krylov subspace by m steps of lanczos
%   and then compute the Ritz value and Ritz vectors 
%
%   BEVecs = B*EVecs  with EVecs =  converged eigenvectors
%

                    % initialization and normalization 
n=size(x,1); 
V = zeros(n,m);
Wr = V; 
Am=zeros(m,m);
beta=norm(x);
V(:,1) = x/beta; 
r=r/beta; 
bx=bx/beta;
Wr(:,1) = r;
if ( ~ isempty(B) ), 
    Bm=zeros(m,m); 
    Wb=V; 
    Wb(:,1) = bx;
end
beta=0;
                    % loop for Lanczos iteration 
for i = 1:(m-2)                          
    r=r-beta*x;     
    x=V(:,i);   
    alpha=x'*r; 
    r=r-alpha*x;
                    % reorthogonalization if m > 6
    if( m > 6 ) 
        for k = 1:i,                
            temp= V(:,k)'*r;
            r = r - temp*V(:,k);
        end
    end
                    % construct projection      
    beta=norm(r);
    Am(i,i)=alpha; 
    Am(i+1, i)=beta;
    Am(i, i+1)=beta;    
    if ( ~ isempty(B) ), Bm(1:i,i)=V(:,1:i)'*bx; end
    
    if (beta == 0)           
             break;
    end 
        V(:,i+1)=r/beta; 
        
                    % generate new basis vector     
    if ( (~isempty(B)) & isnumeric(B) )
        bx=B*V(:,i+1);  
        Wb(:,i+1)=bx;
    elseif ( (~isempty(B)) & isstr(B) )
        bx=feval(B, V(:,i+1));
        Wb(:,i+1)=bx;
    else 
        bx=V(:,i+1);  
    end 
    if (isstr(A))
        r=feval(A, V(:,i+1))-lambda*bx;
    else
        r=A*V(:,i+1)-lambda*bx;
    end
    Wr(:, i+1) = r; 
                    % deflation for converged e-vectors
    if (nargin == 15)
        count=size(BEVecs,2);
        for j=1:count
            sigma=-lambda_1(j) + lambda_max(count);
            r=r+sigma*(BEVecs(:,j)'*V(:,i+1))*BEVecs(:,j);
        end
    end
end 
                    % complete projection for the last vec
if (beta == 0) 
    m=i+1; 
else        
    Am(m-1, m-1) = V(:,m-1)'*r;
    if ( ~ isempty(B) ), Bm(1:m-1,m-1)=V(:,1:m-1)'*bx; end
end 

                    % add the diff vec to the basis
                    % and complete projection 
diffNorm=norm(diff);
for k = 1:(m-1),                         
    temp= V(:,k)'*diff;
    diff = diff - temp*V(:,k);
    rDiff = rDiff - temp*Wr(:,k);
    if ( ~ isempty(B) ), bDiff=bDiff-temp*Wb(:,k);  end
end
if ( m > 6)             % reorthogonalization if m>6 
    for k = 1:(m-1),                         
        temp= V(:,k)'*diff;
        diff = diff - temp*V(:,k);
        rDiff = rDiff - temp*Wr(:,k);
        if (~isempty(B)), bDiff=bDiff-temp*Wb(:,k);  end
    end
end 
temp = norm(diff); 

                    % check and add diff only if it's significant
if ( temp <= 1e-8*diffNorm | temp==0 )
        m=m-1; 
    if ( ~ isempty(B) ), 
            Bm=triu(Bm) + triu(Bm,1)';
        end
elseif ( temp <= 1e-2*diffNorm )
    V(:,m)=diff/temp; 
        if ( (~isempty(B)) & isnumeric(B) )
        bDiff=B*V(:,m);
            Wb(:,m)=bDiff;
    elseif ( (~isempty(B)) & isstr(B) )
        bDiff=feval(B, V(:,m));
            Wb(:,m)=bDiff;
    else 
        bDiff=V(:,m); 
    end
    if (isstr(A))
        rDiff=feval(A, V(:,m))-lambda*bDiff;
    else
        rDiff=A*V(:,m)-lambda*bDiff;
    end    
    Wr(:,m)=rDiff;
    if (nargin == 15)
        count=size(BEVecs,2);
        for j=1:count
            sigma=-lambda_1(j) + lambda_max(count);
            rDiff=rDiff+sigma*(BEVecs(:,j)'*V(:,m))*BEVecs(:,j);
        end
    end     
    Am(1:m,m)=V(:,1:m)'*rDiff;    
    Am(m,1:m)=Am(1:m,m)'; 
    if ( ~ isempty(B) ),  
        Bm(1:m,m)=V(:,1:m)'*Wb(:,m);
        Bm=triu(Bm) + triu(Bm,1)'; 
    end
else    
    V(:,m)=diff/temp; 
    Wr(:,m)=rDiff/temp;
    if ( ~ isempty(B) ), Wb(:,m)=bDiff/temp;   end
    r=Wr(:,m);
    if (nargin == 15)
        count=size(BEVecs,2);
        for j=1:count
            sigma=-lambda_1(j) + lambda_max(count);
            r=r  + sigma*(BEVecs(:,j)'*V(:,m))*BEVecs(:,j);
        end
    end
        
    Am(1:m,m)=V(:,1:m)'*r;
    Am(m,1:m)=Am(1:m,m)'; 
    if ( ~ isempty(B) ),  
        Bm(1:m,m)=V(:,1:m)'*Wb(:,m);
        Bm=triu(Bm) + triu(Bm,1)'; 
    end
end

                    % compute Ritz value and vector and
                    % approx gap d2 used for error estimate 
if (~ isempty(B) )
    [U, D]=eig(Am(1:m,1:m),Bm(1:m,1:m),'chol');
    [delta, Ieig]=sort(diag(D));                    
    EVals2=Bisection(Am(1:m,1:m), 0.1); 
    EVals2=sort(EVals2);
    if ( EVals2(2) <= 0 ) 
        d2=-1; 
    elseif ( EVals2(1) > 0 )
        d2=eps;
    else 
        d2=max(abs(EVals2));
    end     
else
    [U, D]=eig(Am(1:m,1:m));
    [delta, Ieig]=sort(diag(D));                
    EVals2=delta(1:2);
    if ( EVals2(2) <= 0 ) 
        d2=-1; 
    elseif ( EVals2(1) > 0 )
        d2=eps;
    else 
        d2=max(abs(EVals2));
    end                     
end

U(:,Ieig(1))=U(:,Ieig(1))/norm(U(:,Ieig(1))); 

x=V(:,1:m)*U(:,Ieig(1));
if ( ~ isempty(B) ),  bx=Wb(:,1:m)*U(:,Ieig(1)); end    
r=Wr(:,1:m)*U(:,Ieig(1));

if (nargin == 15)                   % post processing Ritz vector
    for j=1:count
            temp=BEVecs(:,j)'*x;
            x=x-EVecs(:,j)*temp;
            if ( ~ isempty(B) ),
                    bx=bx-BEVecs(:,j)*temp;
            end
            r=r-(lambda_1(j)-lambda)*BEVecs(:,j)*temp;
    end
end
if ( isempty(B) ), bx=x; end
sigma=(x'*r)/(x'*bx);
lambda=lambda+sigma; 
r=r-sigma*bx;

                    % update new diff and related vectors
                    % for the next iteration  
U(1,Ieig(1)) = -(U(2:m,Ieig(1))'*U(2:m,Ieig(1)))/U(1,Ieig(1)); 
diff=V(:,1:m)*U(1:m,Ieig(1));
if ( ~ isempty(B) ),  
        bDiff=Wb(:,1:m)*U(1:m,Ieig(1)); 
else
        bDiff=diff;
end 
rDiff=Wr(:,1:m)*U(1:m,Ieig(1)); 
rDiff=rDiff-sigma*bDiff;

temp=norm(x,1); 
res=norm(r,1)/temp; 
lam_max=lambda+delta(size(delta,1));    
tol0=tol(1) + abs(lambda)*tol(2);
if ( res < tol0 | m> 10)        % recompute bx and r if necessary 
    if ( (~isempty(B)) & isnumeric(B) )
        bx=B*x;  
    elseif ( (~isempty(B)) & isstr(B) )
        bx=feval(B, x);
        else
            bx=x; 
    end 
    if (isstr(A))
        r=feval(A, x);
    else
        r=A*x;
    end
    lambda = (x'*r)/(x'*bx); 
    r=r-lambda*bx;
    res=norm(r,1)/temp;
    tol0=tol(1) + abs(lambda)*tol(2);   
end


if (res < tol0 )            % use new diff if converged 
        diff = V(:,1:m)*U(:,Ieig(2));
end


 
function [Lambda_1, EVecs] = ...
    ifree(A, B, n, m, tol, itermax, k, X, eta, L, usePrecon, normA, normB,outputYes);
%
%  The main function that carries out the (outer) iteration. 
% 
t=cputime;
                    % initialization 
conv=ones(itermax,1);
rate=ones(10,1);

Lambda_1=zeros(k,1);
Lambda_Max=zeros(k,1);
EVecs=zeros(n,k);
BEVecs=EVecs;
pertbound=0; 
shift=0; 
if ( isempty(L) )
    regenL=1;
    startL=0;
else
    regenL=0;
    startL=1;
end

                    % set initial random vector if not given
if (~ isempty(X))
    x=X(:,1);
    normX=norm(x);
    if (normX == 0 | size(x,1) ~= n)
        error('Input Error: Invalid input of the initial vector.');
    end
else
    x=rand(n,1)-0.5*ones(n,1);
    normX=norm(x);
end 
x=x/normX;

                    % loop for computing k eigenvalues
for l=1:k

        diff=zeros(n,1);        % initialization for each eigenvalue iteration
        rDiff=diff; 
        bDiff=diff;
    
    if ( (~isempty(B)) & isnumeric(B) )
        bx=B*x; 
    elseif ( (~isempty(B)) & isstr(B) )
        bx=feval(B, x);
    else 
        bx=x;  
    end 

    temp=x'*bx; 
    if (temp <= 0)
        error('Error: B is not positive definite.');
    end

    if (isstr(A))
        r=feval(A, x);
    else
        r=A*x;
    end
        
    lambda=(x'*r)/temp;
    r=r-lambda*bx;
    res=norm(r);
    conv(1)=res;
    initialConv=0; 
    tol0=tol(1) + abs(lambda)*tol(2);
    
    matvecCount=1;
    preconCount=0;
    changeCounter=-itermax;
    mValue=m;
    if (mValue < 1)
        mValue=min(n-1,2);
        if (startL == 1), mValue=1; end
        changeCounter=-10;
    end
    priorM= mValue; 
        priorRate=-eps;
        
        if (outputYes ~=0)
        fprintf(1,'\nComputing Eigenvalue %d:\n',l);
        fprintf(1,'  Iteration\tEigenvalue\tResidual\n');
        fprintf(1,'  ======================================================\n');
        fprintf(1,'  %d\t\t%E\t%E\n',1,lambda,res); 
    end
    
                    % iteration for computing l-th eigenvalue
    for iter=2:itermax
        
        if ( res < tol0 )   % if converged, set parameters for next eigen
            lambda_max = normA; 
            initialConv=1;  
            break;
        end
        projSize=mValue+2;  
                    % compute Ritz value and vector by Krylov projection   
                    % call lancb if no preconditioning 
                    % call parnob if using preconditioning
        if (l>1)
                    % with deflation for l > 1      
          if ( startL == 0 )
                    [lambda, x, r, bx, diff, rDiff, bDiff, d2, lambda_max, res] = ...
              lancb(A, B, lambda, x, r, bx, diff, rDiff, bDiff, projSize,  tol,...
                          EVecs(:,1:l-1), BEVecs(:,1:l-1), Lambda_1(1:l-1), Lambda_Max(1:l-1,:));
          else
                    preconCount=preconCount+mValue;
                    [lambda, x, r, bx, diff, rDiff, bDiff, lambda_max, res] = ...
              parnob(A, B, L, lambda, x, r, bx, diff, rDiff, bDiff, projSize, tol,...
                          EVecs(:,1:l-1), BEVecs(:,1:l-1), Lambda_1(1:l-1), Lambda_Max(1:l-1,:));
          end           
        else
                    % no deflation if l=1. 
          if ( startL == 0 )
                    [lambda, x, r, bx, diff, rDiff, bDiff, d2, lambda_max, res] = ...
              lancb(A, B, lambda, x, r, bx, diff, rDiff, bDiff, projSize, tol);             
          else
            preconCount=preconCount+mValue;
                    [lambda, x, r, bx, diff, rDiff, bDiff, lambda_max, res] = ...
              parnob(A, B, L, lambda, x, r, bx, diff, rDiff, bDiff, projSize, tol);
          end           
        end
             
        conv(iter)=res; 
        matvecCount=matvecCount+mValue;         
        if (outputYes ~= 0)
            fprintf(1,'  %d\t\t%E\t%E\n',iter,lambda,res);
        end
        
                    % update tolerance and check convergence
        tol0=tol(1) + abs(lambda)*tol(2);
        if ( res <= tol0 )
            break;
        end

                    % check on convergence rate and update mValue
        changeCounter=changeCounter+1;
        if ( changeCounter >= 19 )
                    rate=[rate(2:10); aveRate(conv, iter-changeCounter+4, iter)];
                    [mValue, priorM, priorRate, fixM]=updateM(rate,mValue,priorM,priorRate,n);
                    changeCounter=(changeCounter+fixM*itermax)*(1-fixM);  
            elseif ( changeCounter >= 10 )
                    rate=[rate(2:10); aveRate(conv, iter-changeCounter+4, iter)];
            end
 
                    % estimate error 
        if ( (startL==0) & (d2>0) )
            err=min((res^2 / d2), res)/(bx'*x); 
        elseif (startL==0)
            err =abs(lambda);
        end

                    % determine whether to switch to preconditioning            
        if ((startL==0) & (usePrecon ~= 0) & (err <= (0.01*abs(lambda))) & (err < (0.0001*max(abs(lambda),normA/normB)) ) )
            startL=1;
            mValue=1;
            changeCounter=0;
            pertbound=0.0002*max(abs(lambda),normA/normB);          
            if (outputYes ~= 0)
                fprintf(1,'   \n');
                fprintf(1,'    Computing ILU factorization: if it takes too long, try the options:\n'); 
                fprintf(1,'      opt.iluThresh, opt.Preconditioner, or opt.usePrecon. \n');         
                fprintf(1,'   \n');
            end
            shift=lambda-err;
            if ( ~isempty(B) )
                L=ildlte(A-shift*B,eta);
            else
                L=ildlte(A-shift*speye(n),eta);
            end
            if (outputYes ~= 0)
                fprintf(1,'    Done!  Starting preconditioned iterations: \n');
                fprintf(1,'      if convergence is not accelerated, try decreasing opt.iluThresh. \n');
                fprintf(1,'   \n');         
            end     
        end
                        % determine whether to switch off preconditioning
        if ((regenL ~= 0) & (startL==1) & (changeCounter>15) & ( abs(shift-lambda) >= min(0.02*abs(shift), pertbound)))
            if ( sum(rate) > -0.2)  
                startL=0;
                changeCounter=-itermax;
                mValue=m;
                if (mValue < 1)
                    mValue=2;
                    changeCounter=0;
                end         
            end 
        end
    end

                        % store eigenpairs and others
    Lambda_1(l)=lambda;
    Lambda_Max(l)=lambda_max;
    temp=sqrt(bx'*x);
    EVecs(:,l)=x/temp; 
    BEVecs(:,l)=bx/temp;

                        % warn if not converged                         
    if ( res >= tol0 )
        if (isnumeric(A) & nnz(A-A') ~= 0)
            error('Input Error: Matrix A is not symmetric');
        end
        if ( (~isempty(B)) & isnumeric(B) & nnz(B-B') ~= 0)
            error('Input Error: Matrix B is not symmetric');
        end
        fprintf(1,'\n');
        fprintf(1,'Warning: Eigenvalue %d not converged to the tolerance within max iteration\n',l);
        fprintf(1,'         Residual = %E , the set tolerance = %E \n',res, tol0);
        fprintf(1,'         Try setting opt.innerit, a larger tolerance, and/or opt.maxit\n');
    end
    if (outputYes ~= 0)
        fprintf(1,'  ------------\n');
        fprintf(1,'  Eigenvalue %d converged. \n',l);       
        fprintf(1,'  # of multiplications by A (and B):      %d. \n',1+matvecCount);        
        if (preconCount > 0)
            fprintf(1,'  # of multiplications by preconditioner: %d.\n',  preconCount);     
        end 
    end

                        %  set next initial vector  
    if ( l<size(X,2) )
        x=X(:,l+1);
        normX=norm(x);
        if (normX == 0)
            fprintf(1,'Input Error: Invalid input of initial vector.\n');
        end
        x=diff;
    elseif (initialConv==1) 
        x=rand(n,1)-0.5*ones(n,1);  
    else
        x=diff;
    end
    x=x-EVecs(:,1:l)*(BEVecs(:,1:l)'*x);
    normX=norm(x);
    if (normX == 0)
        if (l<k) 
            error('Error: Fail to form an initial vector.');
        end 
        return;
    end
    x=x/normX;              

end
if (outputYes ~= 0)
    fprintf(1,'\n-----\n');
    fprintf(1,'  CPU Time:%f.\n',cputime-t);
end

function rate=aveRate(conv, k, iter)
%
% compute average linear convergence rate of conv over steps k+1 to iter 
%
rate=0; 
if ( iter-k < 2 | k<1 ) 
    return;
end

y=log10(conv((k+1):iter));
xAve=(iter-k+1)/2; 
xyAve=((1:iter-k)*y)/(iter-k)-xAve*sum(y)/(iter-k); 
xAve=((iter-k)^2-1)/12;
rate=xyAve/xAve;

        
  
  
function L=ildlte(A, eta)
%
%ILDLT   Threshold incomplete LDL^T  using luinc.m or ilu 
%

n=size(A, 1);
    
                    % compute L U factors 
options.thresh=0;
options.udiag=1;
if (eta ==2)
    options.milu=1; 
    [L, U] = luinc(A, options);
elseif (eta==1) 
    [L, U] = ilu(A);
else 
    options.droptol=eta;
    [L, U] = luinc(A, options);
end 

                    % scale diagonals to get L 
d = diag (U);
d = sqrt(abs(d));  
for k = 1:n, 
    if ( d(k) < 1e-8) 
        d(k) = 1e-8;
    end
end 
L = L * spdiags(d, 0, n, n);


function [L, U] = ilu(A)
%
% ILU   Incomplete LU factorization with no fill-in.
%   No pivoting is used 
%    

n = size(A, 1); 
M = A; 

for k = 1:(n-1), 
    if ( M(k, k) == 0)
        M(k, k) = 1e-8;
    elseif ( abs(M(k, k)) < 1e-8) 
        M(k, k) = 1e-8*sign(M(k, k));
    end 
    ind = find( M(:, k) ); 
    ind = ind(find(ind>k))';
    for i = ind,
            M(i, k) = M(i, k) / M(k, k); 
    end 

    ind_j = find( M(k, :) ); 
    ind_j = ind_j(find(ind_j>k))';

    for j = ind_j, 
            for i = ind,
                if (M(i, j) ~= 0) 
                    M(i, j) = M(i, j) - M(i, k)*M(k,j); 
                end 
            end 
    end 
end 

L=tril(M, -1)+speye(n); 
U=triu(M, 0);


function [mValue, priorM, priorRate, fixM]=updateM(rate, mValue, priorM, priorRate,n);          
%
%  Adaptive update of mValue: inner iteration 
%
fixM=0;
maxm=min(n-1, 128);

                    % update m when rate stagnates  
if ( (max(rate)-min(rate)) < 0.1*(-rate(10)) | min(rate) > 0 )
    k=2;                % increase m by k times, 
                    % use larger k if slower convergence 
        if ((rate(10) > -0.001) & 8*mValue <= maxm), 
            k=8; 
        elseif ((rate(10) > -0.01) & 4*mValue <= maxm), 
            k=4;
        end
                        % increase m by testing acceleration rate 
        incFlag=(rate(10)/priorRate)*(priorM/mValue); 
        if (incFlag > 1.05) 
            if(2*mValue > maxm )
                fixM=-1;
            else
                priorM=mValue; 
                priorRate=rate(10); 
                mValue=k*mValue;
                fixM=1;
            end  
        elseif ((rate(10) > -0.001 ) & 2*mValue <= maxm), 
            mValue=k*mValue;
            fixM=1;  
        elseif ((rate(10) > -0.01 ) & 2*mValue <= maxm), 
            mValue=k*mValue;
            fixM=1;   
        elseif (incFlag < 0.9)
            mValue=priorM;
            fixM=-1;   
        end
end


function lambda=Bisection(T,tol)
%
% estimate two lowest eigenvalues of tridiagonal T by Bisection 
%

a=diag(T);
b=diag(T,1);
k=2; 
lambda=zeros(k,1); 
xb=norm(T, 1)*1.0000001; 
xa=-xb; 
xc=0; 

p=0;
x1=[];
x2=[];
n1=[];
n2=[];
na=0; 
nc=Negcount(a,b,xc);

[p,x1,x2,n1,n2]=put(xa,na,xc,nc,p,x1,x2,n1,n2);
if( nc < k)
    nb=Negcount(a,b,xb);
    [p,x1,x2,n1,n2]=put(xc,nc,xb,nb,p,x1,x2,n1,n2);
end 
q=0;
while ( (p~=0) & (q < k) )
    [low,nlow,up,nup,p,x1,x2,n1,n2]=remove(p,x1,x2,n1,n2);
    if(nlow <= k) 
        if ( abs(up-low)<tol*abs(up) ) 
            ind=q+nup-nlow; 
            ind=min(ind, k); 
            q=q+1;
            lambda(q:ind)=((up+low)/2)*ones((ind-q+1),1);
        else
            mid=(low+up)/2;
            nmid=Negcount(a,b,mid);  
            if ( (nup>nmid) & (nmid <k) )
                [p,x1,x2,n1,n2]=put(mid,nmid,up,nup,p,x1,x2,n1,n2);
            end
            if nmid>nlow  
                [p,x1,x2,n1,n2]=put(low,nlow,mid,nmid,p,x1,x2,n1,n2);
            end
        end
    end 
end 


function [p,x1,x2,n1,n2]=put(xa,na,xb,nb,p,x1,x2,n1,n2)

p=p+1;
x1(p)=xa;
x2(p)=xb;
n1(p)=na;
n2(p)=nb;


function [xa,na,xb,nb,p,x1,x2,n1,n2]=remove(p,x1,x2,n1,n2)
xa=x1(p);
xb=x2(p);
na=n1(p);
nb=n2(p);
p=p-1;

 
function m=Negcount(a,b,z)
%
% NEGCOUNT: evaluates the number of eigenvalues that are smaller than z
%

m=0;
n=size(a,1);
f=size(b,1);

d(1)=a(1)-z;
if d(1)<0
    m=m+1;
elseif d(1)==0,   
        d(1)=eps^2; 
end
for i=2:n
    d(i)=a(i)-z-b(i-1)^2/d(i-1);
    if d(i)<0
        m=m+1;
        elseif d(i)==0
            d(i)=eps^2; 
    end
end


function x=minusA(x)
global Afunc123;

x=feval(Afunc123,x);
x=-x;


