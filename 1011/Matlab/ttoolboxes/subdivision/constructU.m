function [ U, ev, r, errorflag ] = constructU( varargin  )
% [ U, ev, r, errorflag ] = constructU( T, k, [options] )
% Constructs a basis for the space U as described in Charina, Protasov, 2017.
% This function is not well tested.
%
% Input:
% =====
%  T            flat set of transition matrices
%
%  k            is obligatory, but can be empty
%               empty:      all eigenvectors v of T{N} are taken as starting vector.                                               Returns cell array
%               scalar:     starting vector v = eigenvector of T{idx} (idx=1 default) to the k^th eigenvalue (zeroth eigenvalue==1).   Returns matrix
%               vector:     all eigenvectors v = eigenvector of T{idx} (idx=1 default) to the k{i}^th eigenvalue                       Returns cell array
%
% Options:
% ========
%  'sym'                computation is done symbolically. Output is double.
%  'verbose',val        defines verbose level
%  'maxerror',val       maximal error for automatically computed starting vector, default: 1e-15
%  'vpa',val            computes with variable precision arithmetic with precision val. Prints 'X' to the console, if a new vector is added to U.
%  'nosimplify'         does not remove linear dependent rows after computation.
%  'Omega',val          If Omega is given, the vectors spanning U are tried to matched to the geometry of Omega.
%  'idx',val            Controls which transition matrix is taken to compute the eigenvalue. Default idx=1.
%  'v',vec              Overrides 'k'. Only vec is taken as starting vector. k must be empty.
%  'eig',val            Overrides 'k'. Only the eigenvector to eigenvalue val is taken as starting vectors. 'k' must be empty.
%
% Output:
% =======
%  U            Cell array of matrices of row vectors spanning U{k} corresponding to k
%  ev           Eigenvalue of corresponding starting vector. Is empty if 'v' is given
%  r            rank of U
%  errorflag    1 = Errors did probably occur. 0 = Probably no errors occured.
%
% Output: 
%   [{U1, U2, \sum_{w\in\Omega} <Ui,w>w, %, alpha_i, %, Omega}, {rank(U1),rank(U2)}]
%
% Info:
% =====
%   To compute vectors v one can try the following scheme:
%       v = sym( double(sym(T{2})^1000-sym(T{1})^1000) ); %where T is a cell-array of transition matrices
%       v = v(:,1); 
%   or
%       v = vpa( sym(T{2})^1000-sym(T{1})^1000, 100 ); %where T is a cell-array of transition matrices
%       v = v(:,1); 
%
% See also: transitionmatrix
%
% References:
% M. Charina, V. Yu.Protasov, 
% Regularity of anisotropic refinable functions
% Appl. Comput. Harm. A., (2017).
%
% Written by: tommsch, 2018

% XX ev not set yet
% XX zweites argument wird ein obligater parameter der definiert was ich ausrechnen will.
% XX    'u', val        U_i
% XX    scalar, val     U_1, U_2, etc.
% XX    'eig', val      corresponding to eigenvalue val
% XX simplify funktioniert nicht mehr
% XX Wenn es mehr als einen EW "1" gibt, dann muss ein EV ausgewählt werden, für den JSR(T|_U)<1 ist. Falls so einer existiert ist er eindeutig.
% XX Implement brute-force method

ev = 0; %%XX  % XX ev not set yet
errorflag = 0;
if( parsem('bruteforce',varargin) );  %use bruteforce method
    [U1,U2,D1,D2,al1,al2,Om] = computeU_bruteforce( varargin{1} ); 
    U = {U1,U2,D1,D2,al1,al2,Om}; 
    r = {rank(U1),rank(U2)}; 
    return; 
end;

T = varargin{1};
[symflag,varargin]      = parsem( 'sym', varargin );
[vpaflag,varargin]      = parsem( 'vpa', varargin, 0 );
[verbose,varargin]      = parsem( {'verbose'}, varargin, 1 );
[maxerror,varargin]     = parsem( 'maxerror', varargin, 1e-14 );
[nosimplify,varargin]   = parsem( 'nosimplify', varargin );
[Om,varargin]           = parsem( 'Omega', varargin, [] );
[idx,varargin]          = parsem( 'idx', varargin, 1 );
[vgiven,varargin]       = parsem( 'v', varargin, [] );
[eigenvalue,varargin]   = parsem( 'eig', varargin, [] );
if( size(varargin,2)<=1 ); 
    error( 'constructU: ''k'' is an obligatory parameter. If it is not needed for your purpose, give the empty argument [].' ); end;
k = varargin{2};

N = size( T, 2 ); %number of transition matrices
sizeOm = size( T{1}, 2 );



if( symflag || vpaflag ); %make matrices symbolic if wanted
    for i = 1:N; 
        T{i} = sym( T{i} ); end; end; 


if( ~isempty(vgiven) );
    if( ~isempty(k) ); 
        error( 'constructU: ''k'' must be empty.' ); end;
    v{1} = vgiven;
    if( symflag ); 
        v{1} = sym( v{1} ); end;    
    if( sizeOm>=1 && size(v{1},1)==1 );
        error( 'constructU: ''v'' must be a column vector.' ); end;
    sizek = 1;
else
    if( isempty(k) && isempty(eigenvalue) ); 
        k = 0:sizeOm-1; end; %if 'k' is empty, then take all possible eigenvectors
    if( ~isempty(eigenvalue) ); 
        if( ~isempty(k) ); 
            error( 'constructU: ''k'' must be empty.' ); end;
        k=-inf;  %some stupid scalar which should never occur.
    end;
    sizek=length(k);
    if(vpaflag || symflag); vprintf('Compute starting vector. ','imp',[2,verbose]); end; 
    v=cell(1,sizek);
    d=cell(1,sizek);
    p=cell(1,sizek);
    
    for i=1:sizek
        %compute all starting values corresponding to k
        if(symflag)
            [v{i},d{i},p{i}]=eig(T{idx});
            d{i}=d{i}(p{i},p{i});
        elseif(vpaflag); %compute the eigen vector
            [v{i},d{i}]=eig(vpa(T{idx},vpaflag));  %%XX sollte es nicht vpa(eig) sein?
        else
            [v{i},d{i}]=eig(T{idx});  
        end
        d{i}=diag(d{i});
        
        if(any(abs(d{i})>1+maxerror)); vprintf('constructU: Eigenvalue larger than 1 found','imp',[1,verbose]); end;
        val=sort(double(unique(d{i})),'descend');
        
        if(~isempty(eigenvalue)); val=eigenvalue; end;
        
        if(symflag);  %take all eigenvectors corresponding to the kth eigenvalue (eigenvalues larger than 1 are discarded)
            v{i}=v{i}(:,isAlways(abs(d{i}-val(i))<maxerror));  %vpa
        elseif(vpaflag)
            v{i}=v{i}(:,isAlways(abs(vpa(d{i}-val(i),vpaflag))<maxerror));  %vpa
        else
            v{i}=v{i}(:,abs(d{i}-val(i))<maxerror); %floating point
        end;
        
        %XX Teste ob Wert von k mit Wert des Eigenwertes zusammenpast (i.e. Fuer k==0 must val(i)==1 sein)
        if(isempty(v{i}));  vprintf('computeU: Given eigenvalue is no eigenvalue for the transition matrix %i.\n',idx,'imp',[1,verbose]); U=[]; r=[]; return; end;
        if(size(v{i},2)>1); vprintf('computeU: Eigenvalue is not unique. I take the first eigenvector:\n%v\n',v,'imp',[2,verbose]); end;
        vprintf('Eigenvector: \n%v\n',v{i},'imp',[2,verbose]); 
        if( symflag );
            warning('off','symbolic:sym:isAlways:TruthUnknown') %XX fix isAlways call
            if( ~all(isAlways(in(v{i},'rational'))) )
                vprintf( 'Symbolic computation failed to be exact.\n', 'cpr','err', 'imp',[1 verbose] );
                errorflag = 1; end;
            warning( 'on', 'symbolic:sym:isAlways:TruthUnknown' ); end;
        % XX Test if v/val is a eigenvector/value pair for all matrices
        
        %Test if v/val is really eigenvector/value pair for the matrix T{idx}
        if( vpaflag && vpa( norm(T{idx}*v{i}-val(i)*v{i},1), vpaflag )>maxerror  || ...
            norm( T{idx}*v{i}-val(i)*v{i}, 1 )>maxerror ); 
                vprintf( 'Could not find starting vector for eigenvalue %g. Err= %i > maxerror=%i .\n', val(i),  double(norm(T{idx}*v{i}-val(i)*v{i},1)), maxerror, 'cpr','err', 'imp',[1 verbose] );  
                vprintf( 'Use ''vpa'',val or simple floating point operations or increase <''maxerror'',val> .\n', 'cpr','err', 'imp',[1 verbose] );
                errorflag = 1;
                v{i} = zeros( sizeOm, 1 ); end; end;
end;

%make starting vectors
U = cell(1,sizek);
for i = 1:sizek
    U{i} = zeros(sizeOm,0);
    for j = 1:N
        if( vpaflag )
            W = vpa( T{j}*v{i}-v{i}, vpaflag );
        else
            W = T{j}*v{i}-v{i}; end;

        if( symflag );
            if(~isequal(W,sym(zeros(sizeOm,1)))); 
                U{i}=[U{i} W]; end;
        else
            if( norm(W,1)>=maxerror ); 
                U{i}=[U{i} W]; end; end; end;
    
    jj = 1;
    if( vpaflag || symflag ); 
        vprintf( 'Adding vectors: \n', 'imp',[1 verbose] ); end;
    while( true )
        RANK = rank(U{i});   
        RANKU = RANK;
        sizeU = size( U{i}, 2 );
        if( symflag || vpaflag ); 
            vprintf( [repmat('.',[1 (sizeU-jj+1)*N]) '\n'], 'imp',[1 verbose] ); end;
        for kk = jj:sizeU; %only work with the new vectors %%XX Sind indices (siehe auch XX unten) richtig?
            for iii = 1:N  %cycle through all transition matrices
                if( vpaflag )
                    W = vpa( T{iii}*U{i}(:,kk), vpaflag );
                else
                    W = T{iii}*U{i}(:,kk); end;
                if( rank([U{i} W])>RANKU );
                    U{i} = [U{i} W]; 
                    RANKU = rank( U{i} );
                    if( vpaflag || symflag ); vprintf( 'X', 'imp',[1 verbose] ); end; 
                else
                    if( vpaflag || symflag ); vprintf( '|', 'imp',[1 verbose] ); end; end; end; end;
        jj = kk+1;
        if( vpaflag || symflag ); 
            vprintf( '\n', 'imp',[1 verbose] ); end; 
        if( RANK==rank(U{i}) ); 
            break; end; end;
    
    if( symflag );
        for iii = 1:size(U{i},2);
            [~,DEN] = numden( U{i}(:,iii) ); %Make vectors to integers
            LCM = lcm( DEN );
            if( LCM>1e9 ); 
                vprintf( 'constructU: Computation could be wrong.\n', 'cpr','err', 'imp',[1 verbose] );
                vprintf( ' Use ''vpa'',val or simple floating point operations or increase <''maxerror'',val> .\n', 'cpr','err', 'imp',[1 verbose] );
                errorflag=1; end;
            if( all(isfinite(DEN)) ); 
                U{i}(:,iii) = lcm( DEN )*U{i}(:,iii); end; end; end;
    
    
    if( ~nosimplify );
        U{i} = intersectspace( U{i} ); %also simplifies bases

        %try to replace elements from U by elements from Omega
        if( ~isempty(Om) );
            SZE = size(U{i},2);            
            if( symflag || vpaflag ); 
                vprintf( ['Simplify: \n' repmat('.',[1 SZE]) '\n'], 'imp',[1 verbose] ); end;
            Vsimple = constructVt( Om, i-1 );
            for iii = 1:size(U{i},2)
                if( vpaflag || symflag ); 
                    vprintf( '|', 'imp',[1 verbose] ); end;
                for jjj = 1:size( Vsimple, 2 )
                    v = Vsimple(:,jjj);
                    Ux = U{i}(:,[1:iii-1, iii+1:SZE]);
                    if( rank([ v Ux]) == RANK && rank([ v U{i}]) == RANK ); 
                        U{i}(:,iii) = Vsimple(:,jjj);
                        break; end; end; end; end; end;
    
    for j = 1:size(U{i},2);
        U{i} = unique( U{i}.', 'rows' ).'; end;

    U{i} = double( U{i} ); %cast back to double if symbolic computation was 
    r = RANK;
    if( vpaflag || symflag ); 
        vprintf( '\n', 'imp',[1 verbose] ); end; 
end;


if(isscalar(k) || ~isempty(eigenvalue) || ~isempty(vgiven))
    U=U{1}; end;
end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   