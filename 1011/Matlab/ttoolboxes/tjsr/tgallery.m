function [ M ] = tgallery(varargin)
% [ val ] = tgallery( what, dim, N, [k], [options] )
% Returns sets of matrices, mostly used for tjsr.
%
% Input:
%   what            string which controls the return value
%   dim             dimension of the returned value (if applicable), 
%   N               number of matrices which shall be returned
%   k               argument used in some cases
%
% Options:
%   'sparse',val    Returns matrices with sparsity val
%   'pos'           Returns positive matrices
%   'rho'           Returns normalized matrices s.t. all matrices have spectral radius 1
%   'norm'          Returns normalized matrices s.t. all matrices have 2-norm 1
%   'sparse',val    Returns matrices with sparsity val
%   'int'           Returns integer matrices
%   'bool'          Returns matrices with values 0,1.
%
% Further Options:
%   'seed',val      Seed for random number generator. 
%                   If seed is set, the random number generator has the same state before and after that function 
%   'verbose',val   Verbose level
%   'nocell'        returns only the first matrix. Undefined behaviour if N>1.
%
% Output:    
%   val         Cell array of square matrices.
%
% Note:
% Most of the options preserve no special properties of the matrices. 
% For example: The entries in the matrices in tgallery('rand_gauss',10,2,'pos') do not have a Gaussian distribution anymore.
%   
%
% E.g.:  tgallery('rand_gauss',5,2,'seed',100,'rho')
%
% Written by: tommsch, 2018
% For more information write to: <a href="tommsch@gmx.at">tommsch@gmx.at</a>

%#ok<*ALIGN>
%#ok<*AGROW>

    [verbose,varargin] = parsem( {'verbose','v'},varargin,1);
    [seed,varargin]=parsem('seed',varargin,[]);
    [posflag,varargin]=parsem('pos',varargin);
    [rhoflag,varargin]=parsem('rho',varargin);
    [normflag,varargin]=parsem('norm',varargin);
    [sparsity,varargin]=parsem('sparse',varargin,0);
    [intflag,varargin]=parsem('int',varargin);
    [boolflag,varargin]=parsem('bool',varargin);
    [nocell,varargin]=parsem('nocell',varargin);
    
    if(~isempty(seed))
        OLDRNG=rng();
        if(isnumeric(seed))
            rng(mod(seed,2^30)); 
            rand(1000); %initialize random number generator
        else
            rng(seed); 
        end;
    end

    if(numel(varargin)>=4); 
        k=varargin{4}; 
    else; 
        k=0; 
    end;
    if(numel(varargin)>=3); 
        N=varargin{3}; 
    else; 
        N=1; 
    end;
    if(numel(varargin)>=2); 
        dim=varargin{2}; 
    else; 
        dim=2; 
    end;
    
    what=varargin{1};
    
    M=getmatrix(what, dim, N, k, verbose);
    
    J=numel(M);
    
    if(sparsity)
        idx=randperm(dim^2,round(dim^2*sparsity));
        for j=1:J
            M{j}(idx)=0;
        end
    end
    
    if(boolflag); 
        for j=1:J; M{j}=double(logical(round(M{j}))); end; end;
    if(intflag); 
        for j=1:J; M{j}=round(M{j}); end; end;    
    if(posflag); 
        for j=1:J; M{j}=abs(M{j}); end; end;
    
    if(rhoflag); 
        for j=1:J; M{j}=M{j}/rho(M{j}); end; end;
    if(normflag); 
        for j=1:J; M{j}=M{j}/norm(M{j}); end; end;
    if(intflag && rhoflag);  
        vprintf('''int'' and ''rho'' together does not return integer matrices.\n','imp',[0 verbose],'cpr','err'); end;
    if(intflag && normflag); 
        vprintf('''int'' and ''norm'' together does not return integer matrices.\n','imp',[0 verbose],'cpr','err'); end;
    if(rhoflag && normflag); 
        vprintf('''rho'' and ''norm'' does not return matrices with norm=rho=1.\n','imp',[0 verbose],'cpr','err'); end;
    
    
    if(~isempty(seed))
        rng(OLDRNG);
    end
    if(nocell)
        M=M{1};
    end
end

function [ val ] = getmatrix(what, dim, N, k, verbose)    

    %Random matrices, need two arguments 'rand_XXX',<dim>,<number of matrices>
    if(parsem({'rand_stochastic','r_stoch'},what)); %Stochastic matricesTransition matrices with mask length dim, dilation=N;
        val=matrix_stochastic(dim, N); 
        return; end;       
    if(parsem({'rand_doublestochastic','r_dstoch'},what)); %Double stochastic matrices
        val=matrix_doublestochastic(dim, N); 
        return; end;        
    if(parsem({'rand_stochastic_neg','r_stoch_n'},what)); %Stochastic matricesTransition matrices with mask length dim, dilation=N;
        val=matrix_stochastic_neg(dim, N); 
        return; end;        
    if(parsem({'rand_doublestochastic_neg','r_dstoch_n'},what)); %Double stochastic matrices
        val=matrix_doublestochastic_neg(dim, N); 
        return; end;       
    if(parsem({'rand_pm1','r_pm1'},what)); %random matrix with values -1, 0 1.
        for i=1:N; val{i}=randi(3,dim,dim)-2; end; 
        return; end;        
    if(parsem({'rand_bool','r_b'},what)); %random boolean matrix
        for i=1:N; val{i}=randi(2,dim,dim)-1; end; 
        return; end;        
    if(parsem({'rand_gauss','r_g'},what)); %random matrix with normally distributed values
        for i=1:N; val{i}=randn(dim); end; 
        return; end;        
    if(parsem({'rand_equal','rand','r_e','r'},what)); %random random matrix with equally distributed values in [0 1]          
        for i=1:N; val{i}=rand(dim); end; 
        return; end;        
    if(parsem({'rand_equal_neg','rand_neg','r_e_n','r_n'},what)); %random random matrix with equally distributed values in [-1 1] 
        for i=1:N; val{i}=2*rand(dim)-.5; end; 
        return; end;        
    if(parsem({'rand_unitary','r_u'},what)); %Unitary matrices 
        val=matrix_unitary(dim, N); 
        return; end;        
    if(parsem({'rand_zero','r_z'},what)); %matrices with spectral radius equal zero             
        val=matrix_zerospectralradius(dim, N); 
        return; end;        
    if(parsem({'rand_colu_1','r_colu_1'},what)); %random n-by-n matrix with columns of unit 2-norm, with random singular values whose squares are from a uniform distribution, contains no zeros   
        val{1}=gallery('randcolu',dim,dim,0); 
        val{2}=gallery('randcolu',dim,dim,0); 
        return; end; 
    if(parsem({'rand_colu_0','r_colu_0'},what)); %random n-by-n matrix with columns of unit 2-norm, with random singular values whose squares are from a uniform distribution, may contains zeros   
        val{1}=gallery('randcolu',dim,dim,1); 
        val{2}=gallery('randcolu',dim,dim,1); 
        return; end; 
    if(parsem({'rand_corr_1','r_corr_1'},what)); %random n-by-n correlation matrix with random eigenvalues from a uniform distribution, contains no zeros. Makes maximal trees for high dimensions   
        val{1}=gallery('randcorr',dim,dim,0); 
        val{2}=gallery('randcolu',dim,dim,0); 
        return; end; 
    if(parsem({'rand_corr_0','r_corr_0'},what)); %random n-by-n correlation matrix with random eigenvalues from a uniform distribution, may contains zeros. Makes maximal trees for high dimensions
        val{1}=gallery('randcorr',dim,dim,1); 
        val{2}=gallery('randcolu',dim,dim,1); 
        return; end; 
    if(parsem({'rand_hess','r_hess'},what));  %random n-by-n real, orthogonal upper Hessenberg matrix.
        val{1}=gallery('randhess',dim); 
        val{2}=gallery('randhess',dim); 
        return; end;
    if(parsem({'rand_o_1','r_o_1'},what));          
        val{1}=gallery('rando',dim,1); 
        val{2}=gallery('rando',dim,1); 
        return; end; %random 0-1 matrix. Interesting for balancing.
    if(parsem({'rand_o_2','r_o_2'},what));  %random 0-1 matrix. Interesting for balancing.        
        val{1}=gallery('rando',dim,1); 
        val{2}=gallery('rando',dim,2); 
        return; end; 
    if(parsem({'rand_o_3','r_o_3'},what));  %random 0-1 matrix. Interesting for balancing.
        val{1}=gallery('rando',dim,1); 
        val{2}=gallery('rando',dim,3); 
        return; end; 
    
    %Experimental
    if(parsem({'rand_TU','r_TU'},what)); %Transition matrices with mask length dim, dilation=N; dimension is implicitely defined via dilation            
        val=matrix_TU(dim, N); 
        return; end;        
    if(parsem({'rand_TV0','r_TV0'},what)); %Transition matrices with mask length dim, dilation=N;           
        val=matrix_TV0(dim, N); 
        return; end;        

    %Matrices from matlab gallery
    if(parsem({'orthog_1','o_1'},what)); %complex eigenvalues / maximum number of edges
        val{1}=gallery('orthog',dim,1); 
        val{2}=gallery('orthog',dim,2); 
        return; end; 
    if(parsem({'orthog_2','o_2'},what)); %real eigenvalues, but complex dual eigenvectors
        val{1}=gallery('orthog',4,4); 
        val{2}=circshift(val{1},1); 
        return; end; 
    if(parsem({'orthog_3','o_3'},what)); %Polytope is a sphere
        val{1}=gallery('orthog',3,5); 
        val{2}=circshift(val{1},1); 
        return; end; 

    %Matrices from Examples/Papers
    if(parsem('cex',what)); %jsr finitesness conjecture counterexample      
        val{1}=sym([1 1; 0 1]);  
        val{2}=sym('0.749326546330367557943961948091344672091327370236064317358024')*[1 0; 1 1]; 
        return; end;
    if(parsem('cex2',what)); %jsr finitesness conjecture counterexample
        val{1}=[1 1; 0 1]; 
        val{2}=0.7493265463303675579439619*[1 0; 1 1]; 
        return; end;
    if(parsem('cex3',what)); %jsr finitesness conjecture counterexample
        val{1}=[1 1; 0 1]; 
        val{2}=0.75*[1 0; 1 1]; 
        return; end;
    if(parsem('1',what)); %program works
        val{1}=[-0.0766   -0.7779   -0.4523; 0.5958    0.7105    0.7678; -0.8199    0.6780   -0.5047];
        val{2}=[0.2923    0.7024    0.6463; 0.7111   -0.1468   -0.1636; -0.0213   -0.2592    0.2707];
        val{3}=[0.6170   -0.2134   -0.3799; -0.7410   -0.4941    0.6517; 0.3961   -0.7559    0.5435]; 
        return; end;
    if(parsem('code',what));    
        val=codecapacity(dim,'verbose',verbose); 
        return; end;   
    if(parsem('euler',what));    
        val=matrix_euler(dim); 
        return; end;   
    if(parsem('ones',what)); %matrices of ones - should be the same as euler. But this is not tested yet    
        val=transitionmatrix({ones(1,dim)',N},'Omega',0:dim-1,'nocheck');
        val=cellfun(@transpose,val,'uniformoutput',0); 
        return; end;       
    if(parsem('daub',what)); %Transition matrices of difference scheme for Daubechies wavelets.         
        val=daubechiesmatrix(dim); 
        return; end; 
    if(parsem('binary',what)); %Binary matrices in ascending order  
        val=binarymatrix(dim, N, k, 0); 
        return; end; 
    if(parsem('binary2',what)); %Binary matrices in ascending order  
        val=binarymatrix(dim, N, k, 1); 
        return; end; 
    if(parsem('t',what)); %Example from my thesis.       
        val{1}=[2 -1; 2 0]; val{2}=[0 -1; 2 2]; 
        return; end; 
    if(parsem('DD8',what)); %Dubuc Deslauriers 8-point scheme  %Prot 2016 p29
        a=[-5    30   -56   -14   154   -14   -56    30    -5]; 
        val=transitionmatrix({a',2,[0 1]},'Omega',0:7); 
        val=cellfun(@transpose,val,'uniformoutput',0); 
        return; end; 

    if(parsem('morris_p3',what)); %A RAPIDLY-CONVERGING LOWER BOUND FOR THE JOINT SPECTRAL RADIUS VIA MULTIPLICATIVE ERGODIC THEORY, IAN D. MORRIS
        val{1}=[2 2; 0 0]; 
        val{2}=[1 1; 1 1]; 
        return; end; 
    if(parsem('prot2012_p40',what));%pascal rhombus %jsr=2; Prot2012, p40 %program works
        val{1}=[0 1 0 0 0; 1 0 2 0 0; 0 0 0 0 0; 0 1 0 0 1; 0 0 0 2 1]; 
        val{2}=[1 0 2 0 0; 0 0 0 2 1; 1 1 0 0 0; 0 0 0 0 0; 0 1 0 0 0]; 
        return; end 
    if(parsem('prot2012_p35',what)); %Prot2012, p35 %program works
        val{1}=[1 1; 0 1]; 
        val{2}=9/10*[1 0; 1 1]; 
        return; end;
    if(parsem('prot2012_p43',what)); %Euler binary problem for r=7, %Prot12 p43
        val{1}=[1 1 1 1 0 0;0 1 1 1 0 0; 0 1 1 1 1 0; 0 0 1 1 1 0; 0 0 1 1 1 1; 0 0 0 1 1 1]; 
        val{2}=[1 1 1 0 0 0; 1 1 1 1 0 0;0 1 1 1 0 0; 0 1 1 1 1 0; 0 0 1 1 1 0; 0 0 1 1 1 1]; 
        return; end; 
    if(parsem('prot2012_p44',what)); %Euler ternary problem %Prot12 p44, %algorithm working %tree looks different then in paper
        val{1}=[1 1 1 1 1 0 0; 0 1 1 1 1 0 0; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 0 0 1 1 1 1 0; 0 0 1 1 1 1 1; 0 0 1 1 1 1 1];
        val{2}=[1 1 1 1 1 0 0; 1 1 1 1 1 0 0; 0 1 1 1 1 0 0; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 0 0 1 1 1 1 0; 0 0 1 1 1 1 1];
        val{3}=[1 1 1 1 0 0 0; 1 1 1 1 1 0 0; 1 1 1 1 1 0 0; 0 1 1 1 1 0 0; 0 1 1 1 1 1 0; 0 1 1 1 1 1 0; 0 0 1 1 1 1 0];
    return; 
    end; 

    if(parsem('grip',what)); %Gripenberg
        val{1}=[0 1; 0 0]; 
        val{2}=1/5*[0 0; 1 0]; 
        return; end; 
    if(parsem('grip_p52',what));  %Gripenberg, p52 %long smp=1111111111112 %sic
        val{1}=1/5*[3 0; 1 3]; 
        val{2}=1/5*[3 -3; 0 -1]; 
        return; end;
    if(parsem('prot2016_p35',what)); %Prot2016 p50 p35 p33 %working %matrices have common invariant subspace, algorithm still works
        val{1}=1/48*[48     0     0   -12     0     0;     0     0     0     0    12     0;     0     0     0     0     0    12;     0     0     0   -12     0     0;     0    12     0     0    36     0;     0     0    12     0     0    36;];
        val{2}=1/48*[0    12     0     4     9    -4;    48    36    36   -92    -9    11;     0     0   -12   100     0   -25;     0     0     0    20     0    -8;     0     0   -48    80    48   -20;     0     0     0  -112     0    16];
        val{3}=1/48*[0     0    12     4    -4     9;     0   -12     0   100   -25     0;    48    36    36   -92    11    -9;     0     0     0    20    -8     0;     0     0     0  -112    16     0;     0   -48     0    80   -20    48];
        val{4}=1/48*[0    -4    -4   -12     5     5;     0    80  -100     0    80  -103;     0  -100    80     0  -103    80;    48    16    16   -12    16    16;     0   -80   112     0   -92   100;     0   112   -80     0   100   -92];
        return; end;
    if(strcmp(what,'mejstrik_119_long')); %smp length of 119!. norms converge very wiggly to the spectral radius. Experimentally found numbers. 
        val{1}=[0.163026496094203    -0.92406403398077952;    0.94910922519637453 0.75424733784141262]; 
        val{2}=[-0.95851077103541249    -0.65295899115856648;    0.6731845518678724   -0.58469670981727162]; 
        return; end; 
    if(parsem('mejstrik_119',what)); %smp length of 119!. norms converge very wiggly to the spectral radius. Rational approximations of the long numbers above.
        val{1}=[15/92 -73/79;56/59 89/118]; 
        val{2}=[-231/241 -143/219; 103/153 -38/65]; 
        return; end; 
    if(parsem('mejstrik_Cn',what)); %smp length is equal to dim_2ndargument
        val{1}=[1 1; 0 1]; 
        val{2}=[0 0; 1/dim*exp(1+1/dim) 0]; 
        return; end; 
    if(parsem('mejstrik_longsmp',what)); %smp length goes to infinity as dim_2ndargument goes to zero.
        val{1}=[1 1; 0 1]; 
        val{2}=[0 0; dim 0]; 
        return; end; 
    
    if(parsem('gwz05_ex61',what)); %Guglielmi, Wirth, Zennaro - 2005 - Complex polytope extramality results for families of matrices
        val{1}=[cos(sym(1)) -sin(sym(1)); sin(sym(1)) cos(sym(1))]; 
        val{2}=[1/sym(2) -1/2;0 0]; 
        return; end; 
    if(parsem('gwz05_ex62',what)); %Guglielmi, Wirth, Zennaro - 2005 - Complex polytope extramality results for families of matrices
        val{1}=(3-sqrt(sym(5)))*[2 1;1 1]; 
        val{2}=(3-sqrt(sym(5)))*[1 1;1 2]; 
        return; end; 
    if(parsem('gwz05_ex63',what)); %Guglielmi, Wirth, Zennaro - 2005 - Complex polytope extramality results for families of matrices
        val{1}=[1 0 0;0 1/2 0;0 0 1/4]; 
        val{2}=[1/2 0 0;1/2 0 0;1/2 0 0]; 
        return; end; 
    if(parsem('gwz05_ex64',what)); %Guglielmi, Wirth, Zennaro - 2005 - Complex polytope extramality results for families of matrices
        val{1}=[1 1;0 0]; 
        val{2}=[0 0;1 1]; 
        val{3}=1/2*[1 1;1 1]; 
        val{4}=2/3*[1 0;-1 0]; 
        return; end; 
    if(parsem('gz09_ex31',what)); %Guglielmi, Zennaro - Finding Extremal Complex Polytope Norms for Families of Real Matrices 
        val{1}=[cos(sym(1)) -sin(sym(1)); sin(sym(1)) cos(sym(1))]; 
        val{2}=dim/sqrt(sym(2))*[1 1;0 0]; 
        return; end; 
    if(parsem('gz09_ex32',what)); %Guglielmi, Zennaro - Finding Extremal Complex Polytope Norms for Families of Real Matrices
        val{1}=[-1/2 -sqrt(sym(3))/2; sqrt(sym(3))/2 -1/2]; 
        val{2}=[11/20 11/20 -11/20 -11/20]; 
        return; end; 
    if(parsem('gz09_ex41',what)); %Guglielmi, Zennaro - Finding Extremal Complex Polytope Norms for Families of Real Matrices
        val{1}=[-3 -2 1 2;-2 0 -2 1; 1 3 -1 -5;-3 -3 -2 -1]; 
        val{2}=[1 0 -3 -1;-4 -2 -1 -4; -1 0 -1 2;-1 -2 -1 2]; 
        return; end; 
    if(parsem('mej20_ex31',what)); 
        val{1}=[0 0;1 1];
        val{2}=[1 1;0 1]; 
        return; end; 
    
    %Bad examples where something goes wrong
    
    if(parsem('nondominant',what)); %matrix family has no dominant smp. Algorithm must not converge.
        val{1}=[1 1; 0 1]; 
        val{2}=4/5*[1 0; 1 1]; return; end; 
    if(parsem('nondominant_limit',what)); %matrix family from nondominant.and limit matrix of sequence approaching smp-spectral radius is also included.
        val{1}=[1 1; 0 1]; 
        val{2}=4/5*[1 0; 1 1]; 
        val{3}=[(sqrt(5)+1)/4 1/2; (sqrt(5)-1)/4 (3-sqrt(5))/4]; 
        return; end; 
    if(parsem('multv0',what)); %M2^2=1/3*I
        val{1}=1/4*[2 -1; 0 2]; 
        val{2}=1/3*[1 2; 1 -1]; 
        return; end; 
    if(parsem('5',what)); % _very_ long computation duration
        a=((ones(1,15)*-1).^(1:15).*floor(sqrt(1:15)))'; 
        a(10)=-2;    
        %a(10)=-2.11    
        %a(10)=-2.237884497219498; (smallest?? value for which the matrix has real eigenvector
        val=transitionmatrix(getS({a,2},'nocheck'),'Omega',0:13);
        return; end;
    
    


    %Other matrices
    
    if(parsem('nearlycomplex',what)); %two matrices who are nearly the complex case
        val=nearlycomplex(); 
        return; end; 
    if(parsem('2',what)); %program works
        val{1}=[1 -1; 3 -2]; 
        val{2}=[1 3; -1 -1]; 
        return; end;  
    if(parsem('3',what)); %program works
        val{1}=[1 4; 4 1]; 
        val{2}=[1 1; 0 -1]; 
        val{3}=[1 0; 1 -1]; 
        return; end;  
    if(parsem('lotsmps',what));  %a huge amount of smp's %not working
        val=matrix_3(1, -2, 2 ); 
        return; end;
    if(parsem('complex',what)); %complex eigenvector
        a=(ones(15,1)*-1).^(1:15).*floor(sqrt(1:15)); a(10)=-3; 
        val=transitionmatrix({a,2,[0 1]},'Omega',0:13)'; 
        return; end;  
    if(parsem('8',what)); %very long smp
        val{1}=[1 1 -1 1;1 0 1 0;-1 0 1 0;-1 -1 0 0];
        val{1}=val{1}/trho(val{1}); 
        val{2}=[1 1 0 0;1 1 -1 0;-1 1 0 -1; 1 -1 0 0];val{2}=val{2}/trho(val{2}); 
        return; end; 
    if(parsem('9',what)); %same eigenvectors, every product is smp    
        val{1}=[1 1; 0 2]; 
        val{2}=[2 0; 1 1]; 
        return; end; 
    
    if(parsem('psp',what)); %diffscheme for pseudospline(2,3,2) %program works
        a=1/256*[3 -18 38 -18 3]; 
        val=transitionmatrix({a,2,[0 1]},'Omega',0:3); 
        val=cellfun(@transpose,val,'uniformoutput',0); 
        return; end; 

    vprintf('Error: No matrices returned.','cpr','err'); %when 'what' contains something wrong
    val={};
end

function [T] = matrix_3(c0,c1,c2);
    T{1}=[c0 0; c2 c1];
    T{2}=[c1 c0; 0 c2];
end

function [T] = matrix_euler(dim)
    %length of smp==2 if dim has even number of ones in binary expansion
    %length of smp==1 if dim has odd number of ones in binary expansion
    T{1}=zeros(dim,dim);
    T{2}=T{1};
    for i=1:dim
        for j=1:dim
            T{1}(i,j)=tif(1<= 2*j-i && 2*j-i <= dim+1, 1, 0);
            T{2}(i,j)=tif(0<= 2*j-i && 2*j-i <= dim, 1, 0);
        end
    end
end

function [TU] = matrix_TU(l, M)
    dim=size(M,2);
    while(true);
        a=randi(11,l*ones(1,dim))-6;
        a=sequence(a);
        S=getS('a',a,'M',M,'verbose',0);
        S=normalizeS(S,'verbose',1,'gauss');         
        if(~isempty(S)); 
            break; end;
    end
    T=transitionmatrix(S);
    U=constructU(T,0);
    TU=restrictmatrix(T,U);            
end

function [TV0] = matrix_TV0(l, M)
    dim=size(M,2);
    while(true);
        a=randi(11,l*ones(1,dim))-6;
        a=sequence(a);
        S=getS('a',a,'M',M);
        S=normalizeS(S,'verbose',1,'scale');         
        if(~isempty(S)); 
            break; end;
    end
    [T,Om]=transitionmatrix(S);
    V0=constructV(Om,0);
    TV0=restrictmatrix(T,V0);            
end

function [T] = matrix_stochastic(dim, N)
    T=cell(1,N);
    for i=1:N
        T{i}=normalizematrix(rand(dim),'colsum');
    end
end

function [T] = matrix_doublestochastic(dim, N)
    T=cell(1,N);
    for i=1:N
        T{i} = zeros(dim,dim);
        c=rand(1,dim);c=c/sum(c);
        for k = 1:dim
            [~,idx] = sort(randperm(dim)); 
            idx = idx + (0:dim-1)*dim;     
            T{i}(idx) = T{i}(idx) + c(k);
        end
    end
end
function [T] = matrix_stochastic_neg(dim, N)
    T=cell(1,N);
    for i=1:N
        T{i}=normalizematrix(rand(dim)-.5,'colsum');
    end
end
function [T] = matrix_doublestochastic_neg(dim, N)
    T=cell(1,N);
    for i=1:N
        T{i} = zeros(dim,dim);
        c=rand(1,dim)-.5;c=c/sum(c);
        for k = 1:dim
            [~,idx] = sort(randperm(dim)); 
            idx = idx + (0:dim-1)*dim;     
            T{i}(idx) = T{i}(idx) + c(k);
        end
    end
end


function val = matrix_unitary(dim,J)
    % this function generates a random unitary matrix of order 'n' and verifies
    val=cell(1,J);
    for j=1:J
       X = randn(dim); % generate a random matrix 
       [Q,R] = qr(X); % factorize the matrix
       R = diag(diag(R)./abs(diag(R)));
       val{j} = Q*R; % unitary matrix         
    end
end

function val = matrix_zerospectralradius(dim,J)

if(dim==0); 
    val=repcellm(0,[1 J]); 
    return; end;

    % this function generates a random unitary matrix of order 'n' and verifies
%     val=cell(1,J);
%     for j=1:J
%         v=randn(dim,1);
%         w=randn(dim,1);
%         %U=matrix_unitary(dim,1); U=U{1};
%         A=randn(dim);
%         w(end)=0;
%         v(end)=1;
%         w(end)=-v'*w;
%         val{j}=A\v*w'*A;
%         %val{j}=U'*v*w'*U;
%         %trho(val{j});
%     end
    
val=cell(1,J);   
j=1;
N=3;
while(true)
    A=randi(N,dim)-floor(N/2);
    while(rho(A)~=0)
        idx=find(A);
        A(idx(randi(numel(idx))))=0;
    end
    if(anym(A)~=0);
        val{j}=A;
        j=j+1;
    end
    if(j==J+1); 
        break; end;
end
    
    %%This is slower, and does not have better numerical properties
%       val=cell(1,J);
%     jordan=[[zeros(dim-1,1) eye(dim-1) ];zeros(1,dim)];
%     for j=1:J
%        U=matrix_unitary(dim,1); U=U{1};
%        val{j}=U*jordan*U';
%        trho(val{j})
%     end
end

function val = nearlycomplex()
   val{1}=[
        0.260779273117590   0.321746900633509   0.430143562965100
        0.332728106032852   0.165237727924909   0.246741880562854
        0.435544012664021   0.460173223125738   0.331865771472727
        ];
   
   val{2}=[
        0.682274614607416  -0.373006400071818   0.529326054171797
        0.741899861356640   0.551196958575390  -0.263494595696068
       -0.371216501721239   1.002783855858494   0.524985866808493
       ];
end



function comment; %#ok<DEFNU>

    tjsr(tgallery('lotsmps'),'verbose',1,'notestsmp','nobalancing','maxsmpdepth',2) %not working

    tjsr(tgallery('5'),'verbose','maxsmpdepth',15) %very bad estimate of the JSR

    tjsr(tgallery('ones',2,2),'verbose',1,'proof','clc','maxsmpdepth',5,'balancingvector',ones(1,14),'nosmallrank','noinvariantsubspaces') %not working
    tjsr(tgallery('ones',2,2),'verbose',1,'proof','clc','maxsmpdepth',15,'nobalancing','nosmallrank','noinvariantsubspaces') %not working
    tjsr(tgallery('ones',15,2),'verbose',1,'proof','clc')
    tjsr(tgallery('ones',14,2),'verbose',1,'proof','clc')

    %example where epspolytope works wonder
    tjsr(tgallery('ani',5),'maxsmpdepth',1,'verbose',2,'invariantsubspace','none')
    tjsr(tgallery('vale'),'maxsmpdepth',2,'autoextravertex',1,'plot','norm','clc','verbose',2)

    %nice bad example, funktioniert in Matlab2016a, aber nicht mehr in Matlab2017a
    %funktioniert in L'Aquila auch mit Matlab2017a
    M=[2 1; 0 3]; a0=[1 0 0; 1 1 1; 0 1 1];D=[     0     1     1     1     2     2;     0     0     1     2     1     2 ];
    S=getS({1/6*convm(a0,a0),M,D},'bigcheck');
    [T,~]=transitionmatrix(S);
    U2=double([1 0 0 -1 0 0 0 0; 0 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0 -1; 1 0 -2 0 1 0 0 0; 0 1 0 -2 0 1 0 0; 0 0 1 0 -2 0 1 0; 0 0 0 1 0 -2 0 1]');
    TU2=restrictmatrix(T,U2,'outputdouble'); 
    tjsr(TU2,'clc','verbose',2,'autoextravertex',1e12)

    
    getdefault('rand_rho',15,2,100) %Beispiel wo simplepolytope etwas bewirkt

    tgallery('rand_rho',2,2,1000008); %langer Kandidat: 24
    tgallery('rand_rho',2,2,1000021); %Kandidate hat Laenge 119
    tjsr(tgallery('rand_rho',10,2,1000023),'plot','tree','verbose',1)
    tjsr(tgallery('cex2'),'plot','tree','verbose',2,'autoextravertex',0,'clc')
    tjsr({[2 1; 0 -2],[2 0; -1 -2]},'invariantsubspace','none','maxsmpdepth',2); %gut um selectordering zu testen
    {[2 1; 0 2],[2 0; -1 2]}; %#ok<VUNUS> %alle orderings mit Laenge kleiner gleich 7 haben rho==2, darüber steigt er an auf 2.36007056204988
    tjsr(tgallery('ani',5),'maxsmpdepth',1,'verbose',2,'invariantsubspace','none','plot','norm','autoextravertex',0,'clc','waitafterbalancing'); %Multivariates Subdivision Bsp das nur mit Balancing funktioniert. Gegenbeispiel zur Vermutung Protasovs.


    %tjsr(tgallery('nondominantsmp'),'addtranspose','plot','polytope','nobalancing','verbose',2,'simplepolytope',1e-10,'testeigenplane',0,'testspectralradius',0,'findsmp_N',1000,'maxsmpdepth',500)
    tjsr(tgallery('vale'),'maxsmpdepth',1,'autoextravertex',1,'plot','norm','clc','verbose',2,'simplepolytope',1,'clc','balancingdepth',3,'ordering',{[1]},'waitafterbalancing',1)
    
    %norm estimate stimmt nicht
    tjsr(tgallery('rand_gauss',50,2,'seed',102,'pos','rho','sparse',.95))
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 