function [ T ,Om, Vt] = transitionmatrix(varargin)
% [ T, Om, Vt ] = transitionmatrix( S, [options] )
% Constructs transition matrices.
% T_d(al,be)=(M.al-be+d)_{al,be\in\Omega}, d\in\ZZ^s/M\ZZ^s
%
% Input:  
%   S       cell array of subdivision operators
%
% Options: 
%  'Omega',val          Where <val>=dimxN integer matrix, index set for constructing the transition matrices. 
%                       If no set is given, a set is constructed by calling constructOmega()
%  'V',val              (experimental) constructs set Om such that Xmu is non-empty for all abs(mu)<=val+1, see Mejstrik, PhD Thesis, 2019
%                                      This option is needed if one wants to compute whether a subdivision scheme is C^val or not
%  'onlyindex'          Function instead returns (dim+#Omega)x#Omega matrix with the indices (instead of the values of the mask a at the indices position) entries
%  'noflat'             Function instead returns cell array of cell arrays, each corresponding to one subdivision operator
%  'colsum',val         Tests if the columns sum up to val. If val==zero, then the sum of the first column in the first transition matrix is taken as value.
%  'verbose',val        Verbose level
%
% Output:
%   T           cell array of transition matrices {T_(1,d1),...,T_(1,dj_1), T_(2,d1), T_(2,d2), ... , T_(n,dj_n)}
%   Omega       a set such that \ell(Omega) if invariant for all matrices in T
%   Vt          space of span of differences of shifts of the delta sequence with support in Omega
%               Only returned if 'V' is given.
%
% E.g.: vdisp(transitionmatrix({1/4*[1 4 3]',2}))
%
% See also: getS, constructOmega, restrictmatrix, constructV, constructVt, constructU
%
% Written by: tommsch, 2018

 %#ok<*ALIGN>



    [verbose,varargin]=parsem({'verbose','v'},varargin,1);
    [noflat,varargin]=parsem('noflat',varargin);
    [Om,varargin]=parsem({'Omega','Om','O'},varargin,[]);
    [Vval,varargin]=parsem('V',varargin,[]);
    [colsum,varargin]=parsem({'colsum','col','c'},varargin,[]);
    [onlyindex,varargin]=parsem({'onlyindex','oi'},varargin);

    if(isempty(Om))
        vprintf('You should specify a point which is for sure inside of Omega via <''Omega'',vector>.','imp',[2,verbose]); 
        Om=constructOmega(varargin{:},'V',Vval); end;



    S=getS(varargin{1});
    varargin(1)=[];

    parsem(varargin,'test');
    dim=size(S{1,2},1);    

    if(size(Om,1)~=dim); 
        error( 'transitionmatrix:dimOm', 'Dimension of set Omega is wrong.' ); end;

    sizeS=size(S,1);

    T=cell(sizeS,1);

    vprintf('Omega=\n%v\n',characteristic(Om),'imp',[2,verbose]); 

    for j=1:sizeS
        a=S{j,1}; M=S{j,2};D=S{j,3};
         T{j}=cell(1,size(D,2));  %initialize T

        if(size(D,1)~=dim); 
            error( 'transitionmatrix:dimD', 'Dimension of digit-set D is wrong'); end;

        vprintf('Construct Set %i:\n%v\n%v\n%v',j,a,M,D,'imp',[2,verbose]); 

        for i=1:size(D,2) %i: index of digit

            d=D(:,i);
            if(~onlyindex);
                T{j}{i}=zeros(size(Om,2));
            else
                T{j}{i}=zeros(size(Om,2)*dim,size(Om,2)); end;
            for alpha=1:size(Om,2)
                al=Om(:,alpha);
                for beta=1:size(Om,2)
                    be=Om(:,beta);
                    index=M*al-be+d;
                    if(onlyindex); 
                        T{j}{i}(dim*alpha-dim+1:dim*alpha, beta)=index; 
                    else; 
                        T{j}{i}(alpha,beta)=a.ref(index); end; end; end; end; end;


    %Test transition matrices for correct one-eigenvector
    if(~onlyindex && ~isempty(colsum))
        errorprinted=0;
        if(colsum==0); 
            colsum=sum(T{1}{1}(:,1)); end; %take default value for colsum
        for j=1:sizeS
            for i=1:size(S{j,3},2)
                wrong=0;
                s=sum(T{j}{i});
                if(any(abs(s-colsum)>eps)); 
                    wrong=1; end;
                if(wrong);
                    if(errorprinted==0); 
                        warning( 'transitionmatrix:colsum', 'Error in Transition Matrix (Column sum is wrong)' ); end;
                    errorprinted=1;
                    vprintf('%i/%i ',j,i,'cpr','err','imp',[1,verbose]); end; end; end;
        if(errorprinted==1); 
            vprintf('\n','imp',[1 verbose]); end; end;

    vprintf('\n%v\n',characteristic(Om),'imp',[2,verbose]);
    if(~noflat); 
        T=[T{:}]; end;

    if(~isempty(Vval) && nargin>=3)
        Vt=constructVt(Om,Vval);
        V =constructV(Om,Vval);
        if(size(V,2)~=size(Vt,2));
            warning( 'transitionmatrix:VtV', 'Vt(Omega)~=V(Omega). Vt(Omega) may not be invariant under all transition matrices.' ); end; end;
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 