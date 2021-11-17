function [Sp, mu]=differencescheme(varargin)
% [Sp, mu] = difference scheme(S, k, [options])
% Computes the masks of difference schemes (of order k)
% !! This function is not well tested and should not be used !!
%
% Input:
%   S               Cell array of subdivision schemes (as returned by 'getS') (value of S is passed to getS, return value is used.)
%   k               integer, order of wanted difference schemes
%
% Options:
%   'verbose',val   Verbose level
%   'sym'           Computation is done partly symbolically. Default: numerically
%   'pinv'          Computation uses the pinv function. Default: slash-operator
%   'polyonimal'    (deprecated) Uses polynomial-long division to find difference scheme. Yields wrong results sometimes. No options available for this algorithm.
%   'linear'        (deprecated) Solves a linear system and has some bugs.
%
% Output:
%   Sp              Cell array of vector valued subdivision schemes. The masks are a cell array which contain the difference-masks
%   mu              Ordering of the differences used to compute B

% XX Function untested
% XX 'mu' als option hinzufügen
% XX 'normalize' machen
% XX zweite Art der Vereinfachung bei der Konstruktion des Difference-schemes: b12=b21'


if(parsem('polynomial',varargin));
    Sp=differencescheme_polynomials(varargin{:});
elseif(parsem('linear',varargin));    
    Sp=differencescheme_linearsystem(varargin{:});
else %default
    Sp=differencescheme_mixed(varargin{:});
end

end

function Sp = differencescheme_mixed(varargin)

if(size(varargin,2)<=1); error('differencescheme: Too less arguments. Order of differences ''k'' is missing.'); end;
verbose=parsem({'verbose','v'},varargin,1);
symflag=parsem('sym',varargin);
pinvflag=parsem('pinv',varargin);
extrasize=parsem('extrasize',varargin,0);
S=varargin{1}; S=getS(S);
k=varargin{2}; %order of difference-scheme
J=size(S,1);
if(J>1 && isempty(k) || J>1 && length(k)>1); error('differencescheme: ''k'' must be an integer.'); end;
if(any(k<0)); error('differencescheme: ''k'' must be nonnegative.'); end;
dim=size(S{1,2},1);

Sp=S; %initialize Sp


    for j=1:J 
        
        Sj=S(j,:);
        if(~isequal(size(Sj{1}),[1 1])); error('differencescheme: Works only for scalar-valued schemes.'); end;
        
        lhs=diffsequ(Sj{1},k); %lhs of linear system
        lhs=lhs.equalize;
        
        n=sequence([1;0]); n.dim=dim;
        n=diffsequ(n,k); %nabla
        nup=upsample(n,Sj{2});
        dim_outer=sizem(lhs,[],1);

        B=cell(dim_outer,dim_outer); %difference masks
        
        for i=1:numel(B);
            B{i}=sequence(sym(['B' num2str(i) '_'], size(lhs.M{1}.c)+extrasize),lhs.M{1}.idx);
        end;
        B=svm([dim_outer dim_outer],B{:});
        
        rhs=symbol2mask(B.symbol*nup.symbol,'var','z1','dim',dim,'sym'); %rhs of linear system
        val=[rhs{:}];
        SYMS=symvar(val); %list of symbols
        
        

        EQU=cell(size(lhs.M));
        for i=1:lhs.numel
            val1=lhs.M{i}.c;
            val2=rhs{i};
            szelhs=sizem(val1,[],dim).'; szerhs=sizem(val2,[],dim).'; szemax=max(szelhs,szerhs);
            if(any(szelhs~=szemax)); val1=subsco(val1,szemax,0); end; %make matrices the same size
            if(any(szerhs~=szemax)); val2=subsco(val2,szemax,0); end;
            EQU{i}=val1==val2; %generate symbolic equations
            %EQU{i}=lhs.M{i}.c==svm(size(lhs),rhs{:}); %generate symbolic equations
            
        end
        EQU=[EQU{:}]; EQU=EQU(:);
        
        
        warning('off','MATLAB:rankDeficientMatrix');
        warning('off','symbolic:mldivide:RankDeficientSystem');
        warning('off','symbolic:mldivide:InconsistentSystem');
        for tries=1:3
            if(tries==1);     BADD=nondiag(B);  %try all non-diagonal blocks are zero
            elseif(tries==2); BADD=triu(B,1);   %try all lower upper blocks are zero
            else;             BADD=[];          %everything can be nonzero
            end
            ADD=cell(1,numel(BADD)); 
            for i=1:numel(BADD); 
                val=sym(BADD.M{i}.c)==zeros(size(BADD.M{i})); 
                val=val(:);
                ADD{i}=val;
            end; 
            ADD=vertcat(ADD{:});
            ADD=ADD(:);             
            [AA,BB]=equationsToMatrix([EQU; ADD],SYMS); %make symbolic equations to matrix , vector system

            if(~symflag); AA=double(AA); BB=double(BB); end;  %solve the system
            if(pinvflag); X=pinv(AA)*BB;
            else;         X=AA\BB;
            end;
            
            X=double(X);
            
            ERR=norm(double(AA)*X-double(BB),inf); %check if the solution is really a solution
            if(ERR<1e-13 && ~isnan(ERR)); break; end;
            
        end
        if(ERR>1e-13 || isnan(ERR)); vprintf('differencescheme: Computational error = %i for difference scheme of order %i for subdivision scheme %i, name: %s\n', ERR, k, j,S{j,6},'imp',[1,verbose]); end;
        warning('on','MATLAB:rankDeficientMatrix');
        warning('on','symbolic:mldivide:RankDeficientSystem');
        warning('on','symbolic:mldivide:InconsistentSystem');

        for i=1:numel(B); %replace symbolic variables with their values
            B.M{i}.c=double(subs(B.M{i}.c,SYMS,X.'));
            B.M{i}=B.M{i}.shrink;
        end;
        Sp{j,1}=B; 
        Sp{j,6}=strcat(Sp{j,NAMEIDX},[' (diff ' num2str(k) ')']);
    end

    Sp=getS(Sp,'noprint');



end


function [Sp,mu] = differencescheme_linearsystem(varargin)

if(size(varargin,2)<=1); error('differencescheme: Too less arguments. Order of differences ''k'' is missing.'); end;
noprint=parsem('noprint',varargin);
symflag=parsem('sym',varargin);
slashflag=parsem('slash',varargin);
S=varargin{1}; S=getS(S);
k=varargin{2}; %order of difference-scheme
J=size(S,1);
if(J>1 && isempty(k) || J>1 && length(k)>1); error('differencescheme: ''k'' must be an integer.'); end;
if(any(k<0)); error('differencescheme: ''k'' must be nonnegative.'); end;
dim=size(S{1,2},1);


Sp=S; %initialize Sp

%make vector of differences
mu=constructmu(k,dim); if(dim==1); mu=[mu 0]; end;
sizemu=size(mu,2); if(dim==1); sizemu=1; end;

    for j=1:J 

        Sj=S(j,:);
        if(~isequal(size(Sj{1}),[1 1])); error('differencescheme: Works only for scalar-valued schemes.'); end;
        Mj=Sj{2};
        Sp{j,MASKIDX}=cell(sizemu); %mask
        %{j,MIDX}==M %stays the same
        %{j,DIDX}==D %stays the same
        Sp{j,4}=[]; %Omega
        Sp{j,MASKMINIDX}=unflatten(Sj{MASKMINIDX}{1},cell(sizemu),'assign'); %amin
        Sp{j,NAMEIDX}=strcat(Sp{j,NAMEIDX},[' (diff ' num2str(k) ')']);
        Sp{j,MASKMAXIDX}=[]; %amax
        Sp{j,SUPPMASKIDX}=[]; %suppa
        Sp{j,FLAGSIDX}=int32(0); %flags
        %{j,10}==det %stays the same
        Sp{j,SYMIDX}=[];

        k=1; %order
        
        aj=S{j,1}{1};

        SUPP=supp(aj,2,[1 1]');



        %get biggest d
        BIGD=zeros(dim,sizemu); if(dim==1); BIGD=zeros(2,sizemu); end;
        for l=1:sizemu
            d=upsamplen(diffsequ(1,mu(:,l)),Mj);
            BIGD(:,l)=size(d)';
        end
        BIGD=max(BIGD,[],2).';
        BIGD=BIGD-1;

        SIZE=size(aj)+BIGD; %sizes of the output mask

        A=zeros(SIZE); %I need that matrix


        lhs=cell(sizemu,1);
        for i=1:sizemu        
            lhs{i}=zeros(numel(aj),numel(A));
            d=diffsequ(1,mu(:,i));
            shift=supp(d,2,[0 0]');
            for kk=1:size(SUPP,2)
                s=setplus(SUPP(:,kk),shift,'nounique');
                N=subsco(A,s,d,'save');
                N=N(:).';
                lhs{i}(kk,:)=N;  
            end
            lhs{i}=lhs{i}.';
        end

        RHS=cell(sizemu,sizemu);
        origin=cell(sizemu,sizemu);
        rhs=cell(sizemu,1);
        for i=1:sizemu
            for l=1:sizemu
                d=upsamplen(diffsequ(1,mu(:,l)),Mj);

                shift=supp(ones(size(d)),2,[0 0]');

                RHS{i,l}=zeros(numel(aj),numel(A));
                for kk=1:size(SUPP,2)
                    s=setplus(SUPP(:,kk),shift,'nounique');
                    N=subsco(A,s,d,'save');
                    N=N(:).';
                    RHS{i,l}(kk,:)=N;
                end
                RHS{i,l}=RHS{i,l}.';
            end
            rhs{i}=horzcat(RHS{i,:});
        end

        B=cell(sizemu,1); %\-operator
        for i=1:sizemu
            
            
            if(symflag);    rhs{i}=sym(rhs{i}); end;
            
            if(slashflag);  
                warning('off','MATLAB:rankDeficientMatrix');
                B{i}=rhs{i}\lhs{i}*aj(:); 
                warning('on','MATLAB:rankDeficientMatrix');
            else; 
                B{i}=pinv(rhs{i})*lhs{i}*aj(:); 
            end;
            
            
            %Test result
            ERR=double(norm(rhs{i}*B{i}-lhs{i}*aj(:),1));
            if(ERR>1e-14);
                vprintf('err','differencescheme: Computational error = %i for difference scheme of order %i for subdivision scheme %i, name: %s\n', ERR, k, j,S{j,6},'imp',[1,verbose]); 
            end
            for l=1:sizemu
                RHS{i,l}=B{i}(numel(aj)*(l-1)+1:numel(aj)*l);
                RHS{i,l}=double(RHS{i,l});
                RHS{i,l}=reshape(RHS{i,l},size(aj));
                [RHS{i,l},origin{i,l}]=removezero(RHS{i,l},'border'); %save the change of amin
                origin{i,l}=origin{i,l}-1;
            end
        end
        Sp{j,1}=RHS;
        Sp{j,MASKMINIDX}=nestedcellfun(@plus,Sp{j,MASKMINIDX},origin,'UniformOutput',false);
        

        
    end
    
    Sp=getS(Sp,'verbose',verbose);

end



function [Sp,mu]=differencescheme_polynomials(varargin)
%produces wrong output
cprintf('err','differencescheme ''polynomials'': This algorithm yields wrong results and should not be used.\n');

if(size(varargin,2)<=1); error('differencescheme: Too less arguments. Order of differences ''k'' is missing.'); end;
verbose=parsem({'verbose','v'},varargin); %sic
S=varargin{1}; S=getS(S);
k=varargin{2};
J=size(S,1);
if(J>1 && isempty(k) || J>1 && length(k)>1); error('differencescheme: ''k'' must be an integer.'); end;
if(any(k<0)); error('differencescheme: ''k'' must be nonnegative.'); end;
dim=size(S{1,2},1);
z=sym(zeros(1,0));

for j=1:J; z=union(z,symvar([S{j,11}{:}])); end; %get all variables from the symbols 
if(isrow(z)); z=z.'; end;
%test if number of variables is equals to the dimension
if(size(z,1)~=dim); 
    vprintf('err','Found Variables: \n%v\nDimension: %i\n',z,dim,'cpr','err','imp',[1,verbose]); 
    error('differencescheme: Number of variables in symbols unequal to dimensions.'); 
end

Sp=S; %initialize Sp



%make vector of differences
mu=constructmu(k,dim);
sizemu=size(mu,2);


for j=1:J
    Sj=S(j,:);
    if(~isequal(size(Sj{1}),[1 1])); error('differencescheme: Works only for scalar-valued schemes.'); end;
    Mj=Sj{2};
    Sp{j,1}=cell(sizemu); %mask
    %{j,2}==M %stays the same
    %{j,3}==D %stays the same
    Sp{j,4}=[]; %Omega
    Sp{j,5}=cell(sizemu); %amin
    Sp{j,NAMEIDX}=strcat(Sp{j,NAMEIDX},[' (diff ' num2str(k) ')']);
    Sp{j,7}=[]; %amax
    Sp{j,8}=[]; %suppa
    Sp{j,9}=int32(0); %flags
    %{j,10}==det %stays the same
    Sp{j,11}=cell(sizemu); %symbol
    sya=Sj{11}{1}; % I can put here {1}, because Sj is a scalar valued scheme.

    syright=sym(zeros(sizemu,1));   %syright = (1-z^M)^\mu
    syleft=sym(zeros(sizemu,1));    %syleft  = (1-z)^\mu
    zM=sym(zeros(dim,1));           %z^M

    for i=1:dim; 
        zM(i)=prod(z.^Mj(:,i)); 
    end;
    for i=1:sizemu
        syleft(i)=prod((1-z).^mu(:,i));
        syright(i)=prod((1-zM).^mu(:,i));
        
    end

    [~,corrector2]=numden(sya);
    for m=1:sizemu %rows of difference mask B
        rem=syleft(m).*sya*corrector2;
        for n=1:sizemu %columns of difference mask B
            %[quo2,rem2]=quorem(rem,syright(n),z2);
            
            %Use MuPAD because Matlab cant do it %XX I dont want to use MuPAD here
            %Compute corrector for 'divide'
            
            [~,corrector1]=numden(syright(n));
            val=feval(symengine,'divide',expand(rem*corrector1),expand(syright(n)*corrector1)); 
            if(size(val,2)==1); rem=sym(1); break; end; %sym(1) indicates that the computation failed
            rem=val(2)/corrector1; quo=val(1);
            quo=expand(quo);
            rem=expand(rem);
            [b,bmin]=symbol2mask(quo/corrector2,'dim',dim,'var',z);
            Sp{j,1}{m,n}=b; %a
            Sp{j,5}{m,n}=bmin; %bmin
            Sp{j,11}{m,n}=quo; %symbol
        end
        if(~isequal(rem,sym(0))); %then there we could not find a difference scheme of order k
            vprintf('differencescheme: Failed to find difference scheme of order %i for subdivision scheme %i, name: %s\n', k, j,S{j,6},'imp',[1,verbose]); 
            Sp{j,1}=cell(sizemu);
            Sp{j,4}=[];
            break;
        end;

    end

end

Sp=getS(Sp,'verbose',verbose);


end




function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 


