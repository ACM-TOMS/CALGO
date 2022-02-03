function [probs,funs,bcfuns,dfuns,dbcfuns,esolus,setoutputs,settolerancess]=bvpT26(ExtraArgs)
%   For a description of the input parameters see the report manual_testset
%
%
%   
%
%       Authors:
%
%       Jeff R. Cash 
%            (Department of Mathematics, Imperial College,  London, England.)
%       Davy  Hollevoet 
%            (Vakgroep Toegepaste Wiskunde en Informatica, Universiteit Gent, Belgium.)
%       Francesca Mazzia  
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%       Abdelhameed Nagy Abdo
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%            (Dept. of Mathematics, Faculty of Sciences, Benha  University,Egypt)
%            
%
probs    = @prob;
funs     = @fun;
bcfuns   = @bcfun;
dfuns    = @dfun;
dbcfuns  = @dbcfun;
esolus   = @esolu;
setoutputs = @setoutput;
settolerancess = @settolerances; 
if nargin==1 && ~isempty(ExtraArgs)
    if iscell(ExtraArgs)
       lambda=ExtraArgs{:};
    else
        lambda=ExtraArgs;
    end
end

    function [problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = prob()
        problm = 'bvpT26';
        type   = 'ODEBVP';
        t(1)   = 0;
        t(2)   = 1;
        m   = 2;       
        y0  = 0;        
        Linear  = 'off';
        numjac = 0;
        numbcjac = 0;
        Vectorized = 'on';
        JVectorized = 'on';
        solinit = bvpinit(linspace(t(1),t(2),11),repmat(y0,1,m));
    end

    function F = fun(X,Z,ExtraArgs)
        if nargin==3
           lambda=ExtraArgs;
        end
        F=zeros(size(Z));
        F(1,:) = Z(2,:);
        F(2,:) = (Z(1,:).*(1.e0-Z(2,:)))/lambda;
    end

    function  bc = bcfun(ya,yb,ExtraArgs)
        if nargin==3
           lambda=ExtraArgs;
        end
        C0 = [1,0;0,0];
        C1 = [0,0;1,0];
        Eta = [1;-1/3];
        bc = C0*ya + C1*yb  - Eta;
    end

    function Df = dfun(X,Z,ExtraArgs)
        if nargin==3
           lambda=ExtraArgs;
        end
        ncomp=size(Z,1);
        nmsh=size(Z,2);
        Z=reshape(Z,ncomp,1,nmsh);
        Df=zeros(ncomp,ncomp,nmsh);
        Df(1,2,:) = 1.0e0;
        Df(2,1,:) = (1.e0-Z(2,1,:))/lambda;
        Df(2,2,:) = -Z(1,1,:)/lambda;
    end

    function  [C0,C1] = dbcfun(ya,yb,ExtraArgs)
        if nargin==3
           lambda=ExtraArgs;
        end
        C0 = [1,0;0,0];
        C1 = [0,0;1,0];
    end

    function  tolvec = settolerances(tol)
     tolvec = tol;
    end 

    function [solref,printsolout,nindsol,indsol] = setoutput(plotsol)
      solref = 0; 
       if isempty(plotsol)
          printsolout = 0;
          nindsol = 1;
          indsol = 1;
       else
          printsolout = 1;
          nindsol = length(plotsol);
          indsol = plotsol;  
       end         
    end

    function Exact = esolu(X,ExtraArgs) 
        error('MATLAB:twpbvpc:bvp_examples:noexactsolution', 'No exact solution available');
    end
end