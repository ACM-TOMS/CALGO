function varargout = constructOmega(varargin)
% [ Omega ] = constructOmega( S, [options] )
% Constructs the set Omega as described in Charina, Mejstrik, 2018
%
% Input: 
%       S       cell array of subdivision operators
%
% Options:
%   'Omega',val          Starting set. If not given, a point is computed where to start from.
%                        In some cases, the algorithm may not find a good starting point, and the returned set is not minimal
%   'plot'               plots the set Omega
%   'verbose',val        verbose level
%   'lexicographic'     (experimental) orders the set Omega lexicographically
%   'stable'             does not order the set Omega
%   'legacy'            (experimental) Uses a different algorithm which is slower in Matlab, but which constructs the same set
%   'V',val             (experimental) constructs set Om such that Xmu is non-empty for all abs(mu)<=val+1, see [Mejstrik, PhD Thesis, 2019]
%                        This option is needed if one wants to compute whether a subdivision scheme is C^val or not 
%                        Even without this option, all sets X_mu may be non-empty
%
% Output:
%       Omega       the set Omega
%
% E.g.: constructOmega('2_butterfly','Omega',[2;2],'stable')
%
% See also: getS, transitionmatrix
%
% Reference:
%   Algorithm is described in:
%       Maria Charina, Thomas Mejstrik,
%       Multiple multivariate subdivision schemes: Matrix and operator approaches, 
%       Journal of Computational and Applied Mathematics, 2018
%
% Written by: tommsch, 2017

% XX Option which constructs Omega big enough V=Vt, Xmu non-empty
% XX Option which constructs smallest possible Omega

%#ok<*ALIGN>

[lexico, varargin] = parsem({'lexicographic','lexico'},varargin);
[stable_flag, varargin] = parsem('stable',varargin);
[plotflag, varargin] = parsem('plot',varargin);


[Om]=constructOmega_worker(varargin{:});    
dim=size(Om,1); 

if(lexico); Om=lexicographic(Om,varargin{lexico:end}); 
elseif(~stable_flag); Om=unique(Om','rows')'; end;
varargout{1}=Om;

if(plotflag && dim<=3); 
    plotm(Om,'.'); axis equal; drawnow; 
end;

end

function [Om] = constructOmega_worker(varargin)
% Constructs the set Omega for the multiple subdivision scheme case.
% Iteratively adds points to Omega, till all transition matrices are invariant under the set Omega

S=getS(varargin{1}); varargin(1)=[];
if(~isS(S)); 
    error('First input argument must be a cell array of subdivision schemes.'); end;
[verbose,varargin] = parsem({'verbose','v'},varargin,1);
[Om,varargin] = parsem({'Omega','Om'},varargin,[]);
[legacy,varargin] = parsem('legacy',varargin);
[Vval,varargin] = parsem('V',varargin,[]);

if(~all(iswholenumber(Om))); 
    error('<''Omega'',val> must be integers.'); end;
parsem(varargin,'test');

sizeS=size(S,1);
dim=size(S{1,2},1);

if(isempty(Om)); %compute a point which is in the set K_{supp a}-K_D
    Om=zeros(dim,1);
    for i=1:10 
        SUPP=S{1,1}.supp;
        Om=S{1,2}^(-1)*(SUPP(:,1)-S{1,3}(:,1)+Om); end;
    Om=round(Om);
    if(~isempty(Vval));
        Om=setplus(Om,mixvector(0:Vval+1,dim)); end; %XX very bad computation for non-empty Xmu
    vprintf('You should specify a point which is for sure inside of Omega via <''Omega'',vector>. My guess is: \n%v\n',Om,'imp',[2,verbose]); 
end;
if(size(Om,1)~=dim); 
    error('Dimension of set Omega is wrong.'); end;

vprintf('\nConstruct Omega: \n','imp',[2,verbose]); 
sizebefore=0;
while(true)
    vprintf('Omega: \n%v\n',Om(:, sizebefore+1:end),'imp',[2,verbose]); 
    sizebefore=size(Om,2);
    for j=1:sizeS %iterate through all S
        Mj=S{j,2}; Dj=S{j,3}; detMj=abs(det(Mj));
        suppaj=S{j,1}.supp;
        if(~legacy)
            Omnew=Mj\setplus(suppaj,Om,-Dj);
            Omnew=round(Omnew(:,sum(abs(Omnew-round(Omnew)),1)<1/(2*detMj)));
            Om=unique([Om Omnew]','rows','stable')';
        else
           %this variant seems slower
            for di=1:size(Dj,2); %iterate through all of D
                dj=Dj(:,di);
                Omnew=Mj\setplus(suppaj,Om,-dj);
                Omnew=round(Omnew(:,sum(abs(Omnew-round(Omnew)),1)<1/(2*detMj)));
                Om=unique([Om Omnew]','rows','stable')';
            end
        end
    end
    numadded=size(Om,2)-sizebefore;
    vprintf('Number of added points: %i \n',numadded,'imp',[2,verbose]); 
    if(verbose>=3); 
        plotm(Om,'x','resolution',0); 
        drawnow; 
        if(verbose>=4); 
            pause; end
    end;
    if(numadded==0); 
        break; end;
end
vprintf('\n','imp',[2,verbose]); 


end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 