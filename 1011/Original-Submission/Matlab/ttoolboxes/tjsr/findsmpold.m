function [ cand, nearlycand, info] = findsmp(varargin)
% [ cand, nearlycand, info ] = findsmp( T, [options] )
% Searches for s.m.p.-candidates. Various algorithms available.
%
% Input:    
%   T                       Cell array of square matrices of the same size.
%
% Output: 
%   cand                    Ordering of the products with highest spectral radius.
%   nearlycand              Orderings of products with spectral radius greater than val*bounds(1). Can be empty, depending on the algorithm.
%   info                    (struct) Additional info, depends on the algorithm used.
%                               info.time         time in seconds needed for computation
%                               info.jsrbound     (Interval) estimate on the JSR based on the spectral radius of the found candidates and norm estimates
%                               info.graph        All computed norms and spectral radii as directed graph
%                               info.spectralgap  relative difference between info.jsrbound and second biggest eigenvalue found (either from candidate or from other candidates)
%                               info.count        approximate number of computed matrices
%   
% Algorithms implemented: 'modifiedgripenberg' (default), 'gripenberg', 'bruteforce', 'genetic' (not recommended)
%
% Options for all algorithms:
%   'maxsmpdepth', val              maximal length of products which is searched for. Can be 
%                                   very high (>100) for 'gripenbergmodifed' algorithm (default)
%                                   high (<15) for 'gripenberg'
%                                   small (<12) for 'brute force'
%   'verbose',val                   Defines the verbose level. Default=1
%   'maxtime',val                   Maximal runtime in seconds (is only checked at the beginning of each iteration).
%   'norm',functionhandle           Handle to a norm function. Default: 2-Norm.
%                                      'bruteforce' and 'genetic' do not compute norms.
%   'nosimplify'                    Does not simplify products of c and nc
%   'maxtime'                       Maximal time used for computation
%   'minJSR',val                    Candidates to be found must have at least val spectral radius.
%   'searchonlyonecandidate'        Searches until one candidate is found, maxsmpdepth is ignored
%   'maxeval'                       Maximum number of evaluations (approximate)
%   'testsuite'                     Runs a self-test
%
% Options which apply to only one algorithm: 
% ------------------------------------------
% 'modifiedgripenberg' (default)    Keeps at most 2*N products in each step. Therefore can compute very long products, 
%                                   Gives only good lower bounds.
%                                   Can miss candidates!
%       'N',val                         Number of kept products in each step. If N==inf, then the algorithm behaves like default Gripenberg algorithm.
%       'searchonlyonecandidate'        Searches until one candidate is found, maxsmpdepth is ignored
%       'sufficientbound',val           Algorithm terminates if a candidate with bound larger than val has been found
%       'plot',string                   Some plot output
%                                           'none' or Default       no plot output
%                                           'tree'                  Return the field info.graph and plots the graph
%       'nearlycandidate',val           Nearlycandidates must have spectral radius larger than val*bounds(1), default: val=0.99
%       'shortnearlycandidate',val      Removes all nearlycandidates whose ordering is longer than the val*maximal-length-of-candidates-ordering. Default: true
%       'sparse',val                    Flag if the matrices are sparse. If not given, it is determined automatically.
%       'select',val                    How to select new candidates:
%                                       0: N with biggest norm (bad)
%                                       1: N/2 with biggest and N/2 with smallest norm (default)
%                                       2: N random (not so good)
%
% 'genetic'                         Written by Chia-Tche Chang, 2011. 
%                                   There is a bug and the returned value is sometimes wrongly normalized.
%       'popsize',val                   (default 300, minimum 10) population size.
%       'maxgen',val                    (default 1000) maximum number of generations.
%       'maxstall',val                  (default 15) maximum number of stalling iterations before increasing the maximum product length.
%       'maxtotstall',val               (default 100) maximum number of stalling iterations before terminating the algorithm.
%                                   
%       'mutantprop',val                (default 0.3) mutation probability of a given product.
%       'muteprop',val                  (default 0.2) mutation proportion of a given product.
%
% 'gripenberg'                      Gripenbergs algorithm as implemented by R. Jungers (jsr_prod_gripenberg), modified to have the same frontend as the other functions
%                                   Gives lower and upper bounds.
%                                   Returns only one candidate. Not suitable for finding smp-candidates.
%                                   Gives no nearlycandidates
%       'delta'                         Delta used in the Gripenberg Algorithm. Default: .99
%
% 'bruteforce'                      Tests all possible orderings up to depth 'maxsmpdepth'.
%                                   Gives only good lower bounds (to speep up computation)
%                                   Does not miss candidates
%       'minsmpdepth',val               Minimal length of products
%
% E.g.: [ c, nc, info]=findsmp({[1 -1; 3 -2], [1 3; -1 -1]},'maxsmpdepth',15)
%
% See also: tjsr
%
% Written by: tommsch, 2018

% XX implement spectralgap
% XX implement LSR
% XX Fuer modgrip selbst herausfinden lassen, wann vpa verwendet werden soll
% XX Legende fuer modgrip hinzufuegen
% XX Speichere auch matrizen, sodass nicht immer das gesamte Produkt berechnet werden muss
% XX Implement 'hardworking'
% XX Alle Algorithmen sollten die gleichen Optionen haben, insb. searchonlyonecandidate, maxtime, nearlycandidate, shortnearylcand, JSR, sufficientbound, nosimplify

if(parsem('selftest',varargin)); selftest(); return; end;

verbose=parsem({'verbose','v'},varargin,1);
starttime=clock;

vprintf('Search candidate-smp: ','imp',[1,verbose]);
vprintf('\n','imp',[2,verbose]);

[geneticflag,varargin]=             parsem({'genetic','gen'},varargin);
[bruteforceflag,varargin]=          parsem({'bruteforce','bf'},varargin);
[gripenbergflag,varargin]=          parsem({'gripenberg','grip'},varargin);
[modifiedgripenberg,varargin]=      parsem({'modifiedgripenberg','modgrip'},varargin); %#ok<ASGLU>
[nosimplify,varargin]=              parsem('nosimplify',varargin);

if(geneticflag);
    [cand,nearlycand,info]= tjsr_genetic(varargin{:});
elseif(bruteforceflag);
    [cand, nearlycand, info]= tjsr_bruteForce(varargin{:}); 
elseif(gripenbergflag);
    [cand, nearlycand, info]= tjsr_gripenberg(varargin{:}); 
else %modifiedgripenberg
    [cand, nearlycand, info]=tjsr_gripenberg_modified(varargin{:}); end;

if(~nosimplify)
    %Simplify candidates and nearlycandidates
    vprintf('\n', 'imp',[2,verbose]);             cand=simplify_ordering(cand);
    vprintf('Simplify %i candidates. ', numel(cand), 'imp',[1,verbose]);             cand=simplify_ordering(cand);
    vprintf('Simplify %i nearly candidates. ', numel(nearlycand), 'imp',[1,verbose]);      nearlycand=simplify_ordering(nearlycand);
    vprintf('Remove wrong nearly candidates. ','imp',[1,verbose]);  nearlycand=remove_ordering(nearlycand,cand);
    vprintf('\n','imp',[1,verbose]); end;

vprintf('Length of candidates / nearly-candidates: %r / %r. ', cellfun(@length,cand), cellfun(@length,nearlycand), 'imp',[1,verbose]);
vprintf('\n', 'imp',[2,verbose]);             cand=simplify_ordering(cand);
vprintf('Candidates: \n%v\n', cand, 'imp',[2,verbose]);
if(~isempty(nearlycand));
    vprintf('Nearly candidates: \n%v\n', nearlycand, 'imp',[2,verbose]); end;
vprintf('\n', 'imp',[2,verbose]);             cand=simplify_ordering(cand);
vprintf('Bounds on the jsr : [%.15g, %.15g]\n', info.jsrbound(1), info.jsrbound(2),'imp',[1,verbose]);

info.time=etime(clock,starttime);

end


function [cand, nearlycand, info] = tjsr_gripenberg_modified(varargin)
%this function tries to find candidates in a fast way, but it could miss some potential candidates

M=varargin{1}; varargin(1)=[];
J = length(M);
dim=size(M{1},1);

[normfun,varargin]=                parsem('norm',varargin,[]); 
if(isempty(normfun)); 
    normfun=@(x) norm(x,2); end; %which norm to use
if(isnumeric(normfun)); 
    normfun=@(x) norm(x,normfun); end;
[maxsmpdepth,varargin]=            parsem('maxsmpdepth',varargin,100); %maximum length of products to be computed
    val=ceil(min( sqrt(dim*J)*10, tavailable_memory()*3/4*1/(2*J*(maxsmpdepth+2)*8*feature('numcores')) ));
[N,varargin]=                      parsem('N',varargin,val); %the number of kept products in each step equals dim*J
[minJSR,varargin]=                 parsem('minJSR',varargin,0);    %Spectral radius: current lower bound for JSR
[searchonlyonecandidate,varargin]= parsem('searchonlyonecandidate',varargin); %stops algorithm after one candidate has been found. To be used together with minJSR
[verbose,varargin]=                parsem({'verbose','v'},varargin,1);
[sufficientbound,varargin]=        parsem('sufficientbound',varargin,inf); %stops if (J)SR>val resp. (L)SR<val
[plotflag,varargin]=               parsem('plot',varargin,'none'); %flag that the graph shall be returned
if(isequal(plotflag,'tree')); 
    treeflag=true; 
else; 
    treeflag=false; end;
[nearlycanddelta,varargin]=        parsem('nearlycandidate',varargin,.99); %delta for nearlycandidates

[delta,varargin]=                  parsem('delta',varargin,1); %delta from Gripenberg Algorithm
[vpaflag,varargin]=                parsem('vpa',varargin,0); %uses vpa to compute spectral radii. val=0: do not use vpa (default), val=1: use vpa if difference to current bound is less than 100*eps, val=2: always use vpa.
[epsilon,varargin]=                parsem('epsilon',varargin,10e-11); 
[shortnearlycand,varargin]=        parsem('shortnearlycandidate',varargin,1); %only keeps nearlycandidates whose length is less or equal than val times the length of the shortest candidate
[maxnumnearlycand,varargin]=       parsem('maxnumnearlycandidate',varargin,100); %maximal number of nearlycandidates. If there are more, then delta is set to (delta+1)/2
[selecttype,varargin]=             parsem('select',varargin,1); %How to select new candidates
[maxtime,varargin]=                parsem('maxtime',varargin,inf);
[maxeval,varargin]=                parsem('maxeval',varargin,inf);
[sparse,varargin]=                 parsem('sparse',varargin,-1);

parsem(varargin,'test');

oo=1:J;                     %the orderings of the products to be checked
G=cell(2,0);                %the graph, is transformed to a graph at the end of the function
cand=cell(1,0);             %list of candidates
nearlycand=cell(1,0);       %list of nearly candidates
nearlycandrho=zeros(1,0);   %list of nearly candidates spectral radii
starttime=clock;

opts.disp = 0;
if(sparse==-1)    
    val=cellfun(@(x) nnz(x)/numel(x),M);
    if(all(val<.1)); 
        sparse=1; 
    else; 
        sparse=0; end; end;

if(vpaflag)
    olddigits=digits; %save current value for digits
    Msym=cell(1,J);
    for j=1:J
        Msym{j}=sym(M{j});
    end
else
    for j=1:J
        M{j}=double(M{j}); end; end;

vprintf('Search for candidates: ','imp',[2,verbose]);
counter=0;
while(true)
    vprintf('Time: %s\n',datetime('now'),'imp',[4,verbose]);
    counter=counter+size(oo,2);
    vprintf('.','imp',[1,verbose]);
    sze=size(oo,2);
    NRj=zeros(2,sze);
    loo=size(oo,1);
    vprintf('\bTest %i matrices.\n',sze,'imp',[3,verbose]);
    switch vpaflag
        case 0;
            if(sparse)
                for i=1:sze; 
                    P=tbuildproduct_fast(M,oo(:,i));
                    NRj(:,i)=[normfun(P); eigs(P,1,'LM',opts)]; end;
            else
                parfor i=1:sze; 
                    P=tbuildproduct_fast(M,oo(:,i));
                    NRj(:,i)=[normfun(P); trho(P)]; end; end;
            
            NRj=NRj.^(1/loo); %Norm und Spektralradius der Kandidaten
        case 1;
            digits(dim*loo+32);
            epsilon=vpa(10^(-digits+10));
            if(sparse)
                parfor i=1:sze; 
                    P=tbuildproduct_fast(M,oo(:,i));
                    NRj(:,i)=[normfun(P); eigs(P,1,'LM',opts)]; end;
            else
                parfor i=1:sze; 
                    P=tbuildproduct_fast(M,oo(:,i));
                    NRj(:,i)=[normfun(P); trho(P)]; end; end;
            
            NRj=NRj.^(1/loo); %Norm und Spektralradius der Kandidaten            
            val=NRj(1,:)/max(minJSR,max(NRj(1,:)))-1;
            idx=find(abs(val)<100*eps);
            NRj2=sym(zeros(1,size(idx,2)));
            vprintf('(%i)',size(idx,2),'imp',[1,verbose]);
            parfor i=1:size(idx,2);
                P=tbuildproduct(Msym,oo(:,i));
                NRj2(i)=vpa(trho(P).^(1/loo));
                if(verbose>=1); 
                    fprintf('.'); end; end;
            NRj=vpa(NRj);
            NRj(2,idx)=NRj2;
            
        case 2;
            NRj=sym(NRj);
            digits(dim*loo+32);
            epsilon=vpa(10^(-digits+10));
            if(sparse)
                parfor i=1:sze; 
                    P=tbuildproduct_fast(Msym,oo(:,i));
                    NRj(:,i)=[vpa(normfun(P)); vpa(eigs(P,1,'LM',opts))]; end;
            else
                parfor i=1:sze; 
                    P=tbuildproduct_fast(Msym,oo(:,i));
                    NRj(:,i)=[vpa(normfun(P)); vpa(trho(P))]; end; end;

            NRj=NRj.^(1/loo); %Norm und Spektralradius der Kandidaten            
            
        otherwise;
            error('bad value for ''vpa''.'); end;

    RMAX=max(NRj(2,:)); %biggest spectral radius in this level
    %clear cand if higher spectral radius was found
    %add cand to nearlycand
    if(minJSR<RMAX*(1-epsilon)); 
        nearlycand=[nearlycand cand]; %#ok<AGROW>
        nearlycandrho=[nearlycandrho repmat(minJSR,1,size(cand,2))]; %#ok<AGROW>
        cand={};
        minJSR=RMAX;                %set new value for JSR
        vprintf('Time: %s\n',datetime('now'),'imp',[4,verbose]);
        vprintf('Candidate found: ','imp',[2,verbose]);
        if(vpaflag && verbose>=1)
            RMAX %#ok<NOPRT>
        elseif(verbose>=2);
            fprintf('rho= %15.12g ',double(RMAX)); end;
        vprintf('\n','imp',[2,verbose]); end;
    
    if(treeflag); %store computed values in the graph: Format: First row: N, Rho, Second Row: Ordering
        G(1:2,end+1:end+size(oo,2))=[num2cell(double(NRj(1:2,:)),1); num2cell(oo,1)]; end; 
    
    %add candidates to cand
    idxrho=find(NRj(2,:)>=minJSR*(1-epsilon)); %indices of candidates
    cand=[cand num2cell(oo(:,idxrho),1)]; %#ok<AGROW>
    if(~isempty(idxrho) && verbose >=2); 
        for iii=1:length(idxrho);
            fprintf('%i ',oo(:,idxrho(iii)));
            fprintf('--- '); end
        fprintf('\n'); end;
    
    %add nearlycandidates to nearlycand (and nearlycandrho)
    idxrho=intersect(find(NRj(2,:)<minJSR*(1-epsilon)), find(NRj(2,:)>minJSR*nearlycanddelta)); %indices of candidates
    nearlycand=[nearlycand num2cell(oo(:,idxrho),1)]; %#ok<AGROW>
    nearlycandrho=[nearlycandrho NRj(2,idxrho)]; %#ok<AGROW>
    
    if(isempty(cand)); 
        idxlength=[]; 
    else; 
        idxlength=cellfun(@length,nearlycand)>min(cellfun(@length,cand))*shortnearlycand; end;
        
    nearlycand(idxlength)=[];
    nearlycandrho(idxlength)=[]; 
    if(size(nearlycand,2)>maxnumnearlycand); 
        nearlycanddelta=min(max(nearlycandrho)/minJSR-100*eps,(nearlycanddelta+1)/2);
        if(verbose>=2); 
            vprintf('Too many nearlycandidates (# %i). Set delta to %g. \n',size(nearlycand,2),nearlycanddelta,'imp',[2,verbose]); 
        elseif(verbose>1); 
            vprintf('delta = %i ',nearlycanddelta,'imp',[2,verbose]); end; end; 
    idx=nearlycandrho<minJSR*nearlycanddelta; 
    idx=idx | nearlycandrho>minJSR*(1-epsilon);
    
    nearlycand(idx)=[]; 
    nearlycandrho(idx)=[]; %clean up nearlycand
    
    
    idx=NRj(1,:)<minJSR*(1-epsilon)/delta;
    NRj(:,idx)=[]; %remove all orderings which have norm less than SR
    oo(:,idx)=[];
    
    
    %Select new candidates
    switch selecttype
        case 0;
            [~,idx]=sortrows(NRj'); %sort correspdonding to Norm
            idx=idx.';
            noo=size(oo,2);
            if(noo>N); %keep highest and lowest norms %funktioniert bei "findsmp(tgallery('CH'),'maxsmpdepth',50,'N',2,'vpa','graph')" besser
                oo=oo(:,[idx(noo-ceil(N)+1:noo)]); end;  %#ok<NBRAK>
            noo=size(oo,2);
            oo=[repmat(oo,[1 J]); reshape(repmat(1:J,[noo 1]),1,[])]; %all possible new orderings of products
        case 1;
            [~,idx]=sortrows(NRj'); %sort correspdonding to Norm %XX ist das das richtige?
            idx=idx.';
            noo=size(oo,2);
            if(noo>2*N); %keep highest and lowest norms
                oo=oo(:,[idx(1:ceil(N/2)) idx(noo-ceil(N/2)+1:noo)]); end; 
            noo=size(oo,2);
            oo=[repmat(oo,[1 J]); reshape(repmat(1:J,[noo 1]),1,[])]; %all possible new orderings of products
            
        case 2
            noo=size(oo,2);
            if(noo>N); 
                oo=oo(:,randperm(noo,N)); end;
            noo=size(oo,2);
            oo=[repmat(oo,[1 J]); reshape(repmat(1:J,[noo 1]),1,[])]; %all possible new orderings of products
            
        otherwise
            error('Wrong value for ''select''.'); end;

    
    if(size(oo,1)>maxsmpdepth && ~searchonlyonecandidate); 
        break; end;   %test for maxsmpdepth
    if(searchonlyonecandidate && ~isempty(cand)); 
        break; end;           %test for searchonlyonecandidate
    if(minJSR>sufficientbound); 
        break; end;                                 %test for sufficientbound
    if(etime(clock,starttime)>=maxtime); 
        break; end;
    if(counter>maxeval); 
        break; end;
end
vprintf('Computed matrices in total: %i\n','cpr',[.1,0.5,.1],counter,'imp',[2,verbose]);

info=struct;
info.count=counter;
if(N~=inf)
    info.jsrbound=[minJSR inf];
else
    info.jsrbound=[minJSR max(NRj(1,:))]; end;


%make G to graph
if(treeflag)
    val=cell2mat(G(1,:));
    G=makeorderinggraph(G(2,:), 'value',{'rho',val(2,:).';'norm',val(1,:).'},'verbose',verbose);
    title('Computed orderings');
    %G=makeorderinggraph(G(2,:), 'norm',normval, 'rho',rhoval, 'JSR',minJSR, 'delta',delta, 'movevertex');
    %G=makeorderinggraph(G(2,:), 'norm',normval, 'movevertex');
    info.G=G; end;

if(vpaflag)
    digits(olddigits); end;

end

function [cand, nearlycand, info] = tjsr_gripenberg(varargin)
% JSR_PROD_GRIPENBERG(M) Computes bounds by a branch and bound method on
%                        products of growing length.
%
% [BOUNDS, PRODOPT, ALLPROD, INFO] = JSR_PROD_GRIPENBERG(M)
%     For M a cell array of matrices, computes products of growing length
%     and keeps only the products that could produce an upper bound higher
%     than current_lower + delta. See default values below.
%cand
% [BOUNDS, PRODOPT, ALLPROD, INFO] = JSR_PROD_GRIPENBERG(M,OPTS) or
% [BOUNDS, PRODOPT, ALLPROD, INFO] = JSR_PROD_GRIPENBERG(M,OPTSNAME,OPTSVALUE)
%     Uses the values of the parameters in the structure OPTS generated by
%     jsrsettings. Or the values given by pair of OPTSNAME, OPTSVALUE
%     arguments. OPTSNAME can be 'delta', 'maxEval' or 'normfun' for
%     instance.
%
%   BOUNDS  is a vector with the bounds : [lower, upper]
%
%   PRODOPT contains the indices of a product whose average spectral radius
%           attains the lower bound (use with buildProduct)
%   ALLPROD contains on each row the indices of the products kept
%
%   INFO is a structure containing various data about the iterations :
%         info.status         - 0 if delta was achieved, 1 if maxEval was reached, 2 if strange happens
%         info.allLb          - lower bounds at each depth
%         info.allUb          - upper bounds at each depth
%         info.popHist        - number of candidates along depths
%
%  The field OPTS.grip (generated by jsrsettings) can be used to
%  set parameters and options :
%
%      delta         - aimed difference between bounds, influences the the number of kept products and hence the speed, (1e-2)
% Original-code by  R. Jungers.

starttime=clock;

M=varargin{1};
info=struct;
maxsmpdepth=parsem('maxsmpdepth',varargin,10); %maximum length of products to be computed
verbose=parsem({'verbose','v'},varargin,1);
maxeval=parsem('maxeval',varargin,inf); %maximum number of products to be evaluated
maxtime=parsem('maxtime',varargin,inf);
minJSR=parsem('minJSR',varargin,0);
searchonlyonecandidate= parsem('searchonlyonecandidate',varargin); %stops algorithm after one candidate has been found. To be used together with minJSR
[normfun,varargin]=                parsem('norm',varargin,[]); 
if(isempty(normfun)); 
    normfun=@(x) norm(x,2); end; %which norm to use
if(isnumeric(normfun)); 
    normfun=@(x) norm(x,normfun); end;
delta = parsem('delta',varargin,0.99); delta=1-delta;

% Parameters
m = length(M);

% Initialisation
neval = m;
l = 1;
nt = m;
prods(1:m,1) = (1:m)';
rhoM = trho(M);
[alpha, cand] = max(rhoM);

% Sometimes, scale can be equal to zero
if(alpha < eps)
    scale = eps;
else
    scale = alpha; end;
beta = 0;
P = zeros(1,m);
popHist = zeros(1,50);
popHist(1) = m;
allLb = zeros(1,50);
allUb = zeros(1,50);
% Note: comparison between for and cellfun was performed, for is quicker in
% majority of tests, see testCellfun_for.m
counter=numel(M); %number of computed matrices
for im=1:m
    %normi = feval(normfun,full(M{im}));
    normi = normfun(full(M{im}));
    if (normi>beta)
        beta = normi; end;
    P(im) = normi; end;

M = tcellDivide(M,scale);
T = M;

allLb(1) = alpha;
allUb(1) = beta;

if (m==1)
    info.jsrbound = [alpha, alpha];
    cand = {1};
    nearlycand={};
    info.status = 0;
    info.allLb = allLb(1:l);
    info.allUb = allUb(1:l);
    info.popHist = popHist(1:l);
    return; end;

% Main loop


while (neval<maxeval && beta > alpha + delta)
    
    if(~searchonlyonecandidate && l > maxsmpdepth); 
        break; end;
    if(etime(clock,starttime)>=maxtime); 
        break; end;
    
    
    l = l+1;
    if (l>length(popHist))
        popHist = [popHist zeros(1,50)]; %#ok<AGROW>
        allLb = [allLb zeros(1,50)]; %#ok<AGROW>
        allUb = [allUb zeros(1,50)]; %#ok<AGROW>
    end


    
    newT = cell(1,m*nt);
    newProds = sparse(m*nt,l);
    newP = sparse(1,m*nt);
    newnt = 0;
    maxNewP = alpha;
    vprintf('|','imp',[1, verbose]);
    
    alphaprev=alpha;
    counter=counter+m*nt;

    for it = 1:nt     
        
        for im = 1:m
            MT = M{im}*T{it};
            %normi = feval(normfun,full(MT))^(1/l)*scale; %original code
            normi = normfun(full(MT))^(1/l)*scale; %original code
            rhoMT = trho(MT)^(1/l)*scale;
            
            ilast = find(prods(it,:),1,'last');
            if isempty(ilast)
                ilast = 0;
            end
            prodi = prods(it,1:ilast);
            
            newPi = min(P(it),normi);
            
            % Actualise alpha
            if (rhoMT>alpha);
                alpha = rhoMT;
                cand = [prodi im];        
                vprintf('%v\n',full(prodi),'imp',[2,verbose])
            end
            
            % Keep or not
            if (newPi>alpha+delta)
                newnt = newnt+1;
                newT{newnt} = MT;
                
                % Add to newProds
                newProds(newnt,1:ilast+1) = [prodi im]; %#ok<SPRIX>
                % Actualise newP
                newP(newnt) = newPi; %#ok<SPRIX>
                
                % Keep max newP to actualise beta
                if (newP(newnt)>maxNewP)
                    maxNewP = full(newP(newnt));
                end
                
            end
        end
    end
    %  % Rescale the different sets
    if (alpha>alphaprev)
        newScale = alpha;
        if(newScale < eps)
        % Sometimes scale can be equal to zero
            newScale = beta;
        end
        
        M = tcellDivide(M,newScale/scale);
        %T = tcellDivide(T,(newScale/scale)^(l-1));
        newT = tcellDivide(newT,(newScale/scale)^(l));
        scale = newScale;
    end
    
    %vprintf('Number of products kept : %d\t',newnt,2,verbose);
    neval = neval + m*nt;
    nt = newnt;
    T = newT(1:nt);
    prods = newProds(1:nt,:);
    P = newP(1:nt);
    beta = min(beta,max(alpha+delta,maxNewP));
    
    popHist(l) = nt;
    allLb(l) = alpha;
    allUb(l) = beta;
    if(searchonlyonecandidate && ~isempty(cand) && alpha>minJSR); 
        break; end;
    
    if(abs(round(log2(l))-log2(l))<max(verbose*0.05,.3) || verbose>=5)
        vprintf('Depth %d, starting computation of %d new products. \t JSR =  [%15.12g, %15.12g]\n',l,m*nt,alpha, beta,'imp',[2,verbose]);
    end
end

if (beta-alpha)<=delta
    status = 0;
else
    status=1;
end

% Post-processing
cand = full(cand);
if(isempty(cand)); 
    alpha=0; end;
info.count=counter;
info.jsrbound = [alpha, beta];

info.status = status;
info.allLb = allLb(1:l);
info.allUb = allUb(1:l);
info.popHist = popHist(1:l);
info.count = neval;
nearlycand={};
cand={cand.'};


vprintf('Computed matrices in total (approximate): %i\n','cpr',[.1,0.5,.1],counter,'imp',[2,verbose]);

end

function [cand, nearlycand, info] = tjsr_bruteForce(varargin)
%
M=varargin{1};
minsmpdepth=parsem('minsmpdepth',varargin,1);
maxsmpdepth=parsem('maxsmpdepth',varargin,5);
minJSR=parsem('minJSR',varargin,0);
searchonlyonecandidate=parsem('searchonlyonecandidate',varargin);
verbose=parsem({'verbose','v'},varargin,1);
epsilon=parsem('epsilon',varargin,eps);
maxtime=parsem('maxtime',varargin,inf);
maxeval=parsem('maxeval',varargin,inf);
normfun=parsem('norm',varargin,[]); 
if(~isempty(normfun)); 
    vprintf('''bruteforce'' does not compute norms.\n','imp',[0 verbose],'cpr','err'); end;

% Initialization
JSR=minJSR;
cand={};
nearlycand={};
if(searchonlyonecandidate); 
    maxsmpdepth=inf; 
else
    vprintf([repmat('.',[1 maxsmpdepth-minsmpdepth+1]) '\n\n'],'imp',[1,verbose]); end;
starttime=clock;

j=minsmpdepth-1;
counter=0;
while(true)
    j=j+1; %current length of products
    if(~searchonlyonecandidate && j>maxsmpdepth); 
        break; end;
    if(etime(clock,starttime)>=maxtime); 
        break; end;
    order=tgenNecklaces(j,length(M))'+1; %XX computation of orderings is very expensive and must be accelerated
    for i=1:size(order,2);
        counter=counter+1;
        if(counter>maxeval); 
            break; end;
        o=order(:,i);
        if(~isequal(reducelength(o),o)); 
            continue; end; %%tgenNecklaces should be modified so that is does not return such values XX
        candidate=tbuildproduct(M,o);    
        R=trho(candidate)^(1/length(o));  %r is the spectral radius of the candidate %length_m

        if(JSR<R*(1-epsilon));  %make candidates to nearlycandidates, if their spectral radius is smaller than the maximum
            nearlycand=cand; 
            cand={}; 
            JSR=R; 
            vprintf('New candidate found: %v\n',o,'imp',[2,verbose]); end;

        if(R>=JSR); %JSR is the spectral radius to be outnumbered = current lower/upper bound for JSR/LSR
            JSR=R; end; 
        if(R>=JSR*(1-epsilon)); 
            cand{end+1}=o; %#ok<AGROW>
            if(searchonlyonecandidate);
                info.jsrbound=[JSR inf];
                return; end; end; end;
    if(counter>maxeval); 
        break; end;
    vprintf('|','imp',[1,verbose]); end;


info=struct;
info.jsrbound=[JSR max(estimatejsr(M))];
info.count = counter;
end



function [cand, nearlycand, info] = tjsr_genetic(varargin)
%[cand, nearlycand, info] = tjsr_genetic(M, [options])
%
%Inputs:
%   M      -  a non-empty cell array containing square matrices of identical dimensions.
%
% Options:
%       'verbose',val    (default 1) defines the printing level:
%                        <  1 - no output is printed on the screen;
%                        >= 1 - prints the final lower bound and the associated product at the end of the algorithm;
%                        >= 2 - prints the current lower bound and the associated product at each generation of the algorithm;
%                        >= 3 - prints the current maximum product length and the number of stalling iterations at this current maximum product length at each generation.
%      'popsize',val       (default 150, minimum 10) population size.
%      'maxgen',val        (default 500) maximum number of generations.
%      'maxstall',val      (default 70) maximum number of stalling iterations before increasing the maximum product length.
%      'maxtotstall',val   (default 200) maximum number of stalling iterations before terminating the algorithm.
%                        
%      'mutantprop',val    (default 0.3) mutation probability of a given product.
%      'muteprop',val      (default 0.2) mutation proportion of a given product.
%       
% Output:
%   bound               a lower bound on the jsr of M.
%   c                   a product of matrices corresponding to this lower bound:
%                       the product A1*A2*A3 is expressed as [3 2 1] in this order.
%   nc                  Empty
%   info.elapsedtime    cpu time used to run the algorithm.
%   info.population     the complete population at the end of the algorithm.
%
% Written by Chia-Tche Chang, 2011
% Interface modified by tommsch


M=varargin{1};

searchonlyonecandidate= parsem('searchonlyonecandidate',varargin); %stops algorithm after one candidate has been found. To be used together with minJSR
minJSR=             parsem('minJSR',varargin,0); 
verbose=            parsem({'verbose','v'},varargin,1);
POPSIZE=            parsem('popsize',varargin,500);
MAXGEN=             parsem('maxgen',varargin,200);
MAXSTALLING=        parsem('maxstall',varargin,70);
MAXTOTALSTALLING=   parsem('matotstall',varargin,1000);
MUTANTPROP=         parsem('mutantprop',varargin,0.3);
MUTEPROP=           parsem('muteprop',varargin,0.2);
normfun=            parsem('norm',varargin,[]); 
if(~isempty(normfun)); 
    vprintf('''genetic'' does not compute norms.\n','imp',[0 verbose],'cpr','err'); end;
maxtime=            parsem('maxtime',varargin,inf);
maxeval=            parsem('maxeval',varargin,inf);
starttime=clock; %starttime of tommsch

vprintf('The genetic algorithm has a bug, and sometimes reports spectral radii which are normalized wrongly!\nUse this algorithm with care.\n','cpr','err','imp',[0 verbose]);

% Initialization
STARTTIME = cputime; %starttime of Chiang
m = numel(M);
n = size(M{1}, 1);
k = ceil(log(POPSIZE)/log(m+1));
scaler = ((m+1).^(k-1:-1:0))';

% Cache
cache = tjsr_genetic_genCache(M, k, m, n);
ncache = size(cache, 2);
stalling = 0;
totalstalling = 0;
bound=minJSR^k;
bestpop = [];
for i = 2:ncache;
    value = max(abs(eig(mat(cache(:, i)))));
    if (value > bound);
        bound = value;
        bestpop = zeros(1, k);
        key = i-1;
        for j = 1:k;
            bestpop(j) = floor(key/scaler(j));
            key = key - bestpop(j)*scaler(j);
        end
    end
end
bestpop = bestpop(bestpop ~= 0);
bound = bound^(1/k); %XX Here is a bug. 
                        %T=overlapMat needs this line. 
                        %For a=1/12*[3 3 4 3 3 4 3 3 4 3 3]';M=-3;D=[-2 -1 0];S=getS('a',a,'M',M,'D',D);[T,Om]=transitionmatrix(S);U=constructU(T,1,'sym');T=restrictmatrix(T,U);
                        %the line is wrong.
X = eye(n);
testbestpop = bestpop;
for i = 1:length(testbestpop)-1;
    X = X * M{testbestpop(i)};
    testval = max(abs(eig(X)))^(1/i);
    if (testval >= bound);
        bestpop = testbestpop(1:i);
        bound = testval;
    end
end
if (verbose >= 2);
    fprintf(' \n');
    fprintf('Starting population: init lower bound on the JSR = %.15g with product: %s\n', tif(isempty(bestpop),0,bound), num2str(bestpop));
end

% Initial population
CURLENGTH = 2*k;
POPULATION = floor((m+1)*rand(POPSIZE, CURLENGTH));

% Genetic evolution
gen=0;
counter=0;
while(true)
    counter=counter+gen*MAXSTALLING;
    if(counter>maxeval); 
        break; end;
    gen=gen+1; 
    if(~searchonlyonecandidate && gen>=MAXGEN); 
        break; end;
    % Evaluation
    if(searchonlyonecandidate && ~isempty(bestpop)); 
        break; end;
    if(etime(clock,starttime)>=maxtime); 
        break; end;
    [nbeyes, idx] = sort(sum(POPULATION == 0, 2), 'descend');
    POPULATION = POPULATION(idx, :);
    POPULATION((nbeyes == CURLENGTH), 1) = ceil(m*rand(1));
    nbeyes = min(nbeyes, CURLENGTH-1);
    fitness = zeros(POPSIZE, 1);
    cached = floor(CURLENGTH/k);
    cachekey = zeros(POPSIZE, cached);
    for j = 1:cached;
        cachekey(:, j) = POPULATION(:, k*j-k+1:k*j)*scaler;
    end
    cachekey = cachekey + 1;
    for i = 1:POPSIZE;
        X = eye(n);
        for j = 1:cached;
            X = X * mat(cache(:, cachekey(i, j)));
        end
        for j = k*cached+1:CURLENGTH;
            if (POPULATION(i, j) ~= 0);
                X = X * M{POPULATION(i, j)};
            end
        end
        fitness(i) = max(abs(eig(X)))^(1/(CURLENGTH-nbeyes(i)));
    end
    [fitness, idx] = sort(fitness, 'descend');
    POPULATION = POPULATION(idx, :);
    
    % Local optimization
    if (fitness(1) > bound);
        stalling = 0;
        totalstalling = 0;
        bound = fitness(1);
        bestpop = POPULATION(1, :);
        bestpop = bestpop(bestpop ~= 0);
        testbestpop = bestpop;
        lenbestpop = length(testbestpop);
        CELLX = cell(lenbestpop + 1, 1);
        X = eye(n);
        CELLX{1} = X;
        for i = 1:lenbestpop;
            X = X * M{testbestpop(i)};
            CELLX{i+1} = X;
            testval = max(abs(eig(X)))^(1/i);
            if (testval > bound) || ((testval == bound) && (i < length(bestpop)));
                bestpop = testbestpop(1:i);
                bound = testval;
            end
        end
        lenbestpop = length(testbestpop);
        X = eye(n);
        localimprove = 0;
        for i = lenbestpop:-1:1;
            testval = max(abs(eig(CELLX{i}*X)))^(1/(lenbestpop-1));
            if (testval > bound);
                bestpop = testbestpop([1:i-1 i+1:end]);
                bound = testval;
                localimprove = localimprove + 1;
            end
            for j = 1:m;
                testval = max(abs(eig(CELLX{i}*M{j}*X)))^(1/lenbestpop);
                if (testval > bound);
                    bestpop = testbestpop;
                    bestpop(i) = j;
                    bound = testval;
                    localimprove = localimprove + 1;
                end
            end
            if (lenbestpop < CURLENGTH);
                for j = 1:m;
                    testval = max(abs(eig(CELLX{i+1}*M{j}*X)))^(1/(lenbestpop+1));
                    if (testval > bound);
                        bestpop = [testbestpop(1:i)  j  testbestpop(i+1:end)];
                        bound = testval;
                        localimprove = localimprove + 1;
                    end
                end
            end
            X = M{testbestpop(i)}*X;
        end
        if localimprove;
            POPULATION = [bestpop, zeros(1, CURLENGTH-length(bestpop)); POPULATION(1:end-1, :)];
        end
    else
        stalling = stalling + 1;
        totalstalling = totalstalling + 1;
    end
    
    % Stopping criterion
    if (verbose >= 3);
        fprintf('Generation #%3d: STA = %2d, LEN = %2d, lower bound = %.15g with product: %s\n', gen, stalling, CURLENGTH, tif(isempty(bestpop),0,bound), num2str(bestpop));
    elseif (verbose >= 2);
        fprintf('Generation #%3d:  current lower bound on the JSR = %.15g with product: %s\n', gen, tif(isempty(bestpop),0,bound), num2str(bestpop));
    end
    if (~searchonlyonecandidate && totalstalling >= MAXTOTALSTALLING );
        break; end;
    if (stalling >=  MAXSTALLING);
        stalling = 0;
        CURLENGTH = CURLENGTH + 1;
        POPULATION = [POPULATION, zeros(POPSIZE, 1)]; %#ok<AGROW>
    end
    
    % Selection and crossover
    NB_ELITE = max(2, min(3, floor(POPSIZE/50)));
    NB_SPAWN = max(4, floor(POPSIZE/50));
    NB_SWAP = min(POPSIZE - NB_ELITE - NB_SPAWN, floor(POPSIZE/2));
    NB_MIX = POPSIZE - NB_ELITE - NB_SPAWN - NB_SWAP;
    ID_SWAP = ceil((POPSIZE/2)*rand(NB_SWAP, 2));
    ID_MIX = ceil(POPSIZE*rand(NB_MIX, 2));
    
    POP_ELITE = POPULATION(1:NB_ELITE, :);
    POP_SPAWN = ceil(m*rand(NB_SPAWN, CURLENGTH));
    POP_SWAP = zeros(NB_SWAP, CURLENGTH);
    for i = 1:NB_SWAP;
        cut = ceil(CURLENGTH*rand(1));
        POP_SWAP(i, :) = [POPULATION(ID_SWAP(i, 1), 1:cut) POPULATION(ID_SWAP(i, 2), cut+1:end)];
    end
    POP_MIX = zeros(NB_MIX, CURLENGTH);
    for i = 1:NB_MIX;
        mixer = floor(2*rand(CURLENGTH, 1))';
        POP_MIX(i, :) = POPULATION(ID_MIX(i, 1), :) .* mixer + POPULATION(ID_MIX(i, 2), :) .* (1-mixer);
    end
    
    POPULATION = [POP_ELITE ; POP_SPAWN ; POP_SWAP ; POP_MIX];
    
    % Mutation
    MUTESTR = ceil(CURLENGTH * MUTEPROP);
    for i = 2:POPSIZE;
        if (rand(1) < MUTANTPROP);
            POPULATION(i, ceil(CURLENGTH * rand(1, MUTESTR))) = floor((m+1)*rand(1, MUTESTR));
        end
    end
end

% Termination
elapsedtime = cputime - STARTTIME;
%bestpop = tjsr_genetic_deperiod(bestpop); %we use the function simplify_ordering instead

% POST PROCESSING FOR tommsch-Interface GENERATION  %
if(isempty(bestpop)); 
    bound=0; end;
cand={flip(bestpop).'}; %change order and make it to column vector, since tommsch-programs need that
nearlycand={};
info.time=elapsedtime;
info.population=POPULATION;
info.jsrbound=[bound inf];
info.count=counter;

if (verbose >= 1)
    fprintf('Algorithm terminated with lower bound on the JSR = %.15g with product: %s\n', bound, num2str(bestpop));
end

end

function cache = tjsr_genetic_genCache(M, k, m, n) %  CACHE GENERATION  %
    if (k <= 1);
        cache = zeros(n*n, m+1);
        cache(:, 1) = vec(eye(n));
        for i = 1:m;
            cache(:, i+1) = vec(M{i});
        end;
        return;
    end
    cache = repmat(tjsr_genetic_genCache(M, k-1, m, n), 1, m+1);
    gsize = (m+1)^(k-1);
    for i = 1:m;
        cache(:, gsize*i+1:gsize*(i+1)) = kron(eye(n), M{i}) * cache(:, gsize*i+1:gsize*(i+1));
    end
end

% function [x, k] = tjsr_genetic_deperiod(X) %  DEPERIODIZATION  %
% n = length(X);
% for k = 1:n/2;
%     if (mod(n, k) == 0);
%         ok = 1;
%         for t = 1:n-k;
%             if X(t) ~= X(t+k);
%                 ok = 0;
%                 break;
%             end
%         end
%         if ok;
%             x = X(1:k);
%             return;
%         end
%     end
% end
% x = X;
% k = n;
% end

function cand = simplify_ordering(cand)
    %Simplify candidates and remove duplicates
    if(~isempty(cand));
        cand=cellfun(@reducelength,cand,'UniformOutput',0); %reduce length
        LENGTHMAX=max(cellfun(@(x) max(size(x)),cand))+1;   
        for i=1:size(cand,2); 
            cand{i}(LENGTHMAX,1)=0; end;  
        cand=cell2mat(cand);                                
        cand(end,:)=[];                                     
        cand=unique(cand','rows')';                        
        cand=num2cell(cand,1);
        cand=cellfun(@(x) removezero(x,1),cand,'UniformOutput',0); end;
end

function nc = remove_ordering(nc, c)
    %removes orderings in nc, which are also present in c
    idx=size(nc);
    for i=1:length(nc)
        idx(i)=searchincellarray(nc{i},c,1); end;
    nc(i)=[];
end

function v = vec(M)
% Converts matrices to vectors
% Input: %   M       matrix
% Output: %   v       vector
    v = reshape(M,numel(M),1);
end

function M = mat(v,n)
% Input:    v       vector of length 2^n
%          [n]     (optional) the value of n.
%                  If n is not given, it is computed
% Output:   M       the converted matrix
if(nargin < 2);
    n = floor(sqrt(length(v))); end;
M = reshape(v,n,n);
end



function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   
