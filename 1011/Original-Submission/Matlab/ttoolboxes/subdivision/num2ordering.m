function [do, doint, num, proofnum] = num2ordering(varargin)
% [ do, doint, numout, [proofnum] ] = num2ordering( [oo], S, num, [options] )
% Computes number expansions for multiple multivariate number systems.
% Takes a number in \RR^dim and a subdivision scheme (thus a number-system) and computes the corresponding ordering
%
% Input:
%   oo              Ordering ( {[1xN], [1xM]} ) defining the ordering of the subdivision operators. 
%   S               Cell array of subdivision operators
%   num             dim x 1 vector, number to compute
%
% Options:
%   'check',val     Checks if the given number is inside of the tile. Default=false
%   'plot'          Graphical output. Default=0
%   'l',val         Maximal length of output sequence. Default=automatic
%   'round',val     Parameter used in the computation of the attractor. See 'tile' for explanation: Default: <'round',[.1 1e-8]>  . To increase accuracy change the value of 'round', e.g. to <'round',[.2e-5 1e-12]>  .
%   'proof'         Computes the value of the number right to the radix point using the computed digits
%   'sym'           (experimental) symbolic computation
%   'delta'         Maximum distance of points to attractor until error is reported
%   'numpoint'      The number of points of the more exact attractor.
%
% Output: 
%   do              The digits as index numbers corresponding to the entries in S{:,3} for the number numout Format ordering ( {2xN ,2xM} )
%   doint           dim x l sequence of digits (rigth to the radix point) as integers
%   numout          First column: The number which was used for computation (may be different from the input number). All other columns: intermediate values which occured during the computation
%   [proofnum]      Only set with 'proof'. The numeric value of the computed digit sequence
%
% E.g.: S=getS({[],[1 2; -1 2]}); ordering2num(num2ordering(S,[.24;.3]),S)
% 
% See also: num2ordering>fullhelp, ordering2num, getS, peter, tile, tile
%
% Written by: tommsch, 2018

% XX Implement Symbolic computation
% XX Erkennen ob Output periodisch wird (nur bei symbolischer Berechnung). Rest überprüfen, und merken. Dann direkt Periode setzen
% XX Going back to old digits does not work properly

%#ok<*ALIGN>

parg=varargin;
[l,parg]=parsem('l',parg,0);
[verbose,parg]=parsem({'verbose','v'},parg,1);
[proof,parg]=parsem('proof',parg);
[plotflag,parg]=parsem('plot',parg,0);
[delta,parg]=parsem('delta',parg,.5); %max distance of points to attractor until error is reported
[check,parg]=parsem('check',parg,0);
[roundval,parg]=parsem('round',parg,[.01 1e-8]);
[symflag,parg]=parsem('sym',parg);
[numpoint,parg]=parsem('numpoint',parg,1e5);

if(~isS(parg{1}));
    S=parg{2};
    num=parg{3};
    oo=parg{1}; 
    if(~isordering(oo)); 
        oo=constructordering(oo); end;
    parg(1:3)=[];
else
    S=parg{1};
    if(size(S,1)>1); 
        error( 'num2ordering:argument', 'num2ordering: You must give ''oo'' as first argument, since there is more than one subdivision operator in ''S''.' ); end;
    num=parg{2};
    oo=constructordering([],1);
    parg(1:2)=[]; end;

parsem(parg,'test');

S=getS(S);
M=S(:,2);
D=S(:,3);
dim=size(M{1},2);
sizeD=max(cellfun(@(x) size(x,2),D)); 

if(~isequal(numel(num),dim)); 
    error( 'num2ordering:dim', 'Dimension of point wrong.\n'); end;

if(symflag); 
    num=sym(num); end;
num(1:dim,1) = num;
num(1:dim,2:l+1) = zeros(dim,l);
if(isequal(l,0))
    if(symflag); 
        [~,DENOM]=numden(num); 
    else; 
        [~,DENOM]=rat(num); end;
    l=maxm(max(2*DENOM,15)); end;


do=zeros(sizeD,l);
doint=zeros(dim,l);

%compute tile
vprintf('Compute tile.\n','imp',[2,verbose]);    
Qfine=[]; %initialize Qfine. This set is computed later if necessary.
Q=tile(oo,S,'round',roundval,'plot',0,'verbose',verbose-1); %the tile
if(isempty(Q)); 
    error( 'num2ordering:empty', 'Attractor is empty: Increase <''round'',val>.'); end;

oov=ordering2vector(oo,l);

if(check);
    
    %Split up num into numr (with digits right to the radix point) and numl (with digits left to the radix point)

    DIST=ceil(max(sum(abs(setplus(Q,double(-num(:,1)))),1)));
    vprintf('Find part of the number left to the radix point: ','imp',[2, verbose]);  
    while(true)
        ZZ=mixvector(-DIST:DIST,dim);
        vprintf('.','imp',[2,verbose]);  
        d=distancetoattractor(setplus(ZZ,double(num(:,1))),Q);
        [dmin,idx]=min(d);    
        if(dmin<1); 
            break; end;
        DIST=DIST*2; end;
    vprintf('\n','imp',[2,verbose]);  
    numl=-ZZ(:,idx);
    numr=num(:,1)-numl;
    if( any(numl) && nargout<3 ); 
        warning( 'num2ordering:pointoutside', 'Point is (probably) not in the attractor. Compute digits for the point:\n' );
        disp(numr);
    else
        vprintf('Point is (probably) in the attractor.\n','imp',[2,verbose]); end;
    num(:,1)=numr; end;

SMALLEST=l; %keeps track how far we went back in the computation of the digits

d=inf*ones(sizeD,l);
%digits right to the radix point
vprintf('Compute digits right to the radix point for the number: %v\n',num(:,1),'imp',[2,verbose]);
i=0;
while(true);
    if(i==l); 
        break; end;
    i=i+1;
    DD=D{oov(i)};
    z=setplus(M{oov(i)}*num(:,i),-DD,'nounique'); %number minus all digits. admissible digits are inside the tile
    nDD=size(DD,2);

    if(do(1,i)==0);  %compute new digit, otherwise use existing digit
        d(:,i)=distancetoattractor(double(z),Q); %distance of the digit-many points to the attractor
        if(verbose>=3 && plotflag)
            clf; hold on;
            plotm(Q,'.'); plotm(z,'-o');  plotm(z(:,1),'gx');
            drawnow; end;
        [~, idx] = sort(d(:,i),'ascend'); %find digit which is nearest to the attractor
        if(d(idx(1),i)/d(idx(2),i)>1/4);  %if the digit which shall be chosen is not clear
            vprintf(',','imp',[2,verbose]);
            if(isempty(Qfine));
                vprintf('- Recompute tile -','imp',[2,verbose]);
                Qfine=tile(oo,S,'round',roundval./100,'plot',0,'verbose',verbose-1,'numpoint',numpoint); end; %the tile
            d(:,i)=distancetoattractor(double(z),Qfine); %we compute the distances with the fine attractor
            [~, idx] = sort(d(:,i),'ascend');
        else
            vprintf('.','imp',[2,verbose]); end;
        idx(nDD+1:end)=0; %only keep indices which correspond to the chosen index


    else %try different digits from older times
        idx=do(:,i); idx(idx==0)=[]; end;


    doint(:,i)=DD(:,idx(1)); %integer value of digit

    num(:,i+1)=setplus(M{oov(i)}*num(:,i),-doint(:,i)); %new number
    do(1:length(idx),i)=idx; %First row is the index of the chosen digit. Other rows are wrong digits, 
                   %If the expansion goes wrong, the algorithm goes back digit by digit, and uses the other digits from do.

    if(d(idx(1),i)>delta); %if distance is larger than delta, then the algorithm chose a wrong digit
        if(do(2,i)==0); 
            if(i==1); 
                warning( 'num2ordering:pointoutside', '\nGiven number not in the tile. Set <''check'',1> or increase <''delta'',val>.' );
                break; end;
            vprintf('!','imp',[2,verbose]); 
            do(:,i)=0;
            d(:,i)=inf;
            i=i-1;
            if(i<SMALLEST); 
                SMALLEST=i; vprintf('\nBack to digit %i : ',i,'imp',[3,verbose]);  end; end;

        do(1:end-1,i)=do(2:end,i);
        do(end,i)=0;
        i=i-1;
        %vprintf('%i',do(1,:),2,verbose); vprintf('\n',2,verbose);
        continue; end; end; %XX Is this continue really necessary?

vprintf('\n','imp',[2,verbose]);


%postprocessing
do=do(1,:);

do=removezero(do,'right'); %remove zero indices
doint=doint(:,1:size(do,2)); %make doint and do the same length
num=num(:,1:size(do,2));
oov=oov(1:length(do));
%try to find periodics
do=[oov;do];
vprintf('Search for period.\n','imp',[2,verbose]);
[do,w]=findperiod(do,'everywhere','verbose',verbose-1);
if(w<2); 
    warning( 'num2ordering:period', 'Periodic part of output probably wrong. You should set <''l'',%i> .', 2*l ); end;

if(proof)
    vprintf('Proof number.\n','imp',[2,verbose]);
    proofnum=ordering2num(do,S);
    prooferr=norm(proofnum-num(:,1));
    vprintf('Error: %i\n',prooferr,'imp',[1 verbose]); end;


%plot verbose output
vprintf('\n','imp',[2,verbose]);
if(plotflag);    
    dim=size(D{1},1);
    if(dim==1); 
        val(2,:)=(1:size(num,2))/size(num,2); 
    else; 
        val=num; end;
    clf; hold on;
    if(isempty(Qfine)); 
        plotm(Q,'.','MarkerSize',1); else; plotm(Qfine,'.','MarkerSize',1); end;
    
    plotm(val(:,1),'g.','MarkerSize',20); %mark start point in red and green
    plotm(val(:,1),'go','MarkerSize',20);
    plotm(val,'ro-'); end;

num=num(:,1);

end

function [d]=distancetoattractor(p,Q)
% d = distancetoattractor(p,Q)
% computes the 1-distance from p to the nearest point in Q

d=zeros(1,size(p,2));
for i=1:size(p,2);
    %d(i)=min(sqrt(sum(abs(setplus(Q,-p(:,i))).^2,1))); end; %2-norm
    d(i)=min(sum(abs(setplus(Q,-p(:,i))),1)); end;
end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 