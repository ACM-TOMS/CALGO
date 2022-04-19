function varargout = tile( varargin)
% [ Q, oo ] = tile( [oo], S, [options])
% Plots the tile(attractor) corresponding to a multiple subdivision scheme.
% The scheme is defined by the dilation matrices (S{:,2}) and digit sets (S{:,3}) and the ordering oo
%  
% Input: 
%   oo                              {[1xN],[1xM]}, Defines the ordering in which the subdivision operators are applied.
%                                   If S consists of more than one subdivision operator, and this argument is not given, then a random ordering is used
%   S                               anything which returns subdivision operators when called as getS(S)
%
% Options:  
%       'supertile'                 Computes the supertile instead of the attractor.
%       'verbose',val               default=1, verbose level
%       'peter',peterval:           [Not recommend. Use <'round',roundval> instead. After each step, peter(Q,peterval) is called. peterval -> see roundval
%       'round',roundval            After each step, points which are 1/roundval near to each other are removed.
%                                   roundval is applied before peterval
%                                   'roundval':  1 x size(ordering,2)-vector: gives the round values in each step
%                                       1 x 2 vector: starts with roundval(1), last step is rounded to roundval(end), in between linear interpolation is applied
%                                       scalar: same round value in each iteration
%       'start', array              default=zeros(dim,1), starting set
%       'digit'                     computes digits instead of the attractor
%       'diffdigit'                 computs 'digits minus digits' instead of the attractor
%       'plot',cellarray            default={}, cell array of arguments passed to plotm
%                                   if cellarray==0, nothing is plotted.
%       'iteration',val             default: depends on the input, Makes val iterations.
%       'numpoint',val              default=30000, Estimate of how many points the outcome shall have
%       'maxiteration',val          default=50, Maximum number of iterations. 
%       'interval'                  (Experimental) Uses interval arithmetic. Works only for dim==1
%       'supp'                      (Experimental) Calls getS with option 'supp'. Thus, most likely the support of the blf will be plotted
%       'OmegaRR'                   (Experimental) Calls getS with option 'OmegaRR'. Thus, most likely the set OmegaRR will be plotted
%
% Output: 
%       Q                           computed attractor
%       oo                          used ordering.
%                                                   
% E.g.: tile('2_frayed_squares','round',[.1 1e-2])
%       tile('1_cantor','digit')
%       tile([getS('2_rand'); getS('2_rand')],'supertile','iteration',10)
%
% See also blf, getS
%
% Written by: tommsch, 2018
% For more information write to: <a href="tommsch@gmx.at">tommsch@gmx.at</a>

%#ok<*NOSEM>
% XX 'interval' for arbitrary dimensions

[interval,varargin] = parsem('interval',varargin);
[Q,varargin] = parsem('start',varargin,[]);
[peterval,varargin] = parsem('peter',varargin,0);
[plotval,varargin] = parsem('plot',varargin,{});
[verbose,varargin] = parsem({'verbose','v'},varargin,1);
[digit,varargin] = parsem('digit',varargin);
[diffdigit,varargin] = parsem('diffdigit',varargin); if(diffdigit); digit=1; end;
[roundval,varargin] = parsem('round',varargin,tif(digit || diffdigit,0,[.1 1e-2]));
[iteration,varargin] = parsem('iteration',varargin,[]);
[maxiteration,varargin] = parsem('maxiteration',varargin,50);
[supertile,varargin] = parsem('supertile',varargin);
[suppflag,varargin] = parsem('supp',varargin);
[OmegaRR,varargin] = parsem('OmegaRR',varargin);
[numpoint,varargin] = parsem('numpoint',varargin,30000);

if(size(varargin,2)>=2 && ~isS(varargin{1}));
    if(suppflag); 
        S=getS(varargin{2},'supp'); 
    elseif(OmegaRR); 
        S=getS(varargin{2},'OmegaRR'); 
    else; 
        S=getS(varargin{2}); end;
    oo=varargin{1}; 
%     if(~isordering(oo))
%         oo=constructordering(oo,[]);
%     end
    varargin(1:2)=[];
else
    if(suppflag); 
        S=getS(varargin{1},'supp'); 
    elseif(OmegaRR);  
        S=getS(varargin{1},'OmegaRR'); 
    else; 
        S=getS(varargin{1}); end;
    oo=constructordering(size(S,1),0,'random',1000);
    varargin(1)=[];
end

M=S(:,2);D=S(:,3);
dim=size(M{1},1);
J=size(S,1); %number of subdivision operators
if(isempty(Q)); 
    Q=zeros(dim,1); end;

parsem(varargin,'test');

%compute how often we shall iteratie
if(isempty(iteration)); %compute how many points the outcome will have    
    iteration=1;
    while(true)
        if(supertile); 
            val=abs(prod(cellfun(@(x) det(x.^iteration),S(:,2))))*size(Q,2);;
        else; 
            val=abs(det(tbuildproduct(S(:,2),ordering2vector(oo,iteration))))*size(Q,2); end;
        if(val>numpoint || iteration>maxiteration); 
            iteration=iteration-1; break; 
        else; 
            iteration=iteration+1; end; end; end;
oo=ordering2vector(oo,iteration);

if(size(oo,2)<iteration);
    warning( 'tile:oolength', 'Length of ordering shorter than iteration.' );
    iteration=size(oo,2); end;

if(interval)
    vprintf(['Warning: Experimental option ''interval'' used. Works only for dim==1.\n' ...
             'You MUST set ''start'' and''round'' by hand.\n' ...
             '''start'' should probably be {[0 1]}, ''round'' should probably be 0.\n'],'cpr',[.7 .5 0]','imp',[0 verbose]); end;

if(supertile); 
    if(digit); 
        warning( 'tile:options', 'Options ''digits'' and ''supertile'' do not work together.' ); end;
    oo=1:J; end;
l=iteration;

vprintf('Iterations:: %i\n',l,'imp',[1 verbose]);
vprintf('ordering: %v\n',oo,'imp',[1 verbose]);
vprintf('dim: %i\n',dim,'imp',[2 verbose]);

if(~digit); 
    for i=1:J; 
        M{i}=inv(M{i}); end;
else; 
    oo=fliplr(oo); end;


if(supertile); 
    vprintf([repmat('.',1,l) '\n\n'],'imp',[2 verbose]);
else; 
    vprintf([repmat('.',1,size(oo,2)) '\n\n'],'imp',[2 verbose]); end;

if(isscalar(roundval));         
    roundval=roundval*ones(1,l);
elseif(size(roundval,2)==2);    
    roundval=linspace(roundval(1),roundval(2),l); end;
roundval=fliplr(roundval);

if(isscalar(peterval));         
    peterval=peterval*ones(1,l);
elseif(size(peterval,2)==2);    
    peterval=linspace(peterval(1),peterval(2),l); end;
peterval=fliplr(peterval);


for k=l:-1:1 %iterate iter-times
    tic;
    if(supertile); 
        Qnew=zeros(1,dim);
        for n=oo;
            if n==0; 
                continue; end;
            Qnew=union(Qnew,(M{n}*setplus(Q,D{n})).','rows'); end;
        Q=Qnew';
        clear Qnew;
    elseif(interval)
        if(oo(k)==0); 
            continue; end;
        parfor kk=1:numel(Q)
            Q{kk}=M{oo(k)}*setplus(Q{kk},D{oo(k)},'nounique'); end; %#ok<PFBNS>
        Q=mergeinterval([Q{:}],'output','C'); %merge all intervals
    else; 
        if(oo(k)==0); 
            continue; end;
        if(digit);
            Q=setplus(M{oo(k)}*Q,D{oo(k)});
        else
            Q=M{oo(k)}*setplus(Q,D{oo(k)}); end;
    end;
    if(roundval(k));
        [~,idx,~]=unique((round(Q./roundval(k))).','rows');
        Q=Q(:,idx); end;
    if(peterval(k));
        Q=peter( Q , peterval(k) ); end;
    vprintf('\b|\n','imp',[2 verbose]);
    if(toc>2 && ~roundval(k) && ~interval); 
        vprintf('Use <''round'',val> if computation takes too long.\n','imp',[1 verbose]); end; end;

if(diffdigit); 
    Q=setplus(Q,-Q); end;

%output
if(dim==1 && verbose>=1); 
    if(interval); 
        MIN=minm([Q{:}]); 
        MAX=maxm([Q{:}]);
        val1=sum(diff(vertcat(Q{:}),1,2));
        val2=sum(abs(diff(vertcat(Q{:}),1,2)));
        if(val1~=val2);
            fprintf('Area (oriented/unoriented): %i\n',val1,val2);
        else
            fprintf('Area: %i\n',val1); end;
        fprintf('Tile min:%i, Tile max:%i\n', MIN, MAX);
    else; 
        fprintf('Tile min:%i, Tile max:%i\n', min(Q), max(Q)); end; end;

if(~isequal(plotval,0))
    plotm(Q,'.','MarkerSize',1,plotval{:}); end;

if(nargout>=1); 
    varargout{1}=Q;  end;
if(nargout>=2); 
    varargout{2}=oo; end;

end

% function [y,x]=swap(x,y)
% end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.  