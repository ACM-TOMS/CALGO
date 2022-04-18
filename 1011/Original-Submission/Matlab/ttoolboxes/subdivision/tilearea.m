function [ flag ] = tilearea(varargin)
% [ flag ] = tilearea( Q | S, [options] )
% Good heuristic test if an attractor has area one. 
%
% Q has probably area one if
%   1) the green dots have the same spacing as the red dots, and
%   2) the reported percentage of possibly duplicates is small.
% This is just a heuristic test. The function digittest can answer the question rigorously, but it has memory problems and is broken 
%
% Input:
%   Q               dim x N array, the tile
%                       OR
%   S               subdivision operators
%   
%
% Options:
%   'visual'        Makes plot output to test it visually. Works only for dim=1 and 2. May easily be changed to dimension 3.
%                   This is much faster, but needs human interaction, since it plots only some tiles.
%   'verbose'       Verbose level.
%   'numpoint',val  Approximate number of points (only applicable if S is given)
%
% Output: 
%   flag            1 if area is likely to be one
%                   0 if area is likely to be greater one
%                   If 'visual' is used, this flag always equals -1
%   Graphical and text output to be examined
%
% E.g.: tilearea(tile('2_frayed-squares','round',[.1 1e-2],'plot',0));
%       tilearea(tile(getS('M',2,'D',[0 3]),'plot',0));
%
% See also: checktile
%
% Written by: tommsch, 2016


[visual,varargin]=parsem('visual',varargin);
[numpoint,varargin]=parsem('numpoint',varargin,tif(visual,100000,10000));

if(isS(varargin{1})); 
    varargin{1}=tile(varargin{1},'plot',0,'numpoint',numpoint,'round',0,'verbose',0);
end

if(visual);
    tilearea_visual(varargin{:})
    flag=-1;
else
    flag=tilearea_numerical(varargin{:});
end

end

function [ flag ] = tilearea_numerical(varargin)

[verbose,varargin]=parsem({'verbose','v'},varargin,1);

Q=varargin{1};
sze=size(Q,2);
dim=size(Q,1);
accuracy=10e12;
flag=1;

%compute minimum difference
DIFF=zeros(1,sze);
Q=Q.';
if(sze>10000); 
    vprintf('Computing %i distances. This may take a while.\n', sze,'imp',[1 verbose]); end;
parfor i=1:sze
    val=pdist2(Q,Q(i,:),'minkowski','sm',2);
    DIFF(i)=val(2);
end

QSMALL=Q;
for i=1:dim
    idx=Q(:,i)>=1;
    QSMALL(idx,i)=QSMALL(idx,i)-floor(QSMALL(idx,i));
    
    idx=Q(:,i)<0;
    QSMALL(idx,i)=QSMALL(idx,i)-floor(QSMALL(idx,i));
end

if(sze>10000); 
    vprintf('Computing %i distances again. This may take the same while.\n', sze,'imp',[1 verbose]); end;
DIFFSMALL=zeros(1,sze);
parfor i=1:sze
    val=pdist2(QSMALL,QSMALL(i,:),'minkowski','sm',2);
    DIFFSMALL(i)=val(2);
end


if(median(DIFF)>1.5*median(DIFFSMALL))
    flag=0;
end

Qint=int64(QSMALL*accuracy); %change to int64 to be able to compare the values
Qint=unique(Qint, 'rows')';
percentage=(size(Q,2)-size(Qint,2))/size(Q,2);
if(percentage > 0.25)
    flag=0;
end

if(flag)
    vprintf('Tile may have area equal one.\n','cpr',[0 .4 0],'imp',[1 verbose]);
else
    vprintf('Tile may have area larger than one.\n','cpr','err','imp',[1 verbose]);
end

end

function [ ] = tilearea_visual(varargin)
%tilearea(Q,['accuracy',val])
%Tests if the set Q has area one.
%If yes, then 
% 1) the green dots should have the same spacing as the red dots
% 2) the matlab output should say Q Contains probably 0 duplicates. (or a
%    very small number compared to the size of Q.
%This is just a heuristic test. The function digittest can proof the
%question, but has memory problems.
%figure
clf
hold on
Q=varargin{1};
accuracy=parsem('accuracy',varargin,10*(size(Q,2)));
verbose=parsem({'verbose','v'},varargin,1);

if size(Q,1)==1;
    %plot original data 
    vprintf('Minimum disctance Q=%d\n', min(diff(sort(Q))),'imp',[1 verbose]);
    plot(Q,ones(1,size(Q,2)), 'ro', 'LineWidth',.01);
    %plot(Q, 'ro','LineWidth',.01);
    
    %shift everything into the unit square
    Q(Q(:)>=1)=Q(Q(:)>=1)-floor(Q(Q(:)>=1));
    Q(Q(:)<0)=Q(Q(:)<0)-floor(Q(Q(:)<0));
    plot(Q,ones(1,size(Q,2)), 'b.', 'LineWidth',.01);
    hold off
    %check for duplicates
    
    Qint=int64(Q*accuracy); %change to int64 to be able to compare the values
    Qint=unique(Qint', 'rows')';
    numdupl=size(Q,2)-size(Qint,2);
    vprintf('accuracy=%d, Q Contains probably %d duplicates.\n',accuracy,numdupl,'imp',[1 verbose]);
    vprintf('Minimum distance Q_shifted=%d\n', min(diff(sort(Q))),'imp',[1 verbose]);
    %vprintf('Q Contains points: %i\n',size(Q,2),'imp',[1 verbose]);

elseif size(Q,1)==2;
    %plot original data
    plotm(Q,'ro');
     
    %shift everything into the unit square
    Q(1,Q(1,:)>=1)=Q(1,Q(1,:)>=1)-floor(Q(1,Q(1,:)>=1));
    Q(2,Q(2,:)>=1)=Q(2,Q(2,:)>=1)-floor(Q(2,Q(2,:)>=1));
    Q(1,Q(1,:)<0)=Q(1,Q(1,:)<0)-floor(Q(1,Q(1,:)<0));
    Q(2,Q(2,:)<0)=Q(2,Q(2,:)<0)-floor(Q(2,Q(2,:)<0));
    plotm(Q,'b.');
    axis equal;

    %check for duplicates
    if nargin==1; 
        accuracy=4*(size(Q,2))^(1/2); end;
    Qint=int64(Q*accuracy); %change to int64 to be able to compare the values
    Qint=unique(Qint', 'rows')';
    numdupl=size(Q,2)-size(Qint,2);
    vprintf('accuracy=%d, Q Contains probably %d duplicates.\n',accuracy,numdupl,'imp',[1 verbose]);
    %vprintf('Q Contains points: %i\n',size(Q,2),'imp',[1 verbose]);
end;

hold off

end