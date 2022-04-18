
function [t0,x0,y0,z0,zx0,zy0,zxx0,zxy0,zyy0]=...
                                  plotArgyris(t,x,y,s,varargin)

% PLOTARGYRIS Plot an argyris finite element
% 
%
% [T0,X0,Y0,Z0]=PLOTARGYRIS(T,X,Y,S);   
% 
% Plot the Argyris function defined on the grid {T,X,Y} and given by S.
% T is a nT x 21 matrix (nT is the number of triangles), X,Y are 
% vectors  with the coordinates of the nodes and S is a vector with the 
% coordinates of the Argyris element.
%
% T0,X0,Y0 especify a grid resulting of performing a uniform (red)
% refinment of the initial grid where the Argyris element is defined.
% Z0 has the values  of the Argyris element on the grid {T0, X0,Y0} 
%
% [T0,X0,Y0,Z0]=PLOTARGYRIS(T,X,Y,Z,NREF);   
%        
% NREF is the number of uniform (red) refinements of the initial grid 
% used to construct {T0,X0,Y0,Z0} (NREF=1 if not specified).            
%
%
% [T0,X0,Y0,Z0,ZX0,ZY0,ZXX0,ZXY0,ZYY0]=PLOTARGYRIS(T,X,Y,Z,NREF);
%
% Returns in addition  ZX0, ZY0, ZXX0, ZXY0, ZYY0, 
% the values of the first and second derivatives of the Argyris 
% element on the grid {T0,X0,Y0} and plots the corresponding graphics. 
%
% [T0,X0,Y0,Z0,ZX0,ZY0,ZXX0,ZXY0,ZYY0]= PLOTARGYRIS(T,X,Y,Z,NREF,'N')
%
% Do not plot the pictures. 



s=s(:).';
nref=1;
drawQ='Y';
if nargin>4
    if ~isempty(varargin{1})
        nref=varargin{1};
    end
    if nargin>5 
        if isequal(upper(varargin{2}),'N')
            drawQ='N'; % We plot the result.
        end
    end
end

% nref
% 1 -> 1
% 2 -> 5
% 3 -> 21


nrefaux=(4^(nref)-1)/3; 


nP=max(max(t(:,1:3)));

% construct the grid where the Agyris element is
% goint to be evaluated
nT=size(t,1);
t0=zeros(nT*4^nref,3);
x0=zeros(1,nP*4^nref);
y0=zeros(1,nP*4^nref);
x0(1:nP)=x(1:nP); y0(1:nP)=y(1:nP); 

tSon=zeros(nT,4^nref);

nPnew=nP;
nTaux=4^nref;
adjnew=sparse(4^(nref+1)*nP,4^(nref+1)*nP);

v=[2 3 1 2];
midpoint=zeros(1,3);

for m=1:nT;
    tnew=t(m,1:3);
    % Perform nrefaux refinements of triangle m
    for r=1:nrefaux
        taux=tnew(r,:);
        xnew=x0(taux);
        ynew=y0(taux);
        % mid points of the sides
        px(1)=(xnew(2)+xnew(3)) /2;
        px(2)=(xnew(1)+xnew(3)) /2;
        px(3)=(xnew(1)+xnew(2)) /2;
        py(1)=(ynew(2)+ynew(3)) /2;
        py(2)=(ynew(1)+ynew(3)) /2;
        py(3)=(ynew(1)+ynew(2)) /2;
        % Check if the mid points have been added to the new grid
        for l=1:3 
            if adjnew(taux(v(l)),taux(v(l+1)))==0
                % Mid point have to be added
                midpoint(l)=nPnew+1;
                nPnew=nPnew+1;
                % Write down the new point
                adjnew(taux(v(l)),taux(v(l+1)))=nPnew;
                adjnew(taux(v(l+1)),taux(v(l)))=nPnew;
                x0(nPnew)=px(l);y0(nPnew)=py(l);
            else
                % Mid point have been added in the new grid
                % use the number
                midpoint(l)=adjnew(taux(v(l)),taux(v(l+1)));
            end
        end
        % Add the new four (son) triangles from (father) triangle r
        tnew=[tnew; 
            taux(1,1) midpoint(3)  midpoint(2);... 
            midpoint(3)  taux(1,2) midpoint(1);...            
            midpoint(2)  midpoint(1)  taux(1,3);...  
            midpoint(3)  midpoint(1)  midpoint(2)]; 
    end
    % add the new triangles to the refined grid
    % We use only the triangles result of performing the last
    % refinement
    t0((m-1)*nTaux+1:m*nTaux,1:3)=tnew(end-nTaux+1:end,:);
    
    % Write down that triangle m is the father of the triangles
    % ((m-1)*nTaux+1):(m*nTaux) of the refined grid t0
    tSon(m,1:nTaux)=((m-1)*nTaux+1):(m*nTaux); 
end

x0=x0(1:nPnew);y0=y0(1:nPnew);
z0=x0*0;
zx0=z0; zy0=z0; zxx0=z0; zxy0=z0; zyy0=z0;

% Evaluate the Argyris element on the new grid
for m=1:nT
   % Read the points of the 
   aux=t0(tSon(m,:),:); aux=unique(aux(:));
   coef=s(t(m,:));
   
   % orientation of the edges: needed for the
   % normal derivative associated functions. 
   or=[t(m,1)<t(m,2) t(m,1)< t(m,3) t(m,2)<t(m,3)];
   or=2*or-1; % Now 1 if positive, -1 if negative
   coef(19:21)=coef(19:21).*or;  % 
   
   % Evaluation of the Argyris element
   px=x0(aux); py=y0(aux);
   cx=x(t(m,1:3));cy=y(t(m,1:3));
   auxz=evalArgy(px,py,[cx(:).'; cy(:).']);
   z0(aux)=coef*auxz;
   
   % Evaluation of the first derivatives
   [auxzx,auxzy]=evalGradArgy(px,py,[cx(:).'; cy(:).']);
   zx0(aux)=coef*auxzx;
   zy0(aux)=coef*auxzy;
   
   
   % Evaluation of the second derivatives
   [auxzxx,auxzxy,auxzyy]=evalHessArgy(px,py,[cx(:).'; cy(:).']);
   zxx0(aux)=coef*auxzxx;
   zxy0(aux)=coef*auxzxy;
   zyy0(aux)=coef*auxzyy;
   
end

% if required, figures are plotted
if drawQ=='Y'
    figure(1)
    trisurf(t0,x0,y0,z0,'facecolor','interp')
    title('function')
    figure(2)
    trisurf(t0,x0,y0,zx0,'facecolor','interp')
    title('x-derivative')
    figure(3)
    trisurf(t0,x0,y0,zy0,'facecolor','interp')
    title('y-derivative')
    figure(4)
    trisurf(t0,x0,y0,zxx0,'facecolor','interp')
    title('xx-derivative')
    figure(5)
    trisurf(t0,x0,y0,zxy0,'facecolor','interp')
    title('xy-derivative')
    figure(6)
    trisurf(t0,x0,y0,zyy0,'facecolor','interp')
    title('yy-derivative')
end
return

