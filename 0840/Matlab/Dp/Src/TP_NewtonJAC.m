function xandw=TP_NewtonJac(xandw0,c,Res);
% input:   xandw=vector of (N-1) unknowns, first guess for start
%             of march. This vector stores both grid points and quadrature weights
%         Res = string variable with name of residual function, not
%         including the termination .m.
%                                   (Res must return a vector of
%         length N-1, but as a function of the arguments (xandw,lambda)
%         where xandw is a ROW vector  and c is a scalar

% rhist(1:NewtMax) stores the norms of the residuals at each iteration.
% The numbers in rhist must decrease rapidly if the iteration is converging.
%        rhist is graphed at each Newton iteration so that one can follow the 
%   convergence in ``real time''.
%
% Output: xandw stores the (N-1) unknowns. By symmetry, and using 
%        known values (two of the grid points are always 1 and -1, and x=0
%         is also always a grid point if N is odd), one can deduce the total 
%         of N grid points and N quadrature weights from xandw. This 
%         expansion to N grid points and N weights is done in the 
%         calling program
global exact_integral  psicoeffsymm;

Nmin1=length(xandw0)

N=Nmin1 + 1;    % N is the number of grid pts;  Nmin1 is the number of unknowns

if mod(N,2)==0
	% npts is even
M= N/2 - 1
else
	% npts is ODD
M= (N-3)/2;
end % if

% Parameters controlling the Newton iteration
%        nalpha= number of stages in powers-of-two minimum
%                 residual
%        NewtMax is the maximum number of Newton iterations to try.
%        JacMod: Jacobian matrix is recomputed every JacMod steps
%                JacMod=1 for classical Newton method
%                JacMod > NewtMax computes the Jacobian just once

NewtMax=20;     % maximum number of Newton iterations
nalpha=4;          % smallest value of relaxation parameter gamma tested is
                             %  2**(nalpha)
JacMod=1;         % Jacmod==1 means that the Jacobian is updated at each
                            % iteration

Jacobian=zeros(Nmin1,Nmin1);

xandw=xandw0;

rhist=0;
% Compute Jacobian matrix;

deltanorm=1.E50;          
haha=1.E-5*max(abs(xandw));    if haha==0,   haha=1.E-5,   end % if
anorm=1.E50;   deltanorm=1.E50;  

for inewt=1:NewtMax    % $$$$$$$$$$$$$$$ BEGIN NEWTONS
if deltanorm > 1.E-10*anorm
residual=feval(Res,xandw,c);

%  if rem(inewt-1,JacMod)==0
	
	
	
% Jacobiancheck   .... compute by finite differences in unknowns
%   for jcol=1:Nmin1
%	  atemp=xandw;  % atemp=n-dim column vector
%      atemp(jcol)=xandw(jcol)+haha;
%      Jaccol=(feval(Res,atemp,c)- residual) / haha;
%      Jacobiancheck(1:Nmin1,jcol)=Jaccol;
%   end

	
	
	if mod(N,2)==1
xunk=xandw(1:M);   wunk=xandw(M+1:2*M+2);
else
   M = (N-2)/2;   % N is even
   xunk=xandw(1:M);     wunk=xandw(M+1:2*M+1);
end % if
% M is the number of interior grid points on one-half 
% of the interval, always excluding x=1 and x=0.

% First step: compute grid point values of the Legendre polynomials.
[NEIG,psicoeffcoll]= size(psicoeffsymm);
           if mod(N,2)==1
	             xunkfull = [0   xunk 1];   % 
           else  % N even
	  %   disp('N is even !!!!!!!!!!!!!!!')
                     xunkfull=[xunk   1]  ;     % xunkfull includes the fixed grid pt., x=1
           end % if
  
PN2dquad=PN_Legendre_vectN(xunkfull,2*NEIG);  % evaluate normalized Legendre
                      % at points xunkfull;  P+N2d(ii,j)=P_{j-1}(xunkfull(ii)).
								   
PN2dquadsymm=PN2dquad(:,1:2:2*NEIG);   % grid point values of symetric functions only 
% compute grid point values of the prolate functions at the quadrature points

psi_j = PN2dquadsymm * psicoeffsymm;  %psi_j =   psi( x_i, j )
	
	PNx2dquad=PNx_Legendre_vectN(xunkfull,2*NEIG);  % evaluate normalized Legendre
         % at points xunkfull;  PNx2d(ii,j)=P_{j-1}(xunkfull(ii)).
								   
PNx2dquadsymm=PNx2dquad(:,1:2:2*NEIG);   % grid point values of symetric functions only 
% compute grid point values of the prolate functions at the quadrature points

psix_j = PNx2dquadsymm * psicoeffsymm;  %psix_j =   d psi/dx( x_i, j )


iparity=rem(N,2);
	
for jcol=1:M
	Jacobian(:,jcol) = - wunk(jcol+iparity) * psix_j(jcol+iparity,:)';
end % jcol
for jcol=M+1:Nmin1
	Jacobian(:,jcol) = - psi_j(jcol-M,:)';
end % jcol	  
	  Jacobian = 2*Jacobian; % symmetry correction
if mod(N,2) == 1
	% N odd, avoid double-counting the weight for x=0
	Jacobian(:,M+1)=0.5 * Jacobian(:,M+1);
end % if
	  

% end % if to compute Jacobian

%disp('Jacobian')
%Jacobian
%Jacobiancheck
% Jacdiff=Jacobian - Jacobiancheck


rnorma(inewt)=max(abs(residual));
     rnorm=rnorma(inewt);   anorm=max(abs(xandw));
	 
delta_a=Jacobian\residual;  % delta_a is n-dim column vector

dnorma(inewt)=max(abs(delta_a));
deltanorm=dnorma(inewt);

minrnorm=1.E50;
    for ialph=0:nalpha
    % alphaline=(1+0.005*i)/(2^ialph);
           alphaline= 1 / (2^ialph);
    atemp1=alphaline* delta_a';  % n-dimensional column vectors
    atemp= xandw - atemp1;
	
    residual=feval(Res,atemp,c)';
    rnorma(ialph+1)=max(abs(residual)); 
              if rnorma(ialph+1) < minrnorm
              iminnorm=ialph;
              underrel=alphaline;
              minrnorm=rnorma(ialph+1);
              end % if
    end  % ialpha loop

disp([' inewt and stuff','  c=',num2str(c)])
inewt, minrnorm, underrel
	
rhist(inewt)=minrnorm;    gammaa(inewt)=underrel;
rhist,gammaa

atemp1=underrel* delta_a'; % atemp1 is n-dim column vector
xandw=xandw - atemp1;
iita(inewt)=inewt;

Zinewt_underrel_minrnorm=[inewt underrel minrnorm];
semilogy(1:inewt,rhist,'g-')
set(gca,'FontSize',18)
title(['inewt=',int2str(inewt),' c=',num2str(c)])
ylabel('Residual')
xlabel('Iteration number')

drawnow
                               end % Convergence test if block
             end   % of Newton loop, inewt
