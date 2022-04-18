%
%
%DISODE45  is a program for the numerical integration of ODEs with
%   DISCONTINUITIES.  Filippov type systems are allowed.
%   It is based on the DOPRI5(4) pair of embedded Runge-Kutta methods.
%   The function or functions that define the manifolds of the
%   discontinuities g(x,y) are supposed to be known, provided
%   by the user.
%
%[TOUT,YOUT] = DISODE45(ODEFUN,SWITCHFUN, TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates 
%     the system of differential equations y' = f(t,y) from time T0 to TFINAL 
%     with initial conditions Y0. 
%  Input arguments:
%    ODEFUN is a function handle. For a scalar T and a vector Y, ODEFUN(T,Y) 
%         must return a column vector corresponding to f(t,y).      
%    SWITCHFUN is a function handle. [values, isterminal, direction]=SWITCHFUN(T,Y) must return
%         a column vector VALUES in which the component i contains the value of the 
%         i function g_i(t,y) defining the i-eme discontinuity manifold.
%         It also returns a column vector ISTERMINAL and a column vector
%         DIRECTION.
%  Output arguments:
%     Each row in the solution array YOUT corresponds to a time 
%     returned in the column vector TOUT.  
%                

function [xx,yy,tdis,ydis,idis,stats]=disode45(FUN,switchfun, tspan,Y,options)  %Function disode45
% ----------------------------------------------------------------------
%       disode45 integrates a system of discontinuous IVP
% ----------------------------------------------------------------------
%
%  TO DO!!!!!
%    Permitir tspan un vector con varios puntos de salida
%    Arreglar los datos de estadisticas
%    Arreglar mas opciones tipo ode45
%    Arreglar la salida tipo ode45  (estructura)
%    Calculo automatico del paso inicial tras una discontinuidad
%    Poner un control de maximo numero de pasos o de paso minimo
%    Considerar que en el punto inicial haya varias
%    discontinuidades
%
%  
%  Setting the default working parameters
%
EABS=1.e-4;
EREL=1.e-6;
gradswitchfun=[];
doatswitch=@dummy;
exactgradient=0;
    % 0 means  No external function for the gradient of switchfun is provided
    % 1 means  External function for the gradient of switchfun is provided
H0=0;
eventcontrol=0;
Refine=4;
    % 0 means  minimum level of control of possible discontinuity
if nargin > 4
   if size(options.AbsTol) ~=0;
       EABS=options.AbsTol;
   end
   if size(options.RelTol) ~=0;
       EREL=options.RelTol;
   end
   if size(options.Gradient) ~=0;
       gradswitchfun=options.Gradient;
       exactgradient=1; 
   end
   if size(options.ActionSwitch) ~=0;
       doatswitch=options.ActionSwitch; 
   end
   if size(options.InitialStep) ~=0;
       H0=options.InitialStep;
   end
   if size(options.EventControl) ~=0;
       eventcontrol=options.EventControl;
   end
   if size(options.Refine) ~=0;
       Refine=options.Refine;
   end
end
%Persistents variables
  persistent IFIR;
%Initialization of variables
  % NMAX = 10000;
  IFIR = 1;
        %
        %  stats(1)= Naccpt
        %  stats(2)= Nrejct
        %  stats(3)= Nslideaccpt
        %  stats(4)= Nslidereject
        %  stats(5)=Nstepdis
        %  stats(6)=Number of transversal discontinuities
        %  stats(7)=Number of sliding discontinuities
        %  stats(8)=Number of exits of a sliding region
        %  stats(9)=Number of evaluations of the vector field
        %  stats(10)=Number of evaluations of the manifold function
        %  stats(11)=Number of evaluations of the  gradient function
        %
if (IFIR==1),  % It is the first time the function is called
    IFIR=0;
    disctype=0; %  No discontinuity at this point
    X=tspan(1);
    XEND=tspan(2);
    xx=X;
    yy=Y';
    tdis=[];
    ydis=[];
    idis=[];
    stats=zeros(1,11);
%    [g0,~,direction]=switchfun(X,Y);
    f0 = FUN(X, Y);
    stats(9)=stats(9)+1;
%    stats(10)=stats(10)+1;
    TOL= EABS + EREL* max(abs(Y));
%
%   Checking if the initial point is a discontinuity point
%
    tangentside=1;
    inddis=0;
    [xout,yout,ff,disctype,idis,endslid,stats]=classifypoint(X,Y,X,Y,X,Y,inddis,FUN,switchfun,idis,tangentside,stats,TOL,exactgradient);
    format=' X= %26.20e, Y= (%27.20e';
    for ii=1:size(Y,1)-1
        format=[format ', %27.20e'];
    end
    format=[format ')'];
    if disctype==3 && endslid <0,
%        fprintf(['\n Filippov point  ' format], X, Y);
%        stats(7)=stats(7)+1;
        tdis=[tdis X];
        ydis=[ydis; Y'];
    elseif disctype==3 && endslid >=0,
%        fprintf(['\n Tangent point   '  format], X, Y);
%        stats(7)=stats(7)+1;
        tdis=[tdis X];
        ydis=[ydis; Y'];
    elseif disctype==1,
%        fprintf(['\n Transversal discontinuity, exit at ' format], X, Y);
%        stats(6)=stats(6)+1;
        tdis=[tdis X];
        ydis=[ydis; Y'];
        disctype=0;
    end
    X=xout;
    Y=yout;
    WRK(:,1)=ff;
end
if H0 ==0  %H0= was not set at the user options
   H0= min((TOL/max(1.e-15,norm(WRK(:,1))))^(1./5.),(XEND-X)/10);
end
% ----------------------------------------------------------------------
%      Loop controled by variable disctype
%         disctype=0,  normal step
%         disctype=1,  transversal discontinuity
%         disctype=3,  sliding discontinuity
%         disctype=4,  exit of a sliding region
%         disctype=5,  discontinuity inside a sliding region
%         disctype < 0,  XEND has been reached
% ----------------------------------------------------------------------
    H=H0;
  while disctype >=0,   % main loop that distributes the flow according to the kind of point
    if disctype==0,
      [WRK,xx,yy,xout,yout,HOLD,XOLD,YOLD,WRKOLD,disctype,xdisaprox,stats]=normalintegration(FUN,switchfun,H,X,WRK,Y,EABS,EREL,xx,yy,XEND,stats,Refine,eventcontrol); 
      X=xout;
      Y=yout;
    elseif disctype==1, % Discontinuity detected, it must be accurately located
        tol= EABS + EREL* max(abs(Y));
        [WRKout,xout,yout,H,disctype,endslid,tangentside,tdis,ydis,idis,stats]=FindDisc(FUN,switchfun,HOLD,X,XOLD,WRKOLD,Y,YOLD,xdisaprox,tdis,ydis,idis,stats,exactgradient,tol,gradswitchfun,doatswitch);
        WRK=WRKout;
        X=xout;
        Y=yout;
        xx=[xx X];
        yy=[yy ; Y'];
        if H==HOLD,
           H=H0;  %  this could be improved
        end
    elseif disctype==3,  %Filipov discontinuity; entering a sliding region
        [WRKout,xout,yout,HOLD,XOLD,YOLD,WRKOLD,disctype,xdisaprox,xx,yy,stats]=slide(FUN,switchfun,H,X,WRK,Y,EABS,EREL,xx,yy,XEND,idis,endslid,tangentside,stats,Refine,exactgradient,gradswitchfun); 
        X=xout;
        Y=yout;
        WRK = WRKout;
    elseif disctype==5,  % Discontinuity inside a sliding region detected; it must be accurately located
        tol= EABS + EREL* max(abs(Y));
        [WRKout,xout,yout,H,disctype,tdis,ydis,idis,stats]=FindDiscpro(FUN,switchfun,HOLD,X,XOLD,WRKOLD,Y,YOLD,xdisaprox,tdis,ydis,idis,endslid,tangentside,stats,exactgradient,tol,gradswitchfun,doatswitch);
        WRK=WRKout;
        X=xout;
        Y=yout;
        xx=[xx X];
        yy=[yy; Y'];
        H=H0/4;
    end % end of if disctype='value'
 end %end of  WHILE disctype >= 0
 if xout==tdis(end),
    fprintf('\n Integration ended at the switching point \n');
 end
end %Fin de function INTEGRA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Funcion normalintegration
function [WRKout,xx,yy,xout,yout,HOLD,XOLD,YOLD,WRKOLD,disctype,xdis,stats]=normalintegration(FUN,switchfun,H,X,WRK,Y,EABS,EREL,xx,yy,XEND,stats,Refine,eventcontrol)
% -----------------------------------------------------------------------
% This function integrates the problem along a non sliding region
% 
%       The integration follows until
%               a discontinuity point is detected; in this case it returns
%               with:
%                   A rough approximation of the discontinuity xdis
%                   A last step [x, x+h] such that  xdis-x-h < 0.15*h
%               X reaches XEND or
%               a too short stepsize is needed or a maximum number of
%                   steps is reached.
% ---------------------------------------------------------------------- 
%
% Input arguments
%    FUN, switchfun   function handlers for vercot field and switching surfaces
%    H
%    X, Y,  initial point to start this part of the integration
%    XEND    end point of the integration interval
%    WRK    vector containing the data of the last accepted step  
%    EABS, EREL  absolute and relative error bounds
%    xx, yy   vectors containing the integration points
%    stats    vector containing the statistics of the integration
%    eventcontrol    variable to state the control to see if there is a
%         discontinuity
%
% Output arguments
%    WRKout    vector containing the data of the last accepted step
%    xx, yy    vectors containing the integration points
%    xout,  yout   integration (successful) point
%    HOLD   las succesful step size
%    XOLD, YOLD,  penultimate successful integration point
%    WRKOLD   vector containing the data of the penultimate accepted step
%    disctype   type of discontinuity point if detected (-1 otherwhise)
%    xdis   approximate point where a discontinuity has been located
%    stats  vector containing the statistics of the integration 
%
  advance = true; % set to false when XEND is attained or a discontinuity is detected
  EJECUTAR = true;
  REJECT  = false;
  disctype=0;
  UROUND=eps(1.0);
  xdis=XEND;
  xout=X;
  checkall=eventcontrol;
  discdetected=0;
  HOLD=H;
  XOLD=X;
  YOLD=Y;
  WRKOLD=WRK;
  [gxy,isterminal,direction]=switchfun(X,Y);
   stats(10)=stats(10)+1;
  while (advance), 
%        IF( NSTEP.GT.NMAX)  THEN
%           IERR = -5
%           RETURN
%        END IF
    if (EJECUTAR), % Executed only at accepted steps  
      URO1 = 5.0D0 * UROUND * max(1.0D0, abs(X));
      if ((X-XEND) + URO1 > 0.0D0),
        X = XEND;
        yout=Y;
        xout=X;
        disctype=-1;
        return;  % End of the integration interval reached
      end
      if ((X + H - XEND) > 0.00D0),
        H = XEND - X;
      end  
    end  
    if H < 5*eps(X)
         fprintf('\n Minimum step size  h=%g attained at X= %g \n  Integration stopped \n', H, X);
           yout=Y;
           XOLD=X;
           YOLD=Y;
           WRKOLD=WRK;
           WRKout=WRK;
         disctype=-1;
         return;
    end
    [XPH,Y1,WRKout,ERR,disctype,gxy,xdis,stats]=RK(FUN,switchfun,X,H,Y,WRK, checkall,gxy,stats,EABS,EREL); 
    WRK=WRKout; 
    TOL= EABS + EREL* max(abs(Y1));
    if disctype==0, % There is no discontinuity
       if (ERR <= TOL),
% ----------------------------------------------------------------------
%       Accepted step
%      fprintf('\nPASO ACEPTADO %g   [%g, %g]', H, X, XPH);
% ----------------------------------------------------------------------
           stats(1)=stats(1)+1;  %  Accepted steps counter Naccpt
           [XOLD, HOLD, YOLD, WRKOLD]=STORE(X, H, Y, WRK);
           WRKout=WRK; 
           if Refine>=2
              for irefine=1:Refine-1
                 xrefine=X+H*irefine/Refine;
                 [WRK, Y2]=ESTIRA(FUN,XOLD,YOLD,WRKOLD,HOLD,xrefine,-1); 
                 xx=[xx xrefine];
                 yy=[yy ; Y2'];
              end
           end
% ----------------------------------------------------------------------
%       Update the numerical approximation, the current point, the
%       first evaluation of the RK method at the next step and the
%       next stepsize to be used.
% ----------------------------------------------------------------------
           WRK(:,1) = WRK(:,7);
           Y    = Y1;
           X   = XPH;
           xx=[xx X];
           yy=[yy ; Y'];
           FAC = min(0.9*(TOL/(ERR+1D-17))^(1./5.), 2.0D0);
           if (REJECT), 
             FAC = min(FAC, 1.0D0);
           end               
           H       = FAC*H;
           REJECT  = false;
           EJECUTAR = true;
           checkall=eventcontrol;
           if discdetected==1,  % discon detectada en el paso anterior
              checkall=eventcontrol+1;
              discdetected=0;
           end
      else
% ----------------------------------------------------------------------
%       Rejected step
%        fprintf('\nPASO rechazado %g   [%30.20e, %g]', H, X, X+H);
% ----------------------------------------------------------------------
         FAC = max(0.9*(TOL/(ERR+1D-12))^(1./5.), 0.10D0); 
         REJECT  = true;
         H       = FAC*H;
         if disctype~=0,
             H=min(H,0.5*(xdis-X));
             stats(5)=stats(5)+1;  % Discont. detected counter
         else
             stats(2)=stats(2)+1;  % Normal rejected steps counter Nrejct
         end
         EJECUTAR=false;
       end   
    else  % Discontinuity detected at this step.  
        stats(5)=stats(5)+1;   % Discont. detected counter
%   fprintf('\nPASO rechazado discon %g   [%30.20e, %g] %g', H, X, X+H, xdis);
        if xdis-X >= 0.15*HOLD ||  xout==X,  %Estira can not be directly applied. Repeating a step  
           if xout==X  % && xdis-X < H/4,
               H=0.9*(xdis-X)/30 ;
           else
               H=0.9*(xdis-X);
           end
           REJECT  = true;
           checkall=eventcontrol+1;
           discdetected=1;
           if xdis-X < min(TOL,sqrt(eps(1))/10),  % Discontinuity will be located with Euler method
             disctype=1;
             WRKOLD(:,7)=WRK(:,1);
             HOLD=H;
             XOLD=X;
             YOLD=Y;
             advance=false;
             xout=X;
             yout=Y; 
           end
        else  %Estira can be directly applied
           disctype=1;
           advance=false;
           xout=X;
           yout=Y;
        end
    end
  end
  WRKout=WRK;
end % End of function normalintegration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WRKout,xout,yout,H,disctype,endslid,tangentside,tdis,ydis,idis,stats]=FindDisc(FUN,switchfun,HOLD,X,XOLD,WRKOLD,Y,YOLD,xdisaprox,tdis,ydis,idis,stats,exactgradient,tol,gradswitchfun,doatswitch)
% -----------------------------------------------------------------------
% The discontinuity has been detected in [X, 2*xdis-X]
% Taking the step from XOLD to X with step-size HOLD,
% the discotinuity point is located using the function estira
% -----------------------------------------------------------------------
%  (xplus, yplus) is a point just after the discontinuity point
%  (xminus,yminus) is a point just before the discontinuity point
%  they are used to evaluate the Filipov vector field in case it is needed
% 
%
% Input arguments
%    FUN, switchfun   function handlers for vector field and switching surfaces
%    H
%    X, Y,  initial point to start this part of the integration
%    XOLD, YOLD,  penultimate successful integration point
%    XEND    end point of the integration interval
%    WRK    vector containing the data of the last accepted step  
%    EABS, EREL  absolute and relative error bounds
%    xx, yy   vectors containing the integration points
%    stats    vector containing the statistics of the integration
%    eventcontrol    variable to state the control to see if there is a
%         discontinuity
%    gradswitchfun  Function handler, computes the gradient of switchfun
%    doatswitch    Function handler, to be called at switching points for
%         which isterminal=-1 has been set
%
% Output arguments
%    WRKout    vector containing the data of the last accepted step
%    xx, yy    vectors containing the integration points
%    xout,  yout   integration (successful) point
%    HOLD   las succesful step size
%    XOLD, YOLD,  penultimate successful integration point
%    WRKOLD   vector containing the data of the penultimate accepted step
%    disctype   type of discontinuity point
%    tdis, ydis    point where a discontinuity has been located
%    idis   hypersurface where the discontinuity is located  
%    stats  vector containing the statistics of the integration 
%
%    xout=X;

    format=' X= %26.20e, Y= (%27.20e';
    for ii=1:size(Y,1)-1
       format=[format ', %27.20e'];
    end
    format=[format ')'];
    
    tangentside=1;
    
     if XOLD==X,
  %  if xdisaprox-X < min(sqrt(eps(1))/100,tol/10),
        euler=1;
        XMAX=max(2*xdisaprox-X,X+min(sqrt(eps(1))/10,100*tol));
    else
        euler=0;
        XMAX=max(2*xdisaprox-X,X+0.32*HOLD);
     end
   
%
%   Lower extreme of the interval
%
    X0=X;
    ynew=Y;
    [gg0,~,~]=switchfun(X0, ynew);
    stats(10)=stats(10)+1;
    value0=min(abs(gg0(gg0~=0)));
    if any(abs(gg0)==0),
       xout=X+min(1.e-10, max([tol/100, eps(X), 100*eps(1)]));
       yout=Y+min(1.e-10, max([tol/100, eps(X), 100*eps(1)]))*WRKOLD(:,7);
       WRKout(:,1)=FUN(xout,yout);
       stats(9)=stats(9)+1;
       disctype=0;
       H=HOLD/4;
       endslid=1;
       return
    end
% 
%   Upper extreme of the interval
%
    X1=XMAX;
    if euler,
         ynew=Y+(XMAX-X)*WRKOLD(:,7);
    else
       [WRKout, ynew]=ESTIRA(FUN,XOLD,YOLD,WRKOLD,HOLD,X1, 1);
       WRKOLD=WRKout; % Include in WRKOLD the two additional stages of the continuous
       stats(9)=stats(9)+2;
    end
    [gg1,~,direction]=switchfun(X1, ynew);
    stats(10)=stats(10)+1;
    gcondp=(gg1.*gg0<=0).*(gg1.*direction>=0);
    if any(gcondp)  % Checking if there is a discontinuity in the interval
         value1=-max(abs(gg1(gcondp==1)));
         ind=find(gcondp==1);
    else
%         fprintf('\n Discon fallida !!!! %g %g %g ',X0,X1,euler);
         xout=X;
         yout=Y;
         disctype=0;
         H=HOLD/4;
         endslid=1.0;
         WRKout(:,1)=WRKOLD(:,7);
%          if euler==0,
%             euler=1;
%             X1=max(2*xdisaprox-X,X+min(sqrt(eps(1))/10,100*tol)); 
%             ynew=Y+(XMAX-X)*WRKOLD(:,7);
%             [gg1,~,direction]=switchfun(X1, ynew);
%             stats(10)=stats(10)+1;
%             gcondp=(gg1.*gg0<=0).*(gg1.*direction>=0);
%             if any(gcondp)  % Checking if there is a discontinuity in the interval
%                 value1=-max(abs(gg1(gcondp==1)));
%                 ind=find(gcondp==1);
%             else
% %                fprintf('\n Discon fallida !!!! %g %g %g ',X0,X1,euler);
%                 return
%             end
%          else
             return
%          end
    end  
%
%   Modified Secant Method's loop
% 
    xx1=X0;  % Lower limit of the interval containing the discontinuity
    xx2=X1;  % Upper limit of the interval containing the discontinuity
    value=value0;
    ii=1;
  %  for II=1:60,
    while  abs(X1-X0) >= 10*max(eps(X1),eps(1)) && abs(value)>= 1.e-25  && ii< 60,
       xnew=X1-value1*(X1-X0)/(value1-value0);
       if xnew >= xx2 || xnew<= xx1
           xnew=(xx1+xx2)/2;  % If secant method leads out the interval, use bisection
       end
       if euler,
         ynew=Y+(xnew-X)*WRKOLD(:,7);
       else
         [~, ynew]=ESTIRA(FUN,XOLD,YOLD,WRKOLD,HOLD,xnew,  0);
       end
       [gnew,isterminal,direction]=switchfun(xnew, ynew);
       stats(10)=stats(10)+1;
       gauxb =gnew.*direction>=0;
       gcond =(gnew.*gg0<0).*gauxb;
       gcondp=(gnew.*gg0<=0).*(gnew.*direction>=0);

       if any(gcondp)
          value=-max(abs(gnew(gcondp==1))); % current point is strictly after the discontinuity
          ind=find(gcondp==1);
        else
          value=min(abs(gnew(ind))); % current point is before or at the discontinuity
       end
%
%  Updating the interval where the discontinuity is located
%
       if value <=0
          xx2=xnew;
       else
          xx1=xnew;
       end
%
%  Updating the two last points for the secant method
%
       X0=X1;
       value0=value1;
       X1=xnew;
       value1=value;
       ii=ii+1;
%
%  If convergence is attained, the loop is finished
%
    end %end of while
    if ii>=60,
       ii
       fprintf('\n Secant method attained the maximum of iterations');
    end
   
%     fprintf('\n Secant method  iter=%g %30.20e %30.20e %30.20e',II, X1-X0, gnew);
%     fprintf('\n  discontinuidad en X=%30.20e, Y= (%g, %g %g)', xnew, ynew );
     tdis=[tdis xnew];
     ydis=[ydis; ynew'];
     inddis=ind(1);
%
%  Numerical discontinuity found at (xnew, ynew)
%
%
%  Once the discontinuity point is computed, two points, very close
%  but strictly at the each side of the manifold are computed
%  They are used to evaluate the vector field at both sides
%  and determine if the discontinuity is transversal or Filippov
%
   xminus=xnew;
   yminus=ynew;
   xplus=xnew;
   yplus=ynew;
   gplus=gnew;
   if euler && (xnew-X)^2>tol/100,    
      fprintf('\n Warning!!  Euler used a step xnew-X  %g that can be largs for the tolerance %g', xnew-X, tol);
      X
      euler
   end
   
   if size(ind,1)~=1 || size(ind,2)~=1
       fprintf('\n Warning, multiple dicontinuity  index= !!! %g  %g %g %g  %g %g', ind);
   end
   if any(isterminal(ind)<0) % the switching point requires calling the external function
       xout=xnew;
       yout=doatswitch(xnew,ynew);
       ynew=yout;
   end
   H=HOLD;
   endslid=1;
   if any(isterminal((gplus.*gg0<=0) & (direction.*gg0<=0))==1)  % the switching point ends the integration
       disctype=-4;
       idis=[idis inddis];
       stats(6)=stats(6)+1;
       return
   end
   if all(abs(isterminal(ind))==2),
 %        fprintf(['\n Discontinuity, exit at ' format], xnew, ynew);
         xout=xnew;
         yout=ynew;
         stats(6)=stats(6)+1;
         disctype=0;
         idis=[idis inddis];
         stats(6)=stats(6)+1;
         WRKout(:,1)=FUN(xout,yout);
         stats(9)=stats(9)+1;
         return
   end
   
   if all(isterminal(ind)>=0) % the switching point requires calling the external function
%
%  Compute a point just after the discontinuity if (xnew,ynew) is not
%
    ii=0;
    if all(gcond==0) % means that the numerical discontinuity point is in fact not after it
          ii=1;
          xplus=xnew;
          while all(gcond==0) && ii<10  % the point is not after the discontinuity  Que pasa si ii>10 ??
             xplus=xplus+10^ii*eps(xnew);
             if euler,
                  yplus=Y+(xplus-X)*WRKOLD(:,1);
             else
                  [~, yplus]=ESTIRA(FUN,XOLD,YOLD,WRKOLD,HOLD,xplus, 0);
             end
             [gplus,~,direction]=switchfun(xplus,yplus);
             gcond=(gg0.*gplus<0).*(gg0.*direction<=0);
             stats(10)=stats(10)+1;
             ii=ii+1;
         end
         if ii>=10
             ii
             ' yplus en riesgo '
         end
    end
%
%  Compute a point just before the discontinuity if (xnew,ynew) is not
%
    if any(gcondp) % means that the numerical discontinuity point is in fact not before it
         ii=1;
         xminus=xnew;
         while any(gcondp) && ii<10 % the point is not before the discontinuity  Que pasa si ii>10 ??
            xminus=xminus-10^ii*eps(xnew);
            if euler,
                yminus=Y+(xminus-X)*WRKOLD(:,1);
            else
               [~, yminus]=ESTIRA(FUN,XOLD,YOLD,WRKOLD,HOLD,xminus, 0);
            end
            [gminus,~,direction]=switchfun(xminus,yminus);
            gcondp=(gg0.*gminus<=0).*(gg0.*direction<=0);
            stats(10)=stats(10)+1;
            ii=ii+1;
         end
         if ii>=10
            ii
            ' yminus en riesgo '
         end
    end 
   end
   

   tangentside=-sign(gg0(inddis));
%     if any(isterminal(ind)==-1) % the switching point requires calling the external function
       [xout,yout,ff,disctype,idis,endslid,stats]=classifypoint(xnew,ynew,xminus,yminus,xplus,yplus,inddis,FUN,switchfun,idis,tangentside,stats,tol,exactgradient);
       if disctype==3 && endslid <0,
%         fprintf(['\n Filippov point  ' format], xnew, ynew);
%         stats(7)=stats(7)+1;
       elseif disctype==3 && endslid >=0,
%         fprintf(['\n Tangent point ' format], xnew, ynew);
%         stats(7)=stats(7)+1;
       elseif disctype==1,
%         fprintf(['\n Transversal discontinuity, exit at ' format], xnew, ynew);
%         stats(6)=stats(6)+1;
         disctype=0;
       else
           idis=[idis inddis];
           stats(6)=stats(6)+1;
       end
       WRKout(:,1)=ff;
       H=HOLD;
       return
end   % End of FinDisc function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion RK
function [XPH,Y1,WRKout,ERR,disctype,gxy,xdis,stats]=RK(FUN,switchfun,X,H,Y,WRK,checkall,gxy,stats,EABS,EREL)                   %Funcion RK
%
%  Advance of one step by the DOPRI5(4) pair
%
%
%  TODO!!!!!   Afinar la aproximacion de la discontinuidad,
%   sobre todo si esta cerca del punto X, por si es una
%   discontinuidad muy cercana a la anterior.
%   Optimizar el numero de evaluaciones de switchfun(x,y)
%
%
% Input arguments
%    FUN, switchfun   function handlers for vercot field and switching surfaces
%    H
%    X, Y,  initial point to start this step
%    WRK    vector containing the data of the last accepted step  
%    EABS, EREL  absolute and relative error bounds
%    checkall  how sign change of switchfun is checked
%    stats    vector containing the statistics of the integration
%
% Output arguments
%    WRKout    vector containing the data of the last accepted step
%    XPH, Y1   next integration points
%    ERR   norm of the local error estimation
%    disctype   type of discontinuity point if located (-1 otherwise)
%    xdis    approximate  point where a discontinuity has been located
%    stats  vector containing the statistics of the integration 
%

  A=[0,    0, 0, 0, 0, 0,0; ...
     1/5., 0,0,0,0,0,0; ...
     3/40., 9/40, 0, 0, 0, 0,0; ...
     44/45., -168/45, 160/45, 0, 0, 0,0; ...
     19372/6561,    -76080/6561,  64448/6561,  -1908/6561, 0, 0,0; ...
     477901/167904, -1806240/167904, 1495424/167904, 46746/167904, -45927/167904, 0,0; ...
     12985/142464,    0,   64000/142464,  92750/142464,  -45927/142464,  18656/142464,0];
  
  C=[0,  1/5, 3/10, 4/5, 8/9. , 1.0, 1.0];
  B=[12985/142464,    0,   64000/142464,  92750/142464,  -45927/142464,  18656/142464, 0];
  B1=[1921409/21369600,  0,  9690880/21369600,  13122270/21369600,  -5802111/21369600,  1902912/21369600, 534240/21369600];
  disctype=0;
  xdis=X+H;
  XPH=X;
  WRKout=WRK;
  for K=2:7
    HH =  H*C(K);
    Y1 = Y + (H*WRK(:,1:K-1)*A(K,1:K-1)');
    if checkall>=1,
    [gyx,~,direction]=switchfun(X+HH, Y1);
    stats(10)=stats(10)+1;
    auxa= gyx.*gxy<=0;
    auxb=gyx.*direction>=0;
     if any(auxa.*auxb),
        disctype=1;
        xdis=X+H*(C(K)+min(C(K-1),8/9.))/2;
        teta0=C(K-1);
        teta1=C(K);
        if K>=4
             for it=1:3
             teta=(teta0+teta1)/2;
             acon=[(18*teta - 75*teta^2 + 100*teta^3)/18;-5*(-9*teta^2 + 20*teta^3)/6;10*(-3*teta^2 + 10*teta^3)/9];
             Y1=Y+H*WRK(:,1:3)*acon;
             [estg,isterminal,direction]=switchfun(X+H*(teta0+teta1)/2, Y1);
             stats(10)=stats(10)+1;
             if (any((gxy.*estg<=0))),
                 teta1=teta;
             else
                 teta0=teta;
             end
             end
             xdis=X+H*(teta0+teta1)/2;
        end
        ERR=2.0;
        return
      end
    end
    WRK(:,K) = FUN(X + HH , Y1);
    stats(9)=stats(9)+1;
  end  %end of FOR
  EST= 0.0D0;
  Y1 = Y;
  
  for K=1:7,
      EST = EST+H*WRK(:,K)'*(B(K)-B1(K));
      Y1 = Y1 + (H*WRK(:,K)'*B(K))';
  end
  ERR=max(abs(EST));
  TOL= EABS + EREL* max(abs(Y1));
  stats(10)=stats(10)+1;
  [gyx,isterminal,direction]=switchfun(X+H, Y1);
  gnew=gyx;
  if any((gnew.*gxy<=0).*(gyx.*direction>=0))
      xdis=X+H;
      disctype=1;
      return
  end
  if checkall==0,
     auxa= gyx.*gxy<0;
     auxb=gyx.*direction>=0;
    if ERR > TOL || any(auxa.*auxb),
       for K=2:7
         Y2 = Y + (H*WRK(:,1:K-1)*A(K,1:K-1)');
         [gyx,isterminal,direction]=switchfun(X+H*C(K), Y2);
         stats(10)=stats(10)+1;
         auxa= gyx.*gxy<0;
         auxb=gyx.*direction>=0;
         if any(auxa.*auxb),
            xdis=X+H*(C(K)+min(C(K-1),8/9.))/2;
            teta0=min(C(K-1),8/9.);
            teta1=C(K);
            if K==3
               for it=1:1
                  teta=(teta0+teta1)/2;
                  acon=[(2*teta-5*teta^2)/2;5*teta^2/2];
                  Y2=Y+H*WRK(:,1:2)*acon;
                  [estg,isterminal,direction]=switchfun(X+H*(teta0+teta1)/2, Y2);
                  stats(10)=stats(10)+1;
                  if (any((gxy.*estg<=0))),
                     teta1=teta;
                  else
                     teta0=teta;
                  end
               end
               xdis=X+H*(teta0+teta1)/2;
            elseif K>=4
               for it=1:3
                  teta=(teta0+teta1)/2;
                  acon=[(18*teta - 75*teta^2 + 100*teta^3)/18;-5*(-9*teta^2 + 20*teta^3)/6;10*(-3*teta^2 + 10*teta^3)/9];
                  Y2=Y+H*WRK(:,1:3)*acon;
                  [estg,isterminal,direction]=switchfun(X+H*(teta0+teta1)/2, Y2);
                  stats(10)=stats(10)+1;
                  if (any((gxy.*estg<=0))),
                     teta1=teta;
                  else
                     teta0=teta;
                  end
               end
               xdis=X+H*(teta0+teta1)/2;
            end
            disctype=1;
            return
         end
       end
    end
  end
  
  if checkall>=2 && ERR <= TOL
      for icheck=1:6*checkall
          xcheck=X+H*icheck/(6*checkall+1);
          [WRK, Y2]=ESTIRA(FUN,X,Y,WRK,H,xcheck,-1);
          [gyx,~,direction]=switchfun(xcheck, Y2);
          stats(10)=stats(10)+1;
          auxa= gyx.*gxy<0;
          auxb=gyx.*direction>=0;
          if any(auxa.*auxb),
              disctype=1;
              xdis=xcheck;
              ERR=2.0;
              return
          end
      end
  end
  if ERR <= TOL
     gxy=gnew;
  end
  XPH=X+H;
  WRKout=WRK; %Variable de entrada/salida
end % End of funciton RK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Funcion rkpro
function [XPH,Y1,WRKout,ERR,disctype,xdis,stats]=rkpro(FUN,switchfun,idis,X,H,Y,WRK,stats,exactgradient,tol,endslid0,tangentside,gradswitchfun)                   %Funcion rkpro
% -----------------------------------------------------------------------
% This function advances one step with the DOPRI 5(4) pair
% from X to X+H but  projecting every stage onto the 
% sliding manifold
% -----------------------------------------------------------------------
%
%  TODO!!!!   Hay que mejorarla, ajustando la dirección de projección 
%  w para un caso general.  Ver porque funciona el gradiente peor que [0;1]
%  Controlar posibles discon transversales durante la zona de sliding
%
% Input arguments
%    FUN, switchfun   function handlers for vercot field and switching surfaces
%    H
%    X, Y,  initial point to start this step
%    WRK    vector containing the data of the last accepted step  
%    EABS, EREL  absolute and relative error bounds
%    checkall  how sign change of switchfun is checked
%    endslid0   Indicator of type of sliding
%              endslid < 0   pure sliding
%              endslid >=0  tangent sliding
%    stats    vector containing the statistics of the integration
%
% Output arguments
%    WRKout    vector containing the data of the last accepted step
%    XPH, Y1   next integration points
%    ERR   norm of the local error estimation
%    disctype   type of discontinuity point if located (-1 otherwise)
%    xdis    approximate  point where a discontinuity has been located
%    stats  vector containing the statistics of the integration 
%
  A=[0,    0, 0, 0, 0, 0,0; ...
     1/5., 0,0,0,0,0,0; ...
     3/40., 9/40, 0, 0, 0, 0,0; ...
     44/45., -168/45, 160/45, 0, 0, 0,0; ...
     19372/6561,    -76080/6561,  64448/6561,  -1908/6561, 0, 0,0; ...
     477901/167904, -1806240/167904, 1495424/167904, 46746/167904, -45927/167904, 0,0; ...
     12985/142464,    0,   64000/142464,  92750/142464,  -45927/142464,  18656/142464,0];
  
  C=[0,  1/5, 3/10, 4/5, 8/9. , 1.0, 1.0];
  B=[12985/142464,    0,   64000/142464,  92750/142464,  -45927/142464,  18656/142464, 0];
  B1=[1921409/21369600,  0,  9690880/21369600,  13122270/21369600,  -5802111/21369600,  1902912/21369600, 534240/21369600];
  N=size(Y,1);
%
  XPH=X+H;
  WRKout=WRK;
  ERR=2.0;
  disctype=3;
  xdis=X;
  ll=0;
 
  if exactgradient
      w=feval(gradswitchfun,X,Y,idis);
      w=w/norm(w);
      stats(11)=stats(11)+ 1;
  else
      w=graddif(switchfun,X,Y,idis,tol);
      stats(10)=stats(10)+ size(Y,1);
  end

  [gxy,isterminal,direction]=switchfun(X,Y);   %   Improve !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  stats(10)=stats(10)+ 1;
  gxy(idis)=0;
     
  for K=2:7,
     HH =  H*C(K);    
     Y1 = Y;
     for J=1:K-1,
        Y1 = Y1 + (H*WRK(:,J)'*A(K,J))';
     end
     [ypro,l1,ff,fplus,yplus,endslid1,gyx,isterminal,direction,stats]=proj(FUN,switchfun,idis,X+HH,Y1,w,1.e-12,ll,tangentside,stats,exactgradient,tol,gradswitchfun);
     if endslid1 > 10^5
         ERR=2.0;
         return
     end
     ll=l1;
     WRK(:,K) = ff';
     auxa= gyx.*gxy<0;
     auxb=gyx.*direction>=0;
     if (endslid0>=0 && endslid1 > min(5.e-7, max([tol 9.e-13]))) || (endslid0<0 && endslid1> 0) || any(auxa.*auxb),  
        condtrans=any(auxa.*auxb);
        xdis=X+(HH+H*min(C(K-1),8/9.))/2;
        if K>=4
           teta0=C(K-1);
           teta1=C(K);
           for it=1:2
             teta=(teta0+teta1)/2;
             acon=[(18*teta - 75*teta^2 + 100*teta^3)/18;-5*(-9*teta^2 + 20*teta^3)/6;10*(-3*teta^2 + 10*teta^3)/9];
             Y1=Y+H*WRK(:,1:3)*acon;
             [ypro,l1,ff,fplus,yplus,endslid2,gyx,isterminal,direction,stats]=proj(FUN,switchfun,idis,X+H*(teta0+teta1)/2,Y1,w,1.e-12,ll,tangentside,stats,exactgradient,tol,gradswitchfun);
             auxa= gyx.*gxy<0;
             auxb=gyx.*direction>=0;
             if(endslid0>=0 && endslid2 >  min(5.e-7, max([tol 9.e-13]))  ) || (endslid0<0 && endslid2> 0)  || any(auxa.*auxb), 
                 condtrans=any(auxa.*auxb);
                 teta1=teta;
             else
                 teta0=teta;
             end
           end
           xdis=X+H*(teta0+teta1)/2;
        end
        if condtrans, % A new discontinuity has been crossed
            disctype=5;
        else       % end of sliding found
            disctype=5;
        end
        ERR=3;
        return;
     end
   end  %end del FOR
   
   EST= zeros(N,1);
 %  Y1 = Y;
  
   for K=1:7,
     for I=1:N,
       EST(I) = EST(I)+H*WRK(I,K)*(B(K)-B1(K));
     end
   end
   ERR=max(abs(EST));
   Y1=ypro;
   WRKout=WRK; 
end %  End of function RKpro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Funcion Slide
function [WRKout,xout,yout,HOLD,XOLD,YOLD,WRKOLD,disctype,xdis,xx,yy,stats]=slide(FUN,switchfun,H,X,WRK,Y,EABS,EREL,xx,yy,XEND,idis,endslid,tangentside,stats,Refine,exactgradient,gradswitchfun)   %Funcion slide
% -----------------------------------------------------------------------
% A Fiipov discontinuity has been detected at X
% This function integrates the problem along the sliding region
% 
%       The integration follows until
%               a non sliding point is detected
%               X reaches XEND or
%               a too short stepsize is needed or a maximum number of
%               steps are reached.
% ----------------------------------------------------------------------
%  TODO!!!!   Hay que mejorarla para el caso en que en el primer paso
%  se detecte el punto de salida
%  Podría habeer otra discontinuidad en esta region, y 
%  la salida podría ser distinta  !!!
%
% Input arguments
%    FUN, switchfun   function handlers for vector field and switching surfaces
%    H
%    X, Y,  initial point to start this part of the integration
%    XOLD, YOLD,  penultimate successful integration point
%    XEND    end point of the integration interval
%    WRK    vector containing the data of the last accepted step  
%    EABS, EREL  absolute and relative error bounds
%    xx, yy   vectors containing the integration points
%    stats    vector containing the statistics of the integration
%    endslid    variable to indicate if the point is attratctive  (<0), 
%         tangent (=0) or transversal (=0)
%    exactgradient    to specify if there is a function for the gradient 
%    gradswitchfun  Function handler, computes the gradient of switchfun
%
% Output arguments
%    WRKout    vector containing the data of the last accepted step
%    xx, yy    vectors containing the integration points
%    xout,  yout   integration (successful) point
%    HOLD   las succesful step size
%    XOLD, YOLD,  penultimate successful integration point
%    WRKOLD   vector containing the data of the penultimate accepted step
%    disctype   type of discontinuity point
%    xdis    point where a new discontinuity has been located
%    idis   hypersurface where the discontinuity is located  
%    stats  vector containing the statistics of the integration 
%
  REJECT=false;
  advance = true;
  xout=X;
  yout=Y;
  disctype=-1;
  HOLD=H;
  XOLD=X;
  YOLD=Y;
  WRKOLD=WRK;
  WRKout=WRK;
  jj=size(idis,2);
  while idis(jj) >0
      jj=jj-1;
  end
  if jj>0
       inddis=abs(idis(jj));
  else
       fprintf('\n Wrong  sliding area ');
  end
  while (advance && X < XEND), 
      H=min(H,XEND-X);
     TOL= EABS + EREL* max(abs(Y));
    if H < 5*eps(X)
         fprintf('\n Minimum step size  h=%g when sliding attained at X= %g \n  Integration stopped \n', H, X);
         disctype=-1;
         return;
    end
    [XPH,Y1,WRKout,ERR,disctype,xdis,stats]=rkpro(FUN,switchfun,inddis,X,H,Y,WRK,stats,exactgradient,TOL,endslid,tangentside,gradswitchfun);
    WRK=WRKout;
    if disctype==3 || xout==X,  % There is not an exit point or it is a first step
       TOL= EABS + EREL* max(abs(Y1));
       if (ERR <= TOL),
% ----------------------------------------------------------------------
%       Accepted sliding step
%         fprintf('\n aceptado slide %g   [%30.20e, %g]', H, X, X+H);
% ----------------------------------------------------------------------
         [XOLD, HOLD, YOLD, WRKOLD]=STORE(X, H, Y, WRK);
%            if Refine>=2
%               for irefine=1:Refine-1
%                  xrefine=X+H*irefine/Refine;
%                  [WRK, Y2]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xrefine, 1,stats,exactgradient,TOL,gradswitchfun); 
%                  xx=[xx xrefine];
%                  yy=[yy ; Y2'];
%               end
%            end
% ----------------------------------------------------------------------
%       Update the numerical approximation, the current point, the
%       first evaluation of the RK method at the next step and the
%       next stepsize to be used.
% ----------------------------------------------------------------------
         WRK(:,1) = WRK(:,7);
         FAC = min(0.9*(TOL/(ERR+1D-16))^(1./5.), 2.0D0);
         if (REJECT), 
           FAC = min(FAC, 1.0D0);
         end         
         Y=Y1;
         X   = XPH;
         H       = FAC*H;
         REJECT  = false;
         xx=[xx X];
         yy=[yy ; Y'];
         stats(3)=stats(3)+1; %  Accepted sliding steps counter
       else
% ----------------------------------------------------------------------
%       Rejected sliding step
%        fprintf('\n rechazo slide %g %g %g %g %g\n',X,H, WRK(:,1));
% ----------------------------------------------------------------------
         FAC = max(0.9*(TOL/(ERR+1D-12))^(1./5.), 0.10D0);
         REJECT  = true;
         H       = FAC*H;
         stats(4)=stats(4)+1;
       end
    else  % Exit of sliding detected at this step.
       stats(5)=stats(5)+1;  % Discont. detected counter
%   fprintf('\n rechazo discon slide %g   [%30.20e, %g %g]', H, X, X+H, xdis);
       if xdis-X < 0.15*HOLD, %Estirapro can be directly applied
            xout=X;
            yout=Y;
         %   disctype=4;
            advance=false;
       else  %Estirapro can not be directly applied. Repeating a step
   %        H=0.9*(xdis-X);  % Va mal con el problema 162, tol=1.e-12 y 1.e-9
           H=0.95*(xdis-X); % Va mal con el problema 162, tol=1.e-10 
   %        H=0.97*(xdis-X); % Va mal con el problema 162, tol=1.e-10 
           REJECT  = true;
       end
    end
  end   %  End of While
  if XEND <= X,
      disctype=-1;
  end
  WRKout=WRK;
  xout=X;
end %  End of function slide

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WRKout,xout,yout,H,disctype,tdis, ydis,idis,stats]=FindDiscpro(FUN,switchfun,HOLD,X,XOLD,WRKOLD,Y,YOLD, ...
    xdis,tdis,ydis,idis, endslid,tangentside,stats,exactgradient,tol,gradswitchfun,doatswitch)
% -----------------------------------------------------------------------
% Inside the sliding area, a new discontinuity has been detected in
%  [X, 2*xdis-X]
% Taking the step from XOLD to X with step-size HOLD,
% the discotinuity point is located using the function estirapro
% -----------------------------------------------------------------------
%
%  Utiliza el metodo de la secante.  Afinar la terminacion, teniendo en
%  cuenta la tolerancia del usuario, para reducir iteraciones ???
%
  xout=X;
  xmax=max(2*xdis-X, X+0.3*HOLD);
  
  format=' X= %26.20e, Y= (%27.20e';
  for ii=1:size(Y,1)-1
     format=[format ', %27.20e'];
  end  
  format=[format ')'];
       
  jj=size(idis,2);
  while idis(jj) >0
     jj=jj-1;
  end
  if jj>0
     inddis=abs(idis(jj));   %  sliding surface index
  else
     fprintf('\n Wrong  sliding area ');
  end 
  if exactgradient
     w=feval(gradswitchfun,X,Y,inddis);
     w=w/norm(w);
     stats(11)=stats(11)+ 1;
  else
     w=graddif(switchfun,X,Y,inddis,tol);
     w=w/norm(w);
     stats(10)=stats(10)+size(Y,1);
  end
%
%   Lower extreme of the interval
%
  x0=X;
  [WRKout, Y1,f1,ypro,endslid0,gg0,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,x0,1,tangentside,stats,w,exactgradient,tol,gradswitchfun);
  gg0(inddis)=1000;   %  To eliminate the current sliding surface
  value0=min(abs(gg0(gg0~=0)));
  WRKOLD=WRKout; 
% 
%   Upper extreme of the interval
%
  x1=xmax;
  [WRKout,Y1,f1,ypro,endslid1,gg1,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,x1, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun); 
  gg1(inddis)=1000;   %  To eliminate the current sliding surface
  gcond=(gg1.*gg0<0).*(gg1.*direction>=0);
  endtangent=endslid>=0 && endslid1 > min(5.e-7,max([tol 9.e-13]));
  endfilippov=endslid<0 && endslid1>0;
  newdiscon=any(gcond);
  if (endtangent|| endfilippov) && newdiscon,    % Two kinds of simultaneous switching points
         value1=-max([max(abs(gg1(gcond==1))), abs(endslid1)]); 
         ind=find(gcond==1);
         switchpoint=4;
         bound=max(100*tol, 9.e-13);
  elseif newdiscon,     % new discontinuity point
         value1=-min(abs(gg1(gcond==1)));
         ind=find(gcond==1);
         switchpoint=1;
         bound=max(100*tol, 9.e-13);
  elseif endtangent     % possible end of tangent point
         value1=-abs(endslid1);
         switchpoint=3;
         bound=min(5.e-7,max(tol, 9.e-13));
  elseif endfilippov  % possible end of sliding point
         value1=-abs(endslid1);
         switchpoint=2;
         bound=0;
  else
%       fprintf('\n Discon fallida !!!! %g %g ',x0,x1);
       yout=Y;
       disctype=3;
       H=HOLD/4;
       return
  end 
  if switchpoint==3  && endslid0>=0 
     tdis=[tdis X];
     ydis=[ydis; Y'];
     idis=[idis, -inddis];
     stats(8)=stats(8)+1;
     xplus=X;
     y2=Y;
     f2=WRKOLD(:,7);
     hh=1.e-14;
     for i=1:4
        y2=y2+hh*f2;
        xplus=xplus+hh;
        f2=FUN(xplus,y2);
        stats(9)=stats(9)+1;
     end
     slidp=endslid0;
     xout=xplus;
     yout=y2;
     WRKout(:,1)=f2;
%     fprintf(['\n Exit of tangent sliding ' format], xout, yout);
     disctype=0;
     H=HOLD;
     return
  end
%
%   Secant loop, modified to ensure convergence
%  
  xx1=x0;  % Lower limit of the interval containing the discontinuity
  xx2=x1;  % Upper limit of the interval containing the discontinuity
  ii=1;
  value=value1;
  while abs(x1-x0) >= 10*eps(max(abs(x1),1))  &&  abs(value)>= 1.e-30 && ii< 100,
 % for II=1:100,
     if value1~=value0
       xnew=x1-value1*(x1-x0)/(value1-value0);
     else
       xnew=(xx1+xx2)/2;
     end
     if xnew >=xx2 || xnew <= xx1
       xnew=(xx1+xx2)/2; % If secant method leads out the interval, use bisection
     end
     [WRKout,ynew,f1,yplus,endslid1,gnew,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xnew, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun);
     gnew(inddis)=1000;
     gcond=(gnew.*gg0<0).*(gnew.*direction>=0);
     gcondp=(gnew.*gg0<=0).*(gnew.*direction>=0);
     gcondp(inddis)=0;
     endtangent=endslid>=0 && endslid1 > bound;
     endfilippov=endslid<0 && endslid1>0;
     newdiscon=any(gcond);
     if (endtangent|| endfilippov) && newdiscon, 
         value=-max([max(abs(gnew(gcond==1))), abs(endslid1-bound)]); % current point is after the discontinuity
         switchpoint=4;
     elseif newdiscon
         value=-max(abs(gnew(gcond==1)));
         switchpoint=1;
     elseif endfilippov,
         value=-abs(endslid1-bound);
         switchpoint=2;
     elseif endtangent,
         value=-abs(endslid1-bound);
         switchpoint=3;
     elseif switchpoint==4       % current point is beforer the discontinuity
         value=min(min(abs(gnew(ind))), abs(endslid1-bound));
     elseif switchpoint==3,
         value=abs(endslid1-bound);
     elseif switchpoint==2
         value=abs(endslid1-bound);
     else
         value=min(abs(gnew(ind)));
     end
%
%  Updating the interval where the discontinuity is located
%
     if value <0
        xx2=xnew;
     else
        xx1=xnew;
     end
%
%  Updating the two last points for the secant method
%
     x0=x1;
     value0=value1;
     x1=xnew;
     value1=value;
     ii=ii+1;
  end %end of for II=1:32,
%      fprintf('\n Secant method  iter=%g %30.20e %30.20e %30.20e',II, x1-x0, gnew);
%  fprintf(['\n  discontinuidad en ' format], xnew, ynew );
  tdis=[tdis xnew];
  ydis=[ydis; ynew'];
%
%  Numerical discontinuity found at (xnew, ynew)
%
if ii >=100
    '  se pasa del bucle !!! '
end

%
%  Once the discontinuity point is computed, another point, very close
%  but after it, is computed.
%  Estirpro give us also a point yplus at the side of the manifold where the 
%  integration must continue
%
   xminus=xnew;
   yminus=ynew;
   xplus=xnew;
   yplus=ynew;
   f2=f1;
   gplus=gnew;
   endslid0=endslid1;
   
   if switchpoint==2,    %   Possible end of Filippov region
      fplus=f1;
      xplus=xnew;   
      if (endslid<=0), 
         ii=1;
         while (endslid<=0) && ii<6
           xplus=xnew+8^ii*eps(xnew);
           [WRKout,Y2,fplus,yplus,endslid,gyx,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xplus, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun);
           ii=ii+1;
         end   
      end
      y2=yplus;
      f2=fplus;
      hh=max(1.e-11,10*eps(xplus));
      for j=1:8
         for i=1:2
            y2=y2+hh*f2;
            xplus=xplus+hh;
            f2=FUN(xplus,y2);
            stats(9)=stats(9)+1;
         end
         [gg2,isterminal,direction]=switchfun(xplus,y2);
         stats(10)=stats(10)+1;
         if abs(gg2) > max(10*abs(gyx), 1.e-15),
            break
         end        
      end
      idis=[idis, -inddis];
      stats(8)=stats(8)+1;
      slidp=endslid;
      xout=xplus;
      yout=y2;
      WRKout(:,1)=f2;
%      fprintf(['\n Exit of sliding ' format], xout, yout);
      disctype=0;
      H=HOLD;
      return
   elseif switchpoint==3,   %  Possible end of tangent sliding
      fplus=f1;
      xplus=xnew;   
      if endslid1 <= min(5.e-7, max([tol 9.e-13])), 
         ii=1;
         while endslid1 <= min(5.e-7, max([tol 9.e-13])) && ii<6
           xplus=xnew+8^ii*eps(xnew);
           [WRKout,Y2,fplus,yplus,endslid1,gyx,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xplus, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun);
           ii=ii+1;
         end   
      end
      idis=[idis, -inddis];
      stats(8)=stats(8)+1;
      slidp=endslid;
      xout=xplus;
      yout=yplus;
      WRKout(:,1)=fplus;
%      fprintf(['\n Exit of tangent  sliding ' format], xout, yout);
      disctype=0;
      H=HOLD;
      return
   else
      if switchpoint==1,     %  Possible new discon
%
%  Compute a point just after the discontinuity if (xnew,ynew) is not
%
         if  all(gcond==0), % means that the numerical discontinuity point is in fact before it
            ii=1;
            while all(gcond==0) && ii <8  % the point is not estrictly after the discontinuity
               xplus=xplus+10^ii*eps(xnew);
               [WRKout,Y2,f2,yplus,endslid1,gplus,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xplus, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun);   
               gplus(inddis)=1000;
               gcond=(gg0.*gplus<0).*(gg0.*direction<=0);
               ii=ii+1;
            end
            if ii >=8,
               fprintf('\n\n  Warning !!!   yplus could have failed\n ');
            end
         end
         if  any(gcondp) % means that the numerical discontinuity point is in fact after it
            ii=1;
            while any(gcondp)  && ii <8  % the point is not estrictly before the discontinuity
                xminus=xminus-10^ii*eps(xnew);
                [WRKout,Y2,f1,yminus,endslid0,gminus,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xminus, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun);  
                gminus(inddis)=1000;
                gcondp=(gg0.*gminus<=0).*(gg0.*direction<=0);
                gcondp(inddis)=0;
                ii=ii+1;
            end
            if ii >=8,
               fprintf('\n\n  Warning !!!   yminus could have failedo\n ');
            end
         end    
         %   Comprobar si tambien sale de sliding !!
      else     %  Possible new discon, may be with end of tangent or Filippov sliding

%
%  Compute a point just after the discontinuity if (xnew,ynew) is not
%
         if ((endslid>=0 && endslid1 <= min(5.e-7, max([tol 9.e-13])) ) || (endslid<0 && endslid1<= 0)) && all(gcond==0), % means that the numerical discontinuity point is in fact before it
            ii=1;
            while( ((endslid>=0 && endslid1 <= min(5.e-7, max([tol 9.e-13])) ) || (endslid<0 && endslid1<= 0)) && all(gcond==0)) && ii <10  % the point is not estrictly after the discontinuity
               xplus=xplus+10^ii*eps(xnew);
               [WRKout,Y2,f2,yplus,endslid1,gplus,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xplus, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun);   
               gplus(inddis)=1000;
               gcond=(gg0.*gplus<0).*(gg0.*direction<=0);
               ii=ii+1;
            end
            if ii >=8,
               fprintf('\n\n  Warning !!!   yplus could have failed\n ');
            end   
         end

         if  (endslid>=0 && endslid0 > min(5.e-7,max([tol 9.e-13])) ) || (endslid<0 && endslid0> 0) || any(gcondp) % means that the numerical discontinuity point is in fact after it
            ii=1;
            while ( (endslid>=0 && endslid1 > min(5.e-7,max([tol 9.e-13])) ) || (endslid<0 && endslid1> 0) || any(gcondp))  && ii <10  % the point is not estrictly before the discontinuity
                xminus=xminus-10^ii*eps(xnew);
                [WRKout,Y2,f1,yminus,endslid0,gminus,isterminal,direction,stats]=estirapro(FUN,switchfun,inddis,XOLD,YOLD,WRKOLD,HOLD,xminus, 0,tangentside,stats,w,exactgradient,tol,gradswitchfun);  
                gminus(inddis)=1000;
                gcondp=(gg0.*gminus<=0).*(gg0.*direction<=0);
                gcondp(inddis)=0;
                ii=ii+1;
            end
            if ii >=8,
                fprintf('\n\n  Warning !!!   yminus could have failed\n ');
            end   
         end
      end
%        
% Checking if the discontinuity point is Filipov or transversal
%  
                  
   if switchpoint==1 || switchpoint==4,
        uu=ind(1);
   else
        uu=inddis;
   end
   if exactgradient
      gt=gradt(switchfun,xnew,ynew,uu,tol); 
      gfx=feval(gradswitchfun,xnew,ynew,uu);
      stats(11)=stats(11)+ 1;
      gfxf1=(gt+gfx'*f1)/norm(gfx);
      gfxf2=(gt+gfx'*f2)/norm(gfx);
   else
      gt=gradt(switchfun,xnew,ynew,uu,tol); 
      gfxf1=gt+graddif(switchfun,xnew,ynew,uu,tol,f1);
      gfxf2=gt+graddif(switchfun,xnew,ynew,uu,tol,f2);
      if tol > 1.e-9
          stats(10)=stats(10)+6;
      else
          stats(10)=stats(10)+12;
      end
   end
   if gfxf1*gfxf2<0,
%        
% A Filipov point is detected
%
       idis=[idis -ind(1)];
       fprintf('\n co-dimension 2 Filippov point at X= %g Y= (%30.20e, %g %g %g)', xnew, ynew);
       stats(7)=stats(7)+1;
       disctype=-5;
       alfa=gfxf1/(gfxf1-gfxf2);
       WRKout(:,1)=(1-alfa)*f1+alfa*f2;
       yout=ynew;
       xout=xnew;
   elseif (endslid>=0 && endslid1 > min(5.e-7, max([tol 9.e-13]) )) || (endslid<0 && endslid1>0)
       disctype=0;
       xout=xplus;
       yout=yplus;
       WRKout(:,1)=f2;
%       fprintf(['\n salimos de sliding ' format ' %g'], xout, yout, min(switchfun(xout,yout).*gg0));
       idis=[idis, -inddis];
       stats(8)=stats(8)+1;
   else
%
%  Transversal discontinuity.  The integration proceeds from the point
%  just after it, (xplus,yplus)
%
       idis=[idis uu];
       stats(6)=stats(6)+1;
       xout=xplus;
       yout=yplus;
       WRKout(:,1)=f2;
       disctype=3;
%       fprintf(['\n salimos de ' format ' %g'], xout, yout, min(switchfun(xout,yout).*gg0));
   end
   if any(isterminal((gplus.*gg0<=0) & (direction.*gg0<=0))==1)  % the switching point ends the integration
       disctype=-4;
   end
   if any(isterminal((gplus.*gg0<=0) & (direction.*gg0<=0))==-1) % the switching point requires calling the external function
       xout=xnew;
       yout=doatswitch(xnew,ynew);
       WRKout(:,1)=FUN(xout,yout);
       stats(9)=stats(9)+1;
   end
   H=HOLD/2;
   return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Funcion estirapro
function [WRKout, Y1,f1,yplus,endslid,gxy,isterminal,direction,stats]=estirapro(FUN,switchfun,idis,X,Y,WRK,H,XA, addstages,tangentside,stats,w,exactgradient,tol,gradswitchfun)  %Funcion estirapro
% -----------------------------------------------------------------------
% This function uses the continuous extension  of the step
% [X, X+H] to get an approximation of order five at the
% point X+ H+ TETA*(XMAX-X), with TETA in [0,1],
% but proyecting onto the sliding manifold
% -----------------------------------------------------------------------
%
  ll=0; 
  if (addstages==1), 
%
%  First additional stage for computing the 5th order continuous
%  Saved at WRK(:,8)
%
      Y1=Y+H*(-24018683.D0/8152320000.D0*WRK(:,1)+ ... 
            25144.D0/43425.D0*WRK(:,2)- ...
            76360723.D0/337557000.D0*WRK(:,3)+ ...
            349808429.D0/2445696000.D0*WRK(:,4)- ...
            13643731773.D0/144024320000.D0*WRK(:,5)+ ...
            1.D0/20.D0*WRK(:,6)- ...
            12268567.D0/254760000.D0*WRK(:,7));   
     [ypro,l1,ff,f1,yplus,endslid,gxy,isterminal,direction,stats]=proj(FUN,switchfun,idis,X+0.4d0*H,Y1,w,1.e-10,ll,tangentside,stats,exactgradient,tol,gradswitchfun);
     WRK(:,8) = ff';
%
%  Second additional stage for computing the 5th order continuous
%  Saved at WRK(:,9)
%
      Y1=Y+H*(2104901.D0/23010000.D0*WRK(:,1)+ ...
            54324224.D0/106708875.D0*WRK(:,3)+ ...
            134233.D0/2301000.D0*WRK(:,4)- ...
            13268529.D0/406510000.D0*WRK(:,5)+ ...
            26972.D0/2013375.D0*WRK(:,6)- ...
            6324.D0/479375.D0*WRK(:,7)- ...
            1737.D0/7670.D0*WRK(:,8));
     [ypro,l1,ff,f1,yplus,endslid,gxy,isterminal,direction,stats]=proj(FUN,switchfun,idis,X+0.4d0*H,Y1,w,1.e-10,ll,tangentside,stats,exactgradient,tol,gradswitchfun);
     WRK(:,9) = ff';
  end %end IF
%
  t=(XA-X)/H;
%
%
%  Coefficients dof the 5th order continuous
%
  co5(1) = -1809809.D0/441792.D0*t^2+1150189.D0/147264.D0*t^3-2025025.D0/ ...
           294528.D0*t^4+995255.D0/441792.D0*t^5+t;
  co5(2) = 0.0d0;
  co5(3) = 49017100.D0/853671.D0*t^4-54769600.D0/2561013.D0*t^5+40963600.D0/ ...
           2561013.D0*t^2-14677200.D0/284557.D0*t^3;
  co5(4) = 1189575.D0/49088.D0*t^4-1161175.D0/73632.D0*t^3-834475.D0/ ...
           73632.D0*t^5+259225.D0/73632.D0*t^2;
  co5(5) = 21163599.D0/2601664.D0*t^3-64133775.D0/5203328.D0*t^4+14882535.D0/ ...
           2601664.D0*t^5-4817961.D0/2601664.D0*t^2;
  co5(6) = -53416.D0/16107.D0*t^3+323345.D0/64428.D0*t^4-112475.D0/ ...
           48321.D0*t^5+36542.D0/48321.D0*t^2;
  co5(7) = -14140.D0/2301.D0*t^4+7270.D0/2301.D0*t^5-1901.D0/ ...
           2301.D0*t^2+8771.D0/2301.D0*t^3;
  co5(8) = 120625.D0/6136.D0*t^3-120625.D0/6136.D0*t^4+120625.D0/ ...
           18408.D0*t^5-120625.D0/18408.D0*t^2;
  co5(9) = -125.D0/3.D0*t^4+625.D0/36.D0*t^5-125.D0/18.D0*t^2+125.D0/4.D0*t^3;
  
%
%  Coefficients of the 4th order continuous
%
%      co(1) = -435.D0/384.D0*t^4+1184.D0/384.D0*t^3-1098.D0/384.D0*t^2+t;
%      co(2) = 0.0d0;
%      co(3) = 500.D0*t^2*(6.D0*t^2-14.D0*t+9.D0)/1113.D0;
%      co(4) = -125.D0*t^2*(9.D0*t^2-16.D0*t+6.D0)/192.D0;
%      co(5) = 729.D0*t^2*(35.D0*t^2-64.D0*t+26.D0)/6784.D0;
%      co(6) = -11.D0*t^2*(3.D0*t-2.D0)*(5.D0*t-6.D0)/84.D0;
%      co(7) = t^2*(t-1.D0)*(5.D0*t-3.D0)/2.D0;    
%      cont=WRK(:,1:7)*co(1:7)';
%      
% 
%  Computation of the approximation at X+t*H
% 
  cont=WRK*co5';
  Y1=Y+H*cont;

  [ypro,~,~,f1,yplus,endslid,gxy,isterminal,direction,stats]=proj(FUN,switchfun,idis,XA,Y1,w,1.e-10,ll,tangentside,stats,exactgradient,tol,gradswitchfun);
  Y1=ypro;
  WRKout = WRK;
end %fin de funcion ESTIRApro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Funcion ESTIRA
function [WRKout, Y1]=ESTIRA(FUN,X,Y,WRK,H,XA, addstages)  %Function ESTIRA
% -----------------------------------------------------------------------
% This function uses the continuous extension  of the step
% [X, X+H] to get an approximation of order five at the
% point XA, with XA in [X,X+H+0.3*H]
% -----------------------------------------------------------------------
%
%   FUN  
%   X,Y
%   H
%   WRK,
%   XA
%   addstages  
%       If addstages>0   the fifth order continuous is used
%       If addstages=1   The two additional stages are computed
%       If addstages <0  the fourth order continuous is used
%
%   WRKout   
%   Y1  approximation at XA
%
 if (addstages==1), 
%
%  First additional stage for computing the 5th order continuous
%  Saved at WRK(:,8)
%
     Y1=Y+H*(-24018683.D0/8152320000.D0*WRK(:,1)+ ... 
            25144.D0/43425.D0*WRK(:,2)- ...
            76360723.D0/337557000.D0*WRK(:,3)+ ...
            349808429.D0/2445696000.D0*WRK(:,4)- ...
            13643731773.D0/144024320000.D0*WRK(:,5)+ ...
            1.D0/20.D0*WRK(:,6)- ...
            12268567.D0/254760000.D0*WRK(:,7));
    WRK(:,8)=FUN(X+0.4d0*H, Y1);
%
%  Second additional stage  for computing the 5th order continuous
%  Saved at WRK(:,9)
%
     Y1=Y+H*(2104901.D0/23010000.D0*WRK(:,1)+ ...
            54324224.D0/106708875.D0*WRK(:,3)+ ...
            134233.D0/2301000.D0*WRK(:,4)- ...
            13268529.D0/406510000.D0*WRK(:,5)+ ...
            26972.D0/2013375.D0*WRK(:,6)- ...
            6324.D0/479375.D0*WRK(:,7)- ...
            1737.D0/7670.D0*WRK(:,8));
    WRK(:,9)=FUN(X+0.4d0*H, Y1); 
 end %end IF
%
%
 t=(XA-X)/H;
 if addstages>=0
%
%  Coefficients of the 5th order continuous
%
     co(1) = -1809809.D0/441792.D0*t^2+1150189.D0/147264.D0*t^3-2025025.D0/ ...
           294528.D0*t^4+995255.D0/441792.D0*t^5+t;
     co(2) = 0.0d0;
     co(3) = 49017100.D0/853671.D0*t^4-54769600.D0/2561013.D0*t^5+40963600.D0/ ...
           2561013.D0*t^2-14677200.D0/284557.D0*t^3;
     co(4) = 1189575.D0/49088.D0*t^4-1161175.D0/73632.D0*t^3-834475.D0/ ...
           73632.D0*t^5+259225.D0/73632.D0*t^2;
     co(5) = 21163599.D0/2601664.D0*t^3-64133775.D0/5203328.D0*t^4+14882535.D0/ ...
           2601664.D0*t^5-4817961.D0/2601664.D0*t^2;
     co(6) = -53416.D0/16107.D0*t^3+323345.D0/64428.D0*t^4-112475.D0/ ...
           48321.D0*t^5+36542.D0/48321.D0*t^2;
     co(7) = -14140.D0/2301.D0*t^4+7270.D0/2301.D0*t^5-1901.D0/ ...
           2301.D0*t^2+8771.D0/2301.D0*t^3;
     co(8) = 120625.D0/6136.D0*t^3-120625.D0/6136.D0*t^4+120625.D0/ ...
           18408.D0*t^5-120625.D0/18408.D0*t^2;
     co(9) = -125.D0/3.D0*t^4+625.D0/36.D0*t^5-125.D0/18.D0*t^2+125.D0/4.D0*t^3;
     cont=WRK(:,1:9)*co(1:9)';
 else
%
%  Coefficients of the 4th order continuous
%
     co(1) = -435.D0/384.D0*t^4+1184.D0/384.D0*t^3-1098.D0/384.D0*t^2+t;
     co(2) = 0.0d0;
     co(3) = 500.D0*t^2*(6.D0*t^2-14.D0*t+9.D0)/1113.D0;
     co(4) = -125.D0*t^2*(9.D0*t^2-16.D0*t+6.D0)/192.D0;
     co(5) = 729.D0*t^2*(35.D0*t^2-64.D0*t+26.D0)/6784.D0;
     co(6) = -11.D0*t^2*(3.D0*t-2.D0)*(5.D0*t-6.D0)/84.D0;
     co(7) = t^2*(t-1.D0)*(5.D0*t-3.D0)/2.D0;    
     cont=WRK(:,1:7)*co(1:7)';
 end
% 
%  Computation of the approximation at X+t*H
%
  
  Y1=Y+H*cont;
  WRKout = WRK; 
 %Saved at the output vector
end %End of function ESTIRA

%Funcion STORE
function [XOLD, HOLD, YOLD, WRKOLD]=STORE(X, H, Y, WRK)             %Funcion STORE
% ----------------------------------------------------------------------
%
%       Subroutine to store the vectors needed to calculate the
%       continuous Runge-Kutta method of fifth order and the points
%       taken by the integrator.
%
% ----------------------------------------------------------------------
  XOLD=X;
  HOLD=H;
  YOLD=Y;
  WRKOLD=WRK;
end %fin de funcion STORE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ypro,l1,ff,f1,yplus,endslid,gpro,isterminal,direction,stats]=proj(FUN,switchfun,idis,x,y,w,tol,ll,tangentside,stats,exactgradient,tolrk,gradswitchfun)
% -----------------------------------------------------------------------
% This function projects a point (x,y) of the approximated solution onto 
% the idis component of the switching surface switchfun
% -----------------------------------------------------------------------
%
% Input arguments
%
%   FUN  
%   switchfun
%   idis   component of the switching surface where the solution is sliding
%   x,y    point to be projected
%   w   direction of projection
%   tol   error tolerance for the projection
%   l1   initial value of the parameter
%   stats    statistics vector
%   exactgradient    to specify if there is a function for the gradient
%   tolrk   error tolerance the code is using for the integration
%   gradswitchfun    gradient function, if exactgradient=1
%
% Output arguments
%
%   ypro     projected approximation   
%   l1   parameter finally used for the projection
%   ff  Filippov vector field at the projected point
%   yplus, f1  point and vector field for exiting the sliding region
%   endslid   indicator of exiting sliding region
%   [gpro,isterminal,direction]=   output of switchfun(x,ypro)
%   stats    updated statistics vector
%
  l0=ll;  %  Initial approximation for the parameter
%l0=0;
  ypro=y+l0*w; %   Initial approximation for the projected solution
  er=tol+1;
  ns=0;   %number of iterations
  [gpro,isterminal,direction]=switchfun(x,ypro);
  stats(10)=stats(10)+1;
%   if gpro(idis)==0
%       gproa=gpro;
%       l1=0;
%   else
  if exactgradient  % computing the direction of projection
      gt=gradt(switchfun,x,ypro,idis,tol); 
      grad=feval(gradswitchfun,x,ypro,idis);
      stats(11)=stats(11)+ 1;
      den=grad'*w;
  else
      gt=gradt(switchfun,x,ypro,idis,tol); 
      den=graddif(switchfun,x,ypro,idis,tol,w);
      if tol > 1.e-9
          stats(10)=stats(10)+4;
      else
          stats(10)=stats(10)+8;
      end
  end
%
%   Simplified Newton iteration
%
  while abs(er)>tol && ns < 30
    l1=l0-gpro(idis)/den;
    era=er;
    er=l1-l0;
    l0=l1;
    yproa=ypro;
    gproa=gpro;
    ypro=y+l1*w;
    ns=ns+1;
    [gpro,isterminal,direction]=switchfun(x,ypro);
    stats(10)=stats(10)+1;
    if abs(er)> 5 || abs(er)>abs(era)  % possible divergence
        endslid=10^6;
        ff=y;  %just to have an output (nonsense) value
        f1=y;
        yplus=y;
        return
    end
  end
  if ns >= 29   % no convergence
   fprintf('\n no projection after 30 iter');
    endslid=10^6;
    ff=y;
    f1=y;
    yplus=y;
    return
  end
%  end
  if exactgradient
    aux=feval(gradswitchfun,x,ypro,idis);
    stats(11)=stats(11)+ 1;
  else
    aux=w;
  end
%
%   computation two points at each side of the switching surface
%
 if gpro(idis)==0 && gproa(idis) ==0  % the two last iterations are exact
   yminus=ypro+1.e-14*aux;
   [g1,isterminal,direction]=switchfun(x,yminus);
   gminus=g1(idis);
   stats(10)=stats(10)+1;
   ii=2;
   while g1(idis)==0 && ii<5
     yminus=ypro+10^ii*1.e-14*aux;
     [g1,isterminal,direction]=switchfun(x,yminus);
     gminus=g1(idis);
     stats(10)=stats(10)+1;
     ii=ii+1;
   end
   yplus=ypro-1.e-14*aux;
   [g1,isterminal,direction]=switchfun(x,yplus);
   gplus=g1(idis);
   stats(10)=stats(10)+1;
   ii=2;
   while g1(idis)==0 && ii<5
     yplus=ypro-10^ii*1.e-14*aux;
     [g1,isterminal,direction]=switchfun(x,yplus);
     gplus=g1(idis);
     stats(10)=stats(10)+1;
     ii=ii+1;
   end
 elseif gpro(idis)*gproa(idis) <0, % the two last iterations are at both sides
   if gpro(idis)<0
       yminus=ypro;
       gminus=gpro(idis);
       yplus=yproa;
       gplus=gproa(idis);
   else
       yplus=ypro;
       gplus=gpro(idis);
       yminus=yproa;
       gminus=gproa(idis);
   end
 elseif gpro(idis)==0 % the last iteration is exact. We need one
   yminus=yproa;
   gminus=gproa(idis);
   yplus=ypro+sign(er)*1.e-14*aux;
   [g1,isterminal,direction]=switchfun(x,yplus);
   gplus=g1(idis);
   stats(10)=stats(10)+1;
   ii=2;
   while  g1(idis)*gproa(idis) >=0 && ii<5
      yplus=ypro+sign(er)*10^ii*1.e-14*aux;
      [g1,isterminal,direction]=switchfun(x,yplus);
      gplus=g1(idis);
      stats(10)=stats(10)+1;
      ii=ii+1;
   end   
 else  % the penultimate iteration is exact.  Should be impossible ?
   yminus=ypro;
   gminus=gpro(idis);
   yplus=ypro+sign(er)*1.e-13*aux;
   [g1,isterminal,direction]=switchfun(x,yplus);
   gplus=g1(idis);
   stats(10)=stats(10)+1;
   ii=2;
   while  g1(idis)*gpro(idis) >=0 && ii<5
      yplus=ypro+sign(er)*10^ii*1.e-14*aux;
      [g1,isterminal,direction]=switchfun(x,yplus);
      gplus=g1(idis);
      stats(10)=stats(10)+1;
      ii=ii+1;
   end   
 end
%
%  Computing the Filippov vector field at the projected point
%
 f1=FUN(x,yplus);
 f2=FUN(x,yminus);
 stats(9)=stats(9)+2;
 noraux=norm(aux);
 if exactgradient,
    gt=gradt(switchfun,x,ypro,idis,tolrk); 
    gfxf1=(gt+aux'*f1)/noraux;
    gfxf2=(gt+aux'*f2)/noraux;
 else
    gt=gradt(switchfun,x,ypro,idis,tolrk); 
    gfxf1=(gt+graddif(switchfun,x,ypro,idis,tolrk,f1))/noraux;
    gfxf2=(gt+graddif(switchfun,x,ypro,idis,tolrk,f2))/noraux;
    if tol > 1.e-9
        stats(10)=stats(10)+6;
    else
        stats(10)=stats(10)+12;
    end
 end
 endslid=gfxf1*gfxf2;    %  To check if there is an exit point
    
 if endslid <0,   %  switching surface is attractive
    alfa=gfxf1/(gfxf1-gfxf2);
    ff=(1-alfa)*f1+alfa*f2;  %  Filippov vector field
    if abs(gfxf1)>=abs(gfxf2),
       f1=f2;
       yplus=yminus;
    end
    endslid=-min([abs(gfxf1) abs(gfxf2)]);
 elseif abs(gfxf1)+abs(gfxf2) < min( 1.e-5, 100*tolrk)  % Vector field Tangent at both sides of the surface
     endslid=max([abs(gfxf1) abs(gfxf2)]);
    yplus=ypro;
    if sign(gplus)==tangentside,
      ff=f1;
    else
      ff=f2;
      f1=f2;
    end
 elseif abs(gfxf1)<abs(gfxf2)  % f1 closer to tangent than f2
    [gpru,isterminal,direction]=switchfun(x,ypro+1.e-6*f2);
    stats(10)=stats(10)+1;
    if gpru(idis)*gminus>0  % flow from yplus to yminus
        yplus=yminus;
        f1=f2;
        ff=f2;
        endslid=abs(gfxf2)/norm(f2);
    else      % flow from yminus  to yplus
        ff=f1;
        endslid=abs(gfxf1)/norm(f1);
    end
 else         % f2 closer to tangent than f1
    [gpru,isterminal,direction]=switchfun(x,ypro+1.e-6*f1);
    stats(10)=stats(10)+1;
    if gpru(idis)*gminus>0  % flow from yplus to yminus
        yplus=yminus;
        f1=f2;
        ff=f2;
        endslid=abs(gfxf2)/norm(f2);
    else      % flow from yminus  to yplus
        ff=f1;
        endslid=abs(gfxf1/norm(f1));
    end
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion graddif
function gxd=graddif(switchfun,t,y,ind,tol,v)                                    %Funcion graddif
 if nargin==6
    N=1;
    vv=v;
 elseif nargin==5
    N=size(y,1);
    vv=eye(N);
 end
 gxd=zeros(N,1);
 for it=1:N
   v=vv(:,it);
   if tol > 1.e-7
       e=1.e-8;
      [gy,~,~]=switchfun(t,y);
      z=y+e*v;
      [gz,~,~]=switchfun(t,z);
      gxdit=(gz(ind)-gy(ind))/1.e-8;
   elseif tol > 1.e-9
      e=1.e-5;
      [gz1,~,~]=switchfun(t,y+e*v/2);
      [gz2,~,~]=switchfun(t,y-e*v/2);
      gxdit=(gz1(ind)-gz2(ind))/e;
   else
      n=max([3, floor(14+log10(tol))]);
      e=10^(-n);
      [gz1,~,~]=switchfun(t,y+e*v);
      [gz2,~,~]=switchfun(t,y+e*v/2);
      [gz3,~,~]=switchfun(t,y-e*v/2);
      [gz4,~,~]=switchfun(t,y-e*v);
      gxdit=(-gz1(ind)/6+4*gz2(ind)/3-4*gz3(ind)/3+gz4(ind)/6)/e;
   end
   gxd(it)=gxdit;
  end
end
 
function gt=gradt(switchfun,t,y,ind,tol)                                    %Funcion gradt
 
   if tol > 1.e-7
       e=1.e-8;
      [g1,~,~]=switchfun(t,y);
      [g2,~,~]=switchfun(t+e,y);
      gt=(g2(ind)-g1(ind))/1.e-8;
   elseif tol > 1.e-9
      e=1.e-6;
      [g1,~,~]=switchfun(t+e/2,y);
      [g2,~,~]=switchfun(t-e/2,y);
      gt=(g1(ind)-g2(ind))/e;
   else
      n=max([3, floor(15+log10(tol))]);
      e=10^(-n);
      [g1,~,~]=switchfun(t+e,y);
      [g2,~,~]=switchfun(t+e/2,y);
      [g3,~,~]=switchfun(t-e/2,y);
      [g4,~,~]=switchfun(t-e,y);
      gt=(-g1(ind)/6+4*g2(ind)/3-4*g3(ind)/3+g4(ind)/6)/e;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Funcion dummy
function yout=dummy(t,y) 
   yout=y;
   fprintf('\n\n Warning !!   No function actionatswitch provided \n');
   fprintf('\n Check the options in the call to DISODESET \n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function classifypoint
function [xout,yout,ff,disctype,idis,endslid,stats]=classifypoint(x,y,xminus,yminus,xplus,yplus,inddis,FUN,switchfun,idis,tangentside,stats,tol,exactgradient)
  endslid=1;
  if xplus==xminus,  
   f0 = FUN(x,y);
   stats(9)=stats(9)+1;
   [g0,~,~]=switchfun(x,y);
   stats(10)=stats(10)+1;
   if all(abs(g0) > 100*eps(1))
       ff = f0;
       xout=x;
       yout=y;
       disctype=0;
       return
   end
   ind=find(abs(g0) <= 100*eps(1));
   inddis=ind(1);
   gxd=graddif(switchfun,x,y,inddis,tol);
   gt=gradt(switchfun,x,y,inddis,tol); 
   stats(10)=stats(10)+ size(y,1)+2;
   yminus=y;
   yplus=y;
   xminus=x;
   xplus=x;
%
%  Compute points at both sides of the switching surface
%
   ii=1;
   gaux=abs(g0)<=100*eps(1);
   gaux=gaux(ind);
   while any(gaux>=0)  && ii < 7,   % Should be ii < 4 for example ?
       yplus=y+10^ii*eps(1)*gxd;
       yminus=y-10^ii*eps(1)*gxd;
       xplus=x+10^ii*eps(1)*gt;
       xminus=x-10^ii*eps(1)*gt;
       [gminus,~,~]=switchfun(xminus,yminus);
       [gplus,~,~]=switchfun(xplus,yplus);
       gaux=gminus(ind).*gplus(ind);
       stats(10)=stats(10)+2;
       ii=ii+1;
   end  
   if ii >=7
       fprintf('\n Warning !!! The initial point  is a switching point and code could not determine its type \n\n');
   end 
   ii=0;
   while any(gaux>=0)  && ii < 7,           
       yplus=y+10^ii*eps(1)*f0;
       yminus=x-10^ii*eps(1)*f0;           
       xplus=x+10^ii*eps(1);
       xminus=x-10^ii*eps(1);          
       [gminus,~,~]=switchfun(xminus,yminus);           
       [gplus,~,~]=switchfun(xplus,yplus);        
       gaux=gminus(ind).*gplus(ind);         
       stats(10)=stats(10)+2;          
       ii=ii+1;     
   end
   if ii >=7
       fprintf('\n Warning !!! The initial point  is a switching point and code could not determine its type \n\n');
   end        
%        
% Checking if the discontinuity point is Filipov or transversal
%        
       
   f1=FUN(xminus,yminus);
   f2=FUN(xplus,yplus);    
   stats(9)=stats(9)+2;        
        gfxf1=gt+gxd'*f1;
        gfxf2=gt+gxd'*f2;
        if max(abs(gfxf1), abs(gfxf2)) < min(5.e-7, max(tol, 4000*eps(1))),
            fprintf('\n\n\n Warning! the vector field is tangent at both sides of the switching surface  at t= %g', x); 
            fprintf('\n The integration can not be reliable \n \n'); 
            disctype=3;
            endslid=max([abs(gfxf2) abs(gfxf2) 1.e-30]);
            xout=x;
            yout=y;
            if sign(gplus(ind(1)))==tangentside,
                ff=f2;
            else
                ff=f1;
            end
            idis=[idis -ind(1)];
            stats(7)=stats(7)+1;
        elseif gfxf1*gfxf2<0 && min(abs(gfxf2), abs(gfxf1)) > min([5.e-7 max(20*tol, 4000*eps(1) )]),
%        
% A Filipov point is detected
%
            disctype=3;            
            alfa=gfxf1/(gfxf1-gfxf2);
            ff=(1-alfa)*f1+alfa*f2;
            yout=y;
            xout=x;
            endslid=-min([abs(gfxf1)/norm(f1) abs(gfxf2)/norm(f2)]);
            idis=[idis -ind(1)];
            stats(7)=stats(7)+1;
        elseif  abs(gfxf2)<= min([5.e-7 max(20*tol, 4000*eps(1) )]),
%        
% A one side tangencial point is detected
%
            xout=xminus+max(5.e-7,10*eps(x));
            yout=yminus+max(5.e-7,10*eps(x))*f1;
            [g1,~,~]=switchfun(xout,yout);
            stats(10)=stats(10)+1;
            if sign(gminus(ind(1)))==sign(g1(ind(1))),
                xout=x+min(1.e-9, max(tol/100,eps(x)));
                yout=y+min(1.e-9, max(tol/100,eps(x)))*f1;
                ff=f1;
                disctype=1;
                idis=[idis ind(1)];
                stats(6)=stats(6)+1;
            else
                xout=xplus;
                ff=f2;
                disctype=3;
                endslid=max([min(abs(gfxf2), abs(gfxf1)),  1.e-30]);
                idis=[idis -ind(1)];
                stats(7)=stats(7)+1;
            end
       elseif  abs(gfxf1)<= min([5.e-7 max(20*tol, 5.e-13)]),
%        
% A one side tangencial point is detected
%
            xout=xplus+max(5.e-7,10*eps(x));
            yout=yplus+max(5.e-7,10*eps(x))*f2;
            [g1,~,~]=switchfun(xout,yout);
            stats(10)=stats(10)+1;
            if sign(gplus(ind(1)))==sign(g1(ind(1))),
                xout=x+min(1.e-9, max(tol/100,eps(x)));
                yout=y+min(1.e-9, max(tol/100,eps(x)))*f2;
                ff=f2; 
                disctype=1;
                idis=[idis ind(1)];
                stats(6)=stats(6)+1;
            else
                xout=xminus;
                ff=f1;
                disctype=3;
                endslid=max([min(abs(gfxf2), abs(gfxf1)),  1.e-30]);
                idis=[idis -ind(1)];
                stats(7)=stats(7)+1;
            end
        else
%
%  Transversal discontinuity.  The integration proceeds from the point
%  just after it, (xplus,yplus)
%
            xout=xminus+max(5.e-7,10*eps(x));
            yout=yminus+max(5.e-7,10*eps(x))*f1;
            [g1,~,~]=switchfun(xout,yout);
            stats(10)=stats(10)+1;
            if sign(gminus(ind(1)))==sign(g1(ind(1))),
                xout=x+min(1.e-9, max(tol/100,eps(x)));
                yout=y+min(1.e-9, max(tol/100,eps(x)))*f1;
            else
                xout=x+min(1.e-9, max(tol/100,eps(x)));
                yout=y+min(1.e-9, max(tol/100,eps(x)))*f2;
            end
           ff=FUN(xout,yout);
           stats(9)=stats(9)+1;
           disctype=1;
           endslid=gfxf1*gfxf2;
           idis=[idis ind(1)];
           stats(6)=stats(6)+1;
        end
  else
       f1=FUN(xminus,yminus);
       f2=FUN(xplus,yplus);    
       stats(9)=stats(9)+2;
       if exactgradient
           gt=gradt(switchfun,x,y,inddis,tol); 
           gfx=feval(gradswitchfun,x,y,inddis);
           stats(11)=stats(11)+ 1;
           gfxf1=(gt+gfx'*f1)/norm(gfx);
           gfxf2=(gt+gfx'*f2)/norm(gfx);
        else
           gt=gradt(switchfun,x,y,inddis,tol); 
           gfxf1=gt+graddif(switchfun,x,y,inddis,tol,f1);
           gfxf2=gt+graddif(switchfun,x,y,inddis,tol,f2); 
           if tol > 1.e-9
               stats(10)=stats(10)+6;
           else
               stats(10)=stats(10)+12;
           end

       end 
       if gfxf1*gfxf2<0 && (abs(gfxf2)/norm(f2))/max(1, abs(gfxf1)/norm(f1)) > min([5.e-7 max(20*tol, 4000*eps(1))]),
%        
% A Filipov point is detected
%
            idis=[idis -inddis];
            stats(7)=stats(7)+1;
            disctype=3;            
            alfa=gfxf1/(gfxf1-gfxf2);
            ff=(1-alfa)*f1+alfa*f2;
            yout=y;
            xout=x;
            endslid=-min([abs(gfxf1)/norm(f1) abs(gfxf2)/norm(f2)]);
        elseif  (abs(gfxf2)/norm(f2))/max(1, abs(gfxf1)/norm(f1))<= min([5.e-7 max(20*tol, 4000*eps(1))]),
%        
% A tangencial point is detected
%
            idis=[idis -inddis];
            stats(7)=stats(7)+1;
            yout=yplus;
            xout=xplus;
            ff=f2;
            disctype=3;
            endslid=max([abs(gfxf2)/norm(f2) 1.e-30]);
            if (abs(gfxf1)/norm(f1))/max(1, abs(gfxf2)/norm(f2))<= min([5.e-7 max(20*tol, 4000*eps(1))]),
                 fprintf('\n\n\n Warning! the vector field is tangent at both sides of the switching surface '); 
                 fprintf('\n The integration can not be reliable \n \n'); 
            end
        else
%
%  Transversal discontinuity.  The integration proceeds from the point
%  just after it, (xplus,yplus)
%
           idis=[idis inddis];
           stats(6)=stats(6)+1;
           xout=xplus;
           yout=yplus;
           ff=f2;
           disctype=1;
           endslid=gfxf1*gfxf2;
        end
   end
end
 

        

