%  k-Method

% Script Experiments with HDG3d for Reaction-Diffusion problems
 
clear 

n=input('n: ');
switch n
    case 1
        example=4;
        ctmass=1;
        domain=1;
        what=2;
    case 0
        example=input('Input example (1->P1, 2->P2, 3->P3, 4->Cinfty): ');
        ctmass=input('Coefficients: (1) variable; (0) constant: ');
        domain=input('Domain: (1) chimney (2) cube (3) corner: ');
        what=input('Which triangulation: ');
end

switch domain
    case 1
        load 4Tchimney     % Chimney-shaped domain - Dirichlet = horizontal
        list={T1,T2,T3,T4};
    case 2
        load sixT3dDir     % (0,1)^3 all Dirichlet
        list={T1,T2,T3,T4,T5,T6};
    case 3
        load Corner        % Fichera corner (uniform partition) - mixed BC
        list={T1,T2,T3,T4};
end
T=list{what};

switch ctmass
    case 1
        kappa = @(x,y,z) 2+sin(x).*sin(y).*sin(z);
        kx = @(x,y,z) cos(x).*sin(y).*sin(z);
        ky = @(x,y,z) sin(x).*cos(y).*sin(z);
        kz = @(x,y,z) sin(x).*sin(y).*cos(z);
        c = @(x,y,z) 1.+0.5*(x.^2+y.^2+z.^2);
    case 0
        kappa = @(x,y,z) 2+0.*x;
        kx = @(x,y,z) 0.*x;
        ky = @(x,y,z) 0.*x;
        kz = @(x,y,z) 0.*x;
        c = @(x,y,z) 1+0.*x;
end

switch example
    case 1
        u = @(x,y,z) 2*x+3*y+z;
        ux = @(x,y,z) 2+0.*x;
        uy = @(x,y,z) 3+0.*x;
        uz = @(x,y,z) 1+0.*x;
        uxx = @(x,y,z) 0.*x;
        uyy = @(x,y,z) 0.*x;
        uzz = @(x,y,z) 0.*x;
    case 2
        u = @(x,y,z) x.^2+y.*x+z.*y+2*z.^2;
        ux = @(x,y,z) 2*x + y;
        uy = @(x,y,z) x + z;
        uz = @(x,y,z) y + 4*z;
        uxx = @(x,y,z) 2+0.*x;
        uyy = @(x,y,z) 0.*x;
        uzz = @(x,y,z) 4+0.*x;
    case 3
        u = @(x,y,z) 2*x.^2.*y+4*y.^2.*z+3*x.*z.^2;
        ux = @(x,y,z) 3*z.^2 + 4*x.*y;
        uy = @(x,y,z) 2*x.^2 + 8*y.*z;
        uz = @(x,y,z) 4*y.^2 + 6*x.*z;
        uxx = @(x,y,z) 4*y;
        uyy = @(x,y,z) 8*z;
        uzz = @(x,y,z) 6*x;
    case 4
        u = @(x,y,z) sin(x.*y.*z);
        ux = @(x,y,z) y.*z.*cos(x.*y.*z);
        uy = @(x,y,z) x.*z.*cos(x.*y.*z);
        uz = @(x,y,z) x.*y.*cos(x.*y.*z);
        uxx = @(x,y,z) -y.^2.*z.^2.*sin(x.*y.*z);
        uyy = @(x,y,z) -x.^2.*z.^2.*sin(x.*y.*z);
        uzz = @(x,y,z) -x.^2.*y.^2.*sin(x.*y.*z);
end

uD = u;
km = @(x,y,z) 1./kappa(x,y,z);
qx = @(x,y,z) -kappa(x,y,z).*ux(x,y,z);
qy = @(x,y,z) -kappa(x,y,z).*uy(x,y,z);
qz = @(x,y,z) -kappa(x,y,z).*uz(x,y,z);

f = @(x,y,z) -(kx(x,y,z).*ux(x,y,z)+kappa(x,y,z).*uxx(x,y,z))...
             -(ky(x,y,z).*uy(x,y,z)+kappa(x,y,z).*uyy(x,y,z))...
             -(kz(x,y,z).*uz(x,y,z)+kappa(x,y,z).*uzz(x,y,z))...
             +c(x,y,z).*u(x,y,z);

gx = @(x,y,z) -qx(x,y,z);         
gy = @(x,y,z) -qy(x,y,z);  
gz = @(x,y,z) -qz(x,y,z); 

a = @(x,y,z) 0.*x;

ErrorU=[];
ErrorQ=[];
ErrorUhat=[];
ErrorPu=[];
ErrorPuhat=[];
ErrorUstar=[];
h=[];

% Norms of unknowns for relative error


Nelts=size(T.elements,1);
formulas=checkQuadrature3d(0,ctmass);
normU=errorElem(T,u,zeros(1,Nelts),0,formulas{1});
normQ=errorElem(T,qx,zeros(1,Nelts),0,formulas{1})...
        +errorElem(T,qy,zeros(1,Nelts),0,formulas{1})...
        +errorElem(T,qz,zeros(1,Nelts),0,formulas{1});
    
for k=0:3
    Nelts=size(T.elements,1);
    tau=createTau3d(Nelts,2); 
    formulas=checkQuadrature3d(k,ctmass);
    [Uh,Qxh,Qyh,Qzh,Uhat]=HDG3d(km,c,f,tau,T,k,formulas,uD,gx,gy,gz);
    normUhat=errorFaces(T,u,zeros(size(Uhat)),k,formulas{4});
    
    % Errors with exact solution
    error_uhat=errorFaces(T,u,Uhat,k,formulas{4});
    error_q   =errorElem(T,qx,Qxh,k,formulas{1})...
               +errorElem(T,qy,Qyh,k,formulas{1})...
               +errorElem(T,qz,Qzh,k,formulas{1});
    error_u   =errorElem(T,u,Uh,k,formulas{1});
    ErrorQ=[ErrorQ error_q/normQ]; 
    ErrorU=[ErrorU error_u/normU];
    ErrorUhat=[ErrorUhat error_uhat/normUhat];
    
    % Errors with projections
    Puhat=L2projskeleton3d(u,T,k,formulas{3},formulas{4});
    [Pqx,Pqy,Pqz,Pu]=projectHDG3d(T,{qx,qy,qz,u},k,tau,formulas);
    error_Puhat=errorFaces(T,a,Uhat-Puhat,k,formulas{4});
    error_Pu=errorElem(T,a,Pu-Uh,k,formulas{1});
    ErrorPu=[ErrorPu error_Pu/normU];
    ErrorPuhat=[ErrorPuhat error_Puhat/normUhat];
    
    % Postprocessing
    Uhstar=postprocessing(T,km,Qxh,Qyh,Qzh,Uh,k,formulas{1});
    error_Ustar=errorElem(T,u,Uhstar,k+1,formulas{1});
    ErrorUstar=[ErrorUstar error_Ustar/normU];
end

format shortE
disp('Error |Q-Qh|   Error |U-Uh|     Error |U-Uhat|_h');
disp([ErrorQ' ErrorU' ErrorUhat']);
disp('Error |\Pi u-Uh|  Error |Pu-Uhat|_h Error |U-U*|')
disp([ErrorPu' ErrorPuhat' ErrorUstar']);

logErrorQ=log(ErrorQ);
logErrorU=log(ErrorU);
logErrorUhat=log(ErrorUhat);
logErrorPu=log(ErrorPu);
logErrorPuhat=log(ErrorPuhat);
logErrorUstar=log(ErrorUstar);

k=1:k-1;kp1=k+1;kp2=k+2;
rateQ=(logErrorQ(k)-logErrorQ(kp1))./(logErrorQ(kp1)-logErrorQ(kp2));
rateU=(logErrorU(k)-logErrorU(kp1))./(logErrorU(kp1)-logErrorU(kp2));
rateUhat=(logErrorUhat(k)-logErrorUhat(kp1))./(logErrorUhat(kp1)-logErrorUhat(kp2));
ratePu=(logErrorPu(k)-logErrorPu(kp1))./(logErrorPu(kp1)-logErrorPu(kp2));
ratePuhat=(logErrorPuhat(k)-logErrorPuhat(kp1))./(logErrorPuhat(kp1)-logErrorPuhat(kp2));
rateUstar=(logErrorUstar(k)-logErrorUstar(kp1))./(logErrorUstar(kp1)-logErrorUstar(kp2));

format bank
disp('rateQ    rateU    rateUhat')
disp([rateQ' rateU' rateUhat'])
disp('ratePu ratePuhat  ratePustar')
disp([ratePu' ratePuhat' rateUstar'])