% MAIN.m
%{
/********************              MAIN.m            **********************
 *                                                                        *
 *                                                                        *
 *                      DRIVER PROGRAM FOR Duffy's EXAMPLE                *
 *                             ACCURACY TEST FOR                          *
 *                                                                        *
 *               OpenMP PARALLEL VERSIONS OF TALBOT SUITE DE              *
 *                                                                        *
 *                                                                        *
 *                       AUTHOR: Mariarosaria Rizzardi                    *
 *                                                                        *
 *                  mariarosaria.rizzardi@uniparthenope.it                *
 *                                                                        *
 *                  DiST - Dept. of Science and Technology                *
 *                  "Parthenope" University, Naples (Italy)               *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * M. RIZZARDI: "Algorithm xxx: Talbot Suite DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP ##).            *
 *                                                                        *
 **************************************************************************/
%}

path('../../',path)


%% APPEND the "MinGW-w64\...\bin" folder to your system PATH variable in Windows
append_MinGW_dir


%% MENUS u(eta,t)
Xmin=0;   Xmax=1;
Tmin=0.5;    Tmax=20;

TOLmenu = 2; %TOLmenu = menu('tol = ?','1e-6','1e-8','1e-10','1e-12');
switch TOLmenu
    case 1, tol = 1e-6;
    case 2, tol = 1e-8;
    case 3, tol = 1e-10;
    case 4, tol = 1e-12;
end


%% SINGULARITIES
NOSING=2;   sigma0=0;
SINGS=[-1, 3i];    MULT=[0 2];    
problem_string

%% ACCURACY
if isempty( dir('num_sol.mex*') )
    error('No mex executable file. Run before:  "runme"  script for accuracy test.')
end

 NXval=9; NTval=5;
%NXval=19; NTval=21;

fprintf('\nNTval=%d;  Tmin=%g;  Tmax=%g;  NXval=%d;  Xmin=%g;  Xmax=%g;  tol=%e\n', NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);


%% COMPUTE SOLUTIONS
DEstr = '\_DE';
u1 = num_sol (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); % OMP_Talbot11_DE
u2 = num_sol (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); % OMP_Talbot12_DE
u3 = num_sol (3,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); % OMP_Talbot13_DE

%% COMPUTE ERRORS
ABSERR = true; % absolute or relative error
[t,n]=meshgrid(linspace(Tmin,Tmax,NTval),linspace(Xmin,Xmax,NXval));
ft = zeros(NXval,NTval);
for h=1:NXval
    for k=1:NTval
        ft(h,k) = ILT_fun(n(h,k),t(h,k));
    end
end
ERROR1 = abs(u1 - ft);
ERROR2 = abs(u2 - ft);
ERROR3 = abs(u3 - ft);

%% DISPLAY ERROR MATRICES
display_errors
