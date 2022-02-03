% MAIN.m
%{
/********************              MAIN.m            **********************
 *                                                                        *
 *                                                                        *
 *                      DRIVER PROGRAM FOR Duffy's EXAMPLE                *
 *                             ACCURACY TEST FOR                          *
 *                                                                        *
 *                          SEQUENTIAL VERSIONS OF                        *
 *                    TALBOT SUITE and TALBOT SUITE DE                    *
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

%% check
if isempty( dir('num_sol.mex*') )
    error('No mex executable file. Run before:  "runme"  script for accuracy test.')
end


%% MENU: 'Talbot Suite' or 'Talbot Suite DE'
DEmenu = menu('?','Talbot Suite','Talbot Suite DE');
switch DEmenu
    case 1 % Talbot Suite
        DEstr = '';
    case 2 % Talbot Suite DE
        DEstr = ' DE';
end

%% UNCOMMENT TO SAVE OUTPUT TO A FILE
%{
    FILEname = ['out_accuracy_Tsuite' DEstr '_tol' num2str(TOLmenu) '.txt'];
    diary(FILEname)
%}


%% APPEND the "MinGW-w64\...\bin" folder to your system PATH variable in Windows
append_MinGW_dir

%% DOMAIN (X,t)=(eta,t)
Xmin=0;   Xmax=1;
Tmin=0.5;    Tmax=20;

 NXval=9; NTval=5;
%NXval=19; NTval=21;


%% SINGULARITIES: double poles at +/-3i, branch point at -1
NOSING=2;   sigma0=0;
SINGS=[-1, 3i];    MULT=[0 2];    
problem_string


%% tolerance
TOLmenu = 2; %TOLmenu = menu('tol = ?','1e-6','1e-8','1e-10','1e-12');
switch TOLmenu
    case 1, tol = 1e-6;
    case 2, tol = 1e-8;
    case 3, tol = 1e-10;
    case 4, tol = 1e-12;
end

fprintf('\nNTval=%d;  Tmin=%g;  Tmax=%g;  NXval=%d;  Xmin=%g;  Xmax=%g;  tol=%e\n', NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);

% modified Talbot's method
Tmid = mean([Tmin Tmax]);
[CONLAM, CONSIG, CONNU, NOPTS] = TAPAR (sigma0, Tmid, tol, NOSING, SINGS, MULT);
if DEmenu == 2 % for Talbot Suite DE
    NOPTS1 = Ncorrection(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4); % with N correction
    ZLABEL = 'with N correction';
    fprintf('Apply correction to NOPTS1:\n')
else           % for Talbot Suite
    NOPTS1 = NOPTS;                                                        % without N correction
    ZLABEL = 'without N correction';
end

% classical Talbot's method
NOPTS2=zeros(1,2);
[~, ~, ~, NOPTS2(1)] = TAPAR (sigma0, Tmin, tol, NOSING, SINGS, MULT);
[~, ~, ~, NOPTS2(2)] = TAPAR (sigma0, Tmax, tol, NOSING, SINGS, MULT);

disp(['NOPTS1 = ' num2str(NOPTS1) '  (modified Talbot''s method)'])
disp(['NOPTS2 = [' num2str(NOPTS2) ']  (classical Talbot''s method)'])
fprintf('\n')

%% COMPUTE SOLUTIONS
switch DEmenu
    case 1 % Talbot Suite
        u1 = num_sol (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
        u2 = num_sol (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    case 2 % Talbot Suite DE
        u1 = num_sol (3,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
        u2 = num_sol (4,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
end

%% COMPUTE ERRORS
ABSERR = true; % absolute (true) or relative (false) error

ERROR1=zeros(size(u1));
ERROR2=zeros(size(u2));
[t,n]=meshgrid(linspace(Tmin,Tmax,NTval),linspace(Xmin,Xmax,NXval));
for h=1:NXval
    for k=1:NTval
        ft(h,k) = ILT_fun(n(h,k),t(h,k));
        ERROR1(h,k) = abs(u1(h,k) - ft(h,k));
        ERROR2(h,k) = abs(u2(h,k) - ft(h,k));
    end
end

%% DISPLAY ERROR MATRICES
fprintf('    ***   Talbot Suite%s   ***\n', DEstr)
display_errors
diary off

