% MAIN_accuracy.m
%{
/********************         MAIN_accuracy.m        **********************
 *                                                                        *
 *                                                                        *
 *                      PARALLEL OMP TALBOT SUITE DE                      *
 *                                                                        *
 *                      DRIVER PROGRAM FOR EXAMPLE 5                      *
 *                          SKILL-LEVEL FUNCTION                          *
 *                               ACCURACY TEST                            *
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
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP ##).            *
 *                                                                        *
 **************************************************************************/
%}

THREADS = 4;   % number of parallel threads and workers


%% CHECK
if isempty( dir('OMP_Talbot11_SUM_mex.mex*') )
    error('No mex executable file. Run before:  "runme"  script for accuracy test.')
end

%% CHECK for Parallel Computing Toolbox
fun = 'parpool';      % ==>   distcomp
pat = '(?<=^.+[\\/]toolbox[\\/])[^\\/]+'; % match regular expression (case sensitive)
if isempty(regexp(which(fun), pat, 'match', 'once'))
    error('The Parallel Computing Toolbox is required to run this program')
else
    fprintf('\n*** The Parallel Computing Toolbox is installed. ***\n')
end

%% APPEND the "MinGW-w64\...\bin" folder to your system PATH variable in Windows
append_MinGW_dir

%%
path('../..',path);

%% 
TOLmenu=4; %TOLmenu = menu('tol = ?','1e-6','1e-8','1e-10','1e-12');
switch TOLmenu
    case 1, tol = 1e-6;
    case 2, tol = 1e-8;
    case 3, tol = 1e-10;
    case 4, tol = 1e-12;
end

%% PROBLEM TO SOLVE
[~,Nsings,SINGS,MULT,sigma0] = LTsings2();
problem_string

%% SQUARE SPATIAL DOMAIN
NXYval=9;   XYmin=0;    XYmax=1;
NTval=5;    Tmin=100;   Tmax=500;
Tval = linspace(Tmin,Tmax,NTval)';
fprintf('\n\t\tEx. 5: output from ./2PAR/LTS3_mex_2skillLev/MAIN_accuracy.m\n')
fprintf('\n\tLT samples by solving in MATLAB block tridiagonal systems + PCT\n')
fprintf('\nNTval = %d  t in [%g, %g],  NXYval = %d  x,y in [%g, %g],  tol = %e\n', NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol);
fprintf('====================================================================================\n\n')

%% START PARALLEL POOL
parpool(THREADS);

%% COMPUTE INTERNAL MESH POINTS
XY = internalMeshPts(NXYval, XYmin, XYmax);     % XY: col-wise matrix of size (Urows,2)
Urows = size(XY,1);

%% CALL Talbot Suite DE (MATLAB version for coarse-grain parallelism)
fprintf('\nmodified Talbot''s method   (Talbot Suite DE skill-level functions, coarse-grain parallelism)\n')
[IFAIL_tot,NUMft,IFAIL] = OMP_Talbot11_DE (@LT_samples,sigma0,NXYval,XY,NTval,Tval,tol, ...
                                          Nsings,SINGS,MULT,Tmin,Tmax,THREADS);

%% CLOSE PARALLEL POOL
delete(gcp);

%% COMPUTE AND DISPLAY ERRORS
u = ILTfun2(XY(:,1),XY(:,2),Tval);
RELERR1 = abs(u - NUMft)./abs(u);
fprintf('\nRelative Error in Inverse Laplace Transform u(h,k) = u(X(h),Y(h),t(k)):\n'); disp([(1:Urows)' RELERR1])
