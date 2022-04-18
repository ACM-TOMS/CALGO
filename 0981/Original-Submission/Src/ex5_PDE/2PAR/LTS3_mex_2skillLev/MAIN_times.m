% MAIN_times.m
%{
/********************           MAIN_times.m         **********************
 *                                                                        *
 *                                                                        *
 *                      PARALLEL OMP TALBOT SUITE DE                      *
 *                                                                        *
 *                       DRIVER PROGRAM FOR EXAMPLE 5                     *
 *                          SKILL-LEVEL FUNCTION                          *
 *                                 TIME TEST                              *
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

%% PROBLEM TO SOLVE
[~,Nsings,SINGS,MULT,sigma0] = LTsings2();
problem_string

%% SQUARE SPATIAL DOMAIN
XYmin=0;    XYmax=1;
Tmin=100;   Tmax=500;

%% MENUS
TOLmenu = 4; %TOLmenu = menu('tol = ?','1e-6','1e-8','1e-10','1e-12');
switch TOLmenu
    case 1, tol = 1e-6;
    case 2, tol = 1e-8;
    case 3, tol = 1e-10;
    case 4, tol = 1e-12;
end

%%
%      |  <----- NXYval ----->
% NTval|   5       20      120
%      |----------------------
%   5  |  25      100      600
%  20  | 100      400     2400
% 120  | 600     2400    14400
%
XTmenu = 9;
% XTmenu = menu('CHOOSE','NTval=5,   NXYval=5',   ... % 1 - spatial mesh size 3^2
%                        'NTval=5,   NXYval=20',  ... % 2 - spatial mesh size 18^2
%                        'NTval=5,   NXYval=120', ... % 3 - spatial mesh size 118^2
%                        'NTval=20,  NXYval=5',   ... % 4
%                        'NTval=20,  NXYval=20',  ... % 5
%                        'NTval=20,  NXYval=120', ... % 6
%                        'NTval=120, NXYval=5',   ... % 7
%                        'NTval=120, NXYval=20',  ... % 8
%                        'NTval=120, NXYval=120');    % 9

switch XTmenu
    case 1, NTval=5;    NXYval=5;
    case 2, NTval=5;    NXYval=20;
    case 3, NTval=5;    NXYval=120;
    case 4, NTval=20;   NXYval=5;
    case 5, NTval=20;   NXYval=20;
    case 6, NTval=20;   NXYval=120;
    case 7, NTval=120;  NXYval=5;
    case 8, NTval=120;  NXYval=20;
    case 9, NTval=120;  NXYval=120;
end

%%
[CONLAM,CONSIG,CONNU,NOPTS] = TAPAR (sigma0, mean([Tmin Tmax]), tol, Nsings, SINGS, MULT);
NOPTS = Ncorrection(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4); % correction to NOPTS
Tval = linspace(Tmin,Tmax,NTval)';
fprintf('\n\t   Ex. 5: output from ./2PAR/LTS3_mex_2skillLev/MAIN_times.m\n')
fprintf('\n\tLT samples by solving in MATLAB block tridiagonal systems + PCT\n')
fprintf('\n\t      NOPTS = %u\n', NOPTS)
fprintf('\nNTval = %d  t in [%g, %g],  NXYval = %d  x,y in [%g, %g],  tol = %e\n', NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol);
fprintf('====================================================================================\n\n')

%% START DEFAULT PARALLEL POOL
pool=parpool;
maxTHREADS = pool.NumWorkers; % maximum number of parallel workers

%% COMPUTE INTERNAL MESH POINTS
XY = internalMeshPts(NXYval, XYmin, XYmax);     % XY: col-wise matrix of size (Urows,2)
Urows = size(XY,1);

%% CALL Talbot Suite DE (MATLAB version for coarse-grain parallelism)
TIME = zeros(maxTHREADS,1);

% thrds: number of parallel threads and workers
ch='%';
fprintf('\n\t      Modified Talbot''s method [coarse grain parallelism]')

fprintf('\n\tELAPSED TIME:\n\n')
for thrds = 1:maxTHREADS
    tic
    [IFAIL_tot,NUMft,IFAIL] = OMP_Talbot11_DE (@LT_samples,sigma0,NXYval,XY,NTval,Tval,tol, ...
                                               Nsings,SINGS,MULT,Tmin,Tmax,thrds);
    TIME(thrds)=toc;
    fprintf('\t\t%e\t%c number of threads = %2d\n', TIME(thrds),ch,thrds)
end
fprintf('\n')

%% CLOSE PARALLEL POOL
delete(gcp);
