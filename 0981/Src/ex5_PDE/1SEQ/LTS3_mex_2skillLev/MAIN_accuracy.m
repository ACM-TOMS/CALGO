% MAIN_accuracy.m
%{
/********************         MAIN_accuracy.m        **********************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
 *                                                                        *
 *                      DRIVER PROGRAM FOR EXAMPLE 5                     *
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

%% CHECK
if isempty( dir('SEQ_Talbot1_SUM_mex.mex*') )
    error('No mex executable file. Run before:  "runme"  script for accuracy test.')
end

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

%% SQUARE SPATIAL DOMAIN
problem_string

 NXYval=9;
 NTval =5;

XYmin=0;    XYmax=1;
Tmin=100;   Tmax=500;

Tval = linspace(Tmin,Tmax,NTval)';
fprintf('\nNTval = %d  t in [%g, %g],  NXYval = %d,  x,y in [%g, %g],  tol = %e\n', NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol);

%% COMPUTE POINTS ON TALBOT'S CONTOUR FOR MODIFIED METHOD
[~, Nsings, SINGS, MULT, sigma0] = LTsings2();
[CONLAM, CONSIG, CONNU, NOPTS] = TAPAR (sigma0, mean([Tmin Tmax]), tol, Nsings, SINGS, MULT);
NOPTS = Ncorrection(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);

thetak=linspace(0,pi,NOPTS+1)'; thetak(end)=[]; thetak(1)=[];
S = CONSIG + CONLAM*thetak./tan(thetak) + 1i*CONLAM*CONNU*thetak;
S = [CONSIG + CONLAM; S]; % ass the first for thetak(1)
clear thetak

%% COMPUTE THE LT SAMPLES
% U is the (complex-valued) matrix of LT samples. U is of size (Urows, NOPTS), where
%       Urows = (NXval-2)^2
%       NOPTS = numel(S).
% Each column K contains the solution U(x,y,S(K)) computed at internal mesh points (X,Y):
%           U(:,K) = U(X(:),Y(:),S(K))
% Each row contains the values to be used in the summation step: 
%           U(j,:) = U(X(j),Y(j),S(:))
% U is a MATLAB col-wise matrix.
XY = internalMeshPts(NXYval, XYmin, XYmax); % XY: col-wise matrix of size (Urows,2)
U = LT_samples (NXYval,XY,NOPTS,S,tol);     % U:  col-wise matrix of size (Urows,NOPTS)
Urows = size(U,1);

%% CALL THE C FUNCTION FOR THE TALBOT-CLENSHAW SUMMATION
% we need to transpose the matrix
U = U.';                                    % U:  row-wise matrix (for C language)
%1         2     3                             1      2      3     4     5     6 7     8
[IFAIL_tot,NUMft,IFAIL] = SEQ_Talbot1_SUM_mex (CONLAM,CONSIG,CONNU,NOPTS,Urows,U,NTval,Tval);
% NUMft: array 1D containing a col-wise matrix of size (Urows,NTval)
% IFAIL: array 1D containing a col-wise matrix of size (Urows,NTval)
NUMft = reshape(NUMft,Urows,NTval);
IFAIL = reshape(IFAIL,Urows,NTval);

%% COMPUTE AND DISPLAY ERRORS
u = ILTfun2(XY(:,1),XY(:,2),Tval);
RELERR1 = abs(u - NUMft)./abs(u);
fprintf('\n    num (x,y)   Relative Error in Inverse Laplace Transform u(h,k) = u(X(h),Y(h),t(k)):\n'); disp([(1:Urows)' RELERR1])

