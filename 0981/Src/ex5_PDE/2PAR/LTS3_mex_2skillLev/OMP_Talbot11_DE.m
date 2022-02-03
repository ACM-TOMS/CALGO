function [IFAIL_tot,NUMft,IFAIL] = OMP_Talbot11_DE (LTsamples,sigma0,NXYval,XY,NTval,Tval,tol, ...
                                          Nsings,SINGS,MULT,Tmin,Tmax,THREADS)
%{
/**************************************************************************
 *                                                                        *
 *           OpenMP-based parallel version of TALBOT SUITE DE             *
 *                                                                        *
 **************************************************************************
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
 *             OMP_Talbot11_DE:  DRIVER FUNCTION (user level)             *
 *                                                                        *
 *  IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS *
 *              COARSE-GRAINED PARALLELISM IMPLEMENTATION                 *
 *  MATLAB version of OMP_Talbot11_DE(); see its documentation in:        *
 *          TalbotSuiteDE_src/SRC/FUN_DE/OMP_Talbot_pack_DE.c             *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * M. RIZZARDI: "A modification of Talbot's method for the simultaneous   *
 *               approximation of several values of the Inverse Laplace   *
 *               Transform".                                              *
 *               ACM Trans. Math. Soft., vol. 21, no. 4, Dec. 1995,       *
 *               pp. 347-371.                                             *
 *                                                                        *
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: Application of Modified  *
 *               Talbot's Method to solve differential problems".         *
 *               ACM Trans. Math. Soft., vol. xx, no. x, month year,      *
 *               pp ##.                                                   *
 *                                                                        *
 **************************************************************************/
%}
	%% COMPUTE POINTS ON TALBOT'S CONTOUR FOR MODIFIED METHOD
	[CONLAM, CONSIG, CONNU, NOPTS] = TAPAR (sigma0, mean([Tmin Tmax]), tol, Nsings, SINGS, MULT);
    NOPTS = Ncorrection(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4); % correction to NOPTS
    thetak=linspace(0,pi,NOPTS+1)'; thetak(end)=[]; thetak(1)=[];
    S = CONSIG + CONLAM*thetak./tan(thetak) + 1i*CONLAM*CONNU*thetak;
    S = [CONSIG + CONLAM; S]; % ass the first for thetak(1)

    %% COMPUTE THE LT SAMPLES
    % U is the (complex-valued) matrix of LT samples. U is of size (Urows, NOPTS), where
    %       Urows = (NXval-2)^2
    %       NOPTS = numel(S).
    % Each column K contains the solution U(x,y,S(K)) computed at internal mesh points (X,Y):
    %           U(:,K) = U(X(:),Y(:),S(K))
    % Each row contains the values to be used in the summation step: 
    %           U(j,:) = U(X(j),Y(j),S(:))
    % U is a MATLAB col-wise matrix.
    U = LTsamples (NXYval,XY,NOPTS,S,tol,THREADS); % U:  col-wise matrix of size (Urows,NOPTS)
    Urows = size(U,1);

    %% CALL THE C FUNCTION FOR THE TALBOT-CLENSHAW SUMMATION (coarse-grain parallelism)
    % we need to transpose the matrix
    U = U.';                                    % U:  row-wise matrix (for C language)
    %1         2     3                              1      2      3     4     5     6 7     8
    [IFAIL_tot,NUMft,IFAIL] = OMP_Talbot11_SUM_mex (CONLAM,CONSIG,CONNU,NOPTS,Urows,U,NTval,Tval,THREADS);
    % NUMft: array 1D containing a col-wise matrix of size (Urows,NTval)
    % IFAIL: array 1D containing a col-wise matrix of size (Urows,NTval)
    NUMft = reshape(NUMft,Urows,NTval);
    IFAIL = reshape(IFAIL,Urows,NTval);
end

