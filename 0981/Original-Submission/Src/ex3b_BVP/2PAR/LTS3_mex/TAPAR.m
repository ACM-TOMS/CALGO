function [CONLAM, CONSIG, CONNU, NOPTS] = TAPAR (sigma0, TVALUE, tol, NOSING, SINGS, MULT)

%{
***********************************************************************
   TALBOT'S METHOD IMPLEMENTATION FOR THE NUMERICAL INVERSION OF      *
                         LAPLACE TRANSFORMS                           *
                             (MODULE 1)                               *
                                                                      *
>>>>>>>>>>>>>        VERSION 3.0    FEBRUARY 1, 2014        <<<<<<<<<<*
                                                                      *
                     AUTHOR: MARIAROSARIA RIZZARDI                    *
                                                                      *
               DEPARTMENT OF MATHEMATIC AND APPLICATIONS              *
                      UNIVERSITY OF NAPLES - ITALY                    *
***********************************************************************
   REFERENCES
   ==========
   MURLI A., M. RIZZARDI   -   "ALGORITHM XXX: TALBOT'S METHOD FOR THE
                                LAPLACE INVERSION PROBLEM".
                                ACM TRANS. MATH. SOFTWARE, VOL. ##,
                                NO. #, MONTH YEAR, PP. ##-##.


   PURPOSE
   =======
   THIS SUBROUTINE PROVIDES VALUES OF THE CONTOUR PARAMETERS (LAMBDA, SIGMA,
   NU) AND OF THE ACCURACY PARAMETER (N) ACCORDING TO THE TALBOT METHOD
   FOR THE NUMERICAL INVERSION OF LAPLACE TRANSFORMS.
   AFTER THIS ROUTINE HAS RETURNED ITS RESULTS, THESE WOULD NORMALLY
   BE USED IN SUBROUTINE TSUM (MODULE 2 OF TALBOT'S METHOD IMPLEMENTA-
   TION), WHICH CALCULATES THE INVERSE LAPLACE TRANSFORM     F(TVALUE)
   BY A CONTOUR INTEGRATION.

   CALLING SEQUENCE
   ================
     [CONLAM, CONSIG, CONNU, NOPTS] = TAPAR (sigma0, TVALUE, tol, NOSING,
                                                              SINGS, MULT)

   INPUT PARAMETERS
   ================
   sigma0         - (DOUBLE) THE ABSCISSA OF CONVERGENCE. 
   TVALUE         - (DOUBLE) THE VALUE IN WHICH THE INVERSE LAPLACE
                    TRANSFORM WILL BE COMPUTED.
                    TVALUE MUST BE A POSITIVE NUMBER.

   tol            - (DOUBLE) A TOLERANCE TO THE ABSOLUTE ACCURACY REQUIRED
                    IN THE RESULT.

   NOSING         - (INTEGER) THE NUMBER OF SINGULARITIES OF THE
                    LAPLACE TRANSFORM [NOSING = numel(SINGS)]. THE USER HAS
                    TO PROVIDE ONLY SINGULARITIES WITH NON NEGATIVE
                    IMAGINARY PARTS.

   SINGS         - (COMPLEX ARRAY) THIS ARRAY CONTAINS THE SINGULARITIES OF
                    THE LAPLACE TRANSFORM. IT MUST BE DIMENSIONED AT LEAST
                    "NOSING".

   MULT         - (INTEGER ARRAY) THIS ARRAY CONTAINS THE MULTIPLICITIES OF
                    THOSE SINGULARITIES WHICH ARE POLES, ZERO OTHERWISE.


   OUTPUT PARAMETERS
   =================
   CONLAM         - (DOUBLE)
   CONSIG         -   "
   CONNU          -   "    GEOMETRICAL PARAMETERS FOR THE TALBOT INTE_
                    GRATION CONTOUR (RESPECTIVELY LAMBDA, SIGMA AND NU)
   NOPTS          - (INTEGER) THE NUMBER OF POINTS WHICH WOULD BE
                    REQUIRED BY THE QUADRATURE ROUTINE TO OBTAIN THE
                    ACCURACY INDICATED BY "DECDIG".

***********************************************************************
%}
    IC = floor( -log10(eps/2)*0.75 );
    CASE1=true;
    OMEGA=0.4*(IC+1);

	PMAX = max(real(SINGS));    % max of real part of singularities
	[~,KMAX] = max(MULT);       % KMAX: index of max multiplicity
    if (PMAX > sigma0)
        sigma0=PMAX;
    end

% THE FOLLOWING STATEMENTS COMPUTE THE "DOMINANT SINGULARITY"
	RD=0.; ID=0;
	SID=0.;
	TETD=pi;
	for J=1:NOSING
        if imag(SINGS(J)) > 0.0
            CASE1=false;
            ETAJ=atan2(-real(SINGS(J))+sigma0,imag(SINGS(J)));
            TETJ=ETAJ + pi/2;
            RJ=imag(SINGS(J))/TETJ;
            if RD < RJ
                RD=RJ;
                TETD=TETJ;
                ID=J;
            end
        end
    end

	if CASE1 % (CASE1 == true)
        CONLAM=OMEGA/TVALUE;
        CONSIG=sigma0;
        CONNU=1;
    else     % (CASE1 == false)
        SRD=real(SINGS(ID)) - sigma0;
        SID=imag(SINGS(ID));
        MD =MULT(ID);
        V  =SID*TVALUE;
        OMEGA=min([OMEGA + V/2, 2*(IC+1)/3]);
        if 1.8*V <= OMEGA*TETD
            % GEOMETRICAL PARAMETERS' CHOICE   FOR   "CASE-1"
            CASE1=true;
            CONLAM=OMEGA/TVALUE;
            CONSIG=sigma0;
            CONNU=1;
        else
            % GEOMETRICAL PARAMETERS' CHOICE   FOR   "CASE-2"
            CASE1=false;
            PK=1.6 + 12/(V+25);
            FI=1.05 + 1050/max([553 800-V]);
            PMU=(OMEGA/TVALUE + sigma0 - PMAX)/(PK/FI - 1/tan(FI));
            CONLAM=PK*PMU/FI;
            CONSIG=PMAX-PMU/tan(FI);
            CONNU=SID/PMU;
        end
	end

    % required accurate decimal digits in the result
    DECDIG = ceil(-log10(tol));

% THE FOLLOWING STATEMENTS COMPUTE   NOPTS = MAX (N0, N1, N2)
	E=(2.3*DECDIG+OMEGA)/(CONLAM*TVALUE);
    if E <= 4.4
        UNRO=(16 + 4.3 * E)/(24.8 - 2.5 * E);
    elseif (E <= 10)
        UNRO=(50. + 3. * E)/(129. / E - 4.0);
    else % if E > 10
        UNRO=(44. + 19. * E)/(256. / E + 0.4);
    end
    N1=fix( CONLAM*TVALUE*(UNRO + (CONNU-1)/2) );
	N1=N1+1;
	N0=0;
    if SID ~= 0
        DD=DECDIG+2*min([MD-1 1]) + fix(MD/4);
        if CASE1 && MD ~= 0 % (CASE1 == TRUE && (MD != 0.0))
            P=SRD/CONLAM;
            Q=SID/CONLAM;
            U = INVPC(P,Q,TETD);
            if U > 0
                N0=fix( (2.3*DD+SRD*TVALUE)/U );
                N0=N0+1;
            end
        end
        GAMM=(CONSIG-sigma0)/CONLAM;
        Y=V/1000;
        ETA=min([1.78, 1.236 + 0.0064*1.78^DD]);
        ETA=ETA*(1.09-Y*(0.92-0.8*Y));
        N2=fix( ETA*CONNU*(2.3*DD+OMEGA)/(3+4*GAMM+1/exp(GAMM)) );
        N2=N2+1;
    else
        DD=DECDIG+2*min([MULT(KMAX)-1, 1]) + fix(MULT(KMAX)/4);
        N2=fix( (2.3*DD+OMEGA)/2 );
        N2=N2+1;
        ICM1=IC-1;
        if DD >= ICM1
            Q=0.0;
            for L=1:NOSING
                if MULT(L) ~= 0
                    DD=DECDIG+2*min([MULT(L)-1, 1]) + fix(MULT(L)/4);
                    SRD=real(SINGS(L)) - sigma0;
                    if  DD >= ICM1 && SRD < 0
                        P=SRD/CONLAM;
                        U = INVPC(P,Q,TETD);
                        if U > 0
                            N0L=fix( (2.3*DD+TVALUE*SRD)/U );
                            if N0L > N0
                                N0=N0L;
                            end
                        end
                    end
                end
            end
            N0=N0+1;
        end % if (DD >= ICM1)
    end

    NOPTS = max([N0, N1, N2]);
end % end tapar function



function U = INVPC (P,Q,TETA)
%{
    INTERNAL FUNCTION
***********************************************************************

   PURPOSE
   =======
   FIND THE  PRINCIPAL INVERSE  Z=-U + I*Y  OF  S=P + I*Q=R*EXP(I*TETA)
   IN                  S=Z/(1 - EXP(-Z))
   APPLYING NEWTON'S METHOD FOR REAL ROOTS


   INPUT PARAMETERS
   ================
   P             - (DOUBLE)
   Q             -    "
   TETA          -    "    RESPECTIVELY REAL PART, IMAGINARY PART AND
                   ARGUMENT OF S.
                   IT IS SUPPOSED THAT
                                  P .LE. 0.0
                                  Q .GE. 0.0
                     pi/2 .LE. TETA .LE. pi


   OUTPUT PARAMETERS
   =================
   U             - (DOUBLE) AN APPROXIMATION TO THE REAL PART (SIGN INVER-
                   TED) OF THE  PRINCIPAL INVERSE  OF S.

***********************************************************************
%}
    R=abs(complex(P,Q)); % R=sqrt(P^2+Q^2);

%   SET THE INITIAL GUESS
    X=13/(5.-2*P-Q-0.45*exp(P));
    Y=2*pi - X;

%   START ITERATIONS
    K=1; DY=1;
%       TEST FOR NUMBER OF ITERATION (MAXIMUM IS 100)
%       TEST FOR CONVERGENCE (tolerance is 1.e-4)
    while  K <= 100 && abs(DY) >= 1e-4 
        ANG=Y-TETA;
        S=sin(ANG);
        C=cos(ANG);
        B=Q-Y;
        E=-B/S;
        ARG=E/R;
        if ARG <= 0.0, U = -1; return; end
        U = log(ARG);
        G=P+B*C/S+U;
        DY=B*G/(1+E*(E-2*C));
        K=K+1;
        Y=Y+DY;
    end
end