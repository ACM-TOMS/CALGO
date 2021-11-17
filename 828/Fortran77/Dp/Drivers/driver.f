C
C
C ARRAY length and grid parameters:
C
C   LDA = Row dimension of array ABD used to store D**T*D
C         (and the matrix associated with the cubic spline
C         interpolant) in LAPACK band storage (packed)
C         format for the linear solver DPBSV (and CSPLIN).
C
C   MMAX = Maximum number of grid points.  The array storage
C          requirement is approximately 112*MMAX bytes for
C          8-byte double precision numbers.
C
C   NMAX = Maximum number of data points.
C
      INTEGER LDA, MMAX, NMAX
      PARAMETER (LDA=3, NMAX=100, MMAX=900)
C
C Arrays:
C
      DOUBLE PRECISION D(NMAX), S(NMAX), T(NMAX)
      DOUBLE PRECISION A(MMAX), ABD(LDA,MMAX), D1Y(MMAX),
     .                 D2Y(MMAX), G(MMAX), G0(MMAX),
     .                 P(MMAX), WK(MMAX,3), Y0(MMAX),
     .                 Y(MMAX)
      INTEGER IND(NMAX)
C
C Scalars and functions:
C
      CHARACTER CHR
      CHARACTER*60 FNAME
      DOUBLE PRECISION DDOT, PHI
      DOUBLE PRECISION ALPHA, C2, C3, DD, DG, DGG, DPHI, DT,
     .                 DY, FMTOL, GG, GNRM1, GNRMM, H, PHI0,
     .                 PHIY, PHTOL, PLTSIZ, Q, SSIZ0, SSIZ1,
     .                 YMAX
      INTEGER I, IER, J, J1, J2, LIN, LOUT, LPLT, LPRT, M,
     .        MAXCG, N, NCG, NG, NGMAX, NPHSM, NPHVAL, NPRNT
      LOGICAL L2G, OPT, PLT0, PRNT
C
C Parameters:
C
C   FMTOL = Nonnegative tolerance for the line search:
C           length of the interval of uncertainty.  Refer
C           to Function FMIN.
C
C   L2G = Flag with value TRUE iff the L2-gradient is to be
C         used in place of the Sobolev gradient.
C
C   MAXCG = Number of Polak-Ribiere conjugate gradient steps
C           between restarts with a steepest descent step.
C
C   NGMAX = Maximum number of descent steps (outer itera-
C           tions) before termination without convergence.
C           NGMAX > 0.
C
C   NPRNT = Number of descent steps between print steps.
C           For each print step, several lines of statistics
C           describing the step are written to unit LPRT.
C           The first and last descent steps are necessarily
C           print steps.  NPRNT .GE. 1.
C
C   OPT = Flag with value TRUE iff the optimal step-size
C         (computed by a line search) is to be used at each
C         descent step.
C
C   PHTOL = Nonnegative tolerance which defines convergence
C           of the descent method:  bound on Max(DPHI,DY**2,
C           DG**3), where
C             DPHI = abs(change in phi)/(1+phi),
C             DY = max-norm(change in Y)/(1+max-norm(Y)),
C             DG = max-norm(S-gradient)/(1+max-norm(Y)).
C
C   PLT0 = Flag with value TRUE iff the initial estimate is
C          to be plotted.
C
C   PLTSIZ = Plot size in inches for Postscript plots.
C            1.0 <= PLTSIZ <= 6.5.
C
C   SSIZ0 = Descent step size used if OPT = FALSE.
C
C   SSIZ1 = Initial descent step size used if OPT = TRUE.
C
      DATA FMTOL/1.D-6/,  L2G/.FALSE./,  MAXCG/3/,
     .     NGMAX/200/,    NPRNT/1/,      OPT/.TRUE./,
     .     PHTOL/1.D-3/,  PLT0/.TRUE./,  PLTSIZ/6.5/,
     .     SSIZ0/1.D-3/,  SSIZ1/1.D-5/
C
C Logical unit numbers:
C
C   LIN = Logical unit number for an input data set contain-
C         ing N, (T(I), I = 1,N), (D(I), I = 1,N), and H.
C
C   LOUT = Logical unit number for writing the solution
C          M, (Y(I), I = 1,M), (D1Y(I), I = 1,M), and
C          (D2Y(I), I = 1,M):  Subroutine WRITF.
C
C   LPLT = Logical unit number for files containing
C          Postscript plots:  initial estimate, spline
C          curve, and/or curvature.  Refer to Subroutine
C          PLTCRV.
C
C   LPRT = Logical unit number for output other than the
C          solution.
C
      DATA LIN/5/,  LOUT/6/,  LPLT/3/,  LPRT/4/
C
C Open print file, delete it, and reopen it.  This is
C   necessary on some systems (Microsoft) which ignore
C   end-of-file records.
C
      OPEN (LPRT,FILE='dnsplin1.prt')
      CLOSE (LPRT,STATUS='DELETE')
      OPEN (LPRT,FILE='dnsplin1.prt')
C
C Get an input data set file name.
C
C   1 WRITE (*,200)
C 200 FORMAT (//5X,'Specify the name of a file containing ',
C    .        'an input data set'/
C    .        5X,'(at most 60 characters):'/)
C     READ (*,205,ERR=1) FNAME
C 205 FORMAT (A60)
C     OPEN (LIN,FILE=FNAME,STATUS='OLD',ERR=1)
C
C Read input data, and compute indexes (1 to M=IND(N)) of
C   data points.
C
      CALL READF (LIN,NMAX, N,T,D,H,IND,IER)
      IF (IER .NE. 0) GO TO 21
      M = IND(N)
      IF (M .LT. 6  .OR.  M .GT. MMAX  .OR.  M-N .LT. 1)
     .  GO TO 22
C
C Print a heading with parameter values.
C
      WRITE (LPRT,100)
  100 FORMAT (///27X,'DNSPLIN1 Output'/)
      WRITE (LPRT,110) T(1), T(N), N, M, H, MAXCG, FMTOL,
     .                 PHTOL, NGMAX
  110 FORMAT (//5X,'Domain:  T(1) = ',D10.3,', T(N) = ',
     .             D10.3/
     .        5X,'Number of data points:  N = ',I4//
     .        5X,'Number of grid points:  M = ',I5/
     .        5X,'Mesh width:  H = ',D12.6//1P,
     .        5X,'# CG steps between restarts:  MAXCG = ',
     .           I4/
     .        5X,'Line search tolerance:  FMTOL = ',D9.3/
     .        5X,'Convergence tolerance:  PHTOL = ',D9.3/
     .        5X,'Max number of iterations:  NGMAX = ',I6)
C
C Set Y0 to the gridpoint values of the Hermite cubic
C   interpolant with derivatives S at the data abscissae
C   obtained from the natural cubic spline.
C
      CALL CSPLIN (LDA,N,T,D, ABD, S,IER)
      IF (IER .NE. 0) GO TO 23
C
      J2 = 1
      DO 3 I = 1,N-1
        J1 = J2
        J2 = IND(I+1)
        Y0(J1) = D(I)
        DT = T(I+1)-T(I)
        DD = (D(I+1)-D(I))/DT
        C2 = (3.D0*DD-2.D0*S(I)-S(I+1))/DT
        C3 = (-2.D0*DD+S(I)+S(I+1))/(DT*DT)
        DT = 0.D0
        DO 2 J = J1+1,J2-1
          DT = DT + H
          Y0(J) = ((C3*DT + C2)*DT + S(I))*DT + D(I)
    2     CONTINUE
    3   CONTINUE
      Y0(J2) = D(N)
C
C Copy Y0 into Y.
C
      DO 4 I = 1,M
        Y(I) = Y0(I)
    4   CONTINUE
C
C Prompt for optional plot of initial estimate.
C
C   5 WRITE (*,210)
C 210 FORMAT (//5X,'Create file dnsplin1.ps0 containing a ',
C    .        'plot '/
C    .        5X,'of the initial estimate? (y/n)'/)
C     READ (*,215,ERR=5) CHR
C 215 FORMAT (A1)
C     IF (CHR .NE. 'Y'  .AND.  CHR .NE. 'y'  .AND.
C    .    CHR .NE. 'N'  .AND.  CHR .NE. 'n') GO TO 5
C     IF (CHR .EQ. 'Y'  .OR.  CHR .EQ. 'y') THEN
        OPEN (LPLT,FILE='dnsplin1.ps0')
        CLOSE (LPLT,STATUS='DELETE')
        OPEN (LPLT,FILE='dnsplin1.ps0')
        CALL PLTCRV (LPLT,N,T(1),T(N),T,D,M,Y0,.FALSE.,Y0,
     .               PLTSIZ, IER)
        CLOSE (LPLT)
        IF (IER .NE. 0) GO TO 26
C     ENDIF
C
C Compute PHI0 = phi(Y0).  Grid-point values D1y(i), D2y(i),
C   and a(i) are also computed.
C
      PHI0 = PHI (M,H,Y0, D1Y,D2Y,A)
C
C Initialize iteration counts:
C
C   NCG = Number of conjugate gradient steps since the
C         previous restart with a steepest descent step.
C
C   NG = Number of iterations (gradient evaluations).
C
C   NPHSM = Total number of evaluations of phi in the line
C           searches.
C
C   PRNT = TRUE iff statistics describing the iteration are
C          to be printed.
C
      NCG = MAXCG
      NG = 0
      NPHSM = 0
      PRNT = .TRUE.
C
C Top of loop:  compute the L2-gradient G = grad(phi), its
C               1-norm GNRM1, and its max-norm GNRMM, and
C               overwrite G with the Sobolev gradient unless
C               L2G = TRUE.
C
    7 CALL GRADL2 (N,H,IND,M,D1Y,D2Y,A, WK, G,GNRM1,
     .             GNRMM,IER)
      IF (IER .NE. 0) GO TO 24
      IF (.NOT. L2G) THEN
        CALL GRADS4 (N,H,IND,M,D1Y,D2Y,LDA, ABD,G,
     .               WK, IER)
        IF (IER .GT. 0) WRITE (LPRT,125) IER
  125   FORMAT (//5X,'GRADS4 Failure:  IER = ',I1)
      ENDIF
C
C Compute mean absolute L2-gradient component GNRM1.
C
      GNRM1 = GNRM1/DBLE(M-N)
C
C Update iteration count.
C
      NG = NG + 1
C
C Print statistics associated with the step.
C
      IF (PRNT) THEN
        IF (NCG .GE. MAXCG) THEN
          WRITE (LPRT,130) NG
        ELSE
          WRITE (LPRT,135) NG
        ENDIF
        WRITE (LPRT,140) PHI0, GNRM1, GNRMM
      ENDIF
  130 FORMAT (//15X,'Iteration ',I6,' (Steepest descent)')
  135 FORMAT (//15X,'Iteration ',I6,' (Conjugate gradient)')
  140 FORMAT (5X,'Functional:                  phi(F) = ',
     .           1P,D21.15/
     .        5X,'Mean absolute L2-gradient:    GNRM1 = ',
     .           D8.2/
     .        5X,'Max-norm of L2-gradient:      GNRMM = ',
     .           D8.2)
      PRNT = (NPRNT*(NG/NPRNT) .EQ. NG)  .OR.  (NG .EQ. 1)
C
      IF (NCG .GE. MAXCG) THEN
C
C Steepest descent step.
C
        Q = 0.D0
        NCG = 0
      ELSE
C
C Conjugate gradient step.
C
        GG = DDOT(M,G0,1,G0,1)
        DGG = DDOT(M,G,1,G,1) - DDOT(M,G0,1,G,1)
        Q = DGG/GG
        NCG = NCG + 1
      ENDIF
C
C Compute the search direction p = -G + Q*p.
C
      DO 8 I = 1,M
        P(I) = -G(I) + Q*P(I)
    8   CONTINUE
C
C Copy Y into Y1 (WK column 1), compute the optimal step-
C   size ALPHA, and update Y to Y1+ALPHA*p.  Difference
C   approximations to derivatives and terms defining phi(Y)
C   are returned in D1Y, D2Y, and A.  The Max-norm of Y is
C   returned in YMAX, and DY is the relative change in Y:
C   Max-norm(Y-Y1)/(1+YMAX) = ALPHA*Max-norm(p)/(1+YMAX).
C
C The initial estimate for ALPHA should be small.  Using
C   the step-size computed at the previously iteration can
C   result in an incorrect step-size caused by finding a
C   local minimum which is not a global minimum.
C
      DO 9 I = 1,M
        WK(I,1) = Y(I)
    9   CONTINUE
      IF (OPT) THEN
        ALPHA = SSIZ1
      ELSE
C
C   Use small constant step-size.
C
        ALPHA = SSIZ0
      ENDIF
      CALL LNSRCH (M,H,WK,PHI0,P,FMTOL,OPT, ALPHA,Y,D1Y,
     .             D2Y,A, PHIY,NPHVAL,YMAX,DY,IER)
      IF (IER .NE. 0) GO TO 25
      NPHSM = NPHSM + NPHVAL
C
C Adjust DY to the relative change in Y squared.
C
      DY = DY**2
C
C Compute the max-norm of the Sobolev gradient relative to
C   the solution y.
C
      DG = 0.D0
      DO 10 I = 1,M
        DG = MAX(DG,ABS(G(I)))
   10   CONTINUE
      DG = (DG/(1.D0+YMAX))**3
C
C Save a copy of the Sobolev gradient G in G0 for computing
C   the search direction in the next conjugate gradient
C   step.
C
      IF (NCG .LT. MAXCG) THEN
        DO 11 I = 1,M
          G0(I) = G(I)
   11     CONTINUE
      ENDIF
C
C Compute the relative change DPHI in phi and update PHI0.
C   Terminate on an increase in phi.
C
      DPHI = (PHI0-PHIY)/(1.D0+PHI0)
      PHI0 = PHIY
      IF (DPHI .LT. 0.D0) NGMAX = NG
C
C Print statistics.
C
      IF (PRNT) WRITE (LPRT,150) NPHVAL, ALPHA, PHIY, DPHI,
     .                           DY, DG
  150 FORMAT (/5X,'Line search:     No. phi evals NPHV = ',
     .            I3/
     .         5X,'                Optimal step-size S = ',
     .            1P,D9.3/
     .         5X,'                  Functional phi(Y) = ',
     .            D21.15/
     .         5X,'     Relative change in phi(Y) DPHI = ',
     .            D9.2/
     .         5X,'  Squared relative change in Y:  DY = ',
     .            D9.2/
     .         5X,'    Cubed S-gradient rel. to Y:  DG = ',
     .            D9.2/
     .         1X,'______________________________________',
     .            '___________________________')
C
C Test for termination.
C
      IF ((DPHI .GE. PHTOL  .OR.  DY .GE. PHTOL  .OR.
     .     DG .GE. PHTOL)  .AND.  NG .LT. NGMAX) GO TO 7
      IF (.NOT. PRNT) THEN
C
C Print parameter values.
C
        IF (NCG .EQ. 0) THEN
          WRITE (LPRT,130) NG
        ELSE
          WRITE (LPRT,135) NG
        ENDIF
        WRITE (LPRT,140) PHI0, GNRM1, GNRMM
        WRITE (LPRT,150) NPHVAL, ALPHA, PHIY, DPHI, DY, DG
      ENDIF
C
C Compute the curve length in DY.
C
      C2 = H*H
      DY = 0.D0
      DO 12 I = 2,M
        DY = DY + SQRT(C2 + (Y(I)-Y(I-1))**2)
   12   CONTINUE
C
C Print iteration counts and statistics.
C
      WRITE (LPRT,160) NG
  160 FORMAT (//5X,'Number of descent ',
     .        'iterations: ',I6)
      WRITE (LPRT,170) NPHSM
  170 FORMAT (/5X,'Total number of phi evaluations ',
     .        'in the line searches: ',I7)
      WRITE (LPRT,180) DY
  180 FORMAT (/5X,'Total polygonal curve length: ',D9.3)
C
C Write the solution to disk.
C
C     OPEN (LOUT,FILE='dnsplin1.out')
C     CLOSE (LOUT,STATUS='DELETE')
C     OPEN (LOUT,FILE='dnsplin1.out')
      CALL WRITF (LOUT,M,Y,D1Y,D2Y)
C
C Plot the solution.
C
      OPEN (LPLT,FILE='dnsplin1.ps')
      CLOSE (LPLT,STATUS='DELETE')
      OPEN (LPLT,FILE='dnsplin1.ps')
      CALL PLTCRV (LPLT,N,T(1),T(N),T,D,M,Y,PLT0,Y0,
     .             PLTSIZ, IER)
      CLOSE (LPLT)
      IF (IER .NE. 0) GO TO 26
C
C Prompt for optional plot of curvature.
C
C  13 WRITE (*,220)
C 220 FORMAT (//5X,'Create file dnsplin1.ps2 containing a ',
C    .        'curvature plot? (y/n)'/)
C     READ (*,215,ERR=13) CHR
C     IF (CHR .NE. 'Y'  .AND.  CHR .NE. 'y'  .AND.
C    .    CHR .NE. 'N'  .AND.  CHR .NE. 'n') GO TO 13
C     IF (CHR .EQ. 'Y'  .OR.  CHR .EQ. 'y') THEN
        OPEN (LPLT,FILE='dnsplin1.ps2')
        CLOSE (LPLT,STATUS='DELETE')
        OPEN (LPLT,FILE='dnsplin1.ps2')
        Y(1) = 0.D0
        DO 14 I = 2,M-1
          Y(I) = D2Y(I)/(1.D0+D1Y(I)**2)**1.5D0
   14     CONTINUE
        Y(M) = 0.D0
        DO 15 J = 1,N
          I = IND(J)
          D(J) = Y(I)
   15     CONTINUE
        CALL PLTCRV (LPLT,N,T(1),T(N),T,D,M,Y,.FALSE.,Y0,
     .               PLTSIZ, IER)
        CLOSE (LPLT)
        IF (IER .NE. 0) GO TO 26
C     ENDIF
      STOP
C
C Error encountered in Subroutine READF.
C
   21 WRITE (LPRT,421) IER
  421 FORMAT (///10X,'*** Error in READF:  IER = ',
     .        I1,' ***')
      STOP
C
C Invalid data:  M < 6 or M > MMAX or M-N < 1.
C
   22 WRITE (LPRT,422) M, MMAX, N
  422 FORMAT (///10X,'*** Error in data:  M = ',I5,
     .        ', MMAX = ',I5,', N = ',I4,' ***')
      STOP
C
C Error encountered in Subroutine CSPLIN.
C
   23 WRITE (LPRT,423) IER
  423 FORMAT (///10X,'*** Error in CSPLIN:  IER = ',
     .        I1,' ***')
      STOP
C
C Error encountered in Subroutine GRADL2.
C
   24 WRITE (LPRT,424) IER
  424 FORMAT (///10X,'*** Error in GRADL2:  IER = ',
     .        I1,' ***')
      STOP
C
C Error encountered in Subroutine LNSRCH.
C
   25 WRITE (LPRT,425) IER
  425 FORMAT (///10X,'*** Error in LNSRCH:  IER = ',
     .        I1,' ***')
      STOP
C
C Error encountered in Subroutine PLTCRV.
C
   26 WRITE (LPRT,426) IER
  426 FORMAT (///10X,'*** Error in PLTCRV:  IER = ',
     .        I1,' ***')
      STOP
      END
