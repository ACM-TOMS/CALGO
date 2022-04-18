      REAL             FUNCTION SSINPX (X)
C>> 1996-01-29 SCOSPX WVS JPL Add better acknowledgement for origins.
C>> 1994-10-20 SSINPX Krogh  Changes to use M77CON
C>> 1994-04-22 SSINPX CLL Made SP and DP codes similar.
C>> 1993-05-06 SSINPX WVS JPL Convert from NSWC to Math 77
c--S replaces "?": ?SINPX, ?ERM1
C ----------------------------------------------------------------------
c This procedure was originally procedure SIN1 from the Naval Surface
c Warfare Center library.
C ----------------------------------------------------------------------
C
C                        EVALUATION OF SIN(PI*X)
C
C                             --------------
C
C     THE EXPANSION FOR SIN(PI*A) (ABS(A) .LE. 1/4) USING A1,...,A13
C     IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
C     THE EXPANSION FOR COS(PI*A) (ABS(A) .LE. 1/4) USING B1,...,B13
C     IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.
c
c     The polynomials using coefficients SA0, ... SA6 and SB1, ..., SB6
c     give approximations whose largest observed relative error in the
c     relevant intervals is 0.129e-14.
c     We will use this latter approximation when the machine epsilon
c     is larger than 0.2e-14.
C-----------------------------------------------------------------------
      integer N
      real             R1MACH
      REAL             A, BIG, CUTOFF, EPS, PI, T, W, X
      REAL             A1, A2, A3, A4, A5, A6, A7, A8, A9, A10,
     *                 A11, A12, A13
      REAL             B1, B2, B3, B4, B5, B6, B7, B8, B9, B10,
     *                 B11, B12, B13
      REAL             SA0, SA1, SA2, SA3, SA4, SA5, SA6
      REAL                  SB1, SB2, SB3, SB4, SB5, SB6
      SAVE BIG, EPS
C------------------------
      PARAMETER ( PI = 3.141592653589793238462643383279502884197E+00 )
      parameter ( CUTOFF = 0.2e-14 )
C------------------------
      data EPS / -1.0e0 /
      DATA SA0 /.314159265358979E+01/, SA1 /-.516771278004995E+01/,
     *     SA2 /.255016403987327E+01/, SA3 /-.599264528932149E+00/,
     *     SA4 /.821458689493251E-01/, SA5 /-.737001831310553E-02/,
     *     SA6 /.461514425296398E-03/
      DATA SB1 /-.493480220054460E+01/, SB2 /.405871212639605E+01/,
     *     SB3 /-.133526276691575E+01/, SB4 /.235330543508553E+00/,
     *     SB5 /-.258048861575714E-01/, SB6 /.190653140279462E-02/
      DATA A1  /-.1028083791780141522795259479153765743002E+00/,
     *     A2  / .3170868848763100170457042079710451905600E-02/,
     *     A3  /-.4657026956105571623449026167864697920000E-04/,
     *     A4  / .3989844942879455643410226655783424000000E-06/,
     *     A5  /-.2237397227721999776371894030796800000000E-08/,
     *     A6  / .8847045483056962709715066675200000000000E-11/,
     *     A7  /-.2598715447506450292885585920000000000000E-13/,
     *     A8  / .5893449774331011070033920000000000000000E-16/,
     *     A9  /-.1062975472045522550784000000000000000000E-18/,
     *     A10 / .1561182648301780992000000000000000000000E-21/,
     *     A11 /-.1903193516670976000000000000000000000000E-24/,
     *     A12 / .1956617650176000000000000000000000000000E-27/,
     *     A13 /-.1711276032000000000000000000000000000000E-30/
C------------------------
      DATA B1  /-.3084251375340424568385778437461297229882E+00/,
     *     B2  / .1585434424381550085228521039855226435920E-01/,
     *     B3  /-.3259918869273900136414318317506279360000E-03/,
     *     B4  / .3590860448591510079069203991239232000000E-05/,
     *     B5  /-.2461136950494199754009084061808640000000E-07/,
     *     B6  / .1150115912797405152263195572224000000000E-09/,
     *     B7  /-.3898073171259675439899172864000000000000E-12/,
     *     B8  / .1001886461636271969091584000000000000000E-14/,
     *     B9  /-.2019653396886572027084800000000000000000E-17/,
     *     B10 / .3278483561466560512000000000000000000000E-20/,
     *     B11 /-.4377345082051788800000000000000000000000E-23/,
     *     B12 / .4891532381388800000000000000000000000000E-26/,
     *     B13 /-.4617089843200000000000000000000000000000E-29/
C------------------------
      IF (eps .LT. 0.0e0) then
         eps = R1MACH(3)
         BIG = 1.0e0/eps
      endif
C------------------------
      A = ABS(X)
      IF (A .GE. BIG) THEN
         CALL SERM1 ('SSINPX',1,2,
     1   'No precision because ABS(X) is too large','X',X,'.')
         SSINPX = 0.E0
         RETURN
      END IF
      N = A
      T = N
      A = A - T
      IF (A .GT. 0.75E0) GO TO 20
      IF (A .LT. 0.25E0) GO TO 21
C
C                    0.25 .LE. A .LE. 0.75
C
      A = 0.25E0 + (0.25E0 - A)
      if(eps .lt. cutoff) then
         T = 16.E0*A*A
         SSINPX = (((((((((((((B13*T + B12)*T + B11)*T + B10)*T + B9)*T+
     *                B8)*T + B7)*T + B6)*T + B5)*T + B4)*T + B3)*T +
     *                B2)*T + B1)*T + 0.5E0) + 0.5E0
      else
         T = A*A
         SSINPX = ((((((SB6*T + SB5)*T + SB4)*T + SB3)*T + SB2)*T
     *                  + SB1)*T + 0.5E0) + 0.5E0
      endif
      GO TO 30
C
C                 A .LT. 0.25  OR  A .GT. 0.75
C
   20 A = 0.25E0 + (0.75E0 - A)
   21 continue
      if(eps .lt. cutoff) then
         T = 16.E0*A*A
         W = (((((((((((((A13*T + A12)*T + A11)*T + A10)*T + A9)*T +
     *            A8)*T + A7)*T + A6)*T + A5)*T + A4)*T + A3)*T +
     *            A2)*T + A1)*T + 0.5E0) + 0.5E0
         SSINPX = PI*A*W
      else
         T = A*A
         SSINPX = ((((((SA6*T + SA5)*T + SA4)*T + SA3)*T + SA2)*T
     *                  + SA1)*T + SA0)*A
      endif
C
C                        TERMINATION
C
   30 IF (X .LT. 0.0E0) SSINPX = - SSINPX
      IF (MOD(N,2) .NE. 0) SSINPX = - SSINPX
      RETURN
      END
