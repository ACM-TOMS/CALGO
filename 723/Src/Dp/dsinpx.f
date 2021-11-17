      DOUBLE PRECISION FUNCTION DSINPX (X)
C>> 1996-01-29 SCOSPX WVS JPL Add better acknowledgement for origins.
C>> 1994-10-20 DSINPX Krogh  Changes to use M77CON
C>> 1994-04-22 DSINPX CLL Made SP and DP codes similar.
C>> 1993-05-06 DSINPX WVS JPL Convert from NSWC to Math 77
c--D replaces "?": ?SINPX, ?ERM1
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
c     relevant intervals is 0.129d-14.
c     We will use this latter approximation when the machine epsilon
c     is larger than 0.2d-14.
C-----------------------------------------------------------------------
      integer N
      double precision D1MACH
      DOUBLE PRECISION A, BIG, CUTOFF, EPS, PI, T, W, X
      DOUBLE PRECISION A1, A2, A3, A4, A5, A6, A7, A8, A9, A10,
     *                 A11, A12, A13
      DOUBLE PRECISION B1, B2, B3, B4, B5, B6, B7, B8, B9, B10,
     *                 B11, B12, B13
      DOUBLE PRECISION SA0, SA1, SA2, SA3, SA4, SA5, SA6
      DOUBLE PRECISION      SB1, SB2, SB3, SB4, SB5, SB6
      SAVE BIG, EPS
C------------------------
      PARAMETER ( PI = 3.141592653589793238462643383279502884197D+00 )
      parameter ( CUTOFF = 0.2d-14 )
C------------------------
      data EPS / -1.0d0 /
      DATA SA0 /.314159265358979D+01/, SA1 /-.516771278004995D+01/,
     *     SA2 /.255016403987327D+01/, SA3 /-.599264528932149D+00/,
     *     SA4 /.821458689493251D-01/, SA5 /-.737001831310553D-02/,
     *     SA6 /.461514425296398D-03/
      DATA SB1 /-.493480220054460D+01/, SB2 /.405871212639605D+01/,
     *     SB3 /-.133526276691575D+01/, SB4 /.235330543508553D+00/,
     *     SB5 /-.258048861575714D-01/, SB6 /.190653140279462D-02/
      DATA A1  /-.1028083791780141522795259479153765743002D+00/,
     *     A2  / .3170868848763100170457042079710451905600D-02/,
     *     A3  /-.4657026956105571623449026167864697920000D-04/,
     *     A4  / .3989844942879455643410226655783424000000D-06/,
     *     A5  /-.2237397227721999776371894030796800000000D-08/,
     *     A6  / .8847045483056962709715066675200000000000D-11/,
     *     A7  /-.2598715447506450292885585920000000000000D-13/,
     *     A8  / .5893449774331011070033920000000000000000D-16/,
     *     A9  /-.1062975472045522550784000000000000000000D-18/,
     *     A10 / .1561182648301780992000000000000000000000D-21/,
     *     A11 /-.1903193516670976000000000000000000000000D-24/,
     *     A12 / .1956617650176000000000000000000000000000D-27/,
     *     A13 /-.1711276032000000000000000000000000000000D-30/
C------------------------
      DATA B1  /-.3084251375340424568385778437461297229882D+00/,
     *     B2  / .1585434424381550085228521039855226435920D-01/,
     *     B3  /-.3259918869273900136414318317506279360000D-03/,
     *     B4  / .3590860448591510079069203991239232000000D-05/,
     *     B5  /-.2461136950494199754009084061808640000000D-07/,
     *     B6  / .1150115912797405152263195572224000000000D-09/,
     *     B7  /-.3898073171259675439899172864000000000000D-12/,
     *     B8  / .1001886461636271969091584000000000000000D-14/,
     *     B9  /-.2019653396886572027084800000000000000000D-17/,
     *     B10 / .3278483561466560512000000000000000000000D-20/,
     *     B11 /-.4377345082051788800000000000000000000000D-23/,
     *     B12 / .4891532381388800000000000000000000000000D-26/,
     *     B13 /-.4617089843200000000000000000000000000000D-29/
C------------------------
      IF (eps .LT. 0.0d0) then
         eps = D1MACH(3)
         BIG = 1.0d0/eps
      endif
C------------------------
      A = ABS(X)
      IF (A .GE. BIG) THEN
         CALL DERM1 ('DSINPX',1,2,
     1   'No precision because ABS(X) is too large','X',X,'.')
         DSINPX = 0.D0
         RETURN
      END IF
      N = A
      T = N
      A = A - T
      IF (A .GT. 0.75D0) GO TO 20
      IF (A .LT. 0.25D0) GO TO 21
C
C                    0.25 .LE. A .LE. 0.75
C
      A = 0.25D0 + (0.25D0 - A)
      if(eps .lt. cutoff) then
         T = 16.D0*A*A
         DSINPX = (((((((((((((B13*T + B12)*T + B11)*T + B10)*T + B9)*T+
     *                B8)*T + B7)*T + B6)*T + B5)*T + B4)*T + B3)*T +
     *                B2)*T + B1)*T + 0.5D0) + 0.5D0
      else
         T = A*A
         DSINPX = ((((((SB6*T + SB5)*T + SB4)*T + SB3)*T + SB2)*T
     *                  + SB1)*T + 0.5D0) + 0.5D0
      endif
      GO TO 30
C
C                 A .LT. 0.25  OR  A .GT. 0.75
C
   20 A = 0.25D0 + (0.75D0 - A)
   21 continue
      if(eps .lt. cutoff) then
         T = 16.D0*A*A
         W = (((((((((((((A13*T + A12)*T + A11)*T + A10)*T + A9)*T +
     *            A8)*T + A7)*T + A6)*T + A5)*T + A4)*T + A3)*T +
     *            A2)*T + A1)*T + 0.5D0) + 0.5D0
         DSINPX = PI*A*W
      else
         T = A*A
         DSINPX = ((((((SA6*T + SA5)*T + SA4)*T + SA3)*T + SA2)*T
     *                  + SA1)*T + SA0)*A
      endif
C
C                        TERMINATION
C
   30 IF (X .LT. 0.0D0) DSINPX = - DSINPX
      IF (MOD(N,2) .NE. 0) DSINPX = - DSINPX
      RETURN
      END
