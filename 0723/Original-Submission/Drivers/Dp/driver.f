c     program DRDFRENL
c>> 1996-01-29 DRDFRENL WV Snyder Corrected formats
c>> 1994-10-19 DRDFRENL Krogh  Changes to use M77CON
c>> 1993-02-25 DRDFRENL CLL. Minor edits.  Deleted Format statements.
c>> 1992-03-18 DRDFRENL WV Snyder Corrected failure to call DFRENS
c>> 1992-03-18 DRDFRENL CLL Minor edits.
c>> 1992-02-03 DRDFRENL WV Snyder JPL Original code.
c
c     Demonstration driver for Fresnel Integrals functions.
c
c     ------------------------------------------------------------------
c--D replaces "?": DR?FRENL, ?FRENC, ?FRENS, ?FRENG, ?FRENF
c     ------------------------------------------------------------------
      double precision X, YC, YF, YG, YS
      double precision DFRENC, DFRENF, DFRENG, DFRENS
      external DFRENC, DFRENF, DFRENG, DFRENS
      integer I
c
      print '(11x,a,12x,a,11x,a,11x,a,11x,a)',
     *      'X', 'C(x)', 'S(x)', 'g(x)', 'f(x)'
      do 30 I = -12, 12
         X = 0.5 * I
         YC = DFRENC(X)
         YS = DFRENS(X)
         YG = DFRENG(X)
         YF = DFRENF(X)
         print '(1p5e15.07)', X, YC, YS, YG, YF
30    continue
      stop
      end
