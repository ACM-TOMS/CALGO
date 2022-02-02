c     program DRSFRENL
c>> 1996-01-29 DRSFRENL WV Snyder Corrected formats
c>> 1994-10-19 DRSFRENL Krogh  Changes to use M77CON
c>> 1993-02-25 DRSFRENL CLL. Minor edits.  Deleted Format statements.
c>> 1992-03-18 DRSFRENL WV Snyder Corrected failure to call SFRENS
c>> 1992-03-18 DRSFRENL CLL Minor edits.
c>> 1992-02-03 DRSFRENL WV Snyder JPL Original code.
c
c     Demonstration driver for Fresnel Integrals functions.
c
c     ------------------------------------------------------------------
c--S replaces "?": DR?FRENL, ?FRENC, ?FRENS, ?FRENG, ?FRENF
c     ------------------------------------------------------------------
      real             X, YC, YF, YG, YS
      real             SFRENC, SFRENF, SFRENG, SFRENS
      external SFRENC, SFRENF, SFRENG, SFRENS
      integer I
c
      print '(11x,a,12x,a,11x,a,11x,a,11x,a)',
     *      'X', 'C(x)', 'S(x)', 'g(x)', 'f(x)'
      do 30 I = -12, 12
         X = 0.5 * I
         YC = SFRENC(X)
         YS = SFRENS(X)
         YG = SFRENG(X)
         YF = SFRENF(X)
         print '(1p5e15.07)', X, YC, YS, YG, YF
30    continue
      stop
      end
