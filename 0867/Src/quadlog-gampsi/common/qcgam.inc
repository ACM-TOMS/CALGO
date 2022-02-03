      INTEGER             nc
      PARAMETER           (nc = 30)
*
      REAL*16             c(nc)
*-----------------------------------------------------------------------
*     Coefficients for the asymptotic series of ln(Gamma(x)), for x
*     in (cutoff, INF).  See Abramowitz and Stegen, NBS Handbook #55,
*     p. 257, 6.1.41 (1964).  These can be generated in Maple V5.1 by
*
*         Digits := 75:
*         with(numtheory):
*         bprint := proc(m)
*                  printf("     X      %+.50e,\n", B(2*m)/(2*m*(2*m-1)))
*              end:
*         for k from 1 to 25 do bprint(k) od:
*
*     and then substituting the Maple exponent letter e with Fortran q.
*
*     NB: ACM Algorithm 715, file gamma.tpl, had a peculiar order of the
*     coefficients.  This version instead defines the cofficients in
*     standard order, and the code below in loop 300 has been changed
*     accordingly.
*
*     Because of the Fortran limit of 19 continuation lines, we must
*     split this initializer into two parts.
*-----------------------------------------------------------------------
      DATA (c(i), i = 1,19) /
     X      +8.33333333333333333333333333333333333333333333333333q-02,
     X      -2.77777777777777777777777777777777777777777777777778q-03,
     X      +7.93650793650793650793650793650793650793650793650794q-04,
     X      -5.95238095238095238095238095238095238095238095238095q-04,
     X      +8.41750841750841750841750841750841750841750841750842q-04,
     X      -1.91752691752691752691752691752691752691752691752692q-03,
     X      +6.41025641025641025641025641025641025641025641025641q-03,
     X      -2.95506535947712418300653594771241830065359477124183q-02,
     X      +1.79644372368830573164938490015889396694350254721772q-01,
     X      -1.39243221690590111642743221690590111642743221690590q+00,
     X      +1.34028640441683919944789510006901311249137336093858q+01,
     X      -1.56848284626002017306365132452088973828104262886872q+02,
     X      +2.19310333333333333333333333333333333333333333333333q+03,
     X      -3.61087712537249893571732652192422307364836100468284q+04,
     X      +6.91472268851313067108395250775673467553334071687798q+05,
     X      -1.52382215394074161922833649588867805186590765338393q+07,
     X      +3.82900751391414141414141414141414141414141414141414q+08,
     X      -1.08822660357843910890151491655251053747294348798108q+10,
     X      +3.47320283765002252252252252252252252252252252252252q+11 /
      DATA (c(i), i = 20,30) /
     X      -1.23696021422692744542517103492713248810809786419543q+13,
     X      +4.88788064793079335075815162518022902108470538905674q+14,
     X      -2.13203339609193738969750589821368385574654533198517q+16,
     X      +1.02177529652570007756528762805358550039401103230890q+18,
     X      -5.35754721733002036108277091919692044848490405436588q+19,
     X      +3.06157826370488341504315105132962275819418676561534q+21,
     X      -1.89999174263992040502937142930694290294734245899618q+23,
     X      +1.27633740338288341492349513776978259765416336088299q+25,
     X      -9.25284717612041630723024234834762277951933124346917q+26,
     X      +7.21882259518561029783605018730163792248984042025969q+28,
     X      -6.04518340599585696774314823875454728606614439596720q+30 /