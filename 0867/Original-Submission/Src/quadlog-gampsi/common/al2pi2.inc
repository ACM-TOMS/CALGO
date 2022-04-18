*
*     In ACM Algorithm 715, file gamma.tpl, the variable sqrtpi is
*     misnamed: its value is ln(2*pi)/2, not sqrt(pi), so in this
*     program, we change its name to ln2pi2:
*
*     Maple V6:
*         Digits := 75;
*         printf("%.55f\n", evalf(log(2*Pi)/2)):
*             .9189385332046727417803297364056176398613974736377834128
*
      REAL                ln2pi2
      PARAMETER           (ln2pi2 =
     X    0.9189385332046727417803297364056176398613974736377834128e+00)
