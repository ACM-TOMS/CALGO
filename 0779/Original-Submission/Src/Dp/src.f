      DOUBLE PRECISION FUNCTION FDM0P5(XVALUE)
C
C   DESCRIPTION:
C
C      This function computes the Fermi-Dirac function of
C      order -1/2, defined as
C
C                     Int{0 to inf} t**(-1/2) / (1+exp(t-x)) dt
C         FDM0P5(x) = -----------------------------------------
C                                 Gamma(1/2)
C
C      The function uses Chebyshev expansions which are given to
C      16 decimal places for x <= 2, but only 10 decimal places 
C      for x > 2.
C
C
C   ERROR RETURNS:
C    
C      None.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERMS1 - INTEGER - The number of terms used from the array
C                          ARRFD1. The recommended value is such that
C                               ABS(ARRFD1(NTERMS1)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 14.
C
C      NTERMS2 - INTEGER - The number of terms used from the array
C                          ARRFD2. The recommended value is such that
C                               ABS(ARRFD2(NTERMS2)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 23.
C
C      NTERMS3 - INTEGER - The number of terms used from the array
C                          ARRFD3. The recommended value is such that
C                               ABS(ARRFD3(NTERMS3)) < EPS/10
C                          subject to 1 <= NTERMS3 <= 28.
C
C      XMIN1 - REAL - The value of x below which
C                         FDM0P5(x) = exp(x)
C                     to machine precision. The recommended value
C                     is    LN ( SQRT(2) * EPSNEG )
C
C      XMIN2 - REAL - The value of x below which
C                         FDM0P5(x) = 0.0 
C                     to machine precision. The recommended value
C                     is    LN ( XMIN )
C
C      XHIGH - REAL - The value of x above which
C                         FDM0P5(x) = 2 sqrt (x/pi) 
C                     to machine precision. The recommended value
C                     is    1 / sqrt( 2 * EPSNEG )
C
C      For values of EPS, EPSNEG, and XMIN the user should refer to the
C      paper by Cody in ACM. Trans. Math. Soft. Vol. 14 (1988) p303-311.
C   
C      This code is provided with single and double precision values
C      of the machine-dependent parameters, suitable for machines
C      which satisfy the IEEE floating-point standard.
C
C
C   AUTHOR:
C          DR. ALLAN MACLEOD,
C          DEPT. OF MATHEMATICS AND STATISTICS,
C          UNIVERSITY OF PAISLEY,
C          HIGH ST.,
C          PAISLEY,
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl-ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:
C                 20 NOVEMBER, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION
     1     ARRFD1(0:14),ARRFD2(0:23),ARRFD3(0:58),
     2     CHEVAL,CHV,EXPX,FIFTY,FORTY2,
     3     GAM1P5,ONE,T,THREE,TWO,TWOE,
     4     X,XHIGH,XMIN1,XMIN2,XSQ,XVALUE,ZERO
      DATA ARRFD1/1.7863 5963 8510 2264  D   0,
     1           -0.9993 7200 7632 333   D  -1,
     2            0.6414 4652 2160 54    D  -2,
     3           -0.4356 4153 7134 5     D  -3,
     4            0.3052 1670 0310       D  -4,
     5           -0.2181 0648 110        D  -5,
     6            0.1580 0507 81         D  -6,
     7           -0.1156 2057 0          D  -7,
     8            0.8525 860             D  -9,
     9           -0.6325 29              D -10,
     X            0.4715 9               D -11,
     1           -0.3530                 D -12,
     2            0.265                  D -13,
     3           -0.20                   D -14,
     4            0.2                    D -15/
      DATA ARRFD2( 0)/ 1.6877 1115 2605 2352  D   0/
      DATA ARRFD2( 1)/ 0.5978 3602 2633 6983  D   0/
      DATA ARRFD2( 2)/ 0.3572 2600 4541 669   D  -1/
      DATA ARRFD2( 3)/-0.1321 4478 6506 426   D  -1/
      DATA ARRFD2( 4)/-0.4040 1342 0744 7     D  -3/
      DATA ARRFD2( 5)/ 0.5330 0118 4688 7     D  -3/
      DATA ARRFD2( 6)/-0.1489 2350 4863       D  -4/
      DATA ARRFD2( 7)/-0.2188 6382 2916       D  -4/
      DATA ARRFD2( 8)/ 0.1965 2084 277        D  -5/
      DATA ARRFD2( 9)/ 0.8565 8304 66         D  -6/
      DATA ARRFD2(10)/-0.1407 7231 33         D  -6/
      DATA ARRFD2(11)/-0.3051 7580 3          D  -7/
      DATA ARRFD2(12)/ 0.8352 4532            D  -8/
      DATA ARRFD2(13)/ 0.9025 750             D  -9/
      DATA ARRFD2(14)/-0.4455 471             D  -9/
      DATA ARRFD2(15)/-0.1483 42              D -10/
      DATA ARRFD2(16)/ 0.2192 66              D -10/
      DATA ARRFD2(17)/-0.6579                 D -12/
      DATA ARRFD2(18)/-0.1000 9               D -11/
      DATA ARRFD2(19)/ 0.936                  D -13/
      DATA ARRFD2(20)/ 0.420                  D -13/
      DATA ARRFD2(21)/-0.71                   D -14/
      DATA ARRFD2(22)/-0.16                   D -14/
      DATA ARRFD2(23)/ 0.4                    D -15/
      DATA ARRFD3(0)/  0.8707 1950 2959 0563  D    0/
      DATA ARRFD3(1)/  0.5983 3110 2317 33    D   -2/
      DATA ARRFD3(2)/ -0.4326 7047 0895 746   D   -1/
      DATA ARRFD3(3)/ -0.3930 8368 1608 590   D   -1/
      DATA ARRFD3(4)/ -0.1914 8268 8045 932   D   -1/
      DATA ARRFD3(5)/ -0.6558 2880 9801 58    D   -2/
      DATA ARRFD3(6)/ -0.2227 6691 5163 12    D   -2/
      DATA ARRFD3(7)/ -0.8466 7869 3617 8     D   -3/
      DATA ARRFD3(8)/ -0.2807 4594 8921 9     D   -3/
      DATA ARRFD3(9)/ -0.9555 7502 4348       D   -4/
      DATA ARRFD3(10)/-0.3623 6766 2803       D   -4/
      DATA ARRFD3(11)/-0.1091 5846 8869       D   -4/
      DATA ARRFD3(12)/-0.3935 6701 000        D   -5/
      DATA ARRFD3(13)/-0.1310 8192 725        D   -5/
      DATA ARRFD3(14)/-0.2468 8163 88         D   -6/
      DATA ARRFD3(15)/-0.1048 3803 11         D   -6/
      DATA ARRFD3(16)/ 0.2361 8148 7          D   -7/
      DATA ARRFD3(17)/ 0.2271 4535 9          D   -7/
      DATA ARRFD3(18)/ 0.1457 7517 4          D   -7/
      DATA ARRFD3(19)/ 0.1539 2676 7          D   -7/
      DATA ARRFD3(20)/ 0.5692 4772            D   -8/
      DATA ARRFD3(21)/ 0.5062 3068            D   -8/
      DATA ARRFD3(22)/ 0.2342 6075            D   -8/
      DATA ARRFD3(23)/ 0.1265 2275            D   -8/
      DATA ARRFD3(24)/ 0.8927 773             D   -9/
      DATA ARRFD3(25)/ 0.2994 501             D   -9/
      DATA ARRFD3(26)/ 0.2822 785             D   -9/
      DATA ARRFD3(27)/ 0.9106 85              D  -10/
      DATA ARRFD3(28)/ 0.6962 85              D  -10/
      DATA ARRFD3(29)/ 0.3662 25              D  -10/
      DATA ARRFD3(30)/ 0.1243 51              D  -10/
      DATA ARRFD3(31)/ 0.1450 19              D  -10/
      DATA ARRFD3(32)/ 0.1664 5               D  -11/
      DATA ARRFD3(33)/ 0.4585 6               D  -11/
      DATA ARRFD3(34)/ 0.6092                 D  -12/
      DATA ARRFD3(35)/ 0.9331                 D  -12/
      DATA ARRFD3(36)/ 0.5238                 D  -12/
      DATA ARRFD3(37)/-0.56                   D  -14/
      DATA ARRFD3(38)/ 0.3170                 D  -12/
      DATA ARRFD3(39)/-0.926                  D  -13/
      DATA ARRFD3(40)/ 0.1265                 D  -12/
      DATA ARRFD3(41)/-0.327                  D  -13/
      DATA ARRFD3(42)/ 0.276                  D  -13/
      DATA ARRFD3(43)/ 0.33                   D  -14/
      DATA ARRFD3(44)/-0.42                   D  -14/
      DATA ARRFD3(45)/ 0.101                  D  -13/
      DATA ARRFD3(46)/-0.73                   D  -14/
      DATA ARRFD3(47)/ 0.64                   D  -14/
      DATA ARRFD3(48)/-0.37                   D  -14/
      DATA ARRFD3(49)/ 0.23                   D  -14/
      DATA ARRFD3(50)/-0.9                    D  -15/
      DATA ARRFD3(51)/ 0.2                    D  -15/
      DATA ARRFD3(52)/ 0.2                    D  -15/
      DATA ARRFD3(53)/-0.3                    D  -15/
      DATA ARRFD3(54)/ 0.4                    D  -15/
      DATA ARRFD3(55)/-0.3                    D  -15/
      DATA ARRFD3(56)/ 0.2                    D  -15/
      DATA ARRFD3(57)/-0.1                    D  -15/
      DATA ARRFD3(58)/ 0.1                    D  -15/
      DATA ZERO,ONE,TWO/ 0.0 D 0 , 1.0 D 0 , 2.0 D 0/
      DATA THREE,FORTY2,FIFTY/ 3.0 D 0 , 42.0 D 0 , 50.0 D 0/
      DATA GAM1P5/0.8862 2692 5452 7580 D 0/
      DATA TWOE/5.4365 6365 6918 0905 D 0/
C
C   Machine-dependent constants
C
      DATA NTERM1,NTERM2,NTERM3/14,23,58/
      DATA XMIN1,XMIN2,XHIGH/-36.39023D0,-708.39641D0,67108864.0D0/
C
C   Start calculation
C
      X=XVALUE
C
C   Code for x < -1
C
      IF ( X .LT. -ONE ) THEN
         IF ( X .GT. XMIN1 ) THEN
            EXPX = EXP(X)
            T = TWOE * EXPX - ONE
            FDM0P5 = EXPX * CHEVAL ( NTERM1 , ARRFD1 , T )
         ELSE
            IF ( X .LT. XMIN2 ) THEN
               FDM0P5 = ZERO
            ELSE
               FDM0P5 = EXP(X)
            ENDIF
         ENDIF
      ELSE
C
C   Code for -1 <= x <= 2
C
         IF ( X .LE. TWO ) THEN
            T = ( TWO * X - ONE ) / THREE
            FDM0P5 = CHEVAL ( NTERM2 , ARRFD2 , T )
         ELSE
C
C   Code for x > 2
C
            FDM0P5 = SQRT(X) / GAM1P5
            IF ( X .LE. XHIGH ) THEN 
               XSQ = X * X
               T = ( FIFTY - XSQ ) / ( FORTY2 + XSQ )
               CHV = CHEVAL ( NTERM3 , ARRFD3 , T )
               FDM0P5 = FDM0P5 * ( ONE - CHV / XSQ )
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END



      DOUBLE PRECISION FUNCTION FDP0P5(XVALUE)
C
C   DESCRIPTION:
C
C      This function computes the Fermi-Dirac function of
C      order 1/2, defined as
C
C                     Int{0 to inf} t**(1/2) / (1+exp(t-x)) dt
C         FDP0P5(x) = -----------------------------------------
C                                 Gamma(3/2)
C
C      The function uses Chebyshev expansions which are given to
C      16 decimal places for x <= 2, but only 10 decimal places 
C      for x > 2.
C
C
C   ERROR RETURNS:
C    
C      If XVALUE too large and positive, the function value
C      will overflow. An error message is printed and the function
C      returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERMS1 - INTEGER - The number of terms used from the array
C                          ARRFD1. The recommended value is such that
C                               ABS(ARRFD1(NTERMS1)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 13.
C
C      NTERMS2 - INTEGER - The number of terms used from the array
C                          ARRFD2. The recommended value is such that
C                               ABS(ARRFD2(NTERMS2)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 23.
C
C      NTERMS3 - INTEGER - The number of terms used from the array
C                          ARRFD3. The recommended value is such that
C                               ABS(ARRFD3(NTERMS3)) < EPS/10
C                          subject to 1 <= NTERMS3 <= 32.
C
C      XMIN1 - REAL - The value of x below which
C                         FDP0P5(x) = exp(x)
C                     to machine precision. The recommended value
C                     is   1.5*LN(2) + LN(EPSNEG)
C
C      XMIN2 - REAL - The value of x below which
C                         FDP0P5(x) = 0.0 
C                     to machine precision. The recommended value
C                     is    LN ( XMIN )
C
C      XHIGH1 - REAL - The value of x above which
C                         FDP0P5(x) = x**(3/2)/GAMMA(5/2)
C                     to machine precision. The recommended value
C                     is   pi / SQRT(8*EPS)
C
C      XHIGH2 - REAL - The value of x above which FDP0P5 would 
C                      overflow. The reommended value is
C                              (1.329*XMAX)**(2/3)
C
C      For values of EPS, EPSNEG, and XMIN the user should refer to the
C      paper by Cody in ACM. Trans. Math. Soft. Vol. 14 (1988) p303-311.
C   
C      This code is provided with single and double precision values
C      of the machine-dependent parameters, suitable for machines
C      which satisfy the IEEE floating-point standard.
C
C
C   AUTHOR:
C          DR. ALLAN MACLEOD,
C          DEPT. OF MATHEMATICS AND STATISTICS,
C          UNIVERSITY OF PAISLEY,
C          HIGH ST.,
C          PAISLEY,
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl-ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:
C                 20 NOVEMBER, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION
     1     ARRFD1(0:13),ARRFD2(0:23),ARRFD3(0:53),
     2     CHEVAL,CHV,EXPX,FIFTY,FORTY2,
     3     GAM2P5,ONE,T,THREE,TWO,TWOE,X,XHIGH1,
     4     XHIGH2,XMIN1,XMIN2,XSQ,XVALUE,ZERO
      DATA ARRFD1/1.8862 9683 9273 4597  D   0,
     1           -0.5435 8081 7644 053   D  -1,
     2            0.2364 4975 4397 20    D  -2,
     3           -0.1216 9293 6588 0     D  -3,
     4            0.6869 5130 622        D  -5,
     5           -0.4112 0761 72         D  -6,
     6            0.2563 5162 8          D  -7,
     7           -0.1646 5008            D  -8,
     8            0.1081 948             D  -9,
     9           -0.7239 2               D -11,
     X            0.4915                 D -12,
     1           -0.338                  D -13,
     2            0.23                   D -14,
     3           -0.2                    D -15/
      DATA ARRFD2( 0)/ 2.6982 4927 8817 0612  D   0/
      DATA ARRFD2( 1)/ 1.2389 9141 4113 3012  D   0/
      DATA ARRFD2( 2)/ 0.2291 4393 7981 6278  D   0/
      DATA ARRFD2( 3)/ 0.9031 6534 6872 79    D  -2/
      DATA ARRFD2( 4)/-0.2577 6524 6912 46    D  -2/
      DATA ARRFD2( 5)/-0.5836 8160 5388       D  -4/
      DATA ARRFD2( 6)/ 0.6936 0945 8725       D  -4/
      DATA ARRFD2( 7)/-0.1806 1670 265        D  -5/
      DATA ARRFD2( 8)/-0.2132 1530 005        D  -5/
      DATA ARRFD2( 9)/ 0.1754 9839 51         D  -6/
      DATA ARRFD2(10)/ 0.6653 2547 0          D  -7/
      DATA ARRFD2(11)/-0.1016 7597 7          D  -7/
      DATA ARRFD2(12)/-0.1963 7597            D  -8/
      DATA ARRFD2(13)/ 0.5075 769             D  -9/
      DATA ARRFD2(14)/ 0.4914 69              D -10/
      DATA ARRFD2(15)/-0.2337 37              D -10/
      DATA ARRFD2(16)/-0.6645                 D -12/
      DATA ARRFD2(17)/ 0.1011 5               D -11/
      DATA ARRFD2(18)/-0.313                  D -13/
      DATA ARRFD2(19)/-0.412                  D -13/
      DATA ARRFD2(20)/ 0.38                   D -14/
      DATA ARRFD2(21)/ 0.16                   D -14/
      DATA ARRFD2(22)/-0.3                    D -15/
      DATA ARRFD2(23)/-0.1                    D -15/
      DATA ARRFD3(0)/  2.5484 3841 9800 9122  D    0/
      DATA ARRFD3(1)/  0.5104 3940 8960 652   D   -1/
      DATA ARRFD3(2)/  0.7749 3527 6282 94    D   -2/
      DATA ARRFD3(3)/ -0.7504 1656 5849 53    D   -2/
      DATA ARRFD3(4)/ -0.7754 0826 3202 96    D   -2/
      DATA ARRFD3(5)/ -0.4581 0844 5399 77    D   -2/
      DATA ARRFD3(6)/ -0.2343 1641 5873 63    D   -2/
      DATA ARRFD3(7)/ -0.1178 8049 5135 91    D   -2/
      DATA ARRFD3(8)/ -0.5802 7393 5970 2     D   -3/
      DATA ARRFD3(9)/ -0.2825 3507 0053 7     D   -3/
      DATA ARRFD3(10)/-0.1388 1366 5179 9     D   -3/
      DATA ARRFD3(11)/-0.6806 9508 4875       D   -4/
      DATA ARRFD3(12)/-0.3353 5635 0608       D   -4/
      DATA ARRFD3(13)/-0.1665 3301 8734       D   -4/
      DATA ARRFD3(14)/-0.8271 4908 266        D   -5/
      DATA ARRFD3(15)/-0.4142 5714 409        D   -5/
      DATA ARRFD3(16)/-0.2080 5255 294        D   -5/
      DATA ARRFD3(17)/-0.1047 9767 478        D   -5/
      DATA ARRFD3(18)/-0.5315 2738 02         D   -6/
      DATA ARRFD3(19)/-0.2694 0611 78         D   -6/
      DATA ARRFD3(20)/-0.1374 8787 49         D   -6/
      DATA ARRFD3(21)/-0.7023 0888 7          D   -7/
      DATA ARRFD3(22)/-0.3595 4394 2          D   -7/
      DATA ARRFD3(23)/-0.1851 0612 6          D   -7/
      DATA ARRFD3(24)/-0.9502 3937            D   -8/
      DATA ARRFD3(25)/-0.4918 4811            D   -8/
      DATA ARRFD3(26)/-0.2537 1950            D   -8/
      DATA ARRFD3(27)/-0.1315 1532            D   -8/
      DATA ARRFD3(28)/-0.6835 168             D   -9/
      DATA ARRFD3(29)/-0.3538 244             D   -9/
      DATA ARRFD3(30)/-0.1853 182             D   -9/
      DATA ARRFD3(31)/-0.9589 83              D  -10/
      DATA ARRFD3(32)/-0.5040 83              D  -10/
      DATA ARRFD3(33)/-0.2622 38              D  -10/
      DATA ARRFD3(34)/-0.1372 55              D  -10/
      DATA ARRFD3(35)/-0.7234 0               D  -11/
      DATA ARRFD3(36)/-0.3742 9               D  -11/
      DATA ARRFD3(37)/-0.2005 9               D  -11/
      DATA ARRFD3(38)/-0.1026 9               D  -11/
      DATA ARRFD3(39)/-0.5551                 D  -12/
      DATA ARRFD3(40)/-0.2857                 D  -12/
      DATA ARRFD3(41)/-0.1520                 D  -12/
      DATA ARRFD3(42)/-0.811                  D  -13/
      DATA ARRFD3(43)/-0.410                  D  -13/
      DATA ARRFD3(44)/-0.234                  D  -13/
      DATA ARRFD3(45)/-0.110                  D  -13/
      DATA ARRFD3(46)/-0.67                   D  -14/
      DATA ARRFD3(47)/-0.30                   D  -14/
      DATA ARRFD3(48)/-0.19                   D  -14/
      DATA ARRFD3(49)/-0.9                    D  -15/
      DATA ARRFD3(50)/-0.5                    D  -15/
      DATA ARRFD3(51)/-0.3                    D  -15/
      DATA ARRFD3(52)/-0.1                    D  -15/
      DATA ARRFD3(53)/-0.1                    D  -15/
      DATA ZERO,ONE,TWO/ 0.0 D 0 , 1.0 D 0 , 2.0 D 0/
      DATA THREE,FORTY2,FIFTY/ 3.0 D 0 , 42.0 D 0 , 50.0 D 0/
      DATA GAM2P5/0.1329 3403 8817 9137 D 1/
      DATA TWOE/5.4365 6365 6918 0905 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/13,23,53/
      DATA XMIN1,XMIN2/-35.7D00,-708.394D00/
      DATA XHIGH1,XHIGH2/7.45467D7,3.8392996D205/
C
C   Start calculation
C
      X=XVALUE
C
C   Test for error condition
C
      IF ( X .GT. XHIGH2 ) THEN
         PRINT*,'** ERROR ** - X TOO LARGE FOR FDP0P5'
         FDP0P5 = ZERO
         RETURN
      ENDIF    
C
C   Code for x < -1
C
      IF ( X .LT. -ONE ) THEN
         IF ( X .GT. XMIN1 ) THEN
            EXPX = EXP(X)
            T = TWOE * EXPX - ONE
            FDP0P5 = EXPX * CHEVAL ( NTERM1 , ARRFD1 , T )
         ELSE
            IF ( X .LT. XMIN2 ) THEN
               FDP0P5 = ZERO
            ELSE
               FDP0P5 = EXP(X)
            ENDIF
         ENDIF
      ELSE
C
C   Code for -1 <= x <= 2
C
         IF ( X .LE. TWO ) THEN
            T = ( TWO * X - ONE ) / THREE
            FDP0P5 = CHEVAL ( NTERM2 , ARRFD2 , T )
         ELSE
C
C   Code for x > 2
C
            FDP0P5 = X * SQRT(X) / GAM2P5
            IF ( X .LE. XHIGH1 ) THEN 
               XSQ = X * X
               T = ( FIFTY - XSQ ) / ( FORTY2 + XSQ )
               CHV = CHEVAL ( NTERM3 , ARRFD3 , T )
               FDP0P5 = FDP0P5 * ( ONE + CHV / XSQ )
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END



      DOUBLE PRECISION FUNCTION FDP1P5(XVALUE)
C
C   DESCRIPTION:
C
C      This function computes the Fermi-Dirac function of
C      order 3/2, defined as
C
C                     Int{0 to inf} t**(3/2) / (1+exp(t-x)) dt
C         FDP1P5(x) = -----------------------------------------
C                                 Gamma(5/2)
C
C      The function uses Chebyshev expansions which are given to
C      16 decimal places for x <= 2, but only 10 decimal places 
C      for x > 2.
C
C
C   ERROR RETURNS:
C    
C      If XVALUE too large and positive, the function value
C      will overflow. An error message is printed and the function
C      returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERMS1 - INTEGER - The number of terms used from the array
C                          ARRFD1. The recommended value is such that
C                               ABS(ARRFD1(NTERMS1)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 12.
C
C      NTERMS2 - INTEGER - The number of terms used from the array
C                          ARRFD2. The recommended value is such that
C                               ABS(ARRFD2(NTERMS2)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 22.
C
C      NTERMS3 - INTEGER - The number of terms used from the array
C                          ARRFD3. The recommended value is such that
C                               ABS(ARRFD3(NTERMS3)) < EPS/10
C                          subject to 1 <= NTERMS3 <= 33.
C
C      XMIN1 - REAL - The value of x below which
C                         FDP1P5(x) = exp(x)
C                     to machine precision. The recommended value
C                     is   2.5*LN(2) + LN(EPSNEG)
C
C      XMIN2 - REAL - The value of x below which
C                         FDP1P5(x) = 0.0 
C                     to machine precision. The recommended value
C                     is    LN ( XMIN )
C
C      XHIGH1 - REAL - The value of x above which
C                         FDP1P5(x) = x**(5/2)/GAMMA(7/2) 
C                     to machine precision. The recommended value
C                     is   pi * SQRT(1.6/EPS)
C
C      XHIGH2 - REAL - The value of x above which FDP1P5 would 
C                      overflow. The reommended value is
C                              (3.233509*XMAX)**(2/5)
C
C      For values of EPS, EPSNEG, and XMIN the user should refer to the
C      paper by Cody in ACM. Trans. Math. Soft. Vol. 14 (1988) p303-311.
C   
C      This code is provided with single and double precision values
C      of the machine-dependent parameters, suitable for machines
C      which satisfy the IEEE floating-point standard.
C
C
C   AUTHOR:
C          DR. ALLAN MACLEOD,
C          DEPT. OF MATHEMATICS AND STATISTICS,
C          UNIVERSITY OF PAISLEY,
C          HIGH ST.,
C          PAISLEY,
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:
C                 21 NOVEMBER, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION
     1     ARRFD1(0:12),ARRFD2(0:22),ARRFD3(0:55),
     2     CHEVAL,CHV,EXPX,FIFTY,FORTY2,
     3     GAM3P5,ONE,T,THREE,TWO,TWOE,X,XHIGH1,
     4     XHIGH2,XMIN1,XMIN2,XSQ,XVALUE,ZERO
      DATA ARRFD1/1.9406 5492 1037 8650  D   0,
     1           -0.2878 6747 5518 043   D  -1,
     2            0.8509 1579 5231 3     D  -3,
     3           -0.3327 8452 5669       D  -4,
     4            0.1517 1202 058        D  -5,
     5           -0.7622 0087 4          D  -7,
     6            0.4095 5489            D  -8,
     7           -0.2311 964             D  -9,
     8            0.1355 37              D -10,
     9           -0.8187                 D -12,
     X            0.507                  D -13,
     1           -0.32                   D -14,
     2            0.2                    D -15/
      DATA ARRFD2( 0)/ 3.5862 2516 1563 4306 D   0/
      DATA ARRFD2( 1)/ 1.8518 2900 5626 5751 D   0/
      DATA ARRFD2( 2)/ 0.4612 3491 0241 7150 D   0/
      DATA ARRFD2( 3)/ 0.5793 0397 6126 881  D  -1/
      DATA ARRFD2( 4)/ 0.1704 3790 5548 75   D  -2/
      DATA ARRFD2( 5)/-0.3970 5201 2249 6    D  -3/
      DATA ARRFD2( 6)/-0.7070 2491 890       D  -5/
      DATA ARRFD2( 7)/ 0.7659 9748 792       D  -5/
      DATA ARRFD2( 8)/-0.1857 8113 33        D  -6/
      DATA ARRFD2( 9)/-0.1832 2379 56        D  -6/
      DATA ARRFD2(10)/ 0.1392 4949 5         D  -7/
      DATA ARRFD2(11)/ 0.4670 2027           D  -8/
      DATA ARRFD2(12)/-0.6671 984            D  -9/
      DATA ARRFD2(13)/-0.1161 292            D  -9/
      DATA ARRFD2(14)/ 0.2844 38             D -10/
      DATA ARRFD2(15)/ 0.2490 6              D -11/
      DATA ARRFD2(16)/-0.1143 1              D -11/
      DATA ARRFD2(17)/-0.279                 D -13/
      DATA ARRFD2(18)/ 0.439                 D -13/
      DATA ARRFD2(19)/-0.14                  D -14/
      DATA ARRFD2(20)/-0.16                  D -14/
      DATA ARRFD2(21)/ 0.1                   D -15/
      DATA ARRFD2(22)/ 0.1                   D -15/
      DATA ARRFD3( 0)/12.1307 5817 3688 4627  D   0/
      DATA ARRFD3( 1)/-0.1547 5011 1128 7255  D   0/
      DATA ARRFD3( 2)/-0.7390 0738 8850 999   D  -1/
      DATA ARRFD3( 3)/-0.3072 3537 7959 258   D  -1/
      DATA ARRFD3( 4)/-0.1145 4857 9330 328   D  -1/
      DATA ARRFD3( 5)/-0.4056 7636 8095 39    D  -2/
      DATA ARRFD3( 6)/-0.1398 0158 3732 27    D  -2/
      DATA ARRFD3( 7)/-0.4454 9018 1015 3     D  -3/
      DATA ARRFD3( 8)/-0.1173 9461 1270 4     D  -3/
      DATA ARRFD3( 9)/-0.1484 0898 0093       D  -4/
      DATA ARRFD3(10)/ 0.1188 9515 4223       D  -4/
      DATA ARRFD3(11)/ 0.1464 7695 8178       D  -4/
      DATA ARRFD3(12)/ 0.1132 2874 1730       D  -4/
      DATA ARRFD3(13)/ 0.7576 2292 948        D  -5/
      DATA ARRFD3(14)/ 0.4712 0400 466        D  -5/
      DATA ARRFD3(15)/ 0.2813 2628 202        D  -5/
      DATA ARRFD3(16)/ 0.1637 0517 341        D  -5/
      DATA ARRFD3(17)/ 0.9351 0762 72         D  -6/
      DATA ARRFD3(18)/ 0.5278 6892 10         D  -6/
      DATA ARRFD3(19)/ 0.2951 0798 70         D  -6/
      DATA ARRFD3(20)/ 0.1638 6001 90         D  -6/
      DATA ARRFD3(21)/ 0.9052 0540 9          D  -7/
      DATA ARRFD3(22)/ 0.4977 5697 5          D  -7/
      DATA ARRFD3(23)/ 0.2729 5586 3          D  -7/
      DATA ARRFD3(24)/ 0.1492 1458 5          D  -7/
      DATA ARRFD3(25)/ 0.8142 0359            D  -8/
      DATA ARRFD3(26)/ 0.4434 9200            D  -8/
      DATA ARRFD3(27)/ 0.2411 6032            D  -8/
      DATA ARRFD3(28)/ 0.1310 5018            D  -8/
      DATA ARRFD3(29)/ 0.7109 736             D  -9/
      DATA ARRFD3(30)/ 0.3856 721             D  -9/
      DATA ARRFD3(31)/ 0.2089 529             D  -9/
      DATA ARRFD3(32)/ 0.1131 735             D  -9/
      DATA ARRFD3(33)/ 0.6127 85              D -10/
      DATA ARRFD3(34)/ 0.3314 48              D -10/
      DATA ARRFD3(35)/ 0.1794 19              D -10/
      DATA ARRFD3(36)/ 0.9695 3               D -11/
      DATA ARRFD3(37)/ 0.5246 3               D -11/
      DATA ARRFD3(38)/ 0.2834 3               D -11/
      DATA ARRFD3(39)/ 0.1532 3               D -11/
      DATA ARRFD3(40)/ 0.8284                 D -12/
      DATA ARRFD3(41)/ 0.4472                 D -12/
      DATA ARRFD3(42)/ 0.2421                 D -12/
      DATA ARRFD3(43)/ 0.1304                 D -12/
      DATA ARRFD3(44)/ 0.707                  D -13/
      DATA ARRFD3(45)/ 0.381                  D -13/
      DATA ARRFD3(46)/ 0.206                  D -13/
      DATA ARRFD3(47)/ 0.111                  D -13/
      DATA ARRFD3(48)/ 0.60                   D -14/
      DATA ARRFD3(49)/ 0.33                   D -14/
      DATA ARRFD3(50)/ 0.17                   D -14/
      DATA ARRFD3(51)/ 0.11                   D -14/
      DATA ARRFD3(52)/ 0.5                    D -15/
      DATA ARRFD3(53)/ 0.3                    D -15/
      DATA ARRFD3(54)/ 0.1                    D -15/
      DATA ARRFD3(55)/ 0.1                    D -15/
      DATA ZERO,ONE,TWO/ 0.0 D 0 , 1.0 D 0 , 2.0 D 0/
      DATA THREE,FORTY2,FIFTY/ 3.0 D 0 , 42.0 D 0 , 50.0 D 0/
      DATA GAM3P5/0.3323 3509 7044 7843 D 1/
      DATA TWOE/5.4365 6365 6918 0905 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/12,22,55/
      DATA XMIN1,XMIN2/-35.004D0,-708.396418D0/
      DATA XHIGH1,XHIGH2/166674733.2D0,3.204467D123/
C
C   Start calculation
C
      X=XVALUE
C
C   Test for error condition
C
      IF ( X .GT. XHIGH2 ) THEN
         PRINT*,'** ERROR ** - X TOO LARGE FOR FDP1P5'
         FDP1P5 = ZERO
         RETURN
      ENDIF    
C
C   Code for x < -1
C
      IF ( X .LT. -ONE ) THEN
         IF ( X .GT. XMIN1 ) THEN
            EXPX = EXP(X)
            T = TWOE * EXPX - ONE
            FDP1P5 = EXPX * CHEVAL ( NTERM1 , ARRFD1 , T )
         ELSE
            IF ( X .LT. XMIN2 ) THEN
               FDP1P5 = ZERO
            ELSE
               FDP1P5 = EXP(X)
            ENDIF
         ENDIF
      ELSE
C
C   Code for -1 <= x <= 2
C
         IF ( X .LE. TWO ) THEN
            T = ( TWO * X - ONE ) / THREE
            FDP1P5 = CHEVAL ( NTERM2 , ARRFD2 , T )
         ELSE
C
C   Code for x > 2
C
            FDP1P5 = X * X * SQRT(X) / GAM3P5
            IF ( X .LE. XHIGH1 ) THEN 
               XSQ = X * X
               T = ( FIFTY - XSQ ) / ( FORTY2 + XSQ )
               CHV = CHEVAL ( NTERM3 , ARRFD3 , T )
               FDP1P5 = FDP1P5 * ( ONE + CHV / XSQ )
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END



      DOUBLE PRECISION FUNCTION FDP2P5(XVALUE)
C
C   DESCRIPTION:
C
C      This function computes the Fermi-Dirac function of
C      order 5/2, defined as
C
C                     Int{0 to inf} t**(5/2) / (1+exp(t-x)) dt
C         FDP2P5(x) = -----------------------------------------
C                                 Gamma(7/2)
C
C      The function uses Chebyshev expansions which are given to
C      16 decimal places for x <= 2, but only 10 decimal places 
C      for x > 2.
C
C
C   ERROR RETURNS:
C    
C      If XVALUE too large and positive, the function value
C      will overflow. An error message is printed and the function
C      returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERMS1 - INTEGER - The number of terms used from the array
C                          ARRFD1. The recommended value is such that
C                               ABS(ARRFD1(NTERMS1)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 11.
C
C      NTERMS2 - INTEGER - The number of terms used from the array
C                          ARRFD2. The recommended value is such that
C                               ABS(ARRFD2(NTERMS2)) < EPS/10
C                          subject to 1 <= NTERMS1 <= 21.
C
C      NTERMS3 - INTEGER - The number of terms used from the array
C                          ARRFD3. The recommended value is such that
C                               ABS(ARRFD3(NTERMS3)) < EPS/10
C                          subject to 1 <= NTERMS3 <= 39.
C
C      XMIN1 - REAL - The value of x below which
C                         FDP2P5(x) = exp(x)
C                     to machine precision. The recommended value
C                     is   3.5*LN(2) + LN(EPSNEG)
C
C      XMIN2 - REAL - The value of x below which
C                         FDP2P5(x) = 0.0 
C                     to machine precision. The recommended value
C                     is    LN ( XMIN )
C
C      XHIGH1 - REAL - The value of x above which
C                         FDP2P5(x) = x**(7/2)/GAMMA(9/2) 
C                     to machine precision. The recommended value
C                     is   pi * SQRT(35/(12*EPS))
C
C      XHIGH2 - REAL - The value of x above which FDP2P5 would 
C                      overflow. The reommended value is
C                              (11.6317*XMAX)**(2/7)
C
C      For values of EPS, EPSNEG, and XMIN the user should refer to the
C      paper by Cody in ACM. Trans. Math. Soft. Vol. 14 (1988) p303-311.
C   
C      This code is provided with single and double precision values
C      of the machine-dependent parameters, suitable for machines
C      which satisfy the IEEE floating-point standard.
C
C
C   AUTHOR:
C          DR. ALLAN MACLEOD,
C          DEPT. OF MATHEMATICS AND STATISTICS,
C          UNIVERSITY OF PAISLEY,
C          HIGH ST.,
C          PAISLEY,
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl-ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:
C                 21 NOVEMBER, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION
     1     ARRFD1(0:11),ARRFD2(0:21),ARRFD3(0:61),
     2     CHEVAL,CHV,EXPX,FIFTY,FORTY2,
     3     GAM4P5,ONE,T,THREE,TWO,TWOE,X,XHIGH1,
     4     XHIGH2,XMIN1,XMIN2,XSQ,XVALUE,ZERO
      DATA ARRFD1/1.9694 4166 8589 6693  D   0,
     1           -0.1496 9179 4643 492   D  -1,
     2            0.3006 9558 1662 7     D  -3,
     3           -0.8946 2485 950        D  -5,
     4            0.3298 0720 25         D  -6,
     5           -0.1392 3929 8          D  -7,
     6            0.6455 885             D  -9,
     7           -0.3206 23              D -10,
     8            0.1678 3               D -11,
     9           -0.916                  D -13,
     X            0.52                   D -14,
     1           -0.3                    D -15/
      DATA ARRFD2( 0)/ 4.2642 8383 9865 5301  D   0/
      DATA ARRFD2( 1)/ 2.3437 4268 8491 2867  D   0/
      DATA ARRFD2( 2)/ 0.6727 1197 8005 2076  D   0/
      DATA ARRFD2( 3)/ 0.1148 8263 2796 5569  D   0/
      DATA ARRFD2( 4)/ 0.1093 6396 8046 758   D  -1/
      DATA ARRFD2( 5)/ 0.2567 1739 5701 5     D  -3/
      DATA ARRFD2( 6)/-0.5058 8998 3911       D  -4/
      DATA ARRFD2( 7)/-0.7376 2157 74         D  -6/
      DATA ARRFD2( 8)/ 0.7352 9987 58         D  -6/
      DATA ARRFD2( 9)/-0.1664 2173 6          D  -7/
      DATA ARRFD2(10)/-0.1409 2049 9          D  -7/
      DATA ARRFD2(11)/ 0.9949 192             D  -9/
      DATA ARRFD2(12)/ 0.2991 457             D  -9/
      DATA ARRFD2(13)/-0.4013 32              D -10/
      DATA ARRFD2(14)/-0.6354 6               D -11/
      DATA ARRFD2(15)/ 0.1479 3               D -11/
      DATA ARRFD2(16)/ 0.1181                 D -12/
      DATA ARRFD2(17)/-0.524                  D -13/
      DATA ARRFD2(18)/-0.11                   D -14/
      DATA ARRFD2(19)/ 0.18                   D -14/
      DATA ARRFD2(20)/-0.1                    D -15/
      DATA ARRFD2(21)/-0.1                    D -15/
      DATA ARRFD3( 0)/30.2895 6768 5980 2579  D   0/
      DATA ARRFD3( 1)/ 1.1678 9766 4206 0562  D   0/
      DATA ARRFD3( 2)/ 0.6420 5918 0082 1472  D   0/
      DATA ARRFD3( 3)/ 0.3461 7238 6840 7417  D   0/
      DATA ARRFD3( 4)/ 0.1840 8167 9078 1889  D   0/
      DATA ARRFD3( 5)/ 0.9730 9243 5354 509   D  -1/
      DATA ARRFD3( 6)/ 0.5139 7329 2675 393   D  -1/
      DATA ARRFD3( 7)/ 0.2717 0980 1041 757   D  -1/
      DATA ARRFD3( 8)/ 0.1438 3327 1401 165   D  -1/
      DATA ARRFD3( 9)/ 0.7626 4863 9521 55    D  -2/
      DATA ARRFD3(10)/ 0.4050 3695 7672 02    D  -2/
      DATA ARRFD3(11)/ 0.2154 3961 4641 49    D  -2/
      DATA ARRFD3(12)/ 0.1147 5689 9017 77    D  -2/
      DATA ARRFD3(13)/ 0.6120 6223 6928 2     D  -3/
      DATA ARRFD3(14)/ 0.3268 3403 3785 9     D  -3/
      DATA ARRFD3(15)/ 0.1747 1455 2274 2     D  -3/
      DATA ARRFD3(16)/ 0.9348 7845 7860       D  -4/
      DATA ARRFD3(17)/ 0.5006 9221 2553       D  -4/
      DATA ARRFD3(18)/ 0.2683 7382 1846       D  -4/
      DATA ARRFD3(19)/ 0.1439 5719 1251       D  -4/
      DATA ARRFD3(20)/ 0.7727 2440 700        D  -5/
      DATA ARRFD3(21)/ 0.4150 3820 336        D  -5/
      DATA ARRFD3(22)/ 0.2230 5118 261        D  -5/
      DATA ARRFD3(23)/ 0.1199 3697 093        D  -5/
      DATA ARRFD3(24)/ 0.6452 3443 69         D  -6/
      DATA ARRFD3(25)/ 0.3472 8228 81         D  -6/
      DATA ARRFD3(26)/ 0.1869 9642 15         D  -6/
      DATA ARRFD3(27)/ 0.1007 3002 72         D  -6/
      DATA ARRFD3(28)/ 0.5428 0756 1          D  -7/
      DATA ARRFD3(29)/ 0.2926 0782 9          D  -7/
      DATA ARRFD3(30)/ 0.1577 8591 8          D  -7/
      DATA ARRFD3(31)/ 0.8511 0768            D  -8/
      DATA ARRFD3(32)/ 0.4592 2760            D  -8/
      DATA ARRFD3(33)/ 0.2478 5001            D  -8/
      DATA ARRFD3(34)/ 0.1338 0255            D  -8/
      DATA ARRFD3(35)/ 0.7225 103             D  -9/
      DATA ARRFD3(36)/ 0.3902 350             D  -9/
      DATA ARRFD3(37)/ 0.2108 157             D  -9/
      DATA ARRFD3(38)/ 0.1139 122             D  -9/
      DATA ARRFD3(39)/ 0.6156 38              D -10/
      DATA ARRFD3(40)/ 0.3327 81              D -10/
      DATA ARRFD3(41)/ 0.1799 19              D -10/
      DATA ARRFD3(42)/ 0.9728 8               D -11/
      DATA ARRFD3(43)/ 0.5261 7               D -11/
      DATA ARRFD3(44)/ 0.2846 1               D -11/
      DATA ARRFD3(45)/ 0.1539 7               D -11/
      DATA ARRFD3(46)/ 0.8331                 D -12/
      DATA ARRFD3(47)/ 0.4508                 D -12/
      DATA ARRFD3(48)/ 0.2440                 D -12/
      DATA ARRFD3(49)/ 0.1321                 D -12/
      DATA ARRFD3(50)/ 0.715                  D -13/
      DATA ARRFD3(51)/ 0.387                  D -13/
      DATA ARRFD3(52)/ 0.210                  D -13/
      DATA ARRFD3(53)/ 0.114                  D -13/
      DATA ARRFD3(54)/ 0.61                   D -14/
      DATA ARRFD3(55)/ 0.33                   D -14/
      DATA ARRFD3(56)/ 0.18                   D -14/
      DATA ARRFD3(57)/ 0.11                   D -14/
      DATA ARRFD3(58)/ 0.5                    D -15/
      DATA ARRFD3(59)/ 0.3                    D -15/
      DATA ARRFD3(60)/ 0.2                    D -15/
      DATA ARRFD3(61)/ 0.1                    D -15/
      DATA ZERO,ONE,TWO/ 0.0 D 0 , 1.0 D 0 , 2.0 D 0/
      DATA THREE,FORTY2,FIFTY/ 3.0 D 0 , 42.0 D 0 , 50.0 D 0/
      DATA GAM4P5/0.1163 1728 3965 6745 D 2/
      DATA TWOE/5.4365 6365 6918 0905 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/11,21,61/
      DATA XMIN1,XMIN2/-34.3107854D0,-708.396418D0/
      DATA XHIGH1,XHIGH2/254599860.5D0,2.383665D88/
C
C   Start calculation
C
      X=XVALUE
C
C   Test for error condition
C
      IF ( X .GT. XHIGH2 ) THEN
         PRINT*,'** ERROR ** - X TOO LARGE FOR FDP2P5'
         FDP2P5 = ZERO
         RETURN
      ENDIF    
C
C   Code for x < -1
C
      IF ( X .LT. -ONE ) THEN
         IF ( X .GT. XMIN1 ) THEN
            EXPX = EXP(X)
            T = TWOE * EXPX - ONE
            FDP2P5 = EXPX * CHEVAL ( NTERM1 , ARRFD1 , T )
         ELSE
            IF ( X .LT. XMIN2 ) THEN
               FDP2P5 = ZERO
            ELSE
               FDP2P5 = EXP(X)
            ENDIF
         ENDIF
      ELSE
C
C   Code for -1 <= x <= 2
C
         IF ( X .LE. TWO ) THEN
            T = ( TWO * X - ONE ) / THREE
            FDP2P5 = CHEVAL ( NTERM2 , ARRFD2 , T )
         ELSE
C
C   Code for x > 2
C
            FDP2P5 = X * X * X * SQRT(X) / GAM4P5
            IF ( X .LE. XHIGH1 ) THEN 
               XSQ = X * X
               T = ( FIFTY - XSQ ) / ( FORTY2 + XSQ )
               CHV = CHEVAL ( NTERM3 , ARRFD3 , T )
               FDP2P5 = FDP2P5 * ( ONE + CHV / XSQ )
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END



      DOUBLE PRECISION FUNCTION CHEVAL(N,A,T)
C
C   DESCRIPTION:
C
C      This function evaluates a Chebyshev series, using the
C      Clenshaw method with Reinsch modification, as analysed
C      in the paper by Oliver.
C
C
C   INPUT PARAMETERS
C
C       N - INTEGER - The no. of terms in the sequence
C
C       A - REAL ARRAY, dimension 0 to N - The coefficients of
C           the Chebyshev series
C
C       T - REAL - The value at which the series is to be
C           evaluated
C
C
C   REFERENCES
C
C        "An error analysis of the modified Clenshaw method for
C         evaluating Chebyshev and Fourier series" J. Oliver,
C         J.I.M.A., vol. 20, 1977, pp379-391
C
C
C   MACHINE-DEPENDENT CONSTANTS: NONE
C
C
C   INTRINSIC FUNCTIONS USED;
C
C      ABS
C
C
C    AUTHOR:  Dr. Allan J. MacLeod,
C             Dept. of Mathematics and Statistics,
C             University of Paisley ,
C             High St.,
C             PAISLEY,
C             SCOTLAND
C             ( e-mail:  macl-ms0@paisley.ac.uk )
C
C
C   LATEST MODIFICATION:
C                       21 September , 1995
C
C
      INTEGER I,N
      DOUBLE PRECISION 
     1    A(0:N),D1,D2,HALF,T,TEST,TT,TWO,U0,U1,U2,ZERO
      DATA ZERO,HALF/ 0.0 D 0 , 0.5 D 0 /
      DATA TEST,TWO/ 0.6 D 0 , 2.0 D 0 /
      U1 = ZERO
C
C   If ABS ( T )  < 0.6 use the standard Clenshaw method
C
      IF ( ABS( T ) .LT. TEST ) THEN
         U0 = ZERO
         TT = T + T
         DO 100 I = N , 0 , -1
            U2 = U1
            U1 = U0
            U0 = TT * U1 + A( I ) - U2
 100     CONTINUE
         CHEVAL =  ( U0 - U2 ) / TWO
      ELSE
C
C   If ABS ( T )  > =  0.6 use the Reinsch modification
C
         D1 = ZERO
C
C   T > =  0.6 code
C
         IF ( T .GT. ZERO ) THEN
            TT =  ( T - HALF ) - HALF
            TT = TT + TT
            DO 200 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) + D2
               U1 = D1 + U2
 200        CONTINUE
            CHEVAL =  ( D1 + D2 ) / TWO
         ELSE
C
C   T < =  -0.6 code
C
            TT =  ( T + HALF ) + HALF
            TT = TT + TT
            DO 300 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) - D2
               U1 = D1 - U2
 300        CONTINUE
            CHEVAL =  ( D1 - D2 ) / TWO
         ENDIF
      ENDIF
      RETURN
      END
