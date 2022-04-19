example03
 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(z)[1./(1+z.^2).^2;exp(i*z)./(1+z.^2)]


a =

    -1


b =

    -1


c =

   1.0000 + 0.0000i   0.0000 + 2.0000i


exact =

    1.5708
    1.1557


RES =

   1.5708 + 0.0000i
   1.1557 - 0.0000i


ERR =

   1.0e-08 *

   0.4955 - 0.0000i
   0.1448 + 0.0000i


NSUB =

    13


FL =

     2


ACC =

   1.0e-15 *

  -0.2220 + 0.0000i
   0.0000 - 0.0694i

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)1./sqrt(abs(x))


a =

     0


b =

    10


exact =

    6.3246


RES =

    6.3246


ERR =

   9.1551e-16


NSUB =

     2


FL =

     2


ACC =

     0

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)1./sqrt(abs(x))


a =

   -10


b =

    10


exact =

   12.6491


RES =

   12.6491


ERR =

   1.3279e-08


NSUB =

   101


FL =

     2


ACC =

  -1.2176e-08

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)1./(sqrt(x).*(1+x))


a =

     0


b =

   Inf


exact =

    3.1416


RES =

    3.1416


ERR =

   3.1402e-16


NSUB =

     2


FL =

     2


ACC =

     0

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)log(x)./(1-x.^2)


a =

     0


b =

     1


exact =

   -1.2337


RES =

   -1.2337


ERR =

   9.7126e-09


NSUB =

    10


FL =

     2


ACC =

  -1.7838e-10

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)exp(-x).*x./(1-exp(-2*x))


a =

     0


b =

   Inf


exact =

    1.2337


RES =

    1.2337


ERR =

   1.3871e-08


NSUB =

     6


FL =

     2


ACC =

   2.6601e-13

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)[exp(-x).*x;exp(-x).*x.^2;exp(-x).*x.^3;exp(-x).*x.^4;exp(-x).*x.^5]


a =

     0


b =

   Inf


exact =

     1
     2
     6
    24
   120


RES =

    1.0000
    2.0000
    6.0000
   24.0000
  120.0000


ERR =

   1.0e-08 *

    0.0000
    0.0000
    0.0083
    0.1471
    0.6827


NSUB =

    12


FL =

     2


ACC =

   1.0e-12 *

         0
    0.0004
    0.0036
    0.0249
    0.1279

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)exp(-x.^2)


a =

  -Inf


b =

   Inf


exact =

    1.7725


RES =

    1.7725


ERR =

   1.7675e-09


NSUB =

     9


FL =

     2


ACC =

   4.4409e-16

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)exp(-x.^2).*cos(x)


a =

     0


b =

   Inf


exact =

    0.6902


RES =

    0.6902


ERR =

   4.1347e-09


NSUB =

     5


FL =

     2


ACC =

   1.1102e-16

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)exp(-x.^2)./(1+x.^2)


a =

     0


b =

     1


exact =

    0.6188


RES =

    0.6188


ERR =

   2.2393e-13


NSUB =

     2


FL =

     2


ACC =

   1.1102e-16

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)exp(-x(1,:).^2/2)./(1+x(2,:).^2)


a =

  -Inf  -Inf


b =

   Inf   Inf


exact =

    7.8748


RES =

    7.8748


ERR =

   5.9588e-09


NSUB =

    17


FL =

     2


ACC =

   4.3521e-14

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)exp(-x(1,:).^2/2)./(1+x(2,:).^2)


a =

   -10   -10


b =

    10    10


exact =

    7.3751


RES =

    7.3751


ERR =

   1.3193e-08


NSUB =

    76


FL =

     2


ACC =

  -2.1263e-11

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)[exp(-x(1,:).^2/2);1./(1+x(2,:).^2)]


a =

   -10   -10


b =

    10    10


exact =

   50.1326
   58.8451


RES =

   50.1326
   58.8451


ERR =

   1.0e-07 *

    0.1048
    0.1473


NSUB =

    56


FL =

     2


ACC =

   1.0e-13 *

    0.4974
   -0.1421

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)exp(-x(1,:).^2/2)./(1+x(2,:).^2).*x(3,:).^10.*(1-x(3,:)).^10


a =

   -10   -10     0


b =

    10    10     1


exact =

   1.9009e-06


RES =

   1.9009e-06


ERR =

   9.8659e-10


NSUB =

     8


FL =

     2


ACC =

  -8.9970e-12

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)[exp(-x(1,:).^2/2);1./(1+x(2,:).^2);x(3,:).^10.*(1-x(3,:)).^10]


a =

   -10   -10     0


b =

    10    10     1


exact =

   50.1326
   58.8451
    0.0001


RES =

   50.1326
   58.8451
    0.0001


ERR =

   1.0e-07 *

    0.0742
    0.1473
    0.0002


NSUB =

   104


FL =

     2


ACC =

   1.0e-13 *

    0.9237
   -0.3553
   -0.0000

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)x.^(-1/2).*(1-x).^(-1/2)


a =

     0


b =

     1


exact =

    3.1416


RES =

    3.1416


ERR =

   3.4623e-14


NSUB =

     2


FL =

     2


ACC =

   1.0214e-14

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)x.^(-2/3).*(1-x).^(-2/3)


a =

     0


b =

     1


exact =

    5.2999

[Warning: amgkq:inf-or-nan] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 267)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',267,0)">line 267</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 153)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',153,0)">line 153</a>)] 
[Warning: amgkq:error-exceeds-tolerance] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 358)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',358,0)">line 358</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 153)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',153,0)">line 153</a>)] 

RES =

    5.2999


ERR =

   2.9496e-06


NSUB =

    40


FL =

    -1


ACC =

  -1.3794e-05

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)x.^(-3/4).*(1-x).^(-3/4)


a =

     0


b =

     1


exact =

    7.4163

[Warning: amgkq:inf-or-nan] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 267)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',267,0)">line 267</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 162)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',162,0)">line 162</a>)] 
[Warning: amgkq:error-exceeds-tolerance] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 358)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',358,0)">line 358</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 162)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',162,0)">line 162</a>)] 

RES =

    7.4158


ERR =

   1.4211e-04


NSUB =

    40


FL =

    -1


ACC =

  -4.7584e-04

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)(sin(x)./x).^2


a =

     0


b =

   Inf


exact =

    1.5708

[Warning: amgkq:insufficient-subregions] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 346)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',346,0)">line 346</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 170)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',170,0)">line 170</a>)] 
[Warning: amgkq:error-exceeds-tolerance] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 358)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',358,0)">line 358</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 170)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',170,0)">line 170</a>)] 

RES =

    1.5708


ERR =

   4.4588e-06


NSUB =

   200


FL =

     0


ACC =

   5.9371e-06

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)(sin(x)./x).^3


a =

     0


b =

   Inf


exact =

    1.1781


RES =

    1.1781


ERR =

   1.4654e-08


NSUB =

   146


FL =

     2


ACC =

  -7.5785e-08

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)(sin(x)./x).^4


a =

     0


b =

   Inf


exact =

    1.0472


RES =

    1.0472


ERR =

   1.4862e-08


NSUB =

    35


FL =

     2


ACC =

  -1.9602e-09

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)[(sum(x.^2,1)<1);(sum(x.^2,1)>1)]


a =

    -1    -1


b =

     1     1


exact =

    3.1416
    0.8584

[Warning: amgkq:insufficient-subregions] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 346)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',346,0)">line 346</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 194)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',194,0)">line 194</a>)] 
[Warning: amgkq:error-exceeds-tolerance] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 358)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',358,0)">line 358</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 194)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',194,0)">line 194</a>)] 

RES =

    3.1417
    0.8583


ERR =

   1.0e-04 *

    0.7012
    0.7012


NSUB =

   400


FL =

     0


ACC =

   1.0e-04 *

    0.6943
   -0.6943

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)sin(x)./x


a =

     0


b =

   Inf


exact =

    1.5708


longarg =

  1×6 <a href="matlab:helpPopup cell" style="font-weight:bold">cell</a> array

    {0×0 double}    {0×0 double}    {[1000]}    {0×0 double}    {0×0 double}    {[0]}

[Warning: amgkq:insufficient-subregions] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 346)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',346,0)">line 346</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 203)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',203,0)">line 203</a>)] 
[Warning: amgkq:error-exceeds-tolerance] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 358)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',358,0)">line 358</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 203)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',203,0)">line 203</a>)] 

RES =

    1.4657


ERR =

    0.0623


NSUB =

        1000


FL =

     0


ACC =

   -0.1051

 

fx =

  <a href="matlab:helpPopup function_handle" style="font-weight:bold">function_handle</a> with value:

    @(x)sin(3*x).*cosh(x).*sinh(x)


a =

    10


b =

    15


exact =

   2.5884e+10

[Warning: amgkq:error-exceeds-tolerance] 
[> In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('amgkq', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m', 358)" style="font-weight:bold">amgkq</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Src/amgkq.m',358,0)">line 358</a>)
  In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('example03', '/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m', 211)" style="font-weight:bold">example03</a> (<a href="matlab: opentoline('/home/trh/SVNTREE/Toms/TestBed/Johnson/Driver/example03.m',211,0)">line 211</a>)] 

RES =

   2.5884e+10


ERR =

   6.9474e-06


NSUB =

    33


FL =

     1


RELACC =

  -1.4433e-15

