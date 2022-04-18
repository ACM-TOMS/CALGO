Default line lengths have been changed from 128 to 75.  This change makes
it easier to show some of the features of messy that would not be
illustrated as nicely with the longer line lengths. Below is an example of
a non-stopping error in the integration program diva.

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
tmessy reports error: Stop level=2 Print level=5 Error index=2
At: TN=0., KSTEP=1, with H=2.328306436538696E-10, Error tolerances are too
small.  (Estimated Error) / (Requested Error) for equation 1 is
44408.92097466651.  Tolerance 1 for F(6) = 9.999999999999999E-21. 
Replacing F(6) with 1.421085471189328E-14.
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


You might get the next line, using an i0 format, but there is no similar
format for the following line of reals with 8 sig. digits.  The next output
is for the same data when 4 digits after the decimal are requested.

idat(1:6) is: i1=0 i2=0 i3=17 i4=-17 i5=8922 i6=-8922
rdat(1:6) is: r1=0. r2=0. r3=17.123457 r4=-.17123457 r5=4.1234568E+6
r6=3.1234568E-4
rdat(1:6) is: r1=0. r2=0. r3=17.1235 r4=-0.1712 r5=4123456.7890 r6=.0003
rdat(1:6) is: r1=  0.0000 r2=  0.0000 r3= 17.1235, r4= -0.1712 r5=********
r6=  0.0003

To show how Fortran handles cases with (a possibly) negative 0, here are
results for the same idat(2) and rdat(2) with a Fortran  print statement
  0 -0.00
 messy avoids printing the "-" (unless, printing a vector with entries < 0)
Fortran may not.

rdat as a vector with a starting index of -3:  0.0000000E+0 -0.0000000E+0
-1:  1.7123457E+1 -1.7123457E-1  4.1234568E+6  3.1234568E-4

Next we check the close calls on the rounding.  In order to protect against
a close call giving a bad format we sometimes get more digits and/or a
wider field than would be absolutely necessary.  We are asking for 4
significant figures.  Imagine what you would want to see for the numbers:
9.994, 9.9996, .99994999, .999999999999999, 9.99949999

rdat(1:5) is:  r1=9.994 r2=10.000 r3=0.9999 r4=1.0000 r5= 9.999

An idat as a vector:      0     0    17   -17  8922 -8922     7   123     5
10:   -12    16   -26    13   267    48     9
Then last rdat:  9.9940  9.9996  0.9999  1.0000  9.9995

Integer vector as array:
1:    0   17 8922    7    5   16   13   48
Same with no offset:      0     0    17   -17  8922 -8922     7   123     5
10:   -12    16   -26    13   267    48     9
Same with offset:    0     0    17   -17  8922 -8922     7   123     5
 -1:   -12    16   -26    13   267    48     9

The last vectors could be printed like this as well:
An idat as a vector:
 1:     0     0    17   -17  8922 -8922     7   123     5   -12    16   -26
13:    13   267    48     9
Then last rdat:
1:  9.9940  9.9996  0.9999  1.0000  9.9995

imat: Col  1 Col  2 Col  3 Col  4 Col  5 Col  6 Col  7 Col  8 Col  9 Col 10 
Row 1    -12    -11    -12      3     -2     -6     -9     23      8     -1
Row 2    -12    -12     -1     -4     -7    -10     19      6     -2     -7
Row 3    -13     -5     -6     -8    -10     15      4     -3     -8     35
Row 4     -9     -8     -9    -11     11      2     -4     -8     31     12
Row 5    -10    -10    -11      7      0     -5     -9     27     10      0
      Col 11 Col 12 Col 13 Col 14 Col 15 
Row 1     -7     43     18      4     -4
Row 2     39     16      3     -5     59
Row 3     14      2     -5     55     24
Row 4      1     -6     51     22      6
Row 5     -6     47     20      5     -4

This would be more compact with shorter column labels.  Thus:

imat:  1C  2C  3C  4C  5C  6C  7C  8C  9C 10C 11C 12C 13C 14C 15C 
Row 1 -12 -11 -12   3  -2  -6  -9  23   8  -1  -7  43  18   4  -4
Row 2 -12 -12  -1  -4  -7 -10  19   6  -2  -7  39  16   3  -5  59
Row 3 -13  -5  -6  -8 -10  15   4  -3  -8  35  14   2  -5  55  24
Row 4  -9  -8  -9 -11  11   2  -4  -8  31  12   1  -6  51  22   6
Row 5 -10 -10 -11   7   0  -5  -9  27  10   0  -6  47  20   5  -4

Original, but for the transposed matrix:
       Col 1 Col 2 Col 3 Col 4 Col 5 
Row  1   -12   -12   -13    -9   -10
Row  2   -11   -12    -5    -8   -10
Row  3   -12    -1    -6    -9   -11
Row  4     3    -4    -8   -11     7
Row  5    -2    -7   -10    11     0
Row  6    -6   -10    15     2    -5
Row  7    -9    19     4    -4    -9
Row  8    23     6    -3    -8    27
Row  9     8    -2    -8    31    10
Row 10    -1    -7    35    12     0
Row 11    -7    39    14     1    -6
Row 12    43    16     2    -6    47
Row 13    18     3    -5    51    20
Row 14     4    -5    55    22     5
Row 15    -4    59    24     6    -4

How about just some of this data with fancier labels?

imat:    Earth   Air  Fire Water 
Hydrogen   -12   -11   -12     3
Copper     -12   -12    -1    -4
Iron       -13    -5    -6    -8
Gold        -9    -8    -9   -11
Lead       -10   -10   -11     7

Transpose real matrix printing (6 significant digits).
rmat^T:
          Col 1      Col 2      Col 3      Col 4      Col 5   
Row  1 -12.437500 -12.437500 -12.562500  -8.937500 -10.437500
Row  2 -11.437500 -12.062500  -4.937500  -8.437500 -10.437500
Row  3 -11.562500  -0.937500  -6.437500  -9.437500 -11.062500
Row  4   3.062500  -4.437500  -8.437500 -10.562500   7.062500
Row  5  -2.437500  -7.437500 -10.062500  11.062500  -0.437500
Row  6  -6.437500  -9.562500  15.062500   1.562500  -5.437500
Row  7  -9.062500  19.062500   3.562500  -4.437500  -8.562500
Row  8  23.062500   5.562500  -3.437500  -8.062500  27.062500
Row  9   7.562500  -2.437500  -7.562500  31.062500   9.562500
Row 10  -1.437500  -7.062500  35.062500  11.562500  -0.437500
Row 11  -6.562500  39.062500  13.562500   0.562500  -6.062500
Row 12  43.062500  15.562500   1.562500  -5.562500  47.062500
Row 13  17.562500   2.562500  -5.062500  51.062500  19.562500
Row 14   3.562500  -4.562500  55.062500  21.562500   4.562500
Row 15  -4.062500  59.062500  23.562500   5.562500  -3.562500

Just default real matrix printing (6 significant digits).
rmat:   Col  1     Col  2     Col  3     Col  4     Col  5     Col  6   
Row 1 -12.437500 -11.437500 -11.562500   3.062500  -2.437500  -6.437500
Row 2 -12.437500 -12.062500  -0.937500  -4.437500  -7.437500  -9.562500
Row 3 -12.562500  -4.937500  -6.437500  -8.437500 -10.062500  15.062500
Row 4  -8.937500  -8.437500  -9.437500 -10.562500  11.062500   1.562500
Row 5 -10.437500 -10.437500 -11.062500   7.062500  -0.437500  -5.437500
        Col  7     Col  8     Col  9     Col 10     Col 11     Col 12   
Row 1  -9.062500  23.062500   7.562500  -1.437500  -6.562500  43.062500
Row 2  19.062500   5.562500  -2.437500  -7.062500  39.062500  15.562500
Row 3   3.562500  -3.437500  -7.562500  35.062500  13.562500   1.562500
Row 4  -4.437500  -8.062500  31.062500  11.562500   0.562500  -5.562500
Row 5  -8.562500  27.062500   9.562500  -0.437500  -6.062500  47.062500
        Col 13     Col 14     Col 15   
Row 1  17.562500   3.562500  -4.062500
Row 2   2.562500  -4.562500  59.062500
Row 3  -5.062500  55.062500  23.562500
Row 4  51.062500  21.562500   5.562500
Row 5  19.562500   4.562500  -3.562500

Same, but testing strange labels and indexes.
rmat:
        < -2>      < -1>      < +0>      < +1>      < +2>      < +3>   
-3R: -12.437500 -11.437500 -11.562500   3.062500  -2.437500  -6.437500
-2R: -12.437500 -12.062500  -0.937500  -4.437500  -7.437500  -9.562500
-1R: -12.562500  -4.937500  -6.437500  -8.437500 -10.062500  15.062500
+0R:  -8.937500  -8.437500  -9.437500 -10.562500  11.062500   1.562500
+1R: -10.437500 -10.437500 -11.062500   7.062500  -0.437500  -5.437500
        < +4>      < +5>      < +6>      < +7>      < +8>      < +9>   
-3R:  -9.062500  23.062500   7.562500  -1.437500  -6.562500  43.062500
-2R:  19.062500   5.562500  -2.437500  -7.062500  39.062500  15.562500
-1R:   3.562500  -3.437500  -7.562500  35.062500  13.562500   1.562500
+0R:  -4.437500  -8.062500  31.062500  11.562500   0.562500  -5.562500
+1R:  -8.562500  27.062500   9.562500  -0.437500  -6.062500  47.062500
        <+10>      <+11>      <+12>   
-3R:  17.562500   3.562500  -4.062500
-2R:   2.562500  -4.562500  59.062500
-1R:  -5.062500  55.062500  23.562500
+0R:  51.062500  21.562500   5.562500
+1R:  19.562500   4.562500  -3.562500

Row 2 of rmat: -12.437500 -12.062500  -0.937500  -4.437500  -7.437500
 6:  -9.562500  19.062500   5.562500  -2.437500  -7.062500  39.062500
12:  15.562500   2.562500  -4.562500  59.062500

Using tabs to display portable public integers in messy_ty:
fpprec=16     line_len=75   maxerr=25002  lstop=3       lprint=3     
errcnt=1      dblev=3

Data entered by hand to show what output from the 'DIVA' integration
program would look like.  Even though this output was designed for 128
character lines, it can still be understandable with shorter lines.

T=.06636  H=.02789, NS=13 NSE=0 E=0. ND=15 KQ= 6 7 7
KSTEP=13 T=.0663568416 H=.0278940 LSC=2 EIMIN=.159 EAVE=.0271 KSC=1
SIGMA(8) RQ=79.1000
    I KQ LI    E      EI      EPS          F        
    1  6  1 7.1E-03 1.0E+00 1.0E-09  9.933206800E-01
    2  7  1 2.8E-04 8.1E-03 1.0E-09 -9.382836700E-02
    3  7  1 2.4E-04 6.9E-03 1.0E-09 -4.774820200E-02
    I      High Order Predicted Differences        Rnoise    Stiff  
    1  2.392E-06 -2.223E-08 -8.231E-09 -3.193E-09  1.1E-08  0.00E+00
    2 -1.307E-07  1.103E-09  2.654E-10  2.654E-10  4.4E-09  0.00E+00
    3 -9.342E-08  8.125E-10  2.265E-10  2.265E-10  7.3E-09  0.00E+00
    I     Beta    
    1  4.36246E+00
    2  6.09012E+00
    3  6.09012E+00

T=.09425  H=.03138, NS=14 NSE=0 E=0. ND=16 KQ= 7 7 7
KSTEP=14 T=.0942508693 H=.0313808 LSC=2 EIMIN=.0657 EAVE=.0948 KSC=1
SIGMA(8) RQ=32.9000
The following is just to illustrate $Y...$y actions.
    I KQ LI    E      EI      EPS          F             D6         D7    
    1  7  0 5.1E-03 7.3E-02 1.0E-09  9.881702200E-01 -8.535E-08 -2.466E-08
    2  7  0 2.6E-03 3.7E-02 1.0E-09 -1.246341400E-01 -3.498E-07  6.896E-09
    3  7  0 1.5E-03 2.1E-02 1.0E-09 -6.331711900E-02 -2.503E-07  4.775E-09
    I     D8         D9      Rnoise    Stiff       Beta    
    1  4.624E-09  4.624E-09  9.3E-09  0.00E+00  3.55720E+00
    2  2.973E-09  1.733E-09  3.0E-08  0.00E+00  3.55720E+00
    3  1.885E-09  8.272E-10  3.4E-08  0.00E+00  3.55720E+00

Testing subtle issues with 'F' format.
test1:  .1230 .1240 .1250 .1260 .1270
test2: t=.09424
Test printing column 7 of a sparse matrix.
Col 7: ( 3, 7.4600) ( 9, 0.9460) (23, 1.7800) (40,-9.4000)
Best to put $N before what you want at start of a
line, Strange things happen with long
preamble+short line. Col 7: ( 3, 7.4600)
preamble+short line. Col 7: ( 9, 0.9460)
preamble+short line. Col 7: (23, 1.7800)
preamble+short line. Col 7: (40,-9.4000)

Complex data: zdat(1)=(1.0000,0.), zdat(2)=(.93750,.062500)

Complex vector zdat:
 1: ( 1.0000E+0,0.0000E+0) ( 9.3750E-1,6.2500E-2) ( 8.7500E-1,1.2500E-1)
 4: ( 8.1250E-1,1.8750E-1) ( 7.5000E-1,2.5000E-1) ( 6.8750E-1,3.1250E-1)
 7: ( 6.2500E-1,3.7500E-1) ( 5.6250E-1,4.3750E-1) ( 5.0000E-1,5.0000E-1)
10: ( 4.3750E-1,5.6250E-1) ( 3.7500E-1,6.2500E-1) ( 3.1250E-1,6.8750E-1)
13: ( 2.5000E-1,7.5000E-1) ( 1.8750E-1,8.1250E-1) ( 1.2500E-1,8.7500E-1)
16: ( 6.2500E-2,9.3750E-1) ( 0.0000E+0,1.0000E+0) (-6.2500E-2,1.0625E+0)
19: (-1.2500E-1,1.1250E+0) (-1.8750E-1,1.1875E+0)

Complex matrix zmat:
              Col  1                 Col  2                 Col  3         
Row 1 ( 1.0000E+0,0.0000E+0) ( 6.8750E-1,3.1250E-1) ( 3.7500E-1,6.2500E-1)
Row 2 ( 9.3750E-1,6.2500E-2) ( 6.2500E-1,3.7500E-1) ( 3.1250E-1,6.8750E-1)
Row 3 ( 8.7500E-1,1.2500E-1) ( 5.6250E-1,4.3750E-1) ( 2.5000E-1,7.5000E-1)
Row 4 ( 8.1250E-1,1.8750E-1) ( 5.0000E-1,5.0000E-1) ( 1.8750E-1,8.1250E-1)
Row 5 ( 7.5000E-1,2.5000E-1) ( 4.3750E-1,5.6250E-1) ( 1.2500E-1,8.7500E-1)
              Col  4                 Col  5                 Col  6         
Row 1 ( 6.2500E-2,9.3750E-1) (-2.5000E-1,1.2500E+0) (-5.6250E-1,1.5625E+0)
Row 2 ( 0.0000E+0,1.0000E+0) (-3.1250E-1,1.3125E+0) (-6.2500E-1,1.6250E+0)
Row 3 (-6.2500E-2,1.0625E+0) (-3.7500E-1,1.3750E+0) (-6.8750E-1,1.6875E+0)
Row 4 (-1.2500E-1,1.1250E+0) (-4.3750E-1,1.4375E+0) (-7.5000E-1,1.7500E+0)
Row 5 (-1.8750E-1,1.1875E+0) (-5.0000E-1,1.5000E+0) (-8.1250E-1,1.8125E+0)
              Col  7                 Col  8                 Col  9         
Row 1 (-8.7500E-1,1.8750E+0) (-1.1875E+0,2.1875E+0) (-1.5000E+0,2.5000E+0)
Row 2 (-9.3750E-1,1.9375E+0) (-1.2500E+0,2.2500E+0) (-1.5625E+0,2.5625E+0)
Row 3 (-1.0000E+0,2.0000E+0) (-1.3125E+0,2.3125E+0) (-1.6250E+0,2.6250E+0)
Row 4 (-1.0625E+0,2.0625E+0) (-1.3750E+0,2.3750E+0) (-1.6875E+0,2.6875E+0)
Row 5 (-1.1250E+0,2.1250E+0) (-1.4375E+0,2.4375E+0) (-1.7500E+0,2.7500E+0)
              Col 10         
Row 1 (-1.8125E+0,2.8125E+0)
Row 2 (-1.8750E+0,2.8750E+0)
Row 3 (-1.9375E+0,2.9375E+0)
Row 4 (-2.0000E+0,3.0000E+0)
Row 5 (-2.0625E+0,3.0625E+0)

A user formatted complex number with same format for real and imaginary
parts, then one with both formats specified, then two numbers printed with
a default of 4 significant digits.  z1=( 1.000E+12,-7.300E+00) z2=
(1.000E+12,-7.300E+00) z3=(1.000E+12,-7.300) z4=(7.360,0.)

Just printing part of text.  And some more.

A repeated message, with name="this name" so we can change names.

An example of binary output: nodec=B"10101 01010101 01010101 01010101 &
&        01010101 01010101"
noinc=B"01100 11001100 11001100 11001100 11001100 11001100"
nodecx=B"10010 00110100 01010110 01111000 10011010 10111100"
noincx=B"01110 11101110 11101110 11101110 11101110 11101110"

The same with no leading 0's and hex output:  nodec=Z"1555 55555555"
noinc=Z"CCC CCCCCCCC" nodecx=Z"1234 56789ABC" noincx=Z"EEE EEEEEEEE"

The same with no leading 0's and octal output:  nodec=O"5252525 25252525"
noinc=O"3146314 63146314" nodecx=O"4432126 36115274"
noincx=O"3567356 73567356"

Printing a bit string vector
bitvec:
1: B"10101 01010101 01010101 01010101 01010101 01010101"
2: B"01100 11001100 11001100 11001100 11001100 11001100"
3: B"10010 00110100 01010110 01111000 10011010 10111100"
4: B"01110 11101110 11101110 11101110 11101110 11101110"
Lining up output with vectors longer than the line
length is tricky.
bitvec:
1: B"   10101 01010101 01010101 01010101 01010101 &
&    01010101"
2: B"   01100 11001100 11001100 11001100 11001100 &
&    11001100"
3: B"   10010 00110100 01010110 01111000 10011010 &
&    10111100"
4: B"   01110 11101110 11101110 11101110 11101110 &
&    11101110"
Etc. bitvec:
1: B"   10101 01010101 01010101 01010101 &
&    01010101 01010101"
2: B"   01100 11001100 11001100 11001100 &
&    11001100 11001100"
3: B"   10010 00110100 01010110 01111000 &
&    10011010 10111100"
4: B"   01110 11101110 11101110 11101110 &
&    11101110 11101110"
Etc., without ix (number of bits <= 31) and longer line.
bitvec:
1: B"0000000 00000000 00101010 10101010"
2: B"1010101 01010101 01010101 01010101"
3: B"0000000 00000000 00011001 10011001"
4: B"1001100 11001100 11001100 11001100"
5: B"0000000 00000000 00100100 01101000"
6: B"1010110 01111000 10011010 10111100"
7: B"0000000 00000000 00011101 11011101"
8: B"1101110 11101110 11101110 11101110"

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
sample reports error: Stop level=3 Print level=4 Error index=1
Start of error message from sample what=1.
I want some 3.141592653589793!
Finish with what=2.
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

Finished first error message for sample,  count=1 maxerr=34001

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
sample reports error: Stop level=5 Print level=6 Error index=99
Fatal error called sample with what = 17.
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
tmessy reports error: Stop level=8 Print level=8 Error index=999
This error message illustrated how an error message ends things, tmessy is
done.
Previously there have been 1 error messages.  The most serious had a
stop index of 2, a print index of 5, and an |error index| of 2
$$$$ Fatal error -- Program stopped. $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
