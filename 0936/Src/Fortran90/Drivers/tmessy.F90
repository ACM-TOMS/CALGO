program tmessy
! Time-stamp:  <2015-12-10 19:00:37 m>
! This code by Fred T. Krogh at Math a la Carte.
  use messy_, only: messy_ty, messy, rk, OUTPUT_UNIT
  use sample_,only: sample,sample_ty,setup_sample,partial_message,finish_message

  type(messy_ty) :: e
  type(sample_ty) :: s
  integer :: idat(16), imat(5, 15), i, ix(1)
  integer, allocatable :: iw(:)
  real(rk) :: rdat(40), rmat(5,15), c0
  complex(rk) :: zdat(20), zmat(5,10)


!  ! Lots of character data for next to last case.
  character(len=*), parameter :: diva1=&
    & "$N$D4T=$R  H=$R, NS=$I NSE=$I$D3 E=$R ND=$I KQ=$W$B"
  character(len=*), parameter :: diva2="KSTEP=$I$D9 T=$R $D6H=$R LSC=$I EIMIN=$D3$R&
    & EAVE=$R KSC=$I SIGMA($I) RQ=$D6$R$B"
  character(len=*), parameter :: diva3="$H|5RI|3RKQ|3RLI|8CE|8CEI|8CEPS|CF|&
    &44CHigh Order Predicted Differences|9C Rnoise|10CStiff|13CBeta|"
  character(len=*), parameter :: diva4="$FI5$2FI3$3FE8.1$FD$4FE11.3$FE9.1$FE&
    &10.2$FE13.5"
  character(len=*), parameter :: diva5="$H|5RI|3RKQ|3RLI|8CE|8CEI|8CEPS|CF$Y|&
    &11CD$I$y|9C Rnoise|10CStiff|13CBeta|"
  character(len=*), parameter :: diva6="$FI5$2FI3$3FE8.1$FD$Y$FE11.3$y$FE9.1$FE&
    &10.2$FE13.5"

  integer, parameter :: bitvec(8) = [& ! Note only 31 bits per integer
    int(Z"2AAA"),int(Z"55555555"), int(Z"1999"),int(Z"4CCCCCCC"),&
    int(Z"2468"),int(Z"56789ABC"), int(Z"1DDD"),int(Z"6EEEEEEE")]

! ***********************************************************************
! *********************** Start of examples *****************************
! ***********************************************************************
#ifdef MSWin
! Needed to capture output in a file for comparison using checktmessy.
  OPEN(UNIT=OUTPUT_UNIT,FILE=plet_//'newresult',STATUS='replace') 
#endif
  e%ename = "tmessy" ! Setting the package name for error messages.
  e%fpprec = min(32,e%fpprec) ! So results with NAG compiler will match others
  e%line_len = 75
  s%e=e ! So messages from sample, use the same default line length.
! Illustrates a very simple text message.
  call messy(e, "Default line lengths have been changed from 128 to 75.&
    &  This change makes it easier to show some of the features of messy&
    & that would not be illustrated as nicely with the longer line lengths.&
    & Below is an example of a non-stopping error in the integration program&
    & diva.")
! Illustrates a non-stopping error message.
  call messy(e,"$E25At: TN=$R, KSTEP=$I, with H=$R, Error tolerances are too&
    & small.  (Estimated Error) / (Requested Error) for equation $I is $R.&
    &  Tolerance $I for F($I) = $R.  Replacing F($I) with $R.", idat=&
    [ 2, 1, 1, 1, 6, 6 ], rdat=[ 0._rk, 2.328306436538696E-10_rk,&
    4.440892097466651E+4_rk, 1.e-20_rk, 1.421085471189328E-14_rk ])

  idat(1:6) = [ 0, -0, 17, -17, 8922, -8922 ]
  c0 = real(idat(3) + idat(4), rk) ! Used to show how a -0 will print.
  rdat(1:6) = [ 0.0_rk, -c0, 17.123456789_rk, -.17123456789_rk,&
    4.123456789e6_rk, 3.123456789e-4_rk ]
! Just illustrates how different numbers appear in the output.
  call messy(e,"$NYou might get the next line, using an i0 format, but there is&
    & no similar format for the following line of reals with 8 sig. digits.&
    &  The next output is for the same data when 4 digits after the decimal&
    & are requested.$N$Nidat(1:6) is: i1=$I i2=$I i3=$I i4=$I i5=$I i6=$I$N&
    &$D8rdat(1:6) is: r1=$R r2=$R r3=$R r4=$R r5=$R r6=$R$B", idat, rdat)
  call messy(e,"$D-4rdat(1:6) is: r1=$R r2=$R r3=$R r4=$R r5=$R r6=$R$B",&
    & rdat=rdat)

  call messy(e, "rdat(1:6) is: r1=$FF8.4 r2=$G r3=$G, r4=$G r5=$G r6=$G",&
    & rdat=rdat)
! Fortran output in some cases will be slightly different.
  write(OUTPUT_UNIT, &
    & '(/"To show how Fortran handles cases with (a possibly) negative 0,&
    & here are"/"results for the same idat(2) and rdat(2) with a Fortran&
    &  print statement"/ i3, f6.2/ " messy avoids printing the ""-"" (unless,&
    & printing a vector with entries < 0)"/&
    &"Fortran may not.")') idat(2), rdat(2)
! Illustrates print of a real vector.
  call messy(e,"$N$D8rdat as a vector with a starting index of -3:$V$O-3$B",&
    & idat, rdat(1:6))

  rdat(1:5) = [ 9.994_rk, 9.9996_rk, .99994999_rk, .999999999999999_rk,&
    & 9.99949999_rk ]
! Next is just to test that rounding of numbers is done properly.
  call messy(e,"$NNext we check the close calls on the rounding.  In order to&
    & protect against a close call giving a bad format we sometimes get more&
    & digits and/or a wider field than would be absolutely necessary.&
    &  We are asking for 4 significant figures.  Imagine what you would want&
    & to see for the numbers: 9.994, 9.9996, .99994999, .999999999999999,&
    & 9.99949999$N$N$D4rdat(1:5) is:  r1=$R r2=$R r3=$R r4=$R r5=$R",&
    & idat, rdat)
! Illustrates output of an integer vector.
  idat(7:16) = [ 7, 123, 5, -12, 16, -26, 13, 267, 48, 9  ]
  call messy(e, "$N$D4An idat as a vector:$WThen last rdat:$V$B",&
    & idat(1:16), rdat(1:5))
! Shows how to print a stride of an array.
  allocate (iw(-10:5))
  iw(-10:5) = idat
  call messy(e, "$NInteger vector as array:$N$W", iw(-10:5:2))
! Shows how one can get desired indexes when first index in array is not 1.  
  call messy(e, "Same with no offset:$W",iw)
  call messy(e, "Same with offset:$W$O-10",iw)
! Shows writing both an integer and a real vector at the same time.
  call messy(e,"$NThe last vectors could be printed like this as well:$N&
    &$D4An idat as a vector:$N$WThen last rdat:$N$V$B", idat(1:16), rdat(1:5))

  imat = reshape( [ (nint(-12.9375_rk+real(i,rk)*2.0_rk**(-mod(i,4))),&
    & i = 1, 75)], shape(imat)) !  Numbers should match on different machines.
! Illustrates printing an integer matrix.
  call messy(e, "$Nimat:$M", imat=imat)
! When lines don't hold enough one can tighten things up a bit.
  call messy(e,"$NThis would be more compact with shorter column labels.&
    &  Thus:$N$Nimat:$M$O<1C$O$O", imat=imat)
! Shows printing the transpose of a matrix.
  call messy(e,"$NOriginal, but for the transposed matrix:$M$T", imat=imat)
! Shows how you can make any labels you like.
  call messy(e,"$NHow about just some of this data with fancier labels?$N&
    &$Nimat:$M$O 5Earth  Air FireWater$O 8HydrogenCopper  Iron    Gold    &
    &Lead    $O", imat=imat(1:5,1:4))

  rmat = reshape( [(-12.9375+real(i,rk)*2.0_rk**(-mod(i,4)), i = 1, 75)],&
    & shape(imat))
! And the transpose of a real matrix.
  call messy(e,"$N$D6Transpose real matrix printing (6 significant digits).$N&
    &rmat^T:$A$T$B", rmat=rmat)
! And just printing a real matrix
  call messy(e,"$N$D6Just default real matrix printing (6 significant digits).&
    &$Nrmat:$A$B", rmat=rmat)
! Same but with different labels.
  call messy(e,"$N$D6Same, but testing strange labels and indexes.$Nrmat:&
    &$A$O-2|11<>$O-3<2R:$O$B", rmat=rmat)
! Just a row of a matrix can be printed.
  ix(1) = 2
  call messy(e,"$N$D6Row $X of rmat:$V$B", idat, rmat(ix(1),1:15), ix=ix)
! Illustrates using tabs for printing out a bunch of values.
  call messy(e,"$NUsing tabs to display portable public integers in messy_ty:&
    &$14T$Nfpprec=$I$Tline_len=$I$Tmaxerr=$I$Tlstop=$I$Tlprint=$I&
    &$Terrcnt=$I$Tdblev=$I", [ e%fpprec, e%line_len, e%maxerr, e%lstop,&
    & e%lprint, e%errcnt, e%dblev ])

! This table example is messy as we need to put data in by hand here.

  call messy(e,"$NData entered by hand to show what output from the 'DIVA'&
    & integration program would look like.  Even though this output was&
    & designed for 128 character lines, it can still be understandable with&
    & shorter lines.")

!e%line_len=128 ! Uncomment this to see what output should look like.

  rdat(1:3) = [ 6.636e-2_rk, 2.789e-2_rk, 0.0_rk ]
  idat(1:6) = [ 13, 0, 15, 6, 7, 7  ]
! Illustrates setting up the tables, first a line prior to the table.  
  call messy(e,diva1, idat(1:6), rdat)
  rdat(1:6) = [ 6.63568416E-02_rk, 2.78940E-02_rk, 1.59E-01_rk, 2.71E-02_rk,&
    &7.91E+01_rk, 1.42500_rk ]
  idat(1:4) = [ 13, 2, 1, 8 ]
! Then another line preceding the table.  
  call messy(e,diva2, idat, rdat)
  e%kdf=min(9, precision(1.0_rk)) ! Set so $FD can be used later.
! Setting up the headings for the table.  
  call messy(e,diva3)
  idat(1:3) = [ 1, 6, 1 ]
  rdat(1:11) = [7.1E-03_rk, 1.0E+00_rk, 1.0E-09_rk, 9.9332068E-01_rk,&
    & 2.392E-06_rk, -2.223E-08_rk, -8.231E-09_rk, -3.193E-09_rk, 1.1E-08_rk,&
    & 0.00E+00_rk, 4.36246E+00_rk ]
! Output a line in the table  
  call messy(e, diva4, idat(1:3), rdat(1:11))
  idat(1:3) = [ 2, 7, 1 ]
  rdat(:11) = [ 2.8E-04_rk, 8.1E-03_rk, 1.0E-09_rk, -9.3828367E-02_rk,&
    & -1.307E-07_rk, 1.103E-09_rk, 2.654E-10_rk, 2.654E-10_rk, 4.4E-09_rk,&
    & 0.00E+00_rk, 6.09012E+00_rk ]
! And another line  
  call messy(e, diva4, idat(1:3), rdat(1:11))
  idat(1:3) = [ 3, 7, 1 ]
  rdat(1:11) = [ 2.4E-04_rk, 6.9E-03_rk, 1.0E-09_rk, -4.7748202E-02_rk,&
    & -9.342E-08_rk, 8.125E-10_rk, 2.265E-10_rk, 2.265E-10_rk, 7.3E-09_rk,&
    & 0.00E+00_rk, 6.09012E+00_rk ]
! And another  
  call messy(e, diva4, idat(1:3), rdat(1:11))
! Table must be ended like this with no output except for lines of the table.
  call messy(e, "$B")

  rdat(1:3) = [  9.425E-02_rk, 3.138E-02_rk, 0.00E+00_rk ]
  idat(1:6) = [ 14, 0, 16, 7, 7, 7 ]
! Illustrates another way of setting up things for a table.  
  call messy(e,diva1, idat(1:6), rdat)

  rdat(1:6) = [ 9.42508693E-02_rk, 3.13808E-02_rk, 6.57E-02_rk, 9.48E-02_rk,&
    & 3.29E+01_rk, 1.42500E+00_rk ]
  idat(1:4) = [14, 2, 1, 8 ]
! Etc.  
  call messy(e,diva2, idat, rdat)
  e%kdf=min(9, precision(1.0_rk))
  call messy(e, "The following is just to illustrate $$Y...$$y actions.")
! Like using diva3, but table size could be changed at run time.
  call messy(e,diva5, ix=[4], idat=[6,7,8,9])
  idat(1:3) = [ 1, 7, 0 ]
  rdat(1:11) = [5.1E-03_rk, 7.3E-02_rk, 1.0E-09_rk, 9.8817022E-01_rk,&
    & -8.535E-08_rk, -2.466E-08_rk, 4.624E-09_rk, 4.624E-09_rk, 9.3E-09_rk,&
    & 0.00E+00_rk, 3.55720E+00_rk ]
! Like using diva4, but one could now number of things could change.  
  call messy(e, diva6, idat(1:3), rdat(1:11), ix=[4])
  idat(1:3) = [ 2, 7, 0 ]
  rdat(1:11) = [ 2.6E-03_rk, 3.7E-02_rk, 1.0E-09_rk, -1.2463414E-01_rk,&
    & -3.498E-07_rk, 6.896E-09_rk, 2.973E-09_rk, 1.733E-09_rk, 3.0E-08_rk,&
    & 0.00E+00_rk, 3.55720E+00_rk ]
  call messy(e, diva6, idat(1:3), rdat(1:11), ix=[4])
  idat(1:3) = [3, 7, 0 ]
  rdat(1:11) = [ 1.5E-03_rk, 2.1E-02_rk, 1.0E-09_rk, -6.3317119E-02_rk,&
    & -2.503E-07_rk, 4.775E-09_rk, 1.885E-09_rk, 8.272E-10_rk, 3.4E-08_rk,&
    & 0.00E+00_rk, 3.55720E+00_rk ]
  call messy(e, diva6, idat(1:3), rdat(1:11), ix=[4])
! And the end of this table.  Tables are not easy, but give very nice output.
  call messy(e, "$B")

  rdat(1:5) = [ .123_rk, .124_rk, .125_rk, .126_rk, .127_rk ]
! Checking that 'F' format does the right thing.  
  call messy (e, "$NTesting subtle issues with 'F' format.$Ntest1: $D4$V",&
    &rdat=rdat(1:5))
  rdat(1) = .0942408_rk
  call messy (e, "test2: t=$D4$R", rdat=rdat)
! Illustrates printing a column of a sparse matrix.
  call messy(e, "Test printing column 7 of a sparse matrix.$NCol $X:$S",&
    idat=[3, 9, 23, 40], rdat=[7.46_rk, .946_rk, 1.78_rk, -9.40_rk], ix=[7])
! Illustrates problems with using a long preamble for sparse output.
  call messy(e, "$L50Best to put $$N before what you want at start of a line,&
    & Strange things happen with long preamble+short line. Col $X:$S$B",&
    & idat=[3, 9, 23,40], rdat=[7.46_rk, .946_rk, 1.78_rk, -9.40_rk], ix=[7])

  zdat(1:20) =   [ (cmplx(1.0_rk-.0625_rk*real(i,rk),&
    & .0625_rk*real(i,rk),rk), i = 0, 19) ]
! Illustrates printing of a complex data.
  call messy(e, "$N$D5Complex data: zdat(1)=$ZR, zdat(2)=$ZR", zdat=zdat)
! Illustrates printing of a complex vector.
  call messy(e, "$N$D5Complex vector zdat:$N$ZV", zdat=zdat)
  zmat = reshape( [ (cmplx(1.0_rk-.0625_rk*real(i,rk),&
    & .0625_rk*real(i,rk),rk), i = 0, 49) ], shape(zmat))
! And then a complex matrix.  
  call messy(e,"$N$D5Complex matrix zmat:$N$ZA", zmat = zmat)

  zdat(1) = (1.e12_rk, -7.3_rk)
! Shows how to get the same format for real and imaginary parts.  
  call messy(e, "$NA user formatted complex number with same format for real&
    & and imaginary parts, then one with both formats specified, then two&
    & numbers printed with a default of 4 significant digits.  z1=$ZFE10.3&
    & z2=$ZFE9.3,E10.3 $D4z3=$ZR z4=$ZR",&
    & zdat = [ zdat(1), zdat(1), zdat(1), (7.36_rk, 0.0_rk) ])

! Illustrates how you can just get just part of a message based on e%dblev. 
  call messy(e,"$NJust printing part of text.$K2  And some more.$K4&
    &  But not all of it!")
! Illustrates how a message can be repeated, except for part of text changing.
  call messy(e,"$NA repeated message, with name=$P so we can change names.",&
    & ptext='"this name"')
! Illustrates binary output of a variable with 45 bits.
  call messy(e,"$NAn example of binary output: nodec=$QBI$H noinc=$QBI&
    &$H nodecx=$QBI$H noincx=$QBI$H", idat=bitvec, ix=[45])
! Same but leading 0's are not printed, and output is in hex.
  call messy(e, "$NThe same with no leading 0's and hex output:  nodec=$QZI$H&
    & noinc=$QZI$H nodecx=$QZI$H noincx=$QZI$H", idat=bitvec, ix=[-45])
! And yes Octal output is possible too.
  call messy(e, "$NThe same with no leading 0's and octal output:  nodec=$QOI$H&
    & noinc=$QOI$H nodecx=$QOI$H noincx=$QOI$H", idat=bitvec, ix=[-45])
! Illustrates printing a vector of bit strings.
  call messy(e, "$NPrinting a bit string vector$Nbitvec:$QBV",&
    & idat=bitvec,ix=[45])
! Illustrates how bit strings can be lined up -- not easy.   
  call messy(e, "$L50Lining up output with vectors longer than the line length&
    & is tricky.$Nbitvec:$QBV",idat=bitvec,ix=[45])
  call messy(e, "$L40Etc. bitvec:$QBV$B",idat=bitvec,ix=[45])
! Another binary vector.  
  call messy(e, "Etc., without ix (number of bits <= 31) and longer line.$N&
    &bitvec:$QBV", idat=bitvec)
  s%what = setup_sample  !Code to show various interactions.
! This final section illustrates how different levels of a library package can
! output messages that are separate for different parts of the package.  
  call sample(s)
  s%what = partial_message
  call sample(s) ! Uncomment statement below to see how messages can get mixed.
!call messy(e, "$NStarted error message for sample, count=$I maxerr=$I",&
!  & [ s%e%errcnt, s%e%maxerr ])
  s%what = finish_message
  call sample(s)
  call messy(e,&
    &"Finished first error message for sample,  count=$I maxerr=$I",&
    & [ s%e%errcnt, s%e%maxerr ])
  s%what = 17 ! Used to get fatal errror from sample.

  s%e%lstop = 5 ! Uncomment this to avoid stop below.
  call sample(s)

  call messy(e,"$E88This error message illustrated how an error message ends&
    & things, tmessy is done.", [ 999 ])
end program tmessy
