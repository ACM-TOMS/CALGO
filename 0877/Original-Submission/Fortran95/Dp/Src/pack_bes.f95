! File pack_bes.f95
!
! Author: Masao Kodama
! Address: 21-20, Gakuen 1 chome, Mastue-shi, Shimane-ken, 690-0825 Japan
! Email: mkodama@mable.ne.jp
!
! File pack_bes.f95 contains the present algorithm.
! This algorithm provides a package of subroutines computing the cylindrical
! functions:
! * the Bessel function of the first kind J(x)
! * the Neumann function N(x) {the Bessel function of the second kind Y(x)}
! * the Hankel function of the first kind H(1)(x)
! * the Hankel function of the second kind H(2)(x)
! of a complex order and a nonnegative argument.
! The programs are written in Fortran 95.
!
! Directions for the present algorithm are as follows.
!
! In all the files in archive soft_bes.gz, the following notations are used.
!
! zbessel(znu,xx) =the Bessel function of the first kind J(x).
! zneumann(znu,xx)=the Neumann function N(x), that is, the Bessel function of
!                  the second kind Y(x).
! zhankel1(znu,xx)=the Hankel function of the first kind H(1)(x).
! zhankel2(znu,xx)=the Hankel function of the second kind H(2)(x).
! znu=the complex order
! xx=the nonnegative argument
! zans=the complex value of the invoked cylindrical function.
! kp= the kind type parameter for the real numbers and the complex numbers
!     used in the present algorithm.
! nregion= a region number: an integer which is equal to the suffix n of
!          region Rn defined in Section  2.1.1 of Ref. (1).
!          The case of nregion=0 denotes that the arguments znu, xx are out of
!          range.
! pi= the circle ratio
! zunit= (0,1)
! epsilon0= EPSILON(1._kp)
! epsilon1= MAX(2.22E-16,epsilon0)
! re_znu= REAL(znu,kp)
! aim_znu= AIMAG(znu)
! zconnu= CONJG(znu)
!
! The algorithm consists of the following three Groups (A), (B), (C) of
! subprograms.
!
! Group (A)
!
! SUBROUTINE bessel(znu,xx,zans,info)
! This SUBROUTINE calculates zans.  zans=zbessel(znu,xx).
! This invokes num_region, bes_series, neu_series, bes_recur, han2_temme,
! bes_han_dby, bes_olver, han2_olver, abs2.
!
! SUBROUTINE neumann(znu,xx,zans,info)
! This SUBROUTINE calculates zans.  zans=zneumann(znu,xx).
! This invokes num_region, bes_series, neu_series, bes_recur, han2_temme,
! bes_han_dby, bes_olver, han2_olver, abs2.
!
! SUBROUTINE hankel1(znu,xx,zans,info)
! This SUBROUTINE calculates zans.  zans=zhankel1(znu,xx).
! This invokes hankel2.
!
! SUBROUTINE hankel2(znu,xx,zans,info)
! This SUBROUTINE calculates zans.  zans=zhankel2(znu,xx).
! This invokes num_region, bes_series, neu_series, bes_recur, han2_temme,
! bes_han_dby, han2_olver, abs2.
!
! The arguments znu, xx must be defined before the use of the subroutines.
! zans is the output for the complex value of the invoked function.
! info =  output, an integer; information about the output zans:
!     info= 0: normal output, e_r <= 5E3*epsilon1, e_r=the relative error of the
!              functional value zans.
!     info= 5: low accuracy, 5E3*epsilon1 < e_r <= 2E7*epsilon1.
!     info=10: rough accuracy, e_r > 2E7*epsilon1
!     info=20: (1) There is a possibility of an overflow, and the value of the
!                  function is given HUGE(1._kp).
!              (2) The answer zans is indefinite theoretically.
!                  zbessel(znuit,0) is indefinite for example.
!     info=30  Out of range.  xx < 0 for example.
!
! An increase of a relative error in the neighborhood of a zeros of a
! cylindrical function does not influence the value of argument info [Ref. (1),
! Sec. 3.3].
!
! Group (B)
!
! SUBROUTINE bes_series(znu,xx,zsum,zlogbes,error,info)
! This SUBROUTINE calculates zbessel(znu,xx) using the series expansion when
! nregion=1, 2.
! This is invoked by bessel, neumann, hankel2.
! This invokes abs2, cdlgam.
!
! SUBROUTINE neu_series(znu,xx,zneu,info)
! This SUBROUTINE calculates zneumann(znu,xx) using the series expansion
! when nregion=2.
! This is invoked by bessel, neumann, hankel2.
! This invokes neu_srs_init, abs2.
!
! SUBROUTINE bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
! This SUBROUTINE calculates zhankel1(znu,xx), zhankel2(znu,xx) using Debye's
! expansions when nregion=3.
! This SUBROUTINE calculates zbessel(znu,xx), zhankel2(znu,xx) using Debye's
! expansions when nregion=4.
! This is invoked by bessel, neumann, hankel2.
! This invokes abs2.
!
! SUBROUTINE bes_olver(znu,xx,zbes,info)
! This SUBROUTINE calculates zbessel(znu,xx) using Olver's expansion when
! nregion=5.
! This is invoked by bessel, neumann.
! This invokes abs2, fzeta, aiz, sumaabb.
!
! SUBROUTINE han2_olver(znu,xx,zhan2,info)
! This SUBROUTINE calculates zhankel2(znu,xx) using Olver's expansion when
! nregion=5.
! This is invoked by bessel, neumann, hankel2.
! This invokes abs2, fzeta, aiz, sumaabb.
!
! SUBROUTINE bes_recur(znu,xx,zbes,info)
! This SUBROUTINE calculates zbessel(znu,xx) using the recurrence method
! when nregion=6, 7.
! This is invoked by bessel, neumann, hankel2.
! This invokes cdlgam, abs2.
!
! SUBROUTINE han2_temme(znu,xx,zhan2,info)
! This SUBROUTINE calculates zhankel2(znu,xx) using Temme's algorithm when
! nregion=7.
! This is invoked by bessel, neumann, hankel2.
! This invokes abs2.
!
! Group (C)
!
! INTEGER FUNCTION num_region(znu,xx)
! This SUBROUTINE determines nregion.  nregion=num_region(znu,xx)
! This is invoked by bessel, neumann, hankel2.
! This invokes abs2.
!
! SUBROUTINE neu_srs_init(nn,zeps,xx,zneum,info)
! This is invoked by neu_series.
! This invokes def_bessel, abs2.
!
! SUBROUTINE def_bessel(nn,zeps,xx,bes0,zbes1,info)
! This is invoked by neu_srs_init.
! This invokes abs2.
!
! SUBROUTINE sumaabb(znu2,zeta,zetad,zfacta,zfactb,zsum,info)
! This is invoked by bes_olver, han2_olver.
! This invokes abs2.
!
! SUBROUTINE fzeta(zz,zeta,zetad,info)
! This is invoked by bes_olver, han2_olver.
! This invokes abs2.
!
! SUBROUTINE aiz(ifun,ifac,x0,y0,gair,gaii,ierro)
! This SUBROUTINE calculates the Airy function Ai(za) of the complex argument za
! and its derivative Ai'(za).
! This is invoked by bes_olver, han2_olver.
! This is cited from Algorithm 819 [Ref. (5)].
!
! SUBROUTINE cdlgam(carg,cans,error,lf0)
! This SUBROUTINE calculates the gamma function or the logarithm of the
! gamma function of a complex argument.
! This is invoked by bes_series, bes_recur.
! This is cited from Algorithm 421 [Ref. (6)].
!
! FUNCTION abs2(za)
! This calculates a rough absolute value of the complex argument za.
! This is invoked by bessel, neumann, hankel2, num_region, bes_series,
! neu_series, neu_srs_init, def_bessel, bes_recur, han2_temme, bes_han_dby,
! bes_olver, han2_olver, sumaabb, fzeta.
!
!
! All the above subprograms are the module subprograms of MODULE mod_bes. 
! SUBROUTINEs bessel, neumann, zhankel1, hankel2 of Group (A) can be invoked by
! a user directly. If one of SUBROUTINEs bessel, neumann, hankel2 is invoked by
! the user, this subroutine invokes FUNCTION num_region first and determines 
! nregion and selects one or two of subroutines of Group (B) according to
! nregion. In addition, subprograms of Group (C) help the numerical calculation
! in subroutines of Group (B). Using Eqs. (40)-(42) of Ref. (1) if necessary,
! the subroutine invoked by the user determines its functional value zans.
! On the other hand, zhankel1(znu,xx) is calculated by the formula
! zhankel1(znu,xx)=CONJG(zhankel2(zconnu,xx)), where zhankel2(zconnu,xx) is
! calculated by SUBROUTINE hankel2.
!
!
! Notices
! (1) In this algorithm, the default value of kp is kp=KIND(1D0).
! We can redefine kp in MODULE mod_bes in this file.
! An arbitrary kp can be used if the using processor accepts the kp.
! The type of the arguments znu, xx, zans must be agree with the type
! determined by kp.
!
! (2) The variable epsilon1 decides the relative errors of the cylindrical
! functions computed by the present algorithm, and determines M in Eq. (9),
! K in Eq. (10), M1, M2, M3, M4 in Eq. (14) of Ref. (1), so that the truncation
! relative errors are approximately epsilon1.
! In addition, this algorithm is designed so as to work best when
! epsilon1=2.22E-16. Hence, the precision of the results calculated by this
! algorithm is not always improved even if epsilon0 < 2.22E-16.  Hence it is
! desired to set kp=KIND(1D0) first. If overflows occur when kp=KIND(1D0), it is
! very effective in avoiding the overflows to use quadruple precision, because
! the largest real number of higher precision becomes larger.
!
! (3) Characteristics of the intrinsic functions of Fortran such as EXP, SIN,
! REAL, CMPLX vary from processor to processor, and usually do not so much
! influence the computed results.  However, there are cases where the
! characteristics of the intrinsic functions drastically influence the computed
! results.  For example, am_arg denotes the maximum allowable value of ABS(arg),
! where arg is the argument of SIN(arg).  The value of am_arg depends on the
! processor. The intrinsic functions of Fortran have many characteristic
! constants like the constant am_arg. The present algorithm cannot know the
! values of the characteristic constants  beforehand; therefore there is a
! possibility that the implementation of the  present algorithm stops when
! ABS(arg)>am_arg.
!
! (4) The real numbers in many processors include denormalized numbers
! [Ref. (3), p. 173]. The present algorithm is written under the condition
! that the denormalized numbers are included. If the processor does not have
! the denormalized numbers, there is a possibility that the present algorithm
! gives mistaken functional values even at info=0 when an underflow occurs.
! We can know by executing the following program whether a processor has the
! denormalized numbers.
!   PRINT *,TINY(1D0)/2;   END
! If the output is 0.0, the processor does not use the denormalized numbers.
!
! (5) Even if an answer zans underflows, the present algorithm does not inform
! the user of the underflow.  In the case of the underflow, the value of info is
! indefinite, and the number of the significant figures of the answer zans
! decreases gradually according to the degree of the underflow.
!
!
! The following references are cited in all the files in the archive
! soft_bes.gz.
!
! References
!
! (1) Masao Kodama, "Algorithm XXX: A subroutine package for cylindrical
! functions of complex order and nonnegative argument," ACM Transactions on
! Mathematical Software.
!
! (2) Milton Abramowitz, and Irene A. Stegun, "Handbook of Mathematical
! Functions with Formulas, Graphs, and Mathematical Tables," Dover Publications,
! Inc., New York, November 1970.
!
! (3) Michael Metcalf and John Reid, "Fortran 90/95 Explained (Second Edition),"
! Oxford University Press, 1999.
!
! (4) N. M. Temme, "On the numerical evaluation of the modified Bessel function
! of the third kind," Journal of Computational Physics, vol. 19, pp. 324-337,
! 1975.
!
! (5) Amparo Gil, Javier Segura and Nico M. Temme,
! "Algorithm 819   AIZ, BIZ: two Fortran 77 routines for the computation of
! complex Airy functions," ACM Transactions on Mathematical Software,
! vol. 28, no. 3, September 2002, pp. 325-336.
!
! (6) Hirondo Kuki, "Algorithm 421   Complex gamma function with error control
! [S14]," Communications of the ACM, vol. 15, no. 4, pp. 271-272, April 1972.
!
! (7) D. E. Amos, "Algorithm 644   A portable package for Bessel functions of a
! complex argument and nonnegative order," ACM Transactions on Mathematical
! Software, vol. 12, no. 3, pp. 265-273, Sept. 1986.
!
! (8) Frank W. J. Olver, "Asymptotics and Special Functions," Academic Press,
! New York, chap. 11, 1975.
!
! (9) M. Goldstein and R. M. Thaler, "Recurrence technic for the calculation of
! Bessel functions," Mathematical Tables and Other Aids for Computation, vol.
! 13, pp. 102-108, 1959.
!
! (10) G. N. Watson, "A Treatise on the Theory of Bessel Functions (Second
! Edition)," Cambridge University Press, pp. 235-241, 262-268, March 1944.
!
! (11) M. Kodama, M. Yamasato, and S. Yamashiro, "Numerical calculation of the
! Neumann function Nnu(x) of complex order nu," IEICE Trans. Fundamentals, vol.
! E78-A, no. 6, pp. 727-736, June 1995.
!
! (12) Masao Kodama, "Numerical calculation of the Bessel function of complex
! order using the recurrence method," IEICE Trans. Fundamentals, vol. E78-A,
! no. 4, pp. 506-516, April 1995.
!
! (13) Mohd Abdur Rashid and Masao Kodama, "Numerical calculation of cylindrical
! functions of complex order using Debye's asymptotic series," IEICE Trans.
! Fundamentals, vol. E83-A, no. 12, pp. 2664-2671, Dec. 2000.
!
! (14) Masao Kodama and Kengo Taira, "Polynomials approximating complex
! functions," IEICE Trans. Fundamentals, vol. E80-A, no. 4, pp. 778-781, April
! 1997.
!
! (15) B. R. Fabijonas, "Algorithm 838: Airy functions," ACM Transactions on
! Mathematical Software, vol. 30, no. 4, pp. 491-501, Dec. 2004.
!
! (16) I. J. Thompson and A. R. Barnett, "Coulomb and Bessel functions of
! complex arguments and order," Journal of Computational Physics, vol. 64,
! pp. 490-509, 1986.
!
! (17) I. J. Thompson and A. R. Barnett, "COULCC: A continued-fraction
! algorithm for Coulomb functions of complex order with complex arguments,"
! Computer Physics Communications, vol. 36, pp. 363-372, 1985.
!
! Footnote: The papers appearing in IEICE Transactions can also be inspected
! through the homepage: http://www.ieice.org/eng/index.html

 MODULE mod_bes
! MODULE mod_bes declares the following three groups of variables.
! Group (1) kp, huge1, ihuge, epsilon0, err_range1, err_range2, err_range3,
!           epsilon1, huge1, tiny1, hugelog, sqrt_huge1, theta_lim, alog_base.
!           These are processor dependent-constants. The values of sqrt_huge1,
!           hugelog, epsilon1, err_range1, err_range2, err_range3, alog_base, 
!           theta_lim are defined in FUNCTION num_region.
! Group (2) zunit, pi.
!           zunit=(0,1),  pi=the circular constant.
! Group (3) mmax1, dby1.
!           Bkn in Eqs. (10) and (17) of Ref. (1) is stored in the one-
!           dimensional array dby1 after transformation for saving memory.
!           dby1 is used in SUBROUTINEs sumaabb and  bes_han_dby.
!           The values of DIMENSION dby1 are given in this module.
!           The integer constant mmax1 determines the upper-bound of dimension 
!           dby1.
!
! INTEGER,PARAMETER:: kp=KIND(1.0)                   ! For single precision
 INTEGER,PARAMETER:: kp=KIND(1D0)                   ! For double precision
! INTEGER,PARAMETER:: kp=SELECTED_REAL_KIND(20,550)  ! For quadruple precision
 COMPLEX(kp),PARAMETER:: zunit=(0,1)
 REAL(kp),PARAMETER:: huge1=HUGE(1._kp)/3, tiny1=TINY(1._kp)
 REAL(kp),PARAMETER:: epsilon0=EPSILON(1._kp)
 REAL(kp),PARAMETER:: pi=3.1415926535897932_kp
 INTEGER,PARAMETER:: ihuge=HUGE(1)-1, mmax1=21
 REAL(kp):: err_range1,err_range2,err_range3,theta_lim
 REAL(kp):: epsilon1,sqrt_huge1,alog_base,hugelog=-1 
 REAL(kp):: dby1(((mmax1+2)*(mmax1+1))/2)
 DATA dby1(  1:156)/&
  1.0000000000000000E+00_kp, 1.2500000000000000E-01_kp,-2.0833333333333333E-01_kp, 7.0312500000000000E-02_kp,&
 -4.0104166666666667E-01_kp, 3.3420138888888889E-01_kp, 7.3242187500000000E-02_kp,-8.9121093750000000E-01_kp,&
  1.8464626736111111E+00_kp,-1.0258125964506173E+00_kp, 1.1215209960937500E-01_kp,-2.3640869140625000E+00_kp,&
  8.7891235351562500E+00_kp,-1.1207002616222994E+01_kp, 4.6695844234262474E+00_kp, 2.2710800170898438E-01_kp,&
 -7.3687943594796317E+00_kp, 4.2534998745388455E+01_kp,-9.1818241543240017E+01_kp, 8.4636217674600735E+01_kp,&
 -2.8212072558200245E+01_kp, 5.7250142097473145E-01_kp,-2.6491430486951556E+01_kp, 2.1819051174421159E+02_kp,&
 -6.9957962737613254E+02_kp, 1.0599904525279999E+03_kp,-7.6525246814118164E+02_kp, 2.1257013003921712E+02_kp,&
  1.7277275025844574E+00_kp,-1.0809091978839466E+02_kp, 1.2009029132163525E+03_kp,-5.3056469786134031E+03_kp,&
  1.1655393336864533E+04_kp,-1.3586550006434137E+04_kp, 8.0617221817373094E+03_kp,-1.9194576623184070E+03_kp,&
  6.0740420012734830E+00_kp,-4.9391530477308801E+02_kp, 7.1095143024893637E+03_kp,-4.1192654968897551E+04_kp,&
  1.2220046498301746E+05_kp,-2.0340017728041553E+05_kp, 1.9254700123253153E+05_kp,-9.6980598388637513E+04_kp,&
  2.0204291330966149E+04_kp, 2.4380529699556064E+01_kp,-2.4998304818112096E+03_kp, 4.5218768981362726E+04_kp,&
 -3.3164517248456358E+05_kp, 1.2683652733216248E+06_kp,-2.8135632265865341E+06_kp, 3.7632712976564040E+06_kp,&
 -2.9980159185381068E+06_kp, 1.3117636146629772E+06_kp,-2.4291918790055133E+05_kp, 1.1001714026924674E+02_kp,&
 -1.3886089753717041E+04_kp, 3.0818640461266240E+05_kp,-2.7856181280864547E+06_kp, 1.3288767166421818E+07_kp,&
 -3.7567176660763351E+07_kp, 6.6344512274729027E+07_kp,-7.4105148211532658E+07_kp, 5.0952602492664642E+07_kp,&
 -1.9706819118432227E+07_kp, 3.2844698530720378E+06_kp, 5.5133589612202059E+02_kp,-8.4005433603024085E+04_kp,&
  2.2437681779224494E+06_kp,-2.4474062725738728E+07_kp, 1.4206290779753310E+08_kp,-4.9588978427503031E+08_kp,&
  1.1068428168230145E+09_kp,-1.6210805521083371E+09_kp, 1.5535968995705801E+09_kp,-9.3946235968157840E+08_kp,&
  3.2557307418576575E+08_kp,-4.9329253664509962E+07_kp, 3.0380905109223843E+03_kp,-5.4984232757228869E+05_kp,&
  1.7395107553978165E+07_kp,-2.2510566188941528E+08_kp, 1.5592798648792575E+09_kp,-6.5632937926192843E+09_kp,&
  1.7954213731155600E+10_kp,-3.3026599749800723E+10_kp, 4.1280185579753974E+10_kp,-3.4632043388158778E+10_kp,&
  1.8688207509295825E+10_kp,-5.8664814920518472E+09_kp, 8.1478909611831211E+08_kp, 1.8257755474293175E+04_kp,&
 -3.8718334425726126E+06_kp, 1.4315787671888898E+08_kp,-2.1671649832237951E+09_kp, 1.7634730606834969E+10_kp,&
 -8.7867072178023266E+10_kp, 2.8790064990615059E+11_kp,-6.4536486924537650E+11_kp, 1.0081581068653821E+12_kp,&
 -1.0983751560812233E+12_kp, 8.1921866954857733E+11_kp,-3.9909617522446650E+11_kp, 1.1449823773202581E+11_kp,&
 -1.4679261247695617E+10_kp, 1.1883842625678325E+05_kp,-2.9188388122220813E+07_kp, 1.2470092935127103E+09_kp,&
 -2.1822927757529224E+10_kp, 2.0591450323241002E+11_kp,-1.1965528801961816E+12_kp, 4.6127257808491320E+12_kp,&
 -1.2320491305598287E+13_kp, 2.3348364044581841E+13_kp,-3.1667088584785158E+13_kp, 3.0565125519935321E+13_kp,&
 -2.0516899410934437E+13_kp, 9.1093411852398990E+12_kp,-2.4062979000285040E+12_kp, 2.8646403571767904E+11_kp,&
  8.3285930401628930E+05_kp,-2.3455796352225152E+08_kp, 1.1465754899448237E+10_kp,-2.2961937296824647E+11_kp,&
  2.4850009280340853E+12_kp,-1.6634824724892481E+13_kp, 7.4373122908679145E+13_kp,-2.3260483118893993E+14_kp,&
  5.2305488257844466E+14_kp,-8.5746103298289505E+14_kp, 1.0269551960827625E+15_kp,-8.8949693988102644E+14_kp,&
  5.4273966498765972E+14_kp,-2.2134963870252520E+14_kp, 5.4177510755106049E+13_kp,-6.0197234172340054E+12_kp,&
  6.2529514934347970E+06_kp,-2.0016469281917763E+09_kp, 1.1099740513917901E+11_kp,-2.5215584749128546E+12_kp,&
  3.1007436472896461E+13_kp,-2.3665253045164925E+14_kp, 1.2126758042503474E+15_kp,-4.3793258383640154E+15_kp,&
  1.1486706978449752E+16_kp,-2.2268225133911143E+16_kp, 3.2138275268586241E+16_kp,-3.4447226006485145E+16_kp,&
  2.7054711306197081E+16_kp,-1.5129826322457681E+16_kp, 5.7057821590236708E+15_kp,-1.3010127235496994E+15_kp,&
  1.3552215870309369E+14_kp, 5.0069589531988926E+07_kp,-1.8078220384658064E+10_kp, 1.1287091454108741E+12_kp/
 DATA dby1(157:253)/&
 -2.8863837631414760E+13_kp, 4.0004445704303624E+14_kp,-3.4503855118462725E+15_kp, 2.0064271476309531E+16_kp,&
 -8.2709456515850643E+16_kp, 2.4960365126160426E+17_kp,-5.6263178807463603E+17_kp, 9.5753350981691387E+17_kp,&
 -1.2336116931960695E+18_kp, 1.1961991142756308E+18_kp,-8.5925779803175480E+17_kp, 4.4347954614171904E+17_kp,&
 -1.5552983504313903E+17_kp, 3.3192764720355222E+16_kp,-3.2541926196426688E+15_kp, 4.2593921650476691E+08_kp,&
 -1.7228323871735050E+11_kp, 1.2030115826419192E+13_kp,-3.4396530474307595E+14_kp, 5.3351069787088387E+15_kp,&
 -5.1605093193485227E+16_kp, 3.3766762497906096E+17_kp,-1.5736434765189599E+18_kp, 5.4028948767159819E+18_kp,&
 -1.3970803516443374E+19_kp, 2.7572829816505189E+19_kp,-4.1788614446568389E+19_kp, 4.8599427293248358E+19_kp,&
 -4.3015557038314437E+19_kp, 2.8465212251676571E+19_kp,-1.3639420410571591E+19_kp, 4.4702009640123102E+18_kp,&
 -8.9661142152704633E+17_kp, 8.3019576067319105E+16_kp, 3.8362551802304335E+09_kp,-1.7277040123529995E+12_kp,&
  1.3412416915180639E+14_kp,-4.2619355104268983E+15_kp, 7.3516636109309704E+16_kp,-7.9216511193238321E+17_kp,&
  5.7898876676646531E+18_kp,-3.0255665989903720E+19_kp, 1.1707490535797259E+20_kp,-3.4346213997684169E+20_kp,&
  7.7567049534611368E+20_kp,-1.3602037772849941E+21_kp, 1.8571089321463452E+21_kp,-1.9677247077053125E+21_kp,&
  1.6016898573693597E+21_kp,-9.8244384276898582E+20_kp, 4.3927922008887120E+20_kp,-1.3512175034359961E+20_kp,&
  2.5563802960529235E+19_kp,-2.2424388561867750E+18_kp, 3.6468400807065559E+10_kp,-1.8187262038511037E+13_kp,&
  1.5613123930484673E+15_kp,-5.4840336038832897E+16_kp, 1.0461721131134344E+18_kp,-1.2483700995047233E+19_kp,&
  1.0126774169536592E+20_kp,-5.8917941350694964E+20_kp, 2.5489611146649716E+21_kp,-8.4059158171083504E+21_kp,&
  2.1487414815055883E+22_kp,-4.3025343034823785E+22_kp, 6.7836616429518832E+22_kp,-8.4232227500843226E+22_kp,&
  8.1943310054351296E+22_kp,-6.1732063028844146E+22_kp, 3.5284358439034094E+22_kp,-1.4787743528433614E+22_kp,&
  4.2852960828294940E+21_kp,-7.6719439367290041E+20_kp, 6.3932866139408367E+19_kp, 3.6490108188498336E+11_kp,&
 -2.0052440123627112E+14_kp, 1.8944069842521434E+16_kp,-7.3195014915661331E+17_kp, 1.5365025218443373E+19_kp,&
 -2.0197335419300873E+20_kp, 1.8081594057131944E+21_kp,-1.1640246461465369E+22_kp, 5.5915913803662631E+22_kp,&
 -2.0566149136271543E+23_kp, 5.8965434619782448E+23_kp,-1.3337178907798302E+24_kp, 2.3967237744351683E+24_kp,&
 -3.4308728985157458E+24_kp, 3.9052641035369849E+24_kp,-3.5110965283326441E+24_kp, 2.4615060854038751E+24_kp,&
 -1.3170969618092386E+24_kp, 5.1942890947668122E+23_kp,-1.4228394823321414E+23_kp, 2.4174615008963789E+22_kp,&
 -1.9186202388066499E+21_kp/
 CONTAINS

 INTEGER FUNCTION num_region(znu,xx)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
! This FUNCTION determines nregion.  nregion=num_region(znu,xx)
! This defines the values of sqrt_huge1, hugelog, epsilon1,
! err_range1, err_range2, err_range3, alog_base, theta_lim.
 COMPLEX(kp):: zb,zeps,zgamm,znud
 REAL(kp):: a1,fn,gammr,gammi,re_znu
 REAL(kp):: abszeps,absimznu,absreznu,rp,rr2
 INTEGER:: nn
! Computing  sqrt_huge1, hugelog, epsilon1, err_range1, err_range2,
! err_range3, alog_base, theta_lim.
 IF(hugelog < 0) THEN
   sqrt_huge1=SQRT(huge1)/2
   hugelog=LOG(huge1)-1
   epsilon1=MAX(epsilon0,2.2205E-16_kp)
   err_range1=5E3*epsilon1/epsilon0
   err_range2=2E7*epsilon1/epsilon0
   err_range3=0.1*epsilon1/epsilon0**2
   alog_base=LOG(REAL(RADIX(1._kp)))
   theta_lim=0.08/epsilon1
 ENDIF
! Out of the range.  nregion=0
 IF(xx<0 .OR. xx>theta_lim) THEN;   num_region=0;   RETURN;   ENDIF
 re_znu=REAL(znu,kp)
 absreznu=ABS(re_znu)
 absimznu=ABS(AIMAG(znu))
 IF(absreznu>sqrt_huge1 .OR. absimznu>sqrt_huge1) THEN
   num_region=0;  RETURN;  ENDIF
! By series expansion.  nregion=1,2
 IF(xx < 1) THEN
   IF(absreznu > ihuge) THEN;   num_region=0;   RETURN;   ENDIF
   nn=NINT(re_znu);  zeps=znu-nn
   abszeps=ABS(zeps)
   IF(abszeps > 0.3) THEN;   num_region=1;     RETURN;   ENDIF
   num_region=2;     RETURN
 ENDIF
! By the recurrence formula.  nregion=6,7
 IF(absimznu<20 .AND. ((xx<23 .AND. absreznu<20) &
          .OR. (23<=xx .AND. xx<28 .AND. absreznu<20-4*(xx-23)))) THEN
   IF(absimznu > 0.3) THEN;     num_region=6;     RETURN;   ENDIF
   num_region=7;     RETURN
 ENDIF
! By Debye's method.  nregion=3,4
 rr2=(absreznu-xx)**2+absimznu**2
 IF(xx < 20) THEN
   rp=0.512900114*xx+12
 ELSE
   rp=8.2*xx**0.33333
 ENDIF
 IF(rr2 > rp**2) THEN
   znud=CMPLX(absreznu,absimznu,kp)/xx
   zb=znud+SQRT(znud*znud-1)
   zgamm=LOG(zb)
   gammr=REAL(zgamm,kp);  gammi=AIMAG(zgamm)
   IF(ABS(gammi) <= .001) THEN
     a1=1-gammi**2/3
   ELSE
     a1=gammi/TAN(gammi)
   ENDIF
   fn=gammr*TANH(gammr)+a1-1
   IF(fn <= 0) THEN
     num_region=3;    RETURN      !     |znu/xx|<=1
   ENDIF
   num_region=4;      RETURN      !     |znu/xx|>1
 ENDIF
! Method by Olver's asymptotic series.  nregion=5
 num_region=5
 END FUNCTION num_region

 SUBROUTINE bes_series(znu,xx,zsum,zlogbes,error,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zsum,zlogbes
 REAL(kp),INTENT(OUT):: error
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zsum, zlogbes using series expansion in Section
! 2.2.1 and Eq. (47) of Ref. (1) when nregion=1, 2.
! zbessel(znu,xx)=zsum*EXP(zlogbes)
! This is invoked by bessel, neumann, hankel2.
! error*epsilon0 = a rough estimation of the relative error of zbessel(znu,xx).
! The argument info denotes the output information.
 REAL(kp),DIMENSION(2):: carg,cans
 COMPLEX(kp):: za,zloggam,z1
 REAL(kp):: xx2
 INTEGER:: i,lf0
 zsum=0;  info=0;  zlogbes=0
 IF(ABS(REAL(znu,kp)) > ihuge) THEN;   info=30;    RETURN;   ENDIF
 xx2=(xx/2)**2
 za=1;  zsum=za
 DO i=1,500
   za=-za*xx2/(i*(znu+i))
   zsum=za+zsum
   IF(abs2(za) < epsilon1) EXIT
 ENDDO
 error=epsilon0*abs2(znu+1)
 lf0=0
 carg(1)=REAL(znu+1,kp);     carg(2)=AIMAG(znu)
 CALL cdlgam(carg,cans,error,lf0)
 zloggam=CMPLX(cans(1),cans(2),kp)
 IF(xx < tiny1) THEN
   IF(REAL(znu,kp) > 0) THEN;   zlogbes=-2E20_kp;   RETURN;   ENDIF
   IF(REAL(znu,kp) == 0) THEN
     IF(AIMAG(znu) /= 0) info=20
     zlogbes=0;     RETURN
   ENDIF
   info=20;    zlogbes=2E20_kp;     RETURN    !  REAL(znu) < 0
 ENDIF
 z1=LOG(xx/2)*znu 
 zlogbes=z1-zloggam
 error=abs2(z1)+error/epsilon0
 END SUBROUTINE bes_series

 SUBROUTINE neu_series(znu,xx,zneu,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zneu
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zneu using the series expansion when nregion=2.
! zneu=neumann(znu,xx).  Refer to Ref. (11).
! This is available only for REAL(znu)>=0.
! The argument info denotes the output information.
! This is invoked by bessel, neumann, hankel2.
 COMPLEX(kp):: zeps,z1,z2,z3
 REAL(kp):: a1,re_znu
 INTEGER:: i,info1,nn
 zneu=0;   info=0
 IF(xx < tiny1) THEN;   zneu=huge1;   info=20;   RETURN;  ENDIF
 re_znu=znu;  nn=NINT(re_znu)
 zeps=znu-nn
 IF(nn <= 1) THEN
   CALL neu_srs_init(nn,zeps,xx,zneu,info1)
   info=MAX(info,info1)
   IF(info > 17) THEN;   zneu=huge1;    RETURN;   ENDIF
 ELSE   !  nn >= 2
! When nn >= 2, zneumann(zeps,xx), zneumann(1+zeps,xx) are determined first,
! and zneumann(nn+zeps,xx) is determined by the forward recurrence method.
! nn+zeps=znu.  zeps is epsilon in Eq. (2) of Ref. (11).
   CALL neu_srs_init(0,zeps,xx,z1,info1)
   info=MAX(info,info1)
   CALL neu_srs_init(1,zeps,xx,z2,info1)
   info=MAX(info,info1)
   IF(info > 17) THEN;     zneu=huge1;     RETURN;   ENDIF
   IF(nn > huge1*xx/10) THEN;     info=30;   RETURN;   ENDIF
   a1=abs2(2*((nn-1)+zeps)/xx)*5
   z2=z2/a1;   z1=z1/a1
   DO i=2,nn
     z3=(2*((i-1)+zeps)/xx)*z2-z1
     IF(abs2(z3) > huge1/a1) THEN;   zneu=huge1;   info=20;   RETURN;   ENDIF
     z1=z2;    z2=z3
   ENDDO
   zneu=z3*a1
 ENDIF
 END SUBROUTINE neu_series

 SUBROUTINE neu_srs_init(nn,zeps,xx,zneu,info)
 IMPLICIT NONE
 INTEGER,INTENT(IN):: nn
 COMPLEX(kp),INTENT(IN):: zeps
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zneu
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE is a subroutine subordinate to SUBROUTINE neu_series.
! This SUBROUTINE calculates zneu for nn=0,1 using the series expansion.
! znu=nn+zeps.  nn, zeps are defined by Eq. (8) of Ref. (11).
! zneu=neumann(znu,xx).
! This is invoked by neu_series.
 COMPLEX(kp):: zbes0,zbesp,zbesm,zbes2,zcsbes,zepspi,zepsm,zeps2,zdcos,zsin1, &
   z0,z1
 REAL(kp):: a1,bes0
 INTEGER:: i,nnm
 zneu=0;  info=0
 zepsm=-zeps;  nnm=-nn
 CALL def_bessel(nn,zeps,xx,bes0,zbesp,info)
 IF(info > 17) THEN;   zneu=huge1;     RETURN;   ENDIF
 CALL def_bessel(nnm,zepsm,xx,a1,zbesm,info)
 IF(info > 17) THEN;   zneu=huge1;   RETURN;   ENDIF
 IF(nn >= 1) zbesm=-zbesm
 zbes2=zbesp+zbesm
 zbes0=bes0+zeps*zbesp    ! zbes0=zbessel(nn+zeps,xx)

! Computation of zdcos.
! zdcos is Df{cos(zeps*pi)} defined by Eq. (10) of Ref. (11).
 zepspi=pi*zeps
 zeps2=zepspi**2
 z0=-.5_kp;  z1=z0
 DO i=3,51,2
   z0=-z0*zeps2/(i*(i+1))
   z1=z1+z0
   IF(abs2(z0) < epsilon1) EXIT
   IF(i > 50) THEN;   info=35;   RETURN;   ENDIF
 ENDDO
 zdcos=pi*zepspi*z1

 zcsbes=zdcos*zbes0
! zcsbes=zbessel(nn+zeps,xx)*Df{cos(zeps*pi)}

 IF(abs2(zepspi) < 1E-5) THEN
   zsin1=pi*(1-zepspi**2/6)
 ELSE
   zsin1=SIN(zepspi)/zeps
 ENDIF
! zsin1=sin(pi*zeps)/zeps.  zsin1 is defined by Eq. (12) of Ref. (11).
! The following line depends on Eq. (9) of Ref. (11).
 zneu=(zbes2+zcsbes)/zsin1
 END SUBROUTINE neu_srs_init

 SUBROUTINE def_bessel(nn,zeps,xx,bes0,zbes1,info)
 IMPLICIT NONE
 INTEGER,INTENT(IN):: nn
 COMPLEX(kp),INTENT(IN):: zeps
 REAL(kp),INTENT(IN):: xx
 REAL(kp),INTENT(OUT):: bes0
 COMPLEX(kp),INTENT(OUT):: zbes1
 INTEGER,INTENT(OUT):: info
! This is invoked by neu_srs_init.
! znu=nn+zeps.  nn, zeps are defined by Eq. (8) of Ref. (11).
! This subroutine calculates bes0, zbes1.
! zbessel(nn+zeps,xx)=zbes0+zeps*zbes1,   nn=-1,0,1
! zbessel(nn     ,xx)=zbes0
! Df{zbessel(nn+zeps,xx)}=zbes1
! Df{zbessel(nn+zeps,xx)} is shown in Eq. (9) of Ref. (11).
! The operator Df is defined in Eq. (3) of Ref. (11).
! The argument info denotes the output information.
 COMPLEX(kp):: za,zexpgam1,zexp1,zgamma1,zsig1,z1,z3
 COMPLEX(kp):: zd1,ze1,za1,zs,zss,zeps3
 REAL(kp):: amax,alogxx,a0,a1,a2,expgam0,exp0,gamma0,sig0,bmax,b2,d0,e0,xx2
 INTEGER:: i,i5,n1
 INTEGER, PARAMETER:: m3=22
 REAL(kp):: da(m3)
 DATA da/ &
  0.57721566490153286_kp, -0.65587807152025388_kp,&
 -0.04200263503409524_kp,  0.16653861138229149_kp,&
 -0.04219773455554434_kp, -0.00962197152787697_kp,&
  0.00721894324666310_kp, -0.00116516759185907_kp,&
 -0.00021524167411495_kp,  0.00012805028238812_kp,&
 -0.00002013485478079_kp, -0.00000125049348214_kp,&
  0.00000113302723198_kp, -0.00000020563384170_kp,&
  0.00000000611609510_kp,  0.00000000500200764_kp,&
 -0.00000000118127457_kp,  0.00000000010434267_kp,&
  0.00000000000778226_kp, -0.00000000000369681_kp,&
  0.00000000000051004_kp, -0.00000000000002058_kp/
 info=0;   bes0=0;   zbes1=0

! Calculation of exp0, zexp1
! exp0+zeps*zexp1=(xx/2)**(nn+zeps)
! Refer to Eqs. (20), (21) of Ref. (11).
 z3=(xx/2)**zeps-1
 alogxx=LOG(xx/2)
 IF(abs2(zeps)*alogxx < -50) THEN; info=30; RETURN; ENDIF ! To avoid cancellation
 IF(abs2(z3) < 0.5) THEN
   za=1
   zss=zeps*alogxx
   i5=NINT(AIMAG(zss)/(2*pi))
   zss=CMPLX(REAL(zss,kp),AIMAG(zss)-i5*2*pi,kp)
   zexp1=1
   DO i=2,50
     za=za*zss/i
     zexp1=zexp1+za
     IF(abs2(za) < epsilon1) EXIT
     IF(i == 50) THEN;   info=35;   RETURN;    ENDIF
   ENDDO
   z1=0
   IF(i5 /= 0) z1=CMPLX(0,i5*2*pi,kp)/zeps
   zexp1=zexp1*(alogxx-z1)
 ELSE
   zexp1=z3/zeps
 ENDIF
 IF(alogxx*(nn-abs2(zeps)) > hugelog) THEN;   info=30;   RETURN;  ENDIF
 exp0=EXP(alogxx*nn)
 zexp1=zexp1*exp0

! Calculation of 1/Gamma(1+nn1+zeps)         nn1=0,1
! gamma0+zeps*zgamma1=1/Gamma(1+nn1+zeps)    nn1=MAX(nn,0)   nn=-1,0,1
! Refer to Eqs. (22)-(25) of Ref. (11) and Ref. (14).
 gamma0=1;  zgamma1=0
 zeps3=1;   amax=0
 DO i=1,m3
   zs=da(i)*zeps3
   zgamma1=zgamma1+zs
   a1=abs2(zs)
   amax=MAX(amax,a1)
   IF(a1/amax < epsilon1) EXIT
   IF(i == m3) THEN;   info=35;    RETURN;   ENDIF
   zeps3=zeps3*zeps
 ENDDO
 IF(nn >= 1) THEN       !  nn=nn1=1
   a0=gamma0;   z1=zgamma1
   CALL divis(zeps,a0,z1,1._kp,(1._kp,0._kp),gamma0,zgamma1)
 ENDIF

 CALL multi(zeps,exp0,zexp1,gamma0,zgamma1,expgam0,zexpgam1)
! expgam0+zeps*zexpgam1=(xx/2)**(nn+zeps)/Gamma(1+nn1+zeps)    nn1=MAX(nn,0)

! Calculation of zsigma
! zsigma equals to sigma defined by Eqs. (17), (18a) of Ref. (11).
 xx2=(xx/2)**2
 IF(nn >= 0) THEN   !  nn=0,1
   sig0=1;  zsig1=0
   e0=1;    n1=0
 ELSE               !  nn=-1
   sig0=-xx2;   zsig1=1
   e0=-xx2;     n1=1
 ENDIF
 ze1=0
 amax=MAX(ABS(sig0),tiny1)
 bmax=MAX(abs2(zsig1),tiny1)
 DO i=1,200
   n1=n1+1
   a0=nn+n1
   za1=1
   d0=-e0*xx2/n1
   zd1=-ze1*xx2/n1
   CALL divis(zeps,d0,zd1,a0,za1,e0,ze1)
   sig0=sig0+e0;     zsig1=zsig1+ze1
   a2=ABS(e0);       amax=MAX(a2,amax)
   b2=abs2(ze1);     bmax=MAX(b2,bmax)
   IF(a2/amax<epsilon1 .AND. b2/bmax<epsilon1) EXIT
   IF(i >= 12) THEN;   info=35;     RETURN;   ENDIF
 ENDDO
! zsigma=sig0+zeps*zsig1
 CALL multi(zeps,sig0,zsig1,expgam0,zexpgam1,bes0,zbes1)
 CONTAINS

 SUBROUTINE multi(zeps,ain1,zin1,ain2,zin2,aout,zout)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: zeps,zin1,zin2
 REAL(kp),INTENT(IN):: ain1,ain2
 REAL(kp),INTENT(OUT):: aout
 COMPLEX(kp),INTENT(OUT):: zout
! This SUBROUTINE calculates aout, zout.
! aout+zeps*zout=(ain1+zeps*zin1)*(ain2+zeps*zin2)
! Refer to Eqs. (6c), (7c) of Ref. (11).
 aout=ain1*ain2
 zout=ain2*zin1+ain1*zin2+zeps*zin1*zin2
 END SUBROUTINE multi

 SUBROUTINE divis(zeps,ain1,zin1,ain2,zin2,aout,zout)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: zeps,zin1,zin2
 REAL(kp),INTENT(IN):: ain1,ain2
 REAL(kp),INTENT(OUT):: aout
 COMPLEX(kp),INTENT(OUT):: zout
! This SUBROUTINE calculates aout, zout.
! aout+zeps*zout=(ain1+zeps*zin1)/(ain2+zeps*zin2)
! Refer to Eqs. (6d), (7d) of Ref. (11).
 aout=ain1/ain2
 zout=(zin1-(ain1/ain2)*zin2)/(ain2+zeps*zin2)
 END SUBROUTINE divis
 END SUBROUTINE def_bessel

 SUBROUTINE bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zss1,zss2,zlogmu
 REAL(kp),INTENT(OUT):: error
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zss1, zss2, zlogmu using Debye's expansion when
! nregion=3, 4.
! When nregion=3, 
! zhankel1(znu,xx)=zss1*EXP(zlogmu), 
! zhankel2(znu,xx)=zss2*EXP(-zlogmu).
! When nregion=4, 
! zbessel (znu,xx)=zss1*EXP(zlogmu)/2 for REAL(znu)>=0 and AIMAG(znu)>=0,
! zhankel2(znu,xx)=zss2*EXP(-zlogmu)  for REAL(znu)>=0 and AIMAG(znu)>=0.
! This is invoked by bessel, neumann, hankel2.
! error*epsilon0= a rough estimation of the relative error of the output
! function.
! The argument info denotes the output information.
! This subroutine is based on Refs. (10), (13) and Section 2.2.3 of Ref. (1).
 COMPLEX(kp), DIMENSION(0:mmax1):: zt2n
 COMPLEX(kp):: za1,zbes1,zbes2,zc,zgamm,znu2,zmu,ztnu,ztnuk,ztt
 REAL(kp):: a1,abzb1,abzb2,abzb3
 INTEGER:: i1,k,n
 zss1=0;   zss2=0;   info=0
 IF(xx*epsilon0 > 0.5) THEN;   info=30;   RETURN;  ENDIF
 abzb2=10;   abzb3=1
! zmu is mu in Eqs. (11) of Ref. (1).
 IF(REAL(znu) > xx) THEN
   zmu=SQRT(znu-xx)*SQRT(znu+xx)
 ELSE
   zmu=zunit*SQRT(xx-znu)*SQRT(znu+xx)
 ENDIF
 znu2=znu**2
 ztt=znu2/(znu2-xx**2)
! The double summations of Eqs. (10) of Ref. (1) are computed below.
! The constants Bkn in Eqs. (10) are stored in the one-dimensional array
! dby1(i1).
 zt2n(0)=1;    ztnu=1/zmu
 ztnuk=ztnu;   zbes1=1;   zbes2=1
 a1=-1;        i1=2
 DO k=1,mmax1
   zt2n(k)=zt2n(k-1)*ztt
   zc=0
   DO n=0,k
     zc=zc+dby1(i1)*zt2n(n)
     i1=i1+1
   ENDDO
   za1=zc*ztnuk
   ztnuk=ztnuk*ztnu
   zbes1=zbes1+za1;  zbes2=zbes2+za1*a1
   a1=-a1;  abzb1=abzb2;  abzb2=abzb3
   abzb3=abs2(za1)
   IF(abzb2<=10*epsilon1 .AND. abzb3<=epsilon1) EXIT
   IF(epsilon1<abzb1 .AND. abzb1<abzb2 .AND. abzb2<abzb3) THEN
     info=35;  RETURN;  ENDIF
   IF(k >= mmax1) THEN;   info=35;   RETURN;   ENDIF
 ENDDO
! zgamm is gamma in Eqs. (10) of Ref. (1).
 zgamm=LOG((znu+zmu)/xx)
 zlogmu=zmu-znu*zgamm
 error=abs2(zmu)+abs2(znu*zgamm)
! zss1*EXP( zlogmu) is S(1)(x) of Eq. (10a) of Ref. (1).
! zss2*EXP(-zlogmu) is S(2)(x) of Eq. (10b) of Ref. (1).
 zss1=SQRT(2/(zmu*pi))*zbes1
 zss2=zunit*SQRT(2/(zmu*pi))*zbes2
 END SUBROUTINE bes_han_dby

 SUBROUTINE bes_olver(znu,xx,zbes,error,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zbes
 REAL(kp),INTENT(OUT):: error
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zbes using Olver's expansion when nregion=5.
! zbes=zbessel(znu,xx).
! This is available for REAL(znu) >=0.
! error*epsilon0 = a rough estimation of the relative error of the output zbes.
! This is invoked by bessel, neumann.
! The argument info denotes the output information.
! This subroutine is based on Ref. (8), and Section 2.2.4 of Ref. (1).
 COMPLEX(kp):: za,zarg,zeta,zetad,zfacta,zfactb,znu2,zsum,zz
 REAL(kp):: argr,argi,air1,aii1,air2,aii2
 INTEGER:: ierro
 zbes=0;  error=10
 IF(abs2(znu) > sqrt_huge1) THEN;   info=30;    RETURN;   ENDIF
 znu2=znu**2
 zz=xx/znu
! zz is z in Eq. (15) of Ref. (1).
! zeta is zeta in Eqs. (15) of Ref. (1).
! zetad is zeta/t in Eqs. (14) of Ref. (1).
 CALL fzeta(zz,zeta,zetad,info)
 IF(info > 17) RETURN
 zarg=zeta*znu**(2/3._kp)
 argr=REAL(zarg,kp);     argi=AIMAG(zarg)
 CALL aiz(1,1,argr,argi,air1,aii1,ierro)
 IF(ierro >= 1) THEN;  info=30;   RETURN;   ENDIF
! CMPLX(air1,aii1)=Ai(zeta)
 zfacta=CMPLX(air1,aii1,kp)
 CALL aiz(2,1,argr,argi,air2,aii2,ierro)
 IF(ierro >= 1) THEN;   info=30;   RETURN;  ENDIF
! CMPLX(air2,aii2)=Ai'(zeta)
 za=znu**(1/3._kp)
 zfactb=CMPLX(air2,aii2,kp)/(za*znu)
 CALL sumaabb(znu2,zeta,zetad,zfacta,zfactb,zsum,info)
 IF(info > 17) THEN;   zbes=huge1;    RETURN;   ENDIF
 zbes=zsum*(4*zetad)**.25_kp/za
 IF(ABS(argr)+ABS(argi) > 1) THEN;  info=0;  RETURN;  ENDIF
 error=4*(abs2(znu)/2)**(2./3)*(ABS(air2)+ABS(aii2))/(ABS(air1)+ABS(aii1))
 END SUBROUTINE bes_olver

 SUBROUTINE han2_olver(znu,xx,zhan2,error,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zhan2
 REAL(kp),INTENT(OUT):: error
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zhan2 using Olver's expansion when nregion=5.
! zhan2=zhankel2(znu,xx).
! This is available for REAL(znu) >=0.
! error*epsilon0= a rough estimation of the relative error of the output zhan2.
! This is invoked by bessel, neumann, hankel2.
! The argument info denotes the output information.
! This subroutine is based on Ref. (8), and Section 2.2.4 of Ref. (1).
 COMPLEX(kp):: za,zarg,zeta,zetad,zfacta,zfactb,znu2,zsum,zz
 REAL(kp):: argr,argi,air1,aii1,air2,aii2
 INTEGER:: ierro
 zhan2=0;  error=10
 IF(abs2(znu) > sqrt_huge1) THEN;  info=30;   RETURN;  ENDIF
 znu2=znu**2;  zz=xx/znu
! zz is z in Eq. (15) of Ref. (1).
! zeta is zeta in Eqs. (15) of Ref. (1).
! zetad is zeta/t in Eqs. (14) of Ref. (1).
 CALL fzeta(zz,zeta,zetad,info)
 IF(info > 17) RETURN
 zarg=zeta*znu**(2/3._kp)*exp((-2*pi)*zunit/3)
 argr=REAL(zarg,kp);     argi=AIMAG(zarg)
 CALL aiz(1,1,argr,argi,air1,aii1,ierro)
 IF(ierro >= 1) THEN;  info=30;  RETURN;  ENDIF
! CMPLX(air1,aii1)=Ai(zeta)
 zfacta=EXP(pi*zunit/3)*CMPLX(air1,aii1,kp)
 CALL aiz(2,1,argr,argi,air2,aii2,ierro)
 IF(ierro >= 1) THEN;  info=30;  RETURN;  ENDIF
! CMPLX(air2,aii2)=Ai'(zeta)
 za=znu**(1/3._kp)
 zfactb=EXP(-pi*zunit/3)*CMPLX(air2,aii2,kp)/(za*znu)
 CALL sumaabb(znu2,zeta,zetad,zfacta,zfactb,zsum,info)
 IF(info > 17) RETURN
 zhan2=2*zsum*(4*zetad)**.25_kp/za
 IF(ABS(argr)+ABS(argi) > 1) THEN;  info=0;  RETURN;  ENDIF
 error=4*(abs2(znu)/2)**(2./3)*(ABS(air2)+ABS(aii2))/(ABS(air1)+ABS(aii1))
 END SUBROUTINE han2_olver

 SUBROUTINE sumaabb(znu2,zeta,zetad,zfacta,zfactb,zsum,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu2,zeta,zetad,zfacta,zfactb
 COMPLEX(kp),INTENT(OUT):: zsum
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE is a subordinate subroutine to SUBROUTINEs olver_bes,
! olver_han2.
! This is invoked by bes_olver, han2_olver.
! This subroutine calculates the insides of the braces {...} of Eqs. (14) of 
! Ref. (1).
! The argument info denotes the output information.
 COMPLEX(kp):: zai,zetam,zetamk,zetamkn,zetadsq,zetadk,zetadkn,zeta3, &
   zfact,zn,ztt,z1,zsum0
 REAL(kp):: amsum,amsum0,fa0,fa1,fa2
 INTEGER:: k,k1,mm,m2,n
 INTEGER, PARAMETER:: kt2=31, mmax=6
 REAL(kp), DIMENSION(0:mmax,0:kt2):: tas,tbs
 REAL(kp), DIMENSION(0:mmax1):: duu,dvv
 DATA duu/&
  1.0000000000000000E+00_kp, 1.0416666666666667E-01_kp, 8.3550347222222222E-02_kp, 1.2822657455632716E-01_kp,&
  2.9184902646414046E-01_kp, 8.8162726744375765E-01_kp, 3.3214082818627675E+00_kp, 1.4995762986862555E+01_kp,&
  7.8923013011586518E+01_kp, 4.7445153886826432E+02_kp, 3.2074900908906619E+03_kp, 2.4086549640874005E+04_kp,&
  1.9892311916950979E+05_kp, 1.7919020077753438E+06_kp, 1.7484377180034121E+07_kp, 1.8370737967633073E+08_kp,&
  2.0679040329451552E+09_kp, 2.4827519375935889E+10_kp, 3.1669454981734888E+11_kp, 4.2771126865134717E+12_kp,&
  6.0971132411392562E+13_kp, 9.1486942234356399E+14_kp/
 DATA dvv/&
  1.0000000000000000E+00_kp,-1.4583333333333333E-01_kp,-9.8741319444444444E-02_kp,-1.4331205391589506E-01_kp,&
 -3.1722720267841355E-01_kp,-9.4242914795712025E-01_kp,-3.5112030408263543E+00_kp,-1.5727263620368045E+01_kp,&
 -8.2281439097185944E+01_kp,-4.9235537052367052E+02_kp,-3.3162185685479725E+03_kp,-2.4827674245208590E+04_kp,&
 -2.0452658731512979E+05_kp,-1.8384449170682099E+06_kp,-1.7905687473528919E+07_kp,-1.8783563539939434E+08_kp,&
 -2.1114388546913690E+09_kp,-2.5319153422984124E+10_kp,-3.2261407411300026E+11_kp,-4.3528137960092853E+12_kp,&
 -6.1995857325869748E+13_kp,-9.2950733310106101E+14_kp/
 DATA tas( 1,:)/&
 -4.4444444444444444E-03_kp,-1.4637074635031450E-03_kp, 7.0641727241968957E-04_kp, 6.7288760622093955E-04_kp,&
  1.5400276720923508E-04_kp,-5.7663018476394251E-05_kp,-4.9886522195168320E-05_kp,-1.0429604367829555E-05_kp,&
  3.8752331198978751E-06_kp, 3.1490584761556767E-06_kp, 6.2832879261181444E-07_kp,-2.3288740817602857E-07_kp,&
 -1.8282849503530238E-07_kp,-3.5516623299032355E-08_kp, 1.3117001514284343E-08_kp, 1.0086470236538468E-08_kp,&
  1.9266409049778189E-09_kp,-7.0889344438515159E-10_kp,-5.3776944516159788E-10_kp,-1.0154647136771472E-10_kp,&
  3.7238034463444957E-11_kp, 2.7981059993382281E-11_kp, 5.2398761044890493E-12_kp,-1.9160095921741076E-12_kp,&
 -1.4295821562560585E-12_kp,-2.6603126088049829E-13_kp, 9.7041938442309912E-14_kp, 7.2012262732595635E-14_kp,&
  1.3334884232452868E-14_kp,-4.8543255306589147E-15_kp,-3.5867073033082717E-15_kp,-6.6154066677214473E-16_kp/
 DATA tas( 2,:)/&
  6.9373554135458897E-04_kp, 3.6866079061430036E-04_kp,-2.6986330970626881E-04_kp,-3.5133514343855666E-04_kp,&
 -1.0447400839117945E-04_kp, 5.2408106452547418E-05_kp, 5.5302192195464582E-05_kp, 1.3930130018693300E-05_kp,&
 -6.3002695153511125E-06_kp,-5.9829062080674523E-06_kp,-1.3836191775567955E-06_kp, 5.9343750764257870E-07_kp,&
  5.2990638506672884E-07_kp, 1.1634435490271425E-07_kp,-4.8268492593241100E-08_kp,-4.1408823231282093E-08_kp,&
 -8.7801999003860413E-09_kp, 3.5598512524894753E-09_kp, 2.9693492425525876E-09_kp, 6.1407453218361407E-10_kp,&
 -2.4479643860290886E-10_kp,-1.9999003684120806E-10_kp,-4.0588928059111235E-11_kp, 1.5972063020466106E-11_kp,&
  1.2841509481531001E-11_kp, 2.5683574294034407E-12_kp,-1.0003494294794519E-12_kp,-7.9413286867026562E-13_kp,&
 -1.5697735144899970E-13_kp, 6.0633875508439686E-14_kp, 4.7640389172139987E-14_kp, 9.3270719209874801E-15_kp/
 DATA tas( 3,:)/&
 -3.5421197145774384E-04_kp,-2.4789055466322972E-04_kp, 2.3412119028737724E-04_kp, 3.7696345779888661E-04_kp,&
  1.3525847749464630E-04_kp,-8.2996296644837396E-05_kp,-1.0223189316621075E-04_kp,-2.9785770300621436E-05_kp,&
  1.5692340623662484E-05_kp, 1.6906161946891032E-05_kp, 4.4081304614718038E-06_kp,-2.1365151142413299E-06_kp,&
 -2.1229981061074121E-06_kp,-5.1625698360220605E-07_kp, 2.3720529926670604E-07_kp, 2.2324748774521520E-07_kp,&
  5.1739600581013198E-08_kp,-2.2901873608242246E-08_kp,-2.0728767737587459E-08_kp,-4.6378649944424345E-09_kp,&
  1.9969465687608447E-09_kp, 1.7549815726690861E-09_kp, 3.8222920378669453E-10_kp,-1.6110819379525245E-10_kp,&
 -1.3836296907480886E-10_kp,-2.9500901999602384E-11_kp, 1.2225507750063418E-11_kp, 1.0306976172232469E-11_kp,&
  2.1600570870948614E-12_kp,-8.8286866252751447E-13_kp,-7.3309439954446803E-13_kp,-1.5146399377313705E-13_kp/
 DATA tas( 4,:)/&
  3.7819419920177291E-04_kp, 3.2140419080816258E-04_kp,-3.6482937076827173E-04_kp,-6.9008950585504798E-04_kp,&
 -2.8673227455021060E-04_kp, 2.0513591169330640E-04_kp, 2.8668640832097385E-04_kp, 9.4137483729632422E-05_kp,&
 -5.6144471539415203E-05_kp,-6.7269723879389750E-05_kp,-1.9414648112033552E-05_kp, 1.0438723491294786E-05_kp,&
  1.1371828516705664E-05_kp, 3.0206762376920071E-06_kp,-1.5174520466765040E-06_kp,-1.5487362423995880E-06_kp,&
 -3.8811013921480486E-07_kp, 1.8578505255512123E-07_kp, 1.8078945523008760E-07_kp, 4.3386285194114932E-08_kp,&
 -2.0031362740868734E-08_kp,-1.8795300471185762E-08_kp,-4.3619799376188064E-09_kp, 1.9580390678553789E-09_kp,&
  1.7850769075665493E-09_kp, 4.0335811117396544E-10_kp,-1.7702465954576905E-10_kp,-1.5766439175883287E-10_kp,&
 -3.4857286198319262E-11_kp, 1.5017714449600826E-11_kp, 1.3119509017409341E-11_kp, 2.8483376376472697E-12_kp/
 DATA tas( 5,:)/&
 -6.9114139728829417E-04_kp,-6.8254535963939778E-04_kp, 8.9469107699873365E-04_kp, 1.9267713497801526E-03_kp,&
  9.0294592842939848E-04_kp,-7.3137068939721692E-04_kp,-1.1369660848608659E-03_kp,-4.1320906099652086E-04_kp,&
  2.7356734379850644E-04_kp, 3.5928055401491581E-04_kp, 1.1324042441836278E-04_kp,-6.6614548927214472E-05_kp,&
 -7.8675213133445148E-05_kp,-2.2591357144841771E-05_kp, 1.2281047860783745E-05_kp, 1.3472319402575890E-05_kp,&
  3.6203270420387075E-06_kp,-1.8593422649530313E-06_kp,-1.9312841609483167E-06_kp,-4.9375273357334154E-07_kp,&
  2.4290180978746878E-07_kp, 2.4188218622581626E-07_kp, 5.9479329724974520E-08_kp,-2.8288556570146189E-08_kp,&
 -2.7239049609167991E-08_kp,-6.4918616555387668E-09_kp, 3.0044634200148542E-09_kp, 2.8147196532111456E-09_kp,&
  6.5379907239648796E-10_kp,-2.9584796811080911E-10_kp,-2.7090699820960481E-10_kp,-6.1585693576754200E-11_kp/
 DATA tas( 6,:)/&
  1.9282196424877570E-03_kp, 2.1523979826090860E-03_kp,-3.1755765666574350E-03_kp,-7.6233425583651928E-03_kp,&
 -3.9561737907850920E-03_kp, 3.5566325769629441E-03_kp, 6.0617349749524750E-03_kp, 2.4059483126929379E-03_kp,&
 -1.7431643937939284E-03_kp,-2.4818327268185703E-03_kp,-8.4553509446649723E-04_kp, 5.3841770066316224E-04_kp,&
  6.8341301619312359E-04_kp, 2.1040850952445949E-04_kp,-1.2275813852335326E-04_kp,-1.4372272623730370E-04_kp,&
 -4.1139619424357546E-05_kp, 2.2519837399686868E-05_kp, 2.4821626240847324E-05_kp, 6.7230804658240663E-06_kp,&
 -3.5052321607772612E-06_kp,-3.6862255733981593E-06_kp,-9.5595492348204328E-07_kp, 4.7956442466429617E-07_kp,&
  4.8567567185673600E-07_kp, 1.2159790375358541E-07_kp,-5.9120554485818706E-08_kp,-5.8049445905117900E-08_kp,&
 -1.4117231579053920E-08_kp, 6.6878261599009779E-09_kp, 6.3988169888889280E-09_kp, 1.5185458561661449E-09_kp/
 DATA tbs( 0,:)/&
  1.7998872141355331E-02_kp, 8.8888888888888889E-03_kp, 1.6256871626835735E-03_kp,-3.6428486521990960E-04_kp,&
 -3.0206044899922451E-04_kp,-5.8443572545668709E-05_kp, 1.6769870920170090E-05_kp, 1.3016402516458539E-05_kp,&
  2.4468101612355579E-06_kp,-7.7263598925560737E-07_kp,-5.7902887339204372E-07_kp,-1.0686924823038650E-07_kp,&
  3.5246007722679216E-08_kp, 2.5953663677903903E-08_kp, 4.7402867497067398E-09_kp,-1.5987607555792105E-09_kp,&
 -1.1660522462464008E-09_kp,-2.1163703685703503E-10_kp, 7.2307996005079625E-11_kp, 5.2436882199758380E-11_kp,&
  9.4786533149093573E-12_kp,-3.2652365442990429E-12_kp,-2.3591165030697659E-12_kp,-4.2524948566080466E-13_kp,&
  1.4732015685168165E-13_kp, 1.0616288363349885E-13_kp, 1.9097857813539390E-14_kp,-6.6432344231735149E-15_kp,&
 -4.7782264254825092E-15_kp,-8.5824525371945147E-16_kp, 2.9946699775527868E-16_kp, 2.1508448596067800E-16_kp/
 DATA tbs( 1,:)/&
 -1.4928295321342917E-03_kp,-1.3940630797773655E-03_kp,-3.8209541455316256E-04_kp, 1.6909214802859955E-04_kp,&
  1.7098534913549512E-04_kp, 4.1056073909885070E-05_kp,-1.7066235326534381E-05_kp,-1.5505462076725412E-05_kp,&
 -3.4226070875631647E-06_kp, 1.3772001697435936E-06_kp, 1.1775855270226161E-06_kp, 2.4752762408148760E-07_kp,&
 -9.7522504418527907E-08_kp,-8.0341357113110542E-08_kp,-1.6368639044946662E-08_kp, 6.3525262873386209E-09_kp,&
  5.1050724825074206E-09_kp, 1.0179836640230937E-09_kp,-3.9063793926451710E-10_kp,-3.0842939730776675E-10_kp,&
 -6.0546884180214584E-11_kp, 2.3032399617512805E-11_kp, 1.7946616754605774E-11_kp, 3.4813218340384920E-12_kp,&
 -1.3151890490648631E-12_kp,-1.0143302387096685E-12_kp,-1.9492786982240618E-13_kp, 7.3228717904481976E-14_kp,&
  5.6016705276044817E-14_kp, 1.0683955400531417E-14_kp,-3.9950724984202175E-15_kp,-3.0356851782202066E-15_kp/
 DATA tbs( 2,:)/&
  5.5221307672129279E-04_kp, 7.1104865116708669E-04_kp, 2.5286016094457521E-04_kp,-1.5149350089082804E-04_kp,&
 -1.8614830193107675E-04_kp,-5.3684001061355784E-05_kp, 2.7377121748556901E-05_kp, 2.8968768839784410E-05_kp,&
  7.3912685405114352E-06_kp,-3.4621605971617030E-06_kp,-3.3580620423380641E-06_kp,-7.9598852768413146E-07_kp,&
  3.5393990279009281E-07_kp, 3.2459404774100086E-07_kp, 7.3276227886107934E-08_kp,-3.1453904702404814E-08_kp,&
 -2.7744005279993813E-08_kp,-6.0496501488122393E-09_kp, 2.5314518709133384E-09_kp, 2.1698556404223532E-09_kp,&
  4.6104378160235099E-10_kp,-1.8923284976924440E-10_kp,-1.5868495535973560E-10_kp,-3.3046509649612472E-11_kp,&
  1.3359988711934328E-11_kp, 1.1010832941954977E-11_kp, 2.2565824876368212E-12_kp,-9.0123240099684949E-13_kp,&
 -7.3241059407963607E-13_kp,-1.4815104280203154E-13_kp, 5.8577691212291972E-14_kp, 4.7055618952015272E-14_kp/
 DATA tbs( 3,:)/&
 -4.7461779655995981E-04_kp,-7.5856271658798642E-04_kp,-3.2567548332630983E-04_kp, 2.3883462252518139E-04_kp,&
  3.4254908369517226E-04_kp, 1.1422583074440972E-04_kp,-6.7941577632232691E-05_kp,-8.1521599784337479E-05_kp,&
 -2.3440298279944718E-05_kp, 1.2422076374150894E-05_kp, 1.3401022917758554E-05_kp, 3.5166838941857630E-06_kp,&
 -1.7335316097374033E-06_kp,-1.7433352247630095E-06_kp,-4.2998890758259198E-07_kp, 2.0167782852227561E-07_kp,&
  1.9295507274489267E-07_kp, 4.5504304931699722E-08_kp,-2.0581190241604531E-08_kp,-1.8970882163522409E-08_kp,&
 -4.3243434368499613E-09_kp, 1.9023269668366284E-09_kp, 1.7035057999149978E-09_kp, 3.7808953674304107E-10_kp,&
 -1.6272167712626957E-10_kp,-1.4238696150274311E-10_kp,-3.0930708461281140E-11_kp, 1.3077784859707025E-11_kp,&
  1.1229403531422386E-11_kp, 2.3966017353053719E-12_kp,-9.9854496562904963E-13_kp,-8.4403192768453485E-13_kp/
 DATA tbs( 4,:)/&
  7.3646581057257844E-04_kp, 1.3854690422372401E-03_kp, 6.8917313661988550E-04_kp,-5.8903138992766084E-04_kp,&
 -9.5877430027627460E-04_kp,-3.6032038325291893E-04_kp, 2.4269364604130233E-04_kp, 3.2382750320639959E-04_kp,&
  1.0304984879059063E-04_kp,-6.0606347559499350E-05_kp,-7.1669241373765596E-05_kp,-2.0540810079508347E-05_kp,&
  1.1074795701311619E-05_kp, 1.2075727087270411E-05_kp, 3.2200588553818763E-06_kp,-1.6338956400859472E-06_kp,&
 -1.6804064381103840E-06_kp,-4.2499449686479552E-07_kp, 2.0617990105972662E-07_kp, 2.0287756396323417E-07_kp,&
  4.9271208726441855E-08_kp,-2.3089923952311317E-08_kp,-2.1946145171702701E-08_kp,-5.1615069828158513E-09_kp,&
  2.3531303014594433E-09_kp, 2.1749822913880591E-09_kp, 4.9837356404788743E-10_kp,-2.2216460265718978E-10_kp,&
 -2.0067886737798756E-10_kp,-4.5001100347099178E-11_kp, 1.9689635133257886E-11_kp, 1.7446148243728864E-11_kp/
 DATA tbs( 5,:)/&
 -1.8018219196388570E-03_kp,-3.8637811942002539E-03_kp,-2.1686133637856216E-03_kp, 2.0974523297349090E-03_kp,&
  3.7986118830443898E-03_kp, 1.5801250275318573E-03_kp,-1.1814966487037912E-03_kp,-1.7280200415008586E-03_kp,&
 -6.0051590056809110E-04_kp, 3.8646773366413315E-04_kp, 4.9543877439994537E-04_kp, 1.5348884290142592E-04_kp,&
 -8.9568324917642031E-05_kp,-1.0496507357457789E-04_kp,-3.0011773737039005E-05_kp, 1.6341178768339973E-05_kp,&
  1.7937628794759598E-05_kp, 4.8326649290719236E-06_kp,-2.4985235316640218E-06_kp,-2.6089980116322078E-06_kp,&
 -6.7132220179335078E-07_kp, 3.3337435182826707E-07_kp, 3.3464452136764346E-07_kp, 8.3007566893526022E-08_kp,&
 -3.9911783746651069E-08_kp,-3.8801778371251676E-08_kp,-9.3405917236624775E-09_kp, 4.3738356662889561E-09_kp,&
  4.1409656557539398E-09_kp, 9.7226595671541941E-10_kp,-4.4531582375283522E-10_kp,-4.1229934649535208E-10_kp/
 DATA tbs( 6,:)/&
  6.3858589121204571E-03_kp, 1.5276738425581263E-02_kp, 9.4979998763914257E-03_kp,-1.0191709495545094E-02_kp,&
 -2.0240343119118639E-02_kp,-9.1956609887590049E-03_kp, 7.5243158103600651E-03_kp, 1.1930634443673164E-02_kp,&
  4.4815652491166042E-03_kp,-3.1222288792176478E-03_kp,-4.3015965235668933E-03_kp,-1.4288282009237107E-03_kp,&
  8.9493143585581523E-04_kp, 1.1192667409761611E-03_kp, 3.4087385464895653E-04_kp,-1.9784262236262805E-04_kp,&
 -2.3044235533266587E-04_kp,-6.5772055496720477E-05_kp, 3.6041743154019761E-05_kp, 3.9743877072410871E-05_kp,&
  1.0784589624636757E-05_kp,-5.6494744581294096E-06_kp,-5.9643123900763636E-06_kp,-1.5541044892483991E-06_kp,&
  7.8508098072183484E-07_kp, 7.9990732801321948E-07_kp, 2.0159930364650378E-07_kp,-9.8837517211314685E-08_kp,&
 -9.7770624895054503E-08_kp,-2.3963246861091154E-08_kp, 1.1456852863931195E-08_kp, 1.1053304383088450E-08_kp/
 info = 0
 IF(ABS(zeta) < .86) THEN
! Calculation by series
! Calculation of Am(zeta) by Eq. (18a) of Ref. (1).
! tas(mm,k) means the constants Tamk in Eq. (18a) of Ref. (1).
! mm is m in Eq. (18a) of Ref. (1).
   zfact=zfacta
   fa0=huge1;   fa1=huge1
   ztt=zfact;  zai=ztt
   zsum=zai;   fa2=abs2(zai)
   amsum=abs2(zsum);   zfact=zfact/znu2
   fa0=fa1;    fa1=fa2
   DO mm=1,mmax
     ztt=zfact
     zai=0
     DO k=0,kt2
       z1=tas(mm,k)*ztt
       zai=zai+z1
       IF(abs2(z1) < epsilon1*amsum)  EXIT
       IF(k == kt2) THEN;  info=35;  RETURN;  ENDIF
       ztt=ztt*zeta
     ENDDO
     zsum=zsum+zai
     fa2=abs2(zai)
     amsum=abs2(zsum)
     IF(fa2 < epsilon1*amsum) EXIT
     zfact=zfact/znu2
     IF(fa1>fa0 .AND. fa2>fa1) THEN;  info=35;  RETURN;   ENDIF
     IF(mm == mmax) THEN;  info=35;  RETURN;   ENDIF
     fa0=fa1;   fa1=fa2
   ENDDO
! Calculation of Bm(zeta) by Eq. (18b) of Ref. (1).
! tbs(mm,k) means the constants Tbmk in Eq. (18b) of Ref. (1).
! mm is m in Eq. (18b) of Ref. (1).
   zfact=zfactb;  zsum0=0;  amsum0=0
   fa0=huge1;   fa1=huge1
   DO mm=0,mmax
     ztt=zfact
     zai=0
     DO k=0,kt2
       z1=tbs(mm,k)*ztt
       zai=zai+z1
       IF(abs2(z1) < epsilon1*(amsum+amsum0)) EXIT
       IF(k == kt2) THEN;  info=35;  RETURN;  ENDIF
       ztt=ztt*zeta
     ENDDO
     zsum0=zsum0+zai
     fa2=abs2(zai)
     amsum0=abs2(zsum0)
     IF(fa2 < epsilon1*(amsum+amsum0)) EXIT
     zfact=zfact/znu2
     IF(fa1>fa0 .AND. fa2>fa1) THEN;  info=35;  RETURN;  ENDIF
     IF(mm == mmax) THEN;  info=35;  RETURN;  ENDIF
     fa0=fa1;    fa1=fa2
   ENDDO
   zsum=zsum+zsum0
   RETURN
 ENDIF
! Rigorous calculation by Eqs. (17) of Ref. (1).
! Calculation of Am(zeta) by Eq. (17a) of Ref. (1).
! The constants Bkn in Eqs. (17) are stored in the one-dimensional array 
! dby1(k1).
! The constants v(2m-k) in Eq. (17a) is stored in the array dvv(m2-k).
! mm is m in Eq. (17a).
 info=0;   zsum=0
 zetadsq=SQRT(zetad)
 zeta3=zeta**(-3)
 zfact=zfacta
 fa0=huge1;   fa1=huge1
 zetam=1
 DO mm=0,80
   m2=2*mm
   IF(m2 > mmax1) THEN;   info=35;    RETURN;   ENDIF
   k1=1;  zai=0
   zetamk=zetam       !   zetamk=zeta**(0*k-3*m)
   zetadk=1           !   zetadk=zetad**(0*k/2)
   DO k=0,m2
     zetamkn=zetamk   !  zetamkn=zeta**(k-0*n-3*m)
     zetadkn=zetadk   !  zetadkn=zetad**(0*n+k/2)
     zn=0
     DO n=0,k                              ! zetamkn=zeta**(k-n-3*m)
       zn=dby1(k1)*zetamkn*zetadkn+zn      ! zetadkn=zetad**(n+k/2)
       zetamkn=zetamkn/zeta
       zetadkn=zetadkn*zetad
       k1=k1+1
     ENDDO
     zai=zn*dvv(m2-k)+zai
     zetamk=zetamk*zeta                    !  zetamk=zeta**(k+1-3*m)
     zetadk=zetadk*zetadsq
   ENDDO
   z1=zai*zfact
   zsum=zsum+z1
   fa2=abs2(z1)
   amsum=abs2(zsum)
   IF(fa2 < epsilon1*amsum) EXIT
   zfact=zfact/znu2
   IF(fa1>fa0 .AND. fa2>fa1) THEN;   info=35;   RETURN;   ENDIF
   fa0=fa1;   fa1=fa2
   zetam=zetam*zeta3
 ENDDO
! Calculation of Bm(zeta) by Eq. (17b) of Ref. (1).
! The constants u(2m-k+1) in Eq. (17b) is stored in the array duu(m2-k+1).
! mm is m in Eq. (17b).
 zfact=zfactb
 fa0=huge1;   fa1=huge1
 zetam=zeta**(-2)
 DO mm=0,80
   m2=2*mm+1
   IF(m2 > mmax1) THEN;  info=35;   RETURN;   ENDIF
   k1=1;   zai=0
   zetamk=zetam       !   zetamk=zeta**(0*k-3*m-2)
   zetadk=1           !   zetadk=zetad**(0*k/2)
   DO k=0,m2
     zetamkn=zetamk   !  zetamkn=zeta**(k-0*n-3*m-2)
     zetadkn=zetadk   !  zetadkn=zetad**(0*n+k/2)
     zn=0
     DO n=0,k                              ! zetamkn=zeta**(k-n-3*m-2)
       zn=dby1(k1)*zetamkn*zetadkn+zn      ! zetadkn=zetad**(n+m/2)
       zetamkn=zetamkn/zeta
       zetadkn=zetadkn*zetad
       k1=k1+1
     ENDDO
     zai=zn*duu(m2-k)+zai
     zetamk=zetamk*zeta                    !  zetamk=zeta**(k+1-3*m-2)
     zetadk=zetadk*zetadsq
   ENDDO
   z1=-zai*zfact
   zsum=zsum+z1
   fa2=abs2(z1)
   amsum=abs2(zsum)
   IF(fa2 < epsilon1*amsum) EXIT
   zfact=zfact/znu2
   IF(fa1>fa0 .AND. fa2>fa1) THEN;   info=35;    RETURN;   ENDIF
   fa0=fa1;   fa1=fa2
   zetam=zetam*zeta3
 ENDDO
 END SUBROUTINE sumaabb

 SUBROUTINE fzeta(zz, zeta, zetad, info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: zz
 COMPLEX(kp),INTENT(OUT):: zeta,zetad
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE is a subordinate subroutine to SUBROUTINEs olver_bes, 
! olver_han2.
! This is invoked by olver_bes, olver_han2.
! zeta is defined in Eq. (15) of Ref. (1).
! zz=xx/znu,  ztt=1-zz**2,  zetad=zeta/ztt.
! This subroutine calculates zeta, zetad.
! The argument info denotes the output information.
 INTEGER, PARAMETER:: kdss3=60
 COMPLEX(kp):: za,zeta0,zeta2,ztt,ztta,ztts,zzd,zz2
 INTEGER:: i,imagn
 info=0;  zz2=zz**2;   ztt=1-zz2
 IF(ABS(ztt) < 0.35) THEN
! Calculation of zeta, zetad by the series of Eq. (16) of Ref. (1).
   zeta=0;   ztts=1
   DO i=0,kdss3
     za=ztts/(2*i+3)
     zeta=zeta+za
     IF(abs2(za/zeta) < epsilon1)  EXIT
     ztts=ztts*ztt
     IF(i == kdss3) THEN;   info=35;   RETURN;   ENDIF
   ENDDO
   zetad=(zeta*1.5_kp)**(2/3._kp)
   zeta=ztt*zetad
   RETURN
 ENDIF
! Calculation of zeta, zetad by Eq. (15) of Ref. (1).
 zzd=zz;  imagn=-1
 IF(AIMAG(zz) > 0) THEN;   imagn=1;   zzd=CONJG(zz);  ENDIF
 ztt=1-zzd**2
 ztta=SQRT(ztt)
 za=LOG(1+ztta)-LOG(zzd)-ztta
 zeta0=LOG(za*1.5_kp)
 zeta2=zeta0*(2._kp/3)
 IF(AIMAG(zeta0) < 0) THEN
   zeta2=zeta2-CMPLX(0,pi/6,kp)
   zeta=-zunit*EXP(zeta2)
 ELSE
   zeta=EXP(zeta2)
 ENDIF
 zetad=zeta/ztt
 IF(imagn == 1) THEN
   zeta=CONJG(zeta)
   zetad=CONJG(zetad)
 ENDIF
 END SUBROUTINE fzeta

 SUBROUTINE bes_recur(znu,xx,ifc,zbes,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 INTEGER,INTENT(IN):: ifc
 COMPLEX(kp),INTENT(OUT):: zbes
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zbes using the recurrence method when nregion=6,7.
! zbes=zbessel(znu,xx)*EXP(ifc*zunit*pi*znu), where ifc=-1,0,1.
! This is invoked by bessel, neumann, hankel2.
! The argument info denotes the output information.
! This subroutine is based on Refs. (9), (12) and Section 2.2.5 of Ref. (1)
 COMPLEX(kp):: zalpha,zf(3),zfact,zloggam,zmu,z1
 REAL(kp), DIMENSION(2):: carg,cans
 REAL(kp):: error,re_znu,rnu1,xxd
 INTEGER:: ii,nn,n,n2
 zbes=0;  info=0;  re_znu=znu
 rnu1=MAX(xx-ABS(AIMAG(znu))/1.5089,0._kp)
 nn=MAX(NINT(rnu1-re_znu),0)
 zmu=znu+nn;    xxd=2/xx
 zf(2)=tiny1;   zf(3)=0
 zfact=1;       zalpha=0
! ii equals I defined in Eq. (35) of Ref. (1).
 ii=0.7407*xx+8.359
! Computation of Eqs. (23), (24) of Ref. (1).
! zalpha is alphaI in Eq. (24) of Ref. (1).
 DO n=ii,1,-1
! zf(2)=f(zmu+n2,x)    zf(3)=f(zmu+n2+1,x)    n2=2*n
! f(zmu+n2,x), f(zmu+n2+1,x) are defined by Eqs. (23) of Ref. (1).
   n2=2*n
   zalpha=zalpha+zfact*zf(2)*(zmu+n2)
   zfact=zfact*n/(zmu+(n-1))
   zf(1)=xxd*(zmu+n2    )*zf(2)-zf(3)
   zf(3)=zf(2);   zf(2)=zf(1)
   zf(1)=xxd*(zmu+(n2-1))*zf(2)-zf(3)
   zf(3)=zf(2);   zf(2)=zf(1)
 ENDDO
 zalpha=zalpha+zfact*zf(2)*zmu
 error=epsilon0*abs2(zmu)
 carg(1)=REAL(zmu,kp);     carg(2)=AIMAG(zmu)
 CALL cdlgam(carg,cans,error,0)
 zloggam=CMPLX(cans(1),cans(2),kp)
 z1=ifc*zunit*pi*znu-zloggam-zmu*LOG(xxd)
 IF(ABS(REAL(z1,kp)) > hugelog) THEN;   info=30;   RETURN;   ENDIF
! Computation of Eqs. (25) of Ref. (1).
 zf(2)=zf(2)*zfact/zalpha
 zf(3)=zf(3)*zfact/zalpha
! Computation of Eqs. (26) of Ref. (1).
! zf(2)=zbessel(zmu,x)     zf(3)=zbessel(zmu+1,x)     zmu=znu+nn
 DO n=0,-nn+1,-1
   zf(1)=xxd*(zmu+n)*zf(2)-zf(3)
   zf(3)=zf(2);   zf(2)=zf(1)
 ENDDO
! zf(2)=zbessel(znu,x)     zf(3)=zbessel(znu+1,x)
 IF(EXPONENT(abs2(zf(2)))*alog_base+REAL(z1) > hugelog) THEN
   info=20;   zbes=huge1;   RETURN
 ENDIF
 zbes=zf(2)*EXP(z1)
 END SUBROUTINE bes_recur

 SUBROUTINE han2_temme(znu,xx,zhan2,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(in):: znu
 REAL(kp),INTENT(in):: xx
 COMPLEX(kp),INTENT(out):: zhan2
 INTEGER,INTENT(out):: info
! This SUBROUTINE calculates zhan2 by the use of Temme's algorithm when
! nregion=7.ü@ zhan2=hankel2(znu,xx).
! The argument info denotes the output information.
! This is available for REAL(znu) >=0.
! This is invoked by bessel, neumann, hankel2.
! This subroutine is based on Refs. (4), (7) and Section 2.2.6 of Ref. (1)
 COMPLEX(kp):: zc,zf,zg,zh,zeps,zp,zq,zxx,zxxi
 REAL(kp):: aloge,a1,a2,renu
 INTEGER:: ll,n,nn
 INTEGER, PARAMETER:: ik1=4,ik2=3
 REAL(kp):: dc1(ik1),dc2(ik2)
 DATA dc1/ 7.245428E-1_kp, 8.664168E-1_kp, 1.991080E-1_kp, 1.607260E-3_kp/
 DATA dc2/ 2.260379E+0_kp,-8.945939E-2_kp, 5.189255E-3_kp/
 info=0;  zhan2=0
 renu=REAL(znu,kp)
 nn=NINT(renu);  zeps=znu-nn
 zxx=xx*zunit
 zc=0.25_kp-zeps**2
! Computation of ll.   ll is L defined by Eq. (36a) of Ref. (1).
 aloge=LOG(abs2(COS(zeps*pi))/(pi*epsilon1))
 a1=1;  a2=0
 DO n=1,ik1
   a2=a2+dc1(n)*a1
   a1=a1*aloge
 ENDDO
 a1=1;  a2=a2/xx
 DO n=1,ik2
   a2=a2+dc2(n)*a1
   a1=a1*xx
 ENDDO
 ll=a2
! The backward recurrences of k(n)
 zp=0
 zq=zp         !  zp=zq=k(ll+1)/k(ll)       k(ll+1)=0
 DO n=ll,1,-1
   zp=((n-1)+zc/n)/(2*(zxx+n)-(n+1)*zp)
   zq=zp*(zq+1)
! zp=k(n)/k(n-1)    zq=[k(n)+k(n+1)+....+k(ll)+k(ll+1)]/k(n-1)
 ENDDO
! zp=k(1)/k(0)    1/(zq+1)=k(0)/[k(0)+k(1)+....+k(ll+1)]
 zf=SQRT(pi/(2*zxx))*EXP(-zxx)/(zq+1)
 zg=zf*(zeps+zxx+.5_kp-zp)/zxx
 zxxi=2/zxx
! The forward recurrences of K(zeps+n,zxx)
 DO n=1,nn
! zh, zg, zf= Modified Bessel functions of the third kind.
! zh=K(zeps+n+1,zxx)   zg=K(zeps+n,zxx)   zf=K(zeps+n-1,zxx)
   zh=(zeps+n)*zxxi*zg+zf
   zf=zg;   zg=zh
 ENDDO
! zf=K(zeps+nn,zxx)=K(znu,zxx)=the modified Bessel function of the third kind.
 zhan2=(2*zunit/pi)*EXP(zunit*pi*znu/2)*zf
 END SUBROUTINE han2_temme

! Algorithm 819
! Algorithm 819 [Ref. (5)] was downloaded from the homepage of ACM.
! Algorithm 819 is written in FORTRAN 77, but the present algorithm is written
! in Fortran 95. For convenience, Algorithm 819 used here was translated into
! Fortran 95 by the use of the transformer program stored in file
! transformer.f95.
! In addition, some statements are modified in order that the kind type 
! parameter is variable.
! The modified statements are marked with !!.
!
! Algorithms 644, 838 [Refs. (7), (15)] have programs for the complex Airy 
! functions, which can also be applied to the present algorithm with a few 
! modifications.
!                                                        M. Kodama
!
!**************************************************************
!      ALGORITHM 819, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 28,NO. 3, September, 2002, P.  325--336.

      REAL(kp)          FUNCTION D1MACH(I)
      INTEGER I
!***PURPOSE  RETURNS DOUBLE PRECISION MACHINE DEPENDENT CONSTANTS
!***DESCRIPTION
!
!     D1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
!     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
!     SUBPROGRAM WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
!     AS FOLLOWS, FOR EXAMPLE
!
!          D = D1MACH(I)
!
!     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF D ABOVE IS
!     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
!     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
!
!  DOUBLE-PRECISION MACHINE CONSTANTS
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  D1MACH( 5) = LOG10(B)
!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  D1MACH
!
!***FIRST EXECUTABLE STATEMENT  D1MACH
      REAL(KIND=kp), PARAMETER:: ZERO = 0.0_kp, BASE = 2.0_kp       !!
      SELECT CASE (I)
      CASE(1)
        D1MACH = TINY(ZERO)
      CASE(2)
        D1MACH = HUGE(ZERO)
      CASE(3)
        D1MACH = EPSILON(ZERO)
      CASE(4)
        D1MACH = BASE * EPSILON(ZERO)
      CASE(5)
        D1MACH = LOG10(BASE)
      END SELECT
!
      END FUNCTION
      REAL FUNCTION R1MACH(I)
      INTEGER I
!***PURPOSE  RETURNS SINGLE PRECISION MACHINE DEPENDENT CONSTANTS
!***DESCRIPTION
!
!     R1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
!     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
!     SUBROUTINE WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
!     AS FOLLOWS, FOR EXAMPLE
!
!          A = R1MACH(I)
!
!     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF A ABOVE IS
!     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
!     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
!
!  SINGLE-PRECISION MACHINE CONSTANTS
!  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  R1MACH(5) = LOG10(B)
!***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
!                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
!                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
!                 PP. 177-188.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  R1MACH
!***FIRST EXECUTABLE STATEMENT  R1MACH
      REAL, PARAMETER:: ZERO = 0.0e0, BASE = 2.0e0
      SELECT CASE (I)
      CASE(1)
        R1MACH = TINY(ZERO)
      CASE(2)
        R1MACH = HUGE(ZERO)
      CASE(3)
        R1MACH = EPSILON(ZERO)
      CASE(4)
        R1MACH = BASE * EPSILON(ZERO)
      CASE(5)
        R1MACH = LOG10(BASE)
      END SELECT
!
      END FUNCTION
      INTEGER FUNCTION I1MACH(I)
!***PURPOSE  RETURN INTEGER MACHINE DEPENDENT CONSTANTS.
!***DESCRIPTION
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   THESE MACHINE CONSTANT ROUTINES MUST BE ACTIVATED FOR
!   A PARTICULAR ENVIRONMENT.
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     I1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
!     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
!     SUBROUTINE WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
!     AS FOLLOWS, FOR EXAMPLE
!
!          K = I1MACH(I)
!
!     WHERE I=1,...,16.  THE (OUTPUT) VALUE OF K ABOVE IS
!     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
!     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
!
!  I/O UNIT NUMBERS.
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!
!  WORDS.
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
!
!  INTEGERS.
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
!
!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
!    I1MACH( 7) = A, THE BASE.
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS.
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM
!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, THE BASE.
!
!  SINGLE-PRECISION
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
       INTEGER I
       REAL R
       REAL(kp)          D

       SELECT CASE (I)
       CASE(1)
         I1MACH = 5  ! Standard input
       CASE(2)
         I1MACH = 6  ! Standard output
       CASE(3)
         I1MACH = 6  ! Standard punch :-)
       CASE(4)
         I1MACH = 6  ! Standard error
       CASE(5)
         I1MACH = DIGITS(I) + 1 !  Number of bits /integer (+1 for the s
       CASE(6)
         I1MACH = 4  ! Number of characters / integer :-)
       CASE(7)
         I1MACH = RADIX(I) ! base of integers
       CASE(8)
         I1MACH = DIGITS(I) ! number of base radix digits in integer
       CASE(9)
         I1MACH = HUGE(I) ! Maximum integer
       CASE(10)
         I1MACH = RADIX(R) ! base of floating point
       CASE(11)
         I1MACH = DIGITS(R) ! number of base radix digits in sp
       CASE(12)
         I1MACH = MINEXPONENT(R) ! minimun sp exponent
       CASE(13)
         I1MACH = MAXEXPONENT(R) ! maximum sp exponent
       CASE(14)
         I1MACH = DIGITS(D) ! number of base radix digits in dp
       CASE(15)
         I1MACH = MINEXPONENT(D) ! minimun dp exponent
       CASE(16)
         I1MACH = MAXEXPONENT(D) ! maximum dp exponent
       END SELECT
       END FUNCTION

      SUBROUTINE AIZ(IFUN,IFAC,X0,Y0,GAIR,GAII,IERRO)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COMPUTATION OF THE AIRY FUNCTION AI(Z) OR ITS DERIVATIVE AI'(Z)
! THE CODE USES:
!      1. MACLAURIN SERIES FOR |Y|<3 AND -2.6<X<1.3 (Z=X+I*Y)
!      2. GAUSS-LAGUERRE QUADRATURE  FOR |Z|<15 AND  WHEN
!         MACLAURIN SERIES ARE NOT USED.
!      3. ASYMPTOTIC EXPANSION FOR |Z|>15.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  INPUTS:
!    IFUN:
!         * IFUN=1, THE CODE COMPUTES AI(Z)
!         * IFUN=2, THE CODE COMPUTES AI'(Z)
!    IFAC:
!         * IFAC=1, THE CODE COMPUTES  AI(Z) OR AI'(Z)
!         * IFAC=2, THE CODE COMPUTES NORMALIZED AI(Z) OR AI'(Z)
!    X0:   REAL PART OF THE ARGUMENT Z
!    Y0:   IMAGINARY PART OF THE ARGUMENT  Z
!
!  OUTPUTS:
!    GAIR: REAL PART OF AI(Z) OR AI'(Z)
!    GAII: IMAGINARY PART OF AI(Z) OR AI'(Z)
!
!    IERRO: ERROR FLAG
!          * IERRO=0, SUCCESSFUL COMPUTATION
!          * IERRO=1, COMPUTATION OUT OF RANGE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   ACCURACY:
!
!     1) SCALED AIRY FUNCTIONS:
!        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO
!        THE ZEROS, WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
!        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z|>1000
!        (REACHING 10**(-8) ABSOLUTE ACCURACY FOR |Z| CLOSE
!        TO 10**(6)) IN THE CASE OF PHASE(Z) CLOSE TO PI OR -PI.
!     2) UNSCALED AIRY FUNCTIONS:
!        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR
!        3/2*|Z|**(3/2)>LOG(OVER).
!        FOR |Z|<30:
!        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE
!           ZEROS) BETTER THAN 10**(-13).
!        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER
!           THAN 10**(-13), WHERE R(Z)=REAL(AI)/IMAG(AI)
!           OR R(Z)=REAL(AI')/IMAG(AI').
!        FOR |Z|>30, GRADUAL LOSS OF PRECISION TAKES PLACE
!        AS |Z| INCREASES.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     AUTHORS:
!        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN).
!                      E-MAIL: AMPARO.GIL@UAM.ES
!        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
!                      E-MAIL: JSEGURA@MATH.UC3M.ES
!        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
!                      E-MAIL: NICO.TEMME@CWI.NL
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    REFERENCES:
!         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
!         NUMERICAL ALGORITHMS (2002).
!         A. GIL, J. SEGURA, N.M. TEMME
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL(kp)          X0,Y0,GAIR,GAII,X,W,XD,WD
      REAL(kp)          OVER,UNDER,DL1,DL2,COVER   !,D1MACH
      REAL(kp)          PI,PIHAL,PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,&
      &R32,FACTO,TH025,S3,C3,F23,PI23,SQRT3,XA,YA,F23R,DF1,DF2,&
      &S11,C11,DEX,DRE,DIMA,GAR,GAI,C,S,U,V,V0,AR,AI,AR1,AI1,&
      &RO,COE1,COE2,REX,DFR,DFI,AR11,AI11   !,PHASE
      INTEGER IFUN,IFAC,IERRO,IEXPF,IEXPF2,N
      DIMENSION X(25),W(25)
      DIMENSION XD(25),WD(25)
      COMMON/PARAM1/PI,PIHAL
      COMMON/PARAM2/PIH3,PISR,A,ALF
      COMMON/PARAM3/THET,R,TH15,S1,C1,R32
      COMMON/PARAM4/FACTO,TH025,S3,C3
      SAVE X,W
      SAVE XD,WD
      DATA X,W/.283891417994567679D-1,.170985378860034935D0,&
     &.435871678341770460D0,.823518257913030858D0,1.33452543254227372D0,&
     &1.96968293206435071D0,2.72998134002859938D0,3.61662161916100897D0,&
     &4.63102611052654146D0,5.77485171830547694D0,7.05000568630218682D0,&
     &8.45866437513237792D0,10.0032955242749393D0,11.6866845947722423D0,&
     &13.5119659344693551D0,15.4826596959377140D0,17.6027156808069112D0,&
     &19.8765656022785451D0,22.3091856773962780D0,24.9061720212974207D0,&
     &27.6738320739497190D0,30.6192963295084111D0,33.7506560850239946D0,&
     &37.0771349708391198D0,40.6093049694341322D0,.143720408803313866D0,&
     &.230407559241880881D0,.242253045521327626D0,.203636639103440807D0,&
     &.143760630622921410D0,.869128834706078120D-1,.4541750018329&
      &15883D-1,.206118031206069497D-1,.814278821268606972D-2,.280266&
     &075663377634D-2,.840337441621719716D-3,.219303732907765020D-3,&
     &.497401659009257760D-4,.978508095920717661D-5,.166542824603725&
     &563D-5,.244502736801316287D-6,.308537034236207072D-7,.3332960&
     &72940112245D-8,.306781892316295828D-9,.239331309885375719D-10,&
     &.157294707710054952D-11,.864936011664392267D-13,.394819815&
     &638647111D-14,.148271173082850884D-15,.453390377327054458D-17/
      DATA XD,WD/.435079659953445D-1,.205779160144678D0,&
     &.489916161318751D0,.896390483211727D0,1.42582496737580D0,&
     &2.07903190767599D0,2.85702335104978D0,3.76102058198275D0,&
     &4.79246521225895D0,5.95303247470003D0,7.24464710774066D0,&
     &8.66950223642504D0,10.2300817341775D0,11.9291866622602D0,&
     &13.7699665302828D0,15.7559563095946D0,17.8911203751898D0,&
     &20.1799048700978D0,22.6273004064466D0,25.2389175786164D0,&
     &28.0210785229929D0,30.9809287996116D0,34.1265753192057D0,&
     &37.4672580871163D0,41.0135664833476D0,.576354557898966D-1,&
     &.139560003272262D0,.187792315011311D0,.187446935256946D0,&
     &.150716717316301D0,.101069904453380D0,.575274105486025D-1,&
     &.280625783448681D-1,.117972164134041D-1,.428701743297432D-2,&
     &.134857915232883D-2,.367337337105948D-3,.865882267841931D-4,&
     &.176391622890609D-4,.309929190938078D-5,.468479653648208D-6,&
     &.607273267228907D-7,.672514812555074D-8,.633469931761606D-9,&
     &.504938861248542D-10,.338602527895834D-11,.189738532450555D-12,&
     &.881618802142698D-14,.336676636121976D-15,.104594827170761D-16/
!C CONSTANTS CCCCCCCCCCCCCCCCCCCCCCC
      PI=3.1415926535897932385D0
      PIHAL=1.5707963267948966192D0
      PIH3=4.71238898038469D0
      F23=.6666666666666666D0
      PI23=2.09439510239320D0
      PISR=1.77245385090552D0
      SQRT3=1.7320508075688772935D0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      YA=Y0
      XA=X0
      IERRO=0
      IEXPF=0
      IEXPF2=0
      IF (YA.LT.0.D0) YA=-YA
      R=SQRT(XA*XA+YA*YA)
      R32=R*SQRT(R)
      THET=PHASE(XA,YA)
      COVER=2.D0/3.D0*R32*ABS(COS(1.5D0*THET))
!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
      OVER=D1MACH(2)*1.D-4
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
      UNDER=D1MACH(1)*1.D+4
      DL1=LOG(OVER)
      DL2=-LOG(UNDER)
      IF (DL1.GT.DL2) OVER=1/UNDER
      IF (IFAC.EQ.1) THEN
        IF (COVER.GE.LOG(OVER)) THEN
!CC OVERFLOW/UNDERFLOW PROBLEMS.
!CC   CALCULATION ABORTED
          IERRO=1
          GAIR=0
          GAII=0
        ENDIF
        IF (COVER.GE.(LOG(OVER)*0.2)) IEXPF2=1
      ELSE
        IF (COVER.GE.(LOG(OVER)*0.2)) IEXPF=1
      ENDIF
      IF (IERRO.EQ.0) THEN
        IF (IFUN.EQ.1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF AI(Z) CCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCC SERIES, INTEGRALS OR EXPANSIONS CCCCCCCCCCCCCCCCCCCCCCCC
          IF ((YA.LT.3.D0).AND.(XA.LT.1.3D0).AND.(XA.GT.-2.6D0)) THEN
!CC SERIES CCC
            CALL SERAI(XA,YA,GAR,GAI)
            IF (IFAC.EQ.2) THEN
              THET=PHASE(XA,YA)
              TH15=1.5D0*THET
              S1=SIN(TH15)
              C1=COS(TH15)
              F23R=F23*R32
              DF1=F23R*C1
              DF2=F23R*S1
              S11=SIN(DF2)
              C11=COS(DF2)
              DEX=EXP(DF1)
              DRE=DEX*C11
              DIMA=DEX*S11
              GAIR=DRE*GAR-DIMA*GAI
              GAII=DRE*GAI+DIMA*GAR
            ELSE
              GAIR=GAR
              GAII=GAI
              IF (Y0.EQ.0.) GAII=0.D0
            ENDIF
          ELSE
            IF (R.GT.15.D0) THEN
!CC ASYMPTOTIC EXPANSIONS CCC
              THET=PHASE(XA,YA)
              FACTO=0.5D0/PISR*R**(-0.25D0)
              IF (THET.GT.PI23) THEN
!CCCCCCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
                N=1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL EXPAI(AR1,AI1)
                IF (V0.LT.0.D0) AI1=-AI1
                AR=-(C*AR1-S*AI1)
                AI=-(S*AR1+C*AI1)
                IF (IEXPF.EQ.0) THEN
                  IF (IEXPF2.EQ.0) THEN
!CC     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
                    N=-1
                    C=-0.5D0
                    S=N*0.5*SQRT3
                    U=XA*C-YA*S
                    V=XA*S+YA*C
                    V0=V
                    IF (V.LT.0.D0) V=-V
                    THET=PHASE(U,V)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    TH025=THET*0.25D0
                    S3=SIN(TH025)
                    C3=COS(TH025)
                    CALL EXPAI(AR1,AI1)
                    IF (V0.LT.0.D0) AI1=-AI1
                    THET=PHASE(XA,YA)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    RO=1.333333333333333D0*R32
                    COE1=RO*C1
                    COE2=RO*S1
                    REX=EXP(COE1)
                    DFR=REX*COS(COE2)
                    DFI=REX*SIN(COE2)
                    AR11=DFR*AR1-DFI*AI1
                    AI11=DFR*AI1+DFI*AR1
                    GAIR=AR-(C*AR11-S*AI11)
                    GAII=AI-(S*AR11+C*AI11)
                  ELSE
                    THET=PHASE(XA,YA)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    GAIR=AR
                    GAII=AI
                  ENDIF
                ELSE
                  GAIR=AR
                  GAII=AI
                ENDIF
              ELSE
!CCCCCC  ASYMPTOTIC EXPANSION CCCCCCCCCCCCCCC    	
                THET=PHASE(XA,YA)
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL EXPAI(GAIR,GAII)
              ENDIF
            ELSE
!CC INTEGRALS
              A=0.1666666666666666D0
              ALF=-A
              FACTO=0.280514117723058D0*R**(-0.25D0)
              THET=PHASE(XA,YA)
              IF (THET.LE.PIHAL) THEN
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL AIRY1(X,W,GAIR,GAII)
              ENDIF
              IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL AIRY2(X,W,GAIR,GAII)
              ENDIF
              IF (THET.GT.PI23) THEN
                N=1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                IF (THET.LE.PIHAL) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY1(X,W,AR1,AI1)
                ENDIF
                IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY2(X,W,AR1,AI1)
                ENDIF
                IF (V0.LT.0.D0) AI1=-AI1
                AR=-(C*AR1-S*AI1)
                AI=-(S*AR1+C*AI1)
                N=-1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                IF (THET.LE.PIHAL) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY1(X,W,AR1,AI1)
                ENDIF
                IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY2(X,W,AR1,AI1)
                ENDIF
                IF (V0.LT.0.D0) AI1=-AI1
                THET=PHASE(XA,YA)
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                RO=1.333333333333333D0*R32
                COE1=RO*C1
                COE2=RO*S1
                REX=EXP(COE1)
                DFR=REX*COS(COE2)
                DFI=REX*SIN(COE2)
                AR11=DFR*AR1-DFI*AI1
                AI11=DFR*AI1+DFI*AR1
                GAIR=AR-(C*AR11-S*AI11)
                GAII=AI-(S*AR11+C*AI11)
              ENDIF
            ENDIF
            IF (IFAC.EQ.1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF THE UNSCALED AI(Z) CCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              F23R=F23*R32
              DF1=F23R*C1
              DF2=F23R*S1
              S11=SIN(DF2)
              C11=COS(DF2)
              DEX=EXP(-DF1)
              DRE=DEX*C11
              DIMA=-DEX*S11
              GAR=DRE*GAIR-DIMA*GAII
              GAI=DRE*GAII+DIMA*GAIR
              GAIR=GAR
              GAII=GAI
              IF (Y0.EQ.0.) GAII=0.D0
            ENDIF
          ENDIF
        ELSE
!CCC CALCULATION OF AI³(Z) CCCCCCCCCCC
          ALF=0.1666666666666666D0
          FACTO=-0.270898621247918D0*R**0.25D0
!CCCCCCCCCCCCCC SERIES OR INTEGRALS CCCCCCCCCCCCCCCCCCCCCCCCCC
          IF ((YA.LT.3.D0).AND.(XA.LT.1.3D0).AND.(XA.GT.-2.6D0)) THEN
!CC SERIES
            CALL SERAID(XA,YA,GAR,GAI)
            IF (IFAC.EQ.2) THEN
              THET=PHASE(XA,YA)
              TH15=1.5D0*THET
              S1=SIN(TH15)
              C1=COS(TH15)
              F23R=F23*R32
              DF1=F23R*C1
              DF2=F23R*S1
              S11=SIN(DF2)
              C11=COS(DF2)
              DEX=EXP(DF1)
              DRE=DEX*C11
              DIMA=DEX*S11
              GAIR=DRE*GAR-DIMA*GAI
              GAII=DRE*GAI+DIMA*GAR
            ELSE
              GAIR=GAR
              GAII=GAI
              IF (Y0.EQ.0.) GAII=0.D0
            ENDIF
          ELSE
            IF (R.GT.15.D0) THEN
!CC  ASYMPTOTIC EXPANSIONS CCCCCCCCCCCCC
              THET=PHASE(XA,YA)
              FACTO=0.5D0/PISR*R**0.25D0
              IF (THET.GT.PI23) THEN
!CCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
                N=1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL EXPAID(AR1,AI1)
                IF (V0.LT.0.D0) AI1=-AI1
                AR=-(C*AR1+S*AI1)
                AI=-(-S*AR1+C*AI1)
                IF (IEXPF.EQ.0) THEN
                  IF (IEXPF2.EQ.0) THEN
!CC     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
                    N=-1
                    C=-0.5D0
                    S=N*0.5*SQRT3
                    U=XA*C-YA*S
                    V=XA*S+YA*C
                    V0=V
                    IF (V.LT.0.D0) V=-V
                    THET=PHASE(U,V)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    TH025=THET*0.25D0
                    S3=SIN(TH025)
                    C3=COS(TH025)
                    CALL EXPAID(AR1,AI1)
                    IF (V0.LT.0.D0) AI1=-AI1
                    THET=PHASE(XA,YA)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    RO=1.333333333333333D0*R32
                    COE1=RO*C1
                    COE2=RO*S1
                    REX=EXP(COE1)
                    DFR=REX*COS(COE2)
                    DFI=REX*SIN(COE2)
                    AR11=DFR*AR1-DFI*AI1
                    AI11=DFR*AI1+DFI*AR1
                    GAIR=AR-(C*AR11+S*AI11)
                    GAII=AI-(-S*AR11+C*AI11)
                  ELSE
                    THET=PHASE(XA,YA)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    GAIR=AR
                    GAII=AI
                  ENDIF
                ELSE
                  GAIR=AR
                  GAII=AI
                ENDIF
              ELSE
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL EXPAID(GAIR,GAII)
              ENDIF
            ELSE
!CC INTEGRALS CCCCCCCCCCCCCCCC
              THET=PHASE(XA,YA)
              IF (THET.LE.PIHAL) THEN
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL AIRY1D(XD,WD,GAIR,GAII)
              ENDIF
              IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL AIRY2D(XD,WD,GAIR,GAII)
              ENDIF
              IF (THET.GT.PI23) THEN
                N=1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                IF (THET.LE.PIHAL) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY1D(XD,WD,AR1,AI1)
                ENDIF
                IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY2D(XD,WD,AR1,AI1)
                ENDIF
                IF (V0.LT.0.D0) AI1=-AI1
                  AR=-(C*AR1+S*AI1)
                  AI=-(-S*AR1+C*AI1)
                  N=-1
                  C=-0.5D0
                  S=N*0.5*SQRT3
                  U=XA*C-YA*S
                  V=XA*S+YA*C
                  V0=V
                  IF (V.LT.0.D0) V=-V
                  THET=PHASE(U,V)
                  IF (THET.LE.PIHAL) THEN
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    TH025=THET*0.25D0
                    S3=SIN(TH025)
                    C3=COS(TH025)
                    CALL AIRY1D(XD,WD,AR1,AI1)
                  ENDIF
                  IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    TH025=THET*0.25D0
                    S3=SIN(TH025)
                    C3=COS(TH025)
                    CALL AIRY2D(XD,WD,AR1,AI1)
                  ENDIF
                  IF (V0.LT.0.D0) AI1=-AI1
                  THET=PHASE(XA,YA)
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  RO=1.333333333333333D0*R32
                  COE1=RO*C1
                  COE2=RO*S1
                  REX=EXP(COE1)
                  DFR=REX*COS(COE2)
                  DFI=REX*SIN(COE2)
                  AR11=DFR*AR1-DFI*AI1
                  AI11=DFR*AI1+DFI*AR1
                  GAIR=AR-(C*AR11+S*AI11)
                  GAII=AI-(-S*AR11+C*AI11)
                ENDIF
              ENDIF
              IF (IFAC.EQ.1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF THE UNSCALED AI'(z) CCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                F23R=F23*R32
                DF1=F23R*C1
                DF2=F23R*S1
                S11=SIN(DF2)
                C11=COS(DF2)
                DEX=EXP(-DF1)
                DRE=DEX*C11
                DIMA=-DEX*S11
                GAR=DRE*GAIR-DIMA*GAII
                GAI=DRE*GAII+DIMA*GAIR
                GAIR=GAR
                GAII=GAI
                IF (Y0.EQ.0) GAII=0.D0
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF (Y0.LT.0.D0) GAII=-GAII
        RETURN
        END SUBROUTINE
        SUBROUTINE  AIRY1(X,W,GAIR,GAII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              0 <= PHASE(Z) <= PI/2                        C
!CC                                                           C
!CC INPUTS:                                                   C
!CC      X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN         C
!CC                QUADRATURE                                 C
!CC OUTPUTS:                                                  C
!CC      GAIR, GAII,  REAL AND IMAGINARY PARTS OF AI(Z)       C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          X,W,GAIR,GAII
        REAL(kp)          PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,&
       &R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,&
       &S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2  !,PHASE
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SUMAR=0.D0
        SUMAI=0.D0
        DO 1 I=1,25
          DF1=1.5D0*X(I)/R32
          DF1C1=DF1*C1
          PHI=PHASE(2._kp+DF1C1,DF1*S1)
          PHI6=PHI/6.D0
          S2=SIN(PHI6)
          C2=COS(PHI6)
          DMODU=SQRT(4.D0+DF1*DF1+4.D0*DF1C1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMAR=SUMAR+W(I)*FUNR
          SUMAI=SUMAI+W(I)*FUNI
 1      CONTINUE
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR+FAC2*SUMAI
        GAII=FAC1*SUMAI-FAC2*SUMAR
        RETURN
        END SUBROUTINE
        SUBROUTINE  AIRY2(X,W,GAIR,GAII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              PI/2 < PHASE(Z) <= 2PI/3                     C
!CC                                                           C
!CC INPUTS:                                                   C
!CC      X,W,        NODES AND WEIGHTS FOR THE GAUSSIAN       C
!CC                  QUADRATURE                               C
!CC OUTPUTS:                                                  C
!CC      GAIR, GAII, REAL AND IMAGINARY PARTS OF AI(Z)        C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          X,W,GAIR,GAII
        REAL(kp)          PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,&
       &R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,&
       &S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2   !,PHASE
        REAL(kp)          SQR2,SQR2I,TAU,TGTAU,B,ANG,CTAU,CFAC,CT,ST,&
       &SUMR,SUMI,TTAU,BETA
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SQR2=1.41421356237310D0
        SQR2I=0.707106781186548D0
        TAU=TH15-PIH3*0.5D0
        TGTAU= TAN(TAU)
        B=5.D0*A
        ANG=TAU*B
        CTAU=COS(TAU)
        CFAC=CTAU**(-B)
        CT=COS(ANG)
        ST=SIN(ANG)
        SUMR=0.D0
        SUMI=0.D0
        DO 2 I=1,25
          DF1=3.D0*X(I)/(CTAU*R32)
          DF1C1=DF1*SQR2I*0.5D0
          PHI=PHASE(2._kp-DF1C1,DF1C1)
          PHI6=PHI/6.D0
          TTAU=X(I)*TGTAU
          BETA=PHI6-TTAU
          S2=SIN(BETA)
          C2=COS(BETA)
          DMODU=SQRT(4.D0+DF1*DF1*0.25D0-SQR2*DF1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMR=SUMR+W(I)*FUNR
          SUMI=SUMI+W(I)*FUNI
 2      CONTINUE
        SUMAR=CFAC*(CT*SUMR-ST*SUMI)
        SUMAI=CFAC*(CT*SUMI+ST*SUMR)
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR+FAC2*SUMAI
        GAII=FAC1*SUMAI-FAC2*SUMAR
        RETURN
        END SUBROUTINE
        SUBROUTINE AIRY1D(X,W,GAIR,GAII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              0 <= PHASE(Z) <= PI/2                         C
!CC                                                            C
!CC INPUTS:                                                    C
!CC       X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN         C
!CC                 QUADRATURE                                 C
!CC OUTPUTS:                                                   C
!CC       GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)        C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          X,W,GAIR,GAII
        REAL(kp)          PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,&
       &R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,&
       &S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2  !,PHASE
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SUMAR=0.D0
        SUMAI=0.D0
        DO 3 I=1,25
          DF1=1.5D0*X(I)/R32
          DF1C1=DF1*C1
          PHI=PHASE(2._kp+DF1C1,DF1*S1)
          PHI6=-PHI*ALF
          S2=SIN(PHI6)
          C2=COS(PHI6)
          DMODU=SQRT(4.D0+DF1*DF1+4.D0*DF1C1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMAR=SUMAR+W(I)*FUNR
          SUMAI=SUMAI+W(I)*FUNI
 3      CONTINUE
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR-FAC2*SUMAI
        GAII=FAC1*SUMAI+FAC2*SUMAR
        RETURN
        END SUBROUTINE
        SUBROUTINE  AIRY2D(X,W,GAIR,GAII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              PI/2 < PHASE(Z) <= 3PI/2                      C
!CC                                                            C
!CC INPUTS:                                                    C
!CC      X,W,   NODES AND WEIGHTS FOR THE GAUSSIAN             C
!CC                   QUADRATURE                               C
!CC OUTPUTS:                                                   C
!CC      GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          X,W,GAIR,GAII
        REAL(kp)          PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,&
       &R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,&
       &S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2  !,PHASE
        REAL(kp)          SQR2,SQR2I,TAU,TGTAU,B,ANG,CTAU,CFAC,CT,ST,&
       &SUMR,SUMI,TTAU,BETA
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SQR2=1.41421356237310D0
        SQR2I=0.707106781186548D0
        TAU=TH15-PIH3*0.5D0
        TGTAU= TAN(TAU)
        B=7.D0*ALF
        ANG=TAU*B
        CTAU=COS(TAU)
        CFAC=CTAU**(-B)
        CT=COS(ANG)
        ST=SIN(ANG)
        SUMR=0.D0
        SUMI=0.D0
        DO 4 I=1,25
          DF1=3.D0*X(I)/(CTAU*R32)
          DF1C1=DF1*SQR2I*0.5D0
          PHI=PHASE(2._kp-DF1C1,DF1C1)
          PHI6=-PHI/6.D0
          TTAU=X(I)*TGTAU
          BETA=PHI6-TTAU
          S2=SIN(BETA)
          C2=COS(BETA)
          DMODU=SQRT(4.D0+DF1*DF1*0.25D0-SQR2*DF1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMR=SUMR+W(I)*FUNR
          SUMI=SUMI+W(I)*FUNI
 4      CONTINUE
        SUMAR=CFAC*(CT*SUMR-ST*SUMI)
        SUMAI=CFAC*(CT*SUMI+ST*SUMR)
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR-FAC2*SUMAI
        GAII=FAC1*SUMAI+FAC2*SUMAR
        RETURN
        END SUBROUTINE
        REAL(kp)          FUNCTION PHASE(X,Y)
        REAL(kp)          PI,PIHAL,X,Y,AY,P
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  COMPUTES THE PHASE OF Z = X + IY, IN (-PI,PI]
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMMON/PARAM1/PI,PIHAL
        IF ((X.EQ.0).AND.(Y.EQ.0)) THEN
          P=0.D0
        ELSE
          AY=ABS(Y)
          IF (X.GE.AY) THEN
            P=ATAN(AY/X)
          ELSEIF ((X+AY).GE.0.D0) THEN
            P=PIHAL-ATAN(X/AY)
          ELSE
            P=PI+ATAN(AY/X)
          ENDIF
          IF (Y.LT.0.D0) P=-P
        ENDIF
        PHASE=P
        END FUNCTION
        SUBROUTINE FGP(X,Y,EPS,FR,FI,GR,GI)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    COMPUTES THE FUNCTIONS F AND G FOR THE SERIES  C
!    OF AI'(Z).                                     C
!    THIS ROUTINE IS CALLED BY SERAID.              C
!                                                   C
!    INPUTS:                                        C
!         X,Y,  REAL AND IMAGINARY PARTS OF Z       C
!         EPS,  PRECISION FOR THE COMPUTATION OF    C
!               THE SERIES                          C
!    OUTPUTS:                                       C
!         FR,FI, REAL AND IMAGINARY PARTS OF F      C
!         GR,GI, REAL AND IMAGINARY PARTS OF G      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          X,Y,EPS,FR,FI,GR,GI
        INTEGER A,B,K3
        REAL(kp)          X2,Y2,U,V,P,Q,CR,CI,DR,DI
        X2=X*X
        Y2=Y*Y
        K3=0
        U=X*(X2-3*Y2)
        V=Y*(3*X2-Y2)
        CR=0.5D0
        CI=0.D0
        DR=1.D0
        DI=0.D0
        FR=0.5D0
        FI=0.D0
        GR=1.D0
        GI=0.D0
 70     A=(K3+5)*(K3+3)
        B=(K3+1)*(K3+3)
        P=(U*CR-V*CI)/A
        Q=(V*CR+U*CI)/A
        CR=P
        CI=Q
        P=(U*DR-V*DI)/B
        Q=(V*DR+U*DI)/B
        DR=P
        DI=Q
        FR=FR+CR
        FI=FI+CI
        GR=GR+DR
        GI=GI+DI
        K3=K3+3
        IF ((ABS(CR)+ABS(DR)+ABS(CI)+ABS(DI)).GE.EPS) GOTO 70
        U=X2-Y2
        V=2.D0*X*Y
        P=U*FR-V*FI
        Q=U*FI+V*FR
        FR=P
        FI=Q
        RETURN
        END SUBROUTINE
        SUBROUTINE FG(X,Y,EPS,FR,FI,GR,GI)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    COMPUTES THE FUNCTIONS F AND G IN EXPRESSION   C
!    10.4.2 OF ABRAMOWITZ & STEGUN FOR THE SERIES   C
!    OF AI(Z).                                      C
!    THIS ROUTINE IS CALLED BY SERAI.               C
!                                                   C
!    INPUTS:                                        C
!          X,Y,  REAL AND IMAGINARY PARTS OF Z      C
!          EPS,  PRECISION FOR THE COMPUTATION      C
!                OF THE SERIES.                     C
!    OUTPUTS:                                       C
!          FR,FI, REAL AND IMAGINARY PARTS OF F     C
!          GR,GI, REAL AND IMAGINARY PARTS OF G     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        INTEGER A,B,K3
        REAL(kp)          X2,Y2,U,V,P,Q,CR,CI,DR,DI
        REAL(kp)          X,Y,EPS,FR,FI,GR,GI
        X2=X*X
        Y2=Y*Y
        K3=0
        U=X*(X2-3.D0*Y2)
        V=Y*(3.D0*X2-Y2)
        CR=1.D0
        CI=0.D0
        DR=1.D0
        DI=0.D0
        FR=1.D0
        FI=0.D0
        GR=1.D0
        GI=0.D0
 71     A=(K3+2)*(K3+3)
        B=(K3+4)*(K3+3)
        P=(U*CR-V*CI)/A
        Q=(V*CR+U*CI)/A
        CR=P
        CI=Q
        P=(U*DR-V*DI)/B
        Q=(V*DR+U*DI)/B
        DR=P
        DI=Q
        FR=FR+CR
        FI=FI+CI
        GR=GR+DR
        GI=GI+DI
        K3=K3+3
        IF ((ABS(CR)+ABS(DR)+ABS(CI)+ABS(DI)).GE.EPS) GOTO 71
        P=X*GR-Y*GI
        Q=X*GI+Y*GR
        GR=P
        GI=Q
        RETURN
        END SUBROUTINE
        SUBROUTINE SERAI(X,Y,AIR,AII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI(Z), TAYLOR, COMPLEX Z      CCC
!CC                                      CCC
!CC   INPUTS:                            CCC
!CC        X,Y,    REAL AND IMAGINARY    CCC
!CC                PARTS OF Z            CCC
!CC   OUTPUTS:                           CCC
!CC        AIR,AII, REAL AND IMAGINARY   CCC
!CC                 PARTS OF AI(Z)       CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          X,Y,EPS,AIR,AII
        REAL(kp)          FZR,FZI,GZR,GZI,CONS1,CONS2
 !       REAL(kp)          D1MACH
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15
        CONS1=0.355028053887817239260D0
        CONS2=0.258819403792806798405D0
        CALL FG(X,Y,EPS,FZR,FZI,GZR,GZI)
        AIR=CONS1*FZR-CONS2*GZR
        AII=CONS1*FZI-CONS2*GZI
        RETURN
        END SUBROUTINE
        SUBROUTINE SERAID(X,Y,AIR,AII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI'(Z), TAYLOR, COMPLEX Z     CCC
!CC                                      CCC
!CC   INPUTS:                            CCC
!CC        X,Y,   REAL AND IMAGINARY     CCC
!CC               PARTS OF Z             CCC
!CC   OUTPUTS:                           CCC
!CC        AIR,AII, REAL AND IMAGINARY   CCC
!CC                 PARTS OF AI'(Z)      CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          X,Y,EPS,AIR,AII
        REAL(kp)          FZR,FZI,GZR,GZI,CONS1,CONS2
 !       REAL(kp)          D1MACH
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15
        CONS1=0.355028053887817239260D0
        CONS2=0.258819403792806798405D0
        CALL FGP(X,Y,EPS,FZR,FZI,GZR,GZI)
        AIR=CONS1*FZR-CONS2*GZR
        AII=CONS1*FZI-CONS2*GZI
        RETURN
        END SUBROUTINE
        SUBROUTINE EXPAI(GAIR,GAII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
!CC                                               CCC
!CC   OUTPUTS:                                    CCC
!CC        GAIR, GAII,  REAL AND IMAGINARY        CCC
!CC                     PARTS OF AI(Z)            CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          EPS,GAIR,GAII
        REAL(kp)          THET,R,TH15,S1,&
       &C1,R32,FACTO,TH025,S3,C3
        REAL(kp)          DF1,PSIIR,PSIII,CK,DFRR,DFII,SUMAR,SUMAI,&
       &DFR,DFI,DELTAR,DELTAI,FAC1,FAC2
        REAL(kp)          CO,DF
   !     REAL(kp)          D1MACH
        INTEGER K
        DIMENSION CO(20)
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SAVE CO
        DATA CO/-6.944444444444445D-2,3.713348765432099D-2,&
       &-3.799305912780064D-2,5.764919041266972D-2,-0.116099064025515D0,&
       &0.291591399230751D0,-0.877666969510017D0,3.07945303017317D0,&
       &-12.3415733323452D0,55.6227853659171D0,-278.465080777603D0,&
       &1533.16943201280D0,-9207.20659972641D0,59892.5135658791D0,&
       &-419524.875116551D0,3148257.41786683D0,-25198919.8716024D0,&
       &214288036.963680D0,-1929375549.18249D0,18335766937.8906D0/
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15
        DF1=1.5D0/R32
        PSIIR=C1
        PSIII=-S1
        K=0
        CK=1.D0
        DF=1.D0
        DFRR=1.D0
        DFII=0.D0
        SUMAR=1.D0
        SUMAI=0.D0
80      DF=DF*DF1
        CK=CO(K+1)*DF
        DFR=DFRR
        DFI=DFII
        DFRR=DFR*PSIIR-DFI*PSIII
        DFII=DFR*PSIII+DFI*PSIIR
        DELTAR=DFRR*CK
        DELTAI=DFII*CK
        SUMAR=SUMAR+DELTAR
        SUMAI=SUMAI+DELTAI
        K=K+1
        IF (SUMAR.NE.0) THEN
          IF (ABS(DELTAR/SUMAR).GT.EPS)  GOTO 80
        ENDIF
        IF (SUMAI.NE.0) THEN
          IF (ABS(DELTAI/SUMAI).GT.EPS) GOTO 80
        ENDIF
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR+FAC2*SUMAI
        GAII=FAC1*SUMAI-FAC2*SUMAR
        RETURN
        END SUBROUTINE
        SUBROUTINE EXPAID(GAIR,GAII)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI'(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
!CC                                                CCC
!CC   OUTPUTS:                                     CCC
!CC        GAIR, GAII,  REAL AND IMAGINARY         CCC
!CC                     PARTS OF AI'(Z)            CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        REAL(kp)          EPS,GAIR,GAII
        REAL(kp)          THET,R,TH15,S1,&
       &C1,R32,FACTO,TH025,S3,C3
        REAL(kp)          DF1,PSIIR,PSIII,VK,DFRR,DFII,SUMAR,SUMAI,&
       &DFR,DFI,DELTAR,DELTAI,FAC1,FAC2
        REAL(kp)          CO,DF
   !     REAL(kp)          D1MACH
        INTEGER K
        DIMENSION CO(20)
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SAVE CO
        DATA CO/9.722222222222222D-2,-4.388503086419753D-2,&
       &4.246283078989484D-2,-6.266216349203230D-2,&
       &0.124105896027275D0,-0.308253764901079D0,&
       &0.920479992412945D0,-3.21049358464862D0,&
       &12.8072930807356D0,-57.5083035139143D0,&
       &287.033237109221D0,-1576.35730333710D0,&
       &9446.35482309593D0,-61335.7066638521D0,&
       &428952.400400069D0,-3214536.52140086D0,&
       &25697908.3839113D0,-218293420.832160D0,&
       &1963523788.99103D0,-18643931088.1072D0/
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15
        DF1=1.5D0/R32
        PSIIR=C1
        PSIII=-S1
        K=0
        DF=1.D0
        DFRR=1.D0
        DFII=0.D0
        SUMAR=1.D0
        SUMAI=0.D0
 81     DF=DF*DF1
        VK=CO(K+1)*DF
        DFR=DFRR
        DFI=DFII
        DFRR=DFR*PSIIR-DFI*PSIII
        DFII=DFR*PSIII+DFI*PSIIR
        DELTAR=DFRR*VK
        DELTAI=DFII*VK
        SUMAR=SUMAR+DELTAR
        SUMAI=SUMAI+DELTAI
        K=K+1
        IF (SUMAR.NE.0) THEN
          IF (ABS(DELTAR/SUMAR).GT.EPS)  GOTO 81
        ENDIF
        IF (SUMAI.NE.0) THEN
          IF (ABS(DELTAI/SUMAI).GT.EPS) GOTO 81
        ENDIF
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=-(FAC1*SUMAR-FAC2*SUMAI)
        GAII=-(FAC1*SUMAI+FAC2*SUMAR)
        RETURN
        END SUBROUTINE

       SUBROUTINE BIZ(IFUN,IFAC,X0,Y0,GBIR,GBII,IERRO)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COMPUTATION OF THE AIRY FUNCTION BI(Z) OR ITS DERIVATIVE BI'(Z)
! THE CODE USES THE CONNECTION OF BI(Z) WITH AI(Z).
!                BIZ CALLS THE ROUTINE AIZ
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  INPUTS:
!    IFUN:
!         * IFUN=1, THE CODE COMPUTES BI(Z)
!         * IFUN=2, THE CODE COMPUTES BI'(Z)
!    IFAC:
!         * IFAC=1, THE CODE COMPUTES  BI(Z) OR BI'(Z)
!         * IFAC=2, THE CODE COMPUTES NORMALIZED BI(Z) OR BI'(Z)
!    X0:   REAL PART OF THE ARGUMENT Z
!    Y0:   IMAGINARY PART OF THE ARGUMENT  Z
!
!  OUTPUTS:
!    GBIR: REAL PART OF BI(Z) OR BI'(Z)
!    GBII: IMAGINARY PART OF BI(Z) OR BI'(Z)
!
!    IERRO: ERROR FLAG
!          * IERRO=0, SUCCESSFUL COMPUTATION
!          * IERRO=1, COMPUTATION OUT OF RANGE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   ACCURACY:
!
!     1) SCALED AIRY FUNCTIONS:
!        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO
!        THE ZEROS, WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
!        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z|>1000
!        IN THE CASE OF PHASE(Z) CLOSE TO +3*PI/2 OR -3*PI/2.
!     2) UNSCALED AIRY FUNCTIONS:
!        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR
!        3/2*|Z|**(3/2)>LOG(OVER).
!        FOR |Z|<30:
!        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE
!           ZEROS) BETTER THAN 10**(-13).
!        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER
!           THAN 10**(-13), WHERE R(Z)=REAL(BI)/IMAG(BI)
!           OR R(Z)=REAL(BI')/IMAG(BI').
!        FOR |Z|>30, GRADUAL LOSS OF PRECISION TAKES PLACE
!        AS |Z| INCREASES.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     AUTHORS:
!        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN).
!                      E-MAIL: AMPARO.GIL@UAM.ES
!        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
!                      E-MAIL: JSEGURA@MATH.UC3M.ES
!        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
!                      E-MAIL: NICO.TEMME@CWI.NL
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    REFERENCES:
!         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
!         NUMERICAL ALGORITHMS (2002).
!         A. GIL, J. SEGURA, N.M. TEMME
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       REAL(kp)          X0,Y0,GBIR,GBII
       REAL(kp)          OVER,UNDER,DL1,DL2,COVER  !,D1MACH
       REAL(kp)          PI,PI3,PI23,SQRT3,C,S,C1,S1,U,V,X,Y,AR,AI,&
      &APR,API,BR,BI,BPR,BPI,BBR,BBI,BBPR,BBPI  !,PHASE
       REAL(kp)          THET,R,R32,THET32,A1,B1,DF1,EXPO,EXPOI,&
      &ETAR,ETAI,ETAGR,ETAGI,PIHAL
       INTEGER IFUN,IFAC,IEXPF,IERR,IERRO
       COMMON/PARAM1/PI,PIHAL
       SQRT3=1.7320508075688772935D0
       PI=3.1415926535897932385D0
       PIHAL=1.5707963267948966192D0
       PI3=PI/3.D0
       PI23=2.D0*PI3
       X=X0
       C=0.5D0*SQRT3
       S=0.5D0
       IERRO=0
       IEXPF=0
       IF (Y0.LT.0.D0) THEN
         Y=-Y0
       ELSE
         Y=Y0
       ENDIF
       R=SQRT(X*X+Y*Y)
       R32=R*SQRT(R)
       THET=PHASE(X,Y)
       COVER=2.D0/3.D0*R32*ABS(COS(1.5D0*THET))
!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
       OVER=D1MACH(2)*1.D-4
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
       UNDER=D1MACH(1)*1.D+4
       DL1=LOG(OVER)
       DL2=-LOG(UNDER)
       IF (DL1.GT.DL2) OVER=1/UNDER
       IF (IFAC.EQ.1) THEN
         IF (COVER.GE.LOG(OVER)) THEN
!CC OVERFLOW/UNDERFLOW PROBLEMS.
!CC   CALCULATION ABORTED
           IERRO=1
           GBIR=0
           GBII=0
         ENDIF
       ELSE
         IF (COVER.GE.(LOG(OVER)*0.2)) IEXPF=1
       ENDIF
       IF (IERRO.EQ.0) THEN
         IF (IFAC.EQ.1) THEN
           IF (Y.EQ.0.D0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  TAKE TWICE THE REAL PART OF EXP(-PI I/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             C1=-0.5D0
             S1=-0.5D0*SQRT3
             U=X*C1-Y*S1
             V=X*S1+Y*C1
             CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
             IF (IFUN.EQ.1) THEN
               BR=SQRT3*AR+AI
               BI=0.D0
             ELSE
               U=AR*C1-AI*S1
               V=AR*S1+AI*C1
               APR=U
               API=V
               BPR=SQRT3*APR+API
               BPI=0.D0
             ENDIF
           ELSE
             IF ((X.LT.0.D0).AND.(Y.LT.-X*SQRT3)) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC      2 PI/3  < PHASE(Z) < PI
!CC      BI(Z)=EXP(I PI/6) AI_(-1)(Z) + EXP(-I PI/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               C1=-0.5D0
               S1=0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BR=C*AR-S*AI
                 BI=C*AI+S*AR
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 BPR=C*APR-S*API
                 BPI=C*API+S*APR
               ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   WE NEED ALSO AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               C1=-0.5D0
               S1=-0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 S=-S
                 BR=BR+C*AR-S*AI
                 BI=BI+C*AI+S*AR
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 S=-S
                 BPR=BPR+C*APR-S*API
                 BPI=BPI+C*API+S*APR
               ENDIF
             ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   BI(Z) = I AI(Z) + 2 EXP(-I PI/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               C1=-0.5D0
               S1=-0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BR=SQRT3*AR+AI
                 BI=-AR+SQRT3*AI
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 BPR=SQRT3*APR+API
                 BPI=-APR+SQRT3*API
               ENDIF
               CALL AIZ(IFUN,IFAC,X,Y,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BR=BR-AI
                 BI=BI+AR
               ELSE
                 BPR=BPR-AI
                 BPI=BPI+AR
               ENDIF
             ENDIF
           ENDIF
         ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         SCALED BI AIRY FUNCTIONS         C
!CC   WE USE THE FOLLOWING NORMALIZATION:    C
!CC   LET ARGZ=ARG(Z), THEN:                 C
!CC   A) IF  0 <= ARGZ <= PI/3               C
!CC      BI=EXP(-2/3Z^3/2)BI                 C
!CC   B) IF  PI/3 <= ARGZ <= PI              C
!CC      BI=EXP(2/3Z^3/2)BI                  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           THET=PHASE(X,Y)
           IF (THET.LE.PI3) THEN
             C1=-0.5D0
             S1=-0.5D0*SQRT3
             U=X*C1-Y*S1
             V=X*S1+Y*C1
             CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
             IF (IFUN.EQ.1) THEN
               BR=SQRT3*AR+AI
               BI=-AR+SQRT3*AI
             ELSE
               U=AR*C1-AI*S1
               V=AR*S1+AI*C1
               APR=U
               API=V
               BPR=SQRT3*APR+API
               BPI=-APR+SQRT3*API
             ENDIF
             IF (IEXPF.EQ.0) THEN
               R=SQRT(X*X+Y*Y)
               R32=R*SQRT(R)
               THET32=THET*1.5D0
               A1=COS(THET32)
               B1=SIN(THET32)
               DF1=4.D0/3.D0*R32
               EXPO=EXP(DF1*A1)
               EXPOI=1.D0/EXPO
               ETAR=EXPO*COS(DF1*B1)
               ETAI=EXPO*SIN(DF1*B1)
               ETAGR=EXPOI*COS(-DF1*B1)
               ETAGI=EXPOI*SIN(-DF1*B1)
               CALL AIZ(IFUN,IFAC,X,Y,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BR=BR-AR*ETAGI-ETAGR*AI
                 BI=BI+AR*ETAGR-ETAGI*AI
               ELSE
                 BPR=BPR-AR*ETAGI-ETAGR*AI
                 BPI=BPI+AR*ETAGR-ETAGI*AI
               ENDIF
             ENDIF
           ENDIF
           IF ((THET.GT.PI3).AND.(THET.LE.PI23)) THEN
             IF (IEXPF.EQ.0) THEN
               R=SQRT(X*X+Y*Y)
               R32=R*SQRT(R)
               THET32=THET*1.5D0
               A1=COS(THET32)
               B1=SIN(THET32)
               DF1=4.D0/3.D0*R32
               EXPO=EXP(DF1*A1)
               EXPOI=1.D0/EXPO
               ETAR=EXPO*COS(DF1*B1)
               ETAI=EXPO*SIN(DF1*B1)
               ETAGR=EXPOI*COS(-DF1*B1)
               ETAGI=EXPOI*SIN(-DF1*B1)
               C1=-0.5D0
               S1=-0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BBR=SQRT3*AR+AI
                 BBI=-AR+SQRT3*AI
                 BR=BBR*ETAR-BBI*ETAI
                 BI=BBI*ETAR+BBR*ETAI
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 BBPR=SQRT3*APR+API
                 BBPI=-APR+SQRT3*API
                 BPR=BBPR*ETAR-BBPI*ETAI
                 BPI=BBPI*ETAR+BBPR*ETAI
               ENDIF
             ELSE
               IF (IFUN.EQ.1) THEN
                 BR=0.D0
                 BI=0.D0
               ELSE
                 BPR=0.D0
                 BPI=0.D0
               ENDIF
             ENDIF
             CALL AIZ(IFUN,IFAC,X,Y,AR,AI,IERR)
             IF (IFUN.EQ.1) THEN
               BR=BR-AI
               BI=BI+AR
             ELSE
               BPR=BPR-AI
               BPI=BPI+AR
             ENDIF
           ENDIF
           IF (THET.GT.PI23) THEN
             C1=-0.5D0
             S1=0.5D0*SQRT3
             U=X*C1-Y*S1
             V=X*S1+Y*C1
             CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
             IF (IFUN.EQ.1) THEN
               BR=C*AR-S*AI
               BI=C*AI+S*AR
             ELSE
               U=AR*C1-AI*S1
               V=AR*S1+AI*C1
               APR=U
               API=V
               BPR=C*APR-S*API
               BPI=C*API+S*APR
             ENDIF
             IF (IEXPF.EQ.0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   WE NEED ALSO AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                R=SQRT(X*X+Y*Y)
                R32=R*SQRT(R)
                THET32=THET*1.5D0
                A1=COS(THET32)
                B1=SIN(THET32)
                DF1=4.D0/3.D0*R32
                EXPO=EXP(DF1*A1)
                EXPOI=1.D0/EXPO
                ETAR=EXPO*COS(DF1*B1)
                ETAI=EXPO*SIN(DF1*B1)
                ETAGR=EXPOI*COS(-DF1*B1)
                ETAGI=EXPOI*SIN(-DF1*B1)
                C1=-0.5D0
                S1=-0.5D0*SQRT3
                U=X*C1-Y*S1
                V=X*S1+Y*C1
                CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
                IF (IFUN.EQ.1) THEN
                  S=-S
                  BBR=C*AR-S*AI
                  BBI=C*AI+S*AR
                  BR=BR+ETAR*BBR-ETAI*BBI
                  BI=BI+BBI*ETAR+ETAI*BBR
                ELSE
                  U=AR*C1-AI*S1
                  V=AR*S1+AI*C1
                  APR=U
                  API=V
                  S=-S
                  BBPR=C*APR-S*API
                  BBPI=C*API+S*APR
                  BPR=BPR+ETAR*BBPR-ETAI*BBPI
                  BPI=BPI+BBPI*ETAR+ETAI*BBPR
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          IF (Y0.LT.0) THEN
            BI=-BI
            BPI=-BPI
          ENDIF
          IF (IFUN.EQ.1) THEN
            GBIR=BR
            GBII=BI
          ELSE
            GBIR=BPR
            GBII=BPI
          ENDIF
        ENDIF
        RETURN
        END SUBROUTINE
! The end of Algorithm 819.

! Algorithm 421
! Algorithm 421 was copied from Ref. (6).
! Algorithm 421 is written originally in FORTRAN 77. Therefore, Algorithm 421 
! was translated into Fortran 95 by the use of the transformer program stored
! in file transformer.f95.
! In addition, Algorithm 421 was a little modified for increasing portability of
! the program.
! The lines marked with !! are inserted or modified by M. Kodama.
!                                                        M. Kodama
!
!*************************************************************************
!   Hirondo Kuki, "Algorithm 421
!   Complex gamma function with error control,"
!   Communications of the ACM, April 1972, Volume 15, Number 4, pp. 271-272.
!
      SUBROUTINE CDLGAM(CARG,CANS,ERROR,LF0)
      INTEGER,INTENT(IN):: LF0                      !!
      REAL(kp),INTENT(INOUT):: ERROR                    !!
! COMPLEX GAMMA AND LOGGAMMA FUNCTIONS WITH ERROR ESTIMATE
!
! CARG = A COMPLEX ARGUMENT, GIVEN AS A VECTOR OF 2 DOUBLE
!        PRECISION ELEMENTS CONSISTING OF THE REAL COMPONENT
!        FOLLOWED BY THE IMAGINARY COMPONENT
! CANS = THE COMPLEX ANSWER, OF THE SAME TYPE AS CARG
! ERROR = A REAL VARIABLE.  IT STANDS FOR AN ESTIMATE OF THE
!        ABSOLUTE ERROR OF THE ARGUMENT AS INPUT.  AS OUTPUT
!        IT GIVES AN ESTIMATE OF THE ABSOLUTE (FOR LOGGAMMA)
!        OR THE RELATIVE (FOR GAMMA) ERROR OF THE ANSWER
! LF0  = FLAG.  SET IT TO 0 FOR LOGGAMMA, AND 1 FOR GAMMA
      REAL(kp)          CARG(2),CANS(2),COEF(7),F0,F1,G0,G1,&
         &PI,DPI,HL2P,AL2P,DELTA,DE0,DE1,Z1,Z2,ZZ1,W1,W2,Y1,&
         &A,B,U,U1,U2,UU1,UU2,UUU1,UUU2,V1,V2,VV1,VV2,T1,T2,&
         &H,H1,H2,AL1,AL2,DN,EPS,OMEGA
      DATA COEF(1)/+0.641025641025641026D-2/
      DATA COEF(2)/-0.191752691752691753D-2/
      DATA COEF(3)/+0.841750841750841751D-3/
      DATA COEF(4)/-0.595238095238095238D-3/
      DATA COEF(5)/+0.793650793650793651D-3/
      DATA COEF(6)/-0.277777777777777778D-2/
      DATA COEF(7)/+0.833333333333333333D-1/
      DATA F0/840.07385296052619D0/,F1/20.001230821894200D0/
      DATA G0/1680.1477059210524D0/,G1/180.01477047052042D0/
      DATA PI/3.14159265358979324D0/
      DATA DPI/6.28318530717958648D0/
      DATA HL2P/0.918938533204672742D0/
      DATA AL2P/1.83787706640934548D0/
! CONSTANTS EPS AND OMEGA ARE MACHINE DEPENDENT.
! EPS   IS THE BASIC ROUND-OFF UNIT.  FOR S/360 MACHINES,
! IT IS CHOSEN TO BE 16**-13.  FOR BINARY MACHINE OF N-
! BIT ACCURACY SET IT TO 2**(-N+1), AND INITIALIZE DE0
! AS 5.0 RATHER THAN AS 2.0
! OMEGA   IS THE LARGEST NUMBER REPRESENTABLE BY THE FLOAT
! POINT REPRESENTATION OF THE MACHINE.  FOR S/360
!     MACHINES, IT IS SLIGHTLY LESS THAN 16**63.
!!      DATA EPS/2.20D-16/
!!      DATA OMEGA/7.23700538D75/
      PARAMETER(OMEGA=HUGE(1._kp))              !!
      INTEGER, PARAMETER:: IRADIX=RADIX(1D0)    !!
      EPS=epsilon1                              !!
      Z1 = CARG(1)
      Z2 = CARG(2)
      DELTA = ABS(ERROR)
      DE0 = 2.0D0
      IF(IRADIX == 2) DE0 = 5.0D0               !!
      DE1 = 0.0
! FORCE SIGN OF IMAGINARY PART OF ARG TO NON-NEGATIVE
      LF1 = 0
      IF (Z2 .GE. 0.0) GO TO 20
      LF1 = 1
      Z2 = -Z2
   20 LF2 = 0
      IF (Z1 .GE. 0.0) GO TO 100
! CASE WHEN REAL PART OF ARG IS NEGATIVE
      LF2 = 1
      LF1 = LF1-1
      T1 = AL2P - PI*Z2
      T2 = PI*(0.5D0 - Z1)
      U  = -DPI*Z2
      IF (U .GE. -0.1054D0) GO TO 40
      A  = 0.0D0
! IF E**U .LT. 10**(-17), IGNOR IT TO SAVE TIME AND TO AVOID
! IRRELEVANT UNDERFLOW
      IF (U .LE. -39.15D0) GO TO 30
      A =  EXP(U)
   30 H1 = 1.0D0 - A
      GO TO 50
   40 U2 = U*U
      A = -U*(F1*U2 + F0)
      H1 = (A + A)/((U2 + G1)*U2 + G0 + A)
      A = 1.0D0 - H1
! DINT IS THE DOUBLE PRECISION VERSION OF AINT, INTEGER EX-
!   TRACTION.  THIS FUNCTION IS NOT INCLUDED IN ANSI FORTRAN
!   .  WHEN THIS FUNCTION IS NOT PROVIDED BY THE SYSTEM,
!   EITHER SUPPLY IT AS AN EXTERNAL SUBROUTINE (AND TYPE THE
!   NAME DINT AS DOUBLE PRECISION), OR MODIFY THE NEXT
!   STATEMENT AS THE EXAMPLE FOR S/340 INDICATES.  FOR S/360
!   REPLACE IT WITH
!
!      DOUBLE PRECISION SCALE
!      DATA SCALE/Z4F00000000000000/
!  50 B = Z1 - ((Z1 - 0.5D0) + SCALE)
   50 B = Z1 - AINT(Z1 - 0.5D0)
      H2 = A* SIN(DPI*B)
      B  =  SIN(PI*B)
      H1 = H1 + (B+B)*B*A
      H =  ABS(H2) + H1 - DPI*A*DELTA
      IF (H .LE. 0.0) GO TO 500
      DE0 = DE0 +  ABS(T1) + T2
      DE1 = PI + DPI*A/H
      Z1 = 1.0D0 - Z1
! CASE WHEN NEITHER REAL PART NOR IMAGINARY PART OF ARG IS
! NEGATIVE.  DEFINE THERSHOLD CURVE TO BE THE BROKEN LINES
! CONNECTING POINTS 10F0*I, 10F4.142*I, 0.1F14.042*I,AND
! 0.1FOMEGA*I
  100 LF3 = 0
      Y1 = Z1 - 0.5D0
      W1 = 0.0
      W2 = 0.0
      K  = 0
      B  =   MAX(0.1_kp,   MIN(10.0_kp, 14.142_kp-Z2)) - Z1    !!
      IF (B .LE. 0.0) GO TO 200
! CASE WHEN REAL PART OF ARG IS BETWEEN 0 AND THRESHOLD
      LF3 = 1
      ZZ1 = Z1
      N  = B + 1.0D0
      DN = N
      Z1 = Z1 + DN
      A  = Z1*Z1 + Z2*Z2
      V1 = Z1/A
      V2 = -Z2/A
! INITIALIZE U1+U2*I AS THE RIGHTMOST FACTOR 1-1/(Z+M)
      U1 = 1.0D0 - V1
      U2 = -V2
      K  = 6.0D0 - Z2*0.6D0 - ZZ1
      IF (K .LE. 0) GO TO 120
! FORWORD ASSEMBLY OF FACTORS (Z+J-1)/(Z+M)
      N  = N - K
      UU1 = (ZZ1*Z1 + Z2*Z2) / A
      UU2 = DN*Z2/A
      VV1 = 0.0
      VV2 = 0.0
      DO 110 J = 1,K
        B  = U1*(UU1+VV1) - U2*(UU2+VV2)
        U2 = U1*(UU2+VV2) + U2*(UU1+VV1)
        U1 = B
        VV1 = VV1 + V1
        VV2 = VV2 + V2
  110   CONTINUE
  120 IF (N .LE. 1) GO TO 140
! BACKWARD ASSEMBLY OF FACTORS 1-J/(Z+N)
      VV1 = V1
      VV2 = V2
      DO 130  J = 2,N
        VV1 = VV1 + V1
        VV2 = VV2 + V2
        B  = U1*(1.0D0 - VV1) + U2*VV2
        U2 = -U1*VV2 + U2*(1.0D0 - VV1)
        U1 = B
  130   CONTINUE
  140 U  = U1*U1 + U2*U2
      IF (U .EQ. 0.0) GO TO 500
      IF (LF0 .EQ. 0) GO TO 150
      IF (K .LE. 0) GO TO 200
  150 AL1 =  LOG(U)*0.5D0
      IF (LF0 .NE. 0) GO TO 160
      W1 = AL1
      W2 =  ATAN2(U2,U1)
      IF (W2 .LT. 0.0) W2 = W2 + DPI
      IF (K .LE. 0) GO TO 200
  160 A = ZZ1 + Z2 - DELTA
      IF (A .LT. 0.0) GO TO 500
      DE0 = DE0 - AL1
      DE1 = DE1 + 2.0D0 + 1.0D0/A
! CASE WHEN REAL PART OF ARG IS GREATER THAN THRESHOLD
  200 A = Z1*Z1 + Z2*Z2
      AL1 =  LOG(A)*0.5D0
      AL2 =  ATAN2(Z2,Z1)
      V1 = Y1*AL1 - Z2*AL2
      V2 = Y1*AL2 + Z2*AL1
! EVALUATE ASYMTOTIC TERMS.  IGNORE THIS TERM,IF ABS VAL(ARG) .GT.
! 10**9, TO SAVE TIME AND TO AVOID IRRELEVANT UNDERFLOW
      VV1 = 0.0
      VV2 = 0.0
      IF (A .GT. 1.0D18) GO TO 220
      UU1 = Z1/A
      UU2 = -Z2/A
      UUU1 = UU1*UU1 - UU2*UU2
      UUU2 = UU1*UU2*2.0D0
      VV1 = COEF(1)
      DO 210  J = 2,7
        B  = VV1*UUU1 - VV2*UUU2
        VV2 = VV1*UUU2 + VV2*UUU1
        VV1 = B + COEF(J)
  210   CONTINUE
      B  = VV1*UU1 -VV2*UU2
      VV2 = VV1*UU2 + VV2*UU1
      VV1 = B
  220 W1 = (((VV1 + HL2P) - W1) - Z1) + V1
      W2 = ((VV2 - W2) -Z2) + V2
      DE0 = DE0 +  ABS(V1) +  ABS(V2)
      IF (K .LE. 0) DE1 = DE1 + AL1
! FINAL ASSEMBLY
      IF (LF2 .NE. 0) GO TO 310
      IF (LF0 .EQ. 0) GO TO 400
      A =  EXP(W1)
      W1 = A* COS(W2)
      W2 = A* SIN(W2)
      IF (LF3 .EQ. 0) GO TO 400
      B  = (W1*U1 + W2*U2) / U
      W2 = (W2*U1 - W1*U2) / U
      W1 = B
      GO TO 400
  310 H = H1*H1 + H2*H2
      IF (H .EQ. 0.0) GO TO 500
      IF (LF0 .EQ. 0) GO TO 320
      IF (H .GT. 1.0D-2) GO TO 330
  320 A =  LOG(H)*0.5D0
      IF (H .LE. 1.0D-2) DE0 = DE0 - A
      IF (LF0 .NE. 0) GO TO 330
      W1 = (T1 - A) - W1
      W2 = (T2 -  ATAN2(H2,H1)) - W2
      GO TO 400
  330 T1 = T1 - W1
      T2 = T2 - W2
      A  =  EXP(T1)
      T1 = A* COS(T2)
      T2 = A* SIN(T2)
      W1 = (T1*H1 + T2*H2)/H
      W2 = (T2*H1 - T1*H2)/H
      IF (LF3 .EQ. 0) GO TO 400
      B  = W1*U1 - W2*U2
      W2 = W1*U2 + W2*U1
      W1 = B
  400 IF (LF1 .NE. 0) W2 = -W2
! TRUNCATION ERROR OF STIRLINGS FORMULA IS UP TO 3*10**-17.
      DE1 = DE0*EPS + 3.0D-17 + DE1*DELTA
      GO TO 600
! CASE WHEN ARGUMENT IS TOO CLOSE TO A SINGURARITY
!
  500 W1 = OMEGA
      W2 = OMEGA
      DE1 = OMEGA
!
  600 CANS(1) = W1
      CANS(2) = W2
      ERROR = DE1
      RETURN
      END SUBROUTINE

 REAL(kp) FUNCTION abs2(za)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: za
! FUNCTION abs2 calculates a rough absolute value of the complex number za.
! This is invoked by bessel, neumann, hankel2, num_region, bes_series,
! neu_series, neu_srs_init, def_bessel, bes_recur, han2_temme, bes_han_dby, 
! bes_olver, han2_olver, sumaabb, fzeta.
 abs2=ABS(REAL(za,kp))+ABS(AIMAG(za))
 END FUNCTION abs2

 SUBROUTINE bessel(znu,xx,zans,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zans
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE computes zans.  zans=zbessel(znu,xx).
! This is invoked only by a user's program.
! The dummy arguments znu, xx, zans, info and the outline of the calculation
! in SUBROUTINE bessel are explained in the beginning of this file.
 COMPLEX(kp):: zbesa,zepspi,zhan2a,zlogbes,zlogmu,zneua,znua,zconnu, &
   zconnupi,znupi,zss1,zss2,zsum,zconepspi
 REAL(kp):: aim_znu,re_znu,error,error1,error2
 INTEGER:: info1,info2,nregion,nn
! Determination of nregion
 nregion=num_region(znu,xx)
 zans=0
 IF(nregion == 0) THEN;   info=30;   RETURN;  ENDIF
 re_znu=REAL(znu,kp);  aim_znu=AIMAG(znu)
 znupi=znu*pi;         zconnu=CONJG(znu)
 zconnupi=zconnu*pi;   znua=-znu
 IF(re_znu < -ihuge) THEN;   info=30;     RETURN;  ENDIF
 IF(re_znu < 0) THEN
   nn=NINT(re_znu)
   zepspi=(znu-nn)*pi
   zconepspi=CONJG(zepspi)
 ENDIF
 SELECT CASE(nregion)
 CASE(1)
   CALL bes_series(znu,xx,zsum,zlogbes,error,info)
   IF(info > 17) THEN;    zans=huge1;    RETURN;   ENDIF
   IF(REAL(zlogbes) > hugelog) THEN;   zans=huge1;   info=20;   RETURN;   ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zsum*EXP(zlogbes)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(2)
   IF(re_znu >= 0) THEN
     CALL bes_series(znu,xx,zsum,zlogbes,error,info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(REAL(zlogbes) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=zsum*EXP(zlogbes)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0, zbessel(znu,xx) is calculated using Eq. (40a) of Ref. (1).
   CALL bes_series(znua,xx,zsum,zlogbes,error,info1)
   IF(info1 > 17) THEN;   zans=huge1;   info=info1;   RETURN;   ENDIF
   IF(REAL(zlogbes) > hugelog) THEN;   zans=huge1;   info=20;    RETURN;   ENDIF
   IF(error > err_range3) THEN;   info=30;     RETURN;   ENDIF
   zbesa=zsum*EXP(zlogbes)
   CALL neu_series(znua,xx,zneua,info2)
   IF(info2 > 17) THEN;   zans=huge1;   info=info2;   RETURN;   ENDIF
   info=MAX(info1,info2)
   error=abs2(znupi)+error
   IF(error > err_range3) THEN;  info=30;  RETURN;   ENDIF
   zans=zbesa*COS(znupi)+zneua*SIN(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
   IF(error > err_range1) info=5
 CASE(3)  !  by Debye's method       |znu/xx|<=1
! zbessel(znu,xx) is calculated using Eq. (40b) of Ref. (1).
   CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;     zans=huge1;     RETURN;   ENDIF
   IF(ABS(REAL(zlogmu)) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=(zss1*EXP(zlogmu)+zss2*EXP(-zlogmu))/2
   IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
   IF(error > err_range1) info=5
 CASE(4)        !  by Debye's method       |znu/xx|>1
   IF(re_znu >= 0) THEN
     IF(aim_znu >= 0) THEN
       CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
       IF(info > 17) THEN;     zans=huge1;     RETURN;       ENDIF
       IF(REAL(zlogmu) > hugelog) THEN;  zans=huge1;  info=20;   RETURN;  ENDIF
       IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
       zans=zss1*EXP(zlogmu)/2
       IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
       IF(error > err_range1) info=5
       RETURN
     ENDIF
! When re_znu>=0 and aim_znu<0, zbessel(znu,xx) is calculated 
! using Eq. (40c) of Ref. (1).
     CALL bes_han_dby(zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;   zans=huge1;    RETURN;    ENDIF
     IF(REAL(zlogmu) > hugelog) THEN;   zans=huge1;   info=20;   RETURN;   ENDIF
     IF(error > err_range3) THEN;    info=30;     RETURN;     ENDIF
     zans=CONJG(zss1*EXP(zlogmu))/2
     IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zbessel(znu,xx) is calculated 
! using Eq. (40d) of Ref. (1).
     CALL bes_han_dby(-zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;    zans=huge1;     RETURN;     ENDIF
     IF(AIMAG(zconnupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
     IF(MAX(REAL(zlogmu-zunit*zconnupi), &
              REAL(-zlogmu)+ABS(AIMAG(zconepspi))) > hugelog) THEN
       zans=huge1;    info=20;    RETURN
     ENDIF
     error=error+abs2(znupi)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=EXP(zlogmu-zunit*zconnupi)*zss1/2 &
          +zunit*zss2*EXP(-zlogmu)*SIN(zconnupi)
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zbessel(znu,xx) is calculated by 
! using Eq. (40e) of Ref. (1).
   CALL bes_han_dby(znua,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(MAX(REAL(zlogmu+zunit*zepspi), &
           REAL(-zlogmu)+ABS(AIMAG(zepspi))) > hugelog) THEN
     zans=huge1;     info=20;     RETURN
   ENDIF
   error=error+abs2(znupi)
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=EXP(zlogmu-zunit*znupi)*zss1/2 &
        +zunit*zss2*EXP(-zlogmu)*SIN(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(5)        !  Method by the Olver's asymptotic series
   IF(re_znu >= 0) THEN
     CALL bes_olver(znu,xx,zans,error,info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zbessel(znu,xx) is calculated by 
! using Eq. (40d) of Ref. (1).
     CALL bes_olver(-zconnu,xx,zbesa,error1,info1)
     CALL han2_olver(-zconnu,xx,zhan2a,error2,info2)
     info=MAX(info1,info2)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     error=ABS(REAL(znupi,kp))+error2
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(ABS(AIMAG(zconepspi)) > hugelog) THEN;  info=30;  RETURN;  ENDIF
     zans=zbesa*EXP(-zunit*zconnupi)+zunit*zhan2a*SIN(zconnupi)
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zbessel(znu,xx) is calculated by 
! using Eq. (40e) of Ref. (1).
   CALL bes_olver(znua,xx,zbesa,error1,info1)
   CALL han2_olver(znua,xx,zhan2a,error2,info2)
   info=MAX(info1,info2)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   error=ABS(REAL(znupi,kp))+error2
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=EXP(-zunit*znupi)*zbesa+zunit*zhan2a*SIN(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(6)
   CALL bes_recur(znu,xx,0,zans,info)
 CASE(7)
   IF(re_znu >= 0) THEN
     CALL bes_recur(znu,xx,0,zans,info)
     RETURN
   ENDIF
! When re_znu<0, zbessel(znu,xx) is calculated using Eq. (40e) of Ref. (1).
   CALL bes_recur(znua,xx,1,zbesa,info1)
   CALL han2_temme(znua,xx,zhan2a,info2)
   info=MAX(info1,info2)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zbesa+zunit*SIN(znupi)*zhan2a
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 END SELECT
 END SUBROUTINE bessel

 SUBROUTINE neumann(znu,xx,zans,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zans
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zans.  zans=zneumann(znu,xx).
! This is invoked only by a user's program.
! The dummy arguments znu, xx, zans, info and the outline of the calculation
! in SUBROUTINE neumann are explained in the beginning of this file.
 COMPLEX(kp):: zarg1,zarg2,zbes,zbesa,zhan2,zhan2a,zlogbes1, &
  zlogbes2,zlogmu,zneua,znua,zconnu,zconnupi,znupi,zpart1,zpart2, &
  zss1,zss2,zsum1,zsum2,z1
 REAL(kp):: aim_znu,re_znu,abszpart1,abszpart2
 REAL(kp):: error,error1,error2,abes,ahan2
 INTEGER:: info1,info2,nregion
! Determination of nregion
 nregion=num_region(znu,xx)
 zans=0
 IF(nregion == 0) THEN;  info=30;   RETURN;   ENDIF
 re_znu=REAL(znu,kp);  aim_znu=AIMAG(znu)
 znupi=znu*pi;         zconnu=CONJG(znu)
 zconnupi=zconnu*pi;   znua=-znu
 SELECT CASE(nregion)
 CASE(1)        !  With Taylor's expansion.
   CALL bes_series(znu,xx,zsum1,zlogbes1,error1,info1)
   IF(info1 > 17) THEN;     zans=huge1;     info=info1;     RETURN;   ENDIF
   IF(REAL(zlogbes1) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   CALL bes_series(znua,xx,zsum2,zlogbes2,error2,info2)
   IF(info2 > 17) THEN;  zans=huge1;   info=info2;     RETURN;   ENDIF
   IF(aim_znu >= 0) THEN
! When aim_znu>=0, zneumann(znu,xx) is calculated using Eq. (41b) of Ref. (1).
     zarg1=zlogbes1+2*zunit*znupi
     zarg2=zlogbes2+zunit*znupi
     IF(REAL(zarg1)>hugelog .OR. REAL(zarg2)>hugelog) THEN
       zans=huge1;     info=20;     RETURN;     ENDIF
     IF((ABS(AIMAG(zarg1))>theta_lim) .OR. &
          (ABS(AIMAG(zlogbes1))>theta_lim)) THEN;  info=30;  RETURN;  ENDIF
     zpart1=zsum1*(EXP(zarg1)+EXP(zlogbes1))
     zpart2=-2*zsum2*EXP(zarg2)
     zans=zunit*(zpart1+zpart2)/(EXP(2*zunit*znupi)-1)
   ELSE
! When aim_znu<0, zneumann(znu,xx) is calculated using Eq. (41c) of Ref. (1).
     zarg1=zlogbes1-2*zunit*znupi
     zarg2=zlogbes2-zunit*znupi
     IF(REAL(zarg1)>hugelog .OR. REAL(zarg2)>hugelog) THEN
       zans=huge1;  info=20;   RETURN;   ENDIF
     zpart1=zsum1*(EXP(zlogbes1)+EXP(zarg1))
     zpart2=-2*zsum2*EXP(zarg2)
     zans=zunit*(zpart1+zpart2)/(1-EXP(-2*zunit*znupi))
   ENDIF
   info=MAX(info1,info2)
   error2=error2+abs2(znupi)
   abszpart1=abs2(zpart1)
   abszpart2=abs2(zpart2)
   error=(abszpart1*error1+abszpart2*error2)/(abszpart1+abszpart2)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(2)        !  With Taylor's expansion.   ABS(zeps) < 0.3
   IF(re_znu >= 0) THEN
      CALL neu_series(znu,xx,zans,info)
      RETURN
   ENDIF
! When re_znu<0, zneumann(znu,xx) is calculated using Eq. (41a) of Ref. (1).
   CALL bes_series(znua,xx,zsum1,zlogbes1,error,info1)
   IF(info1 > 17) THEN;   zans=huge1;    info=info1;     RETURN;   ENDIF
   IF(REAL(zlogbes1) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zbesa=zsum1*EXP(zlogbes1)
   CALL neu_series(znua,xx,zneua,info2)
   IF(info2 > 17) THEN;  zans=huge1;   info=info2;   RETURN;   ENDIF
   info=MAX(info1,info2)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=-zbesa*SIN(znupi)+zneua*COS(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(3)        !  With Debye's method         |znu/xx|<=1
! zneumann(znu,xx) is calculated using Eq. (41d) of Ref. (1).
   CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;     zans=huge1;     RETURN;   ENDIF
   IF(ABS(REAL(zlogmu)) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=-zunit*(zss1*EXP(zlogmu)-zss2*EXP(-zlogmu))/2
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info = 5
 CASE(4)        !  With Debye's method         |znu/xx|>1
   IF(re_znu >= 0) THEN
     IF(aim_znu >= 0) THEN
! When re_znu>=0 and aim_znu>=0, zneumann(znu,xx) is calculated 
! using Eq. (41e) of Ref. (1).
       CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
       IF(info > 17) THEN;  zans=huge1;   RETURN;  ENDIF
       IF(ABS(REAL(zlogmu)) > hugelog) THEN
         zans=huge1;  info=20;  RETURN;  ENDIF
       IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
       zans = zunit*(zss2*EXP(-zlogmu)-zss1*EXP(zlogmu)/2)
       IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
       IF(error > err_range1) info=5
       RETURN
     ENDIF
! When re_znu>=0 and aim_znu<0, zneumann(znu,xx) is calculated 
! using Eq. (41f) of Ref. (1).
     CALL bes_han_dby(zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(ABS(REAL(zlogmu,kp)) > hugelog) THEN
       zans=huge1;  info=20;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=-zunit*CONJG(zss2*EXP(-zlogmu)-zss1*EXP(zlogmu)/2)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zneumann(znu,xx) is calculated 
! using Eq. (41g) of Ref. (1).
     CALL bes_han_dby(-zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;   zans=huge1;     RETURN;     ENDIF
     IF(AIMAG(zconnupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
     IF(MAX(REAL(zlogmu-zconnupi*zunit), &
              REAL(-zlogmu)+ABS(AIMAG(znupi))) > hugelog) THEN
       zans=huge1;     info=20;     RETURN;    ENDIF
     error=error+abs2(zconnupi)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=zunit*(-EXP(zlogmu-zunit*zconnupi)*zss1/2 &
           +zss2*EXP(-zlogmu)*COS(zconnupi))
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zneumann(znu,xx) is calculated 
! using Eq. (41h) of Ref. (1).
   CALL bes_han_dby(znua,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;    zans=huge1;     RETURN;   ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(MAX(REAL(zlogmu-znupi*zunit),REAL(-zlogmu)+ABS(AIMAG(znupi))) &
       > hugelog) THEN;  zans=huge1;   info=20;   RETURN;   ENDIF
   error=error+abs2(znupi)
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zunit*(-EXP(zlogmu-zunit*znupi)*zss1/2+zss2*EXP(-zlogmu)*COS(znupi))
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(5)        !  Method by the Olver's asymptotic series
   IF(re_znu >= 0) THEN
! When re_znu>=0, zneumann(znu,xx) is calculated using Eq. (41e) of Ref. (1).
     CALL bes_olver(znu,xx,zbes,error1,info1)
     CALL han2_olver(znu,xx,zhan2,error2,info2)
     info=MAX(info1,info2)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     abes=abs2(zbes);  ahan2=abs2(zhan2)
     error=(abes*error1+ahan2*error2)/(abes+ahan2)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     zans=zunit*(zhan2-zbes)
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zneumann(znu,xx) is calculated 
! using Eq. (41g) of Ref. (1).
     CALL bes_olver(-zconnu,xx,zbesa,error1,info1)
     CALL han2_olver(-zconnu,xx,zhan2a,error2,info2)
     info=MAX(info1,info2)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     error=ABS(REAL(znupi,kp))+error2
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(ABS(AIMAG(zconnupi)) > hugelog) THEN;  info=30;  RETURN;  ENDIF
     zans=zunit*(-EXP(-zunit*zconnupi)*zbesa+COS(zconnupi)*zhan2a)
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zneumann(znu,xx) is calculated 
! using Eq. (41h) of Ref. (1).
   CALL bes_olver(znua,xx,zbesa,error1,info1)
   CALL han2_olver(znua,xx,zhan2a,error2,info2)
   info=MAX(info1,info2)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   error = ABS(REAL(znupi,kp))+error2
   IF(error > err_range3) THEN;  info=30;  RETURN;   ENDIF
   zans=zunit*(-EXP(-zunit*znupi)*zbesa+COS(znupi)*zhan2a)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(6)        !  By recurrence formula.
   CALL bes_recur(znu,xx,0,zbes,info1)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   IF(aim_znu >= 0) THEN
! When aim_znu>=0, zneumann(znu,xx) is calculated using Eq. (41b) of Ref. (1).
     CALL bes_recur(znua,xx,-1,zbesa,info2)
     info=MAX(info1,info2)
     z1=EXP(2*zunit*znupi)
     zans=zunit*(zbes*(z1+1)-2*zbesa)/(z1-1)
   ELSE
! When aim_znu<0, zneumann(znu,xx) is calculated using Eq. (41c) of Ref. (1).
     CALL bes_recur(znua,xx,1,zbesa,info2)
     info=MAX(info1,info2)
     z1=EXP(-2*zunit*znupi)
     zans=zunit*(zbes*(z1+1)-2*zbesa)/(1-z1)
   ENDIF
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(7)        !  By recurrence formula.   ABS(AIMAG(znu)) < 0.5
   IF(re_znu >= 0) THEN
! When re_znu>=0, zneumann(znu,xx) is calculated using Eq. (41e) of Ref. (1).
     CALL bes_recur(znu,xx,0,zbes,info1)
     CALL han2_temme(znu,xx,zhan2,info2)
     info=MAX(info1,info2)
     zans=zunit*(zhan2-zbes)
     RETURN
   ENDIF
! When re_znu<0, zneumann(znu,xx) is calculated using Eq. (41h) of Ref. (1).
   CALL bes_recur(-zconnu,xx,1,zbesa,info1)
   CALL han2_temme(-zconnu,xx,zhan2a,info2)
   info=MAX(info1,info2)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zunit*(-zbesa+COS(zconnupi)*zhan2a)
   zans=CONJG(zans)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 END SELECT
 END SUBROUTINE neumann

 SUBROUTINE hankel1(znu,xx,zans,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zans
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zans.  zans=zhankel1(znu,xx).
! This is invoked only by a user's program.
! The dummy arguments znu, xx, zans, info and the outline of the calculation
! in SUBROUTINE hankel1 are explained in the beginning of this file.
 COMPLEX(kp):: zhan2,zconnu
 zconnu=CONJG(znu)
 CALL hankel2(zconnu,xx,zhan2,info)
 zans=CONJG(zhan2)
 END SUBROUTINE hankel1

 SUBROUTINE hankel2(znu,xx,zans,info)
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zans
 INTEGER,INTENT(OUT):: info
! This SUBROUTINE calculates zans.  zans=zhankel2(znu,xx).
! This is invoked by a user's program and SUBROUTINE hankel1.
! The dummy arguments znu, xx, zans, info and the outline of the calculation
! in SUBROUTINE hankel2 are explained in the beginning of this file.
 COMPLEX(kp):: zarg1,zarg2,zbes,zbes1,zbes2,zlogbes1,zlogbes2, &
  zlogmu,zneu,znua,zconnu,znupi,zpart1,zpart2,zss1,zss2,zsum1,zsum2,z1,z2
 REAL(kp):: aim_znu,re_znu,abszpart1,abszpart2
 REAL(kp):: error,error1,error2
 INTEGER:: info1,info2,nregion
! Determination of nregion
 nregion=num_region(znu,xx)
 zans=0
 IF(nregion == 0) THEN;   info=30;   RETURN;  ENDIF
 re_znu=REAL(znu,kp);  aim_znu=AIMAG(znu)
 znupi=znu*pi;         zconnu=CONJG(znu)
 znua=-znu
 SELECT CASE(nregion)
 CASE(1)
   CALL bes_series(znu,xx,zsum1,zlogbes1,error1,info1)
   IF(info1 > 17) THEN;   zans=huge1;   info=info1;    RETURN;   ENDIF
   IF(error1 > err_range3) THEN;   info=30;     RETURN;    ENDIF
   CALL bes_series(-znu,xx,zsum2,zlogbes2,error2,info2)
   info=MAX(info1,info2)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(error2 > err_range3) THEN;  info=30;  RETURN;  ENDIF
   IF(aim_znu >= 0) THEN
! When aim_znu>0, zhankel2(znu,xx) is calculated using Eq. (42a) of Ref. (1).
     zarg1=zlogbes1+2*zunit*znupi
     zarg2=zlogbes2+zunit*znupi
     IF(REAL(zarg1)>hugelog .OR. REAL(zarg2)>hugelog) THEN
       zans=huge1;     info=20;     RETURN
     ENDIF
     zpart1=zsum1*EXP(zarg1);  zpart2=-zsum2*EXP(zarg2)
     zans=2*(zpart1+zpart2)/(EXP(2*zunit*znupi)-1)
     error1=error1+2*abs2(znupi);  error2=error2+abs2(znupi)
     abszpart1=abs2(zpart1);  abszpart2=abs2(zpart2)
     IF(abszpart1+abszpart2 < tiny1) THEN
       info=20;  RETURN
     ELSE
       error=(abszpart1*error1+abszpart2*error2)/(abszpart1+abszpart2)
     ENDIF
   ELSE
! When aim_znu<=0, zhankel2(znu,xx) is calculated using Eq. (42b) of Ref. (1).
     zarg1=zlogbes1
     zarg2=zlogbes2-zunit*znupi
     IF(REAL(zarg1)>hugelog .OR. REAL(zarg2)>hugelog) THEN
       zans=huge1;     info=20;     RETURN
     ENDIF
     zpart1=zsum1*EXP(zarg1);    zpart2=-zsum2*EXP(zarg2)
     zans=2*(zpart1+zpart2)/(1-EXP(-2*zunit*znupi))
     error2=abs2(znupi)+error2
     abszpart1=abs2(zpart1);     abszpart2=abs2(zpart2)
     error=(abszpart1*error1+abszpart2*error2)/(abszpart1+abszpart2)
   ENDIF
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(2)
   IF(re_znu >= 0) THEN
! When re_znu>=0, zhankel2(znu,xx) is calculated using Eq. (42c) of Ref. (1).
     CALL bes_series(znu,xx,zsum1,zlogbes1,error,info1)
     IF(info1 > 17) THEN;     zans=huge1;     info=info1;     RETURN;    ENDIF
     IF(REAL(zlogbes1) > hugelog) THEN;  zans=huge1;  info=20;   RETURN;   ENDIF
     IF(error > err_range3) THEN;  info=30;   info=info1;   RETURN;   ENDIF
     zbes=zsum1*EXP(zlogbes1)
     IF(error > err_range1) info1=5
     IF(error > err_range2) info1=10
     CALL neu_series(znu,xx,zneu,info2)
     info=MAX(info1,info2)
     zans=zbes-zunit*zneu
     RETURN
   ENDIF
! When re_znu<0, zhankel2(znu,xx) is calculated using Eq. (42d) of Ref. (1).
   CALL bes_series(znua,xx,zsum1,zlogbes1,error1,info1)
   IF(info1 > 17) THEN;  zans=huge1;   info=info1;   RETURN;   ENDIF
   zarg1=zlogbes1+zunit*znupi
   IF(REAL(zarg1) > hugelog) THEN;   zans=huge1;   info=20;    RETURN;   ENDIF
   error=abs2(znupi)+error1
   IF(error > err_range3) THEN;    info=30;     RETURN;   ENDIF
   CALL neu_series(znua,xx,zneu,info2)
   info=MAX(info1,info2)
   zans=EXP(zarg1)*zsum1-EXP(zunit*znupi)*zunit*zneu
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info = 5
 CASE(3)
   CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;     zans=huge1;     RETURN;   ENDIF
   IF(REAL(-zlogmu) > hugelog) THEN;   zans=huge1;   info=20;   RETURN;   ENDIF
   IF(error > err_range3) THEN;   info=30;   RETURN;   ENDIF
   zans=zss2*EXP(-zlogmu)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(4)
   IF(re_znu >= 0) THEN
     IF(aim_znu >= 0) THEN
       CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
       IF(info > 17) THEN;   zans=huge1;     RETURN;    ENDIF
       IF(REAL(-zlogmu) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
       IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
       zans=zss2*EXP(-zlogmu)
       IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
       IF(error > err_range1) info=5
       RETURN
     ENDIF
! When re_znu>=0 and aim_znu<0, zhankel2(znu,xx) is calculated
! using Eq. (42f) of Ref. (1).
     CALL bes_han_dby(zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;   zans=huge1;     RETURN;     ENDIF
     IF(ABS(REAL(zlogmu)) > hugelog) THEN
       zans=huge1;  info=20;   RETURN;   ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=CONJG(zss1*EXP(zlogmu)-zss2*EXP(-zlogmu))
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zhankel2(znu,xx) is calculated 
! using Eq. (42g) of Ref. (1).
     CALL bes_han_dby(-zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;    zans=huge1;     RETURN;     ENDIF
     z1=zlogmu+CONJG(zunit*znupi)
     z2=-zlogmu+CONJG(zunit*znupi)
     IF(MAX(REAL(z1),REAL(z2)) > hugelog) THEN
       zans=huge1;  info=20;  RETURN;  ENDIF
     error=error+abs2(znupi)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=CONJG(zss1*EXP(z1)-zss2*EXP(z2))
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zhankel2(znu,xx) is calculated
! using Eq. (42e) of Ref. (1).
   CALL bes_han_dby(znua,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;   zans=huge1;    RETURN;   ENDIF
   z1=zunit*znupi-zlogmu
   IF(REAL(z1) > hugelog) THEN;  zans=huge1;  info=20;    RETURN;   ENDIF
   error=error+abs2(znupi)
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zss2*EXP(z1)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(5)
   IF(re_znu >= 0) THEN
     CALL han2_olver(znu,xx,zans,error,info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0, zhankel2(znu,xx) is calculated using Eq. (42e) of Ref. (1).
   CALL han2_olver(znua,xx,zans,error,info)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   error=ABS(REAL(znupi,kp))+error
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   IF(REAL(zunit*znupi) > hugelog) THEN;  info=30;  RETURN;  ENDIF
   zans=zans*EXP(zunit*znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(6)
   IF(aim_znu > 0) THEN
! zhankel2(znu,xx) is calculated using Eq. (42a) of Ref. (1).
     CALL bes_recur( znu,xx, 1,zbes1,info1)
     CALL bes_recur(-znu,xx,-1,zbes2,info2)
     info=MAX(info1,info2)
     z1=EXP(zunit*znupi)
     zans=2*(zbes1*z1-zbes2)/(z1*z1-1)
   ELSE
! zhankel2(znu,xx) is calculated using Eq. (42b) of Ref. (1).
     CALL bes_recur( znu,xx,0,zbes1,info1)
     CALL bes_recur(-znu,xx,1,zbes2,info2)
     info=MAX(info1,info2)
     zans=2*(zbes1-zbes2)/(1-EXP(-2*zunit*znupi))
   ENDIF
 CASE(7)
   IF(re_znu >= 0) THEN
     CALL han2_temme(znu,xx,zans,info)
     RETURN
   ENDIF
! When re_znu<0, zhankel2(znu,xx) is calculated using Eq. (42e) of Ref. (1).
   CALL han2_temme(znua,xx,zans,info)
   zans=zans*EXP(zunit*znupi)
 END SELECT
 END SUBROUTINE hankel2

 END MODULE mod_bes

