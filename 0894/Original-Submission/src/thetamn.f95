! This module provides the bounding parameters ($\theta_m^n$)
! for the (modified) scaling and squaring method.
!
! The parameters are defined for each of the four precisions,
! sp, wp, tp, and qp, respectively. It also depends on the maximum index
! of phi-functions to be computed simultaneously.
!
! References:
!
!   Nicholas J. Higham,
!   The scaling and squaring method for the matrix exponential revisited,
!   SIAM Journal on Matrix Analysis and Applications,
!   Volume 26, Number 4, pp. 1179-1193, 2005.
!
!   S. Koikari,
!   An error analysis of the modified scaling and squaring method,
!   Computers and Mathematics with Applications,
!   Volume 53, pp. 1293-1305, 2007.
!
! The module is intended for internal use only.

module thetamn

  use floattypes           ! defines ``sp'', ``wp'', ``tp'', and ``qp''.

  implicit none

  public

contains


!    Name : sp_getthetamn
! Purpose : This function provides the bounding parameters for the (modified)
!         : scaling and squaring method implemented in single precision.
!   Input : ``nmax'' is the maximum index of phi-functions to be computed.
!         : (0 <= nmax <= 5).
!  Output : The return value stores the bounding parameters for [3/3], [5/5],
!         : and [7/7] rational approximations in 3rd, 5th and 7th entries,
!         : respectively.

pure function sp_getthetamn(nmax)          ! for the single precision
  integer,       intent(in)      :: nmax
  real(kind=sp), dimension(3:17) :: sp_getthetamn, ctam

  ctam = 0.0_sp

  select case (nmax)
    case (0)                               ! for phi_{0} = exponential only
      ctam( 3) = 4.2587300e-1_sp
      ctam( 5) = 1.8801526e+0_sp
      ctam( 7) = 3.9257248e+0_sp

    case (1)                               ! for computing phi_{0,1}
      ctam( 3) = 3.0659951e-1_sp
      ctam( 5) = 1.7664200e+0_sp
      ctam( 7) = 3.8470623e+0_sp

    case (2)                               ! for computing phi_{0,1,2}
      ctam( 3) = 1.7404177e-1_sp
      ctam( 5) = 1.6087347e+0_sp
      ctam( 7) = 3.7434055e+0_sp

    case (3)                               ! for computing phi_{0,1,2,3}
      ctam( 3) = 6.2243331e-2_sp
      ctam( 5) = 1.4030294e+0_sp
      ctam( 7) = 3.6074010e+0_sp

    case (4)                               ! for computing phi_{0,1,2,3,4}
      ctam( 3) = 7.0592330e-3_sp
      ctam( 5) = 1.1489968e+0_sp
      ctam( 7) = 3.4322799e+0_sp

    case (5)                               ! for computing phi_{0,1,2,3,4,5}
      ctam( 3) = 0.0000000e+0_sp
      ctam( 5) = 8.5225589e-1_sp
      ctam( 7) = 3.2118216e+0_sp

    case default

  end select

  sp_getthetamn = ctam
end function sp_getthetamn


!    Name : wp_getthetamn
! Purpose : This function provides the bounding parameters for the (modified)
!         : scaling and squaring method implemented in double precision.
!   Input : ``nmax'' is the maximum index of phi-functions to be computed.
!         : (0 <= nmax <= 5).
!  Output : The return value stores the bounding parameters for [3/3], [5/5],
!         : [7/7], [9/9], and [13/13] rational approximations in 3rd, 5th, 7th,
!         : 9th and 13th entries, respectively.

pure function wp_getthetamn(nmax)           ! for the double precision
  integer,       intent(in)      :: nmax
  real(kind=wp), dimension(3:17) :: wp_getthetamn, ctam

  ctam = 0.0_wp

  select case (nmax)

    case (0)                               ! for phi_{0} = exponential only
      ctam( 3) = 1.4955852179582915e-2_wp
      ctam( 5) = 2.5393983300632320e-1_wp
      ctam( 7) = 9.5041789961629318e-1_wp
      ctam( 9) = 2.0978479612570674e+0_wp
      ctam(13) = 5.3719203511481522e+0_wp

    case (1)                               ! for computing phi_{0,1}
      ctam( 3) = 5.6157211029084633e-3_wp
      ctam( 5) = 2.0044664339611825e-1_wp
      ctam( 7) = 8.7908729387086115e-1_wp
      ctam( 9) = 2.0342766605280691e+0_wp
      ctam(13) = 5.3302347837622749e+0_wp

    case (2)                               ! for computing phi_{0,1,2}
      ctam( 3) = 1.1684669784240026e-3_wp
      ctam( 5) = 1.4284050683124090e-1_wp
      ctam( 7) = 7.8683344602375597e-1_wp
      ctam( 9) = 1.9483995540541480e+0_wp
      ctam(13) = 5.2785617757601364e+0_wp

    case (3)                               ! for computing phi_{0,1,2,3}
      ctam( 3) = 7.7543685445764846e-5_wp
      ctam( 5) = 8.8965779243155152e-2_wp
      ctam( 7) = 6.7768547444289960e-1_wp
      ctam( 9) = 1.8389562774597516e+0_wp
      ctam(13) = 5.2137825437140124e+0_wp

    case (4)                               ! for computing phi_{0,1,2,3,4}
      ctam( 3) = 3.0538289913172909e-7_wp
      ctam( 5) = 4.5631120461304921e-2_wp
      ctam( 7) = 5.5679491195623167e-1_wp
      ctam( 9) = 1.7061659666051970e+0_wp
      ctam(13) = 5.1328733502655685e+0_wp

    case (5)                               ! for computing phi_{0,1,2,3,4,5}
      ctam( 3) = 0.0000000000000000e+0_wp
      ctam( 5) = 1.7260590043359373e-2_wp
      ctam( 7) = 4.3055278727389165e-1_wp
      ctam( 9) = 1.5512477914306754e+0_wp
      ctam(13) = 5.0332008644524018e+0_wp

    case default

  end select

  wp_getthetamn = ctam
end function wp_getthetamn

#ifdef __USE_TPREC


!    Name : tp_getthetamn
! Purpose : This function provides the bounding parameters for the (modified)
!         : scaling and squaring method implemented in extended double
!         : precision.
!   Input : ``nmax'' is the maximum index of phi-functions to be computed.
!         : (0 <= nmax <= 5).
!  Output : The return value stores the bounding parameters for [3/3], [5/5],
!         : [7/7], [9/9], and [13/13] rational approximations in 3rd, 5th, 7th,
!         : 9th and 13th entries, respectively.

pure function tp_getthetamn(nmax)          ! for the extended double precision
  integer,       intent(in)      :: nmax
  real(kind=tp), dimension(3:17) :: tp_getthetamn, ctam

  ctam = 0.0_tp

  select case (nmax)

    case (0)                               ! for phi_{0} = exponential only
      ctam( 3) = 4.19684972322669896709e-3_tp
      ctam( 5) = 1.18481167346938230910e-1_tp
      ctam( 7) = 5.51703884806867002738e-1_tp
      ctam( 9) = 1.37598688755878453832e+0_tp
      ctam(13) = 4.02460989066973530633e+0_tp

    case (1)                               ! for computing phi_{0,1}
      ctam( 3) = 1.22255037974879927431e-3_tp
      ctam( 5) = 8.62786462020972603360e-2_tp
      ctam( 7) = 4.93672731345846880760e-1_tp
      ctam( 9) = 1.31542417187649897139e+0_tp
      ctam(13) = 3.98153129997292527556e+0_tp

    case (2)                               ! for computing phi_{0,1,2}
      ctam( 3) = 1.73715255454492039354e-4_tp
      ctam( 5) = 5.53742508758197200192e-2_tp
      ctam( 7) = 4.23052758487723909585e-1_tp
      ctam( 9) = 1.23511197597101207261e+0_tp
      ctam(13) = 3.92584238815642028691e+0_tp

    case (3)                               ! for computing phi_{0,1,2,3}
      ctam( 3) = 6.10626984155541232635e-6_tp
      ctam( 5) = 3.00863262703780426645e-2_tp
      ctam( 7) = 3.45012067470042379961e-1_tp
      ctam( 9) = 1.13651008803443930284e+0_tp
      ctam(13) = 3.85451679027754075822e+0_tp

    case (4)                               ! for computing phi_{0,1,2,3,4}
      ctam( 3) = 6.74807313531318146947e-9_tp
      ctam( 5) = 1.28515538272720943462e-2_tp
      ctam( 7) = 2.64864046775015690720e-1_tp
      ctam( 9) = 1.02199356289452103113e+0_tp
      ctam(13) = 3.76527183595751788698e+0_tp

    case (5)                               ! for computing phi_{0,1,2,3,4,5}
      ctam( 3) = 0.00000000000000000000e+0_tp
      ctam( 5) = 3.76380643481241062467e-3_tp
      ctam( 7) = 1.88141794899770755239e-1_tp
      ctam( 9) = 8.94585265271049610030e-1_tp
      ctam(13) = 3.65657885468975940316e+0_tp

    case default

  end select

  tp_getthetamn = ctam
end function tp_getthetamn

#endif

#ifdef __USE_QPREC


!    Name : qp_getthetamn
! Purpose : This function provides the bounding parameters for the (modified)
!         : scaling and squaring method implemented in quadruple precision.
!   Input : ``nmax'' is the maximum index of phi-functions to be computed.
!         : (0 <= nmax <= 5).
!  Output : The return value stores the bounding parameters for [3/3], [5/5],
!         : [7/7], [9/9], [13/13], and [17/17] rational approximations in
!         : 3rd, 5th, 7th, 9th, 13th, and 17th entries, respectively.

pure function qp_getthetamn(nmax)          ! for the quadruple precision
  integer,       intent(in)      :: nmax
  real(kind=qp), dimension(3:17) :: qp_getthetamn, ctam

  ctam = 0.0_qp

  select case (nmax)
    case (0)                               ! for phi_{0} = exponential only
      ctam( 3) = 1.46053455683310364021066602616115585e-5_qp
      ctam( 5) = 3.96841112814126231019629191039739660e-3_qp
      ctam( 7) = 4.87820434671134757850752555944929064e-2_qp
      ctam( 9) = 2.08803184913573258944286650890651851e-1_qp
      ctam(13) = 1.09577903412722875087933655650390340e+0_qp
      ctam(17) = 2.82622562499379681109119039217579934e+0_qp

    case (1)                               ! for computing phi_{0,1}
      ctam( 3) = 1.37153901327002586257660033527961485e-6_qp
      ctam( 5) = 1.98754311927418042499499932916559914e-3_qp
      ctam( 7) = 3.66264414707061879963862418016414031e-2_qp
      ctam( 9) = 1.82175335104961411724033224307811350e-1_qp
      ctam(13) = 1.05604148374347194047372737194498234e+0_qp
      ctam(17) = 2.79248377638733660534409485893894768e+0_qp

    case (2)                               ! for computing phi_{0,1,2}
      ctam( 3) = 3.56639909546935339860057494575488737e-8_qp
      ctam( 5) = 7.96123287772769672184222929933457400e-4_qp
      ctam( 7) = 2.53707476486966972596906649464272623e-2_qp
      ctam( 9) = 1.52763533866349389909316148481688068e-1_qp
      ctam(13) = 1.00463369396894043809287505954683945e+0_qp
      ctam(17) = 2.74778490088652914956009444801471434e+0_qp

    case (3)                               ! for computing phi_{0,1,2,3}
      ctam( 3) = 0.00000000000000000000000000000000000e+0_qp
      ctam( 5) = 2.35651704789423246593966002923365749e-4_qp
      ctam( 7) = 1.60199477640190548783749722742119286e-2_qp
      ctam( 9) = 1.22900937275393023012949695931443406e-1_qp
      ctam(13) = 9.43202873848938421452521231182399655e-1_qp
      ctam(17) = 2.69092609666437582773060375114078627e+0_qp

    case (4)                               ! for computing phi_{0,1,2,3,4}
      ctam( 3) = 0.00000000000000000000000000000000000e+0_qp
      ctam( 5) = 4.47880432802818328824671226373844945e-5_qp
      ctam( 7) = 9.02373526524527637245354905291214757e-3_qp
      ctam( 9) = 9.43908717588872462209889770342043451e-2_qp
      ctam(13) = 8.73507175928299530556273027518788336e-1_qp
      ctam(17) = 2.62153977168420946059086980261632697e+0_qp

    case (5)                               ! for computing phi_{0,1,2,3,4,5}
      ctam( 3) = 0.00000000000000000000000000000000000e+0_qp
      ctam( 5) = 4.22441512002681475740886617380702540e-6_qp
      ctam( 7) = 4.38414573611837733947981661756452212e-3_qp
      ctam( 9) = 6.86568266566271739806608205871507068e-2_qp
      ctam(13) = 7.97297541003560768542331159245196618e-1_qp
      ctam(17) = 2.53978320334260924747081303674618091e+0_qp

    case default

  end select

  qp_getthetamn = ctam
end function qp_getthetamn

#endif

end module thetamn

