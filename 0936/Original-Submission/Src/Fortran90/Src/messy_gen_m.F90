! Time-stamp:  <2016-08-11 15:03:40 m>
!>> 2015-12-06 Package of messy as suggested by a referee

! This code is free for noncommercial use as long as the following source files
! are included:  messy_m.f90, precision_m.f90, tmessy.f90, sample_m.f90,
! Thrdtmessy.F90, Makefile, messy_doc.tex, and messy_doc.bib.

! Used for printing error messages and for pretty output of other types.
! For documentation of use see the files mentioned above.

! Basic use: (Use _d for double, _s for single, or _q for quad precision.)

! use messy_gen, only: messy_d=>messy, rk_d=>rk ! Change _d for diff. precisions
! type(messy_ty) :: e
! Set data
! call messy(e, "Text, see first big block of comments in messy_bod.F90, ...)
module  messy_m
  use  messy_s_m, only: messy_s => messy, messy_ty_s => messy_ty, rk_s => rk
  use  messy_d_m, only: messy_d => messy, messy_ty_d => messy_ty, rk_d => rk
#ifdef NOQ
  interface  messy
    module procedure messy_s, messy_d
#else
  use  messy_q_m, only: messy_q => messy, messy_ty_q => messy_ty, rk_q => rk
  interface  messy
    module procedure messy_s, messy_d, messy_q
#endif
  end interface messy
end module  messy_m
