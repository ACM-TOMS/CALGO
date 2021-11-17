! This files illustrates how one might setup a generic Fortran Library routine
! Basic use: (Use _d for double, _s for single, or _q for quad precision.)

! use sample_gen, only: sample_d=>sample, rk_d=>rk ! Change _d for diff. precisions
! type(sample_ty) :: s
! Set data
! call sample(s) ! You would probably have more arguments here.
module  sample_m
  use  sample_s_m, only: messy_ty_s => messy_ty, sample_s => sample,&
    & sample_ty_s => sample_ty, fp_s => fp
  use  sample_d_m, only: messy_ty_d => messy_ty, sample_d => sample,&
    & sample_ty_d => sample_ty, fp_d => fp
#ifdef NOQ
  interface  sample
    module procedure sample_s, sample_d
#else  
  use  sample_q_m, only: messy_ty_q => messy_ty, sample_q => sample,&
    & sample_ty_q => sample_ty, fp_q => fp
  interface  sample
    module procedure sample_s, sample_d, sample_q
#endif    
  end interface sample
end module  sample_m
