! This module provides an interface to the BLAS library.
!
! This module is useful only with a FORTRAN compiler that does not invoke
! ``xGEMM'' (x=S,D,C,Z) or its equivalence for the intrinsic function,
! ``matmul''.
!
! This is the case for, e.g., Intel's ``ifort'' command. To enable explicit
! use of BLAS, please add a phrase, ``-Dpure= -D__USE_BLAS'', literally to
! the command line arguments of your compiler. Although the pure attribute
! of all functions and subroutines in this library are lost by this option,
! efficiency improves significantly.
!
! The above option is not necessary for, e.g., Sun's ``f95'' command bundled
! in Sun Studio 11, because ``matmul'' is automatically linked to a GEMM
! equivalence in Sun Performance Library.
!
! The module is intended for internal use only.

module blasinterface

  use floattypes

  implicit none

  integer, parameter :: blasinterface_dummy = 0

#ifdef __USE_BLAS

  ! When ``__USE_BLAS'' is not defined, this module disappears.

  external sgemm, dgemm, cgemm, zgemm     ! are level-3 BLAS routines.

  !    Name : ggemm  ( a generic name )
  !   Usage : call ggemm(A,B,C,alpha,beta)
  ! Purpose : This routine computes ``alpha AB + beta C,''
  !         : where A, B, and C are rectangular matrices,
  !         : and ``alpha'' and ``beta'' are INTEGER scalars.
  !   Input : A, B, and C are rectangular matrices of consistent size.
  !         :``alpha'' and ``beta'' are integer scalars.
  !  Output : C is overwritten by ``alpha AB + beta C''.

  interface ggemm
    module procedure rs_ggemm             ! For Real    in Single    precision
    module procedure rw_ggemm             ! For Real    in Double    precision
    module procedure cs_ggemm             ! For Complex in Single    precision
    module procedure cw_ggemm             ! For Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rt_ggemm             ! For Real    in Triple    precision
    module procedure ct_ggemm             ! For Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rq_ggemm             ! For Real    in Quadruple precision
    module procedure cq_ggemm             ! For Complex in Quadruple precision
#endif
  end interface ggemm

  private
  public  :: ggemm

contains


!    Name : rs_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are real matrices in single precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are real matrices of consistent size.
!         :``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.
!    Note : The BLAS routine ``SGEMM'' is required.

pure subroutine rs_ggemm(a,b,c,alpha,beta)
  integer,       intent(in   ) :: alpha, beta
  real(kind=sp), intent(in   ) :: a(1:,1:), b(1:,1:)
  real(kind=sp), intent(inout) :: c(1:,1:)
  integer :: m, n, k

  if (beta == 0) c = 0.0_sp
  m = size(a,1); n = size(b,2); k = size(a,2)
  call sgemm('N','N',m,n,k,real(alpha,sp),a,m,b,k,real(beta,sp),c,m)
end subroutine rs_ggemm


!    Name : rw_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are real matrices in double precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are real matrices of consistent size.
!         :``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.
!    Note : The BLAS routine ``DGEMM'' is required.

pure subroutine rw_ggemm(a,b,c,alpha,beta)
  integer,       intent(in   ) :: alpha, beta
  real(kind=wp), intent(in   ) :: a(1:,1:), b(1:,1:)
  real(kind=wp), intent(inout) :: c(1:,1:)
  integer :: m, n, k

  if (beta == 0) c = 0.0_wp
  m = size(a,1); n = size(b,2); k = size(a,2)
  call dgemm('N','N',m,n,k,real(alpha,wp),a,m,b,k,real(beta,wp),c,m)
end subroutine rw_ggemm

#ifdef __USE_TPREC


!    Name : rt_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are real matrices in extended double precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are real matrices of consistent size.
!         :``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.

pure subroutine rt_ggemm(a,b,c,alpha,beta)
  integer,       intent(in   ) :: alpha, beta
  real(kind=tp), intent(in   ) :: a(1:,1:), b(1:,1:)
  real(kind=tp), intent(inout) :: c(1:,1:)

  if (beta == 0) c = 0.0_tp
  c = matmul(a,b) * real(alpha,tp) + c * real(beta,tp)
end subroutine rt_ggemm

#endif

#ifdef __USE_QPREC


!    Name : rq_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are real matrices in quadruple precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are real matrices of consistent size.
!         :``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.

pure subroutine rq_ggemm(a,b,c,alpha,beta)
  integer,       intent(in   ) :: alpha, beta
  real(kind=qp), intent(in   ) :: a(1:,1:), b(1:,1:)
  real(kind=qp), intent(inout) :: c(1:,1:)

  if (beta == 0) c = 0.0_qp
  c = matmul(a,b) * real(alpha,qp) + c * real(beta,qp)
end subroutine rq_ggemm

#endif


!    Name : cs_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are complex matrices in single precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are complex matrices of consistent size.
!         : ``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.
!    Note : The BLAS routine ``CGEMM'' is required.

pure subroutine cs_ggemm(a,b,c,alpha,beta)
  integer,          intent(in   ) :: alpha, beta
  complex(kind=sp), intent(in   ) :: a(1:,1:), b(1:,1:)
  complex(kind=sp), intent(inout) :: c(1:,1:)
  integer :: m, n, k

  if (beta == 0) c = cmplx(0.0_sp,0.0_sp,sp)
  m = size(a,1); n = size(b,2); k = size(a,2)
  call cgemm('N','N',m,n,k,cmplx(alpha,0,sp),a,m,b,k,cmplx(beta,0,sp),c,m)
end subroutine cs_ggemm


!    Name : cw_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are complex matrices in double precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are complex matrices of consistent size.
!         : ``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.
!    Note : The BLAS routine ``ZGEMM'' is required.

pure subroutine cw_ggemm(a,b,c,alpha,beta)
  integer,          intent(in   ) :: alpha, beta
  complex(kind=wp), intent(in   ) :: a(1:,1:), b(1:,1:)
  complex(kind=wp), intent(inout) :: c(1:,1:)
  integer :: m, n, k

  if (beta == 0) c = cmplx(0.0_wp,0.0_wp,wp)
  m = size(a,1); n = size(b,2); k = size(a,2)
  call zgemm('N','N',m,n,k,cmplx(alpha,0,wp),a,m,b,k,cmplx(beta,0,wp),c,m)
end subroutine cw_ggemm

#ifdef __USE_TPREC


!    Name : ct_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are complex matrices in extended double precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are complex matrices of consistent size.
!         : ``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.

pure subroutine ct_ggemm(a,b,c,alpha,beta)
  integer,          intent(in   ) :: alpha, beta
  complex(kind=tp), intent(in   ) :: a(1:,1:), b(1:,1:)
  complex(kind=tp), intent(inout) :: c(1:,1:)

  if (beta == 0) c = cmplx(0.0_tp,0.0_tp,tp)
  c = matmul(a,b) * cmplx(alpha,0,tp) + c * cmplx(beta,0,tp)
end subroutine ct_ggemm

#endif

#ifdef __USE_QPREC


!    Name : cq_ggemm
!   Usage : This subroutine is invoked through the generic name ``ggemm''.
! Purpose : This routine computes ``alpha AB + beta C'',
!         : where A, B, and C are complex matrices in quadruple precision,
!         : and ``alpha'' and ``beta'' are INTEGER scalars.
!   Input : A, B, and C are complex matrices of consistent size.
!         : ``alpha'' and ``beta'' are integer scalars.
!  Output : C is overwritten by ``alpha AB + beta C''.

pure subroutine cq_ggemm(a,b,c,alpha,beta)
  integer,          intent(in   ) :: alpha, beta
  complex(kind=qp), intent(in   ) :: a(1:,1:), b(1:,1:)
  complex(kind=qp), intent(inout) :: c(1:,1:)

  if (beta == 0) c = cmplx(0.0_qp,0.0_qp,qp)
  c = matmul(a,b) * cmplx(alpha,0,qp) + c * cmplx(beta,0,qp)
end subroutine cq_ggemm

#endif

#endif
! the end of "#ifdef __USE_BLAS" clause

end module blasinterface

