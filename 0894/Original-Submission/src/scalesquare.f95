! This module computes the matrix values of phi-functions by
! the modified scaling and squaring method (S & S).
!
! Six generic functions and subroutines are provided:
!
!``sasmtrphif'': function  : S & S for a diagonal or a square matrix,
!``sasqtrphif'': function  : S & S for an upper quasi-triangular matrix,
!``sasmtrphi'' : subroutine: a subroutine version of ``sasmtrphif'',
!``sasqtrphi'' : subroutine: a subroutine version of ``sasqtrphif'',
!``sqrmtrphif'': function  : squaring only for a diagonal or a square matrix,
!``sqrqtrphif'': function  : squaring only for an upper quasi-triangular matrix.
!
! The details of usage are described in the module header.
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

module scalesquare

  use floattypes
  use thetamn
  use matrixpwrtag
  use rationalphi
  use mtrcfgphilog

  implicit none


  !    Name : sasmtrphif ( a generic name )
  !         :
  !   Usage : phin = sasmtrphif(upto,scale,A)
  !         :
  ! Purpose : This function computes the matrix values of phi-functions for a
  !         : diagonal or a square matrix by the modified scaling and squaring
  !         : method.
  !         :
  !   Input : - ``upto'' is an integer, and specifies the maximum index of
  !         :   phi-functions to be computed simultaneously (0 <= upto <= 5).
  !         :
  !         : - ``A'' is a real or complex array whose shape is 
  !         :   either A(1:N) or A(1:N,1:N).
  !         :
  !         : - ``scale'' is a real or complex scalar that should multiply the
  !         :   entries of ``A''. The type of ``scale'' must be the same as
  !         :   that of the entries of ``A''.
  !         :
  !  Output : - When the shape of ``A'' is (1:N), it is regarded as a diagonal
  !         :   matrix, and the output ``phin'' has the shape (1:N,0:upto).
  !         :   The subarray, phin(1:N,k), stores the value of phi_k(scale*A)
  !         :   for 0 <= k <= upto.
  !         :
  !         : - When the shape of ``A'' is (1:N,1:N), it is regarded as a square
  !         :   matrix, and the output ``phin'' has the shape (1:N,1:N,0:upto).
  !         :   The subarray, phin(1:N,1:N,k), is the value of phi_k(scale*A)
  !         :   for 0 <= k <= upto.

  interface sasmtrphif
    module procedure rsdg_scsqrf     ! Real    Diagonal in Single    precision
    module procedure csdg_scsqrf     ! Complex Diagonal in Single    precision
    module procedure rwdg_scsqrf     ! Real    Diagonal in Double    precision
    module procedure cwdg_scsqrf     ! Complex Diagonal in Double    precision
#ifdef __USE_TPREC
    module procedure rtdg_scsqrf     ! Real    Diagonal in Triple    precision
    module procedure ctdg_scsqrf     ! Complex Diagonal in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqdg_scsqrf     ! Real    Diagonal in Quadruple precision
    module procedure cqdg_scsqrf     ! Complex Diagonal in Quadruple precision
#endif

    module procedure rssq_scsqrf     ! Real    Square   in Single    precision
    module procedure cssq_scsqrf     ! Complex Square   in Single    precision
    module procedure rwsq_scsqrf     ! Real    Square   in Double    precision
    module procedure cwsq_scsqrf     ! Complex Square   in Double    precision
#ifdef __USE_TPREC
    module procedure rtsq_scsqrf     ! Real    Square   in Triple    precision
    module procedure ctsq_scsqrf     ! Complex Square   in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqsq_scsqrf     ! Real    Square   in Quadruple precision
    module procedure cqsq_scsqrf     ! Complex Square   in Quadruple precision
#endif
  end interface sasmtrphif


  !    Name : sasqtrphif ( a generic name )
  !         :
  !   Usage : phin = sasqtrphif(upto,scale,a)
  !         :
  ! Purpose : This function computes the matrix values of phi-functions
  !         : for an upper quasi-triangular matrix by the modified scaling
  !         : and squaring method.
  !         :
  !   Input : - ``upto'' is an integer, and specifies the maximum index of
  !         :   phi-functions to be computed simultaneously (0 <= upto <= 5).
  !         :
  !         : - ``A'' is a real or complex array whose shape is (1:N,1:N).
  !         :   This array stores an upper quasi-triangular matrix.
  !         :
  !         : - ``scale'' is a real or complex scalar that should multiply the
  !         :   entries of ``A''. The type of ``scale'' must be the same as
  !         :   that of the entries of ``A''.
  !         :
  !  Output : The subarray, phin(1:N,1:N,k), stores the value of
  !         : phi_k(scale*A(1:N,1:N)) for 0 <= k <= upto.
  !         :
  !    Note : In comparison with ``sasmtrphif'', this specialized version 
  !         : is usually three times faster.

  interface sasqtrphif
    module procedure rstr_scsqrf     ! Real    in Single    precision
    module procedure cstr_scsqrf     ! Complex in Single    precision
    module procedure rwtr_scsqrf     ! Real    in Double    precision
    module procedure cwtr_scsqrf     ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rttr_scsqrf     ! Real    in Triple    precision
    module procedure cttr_scsqrf     ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqtr_scsqrf     ! Real    in Quadruple precision
    module procedure cqtr_scsqrf     ! Complex in Quadruple precision
#endif
  end interface sasqtrphif


  !    Name : sasmtrphi ( a generic name )
  !         :
  !   Usage : call sasmtrphi(upto,scale,a,phin,err)
  !         :
  ! Purpose : This subroutine computes the matrix values of phi-functions for a
  !         : diagonal or a square matrix by the modified scaling and squaring
  !         : method.
  !         :
  !   Input : - ``upto'' is an integer, and specifies the maximum index of
  !         :   phi-functions to be computed simultaneously (0 <= upto <= 5).
  !         :
  !         : - ``A'' is a real or complex array whose shape is
  !         :   either A(1:N) or A(1:N,1:N).
  !         :
  !         : - ``scale'' is a real or complex scalar that should multiply the
  !         :   entries of ``A''. The type of ``scale'' must be the same as
  !         :   that of the entries of ``A''.
  !         :
  !  Output : - When the shape of ``A'' is (1:N), it is regarded as a diagonal
  !         :   matrix, and the output ``phin'' has the shape (1:N,0:upto).
  !         :   The subarray, phin(1:N,k), stores the value of phi_k(scale*A)
  !         :   for 0 <= k <= upto.
  !         :
  !         : - When the shape of ``A'' is (1:N,1:N), it is regarded as a square
  !         :   matrix, and the output ``phin'' has the shape (1:N,1:N,0:upto).
  !         :   The subarray, phin(1:N,1:N,k), is the value of phi_k(scale*A)
  !         :   for 0 <= k <= upto.
  !         :
  !         : - ``err'' is of type(mcpsqrlog), and stores error and log
  !         :    informations. The definition of the type is found in the file,
  !         :   ``mtrcfgphilog.f95''.

  interface sasmtrphi
    module procedure rsdg_scsqr      ! Real    Diagonal in Single    precision
    module procedure csdg_scsqr      ! Complex Diagonal in Single    precision
    module procedure rwdg_scsqr      ! Real    Diagonal in Double    precision
    module procedure cwdg_scsqr      ! Complex Diagonal in Double    precision
#ifdef __USE_TPREC
    module procedure rtdg_scsqr      ! Real    Diagonal in Triple    precision
    module procedure ctdg_scsqr      ! Complex Diagonal in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqdg_scsqr      ! Real    Diagonal in Quadruple precision
    module procedure cqdg_scsqr      ! Complex Diagonal in Quadruple precision
#endif

    module procedure rssq_scsqr      ! Real    Square   in Single    precision
    module procedure cssq_scsqr      ! Complex Square   in Single    precision
    module procedure rwsq_scsqr      ! Real    Square   in Double    precision
    module procedure cwsq_scsqr      ! Complex Square   in Double    precision
#ifdef __USE_TPREC
    module procedure rtsq_scsqr      ! Real    Square   in Triple    precision
    module procedure ctsq_scsqr      ! Complex Square   in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqsq_scsqr      ! Real    Square   in Quadruple precision
    module procedure cqsq_scsqr      ! Complex Square   in Quadruple precision
#endif
  end interface sasmtrphi


  !    Name : sasqtrphi ( a generic name )
  !         :
  !   Usage : call sasqtrphi(upto,scale,a,phin,err)
  !         :
  ! Purpose : This subroutine computes the matrix values of phi-functions
  !         : for an upper quasi-triangular matrix by the modified scaling
  !         : and squaring method.
  !         :
  !   Input : - ``upto'' is an integer, and specifies the maximum index of
  !         :   phi-functions to be computed simultaneously (0 <= upto <= 5).
  !         :
  !         : - ``A'' is a real or complex array whose shape is (1:N,1:N).
  !         :   This array stores an upper quasi-triangular matrix.
  !         :
  !         : - ``scale'' is a real or complex scalar that should multiply the
  !         :   entries of ``A''. The type of ``scale'' must be the same as
  !         :   that of the entries of ``A''.
  !         :
  !  Output : - The subarray, phin(1:N,1:N,k), stores the value of
  !         :   phi_k(scale*A(1:N,1:N)) for 0 <= k <= upto.
  !         :
  !         : - ``err'' is of type(mcpsqrlog) and stores error and log
  !         :   informations. The definition of the type is found in the file,
  !         :   ``mtrcfgphilog.f95''.
  !         :
  !    Note : In comparison with ``sasmtrphi'', this specialized version 
  !         : is usually three times faster.

  interface sasqtrphi
    module procedure rstr_scsqr      ! Real    in Single    precision
    module procedure cstr_scsqr      ! Complex in Single    precision
    module procedure rwtr_scsqr      ! Real    in Double    precision
    module procedure cwtr_scsqr      ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rttr_scsqr      ! Real    in Triple    precision
    module procedure cttr_scsqr      ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqtr_scsqr      ! Real    in Quadruple precision
    module procedure cqtr_scsqr      ! Complex in Quadruple precision
#endif
  end interface sasqtrphi


  !    Name : sqrmtrphif ( a generic name )
  !         :
  !   Usage : phiv = sqrmtrphif(upto,phiu)
  !         :
  ! Purpose : This function squares the phi-functions, i.e.,
  !         : { phi_k(2A); 0 <= k <= upto } is computed from
  !         : { phi_k( A); 0 <= k <= upto }.
  !         : The matrix A is either diagonal or square.
  !         :
  !   Input : - ``upto'' is an integer, and specifies the maximum index of
  !         :   phi-functions to be squared simultaneously (0 <= upto <= 5).
  !         :
  !         : - ``phiu'' is a real or complex array, whose shape is either
  !         :   (1:N,0:upto) or (1:N,1:N,0:upto).
  !         :
  !         :   When the shape is (1:N,0:upto), it stores
  !         :     ( phi_0(A), phi_1(A), ..., phi_{upto}(A) )
  !         :   for a diagonal matrix A(1:N). The subarray, phi(1:N,k), stores
  !         :   the value of phi_k(A).
  !         :
  !         :   When the shape is (1:N,1:N,0:upto), it stores
  !         :     ( phi_0(A), phi_1(A), ..., phi_{upto}(A) )
  !         :   for a square matrix A(1:N,1:N). The subarray, phi(1:N,1:N,k),
  !         :   stores the value of phi_k(A).
  !         :
  !  Output : ``phiv'' is a real or complex array, whose shape must be the
  !         : same as ``phiu''. For the diagonal case, phiv(1:N,k) stores
  !         : phi_k(2*A), while for the square case, phiv(1:N,1:N,k) stores
  !         : phi_k(2*A).

  interface sqrmtrphif
    module procedure rsdg_squaref    ! Real    Diagonal in Single    precision
    module procedure csdg_squaref    ! Complex Diagonal in Single    precision
    module procedure rwdg_squaref    ! Real    Diagonal in Double    precision
    module procedure cwdg_squaref    ! Complex Diagonal in Double    precision
#ifdef __USE_TPREC
    module procedure rtdg_squaref    ! Real    Diagonal in Triple    precision
    module procedure ctdg_squaref    ! Complex Diagonal in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqdg_squaref    ! Real    Diagonal in Quadruple precision
    module procedure cqdg_squaref    ! Complex Diagonal in Quadruple precision
#endif

    module procedure rssq_squaref    ! Real    Square   in Single    precision
    module procedure cssq_squaref    ! Complex Square   in Single    precision
    module procedure rwsq_squaref    ! Real    Square   in Double    precision
    module procedure cwsq_squaref    ! Complex Square   in Double    precision
#ifdef __USE_TPREC
    module procedure rtsq_squaref    ! Real    Square   in Triple    precision
    module procedure ctsq_squaref    ! Complex Square   in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqsq_squaref    ! Real    Square   in Quadruple precision
    module procedure cqsq_squaref    ! Complex Square   in Quadruple precision
#endif
  end interface sqrmtrphif


  !    Name : sqrqtrphif ( a generic name )
  !         :
  !   Usage : phiv = sqrqtrphif(upto,phiu)
  !         :
  ! Purpose : This function squares the phi-functions, i.e.,
  !         : { phi_k(2A); 0 <= k <= upto } is computed from
  !         : { phi_k( A); 0 <= k <= upto }.
  !         : The matrix A is upper quasi-triangular.
  !         :
  !   Input : - ``upto'' is an integer, and specifies the maximum index of
  !         :   phi-functions to be squared simultaneously (0 <= upto <= 5).
  !         :
  !         : - ``phiu'' is a real or complex array, whose shape is
  !         :   (1:N,1:N,0:upto). This array stores
  !         :     ( phi_0(A), phi_1(A), ..., phi_{upto}(A) )
  !         :   for an upper quasi-triangular matrix, A(1:N,1:N).
  !         :   The subarray, phi(1:N,1:N,k), stores the value of phi_k(A).
  !         :
  !  Output : ``phiv'' is a real or complex array, whose shape must be the
  !         : same as ``phiu''. The subarray, phiv(1:N,1:N,k), stores
  !         : the value of  phi_k(2*A).

  interface sqrqtrphif
    module procedure rstr_squaref    ! Real    in Single    precision
    module procedure cstr_squaref    ! Complex in Single    precision
    module procedure rwtr_squaref    ! Real    in Double    precision
    module procedure cwtr_squaref    ! Complex in Double    precision
#ifdef __USE_TPREC
    module procedure rttr_squaref    ! Real    in Triple    precision
    module procedure cttr_squaref    ! Complex in Triple    precision
#endif
#ifdef __USE_QPREC
    module procedure rqtr_squaref    ! Real    in Quadruple precision
    module procedure cqtr_squaref    ! Complex in Quadruple precision
#endif
  end interface sqrqtrphif

  ! the highest order of rational approximants used in each precision

  integer, parameter :: sp_hord =  7 ! for single    precision
  integer, parameter :: wp_hord = 13 ! for double    precision
  integer, parameter :: tp_hord = 13 ! for triple    precision
  integer, parameter :: qp_hord = 17 ! for quadruple precision

  private
  public  :: sasmtrphif, sasmtrphi   ! for a square or diagonal matrix
  public  :: sasqtrphif, sasqtrphi   ! for a quasi-upper triangular matrix
  public  :: sqrmtrphif              ! squaring for square or diagonal
  public  :: sqrqtrphif              ! squaring for quasi-upper triangular

contains


!    Name : rsdg_scsqrf
!   Usage : values = rsdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, single precision, diagonal matrix.

pure function rsdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: scale
  real   (kind=sp), intent(in) :: arg(1:)
  real   (kind=sp)             :: rsdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call rsdg_scsqr(upto,scale,arg,rsdg_scsqrf,err)
end function rsdg_scsqrf


!    Name : rsdg_scsqr
!   Usage : call rsdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, single precision, diagonal matrix.

pure subroutine rsdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=sp), intent(in ) :: scale
  real   (kind=sp), intent(in ) :: arg(1:)
  real   (kind=sp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=sp) :: a(size(arg,1))
  real   (kind=sp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_sp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rsdg_normi(a))
  err % ah_norm(1,0) = real(rsdg_norm1(a))

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rsdg_norm1(a),rsdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rsdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rsdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rsdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(rsdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rsdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rsdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(rsdg_norm1(values(:,i)))
    end do
  end if

end subroutine rsdg_scsqr


!    Name : rsdg_squaref
!   Usage : value = rsdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, single precision, diagonal matrix.

pure function rsdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: arg(1:,0:)
  real   (kind=sp)             :: rsdg_squaref(size(arg,1),0:upto)

  rsdg_squaref = arg
  call rsdg_square(upto,rsdg_squaref)
end function rsdg_squaref


!    Name : rsdg_square
!   Usage : call rsdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, single precision, diagonal matrix.

pure subroutine rsdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=sp), intent(inout) :: arg(1:,0:)

  real   (kind=sp) :: Iexp(size(arg,1))

  Iexp = rsdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_sp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_sp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_sp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_sp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_sp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_sp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_sp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_sp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_sp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_sp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_sp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_sp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_sp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_sp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_sp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine rsdg_square


!    Name : rssq_scsqrf
!   Usage : values = rssq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, single precision, square matrix.

pure function rssq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: scale
  real   (kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: rssq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rssq_scsqr(upto,scale,arg,rssq_scsqrf,err)
end function rssq_scsqrf


!    Name : rssq_scsqr
!   Usage : call rssq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, single precision, square matrix.

pure subroutine rssq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=sp), intent(in ) :: scale
  real   (kind=sp), intent(in ) :: arg(1:,1:)
  real   (kind=sp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=sp) :: a(size(arg,1),size(arg,2))
  real   (kind=sp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_sp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rssq_normi(a))
  err % ah_norm(1,0) = real(rssq_norm1(a))

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rssq_norm1(a),rssq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rssq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rssq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rssq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rssq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rssq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rssq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rssq_norm1(values(:,:,i)))
    end do
  end if

end subroutine rssq_scsqr


!    Name : rssq_squaref
!   Usage : value = rssq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, single precision, square matrix.

pure function rssq_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: arg(1:,1:,0:)
  real   (kind=sp)             :: rssq_squaref(size(arg,1),size(arg,2),0:upto)

  rssq_squaref = arg
  call rssq_square(upto,rssq_squaref)
end function rssq_squaref


!    Name : rssq_square
!   Usage : call rssq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, single precision, square matrix.

pure subroutine rssq_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=sp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=sp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = rssq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_sp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_sp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_sp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_sp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_sp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_sp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine rssq_square


!    Name : rstr_scsqrf
!   Usage : values = rstr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, single precision, quasi upper triangular matrix.

pure function rstr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: scale
  real   (kind=sp), intent(in) :: arg(1:,1:)
  real   (kind=sp)             :: rstr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rstr_scsqr(upto,scale,arg,rstr_scsqrf,err)
end function rstr_scsqrf


!    Name : rstr_scsqr
!   Usage : call rstr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, single precision, quasi upper triangular matrix.

pure subroutine rstr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=sp), intent(in ) :: scale
  real   (kind=sp), intent(in ) :: arg(1:,1:)
  real   (kind=sp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=sp) :: a(size(arg,1),size(arg,2))
  real   (kind=sp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_sp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rstr_normi(a))
  err % ah_norm(1,0) = real(rstr_norm1(a))

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rstr_norm1(a),rstr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rstr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rstr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rstr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rstr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rstr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rstr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rstr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_sp)) then
      values(j+1,j,:) = 0.0_sp
    end if
  end do
end subroutine rstr_scsqr


!    Name : rstr_squaref
!   Usage : value = rstr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, single precision, quasi upper triangular matrix.

pure function rstr_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=sp), intent(in) :: arg(1:,1:,0:)
  real   (kind=sp)             :: rstr_squaref(size(arg,1),size(arg,2),0:upto)

  rstr_squaref = arg
  call rstr_square(upto,rstr_squaref)
end function rstr_squaref


!    Name : rstr_square
!   Usage : call rstr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, single precision, quasi upper triangular matrix.

pure subroutine rstr_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=sp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=sp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_sp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = rstr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_sp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_sp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_sp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_sp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_sp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_sp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = 0.0_sp
  end do
end subroutine rstr_square


!    Name : csdg_scsqrf
!   Usage : values = csdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, single precision, diagonal matrix.

pure function csdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: scale
  complex(kind=sp), intent(in) :: arg(1:)
  complex(kind=sp)             :: csdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call csdg_scsqr(upto,scale,arg,csdg_scsqrf,err)
end function csdg_scsqrf


!    Name : csdg_scsqr
!   Usage : call csdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, single precision, diagonal matrix.

pure subroutine csdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=sp), intent(in ) :: scale
  complex(kind=sp), intent(in ) :: arg(1:)
  complex(kind=sp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=sp) :: a(size(arg,1))
  real   (kind=sp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_sp,0.0_sp,sp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(csdg_normi(a))
  err % ah_norm(1,0) = real(csdg_norm1(a))

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(csdg_norm1(a),csdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call csdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call csdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(csdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(csdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call csdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(csdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(csdg_norm1(values(:,i)))
    end do
  end if

end subroutine csdg_scsqr


!    Name : csdg_squaref
!   Usage : value = csdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, single precision, diagonal matrix.

pure function csdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: arg(1:,0:)
  complex(kind=sp)             :: csdg_squaref(size(arg,1),0:upto)

  csdg_squaref = arg
  call csdg_square(upto,csdg_squaref)
end function csdg_squaref


!    Name : csdg_square
!   Usage : call csdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, single precision, diagonal matrix.

pure subroutine csdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=sp), intent(inout) :: arg(1:,0:)

  complex(kind=sp) :: Iexp(size(arg,1))

  Iexp = csdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_sp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_sp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_sp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_sp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_sp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_sp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_sp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_sp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_sp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_sp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_sp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_sp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_sp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_sp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_sp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine csdg_square


!    Name : cssq_scsqrf
!   Usage : values = cssq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, single precision, square matrix.

pure function cssq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: scale
  complex(kind=sp), intent(in) :: arg(1:,1:)
  complex(kind=sp)             :: cssq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cssq_scsqr(upto,scale,arg,cssq_scsqrf,err)
end function cssq_scsqrf


!    Name : cssq_scsqr
!   Usage : call cssq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, single precision, square matrix.

pure subroutine cssq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=sp), intent(in ) :: scale
  complex(kind=sp), intent(in ) :: arg(1:,1:)
  complex(kind=sp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=sp) :: a(size(arg,1),size(arg,2))
  real   (kind=sp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_sp,0.0_sp,sp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cssq_normi(a))
  err % ah_norm(1,0) = real(cssq_norm1(a))

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cssq_norm1(a),cssq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cssq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cssq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cssq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(cssq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cssq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cssq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(cssq_norm1(values(:,:,i)))
    end do
  end if

end subroutine cssq_scsqr


!    Name : cssq_squaref
!   Usage : value = cssq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, single precision, square matrix.

pure function cssq_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: arg(1:,1:,0:)
  complex(kind=sp)             :: cssq_squaref(size(arg,1),size(arg,2),0:upto)

  cssq_squaref = arg
  call cssq_square(upto,cssq_squaref)
end function cssq_squaref


!    Name : cssq_square
!   Usage : call cssq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, single precision, square matrix.

pure subroutine cssq_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=sp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=sp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = cssq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_sp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_sp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_sp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_sp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_sp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_sp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine cssq_square


!    Name : cstr_scsqrf
!   Usage : values = cstr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, single precision, quasi upper triangular matrix.

pure function cstr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: scale
  complex(kind=sp), intent(in) :: arg(1:,1:)
  complex(kind=sp)             :: cstr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cstr_scsqr(upto,scale,arg,cstr_scsqrf,err)
end function cstr_scsqrf


!    Name : cstr_scsqr
!   Usage : call cstr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, single precision, quasi upper triangular matrix.

pure subroutine cstr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=sp), intent(in ) :: scale
  complex(kind=sp), intent(in ) :: arg(1:,1:)
  complex(kind=sp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=sp) :: a(size(arg,1),size(arg,2))
  real   (kind=sp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_sp,0.0_sp,sp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cstr_normi(a))
  err % ah_norm(1,0) = real(cstr_norm1(a))

  cta = sp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cstr_norm1(a),cstr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cstr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = sp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_sp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cstr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cstr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(cstr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cstr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cstr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(cstr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_sp)) then
      values(j+1,j,:) = cmplx(0.0_sp,0.0_sp,sp)
    end if
  end do
end subroutine cstr_scsqr


!    Name : cstr_squaref
!   Usage : value = cstr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, single precision, quasi upper triangular matrix.

pure function cstr_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=sp), intent(in) :: arg(1:,1:,0:)
  complex(kind=sp)             :: cstr_squaref(size(arg,1),size(arg,2),0:upto)

  cstr_squaref = arg
  call cstr_square(upto,cstr_squaref)
end function cstr_squaref


!    Name : cstr_square
!   Usage : call cstr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, single precision, quasi upper triangular matrix.

pure subroutine cstr_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=sp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=sp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_sp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = cstr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_sp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_sp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_sp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_sp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_sp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_sp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_sp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_sp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_sp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = cmplx(0.0_sp,0.0_sp,sp)
  end do
end subroutine cstr_square


!    Name : rwdg_scsqrf
!   Usage : values = rwdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, double precision, diagonal matrix.

pure function rwdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: scale
  real   (kind=wp), intent(in) :: arg(1:)
  real   (kind=wp)             :: rwdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call rwdg_scsqr(upto,scale,arg,rwdg_scsqrf,err)
end function rwdg_scsqrf


!    Name : rwdg_scsqr
!   Usage : call rwdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, double precision, diagonal matrix.

pure subroutine rwdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=wp), intent(in ) :: scale
  real   (kind=wp), intent(in ) :: arg(1:)
  real   (kind=wp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=wp) :: a(size(arg,1))
  real   (kind=wp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_wp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rwdg_normi(a))
  err % ah_norm(1,0) = real(rwdg_norm1(a))

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rwdg_norm1(a),rwdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rwdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rwdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rwdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(rwdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rwdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rwdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(rwdg_norm1(values(:,i)))
    end do
  end if

end subroutine rwdg_scsqr


!    Name : rwdg_squaref
!   Usage : value = rwdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, double precision, diagonal matrix.

pure function rwdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: arg(1:,0:)
  real   (kind=wp)             :: rwdg_squaref(size(arg,1),0:upto)

  rwdg_squaref = arg
  call rwdg_square(upto,rwdg_squaref)
end function rwdg_squaref


!    Name : rwdg_square
!   Usage : call rwdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, double precision, diagonal matrix.

pure subroutine rwdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=wp), intent(inout) :: arg(1:,0:)

  real   (kind=wp) :: Iexp(size(arg,1))

  Iexp = rwdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_wp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_wp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_wp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_wp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_wp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_wp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_wp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_wp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_wp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_wp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_wp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_wp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_wp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_wp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_wp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine rwdg_square


!    Name : rwsq_scsqrf
!   Usage : values = rwsq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, double precision, square matrix.

pure function rwsq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: scale
  real   (kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: rwsq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rwsq_scsqr(upto,scale,arg,rwsq_scsqrf,err)
end function rwsq_scsqrf


!    Name : rwsq_scsqr
!   Usage : call rwsq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, double precision, square matrix.

pure subroutine rwsq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=wp), intent(in ) :: scale
  real   (kind=wp), intent(in ) :: arg(1:,1:)
  real   (kind=wp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=wp) :: a(size(arg,1),size(arg,2))
  real   (kind=wp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_wp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rwsq_normi(a))
  err % ah_norm(1,0) = real(rwsq_norm1(a))

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rwsq_norm1(a),rwsq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rwsq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rwsq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rwsq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rwsq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rwsq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rwsq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rwsq_norm1(values(:,:,i)))
    end do
  end if

end subroutine rwsq_scsqr


!    Name : rwsq_squaref
!   Usage : value = rwsq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, double precision, square matrix.

pure function rwsq_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: arg(1:,1:,0:)
  real   (kind=wp)             :: rwsq_squaref(size(arg,1),size(arg,2),0:upto)

  rwsq_squaref = arg
  call rwsq_square(upto,rwsq_squaref)
end function rwsq_squaref


!    Name : rwsq_square
!   Usage : call rwsq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, double precision, square matrix.

pure subroutine rwsq_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=wp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=wp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = rwsq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_wp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_wp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_wp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_wp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_wp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_wp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine rwsq_square


!    Name : rwtr_scsqrf
!   Usage : values = rwtr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, double precision, quasi upper triangular matrix.

pure function rwtr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: scale
  real   (kind=wp), intent(in) :: arg(1:,1:)
  real   (kind=wp)             :: rwtr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rwtr_scsqr(upto,scale,arg,rwtr_scsqrf,err)
end function rwtr_scsqrf


!    Name : rwtr_scsqr
!   Usage : call rwtr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, double precision, quasi upper triangular matrix.

pure subroutine rwtr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=wp), intent(in ) :: scale
  real   (kind=wp), intent(in ) :: arg(1:,1:)
  real   (kind=wp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=wp) :: a(size(arg,1),size(arg,2))
  real   (kind=wp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_wp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rwtr_normi(a))
  err % ah_norm(1,0) = real(rwtr_norm1(a))

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rwtr_norm1(a),rwtr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rwtr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rwtr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rwtr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rwtr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rwtr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rwtr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rwtr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_wp)) then
      values(j+1,j,:) = 0.0_wp
    end if
  end do
end subroutine rwtr_scsqr


!    Name : rwtr_squaref
!   Usage : value = rwtr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, double precision, quasi upper triangular matrix.

pure function rwtr_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=wp), intent(in) :: arg(1:,1:,0:)
  real   (kind=wp)             :: rwtr_squaref(size(arg,1),size(arg,2),0:upto)

  rwtr_squaref = arg
  call rwtr_square(upto,rwtr_squaref)
end function rwtr_squaref


!    Name : rwtr_square
!   Usage : call rwtr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, double precision, quasi upper triangular matrix.

pure subroutine rwtr_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=wp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=wp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_wp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = rwtr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_wp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_wp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_wp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_wp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_wp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_wp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = 0.0_wp
  end do
end subroutine rwtr_square


!    Name : cwdg_scsqrf
!   Usage : values = cwdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, double precision, diagonal matrix.

pure function cwdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: scale
  complex(kind=wp), intent(in) :: arg(1:)
  complex(kind=wp)             :: cwdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call cwdg_scsqr(upto,scale,arg,cwdg_scsqrf,err)
end function cwdg_scsqrf


!    Name : cwdg_scsqr
!   Usage : call cwdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, double precision, diagonal matrix.

pure subroutine cwdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=wp), intent(in ) :: scale
  complex(kind=wp), intent(in ) :: arg(1:)
  complex(kind=wp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=wp) :: a(size(arg,1))
  real   (kind=wp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_wp,0.0_wp,wp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cwdg_normi(a))
  err % ah_norm(1,0) = real(cwdg_norm1(a))

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cwdg_norm1(a),cwdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cwdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cwdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cwdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(cwdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cwdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cwdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(cwdg_norm1(values(:,i)))
    end do
  end if

end subroutine cwdg_scsqr


!    Name : cwdg_squaref
!   Usage : value = cwdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, double precision, diagonal matrix.

pure function cwdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: arg(1:,0:)
  complex(kind=wp)             :: cwdg_squaref(size(arg,1),0:upto)

  cwdg_squaref = arg
  call cwdg_square(upto,cwdg_squaref)
end function cwdg_squaref


!    Name : cwdg_square
!   Usage : call cwdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, double precision, diagonal matrix.

pure subroutine cwdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=wp), intent(inout) :: arg(1:,0:)

  complex(kind=wp) :: Iexp(size(arg,1))

  Iexp = cwdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_wp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_wp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_wp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_wp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_wp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_wp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_wp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_wp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_wp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_wp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_wp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_wp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_wp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_wp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_wp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine cwdg_square


!    Name : cwsq_scsqrf
!   Usage : values = cwsq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, double precision, square matrix.

pure function cwsq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: scale
  complex(kind=wp), intent(in) :: arg(1:,1:)
  complex(kind=wp)             :: cwsq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cwsq_scsqr(upto,scale,arg,cwsq_scsqrf,err)
end function cwsq_scsqrf


!    Name : cwsq_scsqr
!   Usage : call cwsq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, double precision, square matrix.

pure subroutine cwsq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=wp), intent(in ) :: scale
  complex(kind=wp), intent(in ) :: arg(1:,1:)
  complex(kind=wp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=wp) :: a(size(arg,1),size(arg,2))
  real   (kind=wp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_wp,0.0_wp,wp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cwsq_normi(a))
  err % ah_norm(1,0) = real(cwsq_norm1(a))

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cwsq_norm1(a),cwsq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cwsq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cwsq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cwsq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(cwsq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cwsq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cwsq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(cwsq_norm1(values(:,:,i)))
    end do
  end if

end subroutine cwsq_scsqr


!    Name : cwsq_squaref
!   Usage : value = cwsq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, double precision, square matrix.

pure function cwsq_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: arg(1:,1:,0:)
  complex(kind=wp)             :: cwsq_squaref(size(arg,1),size(arg,2),0:upto)

  cwsq_squaref = arg
  call cwsq_square(upto,cwsq_squaref)
end function cwsq_squaref


!    Name : cwsq_square
!   Usage : call cwsq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, double precision, square matrix.

pure subroutine cwsq_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=wp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=wp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = cwsq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_wp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_wp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_wp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_wp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_wp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_wp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine cwsq_square


!    Name : cwtr_scsqrf
!   Usage : values = cwtr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, double precision, quasi upper triangular matrix.

pure function cwtr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: scale
  complex(kind=wp), intent(in) :: arg(1:,1:)
  complex(kind=wp)             :: cwtr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cwtr_scsqr(upto,scale,arg,cwtr_scsqrf,err)
end function cwtr_scsqrf


!    Name : cwtr_scsqr
!   Usage : call cwtr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, double precision, quasi upper triangular matrix.

pure subroutine cwtr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=wp), intent(in ) :: scale
  complex(kind=wp), intent(in ) :: arg(1:,1:)
  complex(kind=wp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=wp) :: a(size(arg,1),size(arg,2))
  real   (kind=wp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_wp,0.0_wp,wp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cwtr_normi(a))
  err % ah_norm(1,0) = real(cwtr_norm1(a))

  cta = wp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cwtr_norm1(a),cwtr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cwtr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = wp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_wp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cwtr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cwtr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(cwtr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cwtr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cwtr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(cwtr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_wp)) then
      values(j+1,j,:) = cmplx(0.0_wp,0.0_wp,wp)
    end if
  end do
end subroutine cwtr_scsqr


!    Name : cwtr_squaref
!   Usage : value = cwtr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, double precision, quasi upper triangular matrix.

pure function cwtr_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=wp), intent(in) :: arg(1:,1:,0:)
  complex(kind=wp)             :: cwtr_squaref(size(arg,1),size(arg,2),0:upto)

  cwtr_squaref = arg
  call cwtr_square(upto,cwtr_squaref)
end function cwtr_squaref


!    Name : cwtr_square
!   Usage : call cwtr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, double precision, quasi upper triangular matrix.

pure subroutine cwtr_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=wp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=wp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_wp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = cwtr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_wp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_wp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_wp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_wp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_wp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_wp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_wp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_wp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_wp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = cmplx(0.0_wp,0.0_wp,wp)
  end do
end subroutine cwtr_square

#ifdef __USE_TPREC

!    Name : rtdg_scsqrf
!   Usage : values = rtdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, triple precision, diagonal matrix.

pure function rtdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: scale
  real   (kind=tp), intent(in) :: arg(1:)
  real   (kind=tp)             :: rtdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call rtdg_scsqr(upto,scale,arg,rtdg_scsqrf,err)
end function rtdg_scsqrf


!    Name : rtdg_scsqr
!   Usage : call rtdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, triple precision, diagonal matrix.

pure subroutine rtdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=tp), intent(in ) :: scale
  real   (kind=tp), intent(in ) :: arg(1:)
  real   (kind=tp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=tp) :: a(size(arg,1))
  real   (kind=tp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_tp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rtdg_normi(a))
  err % ah_norm(1,0) = real(rtdg_norm1(a))

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rtdg_norm1(a),rtdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rtdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rtdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rtdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(rtdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rtdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rtdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(rtdg_norm1(values(:,i)))
    end do
  end if

end subroutine rtdg_scsqr


!    Name : rtdg_squaref
!   Usage : value = rtdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, triple precision, diagonal matrix.

pure function rtdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: arg(1:,0:)
  real   (kind=tp)             :: rtdg_squaref(size(arg,1),0:upto)

  rtdg_squaref = arg
  call rtdg_square(upto,rtdg_squaref)
end function rtdg_squaref


!    Name : rtdg_square
!   Usage : call rtdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, triple precision, diagonal matrix.

pure subroutine rtdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=tp), intent(inout) :: arg(1:,0:)

  real   (kind=tp) :: Iexp(size(arg,1))

  Iexp = rtdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_tp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_tp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_tp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_tp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_tp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_tp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_tp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_tp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_tp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_tp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_tp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_tp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_tp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_tp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_tp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine rtdg_square


!    Name : rtsq_scsqrf
!   Usage : values = rtsq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, triple precision, square matrix.

pure function rtsq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: scale
  real   (kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: rtsq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rtsq_scsqr(upto,scale,arg,rtsq_scsqrf,err)
end function rtsq_scsqrf


!    Name : rtsq_scsqr
!   Usage : call rtsq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, triple precision, square matrix.

pure subroutine rtsq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=tp), intent(in ) :: scale
  real   (kind=tp), intent(in ) :: arg(1:,1:)
  real   (kind=tp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=tp) :: a(size(arg,1),size(arg,2))
  real   (kind=tp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_tp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rtsq_normi(a))
  err % ah_norm(1,0) = real(rtsq_norm1(a))

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rtsq_norm1(a),rtsq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rtsq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rtsq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rtsq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rtsq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rtsq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rtsq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rtsq_norm1(values(:,:,i)))
    end do
  end if

end subroutine rtsq_scsqr


!    Name : rtsq_squaref
!   Usage : value = rtsq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, triple precision, square matrix.

pure function rtsq_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: arg(1:,1:,0:)
  real   (kind=tp)             :: rtsq_squaref(size(arg,1),size(arg,2),0:upto)

  rtsq_squaref = arg
  call rtsq_square(upto,rtsq_squaref)
end function rtsq_squaref


!    Name : rtsq_square
!   Usage : call rtsq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, triple precision, square matrix.

pure subroutine rtsq_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=tp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=tp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = rtsq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_tp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_tp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_tp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_tp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_tp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_tp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine rtsq_square


!    Name : rttr_scsqrf
!   Usage : values = rttr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, triple precision, quasi upper triangular matrix.

pure function rttr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: scale
  real   (kind=tp), intent(in) :: arg(1:,1:)
  real   (kind=tp)             :: rttr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rttr_scsqr(upto,scale,arg,rttr_scsqrf,err)
end function rttr_scsqrf


!    Name : rttr_scsqr
!   Usage : call rttr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, triple precision, quasi upper triangular matrix.

pure subroutine rttr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=tp), intent(in ) :: scale
  real   (kind=tp), intent(in ) :: arg(1:,1:)
  real   (kind=tp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=tp) :: a(size(arg,1),size(arg,2))
  real   (kind=tp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_tp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rttr_normi(a))
  err % ah_norm(1,0) = real(rttr_norm1(a))

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rttr_norm1(a),rttr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rttr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rttr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rttr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rttr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rttr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rttr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rttr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_tp)) then
      values(j+1,j,:) = 0.0_tp
    end if
  end do
end subroutine rttr_scsqr


!    Name : rttr_squaref
!   Usage : value = rttr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, triple precision, quasi upper triangular matrix.

pure function rttr_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=tp), intent(in) :: arg(1:,1:,0:)
  real   (kind=tp)             :: rttr_squaref(size(arg,1),size(arg,2),0:upto)

  rttr_squaref = arg
  call rttr_square(upto,rttr_squaref)
end function rttr_squaref


!    Name : rttr_square
!   Usage : call rttr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, triple precision, quasi upper triangular matrix.

pure subroutine rttr_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=tp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=tp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_tp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = rttr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_tp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_tp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_tp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_tp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_tp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_tp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = 0.0_tp
  end do
end subroutine rttr_square


!    Name : ctdg_scsqrf
!   Usage : values = ctdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, triple precision, diagonal matrix.

pure function ctdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: scale
  complex(kind=tp), intent(in) :: arg(1:)
  complex(kind=tp)             :: ctdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call ctdg_scsqr(upto,scale,arg,ctdg_scsqrf,err)
end function ctdg_scsqrf


!    Name : ctdg_scsqr
!   Usage : call ctdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, triple precision, diagonal matrix.

pure subroutine ctdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=tp), intent(in ) :: scale
  complex(kind=tp), intent(in ) :: arg(1:)
  complex(kind=tp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=tp) :: a(size(arg,1))
  real   (kind=tp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_tp,0.0_tp,tp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(ctdg_normi(a))
  err % ah_norm(1,0) = real(ctdg_norm1(a))

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(ctdg_norm1(a),ctdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call ctdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call ctdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(ctdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(ctdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call ctdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(ctdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(ctdg_norm1(values(:,i)))
    end do
  end if

end subroutine ctdg_scsqr


!    Name : ctdg_squaref
!   Usage : value = ctdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, triple precision, diagonal matrix.

pure function ctdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: arg(1:,0:)
  complex(kind=tp)             :: ctdg_squaref(size(arg,1),0:upto)

  ctdg_squaref = arg
  call ctdg_square(upto,ctdg_squaref)
end function ctdg_squaref


!    Name : ctdg_square
!   Usage : call ctdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, triple precision, diagonal matrix.

pure subroutine ctdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=tp), intent(inout) :: arg(1:,0:)

  complex(kind=tp) :: Iexp(size(arg,1))

  Iexp = ctdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_tp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_tp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_tp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_tp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_tp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_tp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_tp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_tp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_tp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_tp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_tp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_tp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_tp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_tp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_tp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine ctdg_square


!    Name : ctsq_scsqrf
!   Usage : values = ctsq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, triple precision, square matrix.

pure function ctsq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: scale
  complex(kind=tp), intent(in) :: arg(1:,1:)
  complex(kind=tp)             :: ctsq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call ctsq_scsqr(upto,scale,arg,ctsq_scsqrf,err)
end function ctsq_scsqrf


!    Name : ctsq_scsqr
!   Usage : call ctsq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, triple precision, square matrix.

pure subroutine ctsq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=tp), intent(in ) :: scale
  complex(kind=tp), intent(in ) :: arg(1:,1:)
  complex(kind=tp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=tp) :: a(size(arg,1),size(arg,2))
  real   (kind=tp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_tp,0.0_tp,tp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(ctsq_normi(a))
  err % ah_norm(1,0) = real(ctsq_norm1(a))

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(ctsq_norm1(a),ctsq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call ctsq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call ctsq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(ctsq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(ctsq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call ctsq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(ctsq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(ctsq_norm1(values(:,:,i)))
    end do
  end if

end subroutine ctsq_scsqr


!    Name : ctsq_squaref
!   Usage : value = ctsq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, triple precision, square matrix.

pure function ctsq_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: arg(1:,1:,0:)
  complex(kind=tp)             :: ctsq_squaref(size(arg,1),size(arg,2),0:upto)

  ctsq_squaref = arg
  call ctsq_square(upto,ctsq_squaref)
end function ctsq_squaref


!    Name : ctsq_square
!   Usage : call ctsq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, triple precision, square matrix.

pure subroutine ctsq_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=tp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=tp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = ctsq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_tp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_tp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_tp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_tp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_tp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_tp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine ctsq_square


!    Name : cttr_scsqrf
!   Usage : values = cttr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, triple precision, quasi upper triangular matrix.

pure function cttr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: scale
  complex(kind=tp), intent(in) :: arg(1:,1:)
  complex(kind=tp)             :: cttr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cttr_scsqr(upto,scale,arg,cttr_scsqrf,err)
end function cttr_scsqrf


!    Name : cttr_scsqr
!   Usage : call cttr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, triple precision, quasi upper triangular matrix.

pure subroutine cttr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=tp), intent(in ) :: scale
  complex(kind=tp), intent(in ) :: arg(1:,1:)
  complex(kind=tp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=tp) :: a(size(arg,1),size(arg,2))
  real   (kind=tp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_tp,0.0_tp,tp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cttr_normi(a))
  err % ah_norm(1,0) = real(cttr_norm1(a))

  cta = tp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cttr_norm1(a),cttr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  end if
  if (5 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cttr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = tp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_tp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cttr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cttr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(cttr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cttr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cttr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(cttr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_tp)) then
      values(j+1,j,:) = cmplx(0.0_tp,0.0_tp,tp)
    end if
  end do
end subroutine cttr_scsqr


!    Name : cttr_squaref
!   Usage : value = cttr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, triple precision, quasi upper triangular matrix.

pure function cttr_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=tp), intent(in) :: arg(1:,1:,0:)
  complex(kind=tp)             :: cttr_squaref(size(arg,1),size(arg,2),0:upto)

  cttr_squaref = arg
  call cttr_square(upto,cttr_squaref)
end function cttr_squaref


!    Name : cttr_square
!   Usage : call cttr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, triple precision, quasi upper triangular matrix.

pure subroutine cttr_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=tp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=tp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_tp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = cttr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_tp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_tp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_tp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_tp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_tp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_tp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_tp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_tp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_tp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = cmplx(0.0_tp,0.0_tp,tp)
  end do
end subroutine cttr_square

#endif
#ifdef __USE_QPREC

!    Name : rqdg_scsqrf
!   Usage : values = rqdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, quadruple precision, diagonal matrix.

pure function rqdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: scale
  real   (kind=qp), intent(in) :: arg(1:)
  real   (kind=qp)             :: rqdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call rqdg_scsqr(upto,scale,arg,rqdg_scsqrf,err)
end function rqdg_scsqrf


!    Name : rqdg_scsqr
!   Usage : call rqdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, quadruple precision, diagonal matrix.

pure subroutine rqdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=qp), intent(in ) :: scale
  real   (kind=qp), intent(in ) :: arg(1:)
  real   (kind=qp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=qp) :: a(size(arg,1))
  real   (kind=qp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_qp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rqdg_normi(a))
  err % ah_norm(1,0) = real(rqdg_norm1(a))

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rqdg_norm1(a),rqdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rqdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rqdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rqdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(rqdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rqdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rqdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(rqdg_norm1(values(:,i)))
    end do
  end if

end subroutine rqdg_scsqr


!    Name : rqdg_squaref
!   Usage : value = rqdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, quadruple precision, diagonal matrix.

pure function rqdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: arg(1:,0:)
  real   (kind=qp)             :: rqdg_squaref(size(arg,1),0:upto)

  rqdg_squaref = arg
  call rqdg_square(upto,rqdg_squaref)
end function rqdg_squaref


!    Name : rqdg_square
!   Usage : call rqdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, quadruple precision, diagonal matrix.

pure subroutine rqdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=qp), intent(inout) :: arg(1:,0:)

  real   (kind=qp) :: Iexp(size(arg,1))

  Iexp = rqdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_qp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_qp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_qp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_qp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_qp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_qp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_qp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_qp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_qp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_qp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_qp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_qp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_qp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_qp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_qp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine rqdg_square


!    Name : rqsq_scsqrf
!   Usage : values = rqsq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, quadruple precision, square matrix.

pure function rqsq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: scale
  real   (kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: rqsq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rqsq_scsqr(upto,scale,arg,rqsq_scsqrf,err)
end function rqsq_scsqrf


!    Name : rqsq_scsqr
!   Usage : call rqsq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, quadruple precision, square matrix.

pure subroutine rqsq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=qp), intent(in ) :: scale
  real   (kind=qp), intent(in ) :: arg(1:,1:)
  real   (kind=qp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=qp) :: a(size(arg,1),size(arg,2))
  real   (kind=qp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_qp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rqsq_normi(a))
  err % ah_norm(1,0) = real(rqsq_norm1(a))

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rqsq_norm1(a),rqsq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rqsq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rqsq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rqsq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rqsq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rqsq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rqsq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rqsq_norm1(values(:,:,i)))
    end do
  end if

end subroutine rqsq_scsqr


!    Name : rqsq_squaref
!   Usage : value = rqsq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, quadruple precision, square matrix.

pure function rqsq_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: arg(1:,1:,0:)
  real   (kind=qp)             :: rqsq_squaref(size(arg,1),size(arg,2),0:upto)

  rqsq_squaref = arg
  call rqsq_square(upto,rqsq_squaref)
end function rqsq_squaref


!    Name : rqsq_square
!   Usage : call rqsq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, quadruple precision, square matrix.

pure subroutine rqsq_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=qp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=qp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = rqsq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_qp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_qp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_qp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_qp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_qp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_qp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine rqsq_square


!    Name : rqtr_scsqrf
!   Usage : values = rqtr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a real, quadruple precision, quasi upper triangular matrix.

pure function rqtr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: scale
  real   (kind=qp), intent(in) :: arg(1:,1:)
  real   (kind=qp)             :: rqtr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call rqtr_scsqr(upto,scale,arg,rqtr_scsqrf,err)
end function rqtr_scsqrf


!    Name : rqtr_scsqr
!   Usage : call rqtr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a real, quadruple precision, quasi upper triangular matrix.

pure subroutine rqtr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  real   (kind=qp), intent(in ) :: scale
  real   (kind=qp), intent(in ) :: arg(1:,1:)
  real   (kind=qp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  real   (kind=qp) :: a(size(arg,1),size(arg,2))
  real   (kind=qp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = 0.0_qp
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(rqtr_normi(a))
  err % ah_norm(1,0) = real(rqtr_norm1(a))

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(rqtr_norm1(a),rqtr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call rqtr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call rqtr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(rqtr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(rqtr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call rqtr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(rqtr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(rqtr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_qp)) then
      values(j+1,j,:) = 0.0_qp
    end if
  end do
end subroutine rqtr_scsqr


!    Name : rqtr_squaref
!   Usage : value = rqtr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a real, quadruple precision, quasi upper triangular matrix.

pure function rqtr_squaref(upto,arg)
  integer,          intent(in) :: upto
  real   (kind=qp), intent(in) :: arg(1:,1:,0:)
  real   (kind=qp)             :: rqtr_squaref(size(arg,1),size(arg,2),0:upto)

  rqtr_squaref = arg
  call rqtr_square(upto,rqtr_squaref)
end function rqtr_squaref


!    Name : rqtr_square
!   Usage : call rqtr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a real, quadruple precision, quasi upper triangular matrix.

pure subroutine rqtr_square(upto,arg)
  integer,          intent(in   ) :: upto
  real   (kind=qp), intent(inout) :: arg(1:,1:,0:)

  real   (kind=qp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_qp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = rqtr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_qp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_qp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_qp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_qp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_qp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_qp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = 0.0_qp
  end do
end subroutine rqtr_square


!    Name : cqdg_scsqrf
!   Usage : values = cqdg_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, quadruple precision, diagonal matrix.

pure function cqdg_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: scale
  complex(kind=qp), intent(in) :: arg(1:)
  complex(kind=qp)             :: cqdg_scsqrf(size(arg,1),0:upto)

  type(mcpsqrlog) :: err
  call cqdg_scsqr(upto,scale,arg,cqdg_scsqrf,err)
end function cqdg_scsqrf


!    Name : cqdg_scsqr
!   Usage : call cqdg_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, quadruple precision, diagonal matrix.

pure subroutine cqdg_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=qp), intent(in ) :: scale
  complex(kind=qp), intent(in ) :: arg(1:)
  complex(kind=qp), intent(out) :: values(size(arg,1),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=qp) :: a(size(arg,1))
  real   (kind=qp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_qp,0.0_qp,qp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cqdg_normi(a))
  err % ah_norm(1,0) = real(cqdg_norm1(a))

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cqdg_norm1(a),cqdg_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cqdg_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cqdg_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cqdg_normi(values(:,i)))
      err % sp_norm(1,i) = real(cqdg_norm1(values(:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cqdg_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cqdg_normi(values(:,i)))
      err % sq_norm(1,i) = real(cqdg_norm1(values(:,i)))
    end do
  end if

end subroutine cqdg_scsqr


!    Name : cqdg_squaref
!   Usage : value = cqdg_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, quadruple precision, diagonal matrix.

pure function cqdg_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: arg(1:,0:)
  complex(kind=qp)             :: cqdg_squaref(size(arg,1),0:upto)

  cqdg_squaref = arg
  call cqdg_square(upto,cqdg_squaref)
end function cqdg_squaref


!    Name : cqdg_square
!   Usage : call cqdg_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, quadruple precision, diagonal matrix.

pure subroutine cqdg_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=qp), intent(inout) :: arg(1:,0:)

  complex(kind=qp) :: Iexp(size(arg,1))

  Iexp = cqdg_eye(size(arg,1))
  Iexp = Iexp + arg(:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,5) = arg(:,5) .dgtimes. Iexp  /  32.0_qp
    arg(:,5) = arg(:,5) + arg(:,4)      /  32.0_qp
    arg(:,5) = arg(:,5) + arg(:,3)      /  64.0_qp
    arg(:,5) = arg(:,5) + arg(:,2)      / 192.0_qp
    arg(:,5) = arg(:,5) + arg(:,1)      / 768.0_qp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,4) = arg(:,4) .dgtimes. Iexp  /  16.0_qp
    arg(:,4) = arg(:,4) + arg(:,3)      /  16.0_qp
    arg(:,4) = arg(:,4) + arg(:,2)      /  32.0_qp
    arg(:,4) = arg(:,4) + arg(:,1)      /  96.0_qp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,3) = arg(:,3) .dgtimes. Iexp  /   8.0_qp
    arg(:,3) = arg(:,3) + arg(:,2)      /   8.0_qp
    arg(:,3) = arg(:,3) + arg(:,1)      /  16.0_qp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,2) = arg(:,2) .dgtimes. Iexp  /   4.0_qp
    arg(:,2) = arg(:,2) + arg(:,1)      /   4.0_qp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,1) = arg(:,1) .dgtimes. Iexp  /   2.0_qp
  end if
                                                   ! square the exponential
  arg(:,0) = arg(:,0) .dgtimes. arg(:,0)

end subroutine cqdg_square


!    Name : cqsq_scsqrf
!   Usage : values = cqsq_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, quadruple precision, square matrix.

pure function cqsq_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: scale
  complex(kind=qp), intent(in) :: arg(1:,1:)
  complex(kind=qp)             :: cqsq_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cqsq_scsqr(upto,scale,arg,cqsq_scsqrf,err)
end function cqsq_scsqrf


!    Name : cqsq_scsqr
!   Usage : call cqsq_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, quadruple precision, square matrix.

pure subroutine cqsq_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=qp), intent(in ) :: scale
  complex(kind=qp), intent(in ) :: arg(1:,1:)
  complex(kind=qp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=qp) :: a(size(arg,1),size(arg,2))
  real   (kind=qp) :: cta(3:17), rabs
  integer          :: order, spower, i

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_qp,0.0_qp,qp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cqsq_normi(a))
  err % ah_norm(1,0) = real(cqsq_norm1(a))

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cqsq_norm1(a),cqsq_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cqsq_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cqsq_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cqsq_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(cqsq_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cqsq_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cqsq_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(cqsq_norm1(values(:,:,i)))
    end do
  end if

end subroutine cqsq_scsqr


!    Name : cqsq_squaref
!   Usage : value = cqsq_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, quadruple precision, square matrix.

pure function cqsq_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: arg(1:,1:,0:)
  complex(kind=qp)             :: cqsq_squaref(size(arg,1),size(arg,2),0:upto)

  cqsq_squaref = arg
  call cqsq_square(upto,cqsq_squaref)
end function cqsq_squaref


!    Name : cqsq_square
!   Usage : call cqsq_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, quadruple precision, square matrix.

pure subroutine cqsq_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=qp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=qp) :: Iexp(size(arg,1),size(arg,2))

  Iexp = cqsq_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .sqtimes. Iexp    /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_qp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .sqtimes. Iexp    /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_qp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .sqtimes. Iexp    /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_qp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .sqtimes. Iexp    /   4.0_qp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_qp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .sqtimes. Iexp    /   2.0_qp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .sqtimes. arg(:,:,0)

end subroutine cqsq_square


!    Name : cqtr_scsqrf
!   Usage : values = cqtr_scsqrf(upto,scale,arg)
!         : This function is invoked through the generic name,
!         : ``sasmtrphif'' or ``sasqtrphif''.
! Purpose : This function computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : The subarray, values(*,k), stores the value of phi_k(scale*arg(*)).
!    Note : This function treats a complex, quadruple precision, quasi upper triangular matrix.

pure function cqtr_scsqrf(upto,scale,arg)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: scale
  complex(kind=qp), intent(in) :: arg(1:,1:)
  complex(kind=qp)             :: cqtr_scsqrf(size(arg,1),size(arg,2),0:upto)

  type(mcpsqrlog) :: err
  call cqtr_scsqr(upto,scale,arg,cqtr_scsqrf,err)
end function cqtr_scsqrf


!    Name : cqtr_scsqr
!   Usage : call cqtr_scsqr(upto,scale,arg,values,err)
!         : This subroutine is invoked through the generic name,
!         : ``sasmtrphi'' or ``sasqtrphi''.
! Purpose : This subroutine computes the matrix values of phi-functions
!         : by the modified scaling and squaring method.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be computed simultaneously (0 <= upto <= 5).
!         : - ``arg'' is the matrix argument.
!         : - ``scale'' is a scalar that should multiply the entries of
!         :   the matrix ``arg''.
!  Output : - The subarray, values(*,k), stores phi_k(scale*arg(*)).
!         : - ``err'' is of type(mcpsqrlog) and stores error and log
!         :   informations. The definition of the type is found in
!         :   the file ``mtrcfgphilog.f95''.
!    Note : This subroutine treats a complex, quadruple precision, quasi upper triangular matrix.

pure subroutine cqtr_scsqr(upto,scale,arg,values,err)
  integer,          intent(in ) :: upto
  complex(kind=qp), intent(in ) :: scale
  complex(kind=qp), intent(in ) :: arg(1:,1:)
  complex(kind=qp), intent(out) :: values(size(arg,1),size(arg,2),0:upto)
  type(mcpsqrlog),  intent(out) :: err

  complex(kind=qp) :: a(size(arg,1),size(arg,2))
  real   (kind=qp) :: cta(3:17), rabs
  integer          :: order, spower, i, j, subd(1:size(arg,1))

  ! Checking the arguments

  call init_mcpsqrlog(err)

  values = cmplx(0.0_qp,0.0_qp,qp)
  if (upto < 0 .or. 5 < upto) then; err % error = 1; return; end if
  if (size(arg,1) /= size(arg,2)) then; err % error = 2; return; end if

  subd = 0
  do i=1,size(arg,1)-1
    if ( tiny(abs(arg(i+1,i))) <= abs(arg(i+1,i)) ) subd(i) = 1
  end do
  do i=1,size(arg,1)-2
    if (subd(i) == 1 .and. subd(i+1) == 1) then
      err % error = 3
      return
    end if
  end do
  do i=3,size(arg,1)
    do j=1,i-2
      if ( tiny(abs(arg(i,j))) <= abs(arg(i,j)) ) then
        err % error = 3
        return
      end if
    end do
  end do
  if (size(arg,1) <= 0) then; err % error = 4; return; end if
  if (mpt_blksize <= 3) then; err % error = 5; return; end if

  ! Computations start from here.

  a = arg * scale

  err % upto         = upto
  err % ah_norm(0,0) = real(cqtr_normi(a))
  err % ah_norm(1,0) = real(cqtr_norm1(a))

  cta = qp_getthetamn(upto)                            ! the bounding parameters
  rabs = max(cqtr_norm1(a),cqtr_normi(a))              ! evaluating matrix norm

  order = -1                                           ! order of approximant
  if      (rabs < cta( 3)) then; order =  3
  else if (rabs < cta( 5)) then; order =  5
  else if (rabs < cta( 7)) then; order =  7
  else if (rabs < cta( 9)) then; order =  9
  else if (rabs < cta(13)) then; order = 13
  else if (rabs < cta(17)) then; order = 17
  end if
  if (3 <= upto .and. order == 3) order = 5

  if (0 < order) then                                  ! rational approximation

    ! When the norm of the matrix argument is sufficiently small,
    ! phi-functions are computed by rational approximations
    ! without using the scaling and squaring technique.

    err % rorder = order
    call cqtr_rationalphi(order,upto,a,values,err)

  else                                                 ! scaling and squaring

    ! When the norm of the matrix argument is large,
    ! - the matrix argument is divided by 2^s so that the norm of the
    !   scaled matrix, ``A/2^s'', is sufficiently small,
    ! - rational approximations are computed for { phi_k(A/2^s); 0<=k<=upto },
    ! - required approximations for { phi_k(A); 0 <= k <= upto } are computed
    !   by repeated squarings.

    order = qp_hord; err % rorder = order
    spower = ceiling(log(rabs/cta(order))/log(2.0_qp)) ! log2(|A|/ctam)

    a = a / (2**spower); err % spower = spower         ! scaling

    call cqtr_rationalphi(order,upto,a,values,err)     ! [m/m] approximant

    do i=0,upto
      err % sp_norm(0,i) = real(cqtr_normi(values(:,:,i)))
      err % sp_norm(1,i) = real(cqtr_norm1(values(:,:,i)))
    end do

    do i=1,spower                                      ! repeated squarings
      call cqtr_square(upto,values)
    end do

    do i=0,upto
      err % sq_norm(0,i) = real(cqtr_normi(values(:,:,i)))
      err % sq_norm(1,i) = real(cqtr_norm1(values(:,:,i)))
    end do
  end if

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (abs(arg(j+1,j)) < tiny(1.0_qp)) then
      values(j+1,j,:) = cmplx(0.0_qp,0.0_qp,qp)
    end if
  end do
end subroutine cqtr_scsqr


!    Name : cqtr_squaref
!   Usage : value = cqtr_squaref(upto,arg)
!         : This function is invoked through the generic name,
!         : ``sqrmtrphif'' or ``sqrqtrphif''.
! Purpose : This function squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, value(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This function treats a complex, quadruple precision, quasi upper triangular matrix.

pure function cqtr_squaref(upto,arg)
  integer,          intent(in) :: upto
  complex(kind=qp), intent(in) :: arg(1:,1:,0:)
  complex(kind=qp)             :: cqtr_squaref(size(arg,1),size(arg,2),0:upto)

  cqtr_squaref = arg
  call cqtr_square(upto,cqtr_squaref)
end function cqtr_squaref


!    Name : cqtr_square
!   Usage : call cqtr_square(upto,arg)
! Purpose : This subroutine squares the phi-functions, i.e.,
!         :   { phi_k(2A); 0 <= k <= upto } is computed from
!         :   { phi_k( A); 0 <= k <= upto }.
!   Input : - ``upto'' specifies the maximum index of phi-functions
!         :   to be squared simultaneously (0 <= upto <= 5).
!         : - The subarray, arg(*,k), stores the value of
!         :   phi_k(A) for 0 <= k <= upto.
!  Output : The subarray, arg(*,k), stores the value of phi_k(2A)
!         : for 0 <= k <= upto.
!    Note : This subroutine does not have a generic name.
!         : This subroutine treats a complex, quadruple precision, quasi upper triangular matrix.

pure subroutine cqtr_square(upto,arg)
  integer,          intent(in   ) :: upto
  complex(kind=qp), intent(inout) :: arg(1:,1:,0:)

  complex(kind=qp) :: Iexp(size(arg,1),size(arg,2))
  integer          :: j, subd(1:size(arg,1))

  subd = 0
  do j=1,size(arg,1)-1
    if (tiny(1.0_qp) <= sum(abs(arg(j+1,j,:)))) subd(j) = 1
  end do

  Iexp = cqtr_eye(size(arg,1))
  Iexp = Iexp + arg(:,:,0)

  if (5 <= upto) then                              ! square phi_5
    arg(:,:,5) = arg(:,:,5) .trtimes. Iexp    /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,4)      /  32.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,3)      /  64.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,2)      / 192.0_qp
    arg(:,:,5) = arg(:,:,5) + arg(:,:,1)      / 768.0_qp
  end if

  if (4 <= upto) then                              ! square phi_4
    arg(:,:,4) = arg(:,:,4) .trtimes. Iexp    /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,3)      /  16.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,2)      /  32.0_qp
    arg(:,:,4) = arg(:,:,4) + arg(:,:,1)      /  96.0_qp
  end if

  if (3 <= upto) then                              ! square phi_3
    arg(:,:,3) = arg(:,:,3) .trtimes. Iexp    /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,2)      /   8.0_qp
    arg(:,:,3) = arg(:,:,3) + arg(:,:,1)      /  16.0_qp
  end if

  if (2 <= upto) then                              ! square phi_2
    arg(:,:,2) = arg(:,:,2) .trtimes. Iexp    /   4.0_qp
    arg(:,:,2) = arg(:,:,2) + arg(:,:,1)      /   4.0_qp
  end if

  if (1 <= upto) then                              ! square phi_1
    arg(:,:,1) = arg(:,:,1) .trtimes. Iexp    /   2.0_qp
  end if
                                                   ! square the exponential
  arg(:,:,0) = arg(:,:,0) .trtimes. arg(:,:,0)

  ! for triangular structure to be invariant

  do j=1,size(arg,1)-1
    if (subd(j) == 0) arg(j+1,j,:) = cmplx(0.0_qp,0.0_qp,qp)
  end do
end subroutine cqtr_square

#endif

end module scalesquare

