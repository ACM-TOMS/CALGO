! Log and error information for this library
!

module mtrcfgphilog

  implicit none

  ! A type of error and log information for scaling and squaring methods.

  type mcpsqrlog

    ! The meanings of the ``error'' entry is
    !  0 : No error is detected,
    !  1 : The argument ``upto'' is out of range,
    !  2 : The matrix is not square,
    !  3 : The matrix is not quasi upper triangular.
    !  4 : The size of the matrix is zero.
    !  5 : The block size ``mpt_blksize'' is less then 4.

    integer :: error            ! error code

    integer :: rorder           ! m, the order of rational approximant
    integer :: spower           ! s, the number of squarings in the computation
    integer :: upto             ! the maximal index of phi-function

    ! The following entries store the matrix norm.
    ! (0,:) stores the infinity-norm and (1,:) stores the one-norm.
    ! The second index n of zp_norm(:,n) means that the entry is for |Z^n|,
    ! where Z := hci*A/2^s. The second index of (pm|sp|sq)_norm is the
    ! index of phi-function.

    real :: ah_norm(0:1,0:0)    ! the norm of the matrix hci*A
    real :: zp_norm(0:1,0:8)    ! the norm of the power of the scaled matrix
    real :: qm_norm(0:1,0:0)    ! the denominator of the approximants
    real :: pm_norm(0:1,0:5)    ! the numerators  of the approximants
    real :: sp_norm(0:1,0:5)    ! before squaring, phi_n(hci*A/2^s)
    real :: sq_norm(0:1,0:5)    ! after  squaring, phi_n(hci*A)

  end type mcpsqrlog

  public

contains

! Print the contents of type(mcpsqrlog)

subroutine print_mcpsqrlog(unit,info)
  integer,         intent(in) :: unit
  type(mcpsqrlog), intent(in) :: info

  integer :: i

  write(unit,'(a)'       ) ''
  write(unit,'(a)'       ) '# Status :'
  write(unit,'(a,i3)'    ) 'error         = ', info % error
  write(unit,'(a,i3)'    ) 'order, m,     = ', info % rorder
  write(unit,'(a,i3)'    ) 'scale, s,     = ', info % spower
  write(unit,'(a,i3)'    ) 'upto          = ', info % upto

  write(unit,'(a)'       ) ''
  write(unit,'(a)'       ) '# Argument hci*A :'
  write(unit,'(a,1f10.5)') '|hci*A|_inf   = ', info % ah_norm(0,0)
  write(unit,'(a,1f10.5)') '|hci*A|_1     = ', info % ah_norm(1,0)

  write(unit,'(a)'       ) ''
  write(unit,'(a)'       ) '# After multiplications, (n=2,4,6,8) :'
  write(unit,'(a,4f10.5)') '|A/2^s|^n_inf = ', (info % zp_norm(0,i), i=2,8,2)
  write(unit,'(a,4f10.5)') '|A/2^s|^n_1   = ', (info % zp_norm(1,i), i=2,8,2)

  write(unit,'(a)'       ) ''
  write(unit,'(a,2f10.5)') '# Denominator and numerator of the approximant :'
  write(unit,'(a,2f10.5)') '|qm|_inf,1    = ', (info % qm_norm(i,0), i=0,1,1)
  write(unit,'(a,6f10.5)') '|pm|_inf      = ', (info % pm_norm(0,i), i=0,5,1)
  write(unit,'(a,6f10.5)') '|pm|_1        = ', (info % pm_norm(1,i), i=0,5,1)

  write(unit,'(a)'       ) ''
  write(unit,'(a)'       ) '# Before squaring, (0 <= n <= 5) :'
  write(unit,'(a,6f10.5)') '|phi_n|_inf   = ', (info % sp_norm(0,i), i=0,5,1)
  write(unit,'(a,6f10.5)') '|phi_n|_1     = ', (info % sp_norm(1,i), i=0,5,1)

  write(unit,'(a)'       ) ''
  write(unit,'(a)'       ) '# After  squaring, (0 <= n <= 5) :'
  write(unit,'(a,6f10.5)') '|phi_n|_inf   = ', (info % sq_norm(0,i), i=0,5,1)
  write(unit,'(a,6f10.5)') '|phi_n|_1     = ', (info % sq_norm(1,i), i=0,5,1)

  write(unit,*) ''
end subroutine print_mcpsqrlog

! Initialize the contents of type(mcpsqrlog)

pure subroutine init_mcpsqrlog(info)
  type(mcpsqrlog), intent(out) :: info

  info % error   = 0
  info % rorder  = 0
  info % spower  = 0
  info % zp_norm = 0
  info % qm_norm = 0
  info % pm_norm = 0
  info % sp_norm = 0
  info % sq_norm = 0
end subroutine init_mcpsqrlog

end module mtrcfgphilog

