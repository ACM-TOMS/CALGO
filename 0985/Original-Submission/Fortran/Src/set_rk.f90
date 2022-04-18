! Module set_rk provides the kind type parameter needed to define
! defines the precision of a complete package along
! with values for commonly used standard precisions
module set_rk
    public :: r4,r8,rk

    ! .. Intrinsic Functions ..
    intrinsic kind
    ! .. Parameters ..
    ! Define the standard precisions
    ! For IEEE standard arithmetic we could also use
    !     integer, parameter :: r4 = SELECTED_REAL_KIND(p=6, r=37)
    !     integer, parameter :: r8 = SELECTED_REAL_KIND(p=15, r=307)

    integer, parameter :: r4 = kind(0.0e0)  ! single
    integer, parameter :: r8 = kind(0.0d0)  ! double

    ! Set the precision for the whole package
    integer, parameter :: rk=r8

    ! To change the default package precision to single precision change
    ! the parameter assignment to rk above to
    !     integer, parameter :: rk = r4
    ! and recompile the complete package.

end module set_rk
