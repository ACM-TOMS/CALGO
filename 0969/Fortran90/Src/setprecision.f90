! Module set_precision provides the kind type parameter needed
! to define the precision of a complete package along
! with values for all commonly used precisions
     MODULE set_precision
! ..
! .. Intrinsic Functions ..
      INTRINSIC KIND
! .. Parameters ..
! Define the standard precisions
     INTEGER, PARAMETER :: skind = KIND(0.0E0)
     INTEGER, PARAMETER :: dkind = KIND(0.0D0)
! Set the precision for the whole package
     INTEGER, PARAMETER :: r8 = dkind
! To change the default package precision to single precision change
! the parameter assignment to r8 above to
!    INTEGER, PARAMETER :: r8 = skind
! and recompile the complete package.
!----------------------------------------------------------
! The next statement set quadruple precision.
! This is non-standard and may not be available.
! Different compilers set this value in different ways.
!-----------------------------------------------------------
! For the non-standard quadruple precision:
! IBM and Intel compilers recognize KIND(0.0Q0) 
     INTEGER, PARAMETER:: qkind = KIND(0.0Q0)
! Set the precision for the whole package
     INTEGER, PARAMETER :: r16 = qkind
!-----------------------------------------------------------
! If you are using the NAG compiler then you can
! make use of the NAG-supplied f90_kind module
!     USE, INTRINSIC :: f90_kind
!     INTEGER, PARAMETER:: skind = single ! single precision
!     INTEGER, PARAMETER:: dkind = double ! double precision
!     INTEGER, PARAMETER:: qkind = quad   ! quad   precision
! Then set the precision for the whole package
!     INTEGER, PARAMETER :: r8 = dkind
!     INTEGER, PARAMETER :: r16 = qkind
!-----------------------------------------------------------
    END MODULE set_precision