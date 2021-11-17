! This is a portable version of D1MACH that should work for any compiler
! that complies with the Fortran 90 specifications for the intrinsic
! EPSILON.

      DOUBLE PRECISION FUNCTION D1MACH(I) 
      INTEGER I, I1MACH
!                                                                       
!  DOUBLE-PRECISION MACHINE CONSTANTS                                   
!                                                                       
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.           
!                                                                       
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.           
!                                                                       
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.                 
!                                                                       
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.                 
!                                                                       
!  D1MACH( 5) = LOG10(B)                                                
!                                                                       
      SELECT CASE (I) 
        CASE (1) 
           D1MACH = TINY(1D0) 
        CASE (2) 
           D1MACH = HUGE(1D0) 
        CASE (3) 
           D1MACH = EPSILON(1D0)/RADIX(1D0) 
        CASE (4) 
           D1MACH = EPSILON(1D0) 
        CASE (5) 
           D1MACH = RADIX(1D0) 
           D1MACH = LOG10(D1MACH) 
        CASE DEFAULT 
        WRITE(I1MACH(2),'(A,I10)') ' D1MACH - I OUT OF BOUNDS', I 
        STOP 
      END SELECT 
      END                                           

