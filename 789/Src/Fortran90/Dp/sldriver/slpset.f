C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
      module SLPSET
C The variables which define the current 'settings' of the SLP being
C solved. Put in a module to give more convenient access to the few
C routines that need to access these from 'outside'.
C Done for v4.0. More to add at next revision which is to be more OO!
      double precision A,A1,A2,B,B1,B2,
     +  AORIG,A1ORIG,A2ORIG,BORIG,B1ORIG,B2ORIG

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

        subroutine GTABCU(ACURR,BCURR)
C Access routine for current SLP endpoints
        double precision ACURR,BCURR
        ACURR = A
        BCURR = B
        end subroutine GTABCU

      end module SLPSET
