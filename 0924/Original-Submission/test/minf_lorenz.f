C****************************************************************************
C       
C     This file is part of TIDES.
C       
C     Contributors:
C     
C     A. Abad, R. Barrio, F. Blesa, M. Rodriguez
C     Grupo de Mecanica Espacial
C     University of Zaragoza
C     SPAIN
C       
C     http://gme.unizar.es/software/tides
C     Contact: <tides@unizar.es>
C        
C        
C*****************************************************************************


      SUBROUTINE minfseries(t,v,NVAR,p,NPAR,XVAR,ORDER,MO)
      IMPLICIT   NONE
      INTEGER    i,inext,j,ORDER,NVAR,NPAR,TT,MO
      REAL*8     XVAR(0:MO,0:NVAR)
      REAL*8     XX(0:MO,0:13)
      REAL*8     v(NVAR)
      REAL*8     p(NPAR)
      REAL*8     t
      REAL*8     pr(4)
C-------------------------------------------------------------------------------
      REAL*8     mul_mf
      REAL*8     div_mf
      REAL*8     exp_mf
      REAL*8     pow_mf_c
      REAL*8     sin_mf
      REAL*8     cos_mf
      REAL*8     log_mf
C-------------------------------------------------------------------------------
      DO j = 1, NPAR
          pr(j) = p(j)
      END DO
      pr(4) = -0.1d1 * pr(3)
      TT = 13
      DO j=0, TT
          DO i=0, ORDER
              XX(i,j) = 0.d0
          END DO
      END DO
      XX(0,0) = t
      XX(1,0) = 1.d0
      DO j = 1, NVAR
          XX(0,j) = v(j)
      END DO
C-------------------------------------------------------------------------------
      DO i=0, ORDER-1
          XX(i,4) = XX(i,2)-XX(i,1)
          XX(i,5) = -0.1d1* XX(i,1)
          XX(i,6) = pr(2)* XX(i,1)
          XX(i,7) = mul_mf(1,2,i,XX,TT,MO)
          XX(i,8) = XX(i,6)-XX(i,2)
          XX(i,9) = pr(1)* XX(i,4)
          XX(i,10) = pr(4)* XX(i,3)
          XX(i,11) = mul_mf(5,3,i,XX,TT,MO)
          XX(i,12) = XX(i,7)+XX(i,10)
          XX(i,13) = XX(i,8)+XX(i,11)
          inext = i + 1
          XX(inext,1) = XX(i,9)/inext
          XX(inext,2) = XX(i,13)/inext
          XX(inext,3) = XX(i,12)/inext
      END DO
C-------------------------------------------------------------------------------
      DO j=0, NVAR
          DO i=0, ORDER
              XVAR(i,j) = XX(i,j)
          END DO
      END DO
      RETURN
      END SUBROUTINE
C-------------------------------------------------------------------------------

