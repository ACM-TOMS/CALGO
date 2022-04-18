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
      REAL*8     XX(0:MO,0:14)
      REAL*8     v(NVAR)
      REAL*8     p(NPAR)
      REAL*8     t
C-------------------------------------------------------------------------------
      REAL*8     mul_mf
      REAL*8     div_mf
      REAL*8     exp_mf
      REAL*8     pow_mf_c
      REAL*8     sin_mf
      REAL*8     cos_mf
      REAL*8     log_mf
C-------------------------------------------------------------------------------
      TT = 14
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
          XX(i,5) = -0.1d1* XX(i,1)
          XX(i,6) = -0.1d1* XX(i,2)
          XX(i,7) = mul_mf(1,1,i,XX,TT,MO)
          XX(i,8) = mul_mf(2,2,i,XX,TT,MO)
          XX(i,9) = XX(i,7)+XX(i,8)
          XX(i,10) = p(1)* XX(i,5)
          XX(i,11) = p(1)* XX(i,6)
          XX(i,12) = pow_mf_c(9,-0.3d1/0.2d1,12,i,XX,TT,MO)
          XX(i,13) = mul_mf(10,12,i,XX,TT,MO)
          XX(i,14) = mul_mf(11,12,i,XX,TT,MO)
          inext = i + 1
          XX(inext,1) = XX(i,3)/inext
          XX(inext,2) = XX(i,4)/inext
          XX(inext,3) = XX(i,13)/inext
          XX(inext,4) = XX(i,14)/inext
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