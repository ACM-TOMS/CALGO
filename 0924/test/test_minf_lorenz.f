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


      Program  dr_minf_lorenz
      IMPLICIT NONE
      INTEGER  i,j
C --- NUMBER OF VARIABLES AND PARAMETERS
      INTEGER  NVAR,NPAR
      PARAMETER  (NVAR = 3)
      PARAMETER  (NPAR = 3)
C --- TOLERANCES
      REAL*8 tolabs,tolrel
      REAL*8 aux, error
C --- TIMES: INITIAL, FINAL, INCREMENT
      REAL*8 tini, tend, dt
C --- VARIABLES AND PARAMETERS
      REAL*8 v(NVAR)
      REAL*8 p(NPAR)
      REAL*8 x(NVAR)
C --- FILE NAME AND UNIT NUMBER OF DENSE OUTPUT
      CHARACTER fname*20
      INTEGER   FL
C --- OPTIONS
      LOGICAL dense_output, defect_error_control
C --- COUNTERS
      INTEGER accepted_steps, rejected_steps
C --- CONSTANTS OF THE METHOD (safety factors, maximum order, ...)
      REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
      INTEGER nitermax,nordinc,minord,maxord
C --- GLOBALS
      COMMON /OPT/ dense_output, defect_error_control
      COMMON /ARS/ accepted_steps, rejected_steps
      COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
      COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
      COMMON /FILE/ FL




C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C     INITIAL CONDITIONS,  INTEGRATION TIMES, TOLERANCES
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

C --- PARAMETERS VALUE
      p(1) = 10.d0
      p(2) = 28.d0
      p(3) = 8.d0/3.0d0

C --- INITIAL VALUES
      v(1) = -13.7636106821342005d0
      v(2) = -19.5787519424517955d0
      v(3) = 27.d0

C --- INITIAL INTEGRATION POINT
      tini = 0.d0

C --- ENDPOINT OF INTEGRATION
      tend = 1.558652210716175d0

C --- DELTA t FOR DENSE OUTPUT
      dt   = 1.558652210716175d0


C --- REQUIRED TOLERANCES
      tolrel = 1.d-14
      tolabs = 1.d-14

C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C       DENSE OUTPUT (file , screen or none)
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

      dense_output = .FALSE.

C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C       SHOW INITIAL POINT ON THE SCREEN
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

C      WRITE (*,91) tini,(v(i),i=1,NVAR)

C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C       CALL THE INTEGRATOR
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

      DO 10 i = 1,NVAR
     	x(i) = v(i)
 10   CONTINUE

      CALL minf_tides(v,NVAR,p,NPAR,tini,tend,dt,
     &   tolrel,tolabs)

C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C       SHOW FINAL POINT ON THE SCREEN
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

C      WRITE (*,91) tend,(v(i),i=1,NVAR)

C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C       CHECK PERIODIC ORBIT
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

      aux = 0.d0
      error = -1.d0

      DO 20 i = 1,NVAR
      	aux = ABS(v(i)-x(i))
      	IF (aux .gt. error) error = aux
 20   CONTINUE

      IF (error .lt. 1.d-10) THEN
      	CALL exit(0) 
      	STOP
      ELSE
      	write (*,*) "ERROR IN THE ORBIT = ", error
      	CALL exit(1) 
      	STOP
      ENDIF

91    FORMAT(1X,'t =',E25.16,'    X =',90E25.16)

      STOP
      END




