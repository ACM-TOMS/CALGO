
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

      BLOCKDATA CONSTMETHOD
      REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
      INTEGER nitermax,nordinc,minord,maxord
      INTEGER accepted_steps, rejected_steps
      LOGICAL dense_output, defect_error_control
      COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
      COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
      COMMON /OPT/ dense_output, defect_error_control
      COMMON /ARS/ accepted_steps, rejected_steps
	  DATA fac1,fac2,fac3/0.95d0,10.d0,0.8d0/
	  DATA rminstep,rmaxstep/1.0d2,1.0d-2/
	  DATA nitermax,nordinc/5,5/
	  DATA minord,maxord/6,26/
      DATA dense_output, defect_error_control/.TRUE.,.FALSE./
      DATA accepted_steps, rejected_steps/0,0/
      END


      SUBROUTINE minf_tides(v,numvar,p,numpar,tini,tend,dt,
     &   tolrel,tolabs)

        IMPLICIT NONE
        LOGICAL dense_output, defect_error_control
        CHARACTER fname*20
        INTEGER accepted_steps, rejected_steps
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        INTEGER FL
        INTEGER numvar,numpar
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /VP/ NVAR, NPAR
        COMMON /OPT/ dense_output, defect_error_control
        COMMON /ARS/ accepted_steps, rejected_steps
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /FILE/FL
        REAL*8 v(numvar),p(numpar),XVAR(0:maxord,0:numvar)
        REAL*8 t0,tini,tend,dt,tolrel,tolabs
        REAL*8 tol,tolo, step, temp, nstep, extra
        INTEGER ORDER
        INTEGER tflag

        TZERO = tini
        t0 = tini
        DELTAT = dt
        TOL_REL = tolrel
        TOL_ABS = tolabs
        NVAR = numvar
        NPAR = numpar
        extra = 0.d0

        IF (t0 .LT. tend) THEN
            tflag = 0
            IF (DELTAT .LT. 0.d0) THEN
                DELTAT = -DELTAT
            END IF
        ELSE
            tflag = 1
            IF (DELTAT .GT. 0.d0) THEN
                DELTAT = -DELTAT
            END IF
        END IF

        IF(dense_output) THEN
       		WRITE(FL,'(90E25.16)') TZERO, v
        END IF
           
        DO WHILE(((t0 .LT. tend).AND.(tflag.EQ.0)) .OR.
     &           ((t0 .GT. tend).AND.(tflag.EQ.1)))
            CALL tolerances_mf(v, tol,tolo, ORDER)
            CALL minfseries(t0,v,NVAR,p,NPAR,XVAR,ORDER,maxord)
            CALL steps_mf(tol, XVAR, ORDER, step, tflag)

            IF (defect_error_control) THEN
                CALL steps_DEC_mf(t0,step,tolo,p,XVAR,ORDER)
            END IF
 
            temp = t0
            nstep = step + extra
            t0 = temp + nstep
            extra = (temp-t0)+nstep
            IF(((t0 .GT. tend).AND.(tflag.EQ.0)) .OR.
     &           ((t0 .LT. tend).AND.(tflag.EQ.1))) THEN
                nstep  = (tend-temp)
            END IF
            accepted_steps =  accepted_steps + 1
            IF(dense_output) THEN
                CALL dense_output_mf(temp,nstep,XVAR,ORDER,tflag)
            END IF
            CALL horner_mf(v,XVAR,ORDER,nstep)
                
        END DO
           


        RETURN
      END SUBROUTINE
           
C--------------------------------------------------------------
C--------------------------------------------------------------
C     tolerances_mf
C--------------------------------------------------------------
C--------------------------------------------------------------

      SUBROUTINE tolerances_mf(v, tol, tolo, ORDER)       
        IMPLICIT NONE
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        INTEGER ORDER
        REAL*8 yna,ynb,tol,tolo,v(NVAR),miny
           DATA ynb /0.0d0/
           SAVE ynb
          
        CALL norm_inf_vec_mf(v, yna)
        tol = TOL_ABS + MAX(yna,ynb)*TOL_REL
        miny = MIN(yna,ynb)
        IF(miny .gt. 0.d0) THEN
             tolo = MIN(TOL_ABS/miny, TOL_REL)
        ELSE
             tolo = MIN(TOL_ABS, TOL_REL)
        END IF


        ORDER = MIN(maxord, int(-log(tolo)/2)+nordinc)
        ORDER = MAX(minord,ORDER)
        
        ynb = yna
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     steps_mf
C--------------------------------------------------------------
C--------------------------------------------------------------


      SUBROUTINE steps_mf(tol, XVAR, ORDER, step, tflag)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 XVAR(0:maxord,0:NVAR)
        REAL*8 tol, ynu,ynp,dord,dorda,dordp,sp,su
        INTEGER ORDER, ord, orda, ordp
        REAL*8 step, stepant, rstep
        DATA stepant /0.0d0/
        SAVE stepant
        INTEGER tflag
          
        ord = ORDER+1
        ynu = 0.d0
        DO WHILE((ynu .EQ. 0.d0) .AND. (ord .GT. 0))
           ord = ord-1
           CALL norm_inf_mat_mf(XVAR,ord, ynu)
        END DO
         
        IF(ord .NE. 0)  THEN
             orda = ORDER-1
             ordp = ORDER+1
             dord  = 1.d0/DFLOAT(ord)
             dorda = 1.d0/DFLOAT(orda)
             dordp = 1.d0/DFLOAT(ordp)
             CALL norm_inf_mat_mf(XVAR,orda, ynp)
             IF(ynp .eq. 0.d0) THEN
               step = (tol**dordp) * ((1.d00/ynu)**dord)
             ELSE
              sp = (tol**dord) * ((1.d00/ynp)**dorda)
              su = (tol**dordp) * ((1.d00/ynu)**dord)
              step = MIN(sp,su)
            END IF 
            IF(stepant. NE. 0.0d0) THEN
              rstep = step/stepant
              IF(rstep .GT. rmaxstep) THEN
                 step = rmaxstep*stepant
              ELSE IF(rstep .LT. rminstep) THEN
                 step = rminstep*stepant
              END IF
            END IF
            step = fac1*step
        ELSE
            WRITE(*,*) '*********Error*********'
            STOP
        END IF
        IF (tflag .EQ. 1) THEN
            step = -step
        END IF
        RETURN
      END SUBROUTINE

 

C--------------------------------------------------------------
C--------------------------------------------------------------
C     steps_DEC_mf
C--------------------------------------------------------------
C--------------------------------------------------------------


      SUBROUTINE steps_DEC_mf(t0,step,tolo,p,XVAR,ORDER)      
        IMPLICIT NONE
        INTEGER accepted_steps, rejected_steps
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /VP/ NVAR, NPAR
        COMMON /ARS/ accepted_steps, rejected_steps
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 vh(NVAR),vdh(NVAR), NORM
        REAL*8 t0,step,tolo, p(NPAR)
        REAL*8 XVAR(0:maxord,0:NVAR), XVAR2(0:maxord,0:NVAR)
        INTEGER ITER, ORDER,I
c
        ITER = 1
        NORM = 1.d99
           
        DO WHILE ((NORM .GT. fac2*tolo) .AND. (ITER .LT. nitermax))
             
            CALL horner_mf(vh,XVAR,ORDER,step)
            CALL hornerd_mf(vdh,XVAR,ORDER,step)
            CALL minfseries(t0+step,vh,NVAR,p,NPAR,XVAR2,1,1)
            NORM = ABS(XVAR2(1,1)-vdh(1))
            DO I=2, NVAR
                    NORM = NORM + ABS(XVAR2(1,i)-vdh(i))
            END DO
            IF(ITER .GT.1) THEN
                rejected_steps = rejected_steps+1
                step = fac3*step
            END IF
            ITER = ITER+1
               
        END DO
c
        RETURN
      END SUBROUTINE

         
C--------------------------------------------------------------
C--------------------------------------------------------------
C     norm_inf_vec_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE norm_inf_vec_mf(v, norm)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        COMMON /VP/ NVAR, NPAR
        REAL*8 norm
        REAL*8 v(NVAR)
        INTEGER I
        norm = 0.d0
        DO I =1,NVAR
             norm =MAX(ABS(v(I)), norm)
        END DO
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     norm_inf_mat_mf
C--------------------------------------------------------------
C--------------------------------------------------------------

      SUBROUTINE norm_inf_mat_mf(XVAR,ord, norm)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        COMMON /VP/ NVAR, NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 norm
        REAL*8 XVAR(0:maxord,0:NVAR)
        INTEGER I,ord
        norm = 0.e0
        DO I =1,NVAR
             norm =MAX(abs(XVAR(ord,I)), norm)
        END DO
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     dense_output_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE dense_output_mf(t0,step,XVAR,ORDER,tflag)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER FL
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /VP/ NVAR, NPAR
        COMMON /FILE/ FL
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 v(NVAR), XVAR(0:maxord,0:NVAR)
        REAL*8 t0,ti,tit,tend,step
        INTEGER  IPOS, ORDER
        DATA IPOS /1/
        SAVE IPOS
        INTEGER tflag

        tend = t0 + step
        ti = TZERO + IPOS*DELTAT
        tit= ti    - t0
        DO WHILE(((tit .LE. step).AND.(tflag.EQ.0)) .OR.
     &           ((tit .GE. step).AND.(tflag.EQ.1)))
            CALL horner_mf(v,XVAR,ORDER,tit)
            WRITE(FL,'(90E25.16)') ti, v
            IPOS =IPOS + 1
            ti = TZERO + IPOS*DELTAT
            tit= ti    - t0
        END DO
        RETURN
      END SUBROUTINE
C--------------------------------------------------------------
C--------------------------------------------------------------
C     horner_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE horner_mf(v,XVAR,ORDER,t)       
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        INTEGER ORDER, i, nd
        REAL*8 XVAR(0:maxord,0:NVAR), t, au
        REAL*8 v(NVAR)
        DO nd=1, NVAR
          au = XVAR(ORDER,nd)*t
          DO i = ORDER-1, 1, -1
            au = (au+XVAR(i,nd))*t
          END DO
          v(nd) = au+XVAR(0,nd)
        END DO
        RETURN
      END SUBROUTINE

C--------------------------------------------------------------
C--------------------------------------------------------------
C     hornerd_mf 
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE hornerd_mf(v,XVAR,ORDER,t)     
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        INTEGER ORDER, i, nd
        REAL*8 XVAR(0:maxord,0:NVAR), t, au
        REAL*8 v(NVAR)
        DO nd=1, NVAR
          au = ORDER*XVAR(ORDER,nd)*t
          DO i = ORDER-1, 2, -1
            au = (au+i*XVAR(i,nd))*t
          END DO
          v(nd) = au+XVAR(1,nd)
        END DO
        RETURN
       END SUBROUTINE

C--------------------------------------------------------------
C--------------------------------------------------------------
C     mul_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION mul_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 mul_mf
        au = XX(0,n1)*XX(i,n2)
        DO m = 1, i
          au = au + XX(m,n1)*XX(i-m,n2)
        END DO
        mul_mf = au
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     div_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION div_mf(n1,n2,n3,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,n3,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 div_mf
        IF(XX(0,n2) .EQ. 0.d0)THEN
          WRITE(*,*) ' Function div_mf found division by zero'
          STOP
        ELSE
          au = XX(i,n1)
          DO m = 1, i
            au = au - XX(m,n2)*XX(i-m,n3)
          END DO
          div_mf = au/XX(0,n2)
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     pow_mf_c
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION pow_mf_c(n1,ex,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 ex,au
        REAL*8 pow_mf_c
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .EQ. 0.d0) THEN
            WRITE(*,*) ' Function pow_mf_c found division by zero'
            STOP
          ELSE
            pow_mf_c = XX(0,n1)**ex
          END IF
        ELSE
          IF(XX(0,n1) .NE. 0.d0) THEN
            au = ex*i*XX(0,n2)*XX(i,n1)
            DO m = 1, i - 1
              au = au+(ex*(i-m)-m)*XX(m,n2)*XX(i-m,n1)
            END DO
            pow_mf_c = au/(i*XX(0,n1))
          ELSE
            pow_mf_c = 0.D0
          END IF
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     exp_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION exp_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 exp_mf
        IF(i.EQ.0)THEN
          exp_mf = EXP(XX(0,n1))
        ELSE
          au = i*XX(0,n2)*XX(i,n1)
          DO m = 1, i - 1
            au = au+(i-m)*XX(m,n2)*XX(i-m,n1)
          END DO
          exp_mf = au/i
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     log_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION log_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 log_mf
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .LE. 0.d0)THEN
            WRITE(*,*) 'Function log_mf found log of a negative value'
            STOP
          ELSE
            log_mf = LOG(XX(0,n1))
          END IF
        ELSE
          au = i*XX(i,n1)
          DO m = 1, i - 1
            au = au-(i-m)*XX(m,n1)*XX(i-m,n2)
          END DO
          log_mf = au/i/XX(0,n1)
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     sin_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION sin_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 sin_mf
        IF(i.EQ.0)THEN
          sin_mf = SIN(XX(0,n1))
        ELSE
          au = XX(1,n1)*XX(i-1,n2)
          DO m = 2, i
            au = au+m*XX(m,n1)*XX(i-m,n2)
          END DO
          sin_mf = au/i
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     cos_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION cos_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 cos_mf
        IF(i.EQ.0)THEN
          cos_mf = COS(XX(0,n1))
        ELSE
          au = XX(1,n1)*XX(i-1,n2)
          DO m = 2, i
            au = au+m*XX(m,n1)*XX(i-m,n2)
          END DO
          cos_mf = -au/i
        END IF
        RETURN
      END




