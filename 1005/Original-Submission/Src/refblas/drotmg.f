*> \brief \b DROTMG
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DD1,DD2,DX1,DY1
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DPARAM(5)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
*>    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*>    DY2)**T.
*>    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*>
*>    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
*>
*>      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
*>    H=(          )    (          )    (          )    (          )
*>      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
*>    LOCATIONS 2-4 OF DPARAM CONTAIN DH11,DH21,DH12, AND DH22
*>    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
*>    VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
*>
*>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
*>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
*>    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in,out] DD1
*> \verbatim
*>          DD1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in,out] DD2
*> \verbatim
*>          DD2 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in,out] DX1
*> \verbatim
*>          DX1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] DY1
*> \verbatim
*>          DY1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] DPARAM
*> \verbatim
*>          DPARAM is DOUBLE PRECISION array, dimension (5)
*>     DPARAM(1)=DFLAG
*>     DPARAM(2)=DH11
*>     DPARAM(3)=DH21
*>     DPARAM(4)=DH12
*>     DPARAM(5)=DH22
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup double_blas_level1
*
*  =====================================================================
      SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
*
*  -- Reference BLAS level1 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DD1,DD2,DX1,DY1
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DPARAM(5)
*     ..
*
* =====================================================================
* Corrections by Kristján Jónasson 2018. The subroutines SROTMG and
* DROTMG in the reference BLAS distributed with LAPACK versions 3.2 -3.8
* contain an error, which the original (1979) routines in the reference
* BLAS level 1 (TOMS Algorithm 539) do not have. The original
* subroutines use goto to implement three internal procedures,
* ZERO-H-D-AND-DX1, SCALE-CHECK and FIX-CHECK, assigned goto for the
* last one, The subroutines were translated to be consistent with the
* Fortran 90 standard somtime around the year 2000. When refactoring to
* make the translations goto-less a bug crept in. The bug is filed as
* issue #244 where the current reference BLAS is kept in the github
* repository. The corrections implemented below are commented with
* KJ-2018. They are based on (hopefully careful) compariosn with the
* original BLAS, available on calgo.acm.org.
* =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP,
     $                 DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
*     .. Data statements ..
*
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
*     ..

      IF (DD1.LT.ZERO) THEN
*        PROCEDURE ZERO-H-D-AND-DX1
         DFLAG = -ONE
         DH11 = ZERO
         DH12 = ZERO
         DH21 = ZERO
         DH22 = ZERO
*         
         DD1 = ZERO
         DD2 = ZERO
         DX1 = ZERO
      ELSE
*        CASE-DD1-NONNEGATIVE
         DP2 = DD2*DY1
         IF (DP2.EQ.ZERO) THEN
            DFLAG = -TWO
            DPARAM(1) = DFLAG
            RETURN
         END IF
*        REGULAR-CASE..
         DP1 = DD1*DX1
         DQ2 = DP2*DY1
         DQ1 = DP1*DX1
*
         IF (DABS(DQ1).GT.DABS(DQ2)) THEN
            DH21 = -DY1/DX1
            DH12 = DP2/DP1
*
            DU = ONE - DH12*DH21
*
           IF (DU.GT.ZERO) THEN
             DFLAG = ZERO
             DD1 = DD1/DU
             DD2 = DD2/DU
             DX1 = DX1*DU
          ELSE
*            PROCEDURE ZERO-H-D-AND-DX1
*            KJ-2018. Actually the following section should never be             
*            reached, because with exact arithmetic the preceding if
*            statements ensure that DU > 0, but, nevertheless, the 
*            original from 1979 contains this else clause logically.
             DFLAG = -ONE
             DH11 = ZERO
             DH12 = ZERO
             DH21 = ZERO
             DH22 = ZERO
             DD1 = ZERO
             DD2 = ZERO
             DX1 = ZERO
*            -------------
           END IF
         ELSE

            IF (DQ2.LT.ZERO) THEN
*              PROCEDURE ZERO-H-D-AND-DX1
               DFLAG = -ONE
               DH11 = ZERO
               DH12 = ZERO
               DH21 = ZERO
               DH22 = ZERO
               DD1 = ZERO
               DD2 = ZERO
               DX1 = ZERO
*              ----------
            ELSE
               DFLAG = ONE
               DH11 = DP1/DP2
               DH22 = DX1/DY1
               DU = ONE + DH11*DH22
               DTEMP = DD2/DU
               DD2 = DD1/DU
               DD1 = DTEMP
               DX1 = DY1*DU
            END IF
         END IF
                  
*        KJ-2018: The main refactoring starts here and ends at location (*)         
*        PROCEDURE SCALE-CHECK
         DO WHILE (DD1 .LE. RGAMSQ)
            IF (DD1.EQ.ZERO) EXIT
*           PROCEDURE FIX-H
            IF (DFLAG .EQ. ZERO) THEN
               DH11 = ONE
               DH22 = ONE
            ELSE IF (DFLAG .GT. ZERO) THEN
               DH21 = -ONE
               DH12 = ONE
            ENDIF
            DFLAG = -ONE
*           ------------------
            DD1 = DD1*GAM**2
            DX1 = DX1/GAM
            DH11 = DH11/GAM
            DH12 = DH12/GAM
         ENDDO
         DO WHILE (DD1 .GE. GAMSQ)
*           PROCEDURE FIX-H
            IF (DFLAG .EQ. ZERO) THEN
               DH11 = ONE
               DH22 = ONE
            ELSE IF (DFLAG .GT. ZERO) THEN
               DH21 = -ONE
               DH12 = ONE
            ENDIF
            DFLAG = -ONE
*           --------------------
            DD1 = DD1/GAM**2
            DX1 = DX1*GAM
            DH11 = DH11*GAM
            DH12 = DH12*GAM
         ENDDO
         DO WHILE (DABS(DD2).LE.RGAMSQ)
            IF (DABS(DD2).EQ.ZERO) EXIT
*           PROCEDURE FIX-H
            IF (DFLAG .EQ. ZERO) THEN
               DH11 = ONE
               DH22 = ONE
            ELSE IF (DFLAG .GT. ZERO) THEN
               DH21 = -ONE
               DH12 = ONE
            ENDIF
            DFLAG = -ONE
*           ---------------------------
            DD2 = DD2*GAM**2
            DH21 = DH21/GAM
            DH22 = DH22/GAM
         ENDDO
         DO WHILE (DABS(DD2).GE.GAMSQ)
*           PROCEDURE FIX-H
            IF (DFLAG .EQ. ZERO) THEN
               DH11 = ONE
               DH22 = ONE
            ELSE IF (DFLAG .GT. ZERO) THEN
               DH21 = -ONE
               DH12 = ONE
            ENDIF
            DFLAG = -ONE
*           ---------------------------
            DD2 = DD2/GAM**2
            DH21 = DH21*GAM
            DH22 = DH22*GAM
         END DO
      END IF
*     KJ-2018. Location (*), end of main refactoring section 
      IF (DFLAG.LT.ZERO) THEN
         DPARAM(2) = DH11
         DPARAM(3) = DH21
         DPARAM(4) = DH12
         DPARAM(5) = DH22
      ELSE IF (DFLAG.EQ.ZERO) THEN
         DPARAM(3) = DH21
         DPARAM(4) = DH12
      ELSE
         DPARAM(2) = DH11
         DPARAM(5) = DH22
      END IF

      DPARAM(1) = DFLAG
      RETURN
      END
