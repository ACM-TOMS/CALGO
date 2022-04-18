      SUBROUTINE qtql1(n,d,e,ierr)
************************************************************************
*     (Quadruple-precision symmetric tridiagonal matrix eigenvalues)
*     Given a symmetric tridiagonal matrix of order n stored with its
*     diagonal in d(1..n), and its subdiagonal in e(2..n), with e(1)
*     arbitrary, use the QL method to find the eigenvalues.
*
*     On return, the original contents of d(*) and e(*) will have been
*     destroyed, and d(*) will contain the eigenvalues in ascending
*     order.
*
*     The error indicator, ierr, is normally set to zero on return.
*     However, if an error exit is made, ierr is set to a positive
*     value, and then the eigenvalues are correct and ordered for
*     indices 1, 2, ..., ierr-1, but may not be the smallest
*     eigenvalues.
*
*     This routine is a translation of the Algol procedure tql1() from
*     ``The QR and QL Algorithms for Symmetric Matrices'', Numerische
*     Mathematik 11, 293--306 (1968), by H. J. Bowdler, R. S. Martin,
*     C. Reinsch and J.  H. Wilkinson.  That article was republished
*     in the Handbook for Automatic Computation, Vol. II, Linear
*     Algebra, 227--240 (1971), eds. J. H. Wilkinson and C. Reinsch,
*     Springer-Verlag, ISBN 3-540-05414-6.
*
*     tql1() is part of the EISPACK-1 and EISPACK-2 libraries.
*     (31-Aug-1983)
************************************************************************
*
*     External functions
*
      EXTERNAL            q2norm
*
      REAL*16             q2norm,      qabs,        qsign
*
*     Argument variables
*
      INTEGER             ierr,        n
*
      REAL*16             d(*),        e(*)
*
*     Local variables
*
      INTEGER             i,           ii,          j,           l
      INTEGER             l1,          l2,          m,           mml
*
      REAL*16             c,           c2,          c3,          dl1
      REAL*16             el1,         f,           g,           h
      REAL*16             p,           r,           s,           s2
      REAL*16             tst1,        tst2
*
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL1,
*     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
*     WILKINSON.
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
*
*     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
*     TRIDIAGONAL MATRIX BY THE QL METHOD.
*
*     ON INPUT
*
*        N IS THE ORDER OF THE MATRIX.
*
*        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
*
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
*          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
*
*      ON OUTPUT
*
*        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
*          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
*          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
*          THE SMALLEST EIGENVALUES.
*
*        E HAS BEEN DESTROYED.
*
*        IERR IS SET TO
*          ZERO       FOR NORMAL RETURN,
*          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
*                     DETERMINED AFTER 30 ITERATIONS.
*
*     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
*
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
*     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
*
*     THIS VERSION DATED AUGUST 1983.
*
*     ------------------------------------------------------------------
*
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
*
      DO 100 I = 2, N
          E(I-1) = E(I)
  100 CONTINUE
*
      F = 0.0Q0
      TST1 = 0.0Q0
      E(N) = 0.0Q0
*
      DO 290 L = 1, N
         J = 0
         H = QABS(D(L)) + QABS(E(L))
         IF (TST1 .LT. H) TST1 = H
*     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + QABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
*     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
*                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
*
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
*     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0Q0 * E(L))
         R = Q2NORM(P,1.0Q0)
         D(L) = E(L) / (P + QSIGN(R,P))
         D(L1) = E(L) * (P + QSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
*
         DO 140 I = L2, N
             D(I) = D(I) - H
  140    CONTINUE
*
  145    F = F + H
*     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0Q0
         C2 = C
         EL1 = E(L1)
         S = 0.0Q0
         MML = M - L
*     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = Q2NORM(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
  200    CONTINUE
*
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + QABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  210    P = D(L) + F
*     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
*     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
*
  250    I = 1
  270    D(I) = P
  290 CONTINUE
*
      GO TO 1001
*     .......... SET ERROR -- NO CONVERGENCE TO AN
*                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
