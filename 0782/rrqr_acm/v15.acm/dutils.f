*
*     This code is part of a release of the package for computing
*     rank-revealing QR Factorizations written by:
*     ==================================================================
*     Christian H. Bischof        and   Gregorio Quintana-Orti
*     Math. and Comp. Sci. Div.         Departamento de Informatica
*     Argonne National Lab.             Universidad Jaime I
*     Argonne, IL 60439                 Campus P. Roja, 12071 Castellon
*     USA                               Spain
*     bischof@mcs.anl.gov               gquintan@inf.uji.es
*     ==================================================================
*     $Revision: 1.84 $
*     $Date: 96/12/30 16:59:22 $
*
*********************************************************************
      INTEGER FUNCTION flXGEQPF( m, n )
      INTEGER m, n
*     returns flop count for xgeqpf (lapack routine with column
*     pivoting)
      INTEGER tflops, i
*              initialize column norms
      tflops = 3*n+1
      DO 10 i =1,min(m,n)
*              find pivot column and update partial column norms
         tflops = tflops + 10*(n-i)
*              compute HH vector
         tflops = tflops + 3*(m-i+1)+6
*              update remaining submatrix
         tflops = tflops + 4*(m-i+1)*(n-i)+3*(n-i)
10    CONTINUE
      flXGEQPF = tflops
      RETURN
      END
*********************************************************************
      INTEGER FUNCTION flXGEQR2(m,n)
      INTEGER m, n
*     returns flop count for xgeqr2 (lapack blas 2 routine for
*     QR factorization without column exchanges)
      INTEGER i, tflops
      tflops = 0
      DO 10 i =1,min(m,n)
*              compute HH vector
         tflops = tflops + 3*(m-i+1)+6
*              update remaining submatrix
         tflops = tflops + 4*(m-i+1)*(n-i)+3*(n-i)
10    CONTINUE
      flXGEQR2 = tflops
      RETURN
      END
*********************************************************************
      INTEGER FUNCTION flXLARFT(m,nb)
      INTEGER m, nb
*     returns flop count for slarft (generation of a block
*     reflector)
      INTEGER i, tflops
      tflops = 0
      DO 10 i = 2, nb
*              flops for DGEMV
         tflops = tflops + 2*(m-i+1)*(i-1)+2*(i-1)
*              flops for DTRMV
         tflops = tflops + (i-1)*(i-1)+2*(i-1)
10    CONTINUE
      flXLARFT = tflops
      RETURN
      END
*********************************************************************
      INTEGER FUNCTION flXLARFB(m,n,nb)
      INTEGER m, n, nb
*     returns flop count for applying a m by nb block reflector
*     from the left to a m by n matrix
      INTEGER t
*        XGEMM
      t = nb*(2*m*n + 2*n)
*        XTRMM
      t = t + n*nb*nb
*        XGEMM
      t = t + n*(2*m*nb+nb)
      flXLARFB = t
      RETURN
      END
*******************************************************************
      INTEGER FUNCTION flXGEQRF(m,n,nb)
      INTEGER m, n, nb
*     returns flop count for blocked QR factorization without
*     column exchanges
      INTEGER i, t, kb
      EXTERNAL flXGEQR2, flXLARFT, flXLARFB
      INTEGER  flXGEQR2, flXLARFT, flXLARFB
      t = 0
      DO 10 i = 1,min(m,n),nb
         kb = min(min(m,n)-i+1,nb)
         t = t + flXGEQR2(m-i+1,kb)
     $         + flXLARFT(m-i+1,kb)
     $         + flXLARFB(m-i+1,n-i-kb+1,kb)
10    CONTINUE
      flXGEQRF = t
      RETURN
      END
*********************************************************************
      SUBROUTINE DZLTRI( m, n, a, lda )
*     zeroes lower triangle of m-by-n matrix A
      INTEGER m, n, lda
      DOUBLE PRECISION a( lda, n )
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
      INTEGER i, j
      DO 10 j = 1, n
         DO 20 i = j+1, m
            a(i,j) = ZERO
20       CONTINUE
10    CONTINUE
      RETURN
      END
*********************************************************************
      INTEGER FUNCTION flXLAIC1(j)
      INTEGER j
*     flops for incremental condition estimation excluding
*     construction of nullvector
      flXLAIC1 = 43 + 2*j
      RETURN
      END
*********************************************************************
      LOGICAL FUNCTION iscle(vname,var,bound)
*     INTEGER scalar 'var' less equal 'bound' ?
      CHARACTER*(*) vname
      INTEGER var, bound
      IF( ABS(var) .gt. bound) then
         WRITE(*,1000) vname,var,bound
         iscle = .false.
      ELSE
         iscle = .true.
      END IF
      RETURN
1000  FORMAT(/,1x,a,' = ',i6,' > bound = ',i6)
      END
*********************************************************************
      LOGICAL FUNCTION iarle(vname,var,length,bound)
*     INTEGER array 'var' less equal 'bound' ?
      CHARACTER*(*) vname
      INTEGER length, var(length), bound, i
      DO 10 i = 1,length
         IF( ABS(var(i)) .gt. bound) then
            WRITE(*,1000) vname,i,ABS(var(i)), bound
            iarle = .false.
            RETURN
         END IF
10    CONTINUE
      iarle = .true.
      RETURN
1000  FORMAT(/,1x,a,'(',i3,') = ',i6,' > bound = ',i6)
      END
*********************************************************************
      DOUBLE PRECISION FUNCTION Dckqrf( m, n, qt, ldqt, r, ldr,
     $                                  a, lda, jpvt, work )
      INTEGER m,n,ldqt,ldr,lda
      DOUBLE PRECISION qt(ldqt,m),r(ldr,n),a(lda,n),work(m)
      INTEGER jpvt(n)
*
*     This code computes the frobenius norm of Q'*A*P-R,
*     where permutation matrix P is defined by jpvt.
*
      DOUBLE PRECISION ZERO, ONE
      PARAMETER ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      EXTERNAL Dgemv, Dnrm2
      DOUBLE PRECISION Dnrm2
      INTRINSIC sqrt,min
      INTEGER i, j
      DOUBLE PRECISION    aux
      aux = ZERO
      DO 10 j = 1, n
*        Store j-th column of R in vector work.
         DO 20 i = 1, min(m,j)
            work(i) = r(i,j)
 20      CONTINUE
         DO 30 i = j+1, m
            work(i) = ZERO
 30      CONTINUE
*        Substract j-th column of Q'*A*P and j-th column of R.
         CALL Dgemv('No transpose',m,m,-ONE,qt,ldqt,a(1,jpvt(j)),1,
     $             ONE,work,1)
*        Accumulate the residuals.
         aux = aux + Dnrm2(m,work,1) ** 2
 10   CONTINUE
      Dckqrf = sqrt(aux)
      RETURN
      END
*******************************************************************
      DOUBLE PRECISION FUNCTION Dckort( m, n, q, ldq, work )
      INTEGER m,n,ldq
      DOUBLE PRECISION q(ldq,n), work(n)
*
*     Checks for orthogonality of matrix Q with orthogonal columns.
*     It computes the frobenius norm of Q'*Q-In, where Q is m by n,
*     Q' is the transpose of Q and In is n by n identity matrix.
*
      DOUBLE PRECISION ZERO, ONE
      PARAMETER ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      EXTERNAL Dgemv, Dnrm2
      DOUBLE PRECISION     Dnrm2
      INTRINSIC sqrt
      INTEGER i,j
      DOUBLE PRECISION    aux
      aux = ZERO
      DO 10 j = 1, n
*        Assign j-th column of In (n by n identity matrix) to vector work.
         DO 20 i = 1, n
            work(i) = ZERO
 20      CONTINUE
         work(j) = ONE
*        Compute work:= work - Q'*Q(:,j).
         CALL Dgemv('Transpose',m,n,-ONE,q,ldq,q(1,j),1,
     $              ONE,work,1)
*        Accumulate the residuals.
         aux = aux + Dnrm2(n,work,1) ** 2
 10   CONTINUE
      Dckort = sqrt(aux)
      RETURN
      END
*******************************************************************
      DOUBLE PRECISION FUNCTION Dcksvd( m, n, a, lda, svlues,
     $                                   work1,work2)
      INTEGER m, n, lda
      DOUBLE PRECISION a(lda,*), svlues(*),work1(m,*),
     $     work2(*)
*
*     compares the singular values s of the upper triangle of A
*     with the values in svlues and returns
*            || s - svlues||/||svlues||
*
      DOUBLE PRECISION ZERO, ONE
      PARAMETER ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER i, j, info, mn
      DOUBLE PRECISION dummy, nrmsvl
      EXTERNAL Dgebd2, Dbdsqr, Daxpy, Dnrm2
      DOUBLE PRECISION Dnrm2
      mn = min( m, n )
      nrmsvl = Dnrm2(mn,svlues,1)
      IF( nrmsvl .eq. ZERO) nrmsvl = ONE
*
*     Copy upper triangle of A into work1
*
      DO 10 j = 1, n
         DO 20 i = 1, min( j, m )
            work1( i, j ) = a( i, j )
20       CONTINUE
         DO 30 i = j+1,m
            work1(i,j) = ZERO
30       CONTINUE
10    CONTINUE
*
*     compute SVD of work1
*
      CALL Dgebd2(m,n,work1,m,work2(1),work2(mn+1),work2(2*mn+1),
     $            work2(3*mn+1),work2(4*mn+1),info)
      CALL Dbdsqr('upper',mn,0,0,0,work2(1),work2(mn+1),dummy,mn,
     $            dummy,1,dummy,mn,work1(1,1),info)
*
*     compare svlues and work1
*
      CALL Daxpy(mn,-ONE,svlues,1,work2,1)
      Dcksvd = Dnrm2(mn,work2,1)/nrmsvl
      RETURN
      END
*********************************************************************
      DOUBLE PRECISION FUNCTION Dckpqr(m,n,k,qr,ldq,tau,a,lda,
     $                                     jpvt,work)
      INTEGER m,n,k,lda,ldq
      DOUBLE PRECISION qr(ldq,n),tau(n),a(lda,n),work(*)
      INTEGER jpvt(n)
*
*     Let qr be the (possibly partial) QR-factorization of a matrix B,
*     i.e. the upper triangle of qr(1:k,1:k) is a partial triangular
*     factor and the entries below the diagonal in the first k columns
*     are the Householder vectors. The rest of qr contains a partially
*     updated matrix.
*     The vector tau contain the particulars of the Householder matrices.
*     jpvt contains the pivot inFORMATion
*     The required workspace is: m+n.
*
*     This FUNCTION returns || Q'A*P - R||
*
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
      INTEGER i, j, info
      DOUBLE PRECISION aux
      EXTERNAL Dorm2r, Dcopy, Dnrm2
      DOUBLE PRECISION Dnrm2
      INTRINSIC sqrt
      aux = ZERO
      DO 10 j = 1, n
*        Compute j-th column of Q'*A*P and put into vector work(1:m).
         CALL Dcopy( m,a(1,jpvt(j)),1,work(1),1)
         CALL Dorm2r( 'left','transpose', m, 1, k, qr, ldq, tau,
     $               work( 1 ), m, work( m+1 ), info )
*        Substract j-th column of Q'*A*P and j-th column of R
*        (stored in qr).
         IF( j.gt.k ) then
            DO 20 i = 1, m
               work(i) = work(i) - qr(i,j)
 20         CONTINUE
         ELSE
            DO 30 i = 1, j
               work(i) = work(i) - qr(i,j)
 30         CONTINUE
         END IF
         aux = aux + Dnrm2(m,work,1) ** 2
 10   CONTINUE
      Dckpqr = sqrt(aux)
      RETURN
      END
*******************************************************************
      LOGICAL FUNCTION find(which,an,lan)
      INTEGER which, lan, an(lan)
*     returns TRUE if 'which' is a value in 'an', FALSE otherwise.
      INTEGER i
      find = .false.
      DO 10 i = 1,lan
         IF( an(i) .eq. which) then
            find = .true.
            GOTO 20
         END IF
10    CONTINUE
20    RETURN
      END
*********************************************************************
      SUBROUTINE iZERO(n,x)
      INTEGER n, x(n)
*     ZEROs an n-vector x
      INTEGER i
      DO 10 i = 1,n
         x(i) = 0
10    CONTINUE
      RETURN
      END
*********************************************************************
      SUBROUTINE icopy(n,a,b)
      INTEGER n, a(n), b(n)
*     copies vector a into vector b
      INTEGER i
      DO 10 i = 1,n
         b(i) = a(i)
10    CONTINUE
      RETURN
      END
*********************************************************************
      INTEGER FUNCTION SFRANK(S,N,RCOND)
      INTEGER N
      DOUBLE PRECISION S(N), RCOND
*
*     returns MAX { 1 <= i <= n | s(1)/s(i) < 1/RCOND }
*     The entries of S are assumed to be nonnegative and
*     monotoniCALLy decreasing.
*
      INTEGER I
      SFRANK = 1
      DO 10 I = N,2,-1
         IF( S( 1 )*RCOND.LT.S( I ) ) THEN
            SFRANK = I
            GOTO 20
         END IF
10    CONTINUE
20    RETURN
*
*     END OF SFRANK
*
      END
*********************************************************************
      SUBROUTINE Dsort(n,x,incrx,job)
      INTEGER n, incrx
      CHARACTER*1 job
      DOUBLE PRECISION x(*)
*
*     SUBROUTINE to sort a vector
*
*     On entry:
*     ========
*
*     x    vector of length n to be sorted
*     n    length of vector
*     incrx     element spacing in x
*     job  = 'i' or 'i' sorts in increasing order
*          = 'd' or 'd' sorts in decreasing order
*          otherwise the routine returns without performing
*          any computation
*
*     On exit:
*     ========
*
*     x    sorted in the prescribed order
*
*
*     EXTERNAL entries
*     ================
*
      LOGICAL lsame
      EXTERNAL lsame
*
*     internal variables
*     ==================
*
      INTEGER i, curelt, nextelt, switch,k
      DOUBLE PRECISION temp
      switch = 0
      IF( lsame(job,'i')) switch = 1
      IF( lsame(job,'d')) switch = 2
      IF( switch .eq. 0) RETURN
      GOTO (100,200) switch
*
*     sort in increasing order
*
100   DO 10 i = n-1,1,-1
         k = i
 20      IF( k .eq. n) GOTO 10
            curelt = 1+(k-1)*incrx
            nextelt = 1 + k*incrx
            IF( x(curelt) .le. x(nextelt)) then
               GOTO 10
            ELSE
               temp = x(curelt)
               x(curelt) = x(nextelt)
               x(nextelt) = temp
            END IF
            k = k+1
            GOTO 20
 10   CONTINUE
      RETURN
*
*     sort in decreasing order
*
200   DO 30 i = n-1,1,-1
         k = i
 40      IF( k .eq. n) GOTO 30
            curelt = 1+(k-1)*incrx
            nextelt = 1 + k*incrx
            IF( x(curelt) .ge. x(nextelt)) then
               GOTO 30
            ELSE
               temp = x(curelt)
               x(curelt) = x(nextelt)
               x(nextelt) = temp
            END IF
            k = k+1
            GOTO 40
 30   CONTINUE
      RETURN
*
*     next line is last line of SUBROUTINE Dsort
      END
*********************************************************************
      SUBROUTINE isort(n,ix,job)
      INTEGER n
      CHARACTER*1 job
      INTEGER ix(*)
*
*     SUBROUTINE to sort a vector of INTEGERs
*
*     On entry:
*     ========
*
*     ix   vector of length n to be sorted
*     n    length of vector
*     job  = 'i' or 'i' sorts in increasing order
*          = 'd' or 'd' sorts in decreasing order
*          otherwise the routine returns without performing
*          any computation
*
*     On exit:
*     ========
*
*     ix    sorted in the prescribed order
*
*     EXTERNALs:
*     =========
*
      LOGICAL lsame
      EXTERNAL lsame
*
*     internal variables
*     ==================
*
      INTEGER i, curelt, nextelt, switch, k, temp
      switch = 0
      IF( lsame(job,'i')) switch = 1
      IF( lsame(job,'d')) switch = 2
      IF( switch .eq. 0) RETURN
      GOTO (100,200) switch
*
*     sort in increasing order
*
100   DO 10 i = n-1,1,-1
         k = i
 20      IF( k .eq. n) GOTO 10
            curelt = 1+(k-1)
            nextelt = 1 + k
            IF( ix(curelt) .le. ix(nextelt)) then
               GOTO 10
            ELSE
               temp = ix(curelt)
               ix(curelt) = ix(nextelt)
               ix(nextelt) = temp
            END IF
            k = k+1
            GOTO 20
 10   CONTINUE
      RETURN
*
*     sort in decreasing order
*
200   DO 30 i = n-1,1,-1
         k = i
 40      IF( k .eq. n) GOTO 30
            curelt = 1+(k-1)
            nextelt = 1 + k
            IF( ix(curelt) .ge. ix(nextelt)) then
               GOTO 30
            ELSE
               temp = ix(curelt)
               ix(curelt) = ix(nextelt)
               ix(nextelt) = temp
            END IF
            k = k+1
            GOTO 40
 30   CONTINUE
      RETURN
*
*     next line is last line of SUBROUTINE isort
      END
*********************************************************************
      SUBROUTINE Dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
      INTEGER ldx,n,p,job
      INTEGER jpvt(p)
      DOUBLE PRECISION x(ldx,p),qraux(p),work(p)
c
c     sqrdc uses householder transFORMATions to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     On entry
c
c        x       DOUBLE PRECISION(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     INTEGER.
c                ldx is the leading dimension of the array x.
c
c        n       INTEGER.
c                n is the number of rows of the matrix x.
c
c        p       INTEGER.
c                p is the number of columns of the matrix x.
c
c        jpvt    INTEGER(p).
c                jpvt contains INTEGERs that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in ONE of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the END.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    DOUBLE PRECISION(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     INTEGER.
c                job is an INTEGER that initiates column pivoting.
c                if job .eq. 0, no pivoting is DOne.
c                if job .ne. 0, pivoting is DOne.
c
c     On RETURN
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains inFORMATion from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   DOUBLE PRECISION(p).
c                qraux contains further inFORMATion required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     sqrdc uses the following functions and subprograms.
c
c     blas saxpy,sDOt,sscal,sswap,snrm2
c     fortran ABS,MAX,min0,sqrt
c
c     internal variables
c
      DOUBLE PRECISION ZERO, ONE
      PARAMETER ( ZERO = 0.0D+0, ONE = 1.0D+0 )
c
      INTEGER j,jp,l,lp1,lup,maxj,pl,pu,jj
      DOUBLE PRECISION maxnrm,Dnrm2,tt
      DOUBLE PRECISION DDOt,nrmxl,t
      LOGICAL negj,swapj
c
c
      pl = 1
      pu = 0
      IF( job .eq. 0) GOTO 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         DO 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            IF( negj) jpvt(j) = -j
            IF( .not.swapj) GOTO 10
               IF( j .ne. pl) CALL Dswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       CONTINUE
   20    CONTINUE
         pu = p
         DO 50 jj = 1, p
            j = p - jj + 1
            IF( jpvt(j) .ge. 0) GOTO 40
               jpvt(j) = -jpvt(j)
               IF( j .eq. pu) GOTO 30
                  CALL Dswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          CONTINUE
               pu = pu - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
c
c     compute the norms of the free columns.
c
      IF( pu .lt. pl) GOTO 80
      DO 70 j = pl, pu
         qraux(j) = Dnrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 CONTINUE
   80 CONTINUE
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      DO 200 l = 1, lup
         IF( l .lt. pl .or. l .ge. pu) GOTO 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0
            maxj = l
            DO 100 j = l, pu
               IF( qraux(j) .le. maxnrm) GOTO 90
                  maxnrm = qraux(j)
                  maxj = j
   90          CONTINUE
  100       CONTINUE
            IF( maxj .eq. l) GOTO 110
               CALL Dswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       CONTINUE
  120    CONTINUE
         qraux(l) = ZERO
         IF( l .eq. n) GOTO 190
c
c           compute the householder transFORMATion for column l.
c
            nrmxl = Dnrm2(n-l+1,x(l,l),1)
            IF( nrmxl .eq. ZERO) GOTO 180
               IF( x(l,l) .ne. ZERO) nrmxl = sign(nrmxl,x(l,l))
               CALL Dscal(n-l+1,ONE/nrmxl,x(l,l),1)
               x(l,l) = ONE + x(l,l)
c
c              apply the transFORMATion to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               IF( p .lt. lp1) GOTO 170
               DO 160 j = lp1, p
                  t = - DDOT(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  CALL Daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  IF( j .lt. pl .or. j .gt. pu) GOTO 150
                  IF( qraux(j) .eq. ZERO) GOTO 150
                     tt = ONE - (ABS(x(l,j))/qraux(j))**2
                     tt = MAX(tt,ZERO)
                     t = tt
                     tt = ONE + 0.05*tt*(qraux(j)/work(j))**2
                     IF( tt .eq. ONE) GOTO 130
                        qraux(j) = qraux(j)*sqrt(t)
                     GOTO 140
  130                CONTINUE
                        qraux(j) = Dnrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
c
c              save the transFORMATion.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
**********************************************************************
