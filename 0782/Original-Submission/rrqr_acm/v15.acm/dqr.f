      PROGRAM Dqr
*
*     Test and timing program for the Rank-Revealing QR factorization.
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
*     $Date: 96/12/30 16:59:17 $
*
*     Constants:
*     =========
*
      DOUBLE PRECISION   ZERO, ONE, HUNDRED, MILLION
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   HUNDRED = 1.0D+2, MILLION = 1.0D+6 )
      INTEGER          MMAX, NMAX, MAXTEST, NBMAX, uin, uoutrl, uoutfl
      PARAMETER        ( MMAX=1001, NMAX=MMAX, NBMAX=64, MAXTEST=9,
     $                 uin=15, uoutrl=16, uoutfl=17 )
*
*     Make sure that MAXTEST<=9; otherwise some format statememts bomb.
      CHARACTER*6 infile
      PARAMETER (infile = 'Dqr.in')
*     lda, nmax test matrix array is lda -by- nmax
*     MAXTEST   maximal number of test values for any PARAMETER
*               is assumed to be less than 10 !!!
*     NBMAX     maximal blocksize
*     uin/uout  unit number for input/output file
*     infile    name of input file
      integer flppre, flpice, flppst
      common /CNTFLP/ flppre, flpice, flppst
*     flppre    variable used to accumulate the number of flops
*               performed in the preprocessing
*     flppst    variable used to accumulate the number of flops
*               performed in the postprocessing
*     flpice    variable used to accumulate the number of flops
*               performed in ICE, both in pre and postprocessing
      INTEGER            PERFMSR, TIME, FLOPS, MFLOPS, NORUNS, STFLPS,
     $                   RLFLPS, STMFLP, RLMFLP, TRANK, TRCOND,
     $                   ICFLPS, POFLPS
      PARAMETER          ( PERFMSR = 10, TIME = 1, STFLPS = 2,
     $                   STMFLP = 3, RLFLPS = 4, RLMFLP = 5, NORUNS=6,
     $                   FLOPS = 2, MFLOPS = 3, TRANK = 7, TRCOND = 8,
     $                   ICFLPS = 9, POFLPS = 10 )
*     PERFMSR    number of performance measures taken
*     TIME       total execution time in seconds
*     STFLPS     the no. of flops required by the standard algorithm
*                (xGEQR2: non-block QR factorization with no pivoting)
*     STMFLP     execution rate defined by STFLPS
*     RLFLPS     the number of floating point operations actually
*                performed.
*     RLMFLP     execution rate defined by RLFLPS.
*     NORUNS     total number of runs executed
*     FLOPS      the same as STMFLP for the classical algorithms
*                where STFLPS = RLFLPS
*     MFLOPS     the execution rate induced by FLOPS
*     TRANK      rank as returned by QR routine
*     TRCOND     inverse of estimated condition number
      INTEGER            RELMSR, ACCEPTED, ESTRK, IRCOND, DRCOND,
     $                   ISMAX,DSMAX, ISBEFOR, DSBEFOR,
     $                   ISAFTER, DSAFTER, ISMIN, DSMIN
      PARAMETER          ( RELMSR = 12,
     $                   ACCEPTED = 1, ESTRK = 2,
     $                   IRCOND = 3, DRCOND = 4,
     $                   ISMAX = 5, DSMAX = 6,
     $                   ISBEFOR = 7, DSBEFOR = 8,
     $                   ISAFTER = 9, DSAFTER = 10,
     $                   ISMIN = 11, DSMIN = 12 )
*     Let R1 = R(1:accepted,1:accepted),
*     R2 = R(1:r,1:r), where r = min(mn,accepted+1)
*
*     RELMSR    number of reliability data sampled
*     ACCEPTED  number of columns that was accepted
*     ESTRK     estimated rank (input to xGEQPX and xGEQPY)
*     IRCOND    inverse of condition number of R1
*     DRCOND    the factor rcond_hat/rcond
*     ISMAX     largest singular value of R1
*     DSMAX     the factor smax/smax_hat
*     ISBEFOR   smallest singular value of R1
*     DSBEFOR   the factor sbefor_hat/sbefor
*     ISAFTER   smallest singular value of R2
*     DSAFTER   the factor safter_hat/safter
*     ISMIN     the smallest singular value of R
*     DSMIN     the factor smin_hat/smin
*
*     Indices into the 'svlues' array returned by xGEQPX and xGEQPY
*
*
*     Indices into the 'svlues' array.
*
      INTEGER            IMAX, IBEFOR, IAFTER, IMIN
      PARAMETER          ( IMAX = 1, IBEFOR = 2, IAFTER = 3, IMIN = 4 )
*     ..
      INTEGER            LINRTS, BL1, BL2
      PARAMETER          ( LINRTS = 2, BL1 = 1, BL2 = 2 )
      CHARACTER*6        LINNAMES(LINRTS)
*     LINRTS    number of versions of LINPACK QR routine with column
*               pivoting
*     BL1  blas 1 version
*     BL2  blas 2 version
*     LINNAMES  LINNAMES(i) is the name of routine i, i in {BL1,BL2}
      INTEGER            TOPTMAX, GNTST, QRMTX, LATMS
      PARAMETER          ( TOPTMAX = 3, GNTST = 1, QRMTX = 2,
     $                   LATMS = 3 )
      CHARACTER*6        ctopt( TOPTMAX )
*     TOPTMAX  number of test matrix generators
*     GNTST    testing special cases
*     QRMTX    routine xQRMTX generating strips of dependent and
*              independent columns
*     LATMS    routine xLATMS -- LAPACK test matrix generator
*     ctopt    ctopt(i) is the name or routine i,
*              i in {GNTST,QRMTX,LATMS}
      INTEGER            MODEMAX, FIX1, FIX2, FULLRK, SMFRNT, SMBACK,
     $                   NLLMTX, PETER, BREAK1, GEOM, ARITH
      PARAMETER          ( MODEMAX = 7, FIX1 = 1, FIX2 = 2, FULLRK = 3,
     $                   SMFRNT = 4, SMBACK = 5, PETER = 6, NLLMTX = 7,
     $                   BREAK1 = 2, GEOM = 3, ARITH = 4 )
      CHARACTER*40       cmode( TOPTMAX, -MODEMAX:MODEMAX )
*     MODEMAX  maximal number of different distributions that can
*              be generated by any of the test matrix generators.
*              MAKE SURE that MODEMAX >= 2*(# of options for
*                        distributions DLATMS)
*              otherwise arrays for holding 'mode' are too short.
*
*     options for xGNTST:
*     FIX1     matrix w/ rank MN/2-1
*     FIX2     cols 2:MN full rank, rest dep.
*     FULLRK   full rank
*     SMFRNT   STRIP small cols up front, rest indep.
*     SMBACK   STRIP indep. cols up front, rest dep.
*     PETER    Peter Tang's distribution: a few small singular
*              values that are very close together.
*     NLLMTX   null matrix
*
*     options for xQRMTX and XLATMS:
*     BREAK1   break1 distribution
*     GEOM     geometric distribution
*     ARITH    arithmetic distribution
*
*     cmode    cmode(topt,i) is a description of the matrix that
*              is generated by routine referred to as topt with
*              argument 'mode' set to 'i'.
*
*     INTRINSICS:
*     ===========
*
      INTRINSIC          DBLE, MAX, MIN
*
*     ********************************
*     * Variables read by input file *
*     ********************************
      CHARACTER*40       outfile
      INTEGER            m, nm, am(MAXTEST), im,
     $                   n, nn, an(MAXTEST), in,
     $                   nb, nnb, anb(MAXTEST),inb,
     $                   topt,ntopt,atopt(TOPTMAX),itopt,
     $                   nmod1,amod1(MODEMAX),
     $                   nmod2,amod2(MODEMAX),
     $                   nmod3,amod3(MODEMAX),
     $                   mode,nmode,amode(MODEMAX),imode,
     $                   strip, iseed(4), job, k
      DOUBLE PRECISION   irthresh, orthresh, gap, timmin, scale
*     outfile name of output file
*     m       number of rows of test matrix (m.LE.MMAX)
*     n       number of columns of test matrix (n.LE.NMAX)
*     n       bblock size  (nb. le. NBMAX)
*     topt    different test matrix generators to be CALLed
*                =1: CALL DGNTST
*                =2: CALL DQRMTX
*                =3: CALL DLATMS
*     mod1    distributions to generate for xGNTST
*     mod2    distributions to generate for xQRMTX
*     mod3    distributions to generate for xLATMS
*     mode    distributions to generate for matrix generator
*             actually chosen.
*     dfct     different ways of choosing estimated rank of A.
* For each of <var> in {m,n,nb,topt,mod1,mod2,mod3,mode,dfct},
* n<var> is the number of values for test for <var>, a<var> is
* an array holding these test values, and i<var> is the loop
* variable for the loop stepping through a<var>.
*     strip   choice of strip width for xGNTST and xQRMTX
*     scale   epsilon**scale is taken to multiply dependent
*             columns in xQRMTX.
*     irthresh (input) inverse of threshold for condition number of a
*             matrix
*     orthresh (output) inverse of threshold for condition number of a
*             matrix
*     gap     gap around threshold for generating dependent or
*             full rank matrices
*     timmin  minimum time for a benchmark run
*     iseed   array to initialize the random number generator
*             iseed(4) must be odd
*
*     ********************
*     * Other variables: *
*     ********************
*
*     SCALARS
*     =======
*
      INTEGER            mn, lda, i, j, nobefore, rank,
     $                   runs, ls, info, oiseed(4)
      DOUBLE PRECISION   eps, bstrcond, t1, t2, trt,
     $                   smax, smin, smaxpr, sminpr, mnrm, mnrmpr,
     $                   realsmin, temp
      CHARACTER*1        c1
      CHARACTER*80       fmt
*
*     mn       shorthand for min(m,n)
*     lda      leading dimension of A
*     nobefore the rank of A with respect to the threshold rcond.
*     rank     on output the rank determined by xGEQPB, xGEQPX and xGEQPY.
*     runs     number of runs needed to accumulate timmin seconds
*     iseed    seed for xLATMS
*     bstrcond equal to s(nobefore)/s(1)
*     ls       length of work array for xGEQPB, xGEQPX and xGEQPY
*     info     return PARAMETER of LAPACK routines
*     eps      machine precision
*     trt      total run time
*     smax     estimate for largest singular value
*     smin     estimate for smallest singular value
*
*     FOR PERFORMANCE MEASUREMENT
*     ===========================
      DOUBLE PRECISION   ttime, tircond, tdrcond, tismax,
     $                   tdsmax, tisbefor, tdsbefor, tisafter,
     $                   tdsafter, tismin, tdsmin
      INTEGER            taccptd
      DOUBLE PRECISION   lint(LINRTS,PERFMSR),
     $                   cwytnop(MAXTEST,PERFMSR),
     $                   cwytqpb(MAXTEST,PERFMSR),
     $                   cwytqpx(MAXTEST,PERFMSR),
     $                   cwytqpy(MAXTEST,PERFMSR)
*     lint      performance results for linpack routines
*     linr      reliability data for linpack routines
*               need only one value since routines DO identical operations
*     cwytnop   performance results for xGEQRF
*     cwytqpb   performance results for xGEQPB (preprocessing)
*     cwyrqpb   reliability data for xGEQPB (preprocessing)
*     cwytqpx   performance results for xGEQPX (pre + post(Chandra&Ipsen)
*     cwyrqpx   reliability data for xGEQPX (pre + post(Chandra&Ipsen)
*     cwytqpy   performance results for xGEQPY (pre + post(Pan&Tang)
*     cwyrqpy   reliability data for xGEQPY (pre + post(Pan&Tang)
*     all others
*               temporary variables for time, error and so on.
*     ARRAYS FOR MATRICES
*     ===================
      DOUBLE PRECISION   a(MMAX,NMAX), copya(MMAX,NMAX), s(MMAX),
     $                   copys(MMAX), qraux(NMAX), svlues(4)
      INTEGER            jpvt(NMAX)
*
*     a, copya           matrix to be factored
*     s, copys           singular values of A
*     jpvt               pivot vector
*     WORK SPACE
*     ==========
      DOUBLE PRECISION   work(MMAX*NMAX+4*MMAX+NMAX),wk1(mmax,mmax)
*     we check later on that the length of work is sufficient
*     make sure to keep this check consistent with changes in
*     this declaration
*
*     COMMON BLOCK FOR LAPACK ENVIRONMENT PARAMETERS
*     ==============================================
      INTEGER            nblk, nmnblk, nxover
      common /cenvir/ nblk, nmnblk, nxover
*     nblk is the ideal blocksize
*     mnblk is the minimal blocksize
*     nxover is the crossover point below which an unblocked alg is used
*
*     EXTERNAL ENTRIES
*     ================
*
      EXTERNAL           DSECND, DGNTST, DQRMTX,
     $                   DLATMS, DNRM2, iscle, iarle, find,
     $                   DLAMCH, sfrank, flXGEQPF, flXGEQRF,
     $                   flXGEQR2, DLASMX
      DOUBLE PRECISION   DSECND, DNRM2, DLAMCH, DLASMX
      INTEGER            flXGEQPF, flXGEQRF, flXGEQR2, sfrank
      LOGICAL            iscle, iarle, find
*     iscle    checks INTEGER scalar against bound
*     iarle    checks INTEGER array against bound
*     flops... number of flops of LAPACK routines
*
*     *****************************
*     * start of executable stmts *
*     *****************************
*     Initialize arrays describing testing options
*     ============================================
*
      DATA               lint(BL1,TRANK) /0/, lint(BL1,TRCOND) /0/,
     $                   (cwytnop(i,TRANK),i=1,MAXTEST) /MAXTEST*0/,
     $                   (cwytnop(i,TRCOND),i=1,MAXTEST) /MAXTEST*0/
*     set to zero in case we don't DO error check
      data               LINNAMES /'DQRDC ','DGEQPF'/
      data               ctopt /'DGNTST','DQRMTX','DLATMS'/
      cmode(GNTST,NLLMTX) = 'null matrix'
      cmode(GNTST,FIX1) = 'matrix w/ rank MN/2-1'
      cmode(GNTST,FIX2) = 'cols 2:MN full rank, rest dep.'
      cmode(GNTST,FULLRK) = 'full rank'
      cmode(GNTST,SMFRNT) = 'STRIP small cols up front, rest indep.'
      cmode(GNTST,SMBACK) = 'STRIP indep. cols up front, rest dep.'
      cmode(GNTST,PETER) = 'Peter''s: 5 small close sing. values'
      cmode(QRMTX,BREAK1) = 'break1 distribution'
      cmode(QRMTX,GEOM) = 'geometric distribution'
      cmode(QRMTX,ARITH) = 'arithmetic distribution'
      cmode(QRMTX,-BREAK1) = 'break1 distribution reversed'
      cmode(QRMTX,-GEOM) = 'geometric distribution reversed'
      cmode(QRMTX,-ARITH) = 'arithmetic distribution reversed'
      j = MAX(ARITH,BREAK1,GEOM)
      DO 230 i = -j,j
        cmode(LATMS,i) = cmode(QRMTX,i)
230   CONTINUE
*
*     Envir common block initialization
*     =================================
*
      nmnblk = 1
      nxover = 1
      lda = mmax
      eps = DLAMCH('epsilon')
*     *****************************************************
*     * read data from input file and copy to output file *
*     *****************************************************
*
      OPEN(uin,file=infile)
      REWIND(uin)
*     name of output file
      READ(uin,*) outfile
      OPEN(uoutrl,file='rank.'//outfile)
      OPEN(uoutfl,file='time.'//outfile)
      REWIND(uoutrl)
      REWIND(uoutfl)
      WRITE(uoutrl,1040) outfile
      WRITE(uoutfl,1040) outfile
*     values for m
      READ(uin,*) nm
      IF( .not. iscle('nm',nm,MAXTEST) ) STOP
      READ(uin,*)(am(i),i=1,nm)
      IF( .not. iarle('am',am,nm,lda) ) STOP
      WRITE(c1,'(i1)') nm
      fmt = '(1x,i1,'' nm'',/,1x,'//c1//'(i4,2x),'' m'')'
      WRITE(uoutrl,fmt) nm,(am(i),i=1,nm)
      WRITE(uoutfl,fmt) nm,(am(i),i=1,nm)
*     values for n
      READ(uin,*) nn
      IF( .not. iscle('nn',nn,MAXTEST) ) STOP
      READ(uin,*)(an(i),i=1,nn)
      IF( .not. iarle('an',an,nn,nmax) ) STOP
      WRITE(c1,'(i1)') nn
      fmt = '(1x,i1,'' nn'',/,1x,'//c1//'(i4,2x),'' n'')'
      WRITE(uoutrl,fmt) nn, (an(i),i=1,nn)
      WRITE(uoutfl,fmt) nn, (an(i),i=1,nn)
*     block sizes
      READ(uin,*) nnb
      IF( .not. iscle('nnb',nnb,MAXTEST) ) STOP
      READ(uin,*)(anb(i),i=1,nnb)
      IF( .not. iarle('anb',anb,nnb,NBMAX) ) STOP
      CALL isort(nnb,anb,'i')
      IF( anb(1).NE.1 ) THEN
        WRITE(*,*) '*** ERROR: specify nb = 1 as well ***'
        STOP
      END IF
      WRITE(c1,'(i1)') nnb
      fmt = '(1x,i1,'' nnb'',/,1x,'//c1//'(i2,2x),'' nb'')'
      WRITE(uoutrl,fmt) nnb, (anb(i),i=1,nnb)
      WRITE(uoutfl,fmt) nnb, (anb(i),i=1,nnb)
*     test matrix generation routines
      READ(uin,*) ntopt
      IF( .not. iscle('ntopt',ntopt,TOPTMAX) ) STOP
      READ(uin,*) (atopt(i),i=1,ntopt)
      IF( .not. iarle('atopt',atopt,ntopt,TOPTMAX) ) STOP
      WRITE(c1,'(i1)') ntopt
      fmt = '(1x,i1,'' ntopt'',/,1x,'//c1//'(i3,2x),'' topt'')'
      WRITE(uoutrl,fmt) ntopt, (atopt(i),i=1,ntopt)
      WRITE(uoutfl,fmt) ntopt, (atopt(i),i=1,ntopt)
*     test cases for xGNTST
      READ(uin,*) nmod1
      IF( .not. iscle('nmod1',nmod1,MODEMAX) ) STOP
      READ(uin,*) (amod1(i),i=1,nmod1)
      IF( .not. iarle('amod1',amod1,nmod1,MODEMAX) ) STOP
      IF( find(GNTST,atopt,ntopt) ) THEN
        WRITE(c1,'(i1)') nmod1
        fmt = '(1x,i1,'' nmod1'',/,1x,'//c1//'(i3,2x),'' mod1'')'
        WRITE(uoutrl,fmt) nmod1, (amod1(i),i=1,nmod1)
        WRITE(uoutfl,fmt) nmod1, (amod1(i),i=1,nmod1)
      END IF
      j = MAX(BREAK1,GEOM,ARITH)
*     singular value distributions for xQRMTX
      READ(uin,*) nmod2
      IF( .not. iscle('nmod2',nmod2,MODEMAX) ) STOP
      READ(uin,*) (amod2(i),i=1,nmod2)
      IF( .not. iarle('amod2',amod2,nmod2,j) ) STOP
      IF( find(QRMTX,atopt,ntopt) ) THEN
        WRITE(c1,'(i1)') nmod2
        fmt = '(1x,i1,'' nmod2'',/,1x,'//c1//'(i3,2x),'' mod2'')'
        WRITE(uoutrl,fmt) nmod2, (amod2(i),i=1,nmod2)
        WRITE(uoutfl,fmt) nmod2, (amod2(i),i=1,nmod2)
      END IF
*     singular value distributions for xLATMS
      READ(uin,*) nmod3
      IF( .not. iscle('nmod3',nmod3,MODEMAX) ) STOP
      READ(uin,*) (amod3(i),i=1,nmod3)
      IF( .not. iarle('amod3',amod3,nmod3,j) ) STOP
      IF( find(LATMS,atopt,ntopt) ) THEN
        WRITE(c1,'(i1)') nmod3
        fmt = '(1x,i1,'' nmod3'',/,1x,'//c1//'(i3,2x),'' mod3'')'
        WRITE(uoutrl,fmt) nmod3, (amod3(i),i=1,nmod3)
        WRITE(uoutfl,fmt) nmod3, (amod3(i),i=1,nmod3)
      END IF
*     strip width for xGNTST and xQRMTX
      READ(uin,*) strip
      IF( find(QRMTX,atopt,ntopt).OR.find(GNTST,atopt,ntopt) ) THEN
        WRITE(uoutrl,1050) strip
        WRITE(uoutfl,1050) strip
      END IF
*     scale factor for dependent columns in xQRMTX
      READ(uin,*) scale
      scale = eps**scale
      IF( find(QRMTX,atopt,ntopt) ) THEN
        WRITE(uoutrl,1060) scale
        WRITE(uoutfl,1060) scale
      END IF
*     inverse of acceptance threshold for condition number
      READ(uin,*) irthresh
      WRITE(uoutrl,1070) irthresh
      WRITE(uoutfl,1070) irthresh
*     gap around acceptance threshold
      READ(uin,*) gap
      WRITE(uoutrl,1080) gap
      WRITE(uoutfl,1080) gap
*     minimum time for a benchmark run
      READ(uin,*) timmin
      WRITE(uoutrl,1020) timmin
      WRITE(uoutfl,1020) timmin
*     seed for random number generator
      READ(uin,*) (iseed(i),i=1,4)
      IF( mod(iseed(4),2).EQ.0 ) THEN
        WRITE(*,1090) iseed(4)
        STOP
      END IF
      WRITE(uoutrl,1030) iseed
      WRITE(uoutfl,1030) iseed
      trt = DSECND()
*
*     save values that are overwritten by xGEQPX and xGEQPY
*
*     ***************
*     ***************
*     ** Test loop **
*     ***************
*     ***************
      DO 9001 im = 1,nm
        m = am(im)
        DO 9002 in = 1,nn
          n = an(in)
          mn = min(m,n)
          DO 9003 itopt = 1,ntopt
            topt = atopt(itopt)
            IF( topt.EQ.GNTST ) THEN
              nmode = nmod1
              CALL icopy(nmod1,amod1,amode)
            ELSEIF( topt.EQ.QRMTX ) THEN
              nmode = nmod2
              CALL icopy(nmod2,amod2,amode)
            ELSEIF( topt.EQ.LATMS ) THEN
              nmode = nmod3
              CALL icopy(nmod3,amod3,amode)
            END IF
            DO 9004 imode = 1,nmode
              mode = amode(imode)
              oiseed(1) = iseed(1)
              oiseed(2) = iseed(2)
              oiseed(3) = iseed(3)
              oiseed(4) = iseed(4)
*
*             generate test matrix of size m by n using
*             test matrix generator indicated by 'topt'
*             and singular value distribution by 'mode'.
*             *****************************************
*
              IF( topt.EQ.GNTST ) THEN
                CALL DGNTST(mode,m,n,irthresh,gap,strip,iseed,
     $                      rank,s,a,lda,work)
              ELSEIF( topt.EQ.QRMTX ) THEN
                CALL DQRMTX('all',scale,m,n,irthresh*gap,strip,
     $                      mode,iseed,rank,s,a,lda,work)
                rank = sfrank(s,mn,irthresh)
              ELSEIF( topt.EQ.LATMS ) THEN
                CALL DLATMS(M,N,'Uniform',iseed,'nonsymmetric',
     $                      s,mode,gap/irthresh,ONE,m,n,
     $                      'no packing',a,lda,work,info)
                IF( mode.LT.0 ) THEN
                   CALL Dsort(mn,s,1,'decreasing')
                END IF
                rank = sfrank(s,mn,irthresh)
              END IF
*
*             Save A, its singular values, and the acceptance
*             threshold as well as the best condition number for
*             R that can be achieved. Also save rank with
*             respect to 'irthresh'.
*             ==================================================
*
              CALL DLACPY('all',m,n,a,lda,copya,lda)
              CALL DCOPY(mn,s,1,copys,1)
              IF( rank.GT.0 ) THEN
                nobefore = rank
                bstrcond = s(nobefore)/s(1)
              ELSE
                nobefore = 1
                bstrcond = ZERO
              END IF
              realsmin = s(mn)
*
*             Write info to output file and console
*             =====================================
              WRITE(uoutrl,1200) m,n,ctopt(topt),cmode(topt,mode),
     $                           strip,irthresh,gap,
     $                           rank, nobefore,bstrcond,
     $                           s(1),s(nobefore),
     $                           s(min(mn,nobefore+1)), s(mn),
     $                           oiseed
              WRITE(uoutfl,1200) m,n,ctopt(topt),cmode(topt,mode),
     $                           strip,irthresh,gap,
     $                           rank, nobefore,bstrcond,
     $                           s(1),s(nobefore),
     $                           s(min(mn,nobefore+1)), s(mn),
     $                           oiseed
*             **********************
*             *  LINPACK QR BLAS 1 *
*             **********************
*             DO enough runs for the aggregate runtime to be more
*             than ''timmin''
*
              runs = 0
              t1 = DSECND()
10            CONTINUE
                CALL izero(n,jpvt)
                CALL DLACPY('all',m,n,copya,lda,a,lda)
                CALL DQRDC(a,lda,m,n,qraux,jpvt,work,1)
                t2 = DSECND()
                ttime = t2 - t1
                runs = runs + 1
              IF( ttime.LT.timmin ) GOTO 10
*
*             subtract the time for the DLACPY calls
*
              CALL DLACPY('all',m,n,a,lda,wk1,m)
              t1 = DSECND()
              DO 20 j = 1,runs
                CALL DLACPY('all',m,n,copya,lda,a,lda)
20            CONTINUE
              ttime = (ttime - (DSECND() - t1))/runs
              IF( ttime.LE.ZERO ) ttime = ZERO
              CALL DLACPY('all',m,n,wk1,m,a,lda)
*
*             save test results
*
              lint(BL1,TIME) = ttime
              lint(BL1,FLOPS) = DBLE(flXGEQR2(m,n))
              IF( ttime .EQ. ZERO ) THEN
                lint(BL1,MFLOPS) = ZERO
              ELSE
                lint(BL1,MFLOPS) = lint(BL1,FLOPS)/ttime/MILLION
              END IF
              lint(BL1,NORUNS) = DBLE(runs)
*
*
*             ************************************
*             * LAPACK QR WITH PIVOTING (BLAS-2) *
*             ************************************
*             DO enough runs for the aggregate runtime to be more
*             than ''timmin''
*
              runs = 0
              t1 = DSECND()
100           CONTINUE
                CALL izero(n,jpvt)
                CALL DLACPY('all',m,n,copya,lda,a,lda)
                CALL DGEQPF(m,n,a,lda,jpvt,qraux,work,info)
                t2 = DSECND()
                ttime = t2 - t1
                runs = runs + 1
              IF( ttime.LT.timmin ) GOTO 100
*
*             subtract the time for the DLACPY calls
*
              CALL DLACPY('all',m,n,a,lda,wk1,m)
              t1 = DSECND()
              DO 110 j = 1,runs
                CALL DLACPY('all',m,n,copya,lda,a,lda)
110           CONTINUE
              ttime = (ttime - (DSECND() - t1))/runs
              IF( ttime.LE.ZERO ) ttime = ZERO
              CALL DLACPY('all',m,n,wk1,m,a,lda)
*
*             determine what columns would be accepted
*
              smax = abs(a(1,1))
              smin = abs(a(1,1))
              mnrm = smin
              IF( abs(a(1,1)).LT.irthresh ) THEN
                taccptd = 0
                tircond = abs(a(1,1))
                tdrcond = ONE
                tismax = tircond
                tdsmax = ONE
                tisbefor = tircond
                tdsbefor = ONE
                tisafter = abs(a(min(mn,2),min(mn,2)))
                tdsafter = ONE
                tismin = abs(a(mn,mn))
                IF( realsmin.GT.ZERO ) THEN
                  tdsmin = tismin/realsmin
                ELSE
                  tdsmin = ONE
                END IF
              ELSE
                DO 120 i = 1,mn
                  mnrmpr = MAX(mnrm,DNRM2(i,a(1,i),1))
                  smaxpr = DLASMX(i)*mnrmpr
                  sminpr = min(smin,abs(a(i,i)))
                  IF( smaxpr*irthresh.GT.sminpr ) THEN
                    taccptd = i - 1
                    GOTO 130
                  ELSE
                    smax = smaxpr
                    smin = sminpr
                    mnrm = mnrmpr
                  END IF
120             CONTINUE
                taccptd = mn
130             CONTINUE
              END IF
*
*             save timing results
*
              lint(BL2,TIME) = ttime
              lint(BL2,FLOPS) = lint(BL1,FLOPS)
              IF( ttime.EQ.ZERO ) THEN
                lint(BL2,MFLOPS) = ZERO
              ELSE
                lint(BL2,MFLOPS) = lint(BL1,FLOPS)/ttime/MILLION
              END IF
              lint(BL2,NORUNS) = DBLE(runs)
              lint(BL2,TRANK) = DBLE(taccptd)
              IF( smax.EQ.ZERO ) THEN
                 lint(BL2,TRCOND) = ZERO
              ELSE
                 lint(BL2,TRCOND) = smin/smax
              END IF
              lint(BL1,TRANK) = lint(BL2,TRANK)
              lint(BL1,TRCOND) = lint(BL2,TRCOND)
*
*             Try for all different block sizes
*             ***********************************
              DO 9005 inb = 1,nnb
                nb = anb(inb)
                nblk = nb
*
*               ******************************
*               * LAPACK QR WITHOUT PIVOTING *
*               ******************************
*               DO enough runs for the aggregate runtime to be more
*               than ''timmin''
*
                runs = 0
                t1 = DSECND()
200             CONTINUE
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
                  CALL DGEQRF(m,n,a,lda,qraux,work,n*nb,info)
                  t2 = DSECND()
                  ttime = t2 - t1
                  runs = runs + 1
                IF( ttime.LT.timmin ) GOTO 200
*
*               subtract the time for the DLACPY calls
*
                CALL DLACPY('all',m,n,a,lda,wk1,m)
                t1 = DSECND()
                DO 220 j = 1,runs
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
220             CONTINUE
                ttime = (ttime - (DSECND() - t1))/runs
                IF( ttime.LE.ZERO ) ttime = ZERO
                CALL DLACPY('all',m,n,wk1,m,a,lda)
*
*               save test results
*
                cwytnop(inb,TIME) = ttime
                cwytnop(inb,STFLPS) = DBLE(flXGEQR2(m,n))
                cwytnop(inb,RLFLPS) = DBLE(flXGEQRF(m,n,nb))
                IF( ttime.EQ.ZERO ) THEN
                  cwytnop(inb,STMFLP) = ZERO
                  cwytnop(inb,RLMFLP) = ZERO
                ELSE
                  cwytnop(inb,STMFLP) =
     $                    cwytnop(inb,STFLPS)/ttime/MILLION
                  cwytnop(inb,RLMFLP) =
     $                    cwytnop(inb,RLFLPS)/ttime/MILLION
                END IF
                cwytnop(inb,NORUNS) = DBLE(runs)
*
*               Linpack pivoting strategy is achieved through
*               setting nb = 1
*
*               ********************************
*               * BLOCK QR WITH LOCAL PIVOTING *
*               ********************************
                job = 1
                k = 0
*               Length of work array for DGEQPB
                IF( job.EQ.1 ) THEN
                   IF( nb.EQ.1 ) THEN
                      ls = 2*min(m,n)+3*n
                   ELSE
                      ls = 2*min(m,n)+n*nb
                   END IF
                ELSE
                   IF( nb.EQ.1 ) THEN
                      ls = 2*min(m,n)+2*n+max(k,n)
                   ELSE
                      ls = 2*min(m,n)+nb*nb+nb*max(k,n)
                   END IF
                END IF
                ls = MAX(1,ls)
                IF( ls.GT.MMAX*NMAX+4*MMAX+NMAX ) THEN
                   WRITE(*,*) '**error in xqr.F(DGEQPB):',
     $                        'workspace too short'
                END IF
*               DO enough runs for the aggregate runtime to be more
*               than ''timmin''
                runs = 0
                t1 = DSECND()
300             CONTINUE
                  flppre = 0
                  flpice = 0
                  flppst = 0
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
                  CALL DGEQPB(job,m,n,k,a,lda,wk1,mmax,jpvt,
     $                        irthresh,orthresh,rank,svlues,work,ls,
     $                        info)
                  t2 = DSECND()
                  ttime = t2 - t1
                  runs = runs + 1
                IF( ttime.LT.timmin ) GOTO 300
*
*               subtract the time for the DLACPY calls
*
                CALL DLACPY('all',m,n,a,lda,wk1,m)
                t1 = DSECND()
                DO 320 j = 1,runs
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
320             CONTINUE
                ttime = (ttime - (DSECND() - t1))/runs
                IF( ttime.LE.ZERO ) ttime = ZERO
                CALL DLACPY('all',m,n,wk1,m,a,lda)
*
*               save timing results
*
                cwytqpb(inb,TIME) = ttime
                cwytqpb(inb,STFLPS) = lint(BL1,FLOPS)
                IF( ttime.EQ.ZERO ) THEN
                  cwytqpb(inb,STMFLP) = ZERO
                ELSE
                  cwytqpb(inb,STMFLP) =
     $                    cwytqpb(inb,STFLPS)/ttime/MILLION
                END IF
                cwytqpb(inb,ICFLPS) = DBLE(flpice)
                cwytqpb(inb,POFLPS) = DBLE(flppst)
                cwytqpb(inb,RLFLPS) = DBLE(flppre)+
     $                                DBLE(flpice)+DBLE(flppst)
                IF( ttime.EQ.ZERO ) THEN
                  cwytqpb(inb,RLMFLP) = ZERO
                ELSE
                  cwytqpb(inb,RLMFLP) =
     $                    cwytqpb(inb,RLFLPS)/ttime/MILLION
                END IF
                cwytqpb(inb,NORUNS) = DBLE(runs)
                cwytqpb(inb,TRANK) = DBLE(rank)
                cwytqpb(inb,TRCOND) = orthresh
*
*               ***********************************
*               * xGEQPX: PRE and POSTPROCESSING  *
*               * Modified Chandrasekaran & Ipsen *
*               ***********************************
                job = 1
                k = 0
*               Length of work array for DGEQPB
                IF( job.EQ.1 ) THEN
                   IF( nb.EQ.1 ) THEN
                      ls = 2*min(m,n)+3*n
                   ELSE
                      ls = 2*min(m,n)+n*nb
                   END IF
                ELSE
                   IF( nb.EQ.1 ) THEN
                      ls = 2*min(m,n)+2*n+max(k,n)
                   ELSE
                      ls = 2*min(m,n)+nb*nb+nb*max(k,n)
                   END IF
                END IF
                ls = MAX(1,ls)
                IF( ls.GT.MMAX*NMAX+4*MMAX+NMAX ) THEN
                   WRITE(*,*) '**error in xqr.F (DGEQPX):',
     $                        'workspace too short'
                END IF
*               DO enough runs for the aggregate runtime to be more
*               than ''timmin''
                runs = 0
                t1 = DSECND()
400             CONTINUE
                  flppre = 0
                  flpice = 0
                  flppst = 0
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
                  CALL DGEQPX(job,m,n,k,a,lda,wk1,mmax,jpvt,
     $                        irthresh,orthresh,rank,svlues,work,ls,
     $                        info)
                  t2 = DSECND()
                  ttime = t2 - t1
                  runs = runs + 1
*
*                 Check if xGEQPX was ok.
*
                  IF( info.NE.0 )
     $              WRITE(*,*) 'DGEQPX. Info:',info
                IF( ttime.LT.timmin ) GOTO 400
*
*               subtract the time for the DLACPY calls
*
                CALL DLACPY('all',m,n,a,lda,wk1,m)
                t1 = DSECND()
                DO 420 j = 1,runs
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
420             CONTINUE
                ttime = (ttime - (DSECND() - t1))/runs
                IF( ttime.LE.ZERO ) ttime = ZERO
                CALL DLACPY('all',m,n,wk1,m,a,lda)
*
*               save timing results
*
                cwytqpx(inb,TIME) = ttime
                cwytqpx(inb,STFLPS) = lint(BL1,FLOPS)
                IF( ttime.EQ.ZERO ) THEN
                  cwytqpx(inb,STMFLP) = ZERO
                ELSE
                  cwytqpx(inb,STMFLP) =
     $                    cwytqpx(inb,STFLPS)/ttime/MILLION
                END IF
                cwytqpx(inb,ICFLPS) = DBLE(flpice)
                cwytqpx(inb,POFLPS) = DBLE(flppst)
                cwytqpx(inb,RLFLPS) = DBLE(flppre)+
     $                                 DBLE(flpice)+DBLE(flppst)
                IF( ttime.EQ.ZERO ) THEN
                  cwytqpx(inb,RLMFLP) = ZERO
                ELSE
                  cwytqpx(inb,RLMFLP) =
     $                    cwytqpx(inb,RLFLPS)/ttime/MILLION
                END IF
                cwytqpx(inb,NORUNS) = DBLE(runs)
                cwytqpx(inb,TRANK) = DBLE(rank)
                cwytqpx(inb,TRCOND) = orthresh
*
*               **********************************
*               * xGEQPY: PRE and POSTPROCESSING *
*               * Modified Pan & Tang            *
*               **********************************
                job = 1
                k = 0
*               Length of work array for DGEQPB
                IF( job.EQ.1 ) THEN
                   IF( nb.EQ.1 ) THEN
                      ls = 2*min(m,n)+3*n
                   ELSE
                      ls = 2*min(m,n)+n*nb
                   END IF
                ELSE
                   IF( nb.EQ.1 ) THEN
                      ls = 2*min(m,n)+2*n+max(k,n)
                   ELSE
                      ls = 2*min(m,n)+nb*nb+nb*max(k,n)
                   END IF
                END IF
                ls = MAX(1,ls)
                IF( ls.GT.MMAX*NMAX+4*MMAX+NMAX ) THEN
                   WRITE(*,*) '**error in xqr.F (DGEQPY):',
     $                        'workspace too short'
                END IF
*               DO enough runs for the aggregate runtime to be more
*               than ''timmin''
                runs = 0
                t1 = DSECND()
500             CONTINUE
                  flppre = 0
                  flpice = 0
                  flppst = 0
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
                  CALL DGEQPY(job,m,n,k,a,lda,wk1,mmax,jpvt,
     $                        irthresh,orthresh,rank,svlues,work,ls,
     $                        info)
                  t2 = DSECND()
                  ttime = t2 - t1
                  runs = runs + 1
*
*                 Check if xGEQPY was ok.
*
                  IF( info.NE.0 )
     $              WRITE(*,*) 'DGEQPY. Info:',info
                IF( ttime.LT.timmin ) GOTO 500
*
*               subtract the time for the DLACPY calls
*
                CALL DLACPY('all',m,n,a,lda,wk1,m)
                t1 = DSECND()
                DO 520 j = 1,runs
                  CALL DLACPY('all',m,n,copya,lda,a,lda)
520             CONTINUE
                ttime = (ttime - (DSECND() - t1))/runs
                IF( ttime.LE.ZERO ) ttime = ZERO
                CALL DLACPY('all',m,n,wk1,m,a,lda)
*
*               save timing results
*
                cwytqpy(inb,TIME) = ttime
                cwytqpy(inb,STFLPS) = lint(BL1,FLOPS)
                IF( ttime.EQ.ZERO ) THEN
                  cwytqpy(inb,STMFLP) = ZERO
                ELSE
                  cwytqpy(inb,STMFLP) =
     $                    cwytqpy(inb,STFLPS)/ttime/MILLION
                END IF
                cwytqpy(inb,ICFLPS) = DBLE(flpice)
                cwytqpy(inb,POFLPS) = DBLE(flppst)
                cwytqpy(inb,RLFLPS) = DBLE(flppre)+
     $                                 DBLE(flpice)+DBLE(flppst)
                IF( ttime.EQ.ZERO ) THEN
                  cwytqpy(inb,RLMFLP) = ZERO
                ELSE
                  cwytqpy(inb,RLMFLP) =
     $                    cwytqpy(inb,RLFLPS)/ttime/MILLION
                END IF
                cwytqpy(inb,NORUNS) = DBLE(runs)
                cwytqpy(inb,TRANK) = DBLE(rank)
                cwytqpy(inb,TRCOND) = orthresh
9005          CONTINUE
*             end of inb loop
*
*             print out reliability data
*             ==========================
*
*
*             print out performance numbers
*             =============================
*
              WRITE(uoutfl,1210)
              DO 810 i = 1,LINRTS
                WRITE(uoutfl,1230)
     $                 LINNAMES(i),
     $                 lint(i,TIME),
     $                 lint(i,MFLOPS),
     $                 lint(i,MFLOPS),
     $                 int(lint(i,NORUNS)),
     $                 int(lint(i,TRANK)),
     $                 lint(i,TRCOND)
                WRITE(uoutfl,1220)
810           CONTINUE
              DO 820 inb = 1,nnb
                WRITE(uoutfl,1240)
     $                 anb(inb),
     $                 cwytnop(inb,TIME),
     $                 cwytnop(inb,STMFLP),
     $                 cwytnop(inb,RLMFLP),
     $                 int(cwytnop(inb,NORUNS))
820           CONTINUE
              WRITE(uoutfl,1220)
              DO 830 inb = 1,nnb
                nb = anb(inb)
*
                temp = cwytqpb(inb,RLFLPS)
                IF( temp.GT.0 ) THEN
                   WRITE(uoutfl,1250)'DGEQPB',
     $                    nb,
     $                    cwytqpb(inb,TIME),
     $                    cwytqpb(inb,STMFLP),
     $                    cwytqpb(inb,RLMFLP),
     $                    int(cwytqpb(inb,NORUNS)),
     $                    int(cwytqpb(inb,TRANK)),
     $                    cwytqpb(inb,TRCOND),
     $                    HUNDRED*cwytqpb(inb,ICFLPS)/temp,
     $                    HUNDRED*cwytqpb(inb,POFLPS)/temp
                ELSE
                   WRITE(uoutfl,1260)'DGEQPB',
     $                    nb,
     $                    cwytqpb(inb,TIME),
     $                    cwytqpb(inb,STMFLP),
     $                    cwytqpb(inb,RLMFLP),
     $                    int(cwytqpb(inb,NORUNS)),
     $                    int(cwytqpb(inb,TRANK)),
     $                    cwytqpb(inb,TRCOND)
                END IF
*
                temp = cwytqpx(inb,RLFLPS)
                IF( temp.GT.0 ) THEN
                   WRITE(uoutfl,1250)'DGEQPX',
     $                    nb,
     $                    cwytqpx(inb,TIME),
     $                    cwytqpx(inb,STMFLP),
     $                    cwytqpx(inb,RLMFLP),
     $                    int(cwytqpx(inb,NORUNS)),
     $                    int(cwytqpx(inb,TRANK)),
     $                    cwytqpx(inb,TRCOND),
     $                    HUNDRED*cwytqpx(inb,ICFLPS)/temp,
     $                    HUNDRED*cwytqpx(inb,POFLPS)/temp
                ELSE
                   WRITE(uoutfl,1260)'DGEQPX',
     $                    nb,
     $                    cwytqpx(inb,TIME),
     $                    cwytqpx(inb,STMFLP),
     $                    cwytqpx(inb,RLMFLP),
     $                    int(cwytqpx(inb,NORUNS)),
     $                    int(cwytqpx(inb,TRANK)),
     $                    cwytqpx(inb,TRCOND)
                END IF
*
                temp = cwytqpy(inb,RLFLPS)
                IF( temp.GT.0 ) THEN
                   WRITE(uoutfl,1250)'DGEQPY',
     $                    nb,
     $                    cwytqpy(inb,TIME),
     $                    cwytqpy(inb,STMFLP),
     $                    cwytqpy(inb,RLMFLP),
     $                    int(cwytqpy(inb,NORUNS)),
     $                    int(cwytqpy(inb,TRANK)),
     $                    cwytqpy(inb,TRCOND),
     $                    HUNDRED*cwytqpy(inb,ICFLPS)/temp,
     $                    HUNDRED*cwytqpy(inb,POFLPS)/temp
                ELSE
                   WRITE(uoutfl,1260)'DGEQPY',
     $                    nb,
     $                    cwytqpy(inb,TIME),
     $                    cwytqpy(inb,STMFLP),
     $                    cwytqpy(inb,RLMFLP),
     $                    int(cwytqpy(inb,NORUNS)),
     $                    int(cwytqpy(inb,TRANK)),
     $                    cwytqpy(inb,TRCOND)
                END IF
*
                WRITE(uoutfl,1220)
830           CONTINUE
9004        CONTINUE
*           end of imode loop
9003      CONTINUE
*         end of itopt loop
9002    CONTINUE
*       end of in loop
9001  CONTINUE
*     end of im loop
      trt = DSECND() - trt
      WRITE(uoutfl,1000) trt
      WRITE(uoutrl,1000) trt
      WRITE(*,*) ' End of program'
      CLOSE( uin )
      CLOSE( uoutrl )
      CLOSE( uoutfl )
      STOP
*
1000  FORMAT(/,1x,'total run time: ',f8.2,' seconds')
1010  FORMAT(1x,42('*'),/,1x,'* ','time of run: ',a25,' *',
     $       /,1x,42('*'),/)
1020  FORMAT(1x,f5.3,' minimum time for benchmark run')
1030  FORMAT(1x,4(i5,2x),' seed for RN generator',//)
1040  FORMAT('''',a,'''','  output file')
1050  FORMAT(1x,i3,' strip width')
1060  FORMAT(1x,e8.2,' scale for dependent columns in SQRMTX')
1070  FORMAT(1x,e8.2,' inverse of acceptance threshold')
1080  FORMAT(1x,e8.2,' gap around acceptance threshold')
1090   FORMAT('*** error: iseed(3) = ',i4,' but should be odd')
*
1200  FORMAT(/,1x,74('*'),/,1x,'* m:',i4,'  n:',i4,
     $       '  using ',a6,2x,a40,' *',
     $       /,1x,'* ',16x,'strip: ',i2,
     $       '  rthresh: ',e8.2,'  gap: ',e8.2,11x,' *',
     $       /,1x,'* ','rank: ',i4,5x,'nobefore: ',i4,5x,
     $       ' best rcond: ',e8.2,15x,
     $       ' *',/,1x,'* smax: ',e8.2,'  sbefore: ',e8.2,
     $       '  safter: ',e8.2,'  smin: ',e8.2,'    *',
     $       /,1x,'* seed: ',4(i4,2x),41(' '),'*',
     $       /,1x,74('*'))
1210  FORMAT(/,1x,12x,' | ','time(secs)',' | ','mflops(std)',
     $       ' | ','mflops(real)',' | ','noruns',' | ',
     $       'rank',' | ',' rcond  ',' | ','% ice',' | ','%post',
     $       /,1x,100('='))
1220  FORMAT(1x,100('-'))
1230  FORMAT(1x,a6,'       |  ',e8.2,'  |   ',f7.2,'   |   ',
     $       f7.2,'    | ',1x,i4,1x,' | ',i4,' | ',e8.2)
1240  FORMAT(1x,'DGEQRF',' nb:',i2,' |  ',e8.2,
     $       '  |   ',f7.2,'   |   ',f7.2,'    | ',1x,i4)
1250  FORMAT(1x,a6,' nb:',i2,' |  ',e8.2,'  |   ',
     $       f7.2,'   |   ',f7.2,'    | ',1x,i4,1x,
     $       ' | ',i4,' | ',e8.2,' | ',f5.2,' | ',f5.2)
1260  FORMAT(1x,a6,' nb:',i2,' |  ',e8.2,'  |   ',
     $       f7.2,'   |   ',f7.2,'    | ',1x,i4,1x,
     $       ' | ',i4,' | ',e8.2)
*
1400  FORMAT(/,1x,'DQRDC: ',6x,'  svd:',e8.2)
1410  FORMAT(/,1x,'DGEQPF:',6x,'  svd:',e8.2,'  qrf:',e8.2,
     $       '  ort:',e8.2)
1420  FORMAT(/,1x,'DGEQRF: nb:',i2,'  svd:',e8.2,
     $       '  qrf:',e8.2,'  ort:',e8.2)
1430  FORMAT(1x,'DGEQPB: nb:',i2,
     $       '  svd:',e8.2,'  qrf:',e8.2,'  ort:',e8.2)
1440  FORMAT(1x,'DGEQPX: nb:',i2,
     $       '  svd:',e8.2,'  qrf:',e8.2,'  ort:',e8.2)
1450  FORMAT(1x,'DGEQPY: nb:',i2,
     $       '  svd:',e8.2,'  qrf:',e8.2,'  ort:',e8.2)
*
*
*
       END
