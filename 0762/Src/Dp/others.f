      SUBROUTINE dpsifn(x,n,kode,m,ans)
*include *dpsifndc

C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER kode,m,n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ans(m)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION arg,den,elim,eps,fln,fn,fnp,fns,fx,r1m4,r1m5,rln,
     +                 rxsq,s,slope,t,t1,t2,ta,tk,tol,tols,tss,tst,tt,
     +                 wdtol,xdmln,xdmy,xinc,xln,xm,xmin,xq,yint
      INTEGER i,iflag,iset,j,k,mx,nmax,nn,np,nx
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION b(22),trm(22),trmr(100)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION d1mach
      INTEGER i1mach
      EXTERNAL d1mach,i1mach
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,dmin1,exp,int,min0
C     ..
C     .. Data statements ..
c-----------------------------------------------------------------------
c             bernoulli numbers
c-----------------------------------------------------------------------
      DATA nmax/100/
      DATA b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),b(9),b(10),b(11),
     +     b(12),b(13),b(14),b(15),b(16),b(17),b(18),b(19),b(20),b(21),
     +     b(22)/1.00000000000000000d+00,-5.00000000000000000d-01,
     +     1.66666666666666667d-01,-3.33333333333333333d-02,
     +     2.38095238095238095d-02,-3.33333333333333333d-02,
     +     7.57575757575757576d-02,-2.53113553113553114d-01,
     +     1.16666666666666667d+00,-7.09215686274509804d+00,
     +     5.49711779448621554d+01,-5.29124242424242424d+02,
     +     6.19212318840579710d+03,-8.65802531135531136d+04,
     +     1.42551716666666667d+06,-2.72982310678160920d+07,
     +     6.01580873900642368d+08,-1.51163157670921569d+10,
     +     4.29614643061166667d+11,-1.37116552050883328d+13,
     +     4.88332318973593167d+14,-1.92965793419400681d+16/
C     ..
C     .. Executable Statements ..
c
c
      iflag = 0
      iset = 1
      IF (x.LE.0.0d0) GO TO 300
   10 CONTINUE
      iset = 2
      IF (n.LT.0) GO TO 300
   20 CONTINUE
      iset = 3
      IF (kode.LT.1 .OR. kode.GT.2) GO TO 300
   30 CONTINUE
      iset = 4
      IF (m.LT.1) GO TO 300
   40 CONTINUE
      IF (iflag.NE.0) GO TO 360
      nn = n + m - 1
      fn = dble(nn)
      fnp = fn + 1.0d0
      nx = min0(-i1mach(12),i1mach(13))
      r1m5 = d1mach(5)
      r1m4 = d1mach(4)*0.5d0
      wdtol = dmax1(r1m4,0.5d-18)
c-----------------------------------------------------------------------
c     elim = approximate exponential over and underflow limit
c-----------------------------------------------------------------------
      elim = 2.302d0* (dble(nx)*r1m5-3.0d0)
      xln = dlog(x)
      t = fnp*xln
c-----------------------------------------------------------------------
c     overflow and underflow test for small and large x
c-----------------------------------------------------------------------
      IF (abs(t).GT.elim) GO TO 290
      IF (x.LT.wdtol) GO TO 260
c-----------------------------------------------------------------------
c     compute xmin and the number of terms of the series, fln+1
c-----------------------------------------------------------------------
      rln = r1m5*dble(i1mach(11))
      rln = dmin1(rln,18.06d0)
      fln = dmax1(rln,3.0d0) - 3.0d0
      yint = 3.50d0 + 0.40d0*fln
      slope = 0.21d0 + fln* (0.0006038d0*fln+0.008677d0)
      xm = yint + slope*fn
      mx = int(xm) + 1
      xmin = dble(mx)
      IF (n.EQ.0) GO TO 50
      xm = -2.302d0*rln - dmin1(0.0d0,xln)
      fns = dble(n)
      arg = xm/fns
      arg = dmin1(0.0d0,arg)
      eps = exp(arg)
      xm = 1.0d0 - eps
      IF (abs(arg).LT.1.0d-3) xm = -arg
      fln = x*xm/eps
      xm = xmin - x
      IF (xm.GT.7.0d0 .AND. fln.LT.15.0d0) GO TO 200
   50 CONTINUE
      xdmy = x
      xdmln = xln
      xinc = 0.0d0
      IF (x.GE.xmin) GO TO 60
      nx = int(x)
      xinc = xmin - dble(nx)
      xdmy = x + xinc
      xdmln = dlog(xdmy)
   60 CONTINUE
c-----------------------------------------------------------------------
c     generate w(n+m-1,x) by the asymptotic expansion
c-----------------------------------------------------------------------
      t = fn*xdmln
      t1 = xdmln + xdmln
      t2 = t + xdmln
      tk = dmax1(abs(t),abs(t1),abs(t2))
      IF (tk.GT.elim) GO TO 380
      tss = exp(-t)
      tt = 0.5d0/xdmy
      t1 = tt
      tst = wdtol*tt
      IF (nn.NE.0) t1 = tt + 1.0d0/fn
      rxsq = 1.0d0/ (xdmy*xdmy)
      ta = 0.5d0*rxsq
      t = fnp*ta
      s = t*b(3)
      IF (abs(s).LT.tst) GO TO 80
      tk = 2.0d0
      DO 70 k = 4,22
          t = t* ((tk+fn+1.0d0)/ (tk+1.0d0))* ((tk+fn)/ (tk+2.0d0))*rxsq
          trm(k) = t*b(k)
          IF (abs(trm(k)).LT.tst) GO TO 80
          s = s + trm(k)
          tk = tk + 2.0d0
   70 CONTINUE
   80 CONTINUE
      s = (s+t1)*tss
      IF (xinc.EQ.0.0d0) GO TO 100
c-----------------------------------------------------------------------
c     backward recur from xdmy to x
c-----------------------------------------------------------------------
      nx = int(xinc)
      np = nn + 1
      IF (nx.GT.nmax) GO TO 390
      IF (nn.EQ.0) GO TO 160
      xm = xinc - 1.0d0
      fx = x + xm
c-----------------------------------------------------------------------
c     this loop should not be changed. fx is accurate when x is small
c-----------------------------------------------------------------------
      DO 90 i = 1,nx
          trmr(i) = fx** (-np)
          s = s + trmr(i)
          xm = xm - 1.0d0
          fx = x + xm
   90 CONTINUE
  100 CONTINUE
      ans(m) = s
      IF (fn.EQ.0.0d0) GO TO 180
c-----------------------------------------------------------------------
c     generate lower derivatives, j.lt.n+m-1
c-----------------------------------------------------------------------
      IF (m.EQ.1) RETURN
      DO 150 j = 2,m
          fnp = fn
          fn = fn - 1.0d0
          tss = tss*xdmy
          t1 = tt
          IF (fn.NE.0.0d0) t1 = tt + 1.0d0/fn
          t = fnp*ta
          s = t*b(3)
          IF (abs(s).LT.tst) GO TO 120
          tk = 3.0d0 + fnp
          DO 110 k = 4,22
              trm(k) = trm(k)*fnp/tk
              IF (abs(trm(k)).LT.tst) GO TO 120
              s = s + trm(k)
              tk = tk + 2.0d0
  110     CONTINUE
  120     CONTINUE
          s = (s+t1)*tss
          IF (xinc.EQ.0.0d0) GO TO 140
          IF (fn.EQ.0.0d0) GO TO 160
          xm = xinc - 1.0d0
          fx = x + xm
          DO 130 i = 1,nx
              trmr(i) = trmr(i)*fx
              s = s + trmr(i)
              xm = xm - 1.0d0
              fx = x + xm
  130     CONTINUE
  140     CONTINUE
          mx = m - j + 1
          ans(mx) = s
          IF (fn.EQ.0.0d0) GO TO 180
  150 CONTINUE
      RETURN
c-----------------------------------------------------------------------
c     recursion for n = 0
c-----------------------------------------------------------------------
  160 CONTINUE
      DO 170 i = 1,nx
          s = s + 1.0d0/ (x+dble(nx-i))
  170 CONTINUE
  180 CONTINUE
      IF (kode.EQ.2) GO TO 190
      ans(1) = s - xdmln
      RETURN

  190 CONTINUE
      IF (xdmy.EQ.x) RETURN
      xq = xdmy/x
      ans(1) = s - dlog(xq)
      RETURN
c-----------------------------------------------------------------------
c     compute by series (x+k)**(-(n+1)) , k=0,1,2,...
c-----------------------------------------------------------------------
  200 CONTINUE
      nn = int(fln) + 1
      np = n + 1
      t1 = (fns+1.0d0)*xln
      t = exp(-t1)
      s = t
      den = x
      DO 210 i = 1,nn
          den = den + 1.0d0
          trm(i) = den** (-np)
          s = s + trm(i)
  210 CONTINUE
      ans(1) = s
      IF (n.NE.0) GO TO 220
      IF (kode.EQ.2) ans(1) = s + xln
  220 CONTINUE
      IF (m.EQ.1) RETURN
c-----------------------------------------------------------------------
c     generate higher derivatives, j.gt.n
c-----------------------------------------------------------------------
      tol = wdtol/5.0d0
      DO 250 j = 2,m
          t = t/x
          s = t
          tols = t*tol
          den = x
          DO 230 i = 1,nn
              den = den + 1.0d0
              trm(i) = trm(i)/den
              s = s + trm(i)
              IF (trm(i).LT.tols) GO TO 240
  230     CONTINUE
  240     CONTINUE
          ans(j) = s
  250 CONTINUE
      RETURN
c-----------------------------------------------------------------------
c     small x.lt.unit round off
c-----------------------------------------------------------------------
  260 CONTINUE
      ans(1) = x** (-n-1)
      IF (m.EQ.1) GO TO 280
      k = 1
      DO 270 i = 2,m
          ans(k+1) = ans(k)/x
          k = k + 1
  270 CONTINUE
  280 CONTINUE
      IF (n.NE.0) RETURN
      IF (kode.EQ.2) ans(1) = ans(1) + xln
      RETURN

  290 CONTINUE
      IF (t.GT.0.0d0) GO TO 380
      GO TO 370
c-----------------------------------------------------------------------
c     error messages
c-----------------------------------------------------------------------
  300 IF (iflag.NE.0) GO TO 310
  310 GO TO (320,330,340,350),iset

  320 WRITE (*,*) ' PSIFN, X is not positive.'
      iflag = 1
      GO TO 10

  330 WRITE (*,*) ' PSIFN, N is not greater than or equal to zero.'
      iflag = 1
      GO TO 20

  340 WRITE (*,*) ' PSIFN, KODE is not 1 or 2.'
      iflag = 1
      GO TO 30

  350 WRITE (*,*) ' PSIFN, M is not greater than zero.'
      iflag = 1
      GO TO 40

  360 STOP ' PSIFN, end input errors for PSIFN.'

  370 STOP ' PSIFN, overflow, X too small or N+M-1 too large.'

  380 STOP ' PSIFN, underflow, X too large or N+M-1 too large.'

  390 STOP ' PSIFN, increase the dimension of trmr(NMAX).'

      END
      SUBROUTINE gtcuio(inun,outun)
C**********************************************************************
C
C     SUBROUTINE GTCUIO(INUN,OUTUN)
C
C                         GeT CUrrent Input Output units
C
C
C                              Function
C
C
C          Returns the FORTRAN unit numbers of the current input and
C     output units.  These units are used for the standard reads and
C     writes but may be reassigned from the standard units.  This is
C     convenient on input, for example, when the echo of input from
C     a previous interactive run is used, and on output when a report
C     should be directed to a file rather than a terminal.
C
C
C                              Arguments
C
C
C     INUN <-- The FORTRAN unit number of the current input unit.
C                                                  INUN is INTEGER
C
C     OUTUN <-- The FORTRAN unit number of the current output unit.
C                                                  OUTUN is INTEGER
C
C
C                              Note
C
C
C          If GTCUIO is called before the current units are set by
C     a call to STCUIO, then the standard input and output unit
C     numbers should be returned.
C
C
C----------------------------------------------------------------------
C
C
C     ENTRY STCUIO(INUN,OUTUN)
C
C                    SeT CUrrent Input Output units
C
C
C                              Function
C
C
C          Sets the FORTRAN unit numbers of the current input and
C     output units.  It is the responsibility of the programmer
C     to assure that whatever action is necessary for opening these
C     units for reading and writing has been taken before the
C     units are used.
C
C
C                              Arguments
C
C
C     INUN --> The FORTRAN unit number of the current input unit.
C                                                  INUN is INTEGER
C
C     OUTUN --> The FORTRAN unit number of the current output unit.
C                                                  OUTUN is INTEGER
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      INTEGER inun,outun
C     ..
C     .. Local Scalars ..
      INTEGER in,out
      LOGICAL qset
C     ..
C     .. External Subroutines ..
      EXTERNAL gtstio
C     ..
C     .. Save statement ..
      SAVE in,out,qset
C     ..
C     .. Data statements ..
      DATA qset/.FALSE./
C     ..
C     .. Executable Statements ..
C
      IF (qset) GO TO 10
      CALL gtstio(in,out)
   10 inun = in
C
      outun = out
      RETURN
C
      ENTRY stcuio(inun,outun)

      in = inun
      out = outun
      qset = .TRUE.
      RETURN

      END
      SUBROUTINE gtecun(inun,allun)
C**********************************************************************
C
C     SUBROUTINE GTECUN(INUN,ALLUN)
C
C                    GeT ECho UNits
C
C
C                              Function
C
C
C          Returns the FORTRAN unit numbers of the units to receive
C     the echo of current input and both current input and current
C     output.  A negative value for either unit number indicates
C     that the corresponding echoing function is not desired.
C
C
C                              Note
C
C
C          If GTECUN is called before the units are set with STECUN,
C     then a negative value is returned.
C          An EOF (end of file) from input is echoed as a line
C     consisting of a single blank character.
C
C
C                              Arguments
C
C
C     INUN <-- The FORTRAN unit number of the unit to receive the
C              echo of everything read from current input.  A
C              negative value indicates that this echo function is
C              not to be performed.
C                                                  INUN is INTEGER
C
C     ALLUN <-- The FORTRAN unit number of the unit to receive the
C               echo of everything read from current input or written
C               to current output.  A negative value indicates that
C               this echo function is not to be performed.
C                                                  ALLUN is INTEGER
C
C
C----------------------------------------------------------------------
C
C
C     ENTRY STECUN(INUN,ALLUN)
C
C                    Set ECho Units
C
C                         Function
C
C          Sets the units for the echo of input and echo of input and
C     output functions.  A negative unit number indicates that the
C     function is not to be performed.
C
C
C                              Arguments
C
C
C     INUN --> The FORTRAN unit number of the unit to receive the
C              echo of everything read from current input.  A
C              negative value indicates that this echo function is
C              not to be performed.
C                                                  INUN is INTEGER
C
C     ALLUN --> The FORTRAN unit number of the unit to receive the
C               echo of everything read from current input or written
C               to current output.  A negative value indicates that
C               this echo function is not to be performed.
C                                                  ALLUN is INTEGER
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      INTEGER allun,inun
C     ..
C     .. Local Scalars ..
      INTEGER all,in
      LOGICAL qset
C     ..
C     .. Save statement ..
      SAVE in,all,qset
C     ..
C     .. Data statements ..
      DATA qset/.FALSE./
C     ..
C     .. Executable Statements ..
C
      IF (.NOT. (qset)) GO TO 10
      inun = in
      allun = all
      GO TO 20

   10 inun = -1
      allun = -1
   20 RETURN
C
C
      ENTRY stecun(inun,allun)

      in = inun
      all = allun
      qset = .TRUE.
      RETURN

      END
      SUBROUTINE gtstio(inun,outun)
C**********************************************************************
C
C     SUBROUTINE GTSTIO(INUN,OUTUN)
C
C                    GeT STandard Input Output units
C                      (Machine Dependent Routine)
C
C
C                              Function
C
C
C     Returns the FORTRAN unit numbers corresponding to standard input
C     and output (i.e., the * unit).  These units are usually connected
C     to a terminal for an interactive run.
C
C
C                              Arguments
C
C
C     INUN <-- The unit number corresponding to FORTRAN standard input.
C                                                  INUN is INTEGER
C
C     OUTUN <-- The unit number corresponding to FORTRAN standard output
C                                                  OUTUN is INTEGER
C
C
C                              Note
C
C
C          On the first call, the routine should take whatever action
C     is necessary to associate the unit numbers returned with standard
C     input and output.  The unit numbers should be kept as SAVE
C     variable so that they can be returned on successive calls.
C     The unit numbers returned by this routine will not change
C     within a single run of a program.
C
C**********************************************************************
C**********************************************************************
C
C     Code last modified 87/Feb/20
C
C**********************************************************************
C *** VARIABLES
C
C
C  ** ARGUMENTS
C
C     INUN   <--  UNIT NUMBER CORRESPONDING TO THE STANDARD INPUT UNIT
      INTEGER inun
C
C     OUTUN  <--  UNIT NUMBER CORRESPONDING TO THE STANDARD OUTPUT UNIT
      INTEGER outun
C
C
C
C *** BODY OF ROUTINE
C
C
      inun = 5
      outun = 6
C
      RETURN
C
C
      END
      INTEGER FUNCTION igtfun()
C**********************************************************************
C
C     INTEGER FUNCTION IGTFUN()
C
C                    I GeT Free Unit Number
C
C
C                              Function
C
C
C          Returns a FORTRAN unit number that is not in current use.
C          Returns an invalid unit number (generally less than zero)
C          if no unassigned unit numbers are available.
C
C**********************************************************************
C *** OVERVIEW OF PROGRAM UNIT
C
C
C  ** USAGE
C
C   * OPERATION
C
C     This function is called with no arguments once for each new unit
C     number that is needed by any program unit during execution of
C     a program.
C
C
C
C *** VARIABLES
C
C
C  ** PARAMETERS AND CONSTANTS
C
C   * PROCESSOR-DEPENDENT PARAMETERS
C
C     UNUMMN  --  minimum possible unit number for a file
C
C
C     UNUMMX  --  maximum possible unit number for a file
C
C
C   * CONSTANTS
C
C     ONE     --  mnemonically-named integer constant
C
C
C
C  ** OTHER SIGNIFICANT VARIABLES
C
C     LASTUN  --  last-returned unit number
C
C
C
C  ** MISCELLANEOUS VARIABLES
C     The following variables are used for such purposes as indices and
C     temporary storage.  Each one is intended to be significant only
C     within a well-defined block of code (although it may appear in
C     more than one block), and its meaning and use should be apparent
C     from its name and context.
C
C     QEXIST  --
C
C     QOPEN   --
C
C
C
C
C *** DATA STATEMENTS
C
C     LASTUN  --  initiated to cause first invocation to start at
C                 maximum value
C
C
C
C
C
C     .. Parameters ..
      INTEGER unummn
      PARAMETER (unummn=0)
      INTEGER unummx
      PARAMETER (unummx=99)
      INTEGER one
      PARAMETER (one=1)
C     ..
C     .. Local Scalars ..
      INTEGER lastun
      LOGICAL qexist,qopen
C     ..
C     .. Save statement ..
      SAVE lastun
C     ..
C     .. Data statements ..
      DATA lastun/unummn/
C     ..
C     .. Executable Statements ..
C
      igtfun = lastun
C
   10 igtfun = igtfun - one
      IF (igtfun.LT.unummn) igtfun = unummx
C
      IF (.NOT. (igtfun.EQ.lastun)) GO TO 20
C     *** all available unit numbers are in use
      igtfun = unummn - one
      GO TO 30

   20 INQUIRE (unit=igtfun,exist=qexist,opened=qopen)
C
      IF (qexist .AND. .NOT.qopen) GO TO 30
C
      GO TO 10

   30 lastun = igtfun
C
C
      RETURN
C
C
      END
      INTEGER FUNCTION lens(string)
C
C
C     .. Scalar Arguments ..
      CHARACTER string* (*)
C     ..
C     .. Local Scalars ..
      INTEGER i,lng
C     ..
C     .. Intrinsic Function ..
      INTRINSIC len
C     ..
C     .. Executable Statements ..
C
      lens = len(string)
      IF (.NOT. (lens.GT.0)) GO TO 40
      lng = lens
      DO 30,i = lng,1,-1
          IF (.NOT. (string(i:i).EQ.' ')) GO TO 10
          lens = lens - 1
          GO TO 20

   10     RETURN

   20     CONTINUE
   30 CONTINUE
   40 RETURN

      END
      SUBROUTINE mvlogf(dfn,dfd,mean,var)
      IMPLICIT DOUBLE PRECISION (a-h,o-p,r-z),INTEGER (i-n),LOGICAL (q)
C**********************************************************************
C
C      SUBROUTINE MVLOGF( DFN, DFD, MEAN, VAR )
C      Mean and Variance of the Log F Distribution
C
C
C                              Arguments
C
C
C     DFN --> Numerator degrees of freedom of the log-F
C                    DOUBLE PRECISION DFN
C
C     DFD --> Denominator degrees of freedom of the log-F
C                    DOUBLE PRECISION DFD
C
C     MEAN <-- Mean of log-F distribution
C                    DOUBLE PRECISION MEAN
C
C     VAR <-- Variance of log-F distribution
C                    DOUBLE PRECISION VAR
C
C**********************************************************************

C     .. Parameters ..
      DOUBLE PRECISION half
      PARAMETER (half=0.5D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION dfd,dfn,mean,var
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION psia(2),psib(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL dpsifn
C     ..
C     .. Executable Statements ..

      a = dfn*half
      b = dfd*half

      CALL dpsifn(a,0,2,2,psia)
      CALL dpsifn(b,0,2,2,psib)

C     At this point, PSIA(1) = -psi(A) + ln(A)
C                    PSIA(2) = psi(1,A)
C     Similarly for B

      mean = psib(1) - psia(1)

      var = psia(2) + psib(2)

      RETURN

      END
      INTEGER FUNCTION obfugt(mssg,qfrmtd)
C
C     Calls QGOBFU
C
C     .. Scalar Arguments ..
      LOGICAL qfrmtd
      CHARACTER mssg* (*)
C     ..
C     .. Local Scalars ..
      LOGICAL qok
C     ..
C     .. External Functions ..
      LOGICAL qgobfu
      EXTERNAL qgobfu
C     ..
C     .. Executable Statements ..
      qok = qgobfu(mssg,qfrmtd,obfugt)
      IF (.NOT. (qok)) STOP ' QGOBFU CALLED FROM OBFUGT'
      RETURN

      END
      SUBROUTINE pause
      IMPLICIT INTEGER (a-p,r-z),LOGICAL (q)
C**********************************************************************
C
C     SUBROUTINE PAUSE
C
C
C                              Function
C
C
C     The message
C         'ENTER RETURN/ENTER TO CONTINUE'
C     is displayed at the terminal and the display is stopped until
C     the user hits the return/enter key.
C
C**********************************************************************
C
C     .. Local Scalars ..
      CHARACTER dum*1
C     ..
C     .. External Functions ..
      CHARACTER strgt*1
      EXTERNAL strgt
C     ..
C     .. Executable Statements ..
      dum = strgt('('' Press Return / Enter to continue:'')',.TRUE.)
      RETURN

      END
      SUBROUTINE prompt
C**********************************************************************
C
C     SUBROUTINE PROMPT
C
C
C                              Function
C
C
C          PROMPT issues a line feed and a '?' prompt for FORTRAN system
C     that do not have a built in prompt feature.
C
C**********************************************************************
C**********************************************************************
C
C     This version is specific to D.E.C. VAX running V.M.S. FORTRAN.
C     Code last modified 87/Feb/20
C
C**********************************************************************
      IMPLICIT INTEGER (a-p,r-z),LOGICAL (q)
C
C
C     QECIN   --  Flag as to whether input is echoed
C     QECALL  --  Flag as to whether entire dialog is echoed
C     ECIN    --  Unit to which input is echoed
C     ECALL   --  Unit to which entire dialog is echoed
C
C
C     GTECUN  --  returns unit numbers associated with echo input and
C                 echo all files (negative numbers if none)
C     .. Local Scalars ..
      INTEGER ecall,ecin,inun,outun
      LOGICAL qecall,qecin
C     ..
C     .. External Subroutines ..
      EXTERNAL gtcuio,gtecun
C     ..
C     .. Executable Statements ..
C
      CALL gtcuio(inun,outun)
C
C     Initialize for echo processing
C
      CALL gtecun(ecin,ecall)
      qecin = ecin .GE. 0
      qecall = ecall .GE. 0
      WRITE (outun,'(''$? '')')
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 10
      WRITE (ecall,9000)

 9000 FORMAT (' ? ')

   10 RETURN
C     - === End Echo Code   ===
      END
      LOGICAL FUNCTION qgobfu(mssg,qfrmtd,unit)
C**********************************************************************
C
C     LOGICAL FUNCTION QGOBFU(FORMSG,QFRMTD,UNIT)
C
C                    Q Get Output B File Unit number
C                     (Machine Dependent Routine)
C
C
C                              Function
C
C
C          Conducts a dialog using current input and output to obtain
C     the name of the output file.  If the file is currently open, the
C     unit number with which it is associated is returned.  If the file
C     is not open, then it is opened and associated with a unit number
C     not in current use and this unit number is returned.
C          The function returns .TRUE. if an output file was obtained,
C     .FALSE. if some problem occurred.  An EOF on current input
C     will cause a value of .FALSE. to be returned if current input
C     is not standard input (checked by equality of unit numbers).
C     Three consequtive errors will cause a value of .FALSE. to
C     be returned regardless of the current input file.
C
C     The dialog is sketched below.  Some details may differ from
C     machine to machine as the set of legal file names differs as does
C     the information that the user may need to specify such a name.
C
C     The message specified by FORMSG is presented the user.
C
C     ENTER FILE NAME
C
C     The program then does an INQUIRE on the file.  If it is open,
C     the unit number with which it is associated will be returned
C     (assuming that the user does no specify REDO on the verify).
C     If the file is not open, the program ascertains whether it exists
C     and in either case opens it. Should there be a problem with either
C     the INQUIRE or OPEN (which is probably indicative of an illegal
C     file specification), the following message is presented the user.
C
C     PROBLEM OPENING FILE <filespec> ---
C     --- R(ETRY) FILE SPECIFICATION OR Q(UIT)
C
C     If the user specifies Q the routine returns .FALSE.
C
C     If the user specifies R, the routine starts again from the
C     beginning of the dialog.
C
C     If the file is opened successfully and did not previously exist,
C     the message below is presented.
C
C     OUTPUT FILE CREATED -- <filespec>
C
C     If the file is opened successfully and did previously exist, the
C     user is presented with a choice
C
C     OUTPUT FILE ALREADY EXISTS -- <filespec> --
C     -- O(VERWRITE) OR A(PPEND) TO END OF FILE
C
C     If O is specified, the file is rewound.  If A is specified,
C     the file is positioned immediately prior to the EOF.
C
C     Finally, the user is given the opportunity to change his mind.
C
C     P(ROCEED) OR R(EDO) FILE SPECIFICATION
C
C     If P is specified, QGOBFU returns .TRUE.
C     If R is specified, the dialog begins over.
C
C
C                              Arguments
C
C
C     QGOBFU <-- .TRUE. if all has gone well and a specified old file
C                has been successfully opened, otherwise .FALSE.
C
C     FORMSG --> A FORTRAN format (first and last characters must be
C                left and right parentheses) which, when invoked,
C                will print a message to the user about the
C                output file that is to be specified.
C                If FORMSG is a character string consisting entirely
C                of blanks, then no message is printed to the user.
C                                            FORMSG is CHARACTER*(*)
C
C     QFRMTD --> If .TRUE., the output file will be opened for FORMATTED
C                writing.
C                If .FALSE., the output file will be opened for
C                UNFORMATTED writing.
C
C     UNIT <-- Defined only if QGOBFU is .TRUE.  The unit number
C              associated with the output file.  This unit
C              number will be the one with which the file
C              was associated if it was open at the call, otherwise
C              it will be associated with an unused unit number which
C              is returned.
C                                                  UNIT is INTEGER
C
C
C                                   NOTE
C
C     This routine is very similar to QGTOFU. The difference is that the
C     user may NOT enter * to specify the terminal to be the output file
C
C----------------------------------------------------------------------
C
C
C     ENTRY QNOBFU(FILNAM)
C
C
C                              Function
C
C
C     This entry inherits the characteristic of being a logical function
C     Returns .TRUE. if the previous invocation successfully opened an
C     output file; in this case, FILNAM is the (system dependent) name
C     of the file that was opened.
C     Returns .FALSE. if the previous invocation did not successfully
C     open an output file; in this case, FILNAM is indeterminate.
C
C
C                              Arguments
C
C
C     FILNAM <-- The name of the file opened for output on the previous
C                invocation if it was successful; otherwise
C                indeterminate.
C
C
C----------------------------------------------------------------------
C
C
C                    ===== IBM CMS Specific Information =====
C
C     ENTRY QOPOBF(OPTSTR)
C
C
C                              Function
C
C
C     This entry point allows the specification of the format of the
C     output file.  If the entry is not called before the call to the
C     main routine, the output data file will be opened as a file of
C     fixed-length 80 character records.
C
C     OPTSTR is a string that contains options for a FILEDEF command.
C     The string should not begin with a left parenthesis.  The options
C     specified are in effect for only the next invocation of QGOBFU.
C                                             OPTSTR is CHARACTER*60
C
C     The value returned by QOPOBF is always .TRUE. and should be
C     ignored.
C
C**********************************************************************
C
C *** Variables
C
C
C
C
C  ** Parameters and constants
C
C   * Processor-dependent parameters
C
C     FLNMLM  --  limit of length of a file name
C
C   * CONSTANTS
C
C     ONE     --  mnemonically-named integer constant
C
C     QFALSE  --  mnemonically-named logical constant
C
C     QTRUE   --  mnemonically-named logical constant
C
C     ZERO    --  mnemonically-named integer constant
C
C     Z*****  --  assorted character constants
C
C
C   * Parameters peculiar to this application
C
C     ERRMX   --  maximum number of read errors allowed before automatic
C                 failure occurs
C
C
C  ** Arguments
C
C     MSSG    --> buffer to hold format specification for writing
C                 prompt message to user.  If value is blank,
C                 no prompt message will be written.
C
C     QFRMTD  --> formatted-file flag.  Value is .TRUE. if file to be
C                 opened formatted, .FALSE. if to be opened unformatted
C
C     UNIT   <--  Fortran unit number assigned to file.  Value is
C                 significant only if result of function is .TRUE.
C
C
C  ** Other significant variables
C
C     ECALL   --  Unit to which entire dialog is echoed
C
C     ECIN    --  Unit to which input is echoed
C
C     ERR     --  variable to hold error code returned by I/O statement
C
C     FILEID  --  buffer to hold file id
C
C     FILNAM  --  dummy argument for entry point that returns name
C                 of last file opened by main routine
C
C     FRSTNB  --  position in buffer fileid (q.v.) of first non-blank
C                 character
C
C     INNM    --  Fortran unit number of current input unit
C
C     LASTNB  --  position in buffer fileid (q.v.) of last non-blank
C                 character
C
C     MSG     --> buffer to hold format specification for messages
C                 to be written to user by called routines.
C
C     NEWNUM  --  unassigned unit number
C
C     OUTNM   --  Fortran unit number of current output unit
C
C     QECALL  --  Flag as to whether entire dialog is echoed
C
C     QECIN   --  Flag as to whether input is echoed
C
C     QEXIST  --  does-file-exist flag
C
C     QOPEN   --  is-file-open flag
C
C     QPASS1  --  does-new-unit-number-need-to-be-generated flag
C
C     QPSTAT  --  Value returned in previous invocation of QGOBFU
C
C     QRETRY  --  does-user-want-to-retry-file-specification flag
C
C     STDIFU  --  Fortran unit number of standard input unit
C
C     STDOFU  --  Fortran unit number of standard output unitT
C
C
C  ** Miscellaneous variables
C     The following variables are used for such purposes as indices and
C     temporary storage.  Each one is intended to be significant only
C     within a well-defined block of code (although it may appear in
C     more than one block), and its meaning and use should be apparent
C     from its name and context.
C
C     ANSWER  --
C     FILERR  --
C     I       --
C
C
C
C *** Functions and subroutines
C
C
C  ** Library subprograms
C     Detailed information on each of the following subprograms can be
C     found in the reference for its respective library.
C
C   * COMPLIB
C
C     GTCUIO  --  returns unit numbers associated with current input
C                 and output units
C
C     GTECUN  --  returns unit numbers associated with echo input and
C                 echo all files (negative numbers if none)
C
C     GTSTIO  --  returns unit numbers associated with standard inputT
C                 and output units
C
C     IGTFUN  --  returns a Fortran unit number that is not currently j
C                 in use
C     .. Parameters ..
      INTEGER flnmlm
      PARAMETER (flnmlm=500)
      INTEGER one
      PARAMETER (one=1)
      LOGICAL qfalse
      PARAMETER (qfalse=.FALSE.)
      LOGICAL qtrue
      PARAMETER (qtrue=.TRUE.)
      INTEGER zero
      PARAMETER (zero=0)
      CHARACTER*(1) zblank
      PARAMETER (zblank=' ')
      INTEGER errmx
      PARAMETER (errmx=3)
C     ..
C     .. Scalar Arguments ..
      INTEGER unit
      LOGICAL qfrmtd
      CHARACTER mssg* (*),filnam* (*)
C     ..
C     .. Local Scalars ..
      INTEGER answer,ecall,ecin,err,filerr,frstnb,i,i99954,i99964,
     +        i99967,i99970,i99987,i99992,i99995,i99999,innm,lastnb,
     +        newnum,outnm,stdifu,stdofu
      LOGICAL qecall,qecin,qexist,qok,qopen,qpass1,qpstat,qretry
      CHARACTER msg*100,fileid* (flnmlm)
C     ..
C     .. External Functions ..
      INTEGER igtfun
      LOGICAL qgtchr,qgtstr
      EXTERNAL igtfun,qgtchr,qgtstr
C     ..
C     .. External Subroutines ..
      EXTERNAL gtcuio,gtecun,gtstio
C     ..
C     .. Entry Points ..
      LOGICAL qnobfu
C     ..
C     .. Save statement ..
      SAVE fileid,qpstat
C     ..
C     .. Data statements ..
C
C
C**********************************************************************
      DATA fileid/' '/
      DATA qpstat/qfalse/
C     ..
C     .. Executable Statements ..
C
C
C
C *** Body of routine
C
C
C     GET-OUTPUT-FILE-UNIT-NUMBER
      ASSIGN 10 TO i99999
      GO TO 160
C
   10 RETURN
C
C
C
      ENTRY qnobfu(filnam)
C
C     Entry point to return name of last file opened by main routine
C
      IF (.NOT. (qpstat)) GO TO 20
      filnam = fileid
      GO TO 30

   20 filnam = ' '
   30 qnobfu = qpstat
      RETURN
C
C
      STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO CHECK-OPEN-ERROR
C
C     ERRMX   -->
C     FILEID  -->
C     FRSTNB  -->
C     INNM    -->
C     LASTNB  -->
C     ONE     -->
C     OUTNM   -->
C     QTRUE   -->
C     STDIFU  -->
C     ZERO    -->
C     ZQ      -->
C     ZR      -->
C     *scrt*  --  ANSWER, ERR
C     QRETRY <--
C     *more*  --  other variables may be used in invoked procedures
C
C
   40 WRITE (outnm,FMT='(T2,''Problem opening file '',A,'' --'')')
     +  fileid(frstnb:lastnb)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 50
      WRITE (ecall,FMT='(T2,''PROBLEM OPENING FILE '',A,'' --'')')
     +  fileid(frstnb:lastnb)
   50 msg = '(T5,''-- R(etry file specification) or Q(uit)?'')'
C     - === End Echo Code   ===
C
C
      qok = qgtchr(msg,'RQ',qtrue,answer)
C
      IF (qok) GO TO 70
C     RETURN-FALSE
      ASSIGN 60 TO i99992
      GO TO 470
C
   60 CONTINUE
   70 IF (.NOT. (answer.EQ.1)) GO TO 80
      qretry = qtrue
      GO TO 100
C     RETURN-FALSE
   80 ASSIGN 90 TO i99992
      GO TO 470
C
   90 CONTINUE
  100 GO TO i99995
C     TO FIND-FILE-NAME
C
C     FILEID  -->
C     FLNMLM  -->
C     ONE     -->
C     ZBLANK  -->
C     ZERO    -->
C     *scrt*  --  I
C     FRSTNB <--
C     LASTNB <--
C
C
  110 frstnb = zero
      i = zero
C
  120 IF (i.GE.flnmlm) GO TO 140
C
      i = i + one
C
      IF (.NOT. (fileid(i:i).NE.zblank)) GO TO 130
      frstnb = i
      GO TO 140

  130 GO TO 120
C
  140 lastnb = frstnb
C
      DO 150 i = (frstnb+one),flnmlm
          IF (fileid(i:i).NE.zblank) lastnb = i
  150 CONTINUE
      GO TO i99987
C
C     TO GET-OUTPUT-FILE-UNIT-NUMBER
C
C     ERRMX   -->
C     MSSG    -->
C     ONE     -->
C     QFALSE  -->
C     QTRUE   -->
C     ZBLANK  -->
C     ZERO    -->
C     *scrt*  --  ERR,  FILERR, FILEID, FRSTNB, INNM, LASTNB,
C     *scrt*  --  OUTNM, QPASS1, QRETRY, STDIFU, STDOFU
C     QGOBFU <--
C     UNIT   <--
C     *more*  --  other variables may be used in invoked procedures
C
C
  160 CALL gtcuio(innm,outnm)
C
      CALL gtstio(stdifu,stdofu)
C
C
C     Initialize for echo processing
C
      CALL gtecun(ecin,ecall)
      qecin = ecin .GE. 0
      qecall = ecall .GE. 0
C
      qpass1 = qtrue
C
      filerr = zero
C     *** until get a file or fail or decide to bail out ...
C
  170 IF (.NOT. (filerr.GE.errmx)) GO TO 190
C     RETURN-FALSE
      ASSIGN 180 TO i99992
      GO TO 470
C     *** too many failed attempts at files -- return as .FALSE.
C
C     *** get a file name or fail ...
C
  180 CONTINUE
  190 IF (.NOT. (mssg.NE.zblank)) GO TO 210
      WRITE (outnm,FMT=mssg)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 200
      WRITE (ecall,FMT=mssg)
C     - === End Echo Code   ===
  200 CONTINUE
  210 msg = '('' Specify file name:'')'
      qok = qgtstr(msg,fileid,qfalse)
C
      IF (qok) GO TO 230
C     RETURN-FALSE
      ASSIGN 220 TO i99992
      GO TO 470
C
  220 CONTINUE
C     FIND-FILE-NAME
  230 ASSIGN 240 TO i99987
      GO TO 110
C
C     *** check specified file
  240 INQUIRE (file=fileid(frstnb:lastnb),exist=qexist,opened=qopen,
     +        number=unit,iostat=err)
C
C
      IF (.NOT. (err.NE.zero)) GO TO 280
C     *** file does not exist
C
C     OPEN-FILE
      ASSIGN 250 TO i99970
      GO TO 420

  250 IF (.NOT. (qretry)) GO TO 260
      filerr = filerr + one
      GO TO 170
C     VERIFY-FILE
  260 ASSIGN 270 TO i99967
      GO TO 490
C
C
  270 GO TO 370

  280 IF (.NOT. (qopen)) GO TO 300
C     VERIFY-FILE-EXISTING
      ASSIGN 290 TO i99964
      GO TO 530

  290 GO TO 370
C
C     *** FILE ALREADY OPEN

C
C     OPEN-FILE
  300 ASSIGN 310 TO i99970
      GO TO 420

  310 IF (.NOT. (qretry)) GO TO 320
      filerr = filerr + one
      GO TO 170

  320 IF (.NOT. (qexist)) GO TO 340
C     VERIFY-FILE-EXISTING
      ASSIGN 330 TO i99964
      GO TO 530

  330 GO TO 360
C
C     VERIFY-FILE
  340 ASSIGN 350 TO i99967
      GO TO 490
C
  350 CONTINUE
C
  360 CONTINUE
  370 IF (.NOT. (qretry)) GO TO 380
C
      filerr = filerr + one
      GO TO 170

      GO TO 400
C     RETURN-TRUE
  380 ASSIGN 390 TO i99954
      GO TO 480
C     *** attempt successful
C
  390 CONTINUE
  400 GO TO 170

  410 GO TO i99999
C
C     TO OPEN-FILE
C
C     FILEID  -->
C     FRSTNB  -->
C     LASTNB  -->
C     NEWNUM  -->
C     QFALSE  -->
C     QFRMTD  -->
C     QPASS1  -->
C     UNIT    -->
C     ZERO    -->
C     *scrt*  --  ERR
C     NEWNUM <--
C     QPASS1 <--
C     QRETRY <--
C     UNIT   <--
C     *more*  --  other variables may be used in invoked procedures
C
C
  420 IF (.NOT. (qpass1)) GO TO 430
      newnum = igtfun()
      qpass1 = qfalse
  430 unit = newnum
C
      IF (qfrmtd) THEN
          OPEN (unit=unit,file=fileid(frstnb:lastnb),form='FORMATTED',
     +         iostat=err)

      ELSE
          OPEN (unit=unit,file=fileid(frstnb:lastnb),form='UNFORMATTED',
     +         iostat=err)
      END IF
C
      IF (.NOT. (err.EQ.zero)) GO TO 440
      qretry = qfalse
      GO TO 460
C     CHECK-OPEN-ERROR
  440 ASSIGN 450 TO i99995
      GO TO 40
C
  450 CONTINUE
  460 GO TO i99970
C     TO RETURN-FALSE
  470 qgobfu = qfalse
      qpstat = qfalse
      RETURN

      GO TO i99992
C     TO RETURN-TRUE
  480 qgobfu = qtrue
      qpstat = qtrue
      RETURN

      GO TO i99954
C     TO VERIFY-FILE
C
C     ERRMX   -->
C     FILEID  -->
C     FRSTNB  -->
C     INNM    -->
C     LASTNB  -->
C     ONE     -->
C     OUTNM   -->
C     QFALSE  -->
C     QTRUE   -->
C     UNIT    -->
C     ZERO    -->
C     ZP      -->
C     ZR      -->
C     *scrt*  --  ANSWER, ERR
C     QRETRY <--
C     *more*  --  other variables may be used in invoked procedures
C
C
  490 WRITE (outnm,FMT='(T2,''Created file '',A,'' --'')')
     +  fileid(frstnb:lastnb)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 500
      WRITE (ecall,FMT='(T2,''CREATED FILE '',A,'' --'')')
     +  fileid(frstnb:lastnb)
  500 msg = '(T5,''-- P(roceed) or R(etry file specification)?'')'
C     - === End Echo Code   ===
C
      qok = qgtchr(msg,'PR',qtrue,answer)
C
      IF (qok) GO TO 520
C     RETURN-FALSE
      ASSIGN 510 TO i99992
      GO TO 470
C
  510 CONTINUE
C
  520 IF ((1).EQ. (answer)) THEN
          qretry = qfalse
C
      ELSE IF ((2).EQ. (answer)) THEN
          qretry = qtrue
          CLOSE (unit)
      END IF
C
      GO TO i99967
C
C     TO VERIFY-FILE-EXISTING
C
C     ERRMX   -->
C     FILEID  -->
C     FRSTNB  -->
C     INNM    -->
C     LASTNB  -->
C     ONE     -->
C     OUTNM   -->
C     QFALSE  -->
C     QFRMTD  -->
C     QTRUE   -->
C     UNIT    -->
C     ZA      -->
C     ZERO    -->
C     ZO      -->
C     ZQ      -->
C     ZR      -->
C     *scrt*  --  ANSWER, ERR
C     QRETRY <--
C     *more*  --  other variables may be used in invoked procedures
C
C
  530 WRITE (outnm,FMT='(T2,''File '',A,'' already exists --'')')
     +  fileid(frstnb:lastnb)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 540
      WRITE (ecall,FMT='(T2,''FILE '',A,'' ALREADY EXISTS --'')')
     +  fileid(frstnb:lastnb)
  540 qok = qgtchr(
     +'('' Enter choice of action''/t5,''A(ppend)''/t5,''O(verwrite)''/
     +t5,''R(etry)''/t5,''Q(uit)'')','AORQ',qtrue,answer)
C     - === End Echo Code   ===
C
      IF (qok) GO TO 560
C     RETURN-FALSE
      ASSIGN 550 TO i99992
      GO TO 470
C
  550 CONTINUE
C
  560 IF ((1).NE. (answer)) GO TO 630
C     *** append -- skip to end of file
      IF (qfrmtd) THEN
          GO TO 580

  570     IF (.NOT. (err.EQ.zero)) GO TO 590
  580     READ (unit,FMT='(A)',iostat=err)
          GO TO 570

  590     CONTINUE

      ELSE
          GO TO 610

  600     IF (.NOT. (err.EQ.zero)) GO TO 620
  610     READ (unit,iostat=err)
          GO TO 600

  620     CONTINUE
      END IF

      BACKSPACE unit
      qretry = qfalse
      GO TO 680

  630 IF ((2).NE. (answer)) GO TO 640
C
C     *** overwrite
      REWIND unit
      qretry = qfalse
      GO TO 680

  640 IF ((3).NE. (answer)) GO TO 650
C
C     *** retry
      CLOSE (unit)
      qretry = qtrue
      GO TO 680

  650 IF ((4).NE. (answer)) GO TO 670
C     RETURN-FALSE
      ASSIGN 660 TO i99992
      GO TO 470

  660 GO TO 680
C
C     *** quit
C
  670 CONTINUE
  680 GO TO i99964
C
C
C
      END
      LOGICAL FUNCTION qgtchr(mssg,chrlst,qupcas,iwhich)
C**********************************************************************
C
C     LOGICAL FUNCTION QGTCHR(FORMSG,CHRLST,QUPCAS,IWHICH)
C
C                    Q GeT ChaRacter
C
C
C                              Function
C
C
C          If FORMSG is non-blank, writes it to current output then
C     reads a line from current input and examines the first
C     non-blank character for being a member of CHRLST.  If this
C     character is in CHRLST, returns .TRUE. and sets IWHICH to
C     the ordinal position of the character in CHRLST.
C          If there are three unsuccessful attempts to obtain a
C     character in CHRLST, QGTCHR will return .FALSE. QGTCHR will
C     also return .FALSE. if and EOF is enountered on current input
C     and current input is not standard input (as determined by
C     comparing unit numbers).
C
C
C                              Note
C
C
C     Two error messages can be issued to current output as a result of
C     the use of this routine.
C
C          A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN
C     is issued by QGTSTR, which is called by this routine.
C
C          THE LEGAL FIRST NON-BLANK CHARACTER MUST BE ONE OF
C               THE FOLLOWING - TRY AGAIN
C          (followed by a list of legal characters)
C     is issued if the first non-blank character is not in CHRLST.
C
C     FORMSG is rewritten after the above error messages.
C
C
C                              Arguments
C
C
C     QGTCHR <-- .TRUE. if the first non-blank character obtained is a
C                member of CHRLST (even if it takes three attempts).
C               .FALSE. if an EOF is encountered on current input and
C               current input is not standard input or if three
C               attempts to obtain a character from CHRLST fail.
C                                             QGTCHR is LOGICAL
C
C     FORMSG --> A FORTRAN format (first and last characters must be
C                left and right parentheses) which, when invoked,
C                will print a message to the user about the
C              character that is to be obtained.
C              If FORMSG is a character string consisting entirely of
C              blanks, then no message is printed to the user.
C                                            FORMSG is CHARACTER*(*)
C
C     CHRLST --> A list of characters against which the first non-blank
C               character from current input is matched.
C                                            CHRLST is CHARACTER*(*)
C
C     QUPCAS --> If .TRUE. the character obtained from current input is
C                translated to upper case if it is a lower case letter
C                before the match to CHRLST is done.
C                If .FALSE. no such translation is performed.
C                                                  QUPCAS is LOGICAL
C
C     IWHICH <-- Defined only if QGTCHR is .TRUE.  The ordinal position
C               of the character obtained in CHRLST.
C                                                  IWHICH is INTEGER
C
C**********************************************************************
C
C
C     QECIN   --  Flag as to whether input is echoed
C     QECALL  --  Flag as to whether entire dialog is echoed
C     ECIN    --  Unit to which input is echoed
C     ECALL   --  Unit to which entire dialog is echoed
C
C
C     GTECUN  --  returns unit numbers associated with echo input and
C                 echo all files (negative numbers if none)
C
C
C
C     .. Scalar Arguments ..
      INTEGER iwhich
      LOGICAL qupcas
      CHARACTER mssg* (*),chrlst* (*)
C     ..
C     .. Local Scalars ..
      INTEGER ecall,ecin,i,i99996,icuin,icuot,itries,ivalue,lenstr
      LOGICAL qecall,qecin,qstr
      CHARACTER char*1,string*80
C     ..
C     .. External Functions ..
      LOGICAL qgtstr
      CHARACTER trnchr*1
      EXTERNAL qgtstr,trnchr
C     ..
C     .. External Subroutines ..
      EXTERNAL gtcuio,gtecun
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC index,len
C     ..
C     .. Executable Statements ..
C
      CALL gtcuio(icuin,icuot)
C
C     Initialize for echo processing
C
      CALL gtecun(ecin,ecall)
      qecin = ecin .GE. 0
      qecall = ecall .GE. 0
      itries = 0
      qgtchr = .FALSE.
C
C     loop until valid character read or processing terminated
C
   10 IF (qgtchr) GO TO 130
C
C     increment and test number of tries
C
      itries = itries + 1
      IF (.NOT. (itries.GT.3)) GO TO 30
C     TERMINATE-PROCESSING
      ASSIGN 20 TO i99996
      GO TO 140
C
C     get input string
C
   20 CONTINUE
   30 qstr = qgtstr(mssg,string,.FALSE.)
C
      IF (.NOT. (qstr)) RETURN
C
C     find first non-blank character
C
      lenstr = len(string)
      i = 1
      char = ' '
   40 IF (.NOT. ((i.LE.lenstr).AND. (char.EQ.' '))) GO TO 60
      IF (.NOT. (string(i:i).NE.' ')) GO TO 50
      char = string(i:i)
   50 i = i + 1
      GO TO 40

   60 IF (qupcas) char = trnchr(char)
C
C     check char for membership in character list CHRLST
C
C
      iwhich = index(chrlst,char)
      IF (.NOT. (iwhich.NE.0)) GO TO 70
      qgtchr = .TRUE.
      GO TO 120

   70 WRITE (icuot,*) ' FIRST CHARACTER MUST BE ONE OF THE FOLLOWING:'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 80
      WRITE (ecall,*) ' FIRST CHARACTER MUST BE ONE OF THE FOLLOWING:'
   80 DO 100,i = 1,len(chrlst)
C     - === End Echo Code   ===
          WRITE (icuot,*) '     ',chrlst
C     + === Begin Echo Code ===
          IF (.NOT. (qecall)) GO TO 90
          WRITE (ecall,*) '     ',chrlst
   90     CONTINUE
  100 CONTINUE
C     - === End Echo Code   ===
      WRITE (icuot,*) ' PLEASE TRY AGAIN'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 110
      WRITE (ecall,*) ' PLEASE TRY AGAIN'
C     - === End Echo Code   ===
  110 CONTINUE
  120 GO TO 10
C
  130 RETURN
C
      STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO TERMINATE-PROCESSING
C
C     WHOOPS, LOOKS LIKE WE'VE BOMBED OUT!
C
  140 itries = itries - 1
      qgtchr = .FALSE.
      ivalue = 0
      WRITE (icuot,9000) itries
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 150
      WRITE (ecall,9000) itries
  150 CONTINUE

 9000 FORMAT (' CHARACTER READ WAS UNSUCCESSFUL AFTER ',I2,' TRIES.')
C     - === End Echo Code   ===
      RETURN

      GO TO i99996
C
      END
      LOGICAL FUNCTION qgtrl(mssg,low,hi,value)
C**********************************************************************
C
C     LOGICAL FUNCTION QGTRL(FORMSG,RLOW,RHI,RVALUE)
C
C                         Q GeT ReaL
C
C
C                              Function
C
C
C          If FORMSG is non-blank, prints it to current output, then
C     reads a real value between RLOW and RHI (inclusive) from
C     the current input unit and returns the value in RVALUE.  QGTRL
C     returns .TRUE. if this was successfully accomplished, else
C     it returns .FALSE.
C          If there are three unsuccessful attempts to obtain the
C     number in the range specified, QGTRL returns .FALSE.
C     QGTRL will also return .FALSE. if an EOF is encountered on
C     current input and current input is not standard input
C     (as determined by comparing unit numbers).
C
C
C                              Note
C
C
C     Three error messages can be issued to current output as a result
C     of the use of this routine.
C
C          A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN
C     is issued by this routine if a blank line or EOF is
C     encountered.
C
C          THE STRING ENTERED IS NOT A LEGAL REAL VALUE - TRY AGAIN
C     is issued for the obvious reason.
C
C          NUMBER MUST BE BETWEEN <lo> AND <hi> - TRY AGAIN
C     is issued if a legal but out of range value is encountered.
C
C          FORMSG is re-written to current output after either of the
C     previous error messages.
C
C
C                              Arguments
C
C
C     QGTRL <-- .TRUE. if a legal real in the specified range
C               was obtained from current input, else .FALSE.
C                                                  QGTRL is LOGICAL
C
C     FORMSG --> A FORTRAN format (first and last characters must be
C                left and right parentheses) which, when invoked,
C                will print a message to the user about the
C                real number that is to be obtained.
C                If FORMSG is a character string consisting entirely
C                of blanks, then no message is printed to the user.
C                                       FORMSG is CHARACTER*(*)
C
C     RLOW --> The lowest legal value for the number obtained.
C                                                  RLOW is REAL
C
C     RHI --> The highest legal value for the number obtained.
C                                                  RHI is REAL
C
C     RVALUE <-- Defined only if QGTRL is .TRUE.  The value of the
C               real number obtained from current input.
C                                                  RVALUE is REAL
C
C**********************************************************************
C
C
C     QECIN   --  Flag as to whether input is echoed
C     QECALL  --  Flag as to whether entire dialog is echoed
C     ECIN    --  Unit to which input is echoed
C     ECALL   --  Unit to which entire dialog is echoed
C
C
C     GTECUN  --  returns unit numbers associated with echo input and
C                 echo all files (negative numbers if none)
C
C
C
C     .. Scalar Arguments ..
      REAL hi,low,value
      CHARACTER mssg* (*)
C     ..
C     .. Local Scalars ..
      INTEGER ecall,ecin,i99996,icuin,icuout,ios,istin,istout,itries
      LOGICAL qecall,qecin
C     ..
C     .. External Subroutines ..
      EXTERNAL gtcuio,gtecun,gtstio,prompt
C     ..
C     .. Executable Statements ..
C
C     call routines to get current and standard I/O unit numbers
C
      CALL gtcuio(icuin,icuout)
      CALL gtstio(istin,istout)
C
C     Initialize for echo processing
C
      CALL gtecun(ecin,ecall)
      qecin = ecin .GE. 0
      qecall = ecall .GE. 0
      itries = 0
      qgtrl = .FALSE.
C
C     loop until number read or processing terminated
C
   10 IF (qgtrl) GO TO 190
C
C     increment and test number of tries
C
      itries = itries + 1
      IF (.NOT. (itries.GT.3)) GO TO 30
C     TERMINATE-PROCESSING
      ASSIGN 20 TO i99996
      GO TO 200
C
C     print message if there is one and read input
C
   20 CONTINUE
   30 IF (.NOT. (mssg.NE.' ')) GO TO 50
      WRITE (icuout,mssg)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 40
      WRITE (ecall,mssg)
C     - === End Echo Code   ===
   40 CONTINUE
   50 CALL prompt
      READ (icuin,*,iostat=ios) value
C
C
C     if number read, verify against acceptable range
C
      IF (.NOT. (ios.EQ.0)) GO TO 110
C     + === Begin Echo Code ===
      IF (.NOT. (qecin)) GO TO 60
      WRITE (ecin,*) value
   60 IF (.NOT. (qecall)) GO TO 70
      WRITE (ecall,*) value
   70 IF (.NOT. (value.LT.low.OR.value.GT.hi)) GO TO 90
C     - === End Echo Code   ===
C
C     if number was read correctly but was out of range...
C
      WRITE (icuout,9000) low,hi
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 80
      WRITE (icuout,9000) low,hi
   80 CONTINUE

 9000 FORMAT (' NUMBER MUST BE BETWEEN ',G9.3,' AND ',G9.3,' -',
     +       ' TRY AGAIN')
C     - === End Echo Code   ===
      GO TO 100

   90 qgtrl = .TRUE.
C
  100 GO TO 180

  110 IF (.NOT. (ios.GT.0)) GO TO 130
C
C     if an error is encountered... (IOS > 0)
C
      WRITE (icuout,*)
     +  ' AN ERROR WAS ENCOUNTERED DURING READ - TRY AGAIN'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 120
      WRITE (ecall,*)
     +  ' AN ERROR WAS ENCOUNTERED DURING READ - TRY AGAIN'
  120 GO TO 180
C     - === End Echo Code   ===
  130 IF (.NOT. (icuin.NE.istin)) GO TO 160
C
C     if end-of-file is encountered... (IOS < 0)
C
C     if end-of-file occurred on non-standard input, then give up
C
      WRITE (icuout,*) ' END-OF-FILE ENCOUNTERED ON NON-STANDARD INPUT'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 140
      WRITE (ecall,*) ' END-OF-FILE ENCOUNTERED ON NON-STANDARD INPUT'
C     TERMINATE-PROCESSING
  140 ASSIGN 150 TO i99996
      GO TO 200
C     - === End Echo Code   ===
  150 GO TO 180
C
C     if end-of-file occurred on standard input, then rewind input
C
  160 WRITE (icuout,*) ' A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 170
      WRITE (ecall,*) ' A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN'
  170 REWIND (icuin)
C     - === End Echo Code   ===
  180 GO TO 10

  190 RETURN
C
C     at this point, we know that IOS = 0 (no error on read)
C
      STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO TERMINATE-PROCESSING
C
C     WHOOPS, LOOKS LIKE WE'VE BOMBED OUT!
C
  200 itries = itries - 1
      qgtrl = .FALSE.
      value = 0.0
      WRITE (icuout,9010) itries
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 210
      WRITE (ecall,9010) itries
  210 CONTINUE

 9010 FORMAT (' REAL NUMBER READ WAS UNSUCCESSFUL AFTER ',I2,' TRIES.')
C     - === End Echo Code   ===
      RETURN

      GO TO i99996
C
      END
      LOGICAL FUNCTION qgtstr(mssg,string,qnulok)
C**********************************************************************
C
C     LOGICAL FUNCTION QGTSTR(FORMSG,STRING,QNULOK)
C
C            Q GeT STRing ignoring lines that begin with #
C
C
C
C                              Function
C
C
C          If FORMSG is non-blank, prints it to current output, then
C     reads a line from the current input file and places its
C     contents in STRING.  Returns .TRUE. if the line was successfully
C     read.  It returns .FALSE. if (1) an EOF is encountered and
C     current input is not standard input; (2) an error was encountered
C     in the read; or (3) QNULOK is .FALSE. indicating that a null line
C     is not a legal input, and three successive null lines are obtained
C
C          If QNULOK is .FALSE. indicating that a null line is not a
C     legal input, the message
C          A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN
C     is written to current output and a new line is obtained.
C
C
C                              Note
C
C
C          If an EOF is encountered on standard input, then standard
C     input is rewound.  The EOF is treated as a null
C     line, as is a line consisting only of blank characters, as is
C     a line with '#' as the first non-blank character.
C     If QNULOK is .TRUE. then a null line will cause STRING to
C     be returned filled with all blanks.
C
C
C                              Arguments
C
C
C     FORMSG --> A FORTRAN format (first and last characters must be
C                left and right parentheses) which, when invoked,
C                will print a message to the user about the
C                string that is to be obtained.
C                If FORMSG is a character string consisting entirely
C                of blanks, then no message is printed to the user.
C                                          FORMSG is CHARACTER*(*)
C
C     QGTSTR <-- .TRUE. if QNULOK is .FALSE. and a non-null record was
C                read from current input or if QNULOK is .TRUE. and
C                a null or non-null record was read.  The reasons for
C                returning .FALSE. are listed above.
C                                                QGTSTR is LOGICAL
C
C     STRING <-- Defined only if QGTSTR is .TRUE.  The character string
C                read from current input.
C                                          STRING is CHARACTER*(*)
C
C**********************************************************************
C
C
C     QECIN   --  Flag as to whether input is echoed
C     QECALL  --  Flag as to whether entire dialog is echoed
C     ECIN    --  Unit to which input is echoed
C     ECALL   --  Unit to which entire dialog is echoed
C
C
C     GTECUN  --  returns unit numbers associated with echo input and
C                 echo all files (negative numbers if none)
C
C
C     PROMPT  --  writes prompt character as appropriate for system
C
C
C
C     .. Scalar Arguments ..
      LOGICAL qnulok
      CHARACTER mssg* (*),string* (*)
C     ..
C     .. Local Scalars ..
      INTEGER ecall,ecin,i,i99996,icuin,icuout,ios,istin,istout,itries,
     +        length
      LOGICAL qcomnt,qdone,qecall,qecin
C     ..
C     .. External Subroutines ..
      EXTERNAL gtcuio,gtecun,gtstio,prompt
C     ..
C     .. External Functions ..
      INTEGER lens
      EXTERNAL lens
C     ..
C     .. Executable Statements ..
C
C     call routines to get current and standard I/O unit numbers
C
      CALL gtcuio(icuin,icuout)
      CALL gtstio(istin,istout)
C
C     Initialize for echo processing
C
      CALL gtecun(ecin,ecall)
      qecin = ecin .GE. 0
      qecall = ecall .GE. 0
      itries = 0
      qgtstr = .FALSE.
C
C     loop until string read or processing terminated
C
   10 IF (qgtstr) GO TO 390
C
C     increment and test number of tries
C
      itries = itries + 1
      IF (.NOT. (itries.GT.3)) GO TO 30
C     TERMINATE-PROCESSING
      ASSIGN 20 TO i99996
      GO TO 400

   20 RETURN

   30 IF (.NOT. (mssg.NE.' ')) GO TO 50
C
C     print message if there is one
C
      WRITE (icuout,FMT=mssg)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 40
      WRITE (ecall,FMT=mssg)
C     - === End Echo Code   ===
   40 CONTINUE
   50 IF (icuin.EQ.istin) CALL prompt
C
C     Write prompt character
C
C
C     loop until we get a line without a '#' in the first column
C
      GO TO 70

   60 IF (qdone) GO TO 380
   70 qdone = .TRUE.
C
C     read input
C
      READ (icuin,FMT='(A)',iostat=ios) string
C
C
C     if string read, set QGTSTR
C
      IF (.NOT. (ios.EQ.0)) GO TO 260
C     + === Begin Echo Code ===
      IF (.NOT. (qecin)) GO TO 80
      WRITE (ecin,FMT='(A)') string
   80 IF (.NOT. (qecall)) GO TO 90
      WRITE (ecall,FMT='(A)') string
   90 IF (.NOT. (string(1:1).EQ.'#')) GO TO 100
      qdone = .FALSE.
      GO TO 250
C     - === End Echo Code   ===
C
C     if STRING is not already blank, blank eveything following #
C       including the #, if it already is blank it isn't a comment line
C
  100 IF (.NOT. (string.NE.' ')) GO TO 170
      length = lens(string)
      i = 1
      GO TO 120

  110 i = i + 1
  120 IF (i.GT.length) GO TO 140
      IF (.NOT. (string(i:i).EQ.'#')) GO TO 130
      string(i:length) = ' '
      GO TO 140

  130 GO TO 110

  140 IF (.NOT. (string.EQ.' ')) GO TO 150
      qcomnt = .TRUE.
      GO TO 160

  150 qcomnt = .FALSE.
  160 GO TO 180

  170 qcomnt = .FALSE.
C
C     if string was read correctly but was blank or consisted of nothing
C     but a comment when not allowed. Both cases will make STRING = ' '
C
  180 IF (.NOT. ((.NOT.qnulok).AND. (string.EQ.' '))) GO TO 230
      IF (.NOT. (qcomnt)) GO TO 200
      WRITE (icuout,FMT=9000)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 190
      WRITE (icuout,FMT=9000)
  190 CONTINUE

 9000 FORMAT (' COMMENT LINE NOT ALLOWED - TRY AGAIN')
C     - === End Echo Code   ===
      GO TO 220
C
  200 WRITE (icuout,FMT=9010)
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 210
      WRITE (icuout,FMT=9010)
  210 CONTINUE

 9010 FORMAT (' BLANK LINE NOT ALLOWED - TRY AGAIN')
C     - === End Echo Code   ===
  220 GO TO 240

  230 qgtstr = .TRUE.
C
  240 CONTINUE
  250 GO TO 370

  260 IF (.NOT. (ios.GT.0)) GO TO 290
C
C     if an error is encountered...
C
      WRITE (icuout,FMT=*) ' AN ERROR WAS ENCOUNTERED DURING READ'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 270
      WRITE (ecall,FMT=*) ' AN ERROR WAS ENCOUNTERED DURING READ'
C     TERMINATE-PROCESSING
  270 ASSIGN 280 TO i99996
      GO TO 400
C     - === End Echo Code   ===
  280 RETURN

      GO TO 370

  290 IF (.NOT. (icuin.NE.istin)) GO TO 320
C
C     if end-of-file is encountered... (IOS < 0)
C
C     if end-of-file occurred on non-standard input, then give up
C
      WRITE (icuout,FMT=*)
     +  ' END-OF-FILE ENCOUNTERED ON NON-STANDARD INPUT'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 300
      WRITE (ecall,FMT=*)
     +  ' END-OF-FILE ENCOUNTERED ON NON-STANDARD INPUT'
C     TERMINATE-PROCESSING
  300 ASSIGN 310 TO i99996
      GO TO 400
C     - === End Echo Code   ===
  310 RETURN

      GO TO 370

  320 IF (.NOT. (qnulok)) GO TO 350
C
C     if end-of-file occurred on standard input, then rewind input
C
C     if a null line is acceptable as input
C
      IF (icuin.EQ.istin) REWIND (icuin)
      string = ' '
C     + === Begin Echo Code ===
      IF (.NOT. (qecin)) GO TO 330
      WRITE (ecin,FMT='(A)') string
  330 IF (.NOT. (qecall)) GO TO 340
      WRITE (ecall,FMT='(A)') string
  340 qgtstr = .TRUE.
C     - === End Echo Code   ===
      GO TO 370
C
C     if a null line is not acceptable
C
  350 WRITE (icuout,FMT=*)
     +  ' A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN'
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 360
      WRITE (ecall,FMT=*)
     +  ' A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN'
  360 REWIND (icuin)
C     - === End Echo Code   ===
  370 GO TO 60

  380 GO TO 10

  390 RETURN
C
C     at this point, we know that IOS = 0 (no error on read)
C          and we got a line
C
      STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO TERMINATE-PROCESSING
C
C     WHOOPS, LOOKS LIKE WE'VE BOMBED OUT!
C
  400 itries = itries - 1
      qgtstr = .FALSE.
      WRITE (icuout,FMT=9020) itries
C     + === Begin Echo Code ===
      IF (.NOT. (qecall)) GO TO 410
      WRITE (ecall,FMT=9020) itries
  410 CONTINUE

 9020 FORMAT (' READ WAS UNSUCCESSFUL AFTER ',I2,' TRIES.')
C     - === End Echo Code   ===
      GO TO i99996
C
      END
      LOGICAL FUNCTION qgtyn(mssg,qyn)
C**********************************************************************
C
C     LOGICAL FUNCTION QGTYN(FORMSG,QYN)
C
C                     Q GeT Y-es or N-o
C
C
C                              Function
C
C
C         If FORMSG is not blank, writes it to the current output unit,
C     then reads a line from the current input unit and examines the
C     first non-blank character for being 'Y', 'y', 'N', or 'n'. If
C     the first non-blank character is 'Y' or 'y', QYN is set to
C     .TRUE., if the first non-blank character is 'N' or 'n', QYN
C     is set to .FALSE.
C         If there are three unsuccessful attempts to obtain a
C     yes/no response, QGTYN returns .FALSE., otherwise, QGTYN
C     returns .TRUE.
C
C
C                              Note
C
C
C     Two error messages can be issued to current output as a result of
C     the use of this routine.
C
C          A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN
C     is issued by QGTSTR, which is called indirectly by this routine.
C
C          THE LEGAL FIRST NON-BLANK CHARACTER MUST BE ONE OF
C          THE FOLLOWING:
C          (followed by a list of legal characters)
C     is issued by QGTCHR, which is called by this routine.
C
C     FORMSG is rewritten after the above error messages.
C
C
C                              Arguments
C
C
C     QGTYN <-- .TRUE. if the first non-blank character obtained is a
C               member of 'YyNn' (even if it takes three attempts).
C               .FALSE. if an EOF is encountered on current input and
C               current input is not standard input or if three
C               attempts to obtain a character from 'YyNn' fail.
C                                             QGTYN is LOGICAL
C
C     FORMSG --> A FORTRAN format (first and last characters must be
C                left and right parentheses) which, when invoked,
C                will print a message to the user about the
C                yes or no answer that is to be specified.
C                If FORMSG is a character string consisting entirely
C                of blanks, then no message is printed to the user.
C                                      FORMSG is CHARACTER*(*)
C
C     QYN    <-- Defined only if QGTCHR is .TRUE.  If the response is
C                'Y' or 'y', this is .TRUE., otherwise .FALSE.
C                                               QYN is LOGICAL
C
C**********************************************************************
C
C
C
C     .. Scalar Arguments ..
      LOGICAL qyn
      CHARACTER mssg* (*)
C     ..
C     .. Local Scalars ..
      INTEGER iwhich
C     ..
C     .. External Functions ..
      LOGICAL qgtchr
      EXTERNAL qgtchr
C     ..
C     .. Executable Statements ..
C
      qgtyn = qgtchr(mssg,'YN',.TRUE.,iwhich)
C
      IF (.NOT. (qgtyn)) GO TO 30
      IF (.NOT. (iwhich.EQ.1)) GO TO 10
      qyn = .TRUE.
      GO TO 20

   10 qyn = .FALSE.
   20 CONTINUE
C
   30 END
      LOGICAL FUNCTION qyngt(mssg)
C
C**********************************************************************
C
C     LOGICAL FUNCTION QYNGT(FORMSG)
C
C                     Q Y-es or N-o GeT
C
C
C                              Function
C
C
C         If FORMSG is not blank, writes it to the current output unit,
C     then reads a line from the current input unit and examines the
C     first non-blank character for being 'Y', 'y', 'N', or 'n'. If
C     the first non-blank character is 'Y' or 'y', QYNGT is set to
C     .TRUE., if the first non-blank character is 'N' or 'n', QYN
C     is set to .FALSE.
C         If there are three unsuccessful attempts to obtain a
C     yes/no response, QYNGT STOPs.
C
C
C                              Note
C
C
C     Two error messages can be issued to current output as a result of
C     the use of this routine.
C
C          A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN
C     is issued by QGTSTR, which is called indirectly by this routine.
C
C          THE LEGAL FIRST NON-BLANK CHARACTER MUST BE ONE OF
C          THE FOLLOWING:
C          (followed by a list of legal characters)
C     is issued by QGTCHR, which is called by this routine.
C
C     FORMSG is rewritten after the above error messages.
C
C
C                              Arguments
C
C
C     QYNGT <-- .TRUE. if the first non-blank character obtained is a
C               member of 'Yy' (even if it takes three attempts).
C               .FALSE. if an the character is a member of 'Nn'.
C                                             QYNGT is LOGICAL
C
C     FORMSG --> A FORTRAN format (first and last characters must be
C                left and right parentheses) which, when invoked,
C                will print a message to the user about the
C                yes or no answer that is to be specified.
C                If FORMSG is a character string consisting entirely
C                of blanks, then no message is printed to the user.
C                                        FORMSG is CHARACTER*(*)
C
C
C
C                              Note
C
C
C     This routine merely renames QGTYN.
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER mssg* (*)
C     ..
C     .. Local Scalars ..
      LOGICAL qok
C     ..
C     .. External Functions ..
      LOGICAL qgtyn
      EXTERNAL qgtyn
C     ..
C     .. Executable Statements ..
      qok = qgtyn(mssg,qyngt)
      IF (.NOT. (qok)) STOP ' QGTYN CALLED FROM QYNGT'
      RETURN

      END
      REAL FUNCTION rlgt(mssg,rlow,rhi)
C
C**********************************************************************
C
C     REAL FUNCTION RLGT(FORMSG,RLOW,RHI)
C
C                         ReaL GeT
C
C
C                              Function
C
C
C          If FORMSG is non-blank, prints it to current output, then
C     reads a real value between RLOW and RHI (inclusive) from
C     the current input unit and returns the value.
C     RLGT STOPs if any problem is encountered, as described.
C          If there are three unsuccessful attempts to obtain the
C     number in the range specified, RLGT STOPS.
C     RLGT will also STOP if an EOF is encountered on
C     current input and current input is not standard input
C     (as determined by comparing unit numbers).
C
C
C                              Note
C
C
C     Three error messages can be issued to current output as a result
C     of the use of this routine.
C
C          A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN
C     is issued by this routine if a blank line or EOF is
C     encountered.
C
C          THE STRING ENTERED IS NOT A LEGAL REAL VALUE - TRY AGAIN
C     is issued for the obvious reason.
C
C          NUMBER MUST BE BETWEEN <lo> AND <hi> - TRY AGAIN
C     is issued if a legal but out of range value is encountered.
C
C          FORMSG is re-written to current output after either of the
C     previous error messages.
C
C
C                              Arguments
C
C
C     RLGT <-- The value of the real number obtained from current
C              input.
C                                                  RLGT is REAL
C
C     FORMSG --> A FORTRAN format (first and last characters must be
C                left and right parentheses) which, when invoked,
C                will print a message to the user about the
C                real number that is to be obtained.
C                If FORMSG is a character string consisting entirely of
C                blanks, then no message is printed to the user.
C                                       FORMSG is CHARACTER*(*)
C
C     RLOW --> The lowest legal value for the number obtained.
C                                                  RLOW is REAL
C
C     RHI --> The highest legal value for the number obtained.
C                                                  RHI is REAL
C
C
C
C                              Note
C
C
C     This routine merely renames QGTRL.
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      REAL rhi,rlow
      CHARACTER mssg* (*)
C     ..
C     .. Local Scalars ..
      LOGICAL qok
C     ..
C     .. External Functions ..
      LOGICAL qgtrl
      EXTERNAL qgtrl
C     ..
C     .. Executable Statements ..
      qok = qgtrl(mssg,rlow,rhi,rlgt)
      IF (.NOT. (qok)) STOP ' QGTRL CALLED FROM RLGT'
      RETURN

      END
      CHARACTER*(*) FUNCTION strgt(mssg,qnulok)
C
C**********************************************************************
C
C     CHARACTER*(*) FUNCTION STRGT(MSSG,QNULOK)
C
C                          STRing GeT
C
C
C
C                              Function
C
C
C          If MSSG is non-blank, prints it to current output, then
C     reads a line from the current input file and returns its
C     contents.
C         STRGT returns .FALSE. if (1) an EOF is encountered and
C     current input is not standard input; (2) an error was encountered
C     in the read; or (3) QNULOK is .FALSE. indicating that a null line
C     is not a legal input, and three successive null lines are obtained
C
C          If QNULOK is .FALSE. indicating that a null line is not a
C     legal input, the message
C          A NULL LINE IS NOT ACCEPTABLE HERE - TRY AGAIN
C     is written to current output and a new line is obtained.
C
C
C                              Note
C
C
C          If an EOF is encountered on standard input, then standard
C     input is rewound.  The EOF is treated as a null
C     line as is a line consisting only of blank characters.
C     If QNULOK is .TRUE. then a null line will cause STRING to
C     be returned filled with all blanks.
C
C
C                              Arguments
C
C
C     STRGT <-- The string read from current input.
C                                                STRGT is CHARACTER*(*)
C
C     MSSG --> A FORTRAN format (first and last characters must be
C              left and right parentheses) which, when invoked,
C              will print a message to the user about the
C              string that is to be obtained.
C              If MSSG is a character string consisting entirely of
C              blanks, then no message is printed to the user.
C                                              MSSG is CHARACTER*(*)
C
C     QNULOK --> .TRUE. if a null string is acceptable, else .FALSE.
C                                                QNULOK is LOGICAL
C
C
C                              Note
C
C
C     This routine merely renames QGTSTR.
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      LOGICAL qnulok
      CHARACTER mssg* (*)
C     ..
C     .. Local Scalars ..
      LOGICAL qok
C     ..
C     .. External Functions ..
      LOGICAL qgtstr
      EXTERNAL qgtstr
C     ..
C     .. Executable Statements ..
      qok = qgtstr(mssg,strgt,qnulok)
      IF (.NOT. (qok)) STOP ' QGTSTR CALLED FROM STRGT'
      RETURN

      END
      CHARACTER*1 FUNCTION trnchr(chr)
C
C**********************************************************************
C
C     CHARACTER*1 FUNCTION TRNCHR(CHR)
C                    TRaNslate to upper case one CHaRacter
C
C
C                              Function
C
C
C     If CHR is a lower case letter 'a..z' then returns the upper
C     case letter corresponding.
C
C     If CHR is not a lower case letter, then returns CHR.
C
C
C                              Arguments
C
C
C     TRNCHR <-- CHR or its upper case equivalent if CHR is in 'a..z'
C                                        TRNCHR is CHARACTER*1
C
C     CHR --> Character to be translated to upper case
C                                        CHR is CHARACTER*1
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER chr*1
C     ..
C     .. Local Scalars ..
      INTEGER ix
      LOGICAL qnotlc
      CHARACTER locase*26,upcase*26
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC index
C     ..
C     .. Data statements ..
      DATA locase/'abcdefghijklmnopqrstuvwxyz'/
      DATA upcase/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C     ..
C     .. Executable Statements ..
      qnotlc = chr .LT. 'a' .OR. chr .GT. 'z'
      IF (.NOT. (qnotlc)) GO TO 10
      trnchr = chr
      RETURN

   10 ix = index(locase,chr)
      IF (.NOT. (ix.LE.0)) GO TO 20
      trnchr = chr
      GO TO 30

   20 trnchr = upcase(ix:ix)
   30 RETURN

      END
