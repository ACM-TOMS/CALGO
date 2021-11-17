C*********************************************
C                                            *
C   PROGRAM   : PCOMP                        *
C   MODULE    : ER (ERROR MESSAGES)          *
C   ABSTRACT  : FORTRAN PRECOMPILER          *
C   KEY WORD  : AUTOMATIC DIFFERENTIATION    *
C   SOURCE    : PCOMP 2.3 BY M.LIEPELT       *
C               PCOMP 3.0 BY M.DOBMANN       *
C   COPYRIGHT : C.TRASSL, K.SCHITTKOWSKI     *
C               MATHEMATISCHES INSTITUT,     *
C               UNIVERSITAET BAYREUTH,       *
C               D-95440 BAYREUTH, GERMANY    *
C   DATE      : NOVEMBER 23, 1997            *
C   VERSION   : 5.3                          *
C                                            *
C*********************************************
C
C
C
      SUBROUTINE SYMERR (LNUM,N)
      INTEGER LNUM,N
C
C**********************************************************************
C                                                                      
C   S Y M E R R   -   INDICATE ERROR MESSAGE AND THE LINE NUMBER IN THE
C                     SOURCE CODE, THE ERROR OCCURED.           
C                                                                      
C   PARAMETERS:                                                       
C      LNUM   - LINE NUMBER OF ERROR.
C      N      - ERROR CODE.
C
C**********************************************************************
C
      CHARACTER*10 S
      INTEGER SLEN
C
      CALL ITOA2(LNUM,S,SLEN)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     1       21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     2       38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,
     3       55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71) N
    1 WRITE(*,1000) 'File not found, conpilation denied'
      GOTO 9999
    2 WRITE(*,1000) 'File too long, compilation denied'
      GOTO 9999
    3 WRITE(*,1010) 'Identifier expected',S(1:SLEN)
      GOTO 9999
    4 WRITE(*,1010) 'Identifier redefined',S(1:SLEN)
      GOTO 9999
    5 WRITE(*,1010) ''','' expected',S(1:SLEN)
      GOTO 9999
    6 WRITE(*,1010) '''('' expected',S(1:SLEN)
      GOTO 9999
    7 WRITE(*,1010) 'Identifier not declared',S(1:SLEN)
      GOTO 9999
    8 WRITE(*,1010) 'Type mismatch',S(1:SLEN)
      GOTO 9999
    9 WRITE(*,1010) 'Division by zero',S(1:SLEN)
      GOTO 9999
   10 WRITE(*,1010) 'Constant expected',S(1:SLEN)
      GOTO 9999
   11 WRITE(*,1010) 'Operator expected',S(1:SLEN)
      GOTO 9999
   12 WRITE(*,1010) 'Unexpected end of source file',S(1:SLEN)
      GOTO 9999
   13 WRITE(*,1010) '''.'' expected',S(1:SLEN)
      GOTO 9999
   14 WRITE(*,1010) ''')'' expected',S(1:SLEN)
      GOTO 9999
   15 WRITE(*,1010) 'THEN expected',S(1:SLEN)
      GOTO 9999
   16 WRITE(*,1010) 'ELSE expected',S(1:SLEN)
      GOTO 9999
   17 WRITE(*,1010) 'ENDIF expected',S(1:SLEN)
      GOTO 9999
   18 WRITE(*,1010) 'THEN without IF',S(1:SLEN)
      GOTO 9999
   19 WRITE(*,1010) 'ELSE without IF',S(1:SLEN)
      GOTO 9999
   20 WRITE(*,1010) 'ENDIF without IF',S(1:SLEN)
      GOTO 9999
   21 WRITE(*,1010) '''='' expected',S(1:SLEN)
      GOTO 9999
   22 WRITE(*,1010) 'Bad integer constant',S(1:SLEN)
      GOTO 9999
   23 WRITE(*,1010) 'Bad real constant',S(1:SLEN)
      GOTO 9999
   24 WRITE(*,1010) 'Formula too complex',S(1:SLEN)
      GOTO 9999
   25 WRITE(*,1010) 'Error in expression',S(1:SLEN)
      GOTO 9999
   26 WRITE(*,1010) 'Compiler error',S(1:SLEN)
      GOTO 9999
   27 WRITE(*,1010) 'Identifier not valid',S(1:SLEN)
      GOTO 9999
   28 WRITE(*,1010) 'Unknown type identifier',S(1:SLEN)
      GOTO 9999
   29 WRITE(*,1010) 'Unknown character',S(1:SLEN)
      GOTO 9999
   30 WRITE(*,1010) 'Yacc stack overflow',S(1:SLEN)
      GOTO 9999
   31 WRITE(*,1010) 'Syntax error',S(1:SLEN)
      GOTO 9999
   32 WRITE(*,1010) 'Out of memory',S(1:SLEN)
      GOTO 9999
   33 WRITE(*,1010) 'Bad index',S(1:SLEN)
      GOTO 9999
   34 WRITE(*,1000) 'Internal error of the dynamic parser'   
      GOTO 9999
   35 WRITE(*,1010) 'Wrong number of subscripts',S(1:SLEN)
      GOTO 9999
   36 WRITE(*,1010) 'Wrong number of arguments',S(1:SLEN)
      GOTO 9999
   37 WRITE(*,1010) 'Too many index sets',S(1:SLEN)
      GOTO 9999
   38 WRITE(*,1010) 'Too many integer constants',S(1:SLEN)
      GOTO 9999
   39 WRITE(*,1010) 'Too many real constants',S(1:SLEN)
      GOTO 9999
   40 WRITE(*,1010) 'Too many variables',S(1:SLEN)
      GOTO 9999
   41 WRITE(*,1010) 'Too many functions',S(1:SLEN)
      GOTO 9999
   42 WRITE(*,1010) 'Too many index variables',S(1:SLEN)
      GOTO 9999
   43 WRITE(*,1010) 'Number of variables not consistent',S(1:SLEN)
      GOTO 9999
   44 WRITE(*,1010) 'Number of functions not consistent',S(1:SLEN)
      GOTO 9999
   45 WRITE(*,1010) 'False end symbol',S(1:SLEN)
      GOTO 9999
   46 WRITE(*,1010) 'FORTRAN code exceeds line',S(1:SLEN)
      GOTO 9999
   47 WRITE(*,1010) '**: Domain error',S(1:SLEN)
      GOTO 9999
   48 WRITE(*,1010) 'Bad input format',S(1:SLEN)
      GOTO 9999
   49 WRITE(*,1010) 'Length of working array IWA too small',S(1:SLEN)
      GOTO 9999
   50 WRITE(*,1010) 'Length of working array WA too small',S(1:SLEN)
      GOTO 9999
   51 WRITE(*,1010) 'ATANH: Domain error',S(1:SLEN)
      GOTO 9999
   52 WRITE(*,1010) 'LOG: Domain error',S(1:SLEN)
      GOTO 9999
   53 WRITE(*,1010) 'SQRT: Domain error',S(1:SLEN)
      GOTO 9999
   54 WRITE(*,1010) 'ASIN: Domain error',S(1:SLEN)
      GOTO 9999
   55 WRITE(*,1010) 'ACOS: Domain error',S(1:SLEN)
      GOTO 9999
   56 WRITE(*,1010) 'ACOSH: Domain error',S(1:SLEN)
      GOTO 9999
   57 WRITE(*,1020) 'Label ',S(1:SLEN),' defined mulitply'
      GOTO 9999
   58 WRITE(*,1020) 'Label ',S(1:SLEN),' not found'
      GOTO 9999
   59 WRITE(*,1010) 'Wrong index expression',S(1:SLEN)
      GOTO 9999
   60 WRITE(*,1000) 'Wrong call of subroutine SYMINP'
      GOTO 9999
   61 WRITE(*,1000) 'Wrong call of subroutine SYMPRP'
      GOTO 9999
   62 WRITE(*,1000) 'Compileation of source file in GRAD-mode'
      GOTO 9999
   63 WRITE(*,1010) 'Wrong order of interpolation values',S(1:SLEN)
      GOTO 9999
   64 WRITE(*,1000) 'Insufficient memory for interpolation in '//
     &  'subroutine REVCDE'
      GOTO 9999
   65 WRITE(*,1000) 'Length of working array IWA in subroutine '//
     &  'SYMFOR too small'
      GOTO 9999
   66 WRITE(*,1010) 'Insufficient interpolation values',S(1:SLEN)
      GOTO 9999
   67 WRITE(*,1000) 'Compilation of source file not in GRAD-mode'
      GOTO 9999
   68 WRITE(*,1010) 'Missing macro name',S(1:SLEN)
      GOTO 9999
   69 WRITE(*,1010) 'More than MAXMAC macros defined',S(1:SLEN)
      GOTO 9999
   70 WRITE(*,1010) 'More than MAXBUF lines of macro statements',
     &  S(1:SLEN)
      GOTO 9999
   71 WRITE(*,1010) 'More than MAXBUF statements in function',S(1:SLEN)
      GOTO 9999
 1000 FORMAT(/,1X,'*** Error in PCOMP: ',A,'.')
 1010 FORMAT(/,1X,'*** Error in PCOMP: ',A,' (line ',A,').')
 1020 FORMAT(/,1X,'*** Error in PCOMP: ',3A,'.')
 9999 RETURN
      END
C
C
C
      SUBROUTINE ITOA2 (N,S,SLEN)
      INTEGER N
      CHARACTER*10 S
      INTEGER SLEN
C
      INTEGER I
C
      WRITE(S,'(I10)',ERR=11) N
      DO 10 I=1,10 
        IF (S(I:I) .NE. ' ') THEN
          S=S(I:10)
          SLEN=10-I+1
          RETURN
        ENDIF
 10   CONTINUE
 11   CONTINUE
      SLEN=10
      RETURN
      END
C*********************************************
C                                            *
C   PROGRAM   : PCOMP                        *
C   MODULE    : EV (EVALUATION)              *
C   ABSTRACT  : FORTRAN PRECOMPILER          *
C   KEY WORD  : AUTOMATIC DIFFERENTIATION    *
C   SOURCE    : PCOMP 2.3 by M.LIEPELT       *
C               PCOMP 3.0 by M.DOBMANN       *
C   COPYRIGHT : C.TRASSL, K.SCHITTKOWSKI     *
C               MATHEMATISCHES INSTITUT,     *
C               UNIVERSITAET BAYREUTH,       *
C               D-95440 BAYREUTH, GERMANY    *
C   DATE      : NOVEMBER 23, 1997            *
C   VERSION   : 5.3                          *
C                                            *
C*********************************************
C
C
C
      SUBROUTINE EVAL (START,IVAL,FVAL,MODE,PVVA,PVVAHE,MPIIS,MPVIS,
     1                 MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN,MPVFN,
     1                 MPVHE,MPVPF,MPIV,IINDEX,VINDEX,IICONS,VICONS,
     2                 IRCONS,VRCONS,IVARI,VVARI,IFUNC,VFUNC,VGRAD,
     3                 VHESS,VPFX,VINDVA,GSTACK,HSTACK,DFVAR,MVAR,IERR)
      INTEGER START,IVAL
      DOUBLE PRECISION FVAL
      INTEGER MODE
      INTEGER PVVA,PVVAHE
      INTEGER MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN
      INTEGER MPVFN,MPVHE,MPVPF,MPIV
      INTEGER IINDEX(MPIIS,5),VINDEX(MPVIS),IICONS(MPIIC,4)
      INTEGER VICONS(MPVIC),IRCONS(MPIRC,4),IVARI(MPIVA,3)
      INTEGER IFUNC(MPIFN,7),VPFX(MPVPF),VINDVA(MPIV)
      DOUBLE PRECISION VRCONS(MPVRC),VVARI(MPVVA),VFUNC(MPVFN)
      DOUBLE PRECISION VGRAD(MPVFN,MPVVA),VHESS(MPVHE,PVVAHE*PVVAHE)
      DOUBLE PRECISION GSTACK(MPVVA,10),HSTACK(PVVAHE*PVVAHE,10)
      INTEGER MVAR
      INTEGER DFVAR(*) 
      INTEGER IERR
C
      INTEGER ADD,SUB,MULT,DIV,POWER,LEFT,RIGHT,COMMA,ASSIGN,NLINE
      INTEGER RANGE,RELOP,AND,OR,NOT,INUM,RNUM,ID,SUM,PROD,IN,IF,THEN
      INTEGER ELSE,ENDIF,STDRD,EXTERN,INTERP,PARAM,INDEX,REAL,INT,TABLE
      INTEGER CONINT,LININT,SPLINE,VAR,INFUNC,FUNC,END,GOTO,MARKE,CONTIN
      INTEGER UMINUS,INDVAR,ENDSUM,ENDPRD,BEQ,BRA,LABEL,VECTOR,ACTIVE
      PARAMETER (ADD=43,SUB=45,MULT=42,DIV=47,POWER=94,LEFT=40,RIGHT=41)
      PARAMETER (COMMA=44,ASSIGN=61,NLINE=10,RANGE=257,RELOP=258)
      PARAMETER (AND=259,OR=260,NOT=261,INUM=262,RNUM=263,ID=264)
      PARAMETER (SUM=265,PROD=266,IN=267,IF=268,THEN=269,ELSE=270)
      PARAMETER (ENDIF=271,STDRD=272,EXTERN=273,INTERP=274,PARAM=275)
      PARAMETER (INDEX=276,REAL=277,INT=278,TABLE=279,CONINT=280)
      PARAMETER (LININT=281,SPLINE=282,VAR=283,INFUNC=284,FUNC=285)
      PARAMETER (END=286,GOTO=287,MARKE=288,CONTIN=289,UMINUS=290)
      PARAMETER (INDVAR=291,ENDSUM=292,ENDPRD=293,BEQ=294,BRA=295)
      PARAMETER (LABEL=296,VECTOR=297,ACTIVE=298)
C
      INTEGER MAXEXT
      PARAMETER (MAXEXT=1)
C
      INTEGER EXTTYP(MAXEXT)
      INTEGER EXT
      INTEGER EXTPAR(2)
C
      INTEGER ISMDEP,FSMDEP,GSMDEP,HSMDEP,RSMDEP
      PARAMETER (ISMDEP=10,FSMDEP=10,GSMDEP=10,HSMDEP=10,RSMDEP=40)
      INTEGER ISTACK(ISMDEP),RSTACK(RSMDEP)
      DOUBLE PRECISION FSTACK(FSMDEP)
      INTEGER ITOS,FTOS,GTOS,HTOS,RTOS
C
      INTEGER PC,X,DIM,I,J,HELP
      DOUBLE PRECISION Z,Z1,Z2,STEIG
      DOUBLE PRECISION X0,B,C,D
C
      INTEGER NOGRAD,GRAD,HESS
      PARAMETER (NOGRAD=0,GRAD=1,HESS=2)
C
      DATA EXTTYP /0/
C
      PC=START
      ITOS=0
      FTOS=0
      GTOS=0
      RTOS=0
      HTOS=0
      IF ((MODE .EQ. HESS) .AND. (PVVAHE .NE. PVVA)) THEN
        IERR=59
        RETURN
      ENDIF 
      IF (MVAR .EQ. 0) MODE = NOGRAD
C
      STEIG=0.0D0
      X0=0.0D0
      B=0.0D0
      C=0.0D0
      D=0.0D0
C
 1    X=VPFX(PC)
 10   IF (X .EQ. -1) THEN
        IF (ITOS .GT. 0) IVAL=ISTACK(ITOS)
        IF (FTOS .GT. 0) FVAL=FSTACK(FTOS)
        RETURN
      ENDIF
 20   IF (X .EQ. ADD+128) THEN
        ISTACK(ITOS-1)=ISTACK(ITOS-1)+ISTACK(ITOS)
        ITOS=ITOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 30   IF (X .EQ. SUB+128) THEN
        ISTACK(ITOS-1)=ISTACK(ITOS-1)-ISTACK(ITOS)
        ITOS=ITOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 40   IF (X .EQ. MULT+128) THEN
        ISTACK(ITOS-1)=ISTACK(ITOS-1)*ISTACK(ITOS)
        ITOS=ITOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 50   IF (X .EQ. DIV+128) THEN
        IF (ISTACK(ITOS) .EQ. 0) THEN
          IERR=9
          RETURN
        ENDIF
        ISTACK(ITOS-1)=ISTACK(ITOS-1)/ISTACK(ITOS)
        ITOS=ITOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 60   IF (X .EQ. UMINUS+128) THEN
        ISTACK(ITOS)=-ISTACK(ITOS)
        PC=PC+1
        GO TO 1
      ENDIF
 70   IF (X .EQ. INUM+128) THEN
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=VICONS(VPFX(PC+1))
        PC=PC+2
        GO TO 1
      ENDIF
 80   IF (X .EQ. INDVAR+128) THEN
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=VINDVA(VPFX(PC+1))
        PC=PC+2
        GO TO 1
      ENDIF
 85   IF (X .EQ. FUNC+128) THEN
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=IDNINT(VFUNC(IFUNC(VPFX(PC+1),4)))
        PC=PC+2
        GO TO 1
      ENDIF
 90   IF (X .EQ. ADD) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            DO 92 I=1,MVAR
              DO 91 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)=
     1            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)+
     2            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 91           CONTINUE
 92         CONTINUE
            HTOS=HTOS-1
          ENDIF
          DO 93 I=1,MVAR
            GSTACK(DFVAR(I),GTOS-1)=GSTACK(DFVAR(I),GTOS-1)+
     1                              GSTACK(DFVAR(I),GTOS)
 93       CONTINUE
          GTOS=GTOS-1
        ENDIF
        FSTACK(FTOS-1)=FSTACK(FTOS-1)+FSTACK(FTOS)
        FTOS=FTOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 100  IF (X .EQ. SUB) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            DO 102 I=1,MVAR
              DO 102 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)=
     1            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)-
     2            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 101          CONTINUE
 102        CONTINUE
            HTOS=HTOS-1
          ENDIF
          DO 103 I=1,MVAR
            GSTACK(DFVAR(I),GTOS-1)=GSTACK(DFVAR(I),GTOS-1)-
     1                              GSTACK(DFVAR(I),GTOS)
 103      CONTINUE
          GTOS=GTOS-1
        ENDIF
        FSTACK(FTOS-1)=FSTACK(FTOS-1)-FSTACK(FTOS)
        FTOS=FTOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 110  IF (X .EQ. MULT) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            DO 112 I=1,MVAR
              DO 111 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)=
     1           HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)*FSTACK(FTOS)+
     2           HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)*FSTACK(FTOS-1)+
     3           GSTACK(DFVAR(I),GTOS-1)*GSTACK(DFVAR(J),GTOS)+
     4           GSTACK(DFVAR(J),GTOS-1)*GSTACK(DFVAR(I),GTOS)
 111          CONTINUE
 112        CONTINUE
            HTOS=HTOS-1
          ENDIF
          DO 113 I=1,MVAR
            GSTACK(DFVAR(I),GTOS-1)=FSTACK(FTOS)*GSTACK(DFVAR(I),GTOS-1)
     1                             +FSTACK(FTOS-1)*GSTACK(DFVAR(I),GTOS)
 113      CONTINUE
          GTOS=GTOS-1
        ENDIF
        FSTACK(FTOS-1)=FSTACK(FTOS-1)*FSTACK(FTOS)
        FTOS=FTOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 120  IF (X .EQ. DIV) THEN
        IF (FSTACK(FTOS) .EQ. 0.0D0) THEN
          IERR=9
          RETURN
        ENDIF
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            Z1=FSTACK(FTOS)**2
            Z2=FSTACK(FTOS)**3
            DO 122 I=1,MVAR
              DO 121 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)=
     1            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)/
     2            FSTACK(FTOS)-GSTACK(DFVAR(I),GTOS-1)*
     3            GSTACK(DFVAR(J),GTOS)/Z1-
     4            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)*
     5            FSTACK(FTOS-1)/Z1-GSTACK(DFVAR(J),GTOS-1)*
     6            GSTACK(DFVAR(I),GTOS)/Z1+2*FSTACK(FTOS-1)*
     7            GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)/Z2
 121          CONTINUE
 122        CONTINUE
            HTOS=HTOS-1
          ENDIF
          Z=FSTACK(FTOS)**2
          DO 123 I=1,MVAR
            GSTACK(DFVAR(I),GTOS-1)=GSTACK(DFVAR(I),GTOS-1)/FSTACK(FTOS)
     1                           -FSTACK(FTOS-1)/Z*GSTACK(DFVAR(I),GTOS)
 123      CONTINUE
          GTOS=GTOS-1
        ENDIF
        FSTACK(FTOS-1)=FSTACK(FTOS-1)/FSTACK(FTOS)
        FTOS=FTOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 130  IF (X .EQ. POWER) THEN
        IF (FSTACK(FTOS-1) .LE. 0.0D0) THEN
          IERR=47
          RETURN
        ENDIF
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            Z=DLOG(FSTACK(FTOS-1))*FSTACK(FTOS-1)**FSTACK(FTOS)
            Z1=FSTACK(FTOS-1)**(FSTACK(FTOS)-1.0D0)
            Z2=FSTACK(FTOS-1)**(FSTACK(FTOS)-2.0D0)*FSTACK(FTOS)*
     1         (FSTACK(FTOS)-1.0D0)
            DO 132 I=1,MVAR
              DO 131 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)=
     1           Z*DLOG(FSTACK(FTOS-1))*GSTACK(DFVAR(I),GTOS)*
     2           GSTACK(DFVAR(J),GTOS)+Z1*DLOG(FSTACK(FTOS-1))*
     3           FSTACK(FTOS)*GSTACK(DFVAR(J),GTOS-1)*
     4           GSTACK(DFVAR(I),GTOS)+Z1*GSTACK(DFVAR(I),GTOS)*
     5           GSTACK(DFVAR(J),GTOS-1)+Z*
     6           HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     7           Z1*GSTACK(DFVAR(I),GTOS-1)*GSTACK(DFVAR(J),GTOS)+Z1*
     8           FSTACK(FTOS)*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)+
     9           Z1*FSTACK(FTOS)*DLOG(FSTACK(FTOS-1))*
     O           GSTACK(DFVAR(I),GTOS-1)*GSTACK(DFVAR(J),GTOS)+
     O           Z2*GSTACK(DFVAR(I),GTOS-1)*GSTACK(DFVAR(J),GTOS-1)
 131          CONTINUE
 132        CONTINUE
            HTOS=HTOS-1
          ENDIF
          Z1=FSTACK(FTOS)*FSTACK(FTOS-1)**(FSTACK(FTOS)-1.0D0)
          Z2=FSTACK(FTOS-1)**FSTACK(FTOS)*DLOG(FSTACK(FTOS-1))
          DO 133 I=1,MVAR
            GSTACK(DFVAR(I),GTOS-1)=Z1*GSTACK(DFVAR(I),GTOS-1)+
     1                              Z2*GSTACK(DFVAR(I),GTOS)
 133      CONTINUE
          GTOS=GTOS-1
        ENDIF
        FSTACK(FTOS-1)=FSTACK(FTOS-1)**FSTACK(FTOS)
        FTOS=FTOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 135  IF (X .EQ. POWER+128) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            IF ((VPFX(PC+1) .LT. 2) .AND. 
     1          (FSTACK(FTOS) .EQ. 0.0D0)) THEN
              IERR=47
              RETURN
            ENDIF
            Z1=VPFX(PC+1)*FSTACK(FTOS)**(VPFX(PC+1)-1)
            IF (FSTACK(FTOS).NE.0.0D0) THEN
            Z2=VPFX(PC+1)*(VPFX(PC+1)-1)*FSTACK(FTOS)**(VPFX(PC+1)-2)
            ELSE
            Z2=VPFX(PC+1)*(VPFX(PC+1)-1)
            ENDIF
            DO 137 I=1,MVAR
              DO 136 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1            Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)+
     2            Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 136          CONTINUE
 137        CONTINUE
          ENDIF
          IF ((VPFX(PC+1) .LT. 1) .AND. (FSTACK(FTOS) .EQ. 0.0D0)) THEN
            IERR=47
            RETURN
          ENDIF
          Z=VPFX(PC+1)*FSTACK(FTOS)**(VPFX(PC+1)-1)
          DO 138 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 138      CONTINUE
        ENDIF
        IF ((VPFX(PC+1) .LT. 0) .AND. (FSTACK(FTOS) .EQ. 0.0D0)) THEN
          IERR=47
          RETURN
        ENDIF
        FSTACK(FTOS)=FSTACK(FTOS)**VPFX(PC+1)
        PC=PC+2
        GO TO 1
      ENDIF
 140  IF (X .EQ. UMINUS) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            DO 142 I=1,MVAR
              DO 141 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1            -HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 141          CONTINUE
 142        CONTINUE
          ENDIF
          DO 143 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=-GSTACK(DFVAR(I),GTOS)
 143      CONTINUE
        ENDIF
        FSTACK(FTOS)=-FSTACK(FTOS)
        PC=PC+1
        GO TO 1
      ENDIF
 150  IF (X .EQ. RNUM) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 152 I=1,MVAR
              DO 151 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 151          CONTINUE
 152        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 153 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 153      CONTINUE
        ENDIF
        FTOS=FTOS+1
        IF (FTOS .GT. FSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        FSTACK(FTOS)=VRCONS(VPFX(PC+1))
        PC=PC+2
        GO TO 1
      ENDIF
 160  IF (X .EQ. INUM) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 162 I=1,MVAR
              DO 161 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 161          CONTINUE
 162        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 163 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 163      CONTINUE
        ENDIF
        FTOS=FTOS+1
        IF (FTOS .GT. FSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        FSTACK(FTOS)=DBLE(VICONS(VPFX(PC+1)))
        PC=PC+2
        GO TO 1
      ENDIF
 170  IF (X .EQ. INDVAR) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 172 I=1,MVAR
              DO 171 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 171          CONTINUE
 172        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 173 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 173      CONTINUE
        ENDIF
        FTOS=FTOS+1
        IF (FTOS .GT. FSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        FSTACK(FTOS)=DBLE(VINDVA(VPFX(PC+1)))
        PC=PC+2
        GO TO 1
      ENDIF
 180  IF (X .EQ. REAL) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 182 I=1,MVAR
              DO 181 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 181          CONTINUE
 182        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 183 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 183      CONTINUE
        ENDIF
        DIM=IRCONS(VPFX(PC+1),1)
        IF (DIM .EQ. 0) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          FSTACK(FTOS)=VRCONS(IRCONS(VPFX(PC+1),4))
          PC=PC+2
          GO TO 1
        ELSE IF (DIM .EQ. 1) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          IF ((ISTACK(ITOS) .LT. 0) .OR. (ISTACK(ITOS) .GT.
     1         IRCONS(VPFX(PC+1),2))) THEN
            IERR=33
            RETURN
          ENDIF
          FSTACK(FTOS)=VRCONS(IRCONS(VPFX(PC+1),4)+ISTACK(ITOS)-1)
          ITOS=ITOS-1
          PC=PC+2
          GO TO 1
        ELSE IF (DIM .EQ. 2) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          IF ((ISTACK(ITOS-1) .LE. 0) .OR. (ISTACK(ITOS) .LT. 0) .OR.
     1        (ISTACK(ITOS-1) .GT. IRCONS(VPFX(PC+1),2)) .OR.
     2        (ISTACK(ITOS) .GT. IRCONS(VPFX(PC+1),3))) THEN
            IERR=33
            RETURN
          ENDIF
          FSTACK(FTOS)=VRCONS(IRCONS(VPFX(PC+1),4)+(ISTACK(ITOS-1)-1)*
     1                        IRCONS(VPFX(PC+1),3)+ISTACK(ITOS)-1)
          ITOS=ITOS-2
          PC=PC+2
          GO TO 1
        ENDIF
      ENDIF
 190  IF (X .EQ. INT) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 192 I=1,MVAR
              DO 191 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 191          CONTINUE
 192        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 193 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 193      CONTINUE
        ENDIF
        DIM=IICONS(VPFX(PC+1),1)
        IF (DIM .EQ. 0) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          FSTACK(FTOS)=DBLE(VICONS(IICONS(VPFX(PC+1),4)))
          PC=PC+2
          GO TO 1
        ELSE IF (DIM .EQ. 1) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          IF ((ISTACK(ITOS) .LT. 0) .OR. (ISTACK(ITOS) .GT.
     1        IICONS(VPFX(PC+1),2))) THEN
            IERR=33
            RETURN
          ENDIF
          FSTACK(FTOS)=DBLE(VICONS(IICONS(VPFX(PC+1),4)+ISTACK(ITOS)-1))
          ITOS=ITOS-1
          PC=PC+2
          GO TO 1
        ELSE IF (DIM .EQ. 2) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          IF ((ISTACK(ITOS-1) .LE. 0) .OR. (ISTACK(ITOS) .LT. 0) .OR.
     1        (ISTACK(ITOS-1) .GT. IICONS(VPFX(PC+1),2)) .OR.
     2        (ISTACK(ITOS) .GT. IICONS(VPFX(PC+1),3))) THEN
            IERR=33
            RETURN
          ENDIF
          FSTACK(FTOS)=DBLE(VICONS(IICONS(VPFX(PC+1),4)+(ISTACK(ITOS-1)
     1                 -1)*IICONS(VPFX(PC+1),3)+ISTACK(ITOS)-1))
          ITOS=ITOS-2
          PC=PC+2
          GO TO 1
        ENDIF
      ENDIF
 200  IF (X .EQ. VAR) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 202 I=1,MVAR
              DO 201 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 201          CONTINUE
 202        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 203 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 203      CONTINUE
        ENDIF
        DIM=IVARI(VPFX(PC+1),1)
        IF (DIM .EQ. 0) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          FSTACK(FTOS)=VVARI(IVARI(VPFX(PC+1),3))
          IF (MODE .GE. GRAD) GSTACK(IVARI(VPFX(PC+1),3),GTOS)=1.0D0
          PC=PC+2
          GO TO 1
        ELSE IF (DIM .EQ. 1) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          IF ((ISTACK(ITOS) .LE. 0) .OR. (ISTACK(ITOS) .GT.
     1        IVARI(VPFX(PC+1),2))) THEN
            IERR=33
            RETURN
          ENDIF
          FSTACK(FTOS)=VVARI(IVARI(VPFX(PC+1),3)+ISTACK(ITOS)-1)
          IF (MODE .GE. GRAD) GSTACK(IVARI(VPFX(PC+1),3)
     1                        +ISTACK(ITOS)-1,GTOS)=1.0D0
          ITOS=ITOS-1
          PC=PC+2
          GO TO 1
        ENDIF
      ENDIF
 210  IF (X .EQ. FUNC) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
        ENDIF
        DIM=IFUNC(VPFX(PC+1),1)
        IF (DIM .LE. 0) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              DO 212 I=1,MVAR
                DO 211 J=1,MVAR
                 HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1             VHESS(IFUNC(VPFX(PC+1),7),(DFVAR(I)-1)*PVVA+DFVAR(J))
 211            CONTINUE
 212          CONTINUE 
            ENDIF
            DO 213 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=VGRAD(IFUNC(VPFX(PC+1),5),DFVAR(I))
 213        CONTINUE
          ENDIF
          FSTACK(FTOS)=VFUNC(IFUNC(VPFX(PC+1),4))
          PC=PC+2
          GO TO 1
        ELSE IF (DIM .EQ. 1) THEN
          FTOS=FTOS+1
          IF (FTOS .GT. FSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          IF ((ISTACK(ITOS) .LE. 0) .OR. (ISTACK(ITOS) .GT.
     1        IFUNC(VPFX(PC+1),3))) THEN
            IERR=33
            RETURN
          ENDIF
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              DO 215 I=1,MVAR
                DO 214 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              VHESS(IFUNC(VPFX(PC+1),7)+ISTACK(ITOS)-1,
     2                      (DFVAR(I)-1)*PVVA+DFVAR(J))
 214            CONTINUE
 215          CONTINUE
            ENDIF
            DO 216 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=VGRAD(IFUNC(VPFX(PC+1),5)+
     1                                    ISTACK(ITOS)-1,DFVAR(I))
 216        CONTINUE
          ENDIF
          FSTACK(FTOS)=VFUNC(IFUNC(VPFX(PC+1),4)+ISTACK(ITOS)-1)
          ITOS=ITOS-1
          PC=PC+2
          GO TO 1
        ENDIF
      ENDIF
 220  IF (X .EQ. STDRD) THEN
        IF (VPFX(PC+1) .EQ. 1) THEN
        Z=DSIGN(1.0D0,FSTACK(FTOS))
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              DO 2202 I=1,MVAR
                DO 2201 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 2201           CONTINUE
 2202         CONTINUE
            ENDIF
            DO 2203 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2203       CONTINUE
          ENDIF
          FSTACK(FTOS)=DABS(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 2) THEN
          IF (FSTACK(FTOS) .LT. 0.0D0) THEN
            IERR=53
            RETURN
          ENDIF
          IF (MODE .GE. GRAD) THEN
            IF (FSTACK(FTOS) .EQ. 0.0D0) THEN
              IERR=9
              RETURN
            ENDIF
            IF (MODE .EQ. HESS) THEN
              Z1=-0.25D0/DSQRT((FSTACK(FTOS)**3))
              Z2=0.5D0/DSQRT(FSTACK(FTOS))
              DO 2205 I=1,MVAR
                DO 2204 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z2*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z1*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2204           CONTINUE
 2205         CONTINUE
            ENDIF
            Z=0.5D0/DSQRT(FSTACK(FTOS))
            DO 2206 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2206       CONTINUE
          ENDIF
          FSTACK(FTOS)=DSQRT(FSTACK(FTOS))
          PC=PC+2
        GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 3) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
            Z=DEXP(FSTACK(FTOS))
              DO 2208 I=1,MVAR
                DO 2207 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2207           CONTINUE
 2208         CONTINUE
            ENDIF
            Z=DEXP(FSTACK(FTOS))
            DO 2209 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2209       CONTINUE
          ENDIF
          FSTACK(FTOS)=DEXP(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 4) THEN
          IF (FSTACK(FTOS) .LE. 0.0D0) THEN
            IERR=52
            RETURN
          ENDIF
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/FSTACK(FTOS)
              Z2=1.0D0/(FSTACK(FTOS)**2)
              DO 2211 I=1,MVAR
                DO 2210 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)-
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2210           CONTINUE
 2211         CONTINUE
            ENDIF
            Z=1.0D0/FSTACK(FTOS)
            DO 2212 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2212       CONTINUE
          ENDIF
          FSTACK(FTOS)=DLOG(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 5) THEN
          IF (FSTACK(FTOS) .LE. 0.0D0) THEN
            IERR=52
            RETURN
          ENDIF
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
            Z1=1.0D0/(FSTACK(FTOS)*DLOG(10.0D0))
            Z2=1.0D0/((FSTACK(FTOS)**2)*DLOG(10.0D0))
              DO 2214 I=1,MVAR
                DO 2213 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)-
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2213           CONTINUE
 2214         CONTINUE
            ENDIF
            Z=1.0D0/(FSTACK(FTOS)*DLOG(10.0D0))
            DO 2215 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2215       CONTINUE
          ENDIF
          FSTACK(FTOS)=DLOG10(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 6) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=DSIN(FSTACK(FTOS))
              Z2=DCOS(FSTACK(FTOS))
              DO 2217 I=1,MVAR
                DO 2216 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z2*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)-
     2              Z1*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2216           CONTINUE
 2217         CONTINUE
            ENDIF
            Z=DCOS(FSTACK(FTOS))
            DO 2218 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2218       CONTINUE
          ENDIF
          FSTACK(FTOS)=DSIN(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 7) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=-DSIN(FSTACK(FTOS))
              Z2=-DCOS(FSTACK(FTOS))
              DO 2220 I=1,MVAR
                DO 2219 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2219           CONTINUE
 2220         CONTINUE
            ENDIF
            Z=-DSIN(FSTACK(FTOS))
            DO 2221 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2221       CONTINUE
          ENDIF
          FSTACK(FTOS)=DCOS(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 8) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/DCOS(FSTACK(FTOS))**2
              Z2=2.0D0*DSIN(FSTACK(FTOS))/DCOS(FSTACK(FTOS))**3
              DO 2223 I=1,MVAR
                DO 2222 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2222           CONTINUE
 2223         CONTINUE
            ENDIF
            Z=1.0D0/DCOS(FSTACK(FTOS))**2
            DO 2224 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2224       CONTINUE
          ENDIF 
          FSTACK(FTOS)=DTAN(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 9) THEN
           IF (MODE .GE. GRAD) THEN
            IF (DABS(FSTACK(FTOS)) .GE. 1.0D0) THEN
              IERR = 54
              RETURN
            ENDIF
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/DSQRT(1.0D0-FSTACK(FTOS)**2)
              Z2=FSTACK(FTOS)/DSQRT((1.0D0-FSTACK(FTOS)**2)**3)
              DO 2226 I=1,MVAR
                DO 2225 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2225           CONTINUE
 2226         CONTINUE
            ENDIF
            Z=1.0D0/DSQRT(1.0D0-FSTACK(FTOS)**2)
            DO 2227 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2227       CONTINUE
          ENDIF
          IF (DABS(FSTACK(FTOS)) .GT. 1.0D0) THEN
            IERR = 54
            RETURN
          ENDIF
          FSTACK(FTOS)=DASIN(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 10) THEN
          IF (MODE .GE. GRAD) THEN
            IF (DABS(FSTACK(FTOS)) .GE. 1.0D0) THEN
              IERR = 55
              RETURN
            ENDIF
            IF (MODE .EQ. HESS) THEN
              Z1=-1.0D0/DSQRT(1.0D0-FSTACK(FTOS)**2)
              Z2=-FSTACK(FTOS)/DSQRT((1.0D0-FSTACK(FTOS)**2)**3)
              DO 2229 I=1,MVAR
                DO 2228 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2228           CONTINUE
 2229         CONTINUE
            ENDIF
            Z=-1.0D0/DSQRT(1.0D0-FSTACK(FTOS)**2)
            DO 2230 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2230       CONTINUE
          ENDIF
          IF (DABS(FSTACK(FTOS)) .GT. 1.0D0) THEN
            IERR = 55
            RETURN
          ENDIF
          FSTACK(FTOS)=DACOS(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 11) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/(1.0D0+FSTACK(FTOS)**2)
              Z2=-2.0D0*FSTACK(FTOS)/((1.0D0+FSTACK(FTOS)**2)**2)
              DO 2232 I=1,MVAR
                DO 2231 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2231           CONTINUE
 2232         CONTINUE
            ENDIF
            Z=1.0D0/(1.0D0+FSTACK(FTOS)**2)
            DO 2233 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2233       CONTINUE
          ENDIF
          FSTACK(FTOS)=DATAN(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 12) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=DCOSH(FSTACK(FTOS))
              Z2=DSINH(FSTACK(FTOS))
              DO 2235 I=1,MVAR
                DO 2234 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2234           CONTINUE
 2235         CONTINUE
            ENDIF
            Z=DCOSH(FSTACK(FTOS))
            DO 2236 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2236       CONTINUE
          ENDIF
          FSTACK(FTOS)=DSINH(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 13) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=DSINH(FSTACK(FTOS))
              Z2=DCOSH(FSTACK(FTOS))
              DO 2238 I=1,MVAR
                DO 2237 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2237           CONTINUE
 2238         CONTINUE
            ENDIF
            Z=DSINH(FSTACK(FTOS))
            DO 2239 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2239       CONTINUE
          ENDIF
          FSTACK(FTOS)=DCOSH(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 14) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/DCOSH(FSTACK(FTOS))**2
              Z2=-2.0D0*DSINH(FSTACK(FTOS))/DCOSH(FSTACK(FTOS))**3
              DO 2241 I=1,MVAR
                DO 2240 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2240           CONTINUE
 2241         CONTINUE
            ENDIF
            Z=1.0D0/DCOSH(FSTACK(FTOS))**2
            DO 2242 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2242       CONTINUE
          ENDIF
          FSTACK(FTOS)=DTANH(FSTACK(FTOS))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 15) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/DSQRT(1.0D0+FSTACK(FTOS)**2)
              Z2=-FSTACK(FTOS)/DSQRT((1.0D0+FSTACK(FTOS)**2)**3)
              DO 2244 I=1,MVAR
                DO 2243 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2243           CONTINUE
 2244         CONTINUE
            ENDIF
            Z=1.0D0/DSQRT(1.0D0+FSTACK(FTOS)**2)
            DO 2245 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2245       CONTINUE
          ENDIF
          FSTACK(FTOS)=DLOG(FSTACK(FTOS)+DSQRT(1.0D0+FSTACK(FTOS)**2))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 16) THEN
          IF (MODE .GE. GRAD) THEN
            IF (FSTACK(FTOS) .LE. 1.0D0) THEN
              IERR = 56
              RETURN
            ENDIF
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/DSQRT(FSTACK(FTOS)**2-1.0D0)
              Z2=-FSTACK(FTOS)/DSQRT((FSTACK(FTOS)**2-1.0D0)**3)
              DO 2247 I=1,MVAR
                DO 2246 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2246           CONTINUE
 2247         CONTINUE
            ENDIF
            Z=1.0D0/DSQRT(FSTACK(FTOS)**2-1.0D0)
            DO 2248 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2248       CONTINUE
          ENDIF
          IF (FSTACK(FTOS) .LT. 1.0D0) THEN
            IERR = 56
            RETURN
          ENDIF
          FSTACK(FTOS)=DLOG(FSTACK(FTOS)+DSQRT(FSTACK(FTOS)**2-1.0D0))
          PC=PC+2
          GOTO 1
C
        ELSE IF (VPFX(PC+1) .EQ. 17) THEN
          IF (DABS(FSTACK(FTOS)) .GE. 1.0D0) THEN
            IERR = 51
            RETURN
          ENDIF
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              Z1=1.0D0/(1.0D0-FSTACK(FTOS)**2)
              Z2=2.0D0*FSTACK(FTOS)/((1.0D0-FSTACK(FTOS)**2)**2)
              DO 2250 I=1,MVAR
                DO 2249 J=1,MVAR
                  HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1              Z1*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2              Z2*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)
 2249           CONTINUE
 2250         CONTINUE
            ENDIF
            Z=1.0D0/(1.0D0-FSTACK(FTOS)**2)
            DO 2251 I=1,MVAR
              GSTACK(DFVAR(I),GTOS)=Z*GSTACK(DFVAR(I),GTOS)
 2251       CONTINUE
          ENDIF
          FSTACK(FTOS)=0.5D0*DLOG((1.0D0+FSTACK(FTOS))/
     1                            (1.0D0-FSTACK(FTOS)))
          PC=PC+2
          GOTO 1
        ENDIF
      ENDIF
 230  IF (X .EQ. EXTERN) THEN
        EXT=VPFX(PC+1)
        IF (EXTTYP(EXT) .EQ. 1) THEN
          EXTPAR(1)=ISTACK(ITOS)
          ITOS=ITOS-1
        ELSE IF (EXTTYP(EXT) .EQ. 2) THEN
          EXTPAR(1)=ISTACK(ITOS-1)
          EXTPAR(2)=ISTACK(ITOS)
          ITOS=ITOS-2
        ENDIF
        FTOS=FTOS+1
        IF (FTOS .GT. FSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        CALL EXTFUN(EXT,VVARI,PVVA,FSTACK(FTOS),EXTPAR)
        IF (MODE .GE. GRAD) THEN
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          CALL EXTGRA(EXT,VVARI,PVVA,GSTACK(1,GTOS),EXTPAR)
        ENDIF
        IF (MODE .EQ. HESS) THEN
          HTOS=HTOS+1
          IF (HTOS .GT. HSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          CALL EXTHES(EXT,VVARI,PVVA,HSTACK(1,HTOS),EXTPAR)
        ENDIF
        PC=PC+2
        GO TO 1
      ENDIF
 240  IF (X .EQ. SUM) THEN
        FTOS=FTOS+1
        IF (FTOS .GT. FSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        RTOS=RTOS+4
        IF (RTOS .GT. RSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 242 I=1,MVAR
              DO 241 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 241          CONTINUE
 242        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 243 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 243      CONTINUE
        ENDIF
        FSTACK(FTOS)=0.0D0
        RSTACK(RTOS)=IINDEX(VPFX(PC+2),2)
        RSTACK(RTOS-1)=IINDEX(VPFX(PC+2),5)
        RSTACK(RTOS-2)=VPFX(PC+1)
        RSTACK(RTOS-3)=PC+3
        VINDVA(VPFX(PC+1))=VINDEX(RSTACK(RTOS-1))
        PC=PC+3
        GO TO 1
      ENDIF
 250  IF (X .EQ. ENDSUM) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            DO 252 I=1,MVAR
              DO 251 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)=
     1            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)+
     2            HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 251          CONTINUE
 252        CONTINUE
            HTOS=HTOS-1
          ENDIF
          DO 253 I=1,MVAR
            GSTACK(DFVAR(I),GTOS-1)=GSTACK(DFVAR(I),GTOS-1)+
     1                              GSTACK(DFVAR(I),GTOS)
 253      CONTINUE
          GTOS=GTOS-1
        ENDIF
        FSTACK(FTOS-1)=FSTACK(FTOS-1)+FSTACK(FTOS)
        FTOS=FTOS-1
        RSTACK(RTOS)=RSTACK(RTOS)-1
        IF (RSTACK(RTOS) .GT. 0) THEN
          RSTACK(RTOS-1)=RSTACK(RTOS-1)+1
          VINDVA(RSTACK(RTOS-2))=VINDEX(RSTACK(RTOS-1))
          PC=RSTACK(RTOS-3)
        ELSE
          RTOS=RTOS-4
          PC=PC+1
        ENDIF
        GO TO 1
      ENDIF
 260  IF (X .EQ. PROD) THEN
        FTOS=FTOS+1
        IF (FTOS .GT. FSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        RTOS=RTOS+4
        IF (RTOS .GT. RSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            HTOS=HTOS+1
            IF (HTOS .GT. HSMDEP) THEN
              IERR=24
              RETURN
            ENDIF
            DO 262 I=1,MVAR
              DO 261 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 261          CONTINUE
 262        CONTINUE
          ENDIF
          GTOS=GTOS+1
          IF (GTOS .GT. GSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          DO 263 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 263      CONTINUE
        ENDIF
        FSTACK(FTOS)=1.0D0
        RSTACK(RTOS)=IINDEX(VPFX(PC+2),2)
        RSTACK(RTOS-1)=IINDEX(VPFX(PC+2),5)
        RSTACK(RTOS-2)=VPFX(PC+1)
        RSTACK(RTOS-3)=PC+3
        VINDVA(VPFX(PC+1))=VINDEX(RSTACK(RTOS-1))
        PC=PC+3
        GO TO 1
      ENDIF
 270  IF (X .EQ. ENDPRD) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            DO 272 I=1,MVAR
              DO 271 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)=
     1           HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS-1)*FSTACK(FTOS)+
     2           HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)*FSTACK(FTOS-1)+
     2           GSTACK(DFVAR(I),GTOS-1)*GSTACK(DFVAR(J),GTOS)+
     3           GSTACK(DFVAR(J),GTOS-1)*GSTACK(DFVAR(I),GTOS)
 271          CONTINUE
 272        CONTINUE
            HTOS=HTOS-1
          ENDIF
          DO 273 I=1,MVAR
            GSTACK(DFVAR(I),GTOS-1)=FSTACK(FTOS)*GSTACK(DFVAR(I),GTOS-1)
     1                             +FSTACK(FTOS-1)*GSTACK(DFVAR(I),GTOS)
 273      CONTINUE
          GTOS=GTOS-1
        ENDIF
        FSTACK(FTOS-1)=FSTACK(FTOS-1)*FSTACK(FTOS)
        FTOS=FTOS-1
        RSTACK(RTOS)=RSTACK(RTOS)-1
        IF (RSTACK(RTOS) .GT. 0) THEN
          RSTACK(RTOS-1)=RSTACK(RTOS-1)+1
          VINDVA(RSTACK(RTOS-2))=VINDEX(RSTACK(RTOS-1))
          PC=RSTACK(RTOS-3)
        ELSE
          RTOS=RTOS-4
          PC=PC+1
        ENDIF
        GO TO 1
      ENDIF
 275  IF ((X .EQ. IF+128) .OR. (X .EQ. ENDIF+128)) THEN
        PC=PC+1
        GO TO 1
      ENDIF
 280  IF (X .EQ. RELOP) THEN
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        IF (VPFX(PC+1) .EQ. 1) THEN
          IF (FSTACK(FTOS-1) .EQ. FSTACK(FTOS)) THEN
            ISTACK(ITOS)=1
          ELSE
            ISTACK(ITOS)=0
          ENDIF
        ELSE IF (VPFX(PC+1) .EQ. 2) THEN
          IF (FSTACK(FTOS-1) .NE. FSTACK(FTOS)) THEN
            ISTACK(ITOS)=1
          ELSE
            ISTACK(ITOS)=0
          ENDIF
        ELSE IF (VPFX(PC+1) .EQ. 3) THEN
          IF (FSTACK(FTOS-1) .LT. FSTACK(FTOS)) THEN
            ISTACK(ITOS)=1
          ELSE
            ISTACK(ITOS)=0
          ENDIF
        ELSE IF (VPFX(PC+1) .EQ. 4) THEN
          IF (FSTACK(FTOS-1) .LE. FSTACK(FTOS)) THEN
            ISTACK(ITOS)=1
          ELSE
            ISTACK(ITOS)=0
          ENDIF
        ELSE IF (VPFX(PC+1) .EQ. 5) THEN
          IF (FSTACK(FTOS-1) .GT. FSTACK(FTOS)) THEN
            ISTACK(ITOS)=1
          ELSE
            ISTACK(ITOS)=0
          ENDIF
        ELSE IF (VPFX(PC+1) .EQ. 6) THEN
          IF (FSTACK(FTOS-1) .GE. FSTACK(FTOS)) THEN
            ISTACK(ITOS)=1
          ELSE
            ISTACK(ITOS)=0
          ENDIF
        ENDIF
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) HTOS=HTOS-2
          GTOS=GTOS-2
        ENDIF
        FTOS=FTOS-2
        PC=PC+2
        GO TO 1
      ENDIF
 290  IF (X .EQ. AND) THEN
        IF (ISTACK(ITOS-1)+ISTACK(ITOS) .EQ. 2) THEN
          ISTACK(ITOS-1)=1
        ELSE
          ISTACK(ITOS-1)=0
        ENDIF
        ITOS=ITOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 300  IF (X .EQ. OR) THEN
        IF (ISTACK(ITOS-1)+ISTACK(ITOS) .GT. 0) THEN
          ISTACK(ITOS-1)=1
        ELSE
          ISTACK(ITOS-1)=0
        ENDIF
        ITOS=ITOS-1
        PC=PC+1
        GO TO 1
      ENDIF
 310  IF (X .EQ. NOT) THEN
        ISTACK(ITOS)=1-ISTACK(ITOS)
        PC=PC+1
        GO TO 1
      ENDIF
 320  IF (X .EQ. BEQ) THEN
        IF (ISTACK(ITOS) .EQ. 0) THEN
          ITOS=ITOS-1
          PC=VPFX(PC+1)
        ELSE
          ITOS=ITOS-1
          PC=PC+2
        ENDIF
        GO TO 1
      ENDIF
 330  IF (X .EQ. BRA) THEN
        PC=VPFX(PC+1)
        GO TO 1
      ENDIF
 340  IF (X .EQ. GOTO) THEN
        PC=VPFX(PC+2)
        GO TO 1
      ENDIF
 350  IF (X .EQ. CONTIN) THEN
        PC=PC+2
        GO TO 1
      ENDIF
 360  IF (X .EQ. ASSIGN) THEN
        IF (IFUNC(VPFX(PC+1),1) .LE. 0) THEN
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              DO 402 I=0,PVVA-1
                DO 401 J=1,PVVA
                  VHESS(IFUNC(VPFX(PC+1),7),I*PVVA+J)=0.0D0
 401            CONTINUE
 402          CONTINUE
              DO 362 I=1,MVAR
                DO 361 J=1,MVAR
                  VHESS(IFUNC(VPFX(PC+1),7),(DFVAR(I)-1)*PVVA+DFVAR(J))=
     1              HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 361            CONTINUE
 362          CONTINUE
              HTOS=HTOS-1
            ENDIF
            DO 403 I=1,PVVA
              VGRAD(IFUNC(VPFX(PC+1),5),I)=0.0D0
 403        CONTINUE
            DO 363 I=1,MVAR
              VGRAD(IFUNC(VPFX(PC+1),5),DFVAR(I))=GSTACK(DFVAR(I),GTOS)
 363        CONTINUE
            GTOS=GTOS-1
          ENDIF
          VFUNC(IFUNC(VPFX(PC+1),4))=FSTACK(FTOS)
        ELSE
          IF (MODE .GE. GRAD) THEN
            IF (MODE .EQ. HESS) THEN
              DO 405 I=0,PVVA-1
                DO 404 J=1,PVVA
                  VHESS(IFUNC(VPFX(PC+1),7)+VINDVA(IFUNC(VPFX(PC+1),2))-
     1                        1,I*PVVA+J)=0.0D0
 404            CONTINUE
 405          CONTINUE
              DO 365 I=1,MVAR
                DO 364 J=1,MVAR
                  VHESS(IFUNC(VPFX(PC+1),7)+VINDVA(IFUNC(VPFX(PC+1),
     1              2))-1,(DFVAR(I)-1)*PVVA+DFVAR(J))=
     1              HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 364            CONTINUE
 365          CONTINUE
              HTOS=HTOS-1
            ENDIF
            DO 406 I=1,PVVA
               VGRAD(IFUNC(VPFX(PC+1),5)+VINDVA(IFUNC(VPFX(PC+1),2))-
     1               1,I)=0.0D0
 406        CONTINUE
            DO 366 I=1,MVAR
               VGRAD(IFUNC(VPFX(PC+1),5)+VINDVA(IFUNC(VPFX(PC+1),2))-
     1               1,DFVAR(I))=GSTACK(DFVAR(I),GTOS)
 366        CONTINUE
            GTOS=GTOS-1
          ENDIF
          VFUNC(IFUNC(VPFX(PC+1),4)+VINDVA(IFUNC(VPFX(PC+1),2))-1)=
     1      FSTACK(FTOS)
        ENDIF
        FTOS=FTOS-1
        PC=PC+2
        GO TO 1
      ENDIF
 370  IF ((X .EQ. -ASSIGN) .OR. (X .EQ. -FUNC)) THEN
        IERR=27
        RETURN
      ENDIF
 380  IF (X .EQ. CONINT) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
            DO 382 I=1,MVAR
              DO 381 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=0.0D0
 381          CONTINUE
 382        CONTINUE
          ENDIF
          DO 383 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=0.0D0
 383      CONTINUE
        ENDIF
        HELP=IRCONS(VPFX(PC+1),4)+IRCONS(VPFX(PC+1),2)
        IF (IRCONS(VPFX(PC+1),2) .EQ. 0) THEN
          FSTACK(FTOS)=0.0D0
        ELSE IF (FSTACK(FTOS) .LT. VRCONS(IRCONS(VPFX(PC+1),4))) THEN
          FSTACK(FTOS)=0.0D0
        ELSE IF (FSTACK(FTOS) .GE. VRCONS(HELP-2)) THEN
          FSTACK(FTOS)=VRCONS(HELP-1)
        ELSE
          DO 384 I=IRCONS(VPFX(PC+1),4),HELP-4,2
            IF ((FSTACK(FTOS) .GE. VRCONS(I)) .AND. 
     1          (FSTACK(FTOS) .LT. VRCONS(I+2)) ) THEN
              FSTACK(FTOS)=VRCONS(I+1)
              GO TO 385
            ENDIF
 384      CONTINUE
 385      CONTINUE
        ENDIF
        PC=PC+3
        GO TO 1
      ENDIF
 390  IF (X .EQ. LININT) THEN
        IF (MODE .GE. GRAD) THEN
          HELP=IRCONS(VPFX(PC+1),4)+IRCONS(VPFX(PC+1),2)
          IF (IRCONS(VPFX(PC+1),2) .EQ. 0) THEN
            STEIG=0.0D0
          ELSE IF (FSTACK(FTOS) .LT. VRCONS(IRCONS(VPFX(PC+1),4))) THEN
            STEIG=0.0D0
          ELSE IF (FSTACK(FTOS) .GE. VRCONS(HELP-2)) THEN
            STEIG=0.0D0
          ELSE
            DO 391 I=IRCONS(VPFX(PC+1),4),HELP-4,2
              IF ((FSTACK(FTOS) .GE. VRCONS(I)) .AND. 
     1            (FSTACK(FTOS) .LT. VRCONS(I+2)) ) THEN
                STEIG=(VRCONS(I+3)-VRCONS(I+1))/(VRCONS(I+2)-VRCONS(I))
                GO TO 392
              ENDIF
 391        CONTINUE
 392        CONTINUE
          ENDIF
          IF (MODE .EQ. HESS) THEN
            DO 394 I=1,MVAR
              DO 393 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1            STEIG*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 393          CONTINUE
 394        CONTINUE
          ENDIF
          DO 397 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=STEIG*GSTACK(DFVAR(I),GTOS)
 397      CONTINUE
        ENDIF
        HELP=IRCONS(VPFX(PC+1),4)+IRCONS(VPFX(PC+1),2)
        IF (IRCONS(VPFX(PC+1),2) .EQ. 0) THEN
          FSTACK(FTOS)=0.0D0
        ELSE IF (FSTACK(FTOS) .LT. VRCONS(IRCONS(VPFX(PC+1),4))) THEN
          FSTACK(FTOS)=VRCONS(IRCONS(VPFX(PC+1),4)+1)
        ELSE IF (FSTACK(FTOS) .GE. VRCONS(HELP-2)) THEN
          FSTACK(FTOS)=VRCONS(HELP-1)
        ELSE
          DO 398 I=IRCONS(VPFX(PC+1),4),HELP-4,2
            IF ((FSTACK(FTOS) .GE. VRCONS(I)) .AND. 
     1          (FSTACK(FTOS) .LT. VRCONS(I+2)) ) THEN
              STEIG=(VRCONS(I+3)-VRCONS(I+1))/(VRCONS(I+2)-VRCONS(I))
              FSTACK(FTOS)=VRCONS(I+1)+(FSTACK(FTOS)-VRCONS(I))*STEIG
              GO TO 399
            ENDIF
 398      CONTINUE
 399      CONTINUE
        ENDIF
        PC=PC+3
        GO TO 1
      ENDIF
 410  IF (X .EQ. SPLINE) THEN
        IF (MODE .GE. GRAD) THEN
          IF (MODE .EQ. HESS) THEN
C
          HELP=IRCONS(VPFX(PC+1),4)+IRCONS(VPFX(PC+1),2)
          IF (IRCONS(VPFX(PC+1),2) .EQ. 4*5) THEN
            B=VRCONS(IRCONS(VPFX(PC+1),4)+2)
            C=VRCONS(IRCONS(VPFX(PC+1),4)+3)
            D=VRCONS(IRCONS(VPFX(PC+1),4)+4)
            X0=0.0D0
          ELSE IF (FSTACK(FTOS) .LT. VRCONS(IRCONS(VPFX(PC+1),4)+15)) 
     1    THEN
            B=VRCONS(IRCONS(VPFX(PC+1),4)+2)
            C=VRCONS(IRCONS(VPFX(PC+1),4)+3)
            D=VRCONS(IRCONS(VPFX(PC+1),4)+4)
            X0=0.0D0
          ELSE IF (FSTACK(FTOS) .GE. VRCONS(HELP-5)) THEN
            B=VRCONS(HELP-8)
            C=VRCONS(HELP-7)
            D=VRCONS(HELP-6)
            X0=VRCONS(HELP-10)
          ELSE
            DO 411 I=IRCONS(VPFX(PC+1),4)+15,HELP-10,5
              IF ((FSTACK(FTOS) .GE. VRCONS(I)) .AND.
     1            (FSTACK(FTOS) .LT. VRCONS(I+5)) ) THEN
                B=VRCONS(I+2)
                C=VRCONS(I+3)
                D=VRCONS(I+4)
                X0=VRCONS(I)
                GO TO 412
              ENDIF
 411        CONTINUE
 412        CONTINUE
          ENDIF
C
            DO 414 I=1,MVAR
              DO 413 J=1,MVAR
                HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)=
     1          B*HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     2          2*C*GSTACK(DFVAR(I),GTOS)*GSTACK(DFVAR(J),GTOS)+
     3          2*C*(FSTACK(FTOS)-X0)*
     4          HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)+
     5          6*D*(FSTACK(FTOS)-X0)*GSTACK(DFVAR(I),GTOS)*
     6          GSTACK(DFVAR(J),GTOS)+
     7          3*D*(FSTACK(FTOS)-X0)**2 * 
     8          HSTACK((DFVAR(I)-1)*PVVA+DFVAR(J),HTOS)
 413          CONTINUE
 414        CONTINUE
          ENDIF
          HELP=IRCONS(VPFX(PC+1),4)+IRCONS(VPFX(PC+1),2)
          IF (IRCONS(VPFX(PC+1),2) .EQ. 4*5) THEN
            B=VRCONS(IRCONS(VPFX(PC+1),4)+2)
            C=VRCONS(IRCONS(VPFX(PC+1),4)+3)
            D=VRCONS(IRCONS(VPFX(PC+1),4)+4)
            X0=0.0D0
          ELSE IF (FSTACK(FTOS) .LT. VRCONS(IRCONS(VPFX(PC+1),4)+15)) 
     1    THEN
            B=VRCONS(IRCONS(VPFX(PC+1),4)+2)
            C=VRCONS(IRCONS(VPFX(PC+1),4)+3)
            D=VRCONS(IRCONS(VPFX(PC+1),4)+4)
            X0=0.0D0
          ELSE IF (FSTACK(FTOS) .GE. VRCONS(HELP-5)) THEN
            B=VRCONS(HELP-8)
            C=VRCONS(HELP-7)
            D=VRCONS(HELP-6)
            X0=VRCONS(HELP-10)
          ELSE
            DO 415 I=IRCONS(VPFX(PC+1),4)+15,HELP-10,5
              IF ((FSTACK(FTOS) .GE. VRCONS(I)) .AND.
     1            (FSTACK(FTOS) .LT. VRCONS(I+5)) ) THEN
                B=VRCONS(I+2)
                C=VRCONS(I+3)
                D=VRCONS(I+4)
                X0=VRCONS(I)
                GO TO 416
              ENDIF
 415        CONTINUE
 416        CONTINUE
          ENDIF
          DO 417 I=1,MVAR
            GSTACK(DFVAR(I),GTOS)=
     1        B*GSTACK(DFVAR(I),GTOS)+
     2        2*C*(FSTACK(FTOS)-X0) * GSTACK(DFVAR(I),GTOS)+
     3        3*D*(FSTACK(FTOS)-X0)**2 * GSTACK(DFVAR(I),GTOS)
 417      CONTINUE
        ENDIF
        HELP=IRCONS(VPFX(PC+1),2)+IRCONS(VPFX(PC+1),4)
        IF (IRCONS(VPFX(PC+1),2) .EQ. 4*5) THEN
          FSTACK(FTOS)=VRCONS(IRCONS(VPFX(PC+1),4)+1)+
     1                 VRCONS(IRCONS(VPFX(PC+1),4)+2)*FSTACK(FTOS)+
     2                 VRCONS(IRCONS(VPFX(PC+1),4)+3)*FSTACK(FTOS)**2+
     3                 VRCONS(IRCONS(VPFX(PC+1),4)+4)*FSTACK(FTOS)**3
        ELSE IF (FSTACK(FTOS) .LT. VRCONS(IRCONS(VPFX(PC+1),4)+15)) THEN
          FSTACK(FTOS)=VRCONS(IRCONS(VPFX(PC+1),4)+1)+
     1                 VRCONS(IRCONS(VPFX(PC+1),4)+2)*FSTACK(FTOS)+
     2                 VRCONS(IRCONS(VPFX(PC+1),4)+3)*FSTACK(FTOS)**2+
     3                 VRCONS(IRCONS(VPFX(PC+1),4)+4)*FSTACK(FTOS)**3
        ELSE IF (FSTACK(FTOS) .GE. VRCONS(HELP-5)) THEN
          FSTACK(FTOS)=VRCONS(HELP-9)+
     1                 VRCONS(HELP-8)*(FSTACK(FTOS)-VRCONS(HELP-10))+
     2                 VRCONS(HELP-7)*(FSTACK(FTOS)-VRCONS(HELP-10))**2+
     3                 VRCONS(HELP-6)*(FSTACK(FTOS)-VRCONS(HELP-10))**3
        ELSE
          DO 418 I=IRCONS(VPFX(PC+1),4)+15,HELP-10,5
            IF ((FSTACK(FTOS) .GE. VRCONS(I)) .AND.
     1          (FSTACK(FTOS) .LT. VRCONS(I+5)) ) THEN
              FSTACK(FTOS)=VRCONS(I+1)+
     1                     VRCONS(I+2)*(FSTACK(FTOS)-VRCONS(I))+
     2                     VRCONS(I+3)*(FSTACK(FTOS)-VRCONS(I))**2+
     3                     VRCONS(I+4)*(FSTACK(FTOS)-VRCONS(I))**3
              GO TO 419
            ENDIF
 418      CONTINUE
 419      CONTINUE
        ENDIF
C
        PC=PC+3
        GO TO 1
      ENDIF
C
      IERR=26
      WRITE(*,*) 'EVAL (3901) : unknown opcode ',X
      RETURN
      END
C
C
C
      SUBROUTINE SPLNES (VRCONS,DIM,IERR,LNUM)
      INTEGER DIM,IERR,LNUM
      DOUBLE PRECISION VRCONS((DIM+1)*5)
C
      DOUBLE PRECISION DY0,DYN
      DOUBLE PRECISION HELP,XK,YK,HELP1,HELP2,HELP3,HELP4,Y1,Y2,Y3,Y4
      INTEGER K
C     
      IF (DIM .LE. 2) THEN
        IERR=66
        LNUM=LNUM-1
        GO TO 9999
      ENDIF
      HELP1=VRCONS(2)
      HELP2=(VRCONS(7)-VRCONS(2))/(VRCONS(6)-VRCONS(1))
      HELP3=((VRCONS(12)-VRCONS(2))/
     1      ((VRCONS(11)-VRCONS(1))*(VRCONS(11)-VRCONS(6))))-
     2      (HELP2/(VRCONS(11)-VRCONS(6)))
      HELP4=((VRCONS(17)-VRCONS(2))/
     1      ((VRCONS(16)-VRCONS(1))*(VRCONS(16)-VRCONS(6))))-
     2      ((VRCONS(7)-VRCONS(2))/
     3      ((VRCONS(6)-VRCONS(1))*(VRCONS(16)-VRCONS(6))))-
     4      HELP3
      Y4 = HELP4/(VRCONS(16)-VRCONS(11))
      Y3 = HELP3 - (VRCONS(1) + VRCONS(6) + VRCONS(11)) * Y4
      Y2 = HELP2 - (VRCONS(1) + VRCONS(6)) * Y3 -
     1     (VRCONS(6)**2 + VRCONS(1)*VRCONS(6) + VRCONS(1)**2) * Y4
      Y1 = HELP1 - VRCONS(1) * Y2 - VRCONS(1)**2 * Y3 - 
     1             VRCONS(1)**3 * Y4
      DYN=0.0D0
      DY0=Y2 + 2.0D0*Y3*VRCONS(16) + 3.0D0*Y4*VRCONS(16)**2
      VRCONS(2) = Y1
      VRCONS(3) = Y2
      VRCONS(4) = Y3
      VRCONS(5) = Y4
C
      IF (DIM .GE. 4) THEN
        VRCONS(15+3) = 1.0D0
        HELP = VRCONS(15+6)-VRCONS(15+1)
        VRCONS(15+5) = (6.0D0/HELP)*((VRCONS(15+7)-VRCONS(15+2))/
     1                  HELP - DY0)
C
        DO 10 K=1+(4-1),DIM-1
          VRCONS(5*K+3) = (VRCONS(5*(K+1)+1)-VRCONS(5*K+1))/
     F                    (VRCONS(5*(K+1)+1)-VRCONS(5*(K-1)+1))
          VRCONS(5*K+4) = (VRCONS(5*K+1)-VRCONS(5*(K-1)+1))/
     F                    (VRCONS(5*(K+1)+1)-VRCONS(5*(K-1)+1))
          VRCONS(5*K+5) = 6.0D0/(VRCONS(5*(K+1)+1)-VRCONS(5*(K-1)+1))*
     F                    ( (VRCONS(5*(K+1)+2)-VRCONS(5*K+2))/
     F                      (VRCONS(5*(K+1)+1)-VRCONS(5*K+1)) -
     F                      (VRCONS(5*K+2)-VRCONS(5*(K-1)+2))/
     F                      (VRCONS(5*K+1)-VRCONS(5*(K-1)+1)) )
 10     CONTINUE
        VRCONS(5*DIM+4) = 1.0D0
        HELP = VRCONS(5*DIM+1)-VRCONS(5*(DIM-1)+1)
        VRCONS(5*DIM+5) = (6.0D0/HELP)*(DYN - 
     F                    (VRCONS(5*DIM+2)-VRCONS(5*(DIM-1)+2))/HELP)
C
        VRCONS(3+15) = -0.5D0 * VRCONS(3+15)
        VRCONS(4+15) = VRCONS(5+15)/2.0D0
        VRCONS(5*DIM+3) = 0.0D0
        DO 20 K=1+(4-1),DIM
          HELP = VRCONS(5*K+4)*VRCONS(5*(K-1)+3) + 2.0D0
          VRCONS(5*K+3) = -VRCONS(5*K+3)/HELP
          VRCONS(5*K+4) = (VRCONS(5*K+5)-VRCONS(5*K+4)*
     F                     VRCONS(5*(K-1)+4))/HELP
 20     CONTINUE
C
        VRCONS(5*DIM+5) = VRCONS(5*DIM+4)
        DO 30 K=DIM-1,0+(4-1), -1
          VRCONS(5*K+5) = VRCONS(5*K+3)*VRCONS(5*(K+1)+5)+VRCONS(5*K+4)
 30     CONTINUE
C
        DO 40 K=0+(4-1),DIM-1
          XK = VRCONS(5*(K+1)+1)-VRCONS(5*K+1)
          YK = VRCONS(5*(K+1)+2)-VRCONS(5*K+2)
          VRCONS(5*K+3) = YK/XK - (XK/6) *
     F                    (2*VRCONS(5*K+5)+VRCONS(5*(K+1)+5))
          VRCONS(5*K+4) = VRCONS(5*K+5)/2.0D0
          VRCONS(5*K+5) = (VRCONS(5*(K+1)+5)-VRCONS(5*K+5))/(6.0D0*XK)
 40     CONTINUE
      ENDIF
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      INTEGER FUNCTION GETIWA(FELD,DIM1,DIM2,IWA,LIWA,INFOLI)
C
      INTEGER FELD,DIM1,DIM2,LIWA
      INTEGER IWA(LIWA),INFOLI(15)
C
      INTEGER IIS,VIS,IIC,VIC,IRC,IVA,IFN,VPF,IV
      INTEGER LIIS,LIIC,LIRC,LIVA,LIFN
C
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,IVA=7,IFN=9,VPF=14,IV=15)
C
      LIIS=INFOLI(1)
C      LVIS=INFOLI(2)
      LIIC=INFOLI(3)
C      LVIC=INFOLI(4)
      LIRC=INFOLI(5)
      LIVA=INFOLI(7)
      LIFN=INFOLI(9)
C      LVPF=INFOLI(14)
C      LIV =INFOLI(15)
C
C
      IF (FELD .EQ. IIS) THEN
        IF (IWA(1) .GT. 1) THEN
          GETIWA=IWA(LIIS+IWA(1)*(DIM2-1)+DIM1)
        ELSE
          GETIWA=IWA(LIIS+DIM2)
        ENDIF
C
      ELSE IF (FELD .EQ. VIS) THEN
        GETIWA=IWA(INFOLI(2)+DIM1)
C 
      ELSE IF (FELD .EQ. IIC) THEN
        IF (IWA(3) .GT. 1) THEN
          GETIWA=IWA(LIIC+IWA(3)*(DIM2-1)+DIM1)
        ELSE
          GETIWA=IWA(LIIC+DIM2)
        ENDIF
C
      ELSE IF (FELD .EQ. VIC) THEN
        GETIWA=IWA(INFOLI(4)+DIM1)
C 
      ELSE IF (FELD .EQ. IRC) THEN
        IF (IWA(5) .GT. 1) THEN
          GETIWA=IWA(LIRC+IWA(5)*(DIM2-1)+DIM1)
        ELSE
          GETIWA=IWA(LIRC+DIM2)
        ENDIF
C
      ELSE IF (FELD .EQ. IVA) THEN
        IF (IWA(7) .GT. 1) THEN
          GETIWA=IWA(LIVA+IWA(7)*(DIM2-1)+DIM1)
        ELSE
          GETIWA=IWA(LIVA+DIM2)
        ENDIF
C
      ELSE IF (FELD .EQ. IFN) THEN
        IF (IWA(9) .GT. 1) THEN
          GETIWA=IWA(LIFN+IWA(9)*(DIM2-1)+DIM1)
        ELSE
          GETIWA=IWA(LIFN+DIM2)
        ENDIF
C
      ELSE IF (FELD .EQ. VPF) THEN
        GETIWA=IWA(INFOLI(14)+DIM1)
C 
      ELSE IF (FELD .EQ. IV) THEN
        GETIWA=IWA(INFOLI(15)+DIM1)
C 
      ELSE 
        STOP 'GETIWA: INTERNAL ERROR OF THE DYNAMIC PARSER.'
      ENDIF
      RETURN
      END

C*********************************************
C                                            *
C   PROGRAM   : PCOMP                        *
C   MODULE    : EX (EXTERNAL FUNCTIONS)      *
C   ABSTRACT  : FORTRAN PRECOMPILER          *
C   KEY WORD  : AUTOMATIC DIFFERENTIATION    *
C   SOURCE    : PCOMP 2.3 by M.LIEPELT       *
C               PCOMP 3.0 by M.DOBMANN       *
C   COPYRIGHT : C.TRASSL, K.SCHITTKOWSKI     *
C               MATHEMATISCHES INSTITUT,     *
C               UNIVERSITAET BAYREUTH,       *
C               D-95440 BAYREUTH, GERMANY    *
C   DATE      : NOVEMBER 23, 1997            *
C   VERSION   : 5.3                          *
C                                            *
C*********************************************
C
C
C
      SUBROUTINE EXTFUN (EXT,X,N,F,EXTPAR)
      INTEGER EXT
      INTEGER N
      DOUBLE PRECISION X(N),F
      INTEGER EXTPAR(2)
C
      IF (EXT .EQ. 1) THEN
        CALL EXF001(X,N,F)
        EXTPAR(1)=0
        EXTPAR(2)=0
        RETURN
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE EXTGRA (EXT,X,N,DF,EXTPAR)
      INTEGER EXT
      INTEGER N
      DOUBLE PRECISION X(N),DF(N)
      INTEGER EXTPAR(2)
C
      IF (EXT .EQ. 1) THEN
        CALL EXG001(X,N,DF)
        EXTPAR(1)=0
        EXTPAR(2)=0
        RETURN
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE EXTHES (EXT,X,N,DDF,EXTPAR)
      INTEGER EXT
      INTEGER N
      DOUBLE PRECISION X(N),DDF(N,N)
      INTEGER EXTPAR(2)
C
      IF (EXT .EQ. 1) THEN
        CALL EXH001(X,N,DDF)
        EXTPAR(1)=0
        EXTPAR(2)=0
        RETURN
      ENDIF
      RETURN
      END
C
C***************************************
C
C  EXAMPLE FORMAT FOR EXTERNAL FUNCTIONS
C
C***************************************
C
      SUBROUTINE EXF001 (X,N,F)
      INTEGER N
      DOUBLE PRECISION X(N),F
C
      F=0.0D0*X(N)
      RETURN
      END
C
C
C
      SUBROUTINE EXG001 (X,N,DF)
      INTEGER N
      DOUBLE PRECISION X(N),DF(N)
C
      INTEGER I
C
      DO 10 I=1,N
        DF(I)=0.0D0*X(I)
 10   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE EXH001 (X,N,DDF)
      INTEGER N
      DOUBLE PRECISION X(N),DDF(N,N)
C
      INTEGER I,J
C
      DO 10 I=1,N
        DO 20 J=1,N
          DDF(I,J)=0.0D0*X(I)
 20     CONTINUE
 10   CONTINUE
      RETURN
      END

C*********************************************
C                                            *
C   PROGRAM   : PCOMP                        *
C   MODULE    : G (FORTRAN CODE GENERATOR)   *
C   ABSTRACT  : FORTRAN PRECOMPILER          *
C   KEY WORD  : AUTOMATIC DIFFERENTIATION    *
C   SOURCE    : PCOMP 2.3 by M.LIEPELT       *
C               PCOMP 3.0 by M.DOBMANN       *
C   COPYRIGHT : C.TRASSL, K.SCHITTKOWSKI     *
C               MATHEMATISCHES INSTITUT,     *
C               UNIVERSITAET BAYREUTH,       *
C               D-95440 BAYREUTH, GERMANY    *
C   DATE      : NOVEMBER 23, 1997            *
C   VERSION   : 5.3                          *
C                                            *
C*********************************************
C                                                                      
C
C
      SUBROUTINE SYMPRP (SYMFIL,WA,LWA,IWA,LIWA,UWA,UIWA,IERR,MODE,
     /                   NVAR,NFUNC)
      INTEGER SYMFIL
      INTEGER LWA,LIWA
      DOUBLE PRECISION WA(LWA)
      INTEGER IWA(LIWA)
      INTEGER UWA,UIWA
      INTEGER IERR,MODE,NVAR,NFUNC
C
C**********************************************************************
C
C   S Y M P R P   -   LOAD INTERMEDIATE CODE GENERATED BY SYMINP FROM
C                     SYMFIL INTO WORKING ARRAYS.
C
C   PARAMETERS:
C      SYMFIL    - INPUT DEVICE; THE INTERMEDIATE CODE GENERATED BY
C                  SYMINP WAS WRITTEN TO THIS FILE AND IS NOW LOADED.
C      WA(LWA)   - REAL WORKING ARRAY, REQUIRED BY SYMPRP. ON RETURN,
C                  WA() CONTAINS THE INTERMEDIATE CODE.
C      IWA(LIWA) - INTEGER WORKING ARRAY, CF. WA().
C      UWA,UIWA  - INDICATE THE ACTUAL SPACE OF WA() AND IWA() THAT
C                  HAS BEEN USED BY THE SUBROUTINE.
C      IERR      - THE PARAMETER SHOWS THE REASON FOR TERMINATING THE
C                  SUBROUTINE. ON RETURN IERR COULD CONTAIN THE FOLLOW-
C                  ING VALUES:
C                  IERR = 0 : SUCCESSFUL TERMINATION.
C                  IERR > 0 : AN ERROR HAS BEEN DETECTED. FOR FURTHER
C                             INFORMATION CF. SUBROUTINE SYMERR.
C      MODE      - THE PARAMETER IS USED FOR THE RESERVATION OF SPACE
C                  FOR THE HESSIAN MATRIX
C      NVAR      - ON RETURN, NVAR CONTAINS THE NUMBER OF VARIABLES ON
C                  FUNCTION INPUT FILE
C      NFUNC     - ON RETURN, NFUNC CONTAINS THE NUMBER OF FUNCTIONS ON
C                  INPUT FILE
C
C**********************************************************************
C
      INTEGER I,PWA,PIWA,PX,PIX
C
      INTEGER GSMDEP,GETIWA,HILF
      PARAMETER (GSMDEP=10)
C
      INTEGER HSMDEP
      PARAMETER (HSMDEP=10)
C
      INTEGER GRAD,HESS
      PARAMETER (GRAD=1,HESS=2)
C
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
      INTEGER INFOLI(15)
C
      DO 10 I=1,15
        READ(SYMFIL,'(I6)',ERR=100) IWA(I)
 10   CONTINUE
      PIWA=IWA(1)*5+IWA(2)+IWA(3)*4+IWA(4)+IWA(5)*4+IWA(7)*3+
     1     IWA(9)*7+IWA(14)+15
      PWA=IWA(6)
      IF (MODE .EQ. GRAD) THEN
        PX=IWA(8)+IWA(11)+IWA(12)*IWA(8)+GSMDEP*IWA(8)+HSMDEP+1
      ELSE IF (MODE .EQ. HESS) THEN
        PX=IWA(8)+IWA(11)+IWA(12)*IWA(8)+IWA(13)*IWA(8)*IWA(8)+
     1     GSMDEP*IWA(8)+HSMDEP*IWA(8)*IWA(8)
      ELSE
        IERR=61
        RETURN
      ENDIF
      IF (MODE .NE. HESS) THEN
        IWA(13)=0
      ENDIF
      PIX=IWA(15)
      IF ((PWA+PX .GT. LWA) .OR. (PIWA+PIX .GT. LIWA)) THEN
        IERR=32
        RETURN
      ENDIF
      DO 20 I=16,PIWA
        READ(SYMFIL,'(I6)',ERR=100) IWA(I)
 20   CONTINUE
C
      NVAR=IWA(8)
      NFUNC=IWA(10)
      INFOLI(1)=15
      INFOLI(2)=IWA(1)*5+INFOLI(1)
      INFOLI(3)=IWA(2)+INFOLI(2)
      INFOLI(4)=IWA(3)*4+INFOLI(3)
      INFOLI(5)=IWA(4)+INFOLI(4)
      INFOLI(7)=IWA(5)*4+INFOLI(5)
      INFOLI(9)=IWA(7)*3+INFOLI(7)
      INFOLI(14)=IWA(9)*7+INFOLI(9)
      INFOLI(15)=IWA(14)+INFOLI(14)
C
      INFOLI(6)=0
      INFOLI(8)=IWA(6)+INFOLI(6)
      INFOLI(11)=IWA(8)+INFOLI(8)
      INFOLI(12)=IWA(11)+INFOLI(11)
      INFOLI(13)=IWA(12)+INFOLI(12)
C
      DO 50 I=1,IWA(9)
        HILF=GETIWA(IFN,I,1,IWA,LIWA,INFOLI)
        IF (HILF.EQ.1) THEN
          NFUNC=NFUNC + GETIWA(IFN,I,3,IWA,LIWA,INFOLI)-1
        ENDIF
   50 CONTINUE
C
      DO 30 I=1,PWA
        READ(SYMFIL,'(D24.17)',ERR=100) WA(I)
 30   CONTINUE
      IERR=0
      UWA=PWA+PX
      UIWA=PIWA+PIX
      RETURN
 100  IERR=26
      RETURN
      END
C
C
C
C     INTEGER FUNCTION GETIWA(FELD,DIM1,DIM2,IWA,LIWA,INFOLI)
C
C     INTEGER FELD,DIM1,DIM2,LIWA
C     INTEGER IWA(LIWA),INFOLI(15)
C
C     INTEGER IIS,VIS,IIC,VIC,IRC,IVA,IFN,VPF,IV
C     INTEGER LIIS,LVIS,LIIC,LVIC,LIRC,LIVA,LIFN,LVPF,LIV
C
C     PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,IVA=7,IFN=9,VPF=14,IV=15)
C
C     LIIS=INFOLI(1)
C     LVIS=INFOLI(2)
C     LIIC=INFOLI(3)
C     LVIC=INFOLI(4)
C     LIRC=INFOLI(5)
C     LIVA=INFOLI(7)
C     LIFN=INFOLI(9)
C     LVPF=INFOLI(14)
C     LIV =INFOLI(15)
C
C
C     IF (FELD .EQ. IIS) THEN
C       IF (IWA(1) .GT. 1) THEN
C         GETIWA=IWA(LIIS+IWA(1)*(DIM2-1)+DIM1)
C       ELSE
C         GETIWA=IWA(LIIS+DIM2)
C       ENDIF
C
C     ELSE IF (FELD .EQ. VIS) THEN
C       GETIWA=IWA(INFOLI(2)+DIM1)
C 
C     ELSE IF (FELD .EQ. IIC) THEN
C       IF (IWA(3) .GT. 1) THEN
C         GETIWA=IWA(LIIC+IWA(3)*(DIM2-1)+DIM1)
C       ELSE
C         GETIWA=IWA(LIIC+DIM2)
C       ENDIF
C
C     ELSE IF (FELD .EQ. VIC) THEN
C       GETIWA=IWA(INFOLI(4)+DIM1)
C 
C     ELSE IF (FELD .EQ. IRC) THEN
C       IF (IWA(5) .GT. 1) THEN
C         GETIWA=IWA(LIRC+IWA(5)*(DIM2-1)+DIM1)
C       ELSE
C         GETIWA=IWA(LIRC+DIM2)
C       ENDIF
C
C     ELSE IF (FELD .EQ. IVA) THEN
C       IF (IWA(7) .GT. 1) THEN
C         GETIWA=IWA(LIVA+IWA(7)*(DIM2-1)+DIM1)
C       ELSE
C         GETIWA=IWA(LIVA+DIM2)
C       ENDIF
C
C     ELSE IF (FELD .EQ. IFN) THEN
C       IF (IWA(9) .GT. 1) THEN
C         GETIWA=IWA(LIFN+IWA(9)*(DIM2-1)+DIM1)
C       ELSE
C         GETIWA=IWA(LIFN+DIM2)
C       ENDIF
C
C     ELSE IF (FELD .EQ. VPF) THEN
C       GETIWA=IWA(INFOLI(14)+DIM1)
C 
C     ELSE IF (FELD .EQ. IV) THEN
C       GETIWA=IWA(INFOLI(15)+DIM1)
C 
C     ENDIF
C     RETURN
C     END
C
C
C
      SUBROUTINE SYMFOR (XFIL,WA,LWA,IWA,LIWA,IERR)
      INTEGER XFIL
      INTEGER LWA,LIWA
      DOUBLE PRECISION WA(LWA)
      INTEGER IWA(LIWA)
      INTEGER IERR
C
C**********************************************************************
C
C   S Y M F O R   -   GENERATE EXECUTABLE FORTRAN CODE FOR EVALUATION
C                     OF FUNCTIONS AND GRADIENTS. 
C                                                                      
C   PARAMETERS:                                                         
C      XFIL      - OUTPUT DEVICE; THE EXECUTABLE FORTRAN CODE GENERATED
C                  BY SYMFOR IS WRITTEN TO THIS FILE.
C      WA(LWA)   - REAL WORKING ARRAY, CONTAINS THE INTERMEDIATE CODE
C                  GENERATED BY SYMINP.
C      IWA(LIWA) - INTEGER WORKING ARRAY, CF. WA(LWA). 
C      IERR      - THE PARAMETER SHOWS THE REASON FOR TERMINATING THE 
C                  SUBROUTINE. ON RETURN IERR COULD CONTAIN THE FOLLOW-
C                  ING VALUES:
C                  IERR = 0 : SUCCESSFUL TERMINATION.
C                  IERR > 0 : AN ERROR HAS BEEN DETECTED. FOR FURTHER
C                             INFORMATION CF. SUBROUTINE SYMERR.
C                                                                       
C**********************************************************************
C
C     INTEGER PIIS,PVIS,PIIC,PVIC,PIRC,PVRC,PIVA,PVVA,PIFN,PXFN,PVFN
      INTEGER PIIS,PVIS,PIIC,PVIC,PIRC,PIVA,PIFN
C     INTEGER PVGR,PVPF,PIV,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA,LVVA
      INTEGER PVPF,PIV,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA
C     INTEGER LIFN,LVFN,LVGR,LVPF,LIV
      INTEGER LIFN,LVPF,LIV
      INTEGER PVQD,LVQD
      INTEGER PIP,LIP,MAXVIP
C
      INTEGER MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN
      INTEGER MPXFN,MPVFN,MPVPF,MPIV
C
      INTEGER MAXVQD
C
      INTEGER MAXXC,MAXLC,MAXIC
      LOGICAL EXCALL
C
      PIIS=IWA(1)
      PVIS=IWA(2)
      PIIC=IWA(3)
      PVIC=IWA(4)
      PIRC=IWA(5)
C      PVRC=IWA(6)
      PIVA=IWA(7)
C      PVVA=IWA(8)
      PIFN=IWA(9)
C     PXFN=IWA(10)
C      PVFN=IWA(11)
C     PVGR=IWA(12)
      PVPF=IWA(14)
      PIV=IWA(15)
      PIP=0
      LIIS=16
      LVIS=LIIS+PIIS*5
      LIIC=LVIS+PVIS
      LVIC=LIIC+PIIC*4
      LIRC=LVIC+PVIC
      LVRC=1
      LIVA=LIRC+PIRC*4
C      LVVA=LVRC+PVRC
      LIFN=LIVA+PIVA*3
C      LVFN=LVVA+PVVA
C     LVGR=LVFN+PVFN
      LVPF=LIFN+PIFN*7
      LIV=LVPF+PVPF
      MAXVIP=50
      LIP=LIV+MAXVIP
      LVQD=LIP+PIV
      IF ((LIWA-5+1) .LT. LVQD) THEN
         IERR=65
         RETURN
      ENDIF
      MAXVQD=(LIWA-LVQD+1)/5
C
            MPIIS=MAX(1,IWA(1))
            MPVIS=MAX(1,IWA(2))
            MPIIC=MAX(1,IWA(3))
            MPVIC=MAX(1,IWA(4))
            MPIRC=MAX(1,IWA(5))
            MPVRC=MAX(1,IWA(6))
            MPIVA=MAX(1,IWA(7))
            MPVVA=MAX(1,IWA(8))
            MPIFN=MAX(1,IWA(9))
            MPXFN=MAX(1,IWA(10))
            MPVFN=MAX(1,IWA(11))
            MPVPF=MAX(1,IWA(14))
            MPIV =MAX(1,IWA(15))
      CALL REVCDE(MPIIS,MPIIC,MPVIC,MPIRC,MPIVA,MPVVA,MPIFN,MPXFN,MPVPF,
     1            IWA(LIIS),IWA(LIIC),IWA(LVIC),IWA(LIRC),IWA(LIVA),
     2            IWA(LIFN),IWA(LVPF),IWA(LIP),MAXVIP,PIP,MAXVQD,
     3            IWA(LVQD),PVQD,MAXXC,MAXLC,MAXIC,EXCALL,IERR)
C
c      CALL REVCDE(PIIS,PIIC,PVIC,PIRC,PIVA,PVVA,PIFN,PXFN,PVPF,
c     1            IWA(LIIS),IWA(LIIC),IWA(LVIC),IWA(LIRC),IWA(LIVA),
c     2            IWA(LIFN),IWA(LVPF),IWA(LIP),MAXVIP,PIP,MAXVQD,
c     3            IWA(LVQD),PVQD,MAXXC,MAXLC,MAXIC,EXCALL,IERR)
      IF (IERR .NE. 0) THEN
        RETURN
      ENDIF
C
      CALL FORCDE(XFIL,MPVIS,MPVIC,MPVRC,MPVVA,MPIFN,MPXFN,MPVFN,MPIV,
     1            IWA(LVIS),IWA(LVIC),WA(LVRC),MAXVQD,IWA(LVQD),MAXXC,
     2            MAXLC,MAXIC,EXCALL,IERR,IWA(LIP),MAXVIP,PIP)
C
c      CALL FORCDE(XFIL,PVIS,PVIC,PVRC,PVVA,PIFN,PXFN,PVFN,PIV,
c     1            IWA(LVIS),IWA(LVIC),WA(LVRC),MAXVQD,IWA(LVQD),MAXXC,
c     2            MAXLC,MAXIC,EXCALL,IERR,IWA(LIP),MAXVIP,PIP)
      RETURN
      END
C
C
C
      SUBROUTINE REVCDE (PIIS,PIIC,PVIC,PIRC,PIVA,PVVA,PIFN,PXFN,PVPF,
     1                   IINDEX,IICONS,VICONS,IRCONS,IVARI,IFUNC,VPFX,
     2                   VIP,MAXVIP,PIP,MAXVQD,VQD,PVQD,MAXXC,MAXLC,
     3                   MAXIC,EXCALL,IERR)
      INTEGER PIIS,PIIC,PVIC,PIRC,PIVA,PVVA,PIFN,PXFN,PVPF
      INTEGER IINDEX(PIIS,5),IICONS(PIIC,4),VICONS(PVIC)
      INTEGER IRCONS(PIRC,4),IVARI(PIVA,3),IFUNC(PIFN,7),VPFX(PVPF)
      INTEGER MAXVQD
      INTEGER VQD(MAXVQD,5)
      INTEGER PVQD
      INTEGER MAXXC,MAXLC,MAXIC
      LOGICAL EXCALL
      INTEGER IERR
      INTEGER PIP,MAXVIP
      INTEGER VIP(MAXVIP)
C
      INTEGER ADD,SUB,MULT,DIV,POWER,LEFT,RIGHT,COMMA,ASSIGN,NLINE
      INTEGER RANGE,RELOP,AND,OR,NOT,INUM,RNUM,ID,SUM,PROD,IN,IF,THEN
      INTEGER ELSE,ENDIF,STDRD,EXTERN,INTERP,PARAM,INDEX,REAL,INT,TABLE
      INTEGER CONINT,LININT,SPLINE,VAR,INFUNC,FUNC,END,GOTO,MARKE,CONTIN
      INTEGER UMINUS,INDVAR,ENDSUM,ENDPRD,BEQ,BRA,LABEL,VECTOR,ACTIVE
      PARAMETER (ADD=43,SUB=45,MULT=42,DIV=47,POWER=94,LEFT=40,RIGHT=41)
      PARAMETER (COMMA=44,ASSIGN=61,NLINE=10,RANGE=257,RELOP=258)
      PARAMETER (AND=259,OR=260,NOT=261,INUM=262,RNUM=263,ID=264)
      PARAMETER (SUM=265,PROD=266,IN=267,IF=268,THEN=269,ELSE=270)
      PARAMETER (ENDIF=271,STDRD=272,EXTERN=273,INTERP=274,PARAM=275)
      PARAMETER (INDEX=276,REAL=277,INT=278,TABLE=279,CONINT=280)
      PARAMETER (LININT=281,SPLINE=282,VAR=283,INFUNC=284,FUNC=285)
      PARAMETER (END=286,GOTO=287,MARKE=288,CONTIN=289,UMINUS=290)
      PARAMETER (INDVAR=291,ENDSUM=292,ENDPRD=293,BEQ=294,BRA=295)
      PARAMETER (LABEL=296,VECTOR=297,ACTIVE=298)
C
      INTEGER MAXSTD,MAXEXT
      PARAMETER (MAXSTD=17,MAXEXT=1)
      INTEGER STDTYP(MAXSTD),EXTTYP(MAXEXT)
C
      INTEGER LSMDEP,RSMDEP,XSMDEP,ISMDEP
      PARAMETER (LSMDEP=10,RSMDEP=40,XSMDEP=10,ISMDEP=10)
      INTEGER LSTACK(LSMDEP),RSTACK(RSMDEP)
      INTEGER XSTACK(XSMDEP),ISTACK(ISMDEP)
      INTEGER LTOS,RTOS,XTOS,ITOS
      INTEGER DIM,FC,K,LC,PC,QC,QC1,QC2,X,XC,IC
C
      DATA STDTYP /1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
      DATA EXTTYP /0/
C
      QC=2
      FC=1
      MAXXC=0
      MAXLC=0
      MAXIC=0
      EXCALL=.FALSE.
C
 1    LC=0
      XC=PVVA
      LTOS=0
      RTOS=0
      XTOS=0
      PC=IFUNC(FC,6)
      IC=0
      ITOS=0
C------------------------------------------------------------
      IF (IFUNC(FC,1) .EQ. 1) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        VQD(QC,1)=VECTOR
        VQD(QC,2)=IFUNC(FC,2)
        VQD(QC,3)=IFUNC(FC,3)
        VQD(QC,4)=0
        VQD(QC,5)=0
C*******************************************************
C That's really bad - change the parser YYPAR!
C*******************************************************
        DO 7 K=1,PIIS
          IF ((IINDEX(K,1) .EQ. 1) .AND.
     1        (IINDEX(K,3) .EQ. 1) .AND.
     2        (IINDEX(K,4) .EQ. IFUNC(FC,3))) THEN
            VQD(QC,4)=IINDEX(K,5)
          ENDIF
 7      CONTINUE
C*******************************************************
        RTOS=RTOS+1
        IF (RTOS .GT. RSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        RSTACK(RTOS)=QC
      ENDIF
C------------------------------------------------------------
      QC=QC+1
      IF (QC .GT. MAXVQD) THEN
        IERR=32
        RETURN
      ENDIF
C
      IF (IFUNC(FC,1) .LE. 0) THEN
        VQD(QC,1)=ACTIVE
        VQD(QC,2)=IFUNC(FC,4)
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=0
      ELSE
        VQD(QC,1)=ACTIVE
        VQD(QC,2)=-IFUNC(FC,2)
        VQD(QC,3)=IFUNC(FC,4)-1
        VQD(QC,4)=0
        VQD(QC,5)=0
      ENDIF
C------------------------------------------------------------
C
 2    X=VPFX(PC)
 10   IF (X .EQ. -1) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        VQD(QC,1)=-1
        VQD(QC,2)=0
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=0
C------------------------------------------------------------
        IF (IFUNC(FC,1) .EQ. 1) THEN
          VQD(QC,2)=QC
          VQD(RSTACK(RTOS),5)=QC
          RTOS=RTOS-1
        ENDIF
C------------------------------------------------------------
        FC=FC+1
        IF (FC .GT. PXFN) THEN
C       IF (FC .GT. PIFN) THEN
          PVQD=QC
          RETURN
        ELSE
          GO TO 1
        ENDIF
      ENDIF
C
 20   IF (X .EQ. ADD+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=ADD+128
        VQD(QC,2)=0
        VQD(QC,3)=ISTACK(ITOS-1)
        VQD(QC,4)=ISTACK(ITOS)
        VQD(QC,5)=IC
        ISTACK(ITOS-1)=IC
        ITOS=ITOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 30   IF (X .EQ. SUB+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=SUB+128
        VQD(QC,2)=0
        VQD(QC,3)=ISTACK(ITOS-1)
        VQD(QC,4)=ISTACK(ITOS)
        VQD(QC,5)=IC
        ISTACK(ITOS-1)=IC
        ITOS=ITOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 40   IF (X .EQ. MULT+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=MULT+128
        VQD(QC,2)=0
        VQD(QC,3)=ISTACK(ITOS-1)
        VQD(QC,4)=ISTACK(ITOS)
        VQD(QC,5)=IC
        ISTACK(ITOS-1)=IC
        ITOS=ITOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 50   IF (X .EQ. DIV+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=DIV+128
        VQD(QC,2)=0
        VQD(QC,3)=ISTACK(ITOS-1)
        VQD(QC,4)=ISTACK(ITOS)
        VQD(QC,5)=IC
        ISTACK(ITOS-1)=IC
        ITOS=ITOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 60   IF (X .EQ. UMINUS+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=UMINUS+128
        VQD(QC,2)=0
        VQD(QC,3)=ISTACK(ITOS)
        VQD(QC,4)=0
        VQD(QC,5)=IC
        ISTACK(ITOS)=IC
        PC=PC+1
        GO TO 2
      ENDIF
 70   IF (X .EQ. INUM+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=INUM+128
        VQD(QC,2)=VICONS(VPFX(PC+1))
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=IC
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=IC
        PC=PC+2
        GO TO 2
      ENDIF
 80   IF (X .EQ. FUNC+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=FUNC+128
        VQD(QC,2)=IFUNC(VPFX(PC+1),4)
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=IC
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=IC
        PC=PC+2
        GO TO 2
      ENDIF
 85   IF (X .EQ. INDVAR+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        IC=IC+1
        VQD(QC,1)=INDVAR+128
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=IC
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=IC
        PC=PC+2
        GO TO 2
      ENDIF
 90   IF (X .EQ. ADD) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=ADD
        VQD(QC,2)=0
        VQD(QC,3)=XSTACK(XTOS-1)
        VQD(QC,4)=XSTACK(XTOS)
        VQD(QC,5)=XC
        XSTACK(XTOS-1)=XC
        XTOS=XTOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 100  IF (X .EQ. SUB) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=SUB
        VQD(QC,2)=0
        VQD(QC,3)=XSTACK(XTOS-1)
        VQD(QC,4)=XSTACK(XTOS)
        VQD(QC,5)=XC
        XSTACK(XTOS-1)=XC
        XTOS=XTOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 110  IF (X .EQ. MULT) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=MULT
        VQD(QC,2)=0
        VQD(QC,3)=XSTACK(XTOS-1)
        VQD(QC,4)=XSTACK(XTOS)
        VQD(QC,5)=XC
        XSTACK(XTOS-1)=XC
        XTOS=XTOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 120  IF (X .EQ. DIV) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=DIV
        VQD(QC,2)=0
        VQD(QC,3)=XSTACK(XTOS-1)
        VQD(QC,4)=XSTACK(XTOS)
        VQD(QC,5)=XC
        XSTACK(XTOS-1)=XC
        XTOS=XTOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 130  IF (X .EQ. POWER) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=POWER
        VQD(QC,2)=0
        VQD(QC,3)=XSTACK(XTOS-1)
        VQD(QC,4)=XSTACK(XTOS)
        VQD(QC,5)=XC
        XSTACK(XTOS-1)=XC
        XTOS=XTOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 135  IF (X .EQ. POWER+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=POWER+128
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=XSTACK(XTOS)
        VQD(QC,4)=0
        VQD(QC,5)=XC
        XSTACK(XTOS)=XC
        PC=PC+2
        GO TO 2
      ENDIF
 140  IF (X .EQ. UMINUS) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=UMINUS
        VQD(QC,2)=0
        VQD(QC,3)=XSTACK(XTOS)
        VQD(QC,4)=0
        VQD(QC,5)=XC
        XSTACK(XTOS)=XC
        PC=PC+1
        GO TO 2
      ENDIF
 150  IF (X .EQ. RNUM) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=RNUM
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=XC
        XTOS=XTOS+1
        IF (XTOS .GT. XSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        XSTACK(XTOS)=XC
        PC=PC+2
        GO TO 2
      ENDIF
 160  IF (X .EQ. INUM) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=INUM
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=XC
        XTOS=XTOS+1
        IF (XTOS .GT. XSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        XSTACK(XTOS)=XC
        PC=PC+2
        GO TO 2
      ENDIF
 170  IF (X .EQ. INDVAR) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC,1)=INDVAR
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=XC
        XTOS=XTOS+1
        IF (XTOS .GT. XSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        XSTACK(XTOS)=XC
        PC=PC+2
        GO TO 2
      ENDIF
 180  IF (X .EQ. REAL) THEN
        DIM=IRCONS(VPFX(PC+1),1)
        IF (DIM .EQ. 0) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          VQD(QC,1)=REAL
          VQD(QC,2)=IRCONS(VPFX(PC+1),4)
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 1) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          VQD(QC,1)=REAL
          VQD(QC,2)=IRCONS(VPFX(PC+1),4)-1
          VQD(QC,3)=ISTACK(ITOS)
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          ITOS=ITOS-1
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 2) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=INUM+128
          VQD(QC,2)=1
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=IC
          ITOS=ITOS+1
          IF (ITOS .GT. ISMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=SUB+128
          VQD(QC,2)=0
          VQD(QC,3)=ISTACK(ITOS-2)
          VQD(QC,4)=ISTACK(ITOS)
          VQD(QC,5)=IC
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=INUM+128
          VQD(QC,2)=IRCONS(VPFX(PC+1),3)
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=IC
          ITOS=ITOS+1
          IF (ITOS .GT. ISMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=MULT+128
          VQD(QC,2)=0
          VQD(QC,3)=ISTACK(ITOS-1)
          VQD(QC,4)=ISTACK(ITOS)
          VQD(QC,5)=IC
          ITOS=ITOS-1
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          VQD(QC,1)=REAL
          VQD(QC,2)=IRCONS(VPFX(PC+1),4)-1
          VQD(QC,3)=ISTACK(ITOS-1)
          VQD(QC,4)=ISTACK(ITOS)
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          ITOS=ITOS-3
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ENDIF
      ENDIF
 190  IF (X .EQ. INT) THEN
        DIM=IICONS(VPFX(PC+1),1)
        IF (DIM .EQ. 0) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          VQD(QC,1)=INT
          VQD(QC,2)=IICONS(VPFX(PC+1),4)
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 1) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          VQD(QC,1)=INT
          VQD(QC,2)=IICONS(VPFX(PC+1),4)-1
          VQD(QC,3)=ISTACK(ITOS)
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          ITOS=ITOS-1
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 2) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=INUM+128
          VQD(QC,2)=1
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=IC
          ITOS=ITOS+1
          IF (ITOS .GT. ISMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=SUB+128
          VQD(QC,2)=0
          VQD(QC,3)=ISTACK(ITOS-2)
          VQD(QC,4)=ISTACK(ITOS)
          VQD(QC,5)=IC
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=INUM+128
          VQD(QC,2)=IICONS(VPFX(PC+1),3)
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=IC
          ITOS=ITOS+1
          IF (ITOS .GT. ISMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          IC=IC+1
          VQD(QC,1)=MULT+128
          VQD(QC,2)=0
          VQD(QC,3)=ISTACK(ITOS-1)
          VQD(QC,4)=ISTACK(ITOS)
          VQD(QC,5)=IC
          ITOS=ITOS-1
          ISTACK(ITOS)=IC
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          VQD(QC,1)=INT
          VQD(QC,2)=IICONS(VPFX(PC+1),4)-1
          VQD(QC,3)=ISTACK(ITOS-1)
          VQD(QC,4)=ISTACK(ITOS)
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          ITOS=ITOS-3
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ENDIF
      ENDIF
 200  IF (X .EQ. VAR) THEN
        DIM=IVARI(VPFX(PC+1),1)
        IF (DIM .EQ. 0) THEN
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=IVARI(VPFX(PC+1),3)
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 1) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          IC=IC+1
          VQD(QC,1)=VAR
          VQD(QC,2)=IVARI(VPFX(PC+1),3)-1
          VQD(QC,3)=ISTACK(ITOS)
          VQD(QC,4)=IC
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          ITOS=ITOS-1
          PC=PC+2
          GO TO 2
        ENDIF
      ENDIF
 210  IF (X .EQ. FUNC) THEN
        DIM=IFUNC(VPFX(PC+1),1)
        IF (DIM .LE. 0) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          VQD(QC,1)=FUNC
          VQD(QC,2)=IFUNC(VPFX(PC+1),4)
          VQD(QC,3)=IFUNC(VPFX(PC+1),1)
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 1) THEN
          QC=QC+1
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          XC=XC+1
          IC=IC+1
          VQD(QC,1)=FUNC
          VQD(QC,2)=IFUNC(VPFX(PC+1),4)-1
          VQD(QC,3)=ISTACK(ITOS)
          VQD(QC,4)=IC
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          ITOS=ITOS-1
          PC=PC+2
          GO TO 2
        ENDIF
      ENDIF
 220  IF (X .EQ. STDRD) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        DIM=STDTYP(VPFX(PC+1))
        IF (DIM .EQ. 0) THEN
          VQD(QC,1)=STDRD
          VQD(QC,2)=VPFX(PC+1)
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 1) THEN
          VQD(QC,1)=STDRD
          VQD(QC,2)=VPFX(PC+1)
          VQD(QC,3)=XSTACK(XTOS)
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 2) THEN
          VQD(QC,1)=STDRD
          VQD(QC,2)=VPFX(PC+1)
          VQD(QC,3)=XSTACK(XTOS-1)
          VQD(QC,4)=XSTACK(XTOS)
          VQD(QC,5)=XC
          XTOS=XTOS-1
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ENDIF
      ENDIF
 230  IF (X .EQ. EXTERN) THEN
        EXCALL=.TRUE.
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        DIM=EXTTYP(VPFX(PC+1))
        IF (DIM .EQ. 0) THEN
          VQD(QC,1)=EXTERN
          VQD(QC,2)=VPFX(PC+1)
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 1) THEN
          VQD(QC,1)=EXTERN
          VQD(QC,2)=VPFX(PC+1)
          VQD(QC,3)=ISTACK(ITOS)
          VQD(QC,4)=0
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          ITOS=ITOS-1
          PC=PC+2
          GO TO 2
        ELSE IF (DIM .EQ. 2) THEN
          VQD(QC,1)=EXTERN
          VQD(QC,2)=VPFX(PC+1)
          VQD(QC,3)=ISTACK(ITOS-1)
          VQD(QC,4)=ISTACK(ITOS)
          VQD(QC,5)=XC
          XTOS=XTOS+1
          IF (XTOS .GT. XSMDEP) THEN
            IERR=24
            RETURN
          ENDIF
          XSTACK(XTOS)=XC
          ITOS=ITOS-2
          PC=PC+2
          GO TO 2
        ENDIF
      ENDIF
 240  IF (X .EQ. SUM) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC-1,1)=SUM+128
        VQD(QC-1,2)=0
        VQD(QC-1,3)=0
        VQD(QC-1,4)=0
        VQD(QC-1,5)=XC
        XTOS=XTOS+1
        IF (XTOS .GT. XSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        XSTACK(XTOS)=XC
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=IC
        VQD(QC,1)=SUM
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=IINDEX(VPFX(PC+2),2)
        VQD(QC,4)=IINDEX(VPFX(PC+2),5)
        RTOS=RTOS+1
        IF (RTOS .GT. RSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        RSTACK(RTOS)=QC
        PC=PC+3
        GO TO 2
      ENDIF
 250  IF (X .EQ. ENDSUM) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC-1,1)=ADD
        VQD(QC-1,2)=0
        VQD(QC-1,3)=XSTACK(XTOS-1)
        VQD(QC-1,4)=XSTACK(XTOS)
        VQD(QC-1,5)=XC
        VQD(QC,1)=ENDSUM
        VQD(QC,2)=XC-XSTACK(XTOS-1)
        VQD(QC,3)=IC-ISTACK(ITOS)
        VQD(QC,4)=0
        VQD(QC,5)=RSTACK(RTOS)
        VQD(RSTACK(RTOS),5)=QC
        XC=XC+VQD(QC,2)*(VQD(RSTACK(RTOS),3)-1)
        IC=IC+VQD(QC,3)*(VQD(RSTACK(RTOS),3)-1)
        ITOS=ITOS-1
        XTOS=XTOS-1
        XSTACK(XTOS)=XC
        RTOS=RTOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 260  IF (X .EQ. PROD) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC-1,1)=PROD+128
        VQD(QC-1,2)=0
        VQD(QC-1,3)=0
        VQD(QC-1,4)=0
        VQD(QC-1,5)=XC
        XTOS=XTOS+1
        IF (XTOS .GT. XSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        XSTACK(XTOS)=XC
        ITOS=ITOS+1
        IF (ITOS .GT. ISMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        ISTACK(ITOS)=IC
        VQD(QC,1)=PROD
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=IINDEX(VPFX(PC+2),2)
        VQD(QC,4)=IINDEX(VPFX(PC+2),5)
        RTOS=RTOS+1
        IF (RTOS .GT. RSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        RSTACK(RTOS)=QC
        PC=PC+3
        GO TO 2
      ENDIF
 270  IF (X .EQ. ENDPRD) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC-1,1)=MULT
        VQD(QC-1,2)=0
        VQD(QC-1,3)=XSTACK(XTOS-1)
        VQD(QC-1,4)=XSTACK(XTOS)
        VQD(QC-1,5)=XC
        VQD(QC,1)=ENDPRD
        VQD(QC,2)=XC-XSTACK(XTOS-1)
        VQD(QC,3)=IC-ISTACK(ITOS)
        VQD(QC,4)=0
        VQD(QC,5)=RSTACK(RTOS)
        VQD(RSTACK(RTOS),5)=QC
        XC=XC+VQD(QC,2)*(VQD(RSTACK(RTOS),3)-1)
        IC=IC+(VQD(RSTACK(RTOS),3)-1)*VQD(QC,3)
        ITOS=ITOS-1
        XTOS=XTOS-1
        XSTACK(XTOS)=XC
        RTOS=RTOS-1
        PC=PC+1
        GO TO 2
      ENDIF
 280  IF (X .EQ. RELOP) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        LC=LC+1
        VQD(QC,1)=RELOP
        VQD(QC,2)=VPFX(PC+1)
        VQD(QC,3)=XSTACK(XTOS-1)
        VQD(QC,4)=XSTACK(XTOS)
        VQD(QC,5)=LC
        LTOS=LTOS+1
        IF (LTOS .GT. LSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        LSTACK(LTOS)=LC
        XTOS=XTOS-2
        IF (XC .GT. MAXXC) MAXXC=XC
        IF (IC .GT. MAXIC) MAXIC=IC
        XC=PVVA
        PC=PC+2
        GO TO 2
      ENDIF
 290  IF (X .EQ. AND) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        LC=LC+1
        VQD(QC,1)=AND
        VQD(QC,2)=0
        VQD(QC,3)=LSTACK(LTOS-1)
        VQD(QC,4)=LSTACK(LTOS)
        VQD(QC,5)=LC
        LTOS=LTOS-1
        LSTACK(LTOS)=LC
        PC=PC+1
        GO TO 2
      ENDIF
 300  IF (X .EQ. OR) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        LC=LC+1
        VQD(QC,1)=OR   
        VQD(QC,2)=0          
        VQD(QC,3)=LSTACK(LTOS-1)
        VQD(QC,4)=LSTACK(LTOS)
        VQD(QC,5)=LC
        LTOS=LTOS-1
        LSTACK(LTOS)=LC
        PC=PC+1
        GO TO 2
      ENDIF
 310  IF (X .EQ. NOT) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        LC=LC+1
        VQD(QC,1)=NOT
        VQD(QC,2)=0
        VQD(QC,3)=LSTACK(LTOS)
        VQD(QC,4)=0
        VQD(QC,5)=LC
        LSTACK(LTOS)=LC
        PC=PC+1
        GO TO 2
      ENDIF
 320  IF (X .EQ. BEQ) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        VQD(QC,1)=BEQ
        VQD(QC,2)=0
        VQD(QC,3)=LSTACK(LTOS)
        VQD(QC,4)=0
        VQD(QC,5)=0
        LTOS=LTOS-1
        IF (LC .GT. MAXLC) MAXLC=LC
        LC=0
        PC=PC+2
        GO TO 2
      ENDIF
 330  IF (X .EQ. BRA) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        VQD(QC-1,1)=BRA
        VQD(QC-1,2)=0
        VQD(QC-1,3)=0
        VQD(QC-1,4)=0
        VQD(QC-1,5)=0
        VQD(QC,1)=LABEL
        VQD(QC,2)=QC
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=0
        PC=PC+2
        GO TO 2
      ENDIF
 340  IF (X .EQ. CONTIN) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        VQD(QC,1)=CONTIN
        VQD(QC,2)=0
        VQD(QC,3)=VPFX(PC+1)
        VQD(QC,4)=0
        VQD(QC,5)=0
        PC=PC+2
        GO TO 2
      ENDIF
 350  IF (X .EQ. GOTO) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        VQD(QC,1)=GOTO
        VQD(QC,2)=0
        VQD(QC,3)=VPFX(PC+1)
        VQD(QC,4)=0
        VQD(QC,5)=0
        PC=PC+3
        GO TO 2
      ENDIF
 360  IF (X .EQ. ASSIGN) THEN
        IF (IFUNC(VPFX(PC+1),1) .LE. 0) THEN
          QC=QC+2
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          VQD(QC-1,1)=ASSIGN
          VQD(QC-1,2)=IFUNC(VPFX(PC+1),4)
          VQD(QC-1,3)=XSTACK(XTOS)
          VQD(QC-1,4)=0
          VQD(QC-1,5)=0
          IF (IFUNC(VPFX(PC+1),1) .EQ. -2) VQD(QC-1,4)=-2
          XTOS=XTOS-1
          IF (XC .GT. MAXXC) MAXXC=XC
          IF (IC .GT. MAXIC) MAXIC=IC
          XC=PVVA
          IC=0
          VQD(QC,1)=LABEL
          VQD(QC,2)=QC
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=0
          PC=PC+2
          GO TO 2
        ELSE
          QC=QC+2
          IF (QC .GT. MAXVQD) THEN
            IERR=32
            RETURN
          ENDIF
          VQD(QC-1,1)=ASSIGN
          VQD(QC-1,2)=-IFUNC(VPFX(PC+1),2)
          VQD(QC-1,3)=XSTACK(XTOS)
          VQD(QC-1,4)=IFUNC(VPFX(PC+1),4)-1
          VQD(QC-1,5)=0
          XTOS=XTOS-1
          IF (XC .GT. MAXXC) MAXXC=XC
          IF (IC .GT. MAXIC) MAXIC=IC
          XC=PVVA
          IC=0
          VQD(QC,1)=LABEL
          VQD(QC,2)=QC
          VQD(QC,3)=0
          VQD(QC,4)=0
          VQD(QC,5)=0
          PC=PC+2
          GO TO 2
        ENDIF
      ENDIF
 365  IF (X .EQ. -ASSIGN) THEN
        IERR=27
        RETURN
      ENDIF
 370  IF (X .EQ. IF+128) THEN
        RTOS=RTOS+1
        IF (RTOS .GT. RSMDEP) THEN
          IERR=24
          RETURN
        ENDIF
        RSTACK(RTOS)=QC+1
        PC=PC+1
        GO TO 2
      ENDIF
 380  IF (X .EQ. ENDIF+128) THEN
        QC=QC+1
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        VQD(QC,1)=LABEL
        VQD(QC,2)=QC
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=0
        QC1=RSTACK(RTOS)
        RTOS=RTOS-1
 381    IF (QC1 .LT. QC) THEN
          IF (VQD(QC1,1) .NE. BEQ) THEN
            QC1=QC1+1
            GO TO 381
          ENDIF
        ENDIF
        QC2=QC1+1
 382    IF (QC2 .LT. QC) THEN
          IF (VQD(QC2,1) .NE. BRA) THEN
            QC2=QC2+1
            GO TO 382
          ENDIF
        ENDIF
        IF (QC2 .LT. QC) THEN
          VQD(QC1,2)=QC2+1
          VQD(QC2,2)=QC
          QC1=QC1+1
          GO TO 381
        ELSE
          IF (VQD(QC1,1) .EQ. BEQ) VQD(QC1,2)=QC
          PC=PC+1
          GO TO 2
        ENDIF
      ENDIF
 390  IF (X .EQ. CONINT) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC-1,1)=CONINT
        VQD(QC-1,2)=IRCONS(VPFX(PC+1),4)
        VQD(QC-1,3)=IRCONS(VPFX(PC+1),2)
        VQD(QC-1,4)=VPFX(PC+2)
        VQD(QC-1,5)=XC
        VQD(QC,1)=CONINT+128
        VQD(QC,2)=0
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=XSTACK(XTOS)
        XSTACK(XTOS)=XC
        PIP=PIP+1
        IF (PIP .GT. MAXVIP) THEN
          IERR=64
          RETURN
        ENDIF
        VIP(PIP)=VPFX(PC+2)
        PC=PC+3
        GO TO 2
      ENDIF
 400  IF (X .EQ. LININT) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC-1,1)=LININT
        VQD(QC-1,2)=IRCONS(VPFX(PC+1),4)
        VQD(QC-1,3)=IRCONS(VPFX(PC+1),2)
        VQD(QC-1,4)=VPFX(PC+2)
        VQD(QC-1,5)=XC
        VQD(QC,1)=LININT+128
        VQD(QC,2)=0
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=XSTACK(XTOS)
        XSTACK(XTOS)=XC
        PIP=PIP+1
        IF (PIP .GT. MAXVIP) THEN
          IERR=64
          RETURN
        ENDIF
        VIP(PIP)=VPFX(PC+2)
        PC=PC+3
        GO TO 2
      ENDIF
 410  IF (X .EQ. SPLINE) THEN
        QC=QC+2
        IF (QC .GT. MAXVQD) THEN
          IERR=32
          RETURN
        ENDIF
        XC=XC+1
        VQD(QC-1,1)=SPLINE
        VQD(QC-1,2)=IRCONS(VPFX(PC+1),4)
        VQD(QC-1,3)=IRCONS(VPFX(PC+1),2)
        VQD(QC-1,4)=VPFX(PC+2)
        VQD(QC-1,5)=XC
        VQD(QC,1)=SPLINE+128
        VQD(QC,2)=0
        VQD(QC,3)=0
        VQD(QC,4)=0
        VQD(QC,5)=XSTACK(XTOS)
        XSTACK(XTOS)=XC
        PIP=PIP+1
        IF (PIP .GT. MAXVIP) THEN
          IERR=64
          RETURN
        ENDIF
        VIP(PIP)=VPFX(PC+2)
        PC=PC+3
        GO TO 2
      ENDIF
C
      IERR=26
      WRITE(*,*) 'REVCDE (4866) : unknown token ',X
      RETURN 
      END 
C
C
C
      SUBROUTINE FORCDE (XFIL,PVIS,PVIC,PVRC,PVVA,PIFN,PXFN,PVFN,PIV,
     1                   VINDEX,VICONS,VRCONS,MAXVQD,VQD,MAXXC,MAXLC,
     2                   MAXIC,EXCALL,IERR,VIP,MAXVIP,PIP)
      INTEGER XFIL
      INTEGER PVIS,PVIC,PVRC,PVVA,PIFN,PXFN,PVFN,PIV
      INTEGER VINDEX(PVIS),VICONS(PVIC)
      DOUBLE PRECISION VRCONS(PVRC)
      INTEGER MAXVQD
      INTEGER VQD(MAXVQD,5)
      INTEGER MAXXC,MAXLC,MAXIC
      LOGICAL EXCALL
      INTEGER IERR
      INTEGER MAXVIP,PIP
      INTEGER VIP(MAXVIP)
C
      INTEGER ADD,SUB,MULT,DIV,POWER,LEFT,RIGHT,COMMA,ASSIGN,NLINE
      INTEGER RANGE,RELOP,AND,OR,NOT,INUM,RNUM,ID,SUM,PROD,IN,IF,THEN
      INTEGER ELSE,ENDIF,STDRD,EXTERN,INTERP,PARAM,INDEX,REAL,INT,TABLE
      INTEGER CONINT,LININT,SPLINE,VAR,INFUNC,FUNC,END,GOTO,MARKE,CONTIN
      INTEGER UMINUS,INDVAR,ENDSUM,ENDPRD,BEQ,BRA,LABEL,VECTOR,ACTIVE
      PARAMETER (ADD=43,SUB=45,MULT=42,DIV=47,POWER=94,LEFT=40,RIGHT=41)
      PARAMETER (COMMA=44,ASSIGN=61,NLINE=10,RANGE=257,RELOP=258)
      PARAMETER (AND=259,OR=260,NOT=261,INUM=262,RNUM=263,ID=264)
      PARAMETER (SUM=265,PROD=266,IN=267,IF=268,THEN=269,ELSE=270)
      PARAMETER (ENDIF=271,STDRD=272,EXTERN=273,INTERP=274,PARAM=275)
      PARAMETER (INDEX=276,REAL=277,INT=278,TABLE=279,CONINT=280)
      PARAMETER (LININT=281,SPLINE=282,VAR=283,INFUNC=284,FUNC=285)
      PARAMETER (END=286,GOTO=287,MARKE=288,CONTIN=289,UMINUS=290)
      PARAMETER (INDVAR=291,ENDSUM=292,ENDPRD=293,BEQ=294,BRA=295)
      PARAMETER (LABEL=296,VECTOR=297,ACTIVE=298)
C
      INTEGER MAXSTD,MAXEXT
      PARAMETER (MAXSTD=17,MAXEXT=1)
C
      CHARACTER*6 STDNAM(MAXSTD)
      INTEGER STDTYP(MAXSTD),EXTTYP(MAXEXT)
      INTEGER STDLEN(MAXSTD)
      CHARACTER*4 RELNAM(6)
C
      INTEGER FUNMOD,GRAMOD
      PARAMETER (FUNMOD=0,GRAMOD=1)
C
      INTEGER I,J,K,MODE
      INTEGER AQC,FC,FCOFS,QC,X
      INTEGER DIM,EXT
C     INTEGER B
C
C     INTEGER OFS,IOFS
      INTEGER OFS
C
      CHARACTER*30 S1,S2,S3,S4,S5,S6,S7,S8,S9,S10
      INTEGER SLEN1,SLEN2,SLEN3,SLEN4,SLEN5,SLEN6,SLEN7,SLEN8,SLEN9
      INTEGER SLEN10
C
      DATA STDNAM /'DABS','DSQRT','DEXP','DLOG','DLOG10',
     1             'DSIN','DCOS','DTAN','DASIN','DACOS','DATAN',
     2             'DSINH','DCOSH','DTANH','DASINH','DACOSH','DATANH'/
      DATA STDTYP /1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
      DATA STDLEN /4,5,4,4,6,4,4,4,5,5,5,5,5,5,6,6,6/
      DATA EXTTYP /0/
      DATA RELNAM /'.EQ.','.NE.','.LT.','.LE.','.GT.','.GE.'/
C
      OFS=0
C     IOFS=0
      MODE=FUNMOD
 2    IF (MODE .EQ. FUNMOD) THEN 
        WRITE(XFIL,'(A)') 'C*********************************'
        WRITE(XFIL,'(A)') 'C'
        WRITE(XFIL,'(A)') 'C   P C O M P     (Version 5.3)'
        WRITE(XFIL,'(A)') 'C'
        WRITE(XFIL,'(A)') 'C*********************************'
        WRITE(XFIL,'(A)') 'C'
        WRITE(XFIL,'(6X,A)') 'SUBROUTINE XFUN (X,N,F,M,ACTIVE,IERR)'
        WRITE(XFIL,'(6X,A)') 'INTEGER N,M'
        WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION X(N),F(M)'
        WRITE(XFIL,'(6X,A)') 'LOGICAL ACTIVE(M)' 
        WRITE(XFIL,'(6X,A)') 'INTEGER IERR'
        WRITE(XFIL,'(A)') 'C'
C
        IF (MAXXC .GT. PVVA) THEN
          CALL ITOA(PVVA+1,S1,SLEN1)
          CALL ITOA(MAXXC,S2,SLEN2)
          WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION XAUX('//S1(1:SLEN1)//
     1               ':'//S2(1:SLEN2)//')'
        ENDIF
C
        IF (MAXIC .GT. 0) THEN
          CALL ITOA(MAXIC,S1,SLEN1)
          WRITE(XFIL,'(6X,A)') 'INTEGER IAUX(1:'//S1(1:SLEN1)//')'
        ENDIF
C
        IF (PXFN .LT. PIFN) THEN
          CALL ITOA(PVFN-(PIFN-PXFN)+1,S1,SLEN1)
          CALL ITOA(PVFN,S2,SLEN2)
          WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION VFUNC('//S1(1:SLEN1)//
     1           ':'//S2(1:SLEN2)//')'
        ENDIF
      ELSE 
        WRITE(XFIL,'(A)') 'C'
        WRITE(XFIL,'(A)') 'C'
        WRITE(XFIL,'(A)') 'C'
        WRITE(XFIL,'(6X,A)') 
     1      'SUBROUTINE XGRA (X,N,F,M,DF,MMAX,ACTIVE,IERR)'
        WRITE(XFIL,'(6X,A)') 'INTEGER N,M,MMAX'
        WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION X(N),F(M),DF(MMAX,N)'
        WRITE(XFIL,'(6X,A)') 'LOGICAL ACTIVE(M)' 
        WRITE(XFIL,'(6X,A)') 'INTEGER IERR'
        WRITE(XFIL,'(A)') 'C'
C
        CALL ITOA(PVVA,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION DFHELP('//
     1                        S1(1:SLEN1)//')'
        WRITE(XFIL,'(A)') 'C'
C
        IF (MAXXC .GT. PVVA) THEN
          CALL ITOA(PVVA+1,S1,SLEN1)
          CALL ITOA(MAXXC,S2,SLEN2)
          WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION XAUX('//S1(1:SLEN1)//
     1               ':'//S2(1:SLEN2)//'),YAUX('//S1(1:SLEN1)//
     2               ':'//S2(1:SLEN2)//')'
        ENDIF
        IF (EXCALL) THEN
          CALL ITOA(PVVA,S1,SLEN1)
          WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION Z(1:'//S1(1:SLEN1)//')'
        ENDIF
C
        IF (MAXIC .GT. 0) THEN
          CALL ITOA(MAXIC,S1,SLEN1)
          WRITE(XFIL,'(6X,A)') 'INTEGER IAUX(1:'//S1(1:SLEN1)//')'
        ENDIF
C
        IF (PXFN .LT. PIFN) THEN
          CALL ITOA(PVFN-(PIFN-PXFN)+1,S1,SLEN1)
          CALL ITOA(PVFN,S2,SLEN2)
          CALL ITOA(PVVA,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION VFUNC('//S1(1:SLEN1)//
     1           ':'//S2(1:SLEN2)//'),VGRAD('//S1(1:SLEN1)//
     2           ':'//S2(1:SLEN2)//','//S3(1:SLEN3)//')'
        ENDIF
      ENDIF
C
      IF (MAXLC .GT. 0) THEN
        CALL ITOA(MAXLC,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'LOGICAL LAUX(1:'//S1(1:SLEN1)//')'
      ENDIF
      IF (EXCALL) WRITE(XFIL,'(6X,A)') 'INTEGER EXTPAR(2)'
C---------------------------------------------------------------
      DO 10 I=0,PIV ! -1
        CALL ITOA(I,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'INTEGER I'//S1(1:SLEN1)//',IX'//
     1          S1(1:SLEN1)
 10   CONTINUE
      IF (MODE .EQ. GRAMOD) THEN
        WRITE(XFIL,'(6X,A)') 'INTEGER INIT1,INIT2'
C
        IF (PIP .GT. 0) THEN
          DO 15 I=1,PIP
          CALL ITOA(VIP(I),S1,SLEN1)
          WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION IR'//S1(1:SLEN1)
 15       CONTINUE
        ENDIF
C
      ENDIF
C---------------------------------------------------------------
      WRITE(XFIL,'(6X,A)') 'INTEGER I,OFS,IOFS'
      WRITE(XFIL,'(A)') 'C'
C
      CALL ITOA(PVRC,S3,SLEN3)
      IF (PVIS .GT. 0) THEN
        CALL ITOA(PVIS,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'INTEGER VINDEX('//S1(1:SLEN1)//')'
      ENDIF
      IF (PVIC .GT. 0) THEN
        CALL ITOA(PVIC,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'INTEGER VICONS('//S1(1:SLEN1)//')'
      ENDIF
      IF (PVRC .GT. 0) THEN
        CALL ITOA(PVRC,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'DOUBLE PRECISION VRCONS('//
     1          S1(1:SLEN1)//')'
      ENDIF
      IF (PVIS .GT. 0) CALL IDATA(XFIL,VINDEX,PVIS,'VINDEX')
      IF (PVIC .GT. 0) CALL IDATA(XFIL,VICONS,PVIC,'VICONS')
      IF (PVRC .GT. 0) CALL FDATA(XFIL,VRCONS,PVRC,'VRCONS')
      WRITE(XFIL,'(A)') 'C'
C
C-----------------------------------------------------------------------
C     IF (MODE .EQ. GRAMOD) THEN
C       B=1
C       WRITE(XFIL,'(6X,A)') 'DO 1 INIT1=1,MMAX'
C       WRITE(XFIL,'(6X,A)') 'DO 1 INIT2=1,N'
C       WRITE(XFIL,'(6X,A)') 'DF(INIT1,INIT2)=0.0D0'
C       WRITE(XFIL,'(I5,1X,A)') B,'CONTINUE'
C       IF (PXFN .LT. PIFN) THEN
C         B=2
C         CALL ITOA(PVFN-(PIFN-PXFN)+1,S1,SLEN1)
C         CALL ITOA(PVFN,S2,SLEN2)
C         WRITE(XFIL,'(6X,A)') 'DO 2 INIT1='//S1(1:SLEN1)//','//
C    1                                        S2(1:SLEN2)
C         WRITE(XFIL,'(6X,A)') 'DO 2 INIT2=1,N'
C         WRITE(XFIL,'(6X,A)') 'VGRAD(INIT1,INIT2)=0.0D0'
C         WRITE(XFIL,'(I5,1X,A)') B,'CONTINUE'
C       ENDIF
C     ENDIF
      IF (MODE.EQ.GRAMOD) THEN
        CALL ITOA(PVFN-(PIFN-PXFN),S1,SLEN1)
        CALL ITOA(PVVA,S3,SLEN3)
	WRITE(XFIL,'(6X,A)') 'CALL XINI(DF,1,MMAX,'//
     &    S3(1:SLEN3)//')'
        IF (PXFN.LT.PIFN) THEN
          CALL ITOA(PVFN-(PIFN-PXFN)+1,S1,SLEN1)
          CALL ITOA(PVFN,S2,SLEN2)
	  WRITE(XFIL,'(6X,A)') 'CALL XINI(VGRAD,'//S1(1:SLEN1)//','//
     &      S2(1:SLEN2)//','//S3(1:SLEN3)//')'
	ENDIF
      ENDIF
C-----------------------------------------------------------------------
C    
      DO 20 I=0,PIV ! -1
        CALL ITOA(I,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'IX'//S1(1:SLEN1)//'=0'
 20   CONTINUE
      CALL ITOA(PVVA,S1,SLEN1)
      CALL ITOA(PVFN-(PIFN-PXFN),S2,SLEN2)
      WRITE(XFIL,'(6X,A)') 'IF (N .NE. '//S1(1:SLEN1)//') THEN'
      WRITE(XFIL,'(6X,A)') 'IERR=43'
      WRITE(XFIL,'(6X,A)') 'RETURN'
      WRITE(XFIL,'(6X,A)') 'ENDIF'
      WRITE(XFIL,'(6X,A)') 'IF (M .NE. '//S2(1:SLEN2)//') THEN'
      WRITE(XFIL,'(6X,A)') 'IERR=44'
      WRITE(XFIL,'(6X,A)') 'RETURN'
      WRITE(XFIL,'(6X,A)') 'ENDIF'
      WRITE(XFIL,'(6X,A)') 'OFS=0'
C
      WRITE(XFIL,'(6X,A)') 'IOFS=0'
C----------------------------------------------------------------
      K=1
      QC=3
C----------------------------------------------------------------
 100  X=VQD(QC,1)
C
      IF (X .EQ. ACTIVE) THEN
        IF (VQD(QC,2) .GE. 0) THEN
          CALL ITOA(VQD(QC,2),S1,SLEN1)
          WRITE(XFIL,'(6X,A)') 'IF (ACTIVE('//S1(1:SLEN1)//')) THEN'
        ELSE
          CALL ITOA(VQD(QC,3),S1,SLEN1)
          CALL ITOA(-VQD(QC,2)-1,S2,SLEN2)
          WRITE(XFIL,'(6X,A)') 'IF (ACTIVE('//S1(1:SLEN1)//'+IX'//
     1                         S2(1:SLEN2)//')) THEN'
        ENDIF
        QC=QC+1
        GO TO 100
C
      ELSE IF (X .EQ. -1) THEN
        WRITE(XFIL,'(6X,A)') 'ENDIF'
C----------------------------------------------------------------
        IF (VQD(QC,2) .GT. 0) THEN
          WRITE(XFIL,'(I5,1X,A)') VQD(QC,2),'CONTINUE'
        ENDIF
C----------------------------------------------------------------
        K=K+1
        IF (K .LE. PXFN) THEN 
          QC=QC+1
          GO TO 100
        ELSE
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'END'
          IF (MODE .EQ. FUNMOD) THEN
            MODE=GRAMOD
            GO TO 2
          ELSE
C-----------------------------------------------------------------------
            WRITE(XFIL,'(A)') 'C'
            WRITE(XFIL,'(A)') 'C'
            WRITE(XFIL,'(A)') 'C'
            WRITE(XFIL,'(A)') '      SUBROUTINE XINI (G,ML,MU,N)'
            WRITE(XFIL,'(A)') '      INTEGER ML,MU,N'
            WRITE(XFIL,'(A)') '      DOUBLE PRECISION G(ML:MU,N)'
            WRITE(XFIL,'(A)') 'C'
            WRITE(XFIL,'(A)') '      INTEGER I,J'
            WRITE(XFIL,'(A)') 'C'
            WRITE(XFIL,'(A)') '      DO 20 I=ML,MU'
            WRITE(XFIL,'(A)') '      DO 10 J=1,N'
            WRITE(XFIL,'(A)') '      G(I,J)=0.0D0'
            WRITE(XFIL,'(A)') '   10 CONTINUE'
            WRITE(XFIL,'(A)') '   20 CONTINUE'
            WRITE(XFIL,'(A)') '      RETURN'
            WRITE(XFIL,'(A)') '      END'
C-----------------------------------------------------------------------
            RETURN
          ENDIF
        ENDIF
C
      ELSE IF (X .EQ. ADD+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEI(OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEI(OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)//'+'//
     1        S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. SUB+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEI(OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEI(OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)//'-'//
     1        S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. MULT+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEI(OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEI(OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)//'*'//
     1        S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. DIV+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEI(OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEI(OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)//' .EQ. 0.0D0) THEN'
        WRITE(XFIL,'(6X,A)') 'IERR=9'
        WRITE(XFIL,'(6X,A)') 'RETURN'
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)//'/'//
     1        S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. UMINUS+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEI(OFS,VQD(QC,3),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=-'//S2(1:SLEN2)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. INUM+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL ITOA(VQD(QC,2),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. FUNC+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL ITOA(VQD(QC,2),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=IDNINT(VFUNC('//
     1        S2(1:SLEN2)//'))'
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. INDVAR+128) THEN
        CALL MAKEI(OFS,VQD(QC,5),S1,SLEN1)
        CALL ITOA(VQD(QC,2)-1,S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=IX'//S2(1:SLEN2)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. VAR) THEN
        CALL MAKEI(OFS,VQD(QC,4),S1,SLEN1)
        CALL ITOA(VQD(QC,2),S2,SLEN2)
        CALL MAKEI(OFS,VQD(QC,3),S3,SLEN3)
        IF (VQD(QC,2) .NE. 0) THEN
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)//
     1          '+'//S3(1:SLEN3)
        ELSE
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S3(1:SLEN3)
        ENDIF
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S4,SLEN4)
        CALL ITOA(PVVA,S5,SLEN5)
        WRITE(XFIL,'(6X,A)') 'IF (('//S1(1:SLEN1)//' .GT. '
     1        //S5(1:SLEN5)//') .OR. ('//S3(1:SLEN3)//' .LE. 0)) THEN'
        WRITE(XFIL,'(6X,A)') 'IERR=33'
        WRITE(XFIL,'(6X,A)') 'RETURN'
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(6X,A)') S4(1:SLEN4)//'=X('//S1(1:SLEN1)//')'
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. ADD) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'+'//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. SUB) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'-'//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. MULT) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'*'//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. DIV) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)//' .EQ. 0.0D0) THEN'
        WRITE(XFIL,'(6X,A)') 'IERR=9'
        WRITE(XFIL,'(6X,A)') 'RETURN'
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'/'//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. POWER) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') 'IF ('//S2(1:SLEN2)//' .LE. 0.0D0) THEN'
        WRITE(XFIL,'(6X,A)') 'IERR=9'
        WRITE(XFIL,'(6X,A)') 'RETURN'
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'**'//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. POWER+128) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL ITOA(VQD(QC,2),S3,SLEN3)
        IF (VQD(QC,2) .LT. 0) THEN
          WRITE(XFIL,'(6X,A)') 'IF ('//S2(1:SLEN2)
     1                         //' .EQ. 0.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR=9'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
        ENDIF
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'**'//'('//S3(1:SLEN3)//')'
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. UMINUS) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            '-'//S2(1:SLEN2)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. RNUM) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL FTOA(VRCONS(VQD(QC,2)),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. REAL) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        IF (VQD(QC,3) .EQ. 0) THEN
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//
     1            S2(1:SLEN2)//')'
        ELSE IF (VQD(QC,4) .EQ. 0) THEN
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          CALL MAKEI(OFS,VQD(QC,3),S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//
     1            S2(1:SLEN2)//'+'//S3(1:SLEN3)//')'
        ELSE
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          CALL MAKEI(OFS,VQD(QC,3),S3,SLEN3)
          CALL MAKEI(OFS,VQD(QC,4),S4,SLEN4)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//
     1            S2(1:SLEN2)//'+'//S3(1:SLEN3)//'+'//S4(1:SLEN4)//')'
        ENDIF
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. INUM) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL ITOA(VICONS(VQD(QC,2)),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=DBLE('//
     1            S2(1:SLEN2)//')'
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. INT) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        IF (VQD(QC,3) .EQ. 0) THEN
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=DBLE(VICONS('//
     1            S2(1:SLEN2)//'))'
        ELSE IF (VQD(QC,4) .EQ. 0) THEN
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          CALL MAKEI(OFS,VQD(QC,3),S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=DBLE(VICONS('//
     1            S2(1:SLEN2)//'+'//S3(1:SLEN3)//'))'
        ELSE
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          CALL MAKEI(OFS,VQD(QC,3),S3,SLEN3)
          CALL MAKEI(OFS,VQD(QC,4),S4,SLEN4)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=DBLE(VICONS('//
     1            S2(1:SLEN2)//'+'//S3(1:SLEN3)//'+'//S4(1:SLEN4)//'))'
        ENDIF
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. INDVAR) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL ITOA(VQD(QC,2)-1,S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=DBLE(IX'//
     1            S2(1:SLEN2)//')'
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. FUNC) THEN
        IF (VQD(QC,4) .EQ. 0) THEN
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          CALL MAKEX(PVVA,OFS,VQD(QC,5),S4,SLEN4)
          IF (VQD(QC,3) .EQ. 0) THEN
            WRITE(XFIL,'(6X,A)') S4(1:SLEN4)//'=F('//S2(1:SLEN2)//')'
          ELSE
            WRITE(XFIL,'(6X,A)') S4(1:SLEN4)//'=VFUNC('//
     1                           S2(1:SLEN2)//')'
          ENDIF
        ELSE
          CALL MAKEI(OFS,VQD(QC,4),S1,SLEN1)
          CALL ITOA(VQD(QC,2),S2,SLEN2)
          CALL MAKEI(OFS,VQD(QC,3),S3,SLEN3)
          IF (VQD(QC,2) .NE. 0) THEN
            WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)//
     1            '+'//S3(1:SLEN3)
          ELSE
            WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S3(1:SLEN3)
          ENDIF
          CALL MAKEX(PVVA,OFS,VQD(QC,5),S4,SLEN4)
          CALL ITOA(PVFN-(PIFN-PXFN),S5,SLEN5)
          WRITE(XFIL,'(6X,A)') 'IF (('//S1(1:SLEN1)//' .GT. '//
     1          S5(1:SLEN5)//') .OR. ('//S3(1:SLEN3)//' .LE. 0)) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR=33'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S4(1:SLEN4)//'=F('//S1(1:SLEN1)//')'
        ENDIF
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. STDRD) THEN
        DIM=STDTYP(VQD(QC,2))
C
        IF (DIM .EQ. 0) THEN 
          CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
          S2=STDNAM(VQD(QC,2))
          SLEN2=STDLEN(VQD(QC,2))
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)
C
        ELSE IF (DIM .EQ. 1) THEN
          CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
          S2=STDNAM(VQD(QC,2))
          SLEN2=STDLEN(VQD(QC,2))
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S3,SLEN3)
C----------------------------------------------------------------------
          IF (VQD(QC,2) .EQ. 2) THEN
            WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)
     1                         //' .LT. 0.0D0) THEN'
            WRITE(XFIL,'(6X,A)') 'IERR=53'
            WRITE(XFIL,'(6X,A)') 'RETURN'
            WRITE(XFIL,'(6X,A)') 'ENDIF'
C
          ELSE IF ((VQD(QC,2) .EQ. 4) .OR. (VQD(QC,2) .EQ. 5)) THEN
            WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)
     1                           //' .LE. 0.0D0) THEN'
            WRITE(XFIL,'(6X,A)') 'IERR=52'
            WRITE(XFIL,'(6X,A)') 'RETURN'
            WRITE(XFIL,'(6X,A)') 'ENDIF'
C
          ELSE IF (VQD(QC,2) .EQ. 9) THEN
            WRITE(XFIL,'(6X,A)') 'IF (DABS('//S3(1:SLEN3)//
     1                            ') .GT. 1.0D0) THEN'
            WRITE(XFIL,'(6X,A)') 'IERR = 54'
            WRITE(XFIL,'(6X,A)') 'RETURN'
            WRITE(XFIL,'(6X,A)') 'ENDIF'
C 
          ELSE IF(VQD(QC,2) .EQ. 10) THEN
            WRITE(XFIL,'(6X,A)') 'IF (DABS('//S3(1:SLEN3)//
     1                            ') .GT. 1.0D0) THEN'
            WRITE(XFIL,'(6X,A)') 'IERR = 55'
            WRITE(XFIL,'(6X,A)') 'RETURN'
            WRITE(XFIL,'(6X,A)') 'ENDIF'
C
          ELSE IF(VQD(QC,2) .EQ. 16) THEN
            WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)//
     1                            ' .LT. 1.0D0) THEN'
            WRITE(XFIL,'(6X,A)') 'IERR = 56'
            WRITE(XFIL,'(6X,A)') 'RETURN'
            WRITE(XFIL,'(6X,A)') 'ENDIF'
C
          ELSE IF(VQD(QC,2) .EQ. 17) THEN
            WRITE(XFIL,'(6X,A)') 'IF (DABS('//S3(1:SLEN3)//
     1                            ') .GE. 1.0D0) THEN'
            WRITE(XFIL,'(6X,A)') 'IERR = 51'
            WRITE(XFIL,'(6X,A)') 'RETURN'
            WRITE(XFIL,'(6X,A)') 'ENDIF'
          ENDIF
C----------------------------------------------------------------------
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'('//S3(1:SLEN3)//')'
C
        ELSE IF (DIM .EQ. 2) THEN
          CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
          S2=STDNAM(VQD(QC,2))
          SLEN2=STDLEN(VQD(QC,2))
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S3,SLEN3)
          CALL MAKEX(PVVA,OFS,VQD(QC,4),S4,SLEN4)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'('//S3(1:SLEN3)//
     2            ','//S4(1:SLEN4)//')'
C
        ELSE
          IERR=26
          WRITE(*,*) ' FORCDE (5149) : illegal dimension ',DIM
          RETURN
        ENDIF
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. EXTERN) THEN
        EXT=VQD(QC,2)
        IF (EXTTYP(EXT) .EQ. 1) THEN
          CALL MAKEI(OFS,VQD(QC,3),S1,SLEN1)
          WRITE(XFIL,'(6X,A)') 'EXTPAR(1)='//S1(1:SLEN1)
        ELSE IF (EXTTYP(EXT) .EQ. 2) THEN
          CALL MAKEI(OFS,VQD(QC,3),S1,SLEN1)
          CALL MAKEI(OFS,VQD(QC,4),S2,SLEN2)
          WRITE(XFIL,'(6X,A)') 'EXTPAR(1)='//S1(1:SLEN1)
          WRITE(XFIL,'(6X,A)') 'EXTPAR(2)='//S2(1:SLEN2)
        ENDIF
        CALL ITOA(EXT,S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') 'CALL EXTFUN('//S1(1:SLEN1)//',X,N,'
     1            //S2(1:SLEN2)//',EXTPAR)'
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. SUM+128) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=0.0D0'
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. PROD+128) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=1.0D0'
        QC=QC+1 
        GO TO 100
C
      ELSE IF ((X .EQ. SUM) .OR. (X.EQ. PROD) .OR. (X .EQ. VECTOR)) THEN
        CALL ITOA(VQD(QC,2)-1,S2,SLEN2)
        CALL ITOA(VQD(QC,3)-1,S3,SLEN3)
        CALL ITOA(VQD(QC,4),S4,SLEN4)
        CALL ITOA(VQD(QC,5),S5,SLEN5)
C-----------------------------------------------------
        IF (X .NE. VECTOR) OFS=OFS+1
C-----------------------------------------------------
        WRITE(XFIL,'(6X,A)') 'DO '//S5(1:SLEN5)//' I'//S2(1:SLEN2)//
     1                     '=0,'//S3(1:SLEN3)
        WRITE(XFIL,'(6X,A)') 'IX'//S2(1:SLEN2)//'=VINDEX('//
     1                     S4(1:SLEN4)//'+I'//S2(1:SLEN2)//')'
        QC=QC+1 
        GO TO 100
C
      ELSE IF ((X .EQ. ENDSUM) .OR. (X .EQ. ENDPRD)) THEN
        CALL ITOA(VQD(QC,2),S1,SLEN1)
        CALL ITOA(VQD(QC,2)*VQD(VQD(QC,5),3),S2,SLEN2)
C-----------------------------------------------------
        OFS=OFS-1
C-----------------------------------------------------
        WRITE(XFIL,'(6X,A)') 'OFS=OFS+'//S1(1:SLEN1)
C
        CALL ITOA(VQD(QC,3),S3,SLEN3)
        CALL ITOA(VQD(QC,3)*VQD(VQD(QC,5),3),S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'IOFS=IOFS+'//S3(1:SLEN3)
C
        WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE' 
        WRITE(XFIL,'(6X,A)') 'OFS=OFS-'//S2(1:SLEN2)
C
        WRITE(XFIL,'(6X,A)') 'IOFS=IOFS-'//S4(1:SLEN4)
C
        QC=QC+1 
        GO TO 100
C----------------------------------------------------------------------
      ELSE IF (X .EQ. RELOP) THEN
        CALL MAKEL(VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//RELNAM(VQD(QC,2))//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. AND) THEN
        CALL MAKEL(VQD(QC,5),S1,SLEN1)
        CALL MAKEL(VQD(QC,3),S2,SLEN2)
        CALL MAKEL(VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'.AND.'//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. OR) THEN
        CALL MAKEL(VQD(QC,5),S1,SLEN1)
        CALL MAKEL(VQD(QC,3),S2,SLEN2)
        CALL MAKEL(VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            S2(1:SLEN2)//'.OR.'//S3(1:SLEN3)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. NOT) THEN
        CALL MAKEL(VQD(QC,5),S1,SLEN1)
        CALL MAKEL(VQD(QC,3),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//
     1            '.NOT.'//S2(1:SLEN2)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. BEQ) THEN
        CALL MAKEL(VQD(QC,3),S1,SLEN1)
        CALL ITOA(VQD(QC,2),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') 'IF (.NOT.'//S1(1:SLEN1)//')'//
     1            ' GO TO '//S2(1:SLEN2)
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. BRA) THEN
        CALL ITOA(VQD(QC,2),S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'GO TO '//S1(1:SLEN1) 
        QC=QC+1 
        GO TO 100
C
      ELSE IF ((X .EQ. IF+128) .OR. (X .EQ. ENDIF+128)) THEN
        QC=QC+1 
        GO TO 100
C
      ELSE IF (X .EQ. CONTIN) THEN
        WRITE(XFIL,'(1X,I4,1X,A)') VQD(QC,3),'CONTINUE'
        QC=QC+1
        GO TO 100
C
      ELSE IF (X .EQ. GOTO) THEN
        CALL ITOA(VQD(QC,3),S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'GO TO '//S1(1:SLEN1) 
        QC=QC+1
        GO TO 100
C
      ELSE IF (X .EQ. CONINT) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC+1,5),S2,SLEN2)
        CALL ITOA(VQD(QC,3),S3,SLEN3)
        CALL ITOA(VQD(QC,2),S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)//'.EQ.0) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=0.0D0'
        WRITE(XFIL,'(6X,A)') 'ELSE IF ('//S2(1:SLEN2)//
     1            '.LT.VRCONS('//S4(1:SLEN4)//')) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=0.0D0'
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-2,S3,SLEN3)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-1,S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'ELSE IF ('//S2(1:SLEN2)//
     1            '.GE.VRCONS('//S3(1:SLEN3)//')) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//S4(1:SLEN4)//')'
        WRITE(XFIL,'(6X,A)') 'ELSE'
        CALL ITOA(QC,S5,SLEN5)
        CALL ITOA(VQD(QC,4),S6,SLEN6)
        CALL ITOA(VQD(QC,2),S3,SLEN3)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-4,S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'DO '//S5(1:SLEN5)//' I'//S6(1:SLEN6)//
     1                     '='//S3(1:SLEN3)//','//S4(1:SLEN4)//',2'
        WRITE(XFIL,'(6X,A)') 'IF (('//S2(1:SLEN2)//'.GE.VRCONS(I'//
     1              S6(1:SLEN6)//')) .AND.'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'    ('//S2(1:SLEN2)//'.LT.VRCONS(I'//
     1              S6(1:SLEN6)//'+2))) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS(I'//
     1              S6(1:SLEN6)//'+1)'
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE' 
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        QC=QC+2
        GO TO 100
C
      ELSE IF (X .EQ. LININT) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC+1,5),S2,SLEN2)
        CALL ITOA(VQD(QC,3),S3,SLEN3)
        CALL ITOA(VQD(QC,2),S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)//'.EQ.0) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=0.0D0'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S5,SLEN5)
        WRITE(XFIL,'(6X,A)') 'IR'//S5(1:SLEN5)//'=0.0D0'
        ENDIF
C
        WRITE(XFIL,'(6X,A)') 'ELSE IF ('//S2(1:SLEN2)//
     1            '.LT.VRCONS('//S4(1:SLEN4)//')) THEN'
        CALL ITOA(VQD(QC,2)+1,S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//S3(1:SLEN3)//')'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S5,SLEN5)
        WRITE(XFIL,'(6X,A)') 'IR'//S5(1:SLEN5)//'=0.0D0'
        ENDIF
C
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-2,S3,SLEN3)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-1,S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'ELSE IF ('//S2(1:SLEN2)//
     1            '.GE.VRCONS('//S3(1:SLEN3)//')) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//S4(1:SLEN4)//')'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S5,SLEN5)
        WRITE(XFIL,'(6X,A)') 'IR'//S5(1:SLEN5)//'=0.0D0'
        ENDIF
C
        WRITE(XFIL,'(6X,A)') 'ELSE'
        CALL ITOA(QC,S5,SLEN5)
        CALL ITOA(VQD(QC,4),S6,SLEN6)
        CALL ITOA(VQD(QC,2),S3,SLEN3)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-4,S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'DO '//S5(1:SLEN5)//' I'//S6(1:SLEN6)//
     1                     '='//S3(1:SLEN3)//','//S4(1:SLEN4)//',2'
        WRITE(XFIL,'(6X,A)') 'IF (('//S2(1:SLEN2)//'.GE.VRCONS(I'//
     1              S6(1:SLEN6)//')) .AND.'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'    ('//S2(1:SLEN2)//'.LT.VRCONS(I'//
     1              S6(1:SLEN6)//'+2))) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=((VRCONS(I'//S6(1:SLEN6)//
     1    '+3)-VRCONS(I'//S6(1:SLEN6)//'+1))/'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'      (VRCONS(I'//S6(1:SLEN6)//
     1                          '+2)-VRCONS(I'//S6(1:SLEN6)//')))*'
        J=2
        WRITE(XFIL,'(5X,I1,A)') J,'      ('//S2(1:SLEN2)//'-VRCONS(I'//
     1              S6(1:SLEN6)//'))+VRCONS(I'//S6(1:SLEN6)//'+1)'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S7,SLEN7)
        WRITE(XFIL,'(6X,A)') 'IR'//S7(1:SLEN7)//'=(VRCONS(I'//
     1         S6(1:SLEN6)//'+3)-VRCONS(I'//S6(1:SLEN6)//'+1))/'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'      (VRCONS(I'//S6(1:SLEN6)//
     1                          '+2)-VRCONS(I'//S6(1:SLEN6)//'))'
        ENDIF
C
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE' 
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        QC=QC+2
        GO TO 100
C
      ELSE IF (X .EQ. SPLINE) THEN
        CALL MAKEX(PVVA,OFS,VQD(QC,5),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC+1,5),S2,SLEN2)
        CALL ITOA(VQD(QC,3),S3,SLEN3)
        CALL ITOA(VQD(QC,2),S4,SLEN4)
C
        CALL ITOA(VQD(QC,2)+1,S5,SLEN5)
        CALL ITOA(VQD(QC,2)+2,S6,SLEN6)
        CALL ITOA(VQD(QC,2)+3,S7,SLEN7)
        CALL ITOA(VQD(QC,2)+4,S8,SLEN8)
C
        WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)//'.EQ.20) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//S5(1:SLEN5)//')'//
     1  '+VRCONS('//S6(1:SLEN6)//')*'//S2(1:SLEN2)//'+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS('//S7(1:SLEN7)//')*'//
     1    S2(1:SLEN2)//'**2+VRCONS('//S8(1:SLEN8)//')*'//
     2    S2(1:SLEN2)//'**3'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S10,SLEN10)
        WRITE(XFIL,'(6X,A)') 'IR'//S10(1:SLEN10)//'=VRCONS('//
     1   S6(1:SLEN6)//')+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  2*VRCONS('//S7(1:SLEN7)//')*'//
     1    S2(1:SLEN2)//'+3*VRCONS('//S8(1:SLEN8)//')*'//
     2    S2(1:SLEN2)//'**2'
        ENDIF
C
        CALL ITOA(VQD(QC,2)+15,S9,SLEN9)
        WRITE(XFIL,'(6X,A)') 'ELSE IF ('//S2(1:SLEN2)//
     1            '.LT.VRCONS('//S9(1:SLEN9)//')) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//S5(1:SLEN5)//')'//
     1  '+VRCONS('//S6(1:SLEN6)//')*'//S2(1:SLEN2)//'+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS('//S7(1:SLEN7)//')*'//
     1    S2(1:SLEN2)//'**2+VRCONS('//S8(1:SLEN8)//')*'//
     2    S2(1:SLEN2)//'**3'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S10,SLEN10)
        WRITE(XFIL,'(6X,A)') 'IR'//S10(1:SLEN10)//'=VRCONS('//
     1   S6(1:SLEN6)//')+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  2*VRCONS('//S7(1:SLEN7)//')*'//
     1    S2(1:SLEN2)//'+3*VRCONS('//S8(1:SLEN8)//')*'//
     2    S2(1:SLEN2)//'**2'
        ENDIF
C
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-5,S9,SLEN9)
        WRITE(XFIL,'(6X,A)') 'ELSE IF ('//S2(1:SLEN2)//
     1            '.GE.VRCONS('//S9(1:SLEN9)//')) THEN'
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-9,S5,SLEN5)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-8,S6,SLEN6)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-7,S7,SLEN7)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-6,S8,SLEN8)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-10,S9,SLEN9)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS('//S5(1:SLEN5)//')+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS('//S6(1:SLEN6)//')*('//
     1    S2(1:SLEN2)//'-VRCONS('//S9(1:SLEN9)//'))+'
        J=2
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS('//S7(1:SLEN7)//')*('//
     1    S2(1:SLEN2)//'-VRCONS('//S9(1:SLEN9)//'))**2+'
        J=3
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS('//S8(1:SLEN8)//')*('//
     1    S2(1:SLEN2)//'-VRCONS('//S9(1:SLEN9)//'))**3'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S10,SLEN10)
        WRITE(XFIL,'(6X,A)') 'IR'//S10(1:SLEN10)//'=VRCONS('//
     1    S6(1:SLEN6)//')+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  2*VRCONS('//S7(1:SLEN7)//')*('//
     1    S2(1:SLEN2)//'-VRCONS('//S9(1:SLEN9)//'))+'
        WRITE(XFIL,'(5X,I1,A)') J,'  3*VRCONS('//S8(1:SLEN8)//')*('//
     1    S2(1:SLEN2)//'-VRCONS('//S9(1:SLEN9)//'))**2'
        ENDIF
C
        WRITE(XFIL,'(6X,A)') 'ELSE'
        CALL ITOA(QC,S5,SLEN5)
        CALL ITOA(VQD(QC,4),S6,SLEN6)
        CALL ITOA(VQD(QC,2)+15,S3,SLEN3)
        CALL ITOA(VQD(QC,2)+VQD(QC,3)-10,S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'DO '//S5(1:SLEN5)//' I'//S6(1:SLEN6)//
     1                     '='//S3(1:SLEN3)//','//S4(1:SLEN4)//',5'
        WRITE(XFIL,'(6X,A)') 'IF (('//S2(1:SLEN2)//'.GE.VRCONS(I'//
     1              S6(1:SLEN6)//')) .AND.'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'    ('//S2(1:SLEN2)//'.LT.VRCONS(I'//
     1              S6(1:SLEN6)//'+5))) THEN'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=VRCONS(I'//S6(1:SLEN6)//
     1    '+1)+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS(I'//S6(1:SLEN6)//'+2)*'//
     1     '('//S2(1:SLEN2)//'-VRCONS(I'//S6(1:SLEN6)//'))+'
        J=2
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS(I'//S6(1:SLEN6)//'+3)*'//
     1     '('//S2(1:SLEN2)//'-VRCONS(I'//S6(1:SLEN6)//'))**2+'
        J=3
        WRITE(XFIL,'(5X,I1,A)') J,'  VRCONS(I'//S6(1:SLEN6)//'+4)*'//
     1     '('//S2(1:SLEN2)//'-VRCONS(I'//S6(1:SLEN6)//'))**3'
C
        IF (MODE .EQ. GRAMOD) THEN
        CALL ITOA(VQD(QC,4),S10,SLEN10)
        WRITE(XFIL,'(6X,A)') 'IR'//S10(1:SLEN10)//'=VRCONS(I'//
     1     S6(1:SLEN6)//'+2)+'
        J=1
        WRITE(XFIL,'(5X,I1,A)') J,'  2*VRCONS(I'//S6(1:SLEN6)//'+3)*'//
     1     '('//S2(1:SLEN2)//'-VRCONS(I'//S6(1:SLEN6)//'))+'
        J=2
        WRITE(XFIL,'(5X,I1,A)') J,'  3*VRCONS(I'//S6(1:SLEN6)//'+4)*'//
     1     '('//S2(1:SLEN2)//'-VRCONS(I'//S6(1:SLEN6)//'))**2'
        ENDIF
C
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE' 
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        QC=QC+2
        GO TO 100
C
      ELSE IF (X .EQ. ASSIGN) THEN
        CALL MAKEF(PIFN,PXFN,PVFN,VQD(QC,2),VQD(QC,4),S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S2(1:SLEN2)
        IF (MODE .EQ. GRAMOD) THEN
C
          CALL ITOA(QC,S1,SLEN1)
          CALL ITOA(PVVA,S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,VQD(QC,2),1,VQD(QC,4),
     1               S3,SLEN3)
          S3=S3(1:SLEN3-3)//',I)'
          CALL ITOA(PVVA+1,S4,SLEN4)
          CALL ITOA(VQD(QC,3)-1,S5,SLEN5)
          WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//S2(1:SLEN2)
          WRITE(XFIL,'(6X,A)') 'DFHELP(I)='//S3(1:SLEN3) 
          WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'=0.0D0' 
          WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
C
          IF (PVVA+1 .LT. VQD(QC,3)) THEN
            CALL ITOA(VQD(QC+1,2),S1,SLEN1)
            WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I='
     1            //S4(1:SLEN4)//','//S5(1:SLEN5)
            WRITE(XFIL,'(6X,A)') 'YAUX(I)=0.0D0' 
            WRITE(XFIL,'(I5,1X,A)') VQD(QC+1,2),'CONTINUE'
          ENDIF
C
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,VQD(QC,2),VQD(QC,3),
     1              VQD(QC,4),S1,SLEN1)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'=1.0D0'
C
          AQC=QC
C-----------------------------------------------------------------
          FC=VQD(QC,2)
          FCOFS=VQD(QC,4)
C-----------------------------------------------------------------
          QC=QC-1
 200      IF (QC .GT. 2) THEN 
            IF (VQD(QC,5) .GT. 0) THEN
              CALL FORDF (XFIL,PVVA,PIFN,PXFN,PVFN,MAXVQD,VQD,
     1                    QC,FC,FCOFS,OFS,IERR)
              IF (IERR .NE. 0) RETURN
              QC=QC-1
              GO TO 200
            ENDIF  
          ENDIF  
          QC=AQC+1
C -------------------------
        ELSE
C -------------------------
C   no CONTINUE in XFUN 
C -------------------------
          QC=QC+1
C -------------------------
        ENDIF
        QC=QC+1
        GO TO 100
C
      ELSE IF (X .EQ. LABEL) THEN
        WRITE(XFIL,'(I5,1X,A)') VQD(QC,2),'CONTINUE'
        QC=QC+1
        GO TO 100
C
      ELSE
        IERR=26
        WRITE(*,*) ' FORCDE (5296) : unknown opcode ',X
      ENDIF
      END

      SUBROUTINE FORDF (XFIL,PVVA,PIFN,PXFN,PVFN,MAXVQD,VQD,
     1                  QC,FC,FCOFS,OFS,IERR)
      INTEGER XFIL
      INTEGER PVVA,PIFN,PXFN,PVFN
      INTEGER MAXVQD
      INTEGER VQD(MAXVQD,5)  
      INTEGER QC,FC,FCOFS,OFS
      INTEGER IERR
C
      INTEGER ADD,SUB,MULT,DIV,POWER,LEFT,RIGHT,COMMA,ASSIGN,NLINE
      INTEGER RANGE,RELOP,AND,OR,NOT,INUM,RNUM,ID,SUM,PROD,IN,IF,THEN
      INTEGER ELSE,ENDIF,STDRD,EXTERN,INTERP,PARAM,INDEX,REAL,INT,TABLE
      INTEGER CONINT,LININT,SPLINE,VAR,INFUNC,FUNC,END,GOTO,MARKE,CONTIN
      INTEGER UMINUS,INDVAR,ENDSUM,ENDPRD,BEQ,BRA,LABEL,VECTOR,ACTIVE
      PARAMETER (ADD=43,SUB=45,MULT=42,DIV=47,POWER=94,LEFT=40,RIGHT=41)
      PARAMETER (COMMA=44,ASSIGN=61,NLINE=10,RANGE=257,RELOP=258)
      PARAMETER (AND=259,OR=260,NOT=261,INUM=262,RNUM=263,ID=264)
      PARAMETER (SUM=265,PROD=266,IN=267,IF=268,THEN=269,ELSE=270)
      PARAMETER (ENDIF=271,STDRD=272,EXTERN=273,INTERP=274,PARAM=275)
      PARAMETER (INDEX=276,REAL=277,INT=278,TABLE=279,CONINT=280)
      PARAMETER (LININT=281,SPLINE=282,VAR=283,INFUNC=284,FUNC=285)
      PARAMETER (END=286,GOTO=287,MARKE=288,CONTIN=289,UMINUS=290)
      PARAMETER (INDVAR=291,ENDSUM=292,ENDPRD=293,BEQ=294,BRA=295)
      PARAMETER (LABEL=296,VECTOR=297,ACTIVE=298)
C
      INTEGER MAXSTD,MAXEXT
      PARAMETER (MAXSTD=17,MAXEXT=1)
C
      INTEGER EXTTYP(MAXEXT)
C
      CHARACTER*30 S1,S2,S3,S4,S5,S6,S7,S8
      INTEGER SLEN1,SLEN2,SLEN3,SLEN4,SLEN5,SLEN6,SLEN7,SLEN8
C
      INTEGER EXT,X,Y
C
      DATA EXTTYP /0/
C
      X=VQD(QC,1)
      IF (X .EQ. VAR) THEN
        IF (FC .LT. 0) THEN
          CALL MAKEP(FCOFS,-FC,S1,SLEN1)
        ELSE
          CALL ITOA(FC,S1,SLEN1)
        ENDIF
        CALL ITOA(VQD(QC,4),S2,SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,FC,VQD(QC,5),FCOFS,S3,SLEN3)
        IF (OFS .GT. 0) CALL STRCT(S2,SLEN2,'+IOFS',5)
        IF (FC .LE. PVFN-(PIFN-PXFN)) THEN
          WRITE(XFIL,'(6X,A)') 'DF('//S1(1:SLEN1)//',IAUX('//
     1              S2(1:SLEN2)//'))=DF('//S1(1:SLEN1)//',IAUX('//
     2              S2(1:SLEN2)//'))+'//S3(1:SLEN3)
        ELSE
          WRITE(XFIL,'(6X,A)') 'VGRAD('//S1(1:SLEN1)//',IAUX('//
     1              S2(1:SLEN2)//'))=VGRAD('//S1(1:SLEN1)//',IAUX('//
     1              S2(1:SLEN2)//'))+'//S3(1:SLEN3)
        ENDIF
        RETURN
C 
      ELSE IF (X .EQ. ADD) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='
     1            //S1(1:SLEN1)//'+'//S2(1:SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,4),FCOFS,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='
     1            //S1(1:SLEN1)//'+'//S2(1:SLEN2)
        RETURN  
C
      ELSE IF (X .EQ. SUB) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='
     1            //S1(1:SLEN1)//'+'//S2(1:SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,4),FCOFS,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='
     1            //S1(1:SLEN1)//'-'//S2(1:SLEN2)
        RETURN  
C
      ELSE IF (X .EQ. MULT) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S2,SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S2(1:SLEN2)//'*'//S3(1:SLEN3)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,4),FCOFS,S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S2(1:SLEN2)//'*'//S3(1:SLEN3)
        RETURN  
C
      ELSE IF (X .EQ. DIV) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)//' .EQ. 0.0D0) THEN'
        WRITE(XFIL,'(6X,A)') 'IERR=9'
        WRITE(XFIL,'(6X,A)') 'RETURN'
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S2(1:SLEN2)//'/'//S3(1:SLEN3)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,4),FCOFS,S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S3,SLEN3)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S4,SLEN4)
        CALL FOLD(XFIL,S1(1:SLEN1)//'='//S1(1:SLEN1)//'-'
     1            //S2(1:SLEN2)//'*'//S3(1:SLEN3)
     2            //'/'//S4(1:SLEN4)//'**2',
     3            2*SLEN1+SLEN2+SLEN3+SLEN4+7,IERR)
        RETURN
C
      ELSE IF (X .EQ. POWER) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,4),S3,SLEN3)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'IF ('//S2(1:SLEN2)
     1                         //' .LE. 0.0D0) THEN'
        WRITE(XFIL,'(6X,A)') 'IERR=9'
        WRITE(XFIL,'(6X,A)') 'RETURN'
        WRITE(XFIL,'(6X,A)') 'ENDIF'
        CALL FOLD(XFIL,S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'//S3(1:SLEN3)
     1         //'*'//S2(1:SLEN2)//'**('//S3(1:SLEN3)//'-1.0D0)*'
     2          //S4(1:SLEN4),2*SLEN1+SLEN2+2*SLEN3+SLEN4+14,IERR) 
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,4),FCOFS,S1,SLEN1)
        CALL FOLD(XFIL,S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'//S2(1:SLEN2)
     1          //'**'//S3(1:SLEN3)//'*'//'DLOG('//S2(1:SLEN2)
     2          //')*'//S4(1:SLEN4),2*SLEN1+2*SLEN2+SLEN3+SLEN4+12,IERR)
        RETURN
C
      ELSE IF (X .EQ. POWER+128) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
        CALL ITOA(VQD(QC,2),S2,SLEN2)
        CALL MAKEX(PVVA,OFS,VQD(QC,3),S3,SLEN3)
        CALL ITOA(VQD(QC,2)-1,S4,SLEN4)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S5,SLEN5)
        IF (VQD(QC,2) .GT. 1) THEN
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S2(1:SLEN2)//'*'//S3(1:SLEN3)//'**'
     2            //S4(1:SLEN4)//'*'//S5(1:SLEN5)
        ELSE IF (VQD(QC,2) .EQ. 1) THEN
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1              //S2(1:SLEN2)//'*'//S5(1:SLEN5)
        ELSE
          WRITE(XFIL,'(6X,A)') 'IF ('//S3(1:SLEN3)
     1                         //' .EQ. 0.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR=9'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1               //'('//S2(1:SLEN2)//')'//'*'//S3(1:SLEN3)//'**'
     2               //'('//S4(1:SLEN4)//')'//'*'//S5(1:SLEN5)
        ENDIF
        RETURN
C
      ELSE IF (X .EQ. UMINUS) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S2,SLEN2)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'-'
     1            //S2(1:SLEN2)
        RETURN
C
      ELSE IF ((X .EQ. RNUM) .OR. (X .EQ. REAL)) THEN
        RETURN
C
      ELSE IF ((X .EQ. INUM) .OR. (X .EQ. INT)) THEN
        RETURN
C
      ELSE IF (X .EQ. INDVAR) THEN
        RETURN
C
      ELSE IF (X .EQ. INUM+128) THEN
        RETURN
C
      ELSE IF (X .EQ. ADD+128) THEN
        RETURN
C
      ELSE IF (X .EQ. SUB+128) THEN
        RETURN
C
      ELSE IF (X .EQ. MULT+128) THEN
        RETURN
C
      ELSE IF (X .EQ. DIV+128) THEN
        RETURN
C
      ELSE IF (X .EQ. UMINUS+128) THEN
        RETURN
C
      ELSE IF (X .EQ. FUNC+128) THEN
        RETURN
C
      ELSE IF (X .EQ. INDVAR+128) THEN
        RETURN
C
      ELSE IF (X .EQ. FUNC) THEN
        CALL ITOA(QC,S1,SLEN1)
        CALL ITOA(PVVA,S2,SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,FC,1,FCOFS,S3,SLEN3)
        S3=S3(1:SLEN3-3)//',I)'
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,FC,VQD(QC,5),FCOFS,S5,SLEN5)
        IF (VQD(QC,3) .GT. 0) THEN
          CALL MAKEI(OFS,VQD(QC,4),S4,SLEN4)
C
          IF (FC .GT. 0) THEN
            CALL ITOA(FC,S6,SLEN6)
          ELSE
            CALL MAKEP(FCOFS,-FC,S6,SLEN6)
          ENDIF
          S7='DF('//S4(1:SLEN4)//',I)'
          SLEN7=3+SLEN4+3
          IF (S3 .EQ. S7) THEN
            WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//
     1                            S2(1:SLEN2)
            WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+'
     1                //'DFHELP(I)*'//S5(1:SLEN5)
            WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
          ELSE
            IF ((FCOFS .EQ. VQD(QC,2)) .AND. (VQD(QC,4) .NE. 0)) THEN
              WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//
     1                              S2(1:SLEN2)
              WRITE(XFIL,'(6X,A)') 'IF ('//S6(1:SLEN6)//' .NE. '//
     1                              S4(1:SLEN4)//') THEN'
              WRITE(XFIL,'(6X,A)') 'DFHELP(I)=DF('//S4(1:SLEN4)//',I)'
              WRITE(XFIL,'(6X,A)') 'ENDIF'
              WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+'
     1                  //'DFHELP(I)*'//S5(1:SLEN5)
              WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
            ELSE
              WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//
     1                              S2(1:SLEN2)
              WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+'
     1                  //'DF('//S4(1:SLEN4)//',I)*'//S5(1:SLEN5)
              WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
            ENDIF
          ENDIF
C
          RETURN
        ELSE IF (VQD(QC,3) .EQ. 0) THEN
          CALL ITOA(VQD(QC,2),S4,SLEN4)
C
          IF (FC .GT. 0) THEN
            CALL ITOA(FC,S6,SLEN6)
          ELSE
            CALL MAKEP(FCOFS,-FC,S6,SLEN6)
          ENDIF
C
          S7='DF('//S4(1:SLEN4)//',I)'
          SLEN7=3+SLEN4+3
          IF (S3 .EQ. S7) THEN
            WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//
     1                            S2(1:SLEN2)
            WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+'
     1                //'DFHELP(I)*'//S5(1:SLEN5)
            WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
          ELSE
            WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//
     1                            S2(1:SLEN2)
            WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+'
     1                //'DF('//S4(1:SLEN4)//',I)*'//S5(1:SLEN5)
            WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
          ENDIF
          RETURN 
        ELSE IF (VQD(QC,3) .EQ. -1) THEN
          CALL ITOA(VQD(QC,2),S4,SLEN4)
C
          IF (FC .GT. 0) THEN
            CALL ITOA(FC,S6,SLEN6)
          ELSE
            CALL MAKEP(FCOFS,-FC,S6,SLEN6)
          ENDIF
          S7='VGRAD('//S4(1:SLEN4)//',I)'
          SLEN7=6+SLEN4+3
          IF (S3 .EQ. S7) THEN
            WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//
     1                           S2(1:SLEN2)
            WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+'
     1                //'DFHELP(I)*'//S5(1:SLEN5)
            WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
          ELSE
            WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//
     1                           S2(1:SLEN2)
            WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+'
     1                //'VGRAD('//S4(1:SLEN4)//',I)*'//S5(1:SLEN5)
            WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
          ENDIF
          RETURN
        ELSE
          RETURN
        ENDIF
C-----------------------------------------------------------------
      ELSE IF ((X .EQ. SUM+128) .OR. (X .EQ. PROD+128)) THEN
        RETURN
C
      ELSE IF ((X .EQ. SUM) .OR. (X .EQ. PROD)) THEN
        OFS=OFS-1
        WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE' 
        RETURN
C
      ELSE IF ((X .EQ. ENDSUM) .OR. (X .EQ. ENDPRD)) THEN
        OFS=OFS+1
        CALL ITOA(VQD(QC,2)*VQD(VQD(QC,5),3),S1,SLEN1)
        CALL ITOA(VQD(VQD(QC,5),2)-1,S2,SLEN2)
        CALL ITOA(VQD(VQD(QC,5),3)-1,S3,SLEN3)
        CALL ITOA(VQD(VQD(QC,5),4),S4,SLEN4)
        CALL ITOA(VQD(QC,5),S5,SLEN5)
        CALL ITOA(VQD(QC,2),S6,SLEN6)
        WRITE(XFIL,'(6X,A)') 'OFS=OFS+'//S1(1:SLEN1)
C
        CALL ITOA(VQD(QC,3)*VQD(VQD(QC,5),3),S7,SLEN7)
        CALL ITOA(VQD(QC,3),S8,SLEN8)
        WRITE(XFIL,'(6X,A)') 'IOFS=IOFS+'//S7(1:SLEN7)
C
        WRITE(XFIL,'(6X,A)') 'DO '//S5(1:SLEN5)//' I'//S2(1:SLEN2)//
     1                     '='//S3(1:SLEN3)//',0,-1'
        WRITE(XFIL,'(6X,A)') 'IX'//S2(1:SLEN2)//'=VINDEX('//
     1                     S4(1:SLEN4)//'+I'//S2(1:SLEN2)//')'
        WRITE(XFIL,'(6X,A)') 'OFS=OFS-'//S6(1:SLEN6)
C
        WRITE(XFIL,'(6X,A)') 'IOFS=IOFS-'//S8(1:SLEN8)
C
        RETURN
C-----------------------------------------------------------------
      ELSE IF (X .EQ. STDRD) THEN
        Y=VQD(QC,2)
        IF (Y .EQ. 1) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //'DSIGN(1.0D0,'//S2(1:SLEN2)//')*'//S3(1:SLEN3) 
          RETURN
        ELSE IF (Y .EQ. 2) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'IF ('//S2(1:SLEN2)
     1                         //' .LT. 0.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR=53'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ELSE IF ('//S2(1:SLEN2)
     1                         //' .EQ. 0.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR=9'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //'0.5D0/DSQRT('//S2(1:SLEN2)//')*'//S3(1:SLEN3) 
          RETURN
        ELSE IF (Y .EQ. 3) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //'DEXP('//S2(1:SLEN2)//')*'//S3(1:SLEN3) 
          RETURN
        ELSE IF (Y .EQ. 4) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'IF ('//S2(1:SLEN2)
     1                         //' .LE. 0.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR=52'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/'//S2(1:SLEN2) 
          RETURN
        ELSE IF (Y .EQ. 5) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'IF ('//S2(1:SLEN2)
     1                         //' .LE. 0.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR=52'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          CALL FOLD(XFIL,S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'//S3(1:SLEN3)
     1          //'/('//S2(1:SLEN2)//'*DLOG(10.0D0))',2*SLEN1+SLEN2
     2          +SLEN3+18,IERR)
          RETURN
        ELSE IF (Y .EQ. 6) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //'DCOS('//S2(1:SLEN2)//')*'//S3(1:SLEN3) 
          RETURN
        ELSE IF (Y .EQ. 7) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'-'
     1            //'DSIN('//S2(1:SLEN2)//')*'//S3(1:SLEN3) 
          RETURN
        ELSE IF (Y .EQ. 8) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/DCOS('//S2(1:SLEN2)//')**2'
          RETURN
        ELSE IF (Y .EQ. 9) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'IF (DABS('//S2(1:SLEN2)//
     1                            ') .GE. 1.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR = 54'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/DSQRT(1.0D0-'//S2(1:SLEN2)//'**2)'
          RETURN
        ELSE IF (Y .EQ. 10) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'IF (DABS('//S2(1:SLEN2)//
     1                            ') .GE. 1.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR = 55'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'-'
     1            //S3(1:SLEN3)//'/DSQRT(1.0D0-'//S2(1:SLEN2)//'**2)'
          RETURN
        ELSE IF (Y .EQ. 11) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/(1.0D0+'//S2(1:SLEN2)//'**2)'
          RETURN
        ELSE IF (Y .EQ. 12) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //'DCOSH('//S2(1:SLEN2)//')*'//S3(1:SLEN3)
          RETURN
        ELSE IF (Y .EQ. 13) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //'DSINH('//S2(1:SLEN2)//')*'//S3(1:SLEN3)
          RETURN
        ELSE IF (Y .EQ. 14) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/DCOSH('//S2(1:SLEN2)//')**2'
          RETURN
        ELSE IF (Y .EQ. 15) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/DSQRT(1.0D0+'//S2(1:SLEN2)//'**2)'
          RETURN
        ELSE IF (Y .EQ. 16) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'IF ('//S2(1:SLEN2)//
     1                            ' .LE. 1.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR = 56'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/DSQRT('//S2(1:SLEN2)//'**2-1.0D0)'
          RETURN
        ELSE IF (Y .EQ. 17) THEN
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,3),FCOFS,S1,SLEN1)
          CALL MAKEX(PVVA,OFS,VQD(QC,3),S2,SLEN2)
          CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1           FC,VQD(QC,5),FCOFS,S3,SLEN3)
          WRITE(XFIL,'(6X,A)') 'IF (DABS('//S2(1:SLEN2)//
     1                            ') .GE. 1.0D0) THEN'
          WRITE(XFIL,'(6X,A)') 'IERR = 51'
          WRITE(XFIL,'(6X,A)') 'RETURN'
          WRITE(XFIL,'(6X,A)') 'ENDIF'
          WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1            //S3(1:SLEN3)//'/(1.0D0-'//S2(1:SLEN2)//'**2)'
          RETURN
        ELSE
          IERR=26
          WRITE(*,*) ' FORDF (5504) : unknown standard function ',Y
          RETURN
        ENDIF
C
      ELSE IF (X .EQ. EXTERN) THEN
        EXT=VQD(QC,2)
        IF (EXTTYP(EXT) .EQ. 1) THEN
          CALL MAKEI(OFS,VQD(QC,3),S1,SLEN1)
          WRITE(XFIL,'(6X,A)') 'EXTPAR(1)='//S1(1:SLEN1)
        ELSE IF (EXTTYP(EXT) .EQ. 2) THEN
          CALL MAKEI(OFS,VQD(QC,3),S1,SLEN1)
          CALL MAKEI(OFS,VQD(QC,4),S2,SLEN2)
          WRITE(XFIL,'(6X,A)') 'EXTPAR(1)='//S1(1:SLEN1)
          WRITE(XFIL,'(6X,A)') 'EXTPAR(2)='//S2(1:SLEN2)
        ENDIF
        CALL ITOA(EXT,S1,SLEN1)
        WRITE(XFIL,'(6X,A)') 'CALL EXTGRA('//S1(1:SLEN1)//
     1         ',X,N,Z,EXTPAR)'
C
        CALL ITOA(QC,S1,SLEN1)
        CALL ITOA(PVVA,S2,SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1         FC,1,FCOFS,S3,SLEN3)
        S3=S3(1:SLEN3-3)//',I)'
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1         FC,VQD(QC,5),FCOFS,S4,SLEN4)
        WRITE(XFIL,'(6X,A)') 'DO '//S1(1:SLEN1)//' I=1,'//S2(1:SLEN2)
        WRITE(XFIL,'(6X,A)') S3(1:SLEN3)//'='//S3(1:SLEN3)//'+Z(I)*'
     1          //S4(1:SLEN4)
        WRITE(XFIL,'(I5,1X,A)') QC,'CONTINUE'
        RETURN    
C
      ELSE IF (X .EQ. CONINT) THEN
          RETURN
C
      ELSE IF (X .EQ. CONINT+128) THEN
          RETURN
C
      ELSE IF (X .EQ. LININT) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1         FC,VQD(QC+1,5),FCOFS,S1,SLEN1)
        CALL ITOA(VQD(QC,4),S2,SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1         FC,VQD(QC,5),FCOFS,S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1          //'IR'//S2(1:SLEN2)//'*'//S3(1:SLEN3) 
        RETURN
C
      ELSE IF (X .EQ. LININT+128) THEN
          RETURN
C
      ELSE IF (X .EQ. SPLINE) THEN
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1         FC,VQD(QC+1,5),FCOFS,S1,SLEN1)
        CALL ITOA(VQD(QC,4),S2,SLEN2)
        CALL MAKEY(PIFN,PXFN,PVFN,PVVA,OFS,
     1         FC,VQD(QC,5),FCOFS,S3,SLEN3)
        WRITE(XFIL,'(6X,A)') S1(1:SLEN1)//'='//S1(1:SLEN1)//'+'
     1          //'IR'//S2(1:SLEN2)//'*'//S3(1:SLEN3) 
        RETURN
C
      ELSE IF (X .EQ. SPLINE+128) THEN
          RETURN
C
      ELSE IF ((X .EQ. VECTOR) .OR. (X .EQ. ACTIVE)) THEN
          RETURN
C
      ENDIF
      IERR=26
      WRITE(*,*) ' FORDF (5535) : unknown opcode ',X
      RETURN
      END
C
C
C
      SUBROUTINE MAKEF (PIFN,PXFN,PVFN,I,J,STRING,STRLEN)
      INTEGER PIFN,PXFN,PVFN,I,J
      CHARACTER*(*) STRING
      INTEGER STRLEN
C
      CHARACTER*10 S
      INTEGER SLEN
C
      IF (I .GE. 0) THEN
        CALL ITOA(I,S,SLEN)
      ELSE
        CALL MAKEP(J,-I,S,SLEN)
      ENDIF
      IF (I .LE. PVFN-(PIFN-PXFN)) THEN
        STRING='F('//S(1:SLEN)//')'
        STRLEN=2+SLEN+1
      ELSE
        STRING='VFUNC('//S(1:SLEN)//')'
        STRLEN=6+SLEN+1
      ENDIF
      RETURN
      END
C
C
C 
      SUBROUTINE MAKEL (I,STRING,STRLEN)
      INTEGER I
      CHARACTER*(*) STRING
      INTEGER STRLEN
C
      CHARACTER*10 S
      INTEGER SLEN
C
      CALL ITOA(I,S,SLEN)
      STRING='LAUX('//S(1:SLEN)//')'
      STRLEN=5+SLEN+1
      RETURN
      END
C
C
C 
      SUBROUTINE MAKEI (OFS,I,STRING,STRLEN)
      INTEGER OFS,I
      CHARACTER*(*) STRING
      INTEGER STRLEN
C
      CHARACTER*20 S
      INTEGER SLEN
C
      IF (I .GE. 0) THEN
        CALL ITOA(I,S,SLEN)
        IF (OFS .GT. 0) THEN
          S=S(1:SLEN)//'+IOFS'
          SLEN=SLEN+5
        ENDIF
        STRING='IAUX('//S(1:SLEN)//')'
        STRLEN=5+SLEN+1
        RETURN
      ELSE
        RETURN
      ENDIF
      END
C
C
C 
      SUBROUTINE MAKEX (PVVA,OFS,I,STRING,STRLEN)
      INTEGER PVVA,I,OFS
      CHARACTER*(*) STRING
      INTEGER STRLEN
C
      CHARACTER*20 S
      INTEGER SLEN
C
      IF (I .GE. 0) THEN
        CALL ITOA(I,S,SLEN)
      ENDIF
      IF (I .LE. PVVA) THEN
        STRING='X('//S(1:SLEN)//')'
        STRLEN=2+SLEN+1
      ELSE
        IF (OFS .GT. 0) THEN
          S=S(1:SLEN)//'+OFS'
          SLEN=SLEN+4
        ENDIF
        STRING='XAUX('//S(1:SLEN)//')'
        STRLEN=5+SLEN+1
      ENDIF
      RETURN
      END
C
C
C 
      SUBROUTINE MAKEY (PIFN,PXFN,PVFN,PVVA,OFS,I,J,K,STRING,STRLEN)
      INTEGER PIFN,PXFN,PVFN,PVVA,OFS,I,J
      CHARACTER*(*) STRING
      INTEGER STRLEN
C
      CHARACTER*20 S1,S2
      INTEGER SLEN1,SLEN2
      INTEGER K
C
      IF (I .GE. 0) THEN
        CALL ITOA(I,S1,SLEN1)
      ELSE
        CALL MAKEP(K,-I,S1,SLEN1)
      ENDIF
      IF (J .GE. 0) THEN
        CALL ITOA(J,S2,SLEN2)
      ENDIF
      IF (J .GT. PVVA) THEN
        IF (OFS .GT. 0) THEN
          S2=S2(1:SLEN2)//'+OFS'
          SLEN2=SLEN2+4
        ENDIF
        STRING='YAUX('//S2(1:SLEN2)//')'
        STRLEN=5+SLEN2+1
      ELSE IF (I .LE. PVFN-(PIFN-PXFN)) THEN
        STRING='DF('//S1(1:SLEN1)//','//S2(1:SLEN2)//')'
        STRLEN=3+SLEN1+1+SLEN2+1
      ELSE
        STRING='VGRAD('//S1(1:SLEN1)//','//S2(1:SLEN2)//')'
        STRLEN=6+SLEN1+1+SLEN2+1
      ENDIF
      RETURN
      END
C
C
C     
      SUBROUTINE FTOA (X,S,SLEN)
      DOUBLE PRECISION X
      CHARACTER*24 S
      INTEGER SLEN
C
      INTEGER I
C
      WRITE(S,'(D24.17)',ERR=11) X
      DO 10 I=1,24
        IF (S(I:I) .NE. ' ') THEN
          S=S(I:24)
          SLEN=24-I+1
          RETURN
        ENDIF
 10   CONTINUE
 11   CONTINUE
      SLEN=24
      RETURN
      END
C
C
C     
      SUBROUTINE ITOA (N,S,SLEN)
      INTEGER N
      CHARACTER*10 S
      INTEGER SLEN
C
      INTEGER I
C
      WRITE(S,'(I10)',ERR=11) N
      DO 10 I=1,10 
        IF (S(I:I) .NE. ' ') THEN
          S=S(I:10)
          SLEN=10-I+1
          RETURN
        ENDIF
 10   CONTINUE
 11   CONTINUE
      SLEN=10
      RETURN
      END
C
C
C     
      SUBROUTINE MAKEP (I,J,STRING,STRLEN)
      INTEGER I,J
      CHARACTER*(*) STRING
      INTEGER STRLEN
C
      CHARACTER*10 S1,S2
      INTEGER SLEN1,SLEN2
C
      IF (I .EQ. 0) THEN
        CALL ITOA(J-1,S1,SLEN1)
        STRING='IX'//S1(1:SLEN1)
        STRLEN=2+SLEN1
        RETURN
      ELSE
        CALL ITOA(I,S1,SLEN1)
        CALL ITOA(J-1,S2,SLEN2)
        STRING=S1(1:SLEN1)//'+IX'//S2(1:SLEN2)
        STRLEN=SLEN1+3+SLEN2
        RETURN
      ENDIF 
      END
C
C
C
      SUBROUTINE STRCT (S1,SLEN1,S2,SLEN2)
      CHARACTER*(*) S1,S2
      INTEGER SLEN1,SLEN2
C
      IF (SLEN1 .EQ. 0) THEN
        S1=S2(1:SLEN2)
        SLEN1=SLEN2
      ELSE
        S1=S1(1:SLEN1)//S2(1:SLEN2)
        SLEN1=SLEN1+SLEN2
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE IDATA (XFIL,V,N,ID)
      INTEGER XFIL,N
      INTEGER V(N)
      CHARACTER*6 ID
C
      INTEGER FIRST,LAST,COUNT,ROW,COL
      CHARACTER*80 STRING,S1,S2
      INTEGER STRLEN,SLEN1,SLEN2
C
      CHARACTER CONT(10)
      DATA CONT /'1','2','3','4','5','6','7','8','9','O'/
C
      FIRST=1
 10   IF (N .LT. FIRST+49) THEN
        LAST=N
      ELSE
        LAST=FIRST+49
      ENDIF
      CALL ITOA(FIRST,S1,SLEN1)
      CALL ITOA(LAST,S2,SLEN2)
      WRITE(XFIL,'(6X,A)') 'DATA ('//ID//'(I), I='//S1(1:SLEN1)//
     1                     ','//S2(1:SLEN2)//')'
      COUNT=FIRST
      ROW=1
      COL=1
      STRING='/'
      STRLEN=1
 20   CALL ITOA(V(COUNT),S1,SLEN1)
      CALL STRCT(STRING,STRLEN,S1,SLEN1)
      IF (COUNT .EQ. LAST) THEN
        CALL STRCT(STRING,STRLEN,'/',1)
      ELSE
        CALL STRCT(STRING,STRLEN,',',1)
      ENDIF
      IF ((COUNT .EQ. LAST) .OR. (COL .EQ. 5)) THEN
        WRITE(XFIL,'(5X,A,5X,A)') CONT(ROW),STRING(1:STRLEN)
        STRING=' '
        STRLEN=1
        ROW=ROW+1
        COL=1
      ELSE
        COL=COL+1
      ENDIF
      COUNT=COUNT+1
      IF (COUNT .LE. LAST) GO TO 20
      IF (LAST .LT. N) THEN
        FIRST=LAST+1
        GO TO 10 
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE FDATA (XFIL,V,N,ID)
      INTEGER XFIL,N
      DOUBLE PRECISION V(N)
      CHARACTER*6 ID
C
      INTEGER FIRST,LAST,COUNT,ROW,COL
      CHARACTER*80 STRING,S1,S2
      INTEGER STRLEN,SLEN1,SLEN2
C
      CHARACTER CONT(10)
      DATA CONT /'1','2','3','4','5','6','7','8','9','O'/
C
      FIRST=1
 10   IF (N .LT. FIRST+19) THEN
        LAST=N
      ELSE
        LAST=FIRST+19
      ENDIF
      CALL ITOA(FIRST,S1,SLEN1)
      CALL ITOA(LAST,S2,SLEN2)
      WRITE(XFIL,'(6X,A)') 'DATA ('//ID//'(I), I='//S1(1:SLEN1)//
     1                     ','//S2(1:SLEN2)//')'
      COUNT=FIRST
      ROW=1
      COL=1
      STRING='/'
      STRLEN=1
 20   CALL FTOA(V(COUNT),S1,SLEN1)
      CALL STRCT(STRING,STRLEN,S1,SLEN1)
      IF (COUNT .EQ. LAST) THEN
        CALL STRCT(STRING,STRLEN,'/',1)
      ELSE
        CALL STRCT(STRING,STRLEN,',',1)
      ENDIF
      IF ((COUNT .EQ. LAST) .OR. (COL .EQ. 2)) THEN
        WRITE(XFIL,'(5X,A,5X,A)') CONT(ROW),STRING(1:STRLEN)
        STRING=' '
        STRLEN=1
        ROW=ROW+1
        COL=1
      ELSE 
        COL=COL+1
      ENDIF
      COUNT=COUNT+1
      IF (COUNT .LE. LAST) GO TO 20
      IF (LAST .LT. N) THEN
        FIRST=LAST+1
        GO TO 10 
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE FOLD (XFIL,STRING,STRLEN,IERR)
      CHARACTER*(*) STRING
      INTEGER XFIL,STRLEN,IERR
C
      INTEGER BRACE,BREAK,CLEN,WIDTH
      CHARACTER CONT*4,X,STR*80
      LOGICAL SEP
      INTEGER I
C
      SEP(X)=((X .EQ. '+') .OR. (X .EQ. '-') .OR. (X .EQ. '*') .OR.
     1        (X .EQ. '/') .OR. (X .EQ. '='))
C
      WIDTH=56
      CONT=' '
      CLEN=1
 1    IF (STRLEN .LE. WIDTH) THEN
        STR(1:STRLEN)=STRING(1:STRLEN)
        WRITE(XFIL,'(5X,A)') CONT(1:CLEN)//STR(1:STRLEN)
      ELSE
        BRACE=0
        BREAK=0
        DO 10 I=1,WIDTH
          IF (STRING(I:I) .EQ. '(') BRACE=BRACE+1
          IF (STRING(I:I) .EQ. ')') BRACE=BRACE-1
          IF (SEP(STRING(I:I)) .AND. (BRACE .EQ. 0)) BREAK=I
 10     CONTINUE
        IF (BREAK .EQ. 0) THEN
          IERR=46
          RETURN
        ENDIF
        STR(1:BREAK)=STRING(1:BREAK)
        WRITE(XFIL,'(5X,A)') CONT(1:CLEN)//STR(1:BREAK)
        WIDTH=53
        CONT='/   '
        CLEN=4
        STRING=STRING(BREAK+1:STRLEN)
        STRLEN=STRLEN-BREAK
        IF (STRLEN .GT. 0) GO TO 1
      ENDIF
      RETURN
      END

C*********************************************
C                                            *
C   PROGRAM   : PCOMP                        *
C   MODULE    : PL (LARGE COMPILER)          *
C   ABSTRACT  : FORTRAN PRECOMPILER          *
C   KEY WORD  : AUTOMATIC DIFFERENTIATION    *
C   SOURCE    : PCOMP 2.3 by M.LIEPELT       *
C               PCOMP 3.0 by M.DOBMANN       *
C   COPYRIGHT : C.TRASSL, K.SCHITTKOWSKI     *
C               MATHEMATISCHES INSTITUT,     *
C               UNIVERSITAET BAYREUTH,       *
C               D-95440 BAYREUTH, GERMANY    *
C   DATE      : NOVEMBER 23, 1997            *
C   VERSION   : 5.3                          *
C                                            *
C*********************************************
C
C
C
      SUBROUTINE YYPAR (INPUT,WA,LWA,IWA,LIWA,PLWA,PLIWA,MODE,IERR,LNUM,
     1                  DBGFIL,DBGLVL)
C
      INTEGER INPUT
      INTEGER LWA,LIWA,PLWA,PLIWA,MODE,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
C ------------------------------------------------------
      INTEGER DBGFIL,DBGLVL
      LOGICAL DEBUG
      INTEGER PSLVL
      LOGICAL PSCHK
C ------------------------------------------------------
      INTEGER INFOLI(15)
      INTEGER VEK3(3),VEK4(4),VEK5(5),VEK7(7)
      INTEGER IHELP1,IHELP2,IHELP3,IHELP4,IHELP5,IPOL
C  
      DOUBLE PRECISION HELP3,HELP4,HELP6
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER ADD,SUB,MULT,DIV,POWER,LEFT,RIGHT,COMMA,ASSIGN,NLINE
      INTEGER RANGE,RELOP,AND,OR,NOT,INUM,RNUM,ID,SUM,PROD,IN,IF,THEN
      INTEGER ELSE,ENDIF,STDRD,EXTERN,INTERP,PARAM,INDEX,REAL,INT,TABLE
      INTEGER CONINT,LININT,SPLINE,VAR,INFUNC,FUNC,END,GOTO,MARKE,CONTIN
      INTEGER UMINUS,INDVAR,ENDSUM,ENDPRD,BEQ,BRA,LABEL,VECTOR,ACTIVE
      PARAMETER (ADD=43,SUB=45,MULT=42,DIV=47,POWER=94,LEFT=40,RIGHT=41)
      PARAMETER (COMMA=44,ASSIGN=61,NLINE=10,RANGE=257,RELOP=258)
      PARAMETER (AND=259,OR=260,NOT=261,INUM=262,RNUM=263,ID=264)
      PARAMETER (SUM=265,PROD=266,IN=267,IF=268,THEN=269,ELSE=270)
      PARAMETER (ENDIF=271,STDRD=272,EXTERN=273,INTERP=274,PARAM=275)
      PARAMETER (INDEX=276,REAL=277,INT=278,TABLE=279,CONINT=280)
      PARAMETER (LININT=281,SPLINE=282,VAR=283,INFUNC=284,FUNC=285)
      PARAMETER (END=286,GOTO=287,MARKE=288,CONTIN=289,UMINUS=290)
      PARAMETER (INDVAR=291,ENDSUM=292,ENDPRD=293,BEQ=294,BRA=295)
      PARAMETER (LABEL=296,VECTOR=297,ACTIVE=298)
      INTEGER YYEXCA(0:17),YYACT(0:531),YYPACT(0:387),YYPGO(0:55)
      INTEGER YYR1(0:138),YYR2(0:138),YYCHK(0:387),YYDEF(0:387)
      CHARACTER*10 YYTNAM(0:44)
      INTEGER YYTVAL(0:44)
      CHARACTER*100 YYREDS(0:138)
      INTEGER YYRLEN(0:138)
C
      INTEGER YYMDEP,YYERRC,YYNPRO,YYLAST,YYFLAG,YYACPT
      INTEGER YYABRT
      PARAMETER (YYMDEP=150,YYERRC=256,YYNPRO=139)
      PARAMETER (YYLAST=532,YYFLAG=-3000,YYACPT=0,YYABRT=31)
C
      INTEGER YYV(0:YYMDEP-1),YYS(0:YYMDEP-1)
      INTEGER YYPV,YYPS,YYSTAT,YYTMP,YYNERR,YYERRF,YYCHAR
      INTEGER YYPVT,YYXPV,YYXPS,YYXSTA,YYXN,YYVAL,YYXI,YYXLEN
C
      INTEGER S1,S2,S3,S4,S5,S6,S7,S8,S9,S10
      INTEGER S12,S14
      INTEGER PC,PC1,PC2,DIM
C
      INTEGER IVAL
      DOUBLE PRECISION FVAL
      INTEGER I,J,YYLEX
C
      CHARACTER STAR*5,CONT*1,LINE*66
      INTEGER   LPOS,SPOS
      LOGICAL   EOF
C
      INTEGER MAXSYM
      PARAMETER (MAXSYM=2000)
      INTEGER YYLVAL
      CHARACTER*20 SYMNAM(MAXSYM)
      INTEGER SYMTYP(MAXSYM),SYMREF(MAXSYM),SYMEND
C
      INTEGER OIFLAG,FUFLAG,SPFLAG
      INTEGER MAXSTD,MAXEXT
      PARAMETER (MAXSTD=17,MAXEXT=1)
      INTEGER STDTYP(MAXSTD),EXTTYP(MAXEXT)
C
      INTEGER GRAD,HESS
      PARAMETER (GRAD=1,HESS=2)
C
      INTEGER MAXMAR
      PARAMETER (MAXMAR=100)
      INTEGER MARKST(MAXMAR,2),GOTOST(MAXMAR,2)
      INTEGER MC,GC
C
      INTEGER GETIWA
      DOUBLE PRECISION GETWA
C
      DATA STDTYP /1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
      DATA EXTTYP /0/
C
      DATA (YYEXCA(I),I=0,17)
     1           /-1,1,0,-1,-2,0,-1,109,257,27,
     2           -2,36,-1,111,257,26,-2,37/
      DATA (YYACT(I),I=0,99)
     O           /39,207,210,108,17,18,19,20,21,22,
     1           23,24,25,26,27,5,72,103,366,362,
     2           73,368,385,306,110,377,376,363,356,354,
     3           336,333,332,331,327,304,288,112,286,75,
     4           74,216,113,215,112,134,131,133,130,113,
     5           209,380,379,371,370,116,369,360,359,358,
     6           89,348,347,115,346,345,344,343,116,341,
     7           324,321,318,314,313,309,308,307,279,258,
     8           256,241,229,174,172,237,173,114,175,132,
     9           218,212,204,152,147,146,142,138,106,105/
      DATA (YYACT(I),I=100,199)
     O           /95,92,70,151,68,66,63,58,56,54,
     1           52,50,47,46,45,44,202,199,177,178,
     2           174,172,196,173,144,175,140,136,96,93,
     3           381,329,284,243,249,176,239,167,170,217,
     4           107,104,60,250,312,160,266,174,172,267,
     5           173,176,175,174,172,337,173,174,175,98,
     6           338,101,175,303,214,174,172,334,173,301,
     7           175,299,176,87,297,224,225,226,227,228,
     8           85,261,100,335,373,233,248,235,219,220,
     9           221,222,174,172,86,173,242,175,240,176/
      DATA (YYACT(I),I=200,299)
     O           /238,84,236,83,82,176,232,365,234,176,
     1           252,254,253,310,165,174,172,176,173,166,
     2           175,268,165,163,269,164,192,166,190,193,
     3           99,191,187,262,263,188,317,165,163,271,
     4           164,273,166,275,176,277,185,311,364,186,
     5           350,280,282,283,156,287,247,157,285,168,
     6           213,169,211,161,150,290,109,176,111,148,
     7           294,208,125,124,126,122,123,90,88,367,
     8           352,351,127,128,129,125,124,126,122,123,
     9           293,278,276,274,295,127,128,129,149,255/
      DATA (YYACT(I),I=300,399)
     O           /319,272,320,149,322,305,323,165,163,162,
     1           164,257,166,145,141,137,330,206,94,91,
     2           203,158,302,65,62,231,340,316,174,172,
     3           230,173,184,175,80,298,255,174,172,183,
     4           173,182,175,181,153,97,76,386,250,251,
     5           250,251,250,251,174,172,383,173,34,175,
     6           174,172,300,173,81,175,382,174,172,374,
     7           173,35,175,375,372,378,174,172,296,173,
     8           176,175,355,174,172,292,173,387,175,176,
     9           174,172,291,173,205,175,270,174,172,353/
      DATA (YYACT(I),I=400,499)
     O           /173,189,175,229,174,172,176,173,171,175,
     1           174,172,176,173,342,175,315,165,163,176,
     2           164,339,166,328,326,325,174,172,176,173,
     3           259,175,246,174,172,176,173,245,175,244,
     4           174,172,176,173,201,175,200,198,197,176,
     5           195,264,165,163,265,164,176,166,194,223,
     6           165,163,176,164,159,166,155,154,143,139,
     7           135,79,78,77,49,48,43,42,176,41,
     8           40,28,64,180,179,176,121,120,119,118,
     9           117,384,176,361,357,349,281,102,71,69/
      DATA (YYACT(I),I=500,531)
     O           /67,61,59,57,55,289,260,53,51,16,
     1           38,15,37,14,36,13,12,11,33,10,
     2           32,9,31,8,30,7,29,6,4,3,
     3           2,1/
      DATA (YYPACT(I),I=0,99)
     O           /-3000,-3000,-271,-3000,-3000,471,-3000,-3000,-3000,
     O            -3000,
     1           -3000,-3000,-3000,-3000,-3000,-3000,-3000,470,469,467,
     2           466,-149,-150,-151,-152,465,464,-153,-3000,-154,
     3           -155,-156,-157,-120,61,61,60,-160,-162,-248,
     4           -3000,-3000,-3000,-3000,306,463,462,461,-3000,-3000,
     5           324,-3000,143,-3000,142,-3000,140,-3000,133,-3000,
     6           15,-3000,56,-134,-3000,55,-135,-3000,305,-3000,
     7           -3000,-3000,121,-3000,-272,-121,-165,-3000,-3000,-3000,
     8           -3000,-166,-122,4,23,-216,23,-217,460,-136,
     9           52,459,-137,51,458,-139,50,-170,259,254/
      DATA (YYPACT(I),I=100,199)
     O           /23,-171,304,457,456,213,280,454,-112,219,
     1           265,-3000,-3,-3,398,23,23,-3000,-3000,-3000,
     2           -3000,-3000,-3000,-3000,-3000,-3000,303,301,299,292,
     3           205,191,391,187,185,-3000,448,440,-141,-3000,
     4           438,437,-146,-3000,436,434,-147,279,-3000,-172,
     5           -3000,384,276,10,-3000,-3000,218,-173,216,-3000,
     6           -221,-123,-174,-3,-3,-3,-3,418,-3000,-3000,
     7           -3000,-3000,23,23,23,23,23,362,57,290,
     8           285,-3,23,-3,23,141,-179,139,-126,-3000,
     9           137,-183,135,-129,-3000,-3000,429,-3000,-3000,427/
      DATA (YYPACT(I),I=200,299)
     O           /-3000,-3000,422,212,-3000,-3000,125,93,10,10,
     1           78,-184,270,-185,420,-3000,-3000,-3000,120,172,
     2           172,-3000,-3000,-3000,115,115,57,57,57,-3000,
     3           23,23,410,105,180,355,23,260,23,252,
     4           23,251,23,250,-3000,-3000,-3000,-186,23,-3000,
     5           10,10,-3000,41,91,23,-229,211,-231,-3000,
     6           -3000,-221,348,341,-3000,-3,-3000,23,-3000,-3,
     7           -3000,334,113,325,110,318,108,312,102,-232,
     8           295,-246,-3000,-116,-3000,368,-187,-188,-189,203,
     9           -113,-190,-191,375,286,195,-192,23,-3000,23/
      DATA (YYPACT(I),I=300,387)
     O           /-193,23,-3000,23,-194,-3000,415,414,-233,413,
     1           -3000,-131,-221,-234,-235,-3000,-3000,-3000,-236,123,
     2           173,-237,111,150,411,-3000,-3000,-195,-3000,-3000,
     3           404,-197,-198,-199,-200,-3000,-202,-203,-3000,-3000,
     4           -248,206,-3000,240,239,389,-238,372,-239,-3000,
     5           -205,-3000,-3000,-3000,-206,-3000,-207,-251,-240,204,
     6           163,-253,11,-208,-210,-211,364,144,-3000,363,
     7           -241,-242,-3000,10,-248,-3000,-212,-213,89,356,
     8           346,-3000,-3000,-3000,-247,337,-3000,-248/
      DATA (YYPGO(I),I=0,55)
     O           /0,531,530,529,528,527,526,525,524,523,
     1           522,521,520,519,518,517,358,516,515,514,
     2           513,512,511,510,509,0,508,507,3,506,
     3           505,24,504,2,503,502,501,500,159,499,
     4           498,497,1,496,495,494,493,491,490,489,
     5           488,487,486,484,483,482/
      DATA (YYR1(I),I=0,99)
     O           /0,1,2,2,4,4,4,4,4,4,
     1           4,4,4,4,4,5,6,6,26,7,
     2           8,8,27,29,27,27,28,28,30,30,
     3           31,31,31,31,31,31,31,31,9,10,
     4           10,32,32,32,32,32,11,12,12,34,
     5           34,34,34,34,13,13,14,14,35,35,
     6           35,35,15,17,16,16,36,36,36,36,
     7           20,21,21,37,37,38,38,22,23,23,
     8           39,24,24,25,25,40,40,41,43,44,
     9           40,40,40,47,45,45,46,46,33,33/
      DATA (YYR1(I),I=100,138)
     O           /33,33,33,33,33,33,33,33,33,33,
     1           53,33,54,33,42,42,42,42,42,48,
     2           48,49,49,49,50,50,50,51,51,51,
     3           52,18,19,19,55,55,55,55,3/
      DATA (YYR2(I),I=0,99)
     O           /0,5,4,0,4,4,4,4,4,4,
     1           4,4,4,4,5,4,4,0,9,4,
     2           4,0,13,1,16,21,3,2,7,0,
     3           7,7,7,7,6,5,3,3,5,4,
     4           0,9,23,15,35,19,4,4,0,9,
     5           23,15,35,19,21,33,4,0,7,9,
     6           9,11,7,7,4,0,7,9,9,11,
     7           5,4,0,19,7,7,1,5,4,0,
     8           7,7,21,4,0,9,15,1,1,1,
     9           29,7,7,1,21,0,6,1,7,7/
      DATA (YYR2(I),I=100,138)
     O           /7,7,7,6,5,2,2,2,2,2,
     1           1,19,1,19,7,7,5,6,7,3,
     2           3,3,9,13,3,9,13,3,9,13,
     3           9,7,4,0,7,9,9,11,4/
      DATA (YYCHK(I),I=0,99)
     O           /-3000,-1,-2,-3,-4,286,-5,-7,-9,-11,
     1           -13,-15,-17,-18,-20,-22,-24,275,276,277,
     2           278,279,280,281,282,283,284,285,10,-6,
     3           -8,-10,-12,-14,-16,-16,-19,-21,-23,-25,
     4           10,10,10,10,264,264,264,264,10,10,
     5           264,-26,264,-27,264,-32,264,-34,264,-35,
     6           262,-36,263,45,-55,263,45,-37,264,-39,
     7           264,-40,264,268,288,287,40,10,10,10,
     8           10,40,61,61,61,40,61,40,263,45,
     9           262,263,45,263,263,45,263,40,-38,-38/
      DATA (YYCHK(I),I=100,199)
     O           /61,40,-41,289,262,264,264,262,-28,262,
     1           -31,264,40,45,-33,40,45,-48,-49,-50,
     2           -51,-52,265,266,263,262,264,272,273,274,
     3           264,262,-33,264,262,10,263,263,45,10,
     4           263,263,45,10,263,263,45,264,10,44,
     5           10,-33,264,40,10,10,41,44,41,10,
     6           257,44,44,43,45,42,47,-31,262,264,
     7           -31,10,43,45,42,47,94,-33,-33,-53,
     8           -54,40,40,40,40,41,44,41,44,10,
     9           41,44,41,44,10,10,263,10,10,263/
      DATA (YYCHK(I),I=200,299)
     O           /10,10,263,41,264,10,41,-42,261,40,
     1           -33,44,264,44,-28,264,262,262,264,-31,
     2           -31,-31,-31,41,-33,-33,-33,-33,-33,41,
     3           40,40,-31,-33,-31,-33,61,264,61,262,
     4           61,264,61,262,10,10,10,44,61,41,
     5           259,260,-42,-33,-42,258,264,41,264,10,
     6           -29,61,-33,-33,41,44,41,44,41,44,
     7           41,-33,41,-33,41,-33,41,-33,41,264,
     8           -33,-43,-42,-42,41,-33,267,44,267,-30,
     9           -28,44,44,-31,-33,-31,44,61,10,61/
      DATA (YYCHK(I),I=300,387)
     O           /44,61,10,61,267,10,269,264,264,264,
     1           10,44,257,264,264,41,41,41,264,-33,
     2           -33,264,-33,-33,264,10,10,267,10,262,
     3           -28,267,267,267,44,10,267,44,10,10,
     4           -25,264,10,264,264,264,264,264,264,-44,
     5           44,41,41,10,267,10,267,-45,264,264,
     6           264,-46,270,267,44,44,271,268,10,264,
     7           264,264,10,40,-25,10,267,267,-42,264,
     8           264,41,10,10,-47,269,10,-25/
      DATA (YYDEF(I),I=0,99)
     O           /3,-2,0,1,2,0,17,21,40,48,
     1           57,65,65,133,72,79,84,0,0,0,
     2           0,0,0,0,0,0,0,0,138,4,
     3           5,6,7,8,9,10,11,12,13,14,
     4           15,19,38,46,0,0,0,0,70,77,
     5           0,16,0,20,0,39,0,47,0,56,
     6           0,64,0,0,132,0,0,71,76,78,
     7           76,83,0,87,0,0,0,62,63,131,
     8           81,0,0,0,0,0,0,0,0,0,
     9           0,0,0,0,0,0,0,0,0,0/
      DATA (YYDEF(I),I=100,199)
     O           /0,0,0,0,0,0,0,0,0,-2,
     1           0,-2,0,0,0,0,0,105,106,107,
     2           108,109,110,112,119,120,121,124,127,0,
     3           0,0,0,0,0,58,0,0,0,66,
     4           0,0,0,134,0,0,0,0,74,0,
     5           80,0,0,0,91,92,0,0,0,18,
     6           0,0,0,0,0,0,0,0,36,37,
     7           35,41,0,0,0,0,0,0,104,0,
     8           0,0,0,0,0,0,0,0,0,49,
     9           0,0,0,0,59,60,0,67,68,0/
      DATA (YYDEF(I),I=200,299)
     O           /135,136,0,0,75,85,0,0,0,0,
     1           0,0,0,0,0,26,27,23,0,30,
     2           31,32,33,34,98,99,100,101,102,103,
     3           0,0,0,0,0,0,0,0,0,0,
     4           0,0,0,0,61,69,137,0,0,88,
     5           0,0,116,0,0,0,0,0,0,22,
     6           29,0,0,0,122,0,125,0,128,0,
     7           130,0,0,0,0,0,0,0,0,0,
     8           0,0,114,115,117,118,0,0,0,0,
     9           0,0,0,0,0,0,0,0,43,0/
      DATA (YYDEF(I),I=300,387)
     O           /0,0,51,0,0,86,0,0,0,0,
     1           24,0,0,0,0,123,126,129,0,0,
     2           0,0,0,0,0,84,54,0,82,28,
     3           0,0,0,0,0,45,0,0,53,73,
     4           89,0,25,0,0,0,0,0,0,95,
     5           0,111,113,42,0,50,0,97,0,0,
     6           0,0,0,0,0,0,0,0,84,0,
     7           0,0,90,0,96,55,0,0,0,0,
     8           0,93,44,52,0,0,84,94/
C
      DATA YYV    /YYMDEP*0/
      DATA YYS    /YYMDEP*0/
C
      DATA YYTNAM /
     1  'RANGE','RELOP','AND','OR','NOT','INUM','RNUM','ID',
     2  'SUM','PROD','IN','IF','THEN','ELSE','ENDIF','STANDARD',
     3  'EXTERN','INTERP','PARAM','INDEX','REAL','INT','TABLE',
     4  'CONINT','LININT','SPLINE','VAR','INFUNC','FUNC','END','GOTO',
     5  'LABEL','CONTINUE','+','-','*','/','UMINUS','^','(',')',',','=',
     6  '\n','-unknown-'/
C
      DATA YYTVAL /
     1  257,258,259,260,261,262,263,264,265,266,
     2  267,268,269,270,271,272,273,274,275,276,
     3  277,278,279,280,281,282,283,284,285,286,
     4  287,288,289,43,45,42,47,290,94,40,
     5  41,44,61,10,-1/
C
      DATA (YYRLEN(I),I=0,138) /
     1     19,38,57,32,49,49,47,53,49,69,
     2     69,62,55,51,39,23,57,32,36,23,
     3     57,32,69,40,54,89,20,22,28,22,
     4     32,32,32,32,27,23,15,13,21,54,
     5     31,35,59,48,79,57,23,63,34,38,
     6     62,51,82,60,50,70,57,32,34,38,
     7     39,43,39,39,81,40,42,46,46,50,
     8     24,66,35,54,37,22,20,25,60,33,
     9     35,28,52,18,19,23,34,9,28,44,
     O     76,26,21,52,68,25,26,22,20,20,
     1     20,20,20,19,15,13,17,24,22,29,
     2     10,36,11,37,38,37,27,31,28,13,
     3     13,15,32,45,28,41,50,24,41,54,
     4     44,39,60,33,35,39,39,43,21/
C
      CALL INIRED(YYREDS)
C
      DEBUG=((DBGLVL .EQ. 2) .OR. (DBGLVL .EQ. 3))
      PSLVL=0
      LNUM=0
      LPOS=67
      SPOS=2
      EOF=.FALSE.
      IERR=0
      PC=0
      MC=0
      GC=0
      SYMEND=0
      FUFLAG=0
      SPFLAG=0
C
      PSCHK=.FALSE.
      S10=0
      S14=0
      OIFLAG=0
C
      DO 5 I=1,MAXMAR
        GOTOST(I,1)=0
        MARKST(I,1)=0
        GOTOST(I,2)=0
        MARKST(I,2)=0
 5    CONTINUE   
C
      DO 6 I=1,15
        INFOLI(I)=15
 6    CONTINUE   
      INFOLI(6)=0
      INFOLI(8)=0
      INFOLI(10)=0
      INFOLI(11)=0
      INFOLI(12)=0
      INFOLI(13)=0
C
C ignore leading '\n'!
C
      YYCHAR=YYLEX(INPUT,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF,YYLVAL,
     1             LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,MAXSYM,SYMNAM,
     2             SYMTYP,SYMREF,SYMEND,IERR,DBGFIL,DBGLVL)
C
C Initialize externals - yyparse may be called more than once
C
      YYPV=-1
      YYPS=-1
      YYSTAT=0
      YYTMP=0
      YYVAL=0
      YYNERR=0
      YYERRF=0
      YYCHAR=-1
C
C yystack:
C
 10   YYXPV=YYPV
      YYXPS=YYPS
      YYXSTA=YYSTAT
      IF (DEBUG) THEN
      IF (YYCHAR .EQ. 0) THEN
      WRITE(DBGFIL,'(A,I3,A)') 'YYPAR: State ',YYXSTA,
     1                         ', token end-of-file'
      ELSE IF (YYCHAR .LT. 0) THEN
      WRITE(DBGFIL,'(A,I3,A)') 'YYPAR: State ',YYXSTA,', token -none-'
      ELSE
      YYXI=-1
 101  YYXI=YYXI+1
      IF ((YYTVAL(YYXI) .NE. YYCHAR) .AND.
     1    (YYTVAL(YYXI) .NE. -1)) GO TO 101
      WRITE(DBGFIL,'(A,I3,A,A)') 'YYPAR: State ',YYXSTA,', token ',
     1                           YYTNAM(YYXI)
      ENDIF
      ENDIF
C
C top of for(;;) loop while no reductions done
C
C yy_stack: put a state and value onto the stacks
C
 20   YYXPS=YYXPS+1
      IF (YYXPS .GE. YYMDEP) THEN
        IERR=30
        GO TO 9999
      ENDIF
      YYS(YYXPS)=YYXSTA
      YYXPV=YYXPV+1
      YYV(YYXPV)=YYVAL
C
C yy_newstate: we have a new state - find out what to do
C
 30   YYXN=YYPACT(YYXSTA)
      IF (YYXN .LE. YYFLAG) GO TO 40
      YYTMP=0
      IF (YYCHAR .LT. 0) YYTMP=1
      IF (YYCHAR .LT. 0) THEN
        YYCHAR=YYLEX(INPUT,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF,YYLVAL,
     1               LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,MAXSYM,SYMNAM,
     2               SYMTYP,SYMREF,SYMEND,IERR,DBGFIL,DBGLVL)
        IF (YYCHAR .LT. 0) YYCHAR=0
      ENDIF
      IF (DEBUG .AND. YYTMP.NE.0) THEN
      IF (YYCHAR .EQ. 0) THEN
      WRITE(DBGFIL,'(A)') 'YYPAR: Received token end-of-file'
      ELSE IF (YYCHAR .LT. 0) THEN
      WRITE(DBGFIL,'(A)') 'YYPAR: Received token -none-'
      ELSE
      YYXI=-1
 301  YYXI=YYXI+1
      IF ((YYTVAL(YYXI) .NE. YYCHAR) .AND.
     1    (YYTVAL(YYXI) .NE. -1)) GO TO 301
      WRITE(DBGFIL,'(A,A)') 'YYPAR: Received token ',YYTNAM(YYXI)
      ENDIF
      ENDIF
      YYXN=YYXN+YYCHAR
      IF ((YYXN .LT. 0) .OR. (YYXN .GE. YYLAST)) GO TO 40
      YYXN=YYACT(YYXN)
      IF (YYCHK(YYXN) .EQ. YYCHAR) THEN
        YYCHAR=-1
        YYVAL=YYLVAL
        YYXSTA=YYXN
        IF (YYERRF .GT. 0) YYERRF=YYERRF-1
        GO TO 20
      ENDIF
C
C yydefault:
C
 40   YYXN=YYDEF(YYXSTA)
      IF (YYXN .EQ. -2) THEN
        YYTMP=0
        IF (YYCHAR .LT. 0) YYTMP=1
        IF (YYCHAR .LT. 0) THEN
          YYCHAR=YYLEX(INPUT,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF,YYLVAL,
     1                 LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,MAXSYM,SYMNAM,
     2                 SYMTYP,SYMREF,SYMEND,IERR,DBGFIL,DBGLVL)
          IF (YYCHAR .LT. 0) YYCHAR=0
        ENDIF
      IF (DEBUG .AND. YYTMP.NE.0) THEN
      IF (YYCHAR .EQ. 0) THEN
      WRITE(DBGFIL,'(A)') 'YYPAR: Received token end-of-file'
      ELSE IF (YYCHAR .LT. 0) THEN
      WRITE(DBGFIL,'(A)') 'YYPAR: Received token -none-'
      ELSE
      YYXI=-1
 401  YYXI=YYXI+1
      IF ((YYTVAL(YYXI) .NE. YYCHAR) .AND.
     1    (YYTVAL(YYXI) .NE. -1)) GO TO 401
      WRITE(DBGFIL,'(A,A)') 'YYPAR: Received token ',YYTNAM(YYXI)
      ENDIF
      ENDIF
C
C look through exception table
C
        YYXI=0
 41     IF ((YYEXCA(YYXI) .NE. -1) .OR. (YYEXCA(YYXI+1) .NE.
     1     YYXSTA)) THEN
          YYXI=YYXI+2
          GO TO 41
        ENDIF
 42     YYXI=YYXI+2
        IF ((YYEXCA(YYXI) .GE. 0) .AND. (YYEXCA(YYXI) .NE.
     1     YYCHAR)) GO TO 42
        YYXN=YYEXCA(YYXI+1)
        IF (YYXN .LT. 0) THEN
          IERR=YYACPT
          GO TO 9999
        ENDIF
      ENDIF
C
C check for syntax error
C
      IF (YYXN .EQ. 0) THEN
C
C have an error, switch( yyerrflag )
C
        GO TO (50,51,52,53) YYERRF+1
C
C case 0: new error
C
 50     YYNERR=YYNERR+1
        IF (IERR .EQ. 0) IERR=YYABRT
        GO TO 51
C
C yyerrlab: get globals into registers.
C           we have a user generated syntax type error
C
C501    YYXPV=YYPV
C       YYXPS=YYPS
C       YYXSTA=YYSTAT
C       YYNERR=YYNERR+1
C
C skip_init:
C case 1:
C case 2: incompletely recovered error. try again...
C
 51     CONTINUE
 52     YYERRF=3
C
C find state where "error" is a legal shift action
C
 521    IF (YYXPS .GE. 0) THEN
          YYXN=YYPACT(YYS(YYXPS))+YYERRC
          IF ((YYXN .GE. 0) .AND. (YYXN .LT. YYLAST)) THEN
            IF (YYCHK(YYACT(YYXN)) .EQ. YYERRC) THEN
C
C simulate shift of "error"
C
              YYXSTA=YYACT(YYXN)
              GO TO 20
            ENDIF
          ENDIF
C
C current state has no shift on "error", pop stack
C
          IF (DEBUG) THEN
          IF (YYXPS .GT. 0) THEN
          WRITE(DBGFIL,'(A,I3,A,I3)') 'YYPAR: Error recovery pops '//
     1      'state ',YYS(YYXPS),', uncovers state ',YYS(YYXPS-1)
          ELSE
          WRITE(DBGFIL,'(A,I3)') 'YYPAR: Error recovery pops '//
     1      'state ',YYS(YYXPS)
          ENDIF
          ENDIF
          YYXPS=YYXPS-1
          YYXPV=YYXPV-1
          GO TO 521
        ENDIF
C
C there is no state on stack with "error"
C as a valid shift, give up.
C
        IF (IERR .EQ. 0) IERR=YYABRT
        GO TO 9999
C
C case 3: no shift yet; eat a token
C
 53     IF (DEBUG) THEN
        IF (YYCHAR .EQ. 0) THEN
        WRITE(DBGFIL,'(A,A)') 'YYPAR: Error recovery discards token ',
     1                        'end-of-file'
        ELSE IF (YYCHAR .LT. 0) THEN
        WRITE(DBGFIL,'(A,A)') 'YYPAR: Error recovery discards token ',
     1                        '-none-'
        ELSE
        YYXI=-1
 531    YYXI=YYXI+1
        IF ((YYTVAL(YYXI) .NE. YYCHAR) .AND.
     1      (YYTVAL(YYXI) .NE. -1)) GO TO 531
        WRITE(DBGFIL,'(A,A)') 'YYPAR: Error recovery discards token ',
     1                        YYTNAM(YYXI)
        ENDIF
        ENDIF
        IF (YYCHAR .EQ. 0) THEN
          IF (IERR .EQ. 0) IERR=YYABRT
          GO TO 9999
        ENDIF
        YYCHAR=-1
        GO TO 30
C
C /* end if( yy_n == 0 ) */
C
      ENDIF
C
C reduction by production yy_n
C put stack tops, etc. so things right after switch.
C
      YYTMP=YYXN
      YYPVT=YYXPV
      IF (DEBUG) WRITE(DBGFIL,'(A,I3,A,A,A)') 
     1  'YYPAR: Reduce by (',YYXN,') "',YYREDS(YYXN)(1:YYRLEN(YYXN)),'"'
C
C Look in goto table for next state.
C If yyr2[ yy_n ] doesn't have the low order bit
C set, then there is no action to be done for
C this reduction. So, no saving & unsaving of
C registers done. The only difference between the
C code just after the if and the body of the if is
C the goto yy_stack in the body. This way the test
C can be made before the choice of what to do is needed.
C
C length of production doubled with extra bit.
C
      YYXLEN=YYR2(YYXN)
      IF (MOD(YYXLEN,2) .EQ. 0) THEN
        YYXLEN=YYXLEN/2
        YYXPV=YYXPV-YYXLEN
        YYVAL=YYV(YYXPV+1)
        YYXN=YYR1(YYXN)
        YYXPS=YYXPS-YYXLEN
        YYXSTA=YYPGO(YYXN)+YYS(YYXPS)+1
        IF (YYXSTA .GE. YYLAST) THEN
          YYXSTA=YYACT(YYPGO(YYXN))
        ELSE
          YYXSTA=YYACT(YYXSTA)
          IF (YYCHK(YYXSTA) .NE. -YYXN) YYXSTA=YYACT(YYPGO(YYXN))
        ENDIF
        GO TO 20
      ENDIF
      YYXLEN=YYXLEN/2
      YYXPV=YYXPV-YYXLEN
      YYVAL=YYV(YYXPV+1)
      YYXN=YYR1(YYXN)
      YYXPS=YYXPS-YYXLEN
      YYXSTA=YYPGO(YYXN)+YYS(YYXPS)+1
      IF (YYXSTA .GE. YYLAST) THEN
        YYXSTA=YYACT(YYPGO(YYXN))
      ELSE
        YYXSTA=YYACT(YYXSTA)
        IF (YYCHK(YYXSTA) .NE. -YYXN) YYXSTA=YYACT(YYPGO(YYXN))
      ENDIF
C
C save until reenter driver code
C
      YYSTAT=YYXSTA
      YYPS=YYXPS
      YYPV=YYXPV
C
C code supplied by user is placed in this switch
C
C module : declaration_blocks endmodule {}
C
      IF (YYTMP .EQ. 1) THEN
        IF (GC .NE. 0) THEN
          IF (MC .EQ.0) THEN
            IERR=58
            LNUM=GOTOST(1,1)
            GO TO 9999
          ENDIF
          DO 1010 I=1,GC
            IERR=58
            DO 1000 J=1,MC
              IF (GOTOST(I,1) .EQ. MARKST(J,1)) THEN
                IERR=0
                 CALL PUT1(MARKST(J,2),0.0D0,VPF,GOTOST(I,2)+2,LIWA,
     1               PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
              ENDIF
 1000       CONTINUE
            IF (IERR .NE. 0) THEN
              LNUM=GOTOST(I,1)
              GO TO 9999        
            ENDIF
 1010     CONTINUE
          DO 1020 I=1,GC
            GOTOST(I,1)=0
            GOTOST(I,2)=0
 1020     CONTINUE   
          DO 1030 I=1,MC
            MARKST(I,1)=0
            MARKST(I,2)=0
 1030     CONTINUE   
          MC=0
          GC=0
        ENDIF
        DO 1040 I=1,IWA(14)
          IF (GETIWA(VPF,I,0,IWA,LIWA,INFOLI) .EQ. -ASSIGN) THEN
            CALL PUT1(ASSIGN,0.0D0,VPF,I,LIWA,PLIWA,IWA,LWA,
     1          PLWA,WA,INFOLI,IERR)
            IHELP1=GETIWA(VPF,I+1,0,IWA,LIWA,INFOLI)+IWA(10)
            CALL PUT1(IHELP1,0.0D0,
     1               VPF,I+1,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          ELSE IF (GETIWA(VPF,I,0,IWA,LIWA,INFOLI) .EQ. -FUNC) THEN
            CALL PUT1(FUNC,0.0D0,VPF,I,LIWA,PLIWA,IWA,LWA,
     1          PLWA,WA,INFOLI,IERR)
            IHELP1=GETIWA(VPF,I+1,0,IWA,LIWA,INFOLI)+IWA(10)
            CALL PUT1(IHELP1,0.0D0,
     1          VPF,I+1,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          ELSE IF (GETIWA(VPF,I,0,IWA,LIWA,INFOLI) .EQ. -FUNC-128) THEN
            CALL PUT1(FUNC+128,0.0D0,VPF,I,LIWA,PLIWA,IWA,LWA,
     1          PLWA,WA,INFOLI,IERR)
            IHELP1=GETIWA(VPF,I+1,0,IWA,LIWA,INFOLI)+IWA(10)
            CALL PUT1(IHELP1,0.0D0,
     1          VPF,I+1,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          ENDIF
 1040   CONTINUE
        GO TO 9999
C
C declaration_block : function_head stmts {}
C
      ELSE IF (YYTMP .EQ. 14) THEN
         CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1             IERR)
         IF (IERR .NE. 0) GO TO 9999       
         PC=IWA(14)
C
C param_declaration : ID '=' INUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 18) THEN
        S1=YYV(YYPVT-3)        
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=INT
        VEK4(1)=0
        VEK4(2)=0
        VEK4(3)=0
        VEK4(4)=IWA(4)
        CALL PUT4(VEK4,IIC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999 
        ENDIF
        SYMREF(S1)=IWA(3)
C
C index_declaration : ID '=' index_delimiter RANGE index_delimiter '\n'
C                     {}
C
      ELSE IF (YYTMP .EQ. 22) THEN
        S1=YYV(YYPVT-5)
        S3=YYV(YYPVT-3)
        S5=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=INDEX
        IHELP1=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S5,0,IWA,LIWA,INFOLI)
        IF (IHELP1 .LT. IHELP2) THEN
          VEK5(1)=1
          VEK5(2)=IHELP2-IHELP1+1
          VEK5(3)=IHELP1
          VEK5(4)=IHELP2
          VEK5(5)=IWA(2)+1
          CALL PUT5(VEK5,IIS,0,LIWA,PLIWA,IWA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF 
          DO 1100 I=IHELP1,IHELP2
            CALL PUT1(I,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1                IERR)
 1100     CONTINUE
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF           
        ELSE
          VEK5(1)=1
          VEK5(2)=IHELP1-IHELP2+1
          VEK5(3)=IHELP2
          VEK5(4)=IHELP1
          VEK5(5)=IWA(2)+1
          CALL PUT5(VEK5,IIS,0,LIWA,PLIWA,IWA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF 
          DO 1110 I=IHELP1,IHELP2,-1
            CALL PUT1(I,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1                IERR)
 1110     CONTINUE
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF           
        ENDIF
        SYMREF(S1)=IWA(1)
C
C index_declaration : ID '=' INUM ',' INUM {} opt_inum '\n'
C
      ELSE IF (YYTMP .EQ. 23) THEN
        S1=YYV(YYPVT-4)
        S3=YYV(YYPVT-2)
        S5=YYV(YYPVT)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          GO TO 9999
        ENDIF
        SYMTYP(S1)=INDEX
        IHELP1=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S5,0,IWA,LIWA,INFOLI)
        VEK5(1)=2
        VEK5(2)=2
        VEK5(3)=MIN(IHELP1,IHELP2)
        VEK5(4)=MAX(IHELP1,IHELP2)
        VEK5(5)=IWA(2)+1
        CALL PUT5(VEK5,IIS,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999       
        SYMREF(S1)=IWA(1)
        CALL PUT1(IHELP1,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1            IERR)
        CALL PUT1(IHELP2,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1            IERR)
        IF (IERR .NE. 0) GO TO 9999       
C
C index_declaration : ID '=' ind_expr ',' 
C                     ID '=' index_delimiter RANGE index_delimiter '\n'
C                     {}
C
      ELSE IF (YYTMP .EQ. 25) THEN
        S1=YYV(YYPVT-9)
        S5=YYV(YYPVT-5)
        S7=YYV(YYPVT-3)
        S9=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=INDEX
        VEK5(1)=3
        VEK5(5)=IWA(2)+1
        IF (SYMTYP(S5) .EQ. 0) THEN
          SYMTYP(S5)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S5)=IWA(15)
        ELSE IF (SYMTYP(S5) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(VIC,S7,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S9,0,IWA,LIWA,INFOLI)
        IF (IHELP1 .LT. IHELP2) THEN
          DO 1200 I=IHELP1,IHELP2
            CALL PUT1(I,0.0D0,IV,SYMREF(S5),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1               INFOLI,IERR)
            CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
            IF (IERR .NE. 0) THEN
              LNUM=LNUM-1
              GO TO 9999
            ENDIF
            CALL PUT1(IVAL,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) THEN
              LNUM=LNUM-1
              GO TO 9999
            ENDIF
            IF (I .EQ. IHELP1) THEN
              VEK5(2)=1
              VEK5(3)=IVAL
              VEK5(4)=IVAL
            ELSE
              VEK5(2)=VEK5(2)+1
              VEK5(3)=MIN(IVAL,VEK5(3))
              VEK5(4)=MAX(IVAL,VEK5(4))
            ENDIF
 1200     CONTINUE
          CALL PUT5(VEK5,IIS,0,LIWA,PLIWA,IWA,INFOLI,IERR)  
        ELSE
          DO 1210 I=IHELP1,IHELP2,-1
            CALL PUT1(I,0.0D0,IV,SYMREF(S5),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1          INFOLI,IERR)
            CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
            IF (IERR .NE. 0) THEN
              LNUM=LNUM-1
              GO TO 9999
            ENDIF
            CALL PUT1(IVAL,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) THEN
              LNUM=LNUM-1
              GO TO 9999
            ENDIF
            IF (I .EQ. IHELP1) THEN
              VEK5(2)=1
              VEK5(3)=IVAL
              VEK5(4)=IVAL
            ELSE
              VEK5(2)=VEK5(2)+1
              VEK5(3)=MIN(IVAL,VEK5(3))
              VEK5(4)=MAX(IVAL,VEK5(4))
            ENDIF
 1210     CONTINUE
          CALL PUT5(VEK5,IIS,0,LIWA,PLIWA,IWA,INFOLI,IERR)  
        ENDIF
        SYMREF(S1)=IWA(1)
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
C index_delimiter : ID {}
C
      ELSE IF (YYTMP .EQ. 26) THEN
        S1 = YYV(YYPVT)
        IF (SYMTYP(S1) .NE. INT) THEN
          IERR=8
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIC,SYMREF(S1),4,IWA,LIWA,INFOLI)
        YYVAL = IHELP1
C
C opt_inum : opt_inum ',' INUM {}
C
      ELSE IF (YYTMP .EQ. 28) THEN
        S3=YYV(YYPVT)
        IHELP1=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        CALL PUT1(IHELP1,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1           IERR)  
        IF (IERR .NE. 0) GO TO 9999
        VEK5(1)=GETIWA(IIS,IWA(1),1,IWA,LIWA,INFOLI)
        VEK5(5)=GETIWA(IIS,IWA(1),5,IWA,LIWA,INFOLI)
        VEK5(2)=GETIWA(IIS,IWA(1),2,IWA,LIWA,INFOLI)+1
        IHELP2=GETIWA(IIS,IWA(1),3,IWA,LIWA,INFOLI)
        IHELP3=GETIWA(IIS,IWA(1),4,IWA,LIWA,INFOLI)
        VEK5(3)=MIN(IHELP1,IHELP2)
        VEK5(4)=MAX(IHELP1,IHELP3)
        CALL PUT5(VEK5,IIS,IWA(1),LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C ind_expr : ind_expr '+' ind_expr {}
C
      ELSE IF (YYTMP .EQ. 30) THEN
        CALL PUT1(ADD+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1           IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C ind_expr : ind_expr '-' ind_expr {}
C
      ELSE IF (YYTMP .EQ. 31) THEN
        CALL PUT1(SUB+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1           IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C ind_expr : ind_expr '*' ind_expr {}
C
      ELSE IF (YYTMP .EQ. 32) THEN
        CALL PUT1(MULT+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C ind_expr : ind_expr '/' ind_expr {}
C
      ELSE IF (YYTMP .EQ. 33) THEN
        CALL PUT1(DIV+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1           IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C ind_expr : '-' ind_expr %prec UMINUS {}
C
      ELSE IF (YYTMP .EQ. 35) THEN
        CALL PUT1(UMINUS+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C ind_expr : INUM {}
C
      ELSE IF (YYTMP .EQ. 36) THEN
        CALL PUT1(INUM+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(YYV(YYPVT),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C ind_expr : ID {}
C
      ELSE IF (YYTMP .EQ. 37) THEN
        S1=YYV(YYPVT)
        IF (SYMTYP(S1) .EQ. 0) THEN
          SYMTYP(S1)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMREF(S1)=IWA(15)
        ELSE IF ((SYMTYP(S1) .NE. INDVAR) .AND.
     1           (SYMTYP(S1) .NE. FUNC)) THEN
          IERR=8
          GO TO 9999
        ELSE IF ((SYMTYP(S1) .EQ. FUNC) .AND. 
     1          (GETIWA(IFN,SYMREF(S1),1,IWA,LIWA,INFOLI) .NE. -2)) THEN
          IERR=59
          GO TO 9999
        ENDIF
        IF (SYMTYP(S1) .EQ. INDVAR) THEN
          CALL PUT1(INDVAR+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ELSE IF (SYMTYP(S1) .EQ. FUNC) THEN
          CALL PUT1(-(FUNC+128),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(SYMREF(S1)-IWA(10),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,
     1              PLWA,WA,INFOLI,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 9999
C
C real_head : REAL '\n' {}
C
      ELSE IF (YYTMP .EQ. 38) THEN 
        PSCHK=.FALSE.
C real_declaration : ID '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 41) THEN
        CALL CASE41(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C real_declaration : ID '(' ID ')' '=' expr ',' ID IN ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 42) THEN
        CALL CASE42(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C real_declaration : ID '(' INUM ')' '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 43) THEN
        CALL CASE43(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C real_declaration : ID '(' ID ',' ID ')' '=' expr ','
C                    ID IN ID ',' ID IN ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 44) THEN
        CALL CASE44(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999

C
C real_declaration : ID '(' INUM ',' INUM ')' '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 45) THEN
        CALL CASE45(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C integer_declaration : ID '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 49) THEN
        CALL CASE49(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999

C
C integer_declaration : ID '(' ID ')' '=' expr ',' ID IN ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 50) THEN
        CALL CASE50(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C integer_declaration : ID '(' INUM ')' '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 51) THEN
        CALL CASE51(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C integer_declaration : ID '(' ID ',' ID ')' '=' expr ','
C                       ID IN ID ',' ID IN ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 52) THEN
        CALL CASE52(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C integer_declaration : ID '(' INUM ',' INUM ')' '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 53) THEN
        CALL CASE53(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
        IF (IERR .NE. 0) GO TO 9999
C
C table_head : TABLE ID '(' ID ')' ',' ID IN ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 54) THEN
        CALL CASE54(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,S9)
        IF (IERR .NE. 0) GO TO 9999
C
C table_head : TABLE ID '(' ID ',' ID ')' ',' ID IN ID ',' ID IN ID '\n'
C              {}
C
      ELSE IF (YYTMP .EQ. 55) THEN
        S1=YYV(YYPVT-14)
        S3=YYV(YYPVT-12)
        S5=YYV(YYPVT-10)
        S8=YYV(YYPVT-7)
        S10=YYV(YYPVT-5)
        S12=YYV(YYPVT-3)
        S14=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=REAL
        SYMREF(S1)=IWA(5)+1
        IF (SYMTYP(S3) .EQ. 0) THEN
          SYMTYP(S3)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S3)=IWA(15)
        ELSE IF (SYMTYP(S3) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF (SYMTYP(S5) .EQ. 0) THEN
          SYMTYP(S5)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S5)=IWA(15)
        ELSE IF (SYMTYP(S5) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF ((SYMTYP(S10) .NE. INDEX) .OR. 
     1      (SYMTYP(S14) .NE. INDEX)) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S10),3,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIS,SYMREF(S14),3,IWA,LIWA,INFOLI)
        IF ((S3 .NE. S8) .OR. (S5 .NE. S12) .OR. 
     1      (IHELP1 .LE. 0) .OR. (IHELP2 .LE. 0)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        VEK4(1)=2
        VEK4(2)=GETIWA(IIS,SYMREF(S10),4,IWA,LIWA,INFOLI)
        VEK4(3)=GETIWA(IIS,SYMREF(S14),4,IWA,LIWA,INFOLI)
        VEK4(4)=IWA(6)+1
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IRC,IWA(5),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,IWA(5),3,IWA,LIWA,INFOLI)
        DO 1800 I=1,IHELP1*IHELP2
          CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
 1800   CONTINUE
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
C table_declaration : INUM RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 58) THEN
        S1=YYV(YYPVT-2)
        S2=YYV(YYPVT-1)
        IHELP1=GETIWA(VIC,S1,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIS,SYMREF(S9),4,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 1) .OR.
     1      (IHELP1 .GT. IHELP2)) THEN
          IERR = 33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP3=GETIWA(IRC,IWA(5),4,IWA,LIWA,INFOLI)
        HELP4=GETWA(VRC,S2,WA,LWA,INFOLI)
        CALL PUT1(0,HELP4,VRC,IHELP3+IHELP1-1,LIWA,PLIWA,IWA,LWA,PLWA,
     1            WA,INFOLI,IERR)
C
C table_declaration : INUM '-' RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 59) THEN
        S1=YYV(YYPVT-3)
        S3=YYV(YYPVT-1)
        IHELP1=GETIWA(VIC,S1,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIS,SYMREF(S9),4,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 1) .OR.
     1      (IHELP1 .GT. IHELP2)) THEN
          IERR = 33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP3=GETIWA(IRC,IWA(5),4,IWA,LIWA,INFOLI)
        HELP4=GETWA(VRC,S3,WA,LWA,INFOLI)
        CALL PUT1(0,-HELP4,VRC,IHELP3+IHELP1-1,LIWA,PLIWA,IWA,LWA,PLWA,
     1            WA,INFOLI,IERR)
C
C table_declaration : INUM INUM RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 60) THEN
        S1=YYV(YYPVT-3)
        S2=YYV(YYPVT-2)
        S3=YYV(YYPVT-1)
        IHELP1=GETIWA(VIC,S1,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S2,0,IWA,LIWA,INFOLI)
        IHELP3=GETIWA(IIS,SYMREF(S10),4,IWA,LIWA,INFOLI)
        IHELP4=GETIWA(IIS,SYMREF(S14),4,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 1) .OR. (IHELP1 .GT. IHELP3) .OR.
     1      (IHELP2 .LT. 1) .OR. (IHELP2 .GT. IHELP4)) THEN
          IERR = 33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP5=GETIWA(IRC,IWA(5),4,IWA,LIWA,INFOLI)
        HELP6=GETWA(VRC,S3,WA,LWA,INFOLI)
        CALL PUT1(0,HELP6,VRC,IHELP5+(IHELP1-1)*IHELP4+IHELP2-1,LIWA,
     1            PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
C
C table_declaration : INUM INUM '-' RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 61) THEN
        S1=YYV(YYPVT-4)
        S2=YYV(YYPVT-3)
        S4=YYV(YYPVT-1)
        IHELP1=GETIWA(VIC,S1,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S2,0,IWA,LIWA,INFOLI)
        IHELP3=GETIWA(IIS,SYMREF(S10),4,IWA,LIWA,INFOLI)
        IHELP4=GETIWA(IIS,SYMREF(S14),4,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 1) .OR. (IHELP1 .GT. IHELP3) .OR.
     1      (IHELP2 .LT. 1) .OR. (IHELP2 .GT. IHELP4)) THEN
          IERR = 33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP5=GETIWA(IRC,IWA(5),4,IWA,LIWA,INFOLI)
        HELP6=GETWA(VRC,S4,WA,LWA,INFOLI)
        CALL PUT1(0,-HELP6,VRC,IHELP5+(IHELP1-1)*IHELP4+IHELP2-1,LIWA,
     1            PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
C
C con_interpolation_head : CONINT ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 62) THEN
        IF (SPFLAG .GT. 0) THEN
          IHELP1=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
          DIM=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)/5-1
          CALL SPLNES(WA(IHELP1),DIM,IERR,LNUM)
          IF (IERR .NE. 0) GO TO 9999
          SPFLAG=0
        ENDIF
        S2=YYV(YYPVT-1)
        IPOL=S2
        IF (SYMTYP(S2) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S2)=CONINT
        SYMREF(S2)=IWA(5)+1
        VEK4(1)=1
        VEK4(2)=0
        VEK4(3)=0
        VEK4(4)=IWA(6)+1
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
C lin_interpolation_head : LININT ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 63) THEN
        IF (SPFLAG .GT. 0) THEN
          IHELP1=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
          DIM=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)/5-1
          CALL SPLNES(WA(IHELP1),DIM,IERR,LNUM)
          IF (IERR .NE. 0) GO TO 9999
          SPFLAG=0
        ENDIF
        S2=YYV(YYPVT-1)
        IPOL=S2
        IF (SYMTYP(S2) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S2)=LININT
        SYMREF(S2)=IWA(5)+1
        VEK4(1)=1
        VEK4(2)=0
        VEK4(3)=0
        VEK4(4)=IWA(6)+1
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
C interpolation_declaration : RNUM RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 66) THEN
        S1=YYV(YYPVT-2)
        S2=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S2,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S2,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
        ENDIF
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-3,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=2+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
C interpolation_declaration : RNUM '-' RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 67) THEN
        S1=YYV(YYPVT-3)
        S3=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP4,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ENDIF
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-3,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=2+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
C interpolation_declaration : '-' RNUM RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 68) THEN
        S2=YYV(YYPVT-2)
        S3=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,S2,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ENDIF
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-3,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=2+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
C interpolation_declaration : '-' RNUM '-' RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 69) THEN
        S2=YYV(YYPVT-3)
        S4=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S4,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S4,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S4,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,S2,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,S4,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ENDIF
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-3,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=2+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
C variable_head : VAR '\n' {}
C
      ELSE IF (YYTMP .EQ. 70) THEN
        OIFLAG=1
C
C variable_declaration : ID '(' ID ')' ',' ID IN ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 73) THEN
        S1=YYV(YYPVT-8)
        S3=YYV(YYPVT-6)
        S6=YYV(YYPVT-3)
        S8=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=VAR
        SYMREF(S1)=IWA(7)+1
        IF (SYMTYP(S3) .EQ. 0) THEN
          SYMTYP(S3)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S3)=IWA(15)
        ELSE IF (SYMTYP(S3) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF (SYMTYP(S8) .NE. INDEX) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S8),1,IWA,LIWA,INFOLI)
        IF ((S3 .NE. S6) .OR. (IHELP1 .NE. 1)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        VEK3(1)=1
        VEK3(2)=GETIWA(IIS,SYMREF(S8),2,IWA,LIWA,INFOLI)
        VEK3(3)=IWA(8)+1
        CALL PUT3(VEK3,IVA,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IVA,IWA(7),2,IWA,LIWA,INFOLI)
        DO 2100 I=1,IHELP1
          CALL PUT1(0,0.0D0,VVA,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
 2100   CONTINUE
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
C variable_declaration : ID opt_id '\n' {}
C
      ELSE IF (YYTMP .EQ. 74) THEN
        S1=YYV(YYPVT-2)
        S2=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=VAR
        SYMREF(S1)=S2
        VEK3(1)=0
        VEK3(2)=0
        VEK3(3)=IWA(8)+S2-IWA(7)
        CALL PUT3(VEK3,IVA,S2,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IVA,S2,3,IWA,LIWA,INFOLI)
        CALL PUT1(0,0.0D0,VVA,IHELP1,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1           IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
C opt_id : opt_id ',' ID {}
C
      ELSE IF (YYTMP .EQ. 75) THEN
        S1=YYV(YYPVT-2)
        S3=YYV(YYPVT)
        IF (SYMTYP(S3) .NE. 0) THEN
          IERR=4
          GO TO 9999
        ENDIF
        IF (OIFLAG .EQ. 1) THEN
          CALL PUT1(0,0.0D0,VVA,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1             IERR)
          IF (IERR .NE. 0) GO TO 9999
          VEK3(1)=0
          VEK3(2)=0
          VEK3(3)=IWA(8)
          CALL PUT3(VEK3,IVA,0,LIWA,PLIWA,IWA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMTYP(S3)=VAR
          SYMREF(S3)=IWA(7)
          YYVAL=S1
        ELSE IF (OIFLAG .EQ. 2) THEN
          VEK7(1)=-2
          VEK7(2)=0
          VEK7(3)=0
          VEK7(4)=IWA(11)+1
          VEK7(5)=IWA(12)+1
          VEK7(6)=0
          VEK7(7)=IWA(12)+1
          CALL PUT7(VEK7,IFN,0,LIWA,PLIWA,IWA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMTYP(S3)=FUNC
          SYMREF(S3)=IWA(9)
          IWA(11)=IWA(11)+1
          IF (MODE .GE. GRAD) THEN
            IWA(12)=IWA(12)+1
          ENDIF
          IF (MODE .EQ. HESS) THEN
            IWA(13)=IWA(13)+1
          ENDIF
          IF (IWA(6)+IWA(8)+IWA(11)+IWA(12)*IWA(8)+
     1        IWA(13)*IWA(8)*IWA(8) .GT. LWA) THEN
            IERR=32
            GO TO 9999
          ENDIF
          YYVAL=S1
        ENDIF
C
C opt_id : /* empty */ {}
C
      ELSE IF (YYTMP .EQ. 76) THEN
        IF (OIFLAG .EQ. 1) THEN
          VEK3(1)=0
          VEK3(2)=0
          VEK3(3)=0
          CALL PUT3(VEK3,IVA,0,LIWA,PLIWA,IWA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          CALL PUT1(0,0.0D0,VVA,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1             IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          YYVAL=IWA(7)
        ELSE IF (OIFLAG .EQ. 2) THEN
          VEK7(1)=0
          VEK7(2)=0
          VEK7(3)=0
          VEK7(4)=0
          VEK7(5)=0
          VEK7(6)=0
          VEK7(7)=0
          CALL PUT7(VEK7,IFN,0,LIWA,PLIWA,IWA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          IWA(11)=IWA(11)+1
          IF (MODE .GE. GRAD) THEN
            IWA(12)=IWA(12)+1
          ENDIF
          IF (MODE .EQ. HESS) THEN
            IWA(13)=IWA(13)+1
          ENDIF
          IF (IWA(6)+IWA(8)+IWA(11)+IWA(12)*IWA(8)+
     1        IWA(13)*IWA(8)*IWA(8) .GT. LWA) THEN
            IERR=32
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          YYVAL=IWA(9)
        ENDIF
C
C infunc_head : INFUNC '\n' {}
C
      ELSE IF (YYTMP .EQ. 77) THEN
        OIFLAG=2
C
C infunc_declaration : ID opt_id '\n' {}
C
      ELSE IF (YYTMP .EQ. 80) THEN
        S1=YYV(YYPVT-2)
        S2=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=FUNC
        SYMREF(S1)=S2
        VEK7(1)=-2
        VEK7(2)=0
        VEK7(3)=0
        VEK7(4)=IWA(11)
        VEK7(5)=IWA(12)
        VEK7(6)=0
        VEK7(7)=IWA(12)
        CALL PUT7(VEK7,IFN,S2,LIWA,PLIWA,IWA,INFOLI,IERR)
C
C function_head : FUNC ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 81) THEN
        PSCHK=.TRUE.
        IF (SPFLAG .GT. 0) THEN
          IHELP1=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
          DIM=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)/5-1
          CALL SPLNES(WA(IHELP1),DIM,IERR,LNUM)
          IF (IERR .NE. 0) GO TO 9999
          SPFLAG=0
        ENDIF
        IF (GC .NE. 0) THEN
          IF (MC .EQ.0) THEN
            IERR=58
            LNUM=GOTOST(1,1)
            GO TO 9999
          ENDIF
          DO 2210 I=1,GC
            IERR=58
            DO 2200 J=1,MC
              IF (GOTOST(I,1) .EQ. MARKST(J,1)) THEN
                IERR=0
                 CALL PUT1(MARKST(J,2),0.0D0,VPF,GOTOST(I,2)+2,LIWA,
     1               PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
              ENDIF
 2200       CONTINUE
            IF (IERR .NE. 0) THEN 
              LNUM=GOTOST(I,1)
              GO TO 9999
            ENDIF
 2210     CONTINUE
          DO 2220 I=1,GC
            GOTOST(I,1)=0
            GOTOST(I,2)=0
 2220     CONTINUE   
          DO 2230 I=1,MC
            MARKST(I,1)=0
            MARKST(I,2)=0
 2230     CONTINUE   
          MC=0
          GC=0
        ENDIF
        FUFLAG=0
        S2=YYV(YYPVT-1)
        IF (SYMTYP(S2) .NE. 0) THEN
          IERR=4
          GO TO 9999
        ENDIF
        DO 2240 I=1,SYMEND
          IF ((SYMTYP(I) .EQ. FUNC) .AND. (SYMREF(I) .GT. IWA(10))) THEN
            SYMREF(I)=SYMREF(I)+1
          ENDIF
 2240   CONTINUE
        SYMTYP(S2)=FUNC
        VEK7(1)=0
        VEK7(2)=0
        VEK7(3)=0
        VEK7(4)=0
        VEK7(5)=0
        VEK7(6)=0
        VEK7(7)=0
        CALL PUT7(VEK7,IFN,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        IWA(10)=IWA(10)+1
        SYMREF(S2)=IWA(10)
        DO 2250 I=IWA(9)-1,IWA(10),-1
          VEK7(1)=GETIWA(IFN,I,1,IWA,LIWA,INFOLI)
          VEK7(2)=GETIWA(IFN,I,2,IWA,LIWA,INFOLI)
          VEK7(3)=GETIWA(IFN,I,3,IWA,LIWA,INFOLI)
          VEK7(4)=GETIWA(IFN,I,4,IWA,LIWA,INFOLI)+1
          VEK7(5)=GETIWA(IFN,I,5,IWA,LIWA,INFOLI)+1
          VEK7(6)=GETIWA(IFN,I,6,IWA,LIWA,INFOLI)
          VEK7(7)=GETIWA(IFN,I,7,IWA,LIWA,INFOLI)+1
          CALL PUT7(VEK7,IFN,I+1,LIWA,PLIWA,IWA,INFOLI,IERR)
 2250   CONTINUE
        VEK7(1)=0
        VEK7(2)=0
        VEK7(3)=0
        VEK7(4)=IWA(11)-(IWA(9)-IWA(10))+1
        VEK7(5)=IWA(12)-(IWA(9)-IWA(10))+1
        VEK7(6)=PC+1
        VEK7(7)=IWA(12)-(IWA(9)-IWA(10))+1
        CALL PUT7(VEK7,IFN,IWA(10),LIWA,PLIWA,IWA,INFOLI,IERR)
        IWA(11)=IWA(11)+1
        IF (MODE .GE. GRAD) THEN
          IWA(12)=IWA(12)+1
        ENDIF
        IF (MODE .EQ. HESS) THEN
          IWA(13)=IWA(13)+1
        ENDIF
        IF (IWA(6)+IWA(8)+IWA(11)+IWA(12)*IWA(8)+
     1      IWA(13)*IWA(8)*IWA(8) .GT. LWA) THEN
          IERR=32
          GO TO 9999
        ENDIF
C
C function_head : FUNC ID '(' ID ')' ',' ID IN ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 82) THEN
        PSCHK=.TRUE.
        IF (SPFLAG .GT. 0) THEN
          IHELP1=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
          DIM=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)/5-1
          CALL SPLNES(WA(IHELP1),DIM,IERR,LNUM)
          IF (IERR .NE. 0) GO TO 9999
          SPFLAG=0
        ENDIF
        IF (GC .NE. 0) THEN
          IF (MC .EQ.0) THEN
            IERR=58
            LNUM=GOTOST(1,1)
            GO TO 9999
          ENDIF
          DO 2310 I=1,GC
            IERR=58
            DO 2300 J=1,MC
              IF (GOTOST(I,1) .EQ. MARKST(J,1)) THEN
                IERR=0
                 CALL PUT1(MARKST(J,2),0.0D0,VPF,GOTOST(I,2)+2,LIWA,
     1               PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
              ENDIF
 2300       CONTINUE
            IF (IERR .NE. 0) THEN 
              LNUM=GOTOST(I,1)
              GO TO 9999
            ENDIF
 2310     CONTINUE
          DO 2320 I=1,GC
            GOTOST(I,1)=0
            GOTOST(I,2)=0
 2320     CONTINUE   
          DO 2330 I=1,MC
            MARKST(I,1)=0
            MARKST(I,2)=0
 2330     CONTINUE   
          MC=0
          GC=0
        ENDIF
        FUFLAG=0
        S2=YYV(YYPVT-8)
        S4=YYV(YYPVT-6)
        S7=YYV(YYPVT-3)
        S9=YYV(YYPVT-1)
        IF (SYMTYP(S2) .NE. 0) THEN
          IERR=4
          GO TO 9999
        ENDIF
        DO 2340 I=1,SYMEND
          IF ((SYMTYP(I) .EQ. FUNC) .AND. (SYMREF(I) .GT. IWA(10))) THEN
            SYMREF(I)=SYMREF(I)+1
          ENDIF
 2340   CONTINUE
        IF (SYMTYP(S4) .EQ. 0) THEN
          SYMTYP(S4)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMREF(S4)=IWA(15)
        ELSE IF (SYMTYP(S4) .NE. INDVAR) THEN
          IERR=8
          GO TO 9999
        ENDIF
        IF (SYMTYP(S9) .NE. INDEX) THEN
          IERR=8
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S9),1,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIS,SYMREF(S9),3,IWA,LIWA,INFOLI)
        IF ((S4 .NE. S7) .OR. (IHELP1 .NE. 1) .OR. (IHELP2 .NE. 1)) THEN
          IERR=33
          GO TO 9999
        ENDIF
        SYMTYP(S2)=FUNC
        VEK7(1)=0
        VEK7(2)=0
        VEK7(3)=0
        VEK7(4)=0
        VEK7(5)=0
        VEK7(6)=0
        VEK7(7)=0
        CALL PUT7(VEK7,IFN,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        IWA(10)=IWA(10)+1
        SYMREF(S2)=IWA(10)
        DO 2350 I=IWA(9)-1,IWA(10),-1
          VEK7(1)=GETIWA(IFN,I,1,IWA,LIWA,INFOLI)
          VEK7(2)=GETIWA(IFN,I,2,IWA,LIWA,INFOLI)
          VEK7(3)=GETIWA(IFN,I,3,IWA,LIWA,INFOLI)
          VEK7(4)=GETIWA(IFN,I,4,IWA,LIWA,INFOLI)+
     1            GETIWA(IIS,SYMREF(S9),2,IWA,LIWA,INFOLI)
          VEK7(5)=GETIWA(IFN,I,5,IWA,LIWA,INFOLI)+
     1            GETIWA(IIS,SYMREF(S9),2,IWA,LIWA,INFOLI)
          VEK7(6)=GETIWA(IFN,I,6,IWA,LIWA,INFOLI)
          VEK7(7)=GETIWA(IFN,I,7,IWA,LIWA,INFOLI)+
     1            GETIWA(IIS,SYMREF(S9),2,IWA,LIWA,INFOLI)
          CALL PUT7(VEK7,IFN,I+1,LIWA,PLIWA,IWA,INFOLI,IERR)
 2350   CONTINUE
        VEK7(1)=1
        VEK7(2)=SYMREF(S4)
        VEK7(3)=GETIWA(IIS,SYMREF(S9),2,IWA,LIWA,INFOLI)
        VEK7(4)=IWA(11)-(IWA(9)-IWA(10))+1
        VEK7(5)=IWA(12)-(IWA(9)-IWA(10))+1
        VEK7(6)=PC+1
        VEK7(7)=IWA(12)-(IWA(9)-IWA(10))+1
        CALL PUT7(VEK7,IFN,IWA(10),LIWA,PLIWA,IWA,INFOLI,IERR)
        IWA(11)=IWA(11)+GETIWA(IFN,IWA(10),3,IWA,LIWA,INFOLI)
        IF (MODE .GE. GRAD) THEN
          IWA(12)=IWA(12)+GETIWA(IFN,IWA(10),3,IWA,LIWA,INFOLI)
        ENDIF
        IF (MODE .EQ. HESS) THEN
          IWA(13)=IWA(13)+GETIWA(IFN,IWA(10),3,IWA,LIWA,INFOLI)
        ENDIF
        IF (IWA(6)+IWA(8)+IWA(11)+IWA(12)*IWA(8)+
     1      IWA(13)*IWA(8)*IWA(8) .GT. LWA) THEN
          IERR=32
          GO TO 9999
        ENDIF
C
C
C stmt : ID '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 85) THEN
        S1=YYV(YYPVT-3)
        IF (SYMTYP(S1) .EQ. 0) THEN
          SYMTYP(S1)=FUNC
          VEK7(1)=-1
          VEK7(2)=0
          VEK7(3)=0
          VEK7(4)=IWA(11)+1
          VEK7(5)=IWA(12)+1
          VEK7(6)=0
          VEK7(7)=IWA(12)+1
          CALL PUT7(VEK7,IFN,0,LIWA,PLIWA,IWA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMREF(S1)=IWA(9)
          FUFLAG=0
          IWA(11)=IWA(11)+1
          IF (MODE .GE. GRAD) THEN
            IWA(12)=IWA(12)+1
          ENDIF
          IF (MODE .EQ. HESS) THEN
            IWA(13)=IWA(13)+1
          ENDIF
          IF (IWA(6)+IWA(8)+IWA(11)+IWA(12)*IWA(8)+
     1        IWA(13)*IWA(8)*IWA(8) .GT. LWA) THEN
            IERR=32
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. FUNC) THEN
          IF (GETIWA(IFN,SYMREF(S1),1,IWA,LIWA,INFOLI) .GT. 0) THEN
            IERR=35
            GO TO 9999
          ENDIF
          IF ((GETIWA(IFN,SYMREF(S1),1,IWA,LIWA,INFOLI) .EQ. -2) .AND.
     1        (FUFLAG .EQ. 1)) THEN
            IERR=59
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          FUFLAG=0
        ELSE
          IERR=8
          GO TO 9999
        ENDIF
        IF (GETIWA(IFN,SYMREF(S1),1,IWA,LIWA,INFOLI) .LE. -1) THEN
          CALL PUT1(-ASSIGN,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(SYMREF(S1)-IWA(10),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,
     1              PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ELSE
          CALL PUT1(ASSIGN,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,
     1              PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ENDIF
        PC=IWA(14)
C
C stmt : ID '(' ID ')' '=' expr '\n' {}
C
      ELSE IF (YYTMP .EQ. 86) THEN
        S1=YYV(YYPVT-6)
        S3=YYV(YYPVT-4)
        IF (SYMTYP(S1) .EQ. FUNC) THEN
          IF (GETIWA(IFN,SYMREF(S1),1,IWA,LIWA,INFOLI) .NE. 1) THEN
            IERR=35
            GO TO 9999
          ENDIF
        ELSE
          IERR=8
          GO TO 9999
        ENDIF
        IF (SYMTYP(S3) .NE. INDVAR) THEN
          IERR=8
          GO TO 9999
        ELSE IF (GETIWA(IFN,SYMREF(S1),2,IWA,LIWA,INFOLI) .NE. 
     1           SYMREF(S3)) THEN
          IERR=33
          GO TO 9999
        ENDIF
        CALL PUT1(ASSIGN,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,
     1            PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        PC=IWA(14)
        FUFLAG=0
C
C stmt : IF {} '(' logic_expr ')' THEN '\n' stmts opt_else_if
C        opt_else ENDIF '\n'
C
      ELSE IF (YYTMP .EQ. 87) THEN
        CALL PUT1(IF+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C stmt : IF '(' logic_expr ')' {} THEN '\n' stmts opt_else_if
C        opt_else ENDIF '\n'
C
      ELSE IF (YYTMP .EQ. 88) THEN
        FUFLAG=0
        S1=YYV(YYPVT-3)
        CALL PUT1(BEQ,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(-100,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        YYV(YYPVT-3)=IWA(14)
C
C stmt : IF '(' logic_expr ')' THEN '\n' stmts {} opt_else_if
C        opt_else ENDIF '\n'
C
      ELSE IF (YYTMP .EQ. 89) THEN
        FUFLAG=0
        CALL PUT1(BRA,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(-200,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C stmt : IF '(' logic_expr ')' THEN '\n' stmts opt_else_if
C        opt_else ENDIF '\n' {}
C
      ELSE IF (YYTMP .EQ. 90) THEN
        CALL PUT1(ENDIF+128,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        FUFLAG=0
        S1=YYV(YYPVT-12)
        PC1=S1
 2400   IF ((GETIWA(VPF,PC1,0,IWA,LIWA,INFOLI) .NE. -100) .AND. 
     1       (PC1 .LT. IWA(14))) THEN
          PC1=PC1+1
          GO TO 2400
        ENDIF
        PC2=PC1+1
 2410   IF ((GETIWA(VPF,PC2,0,IWA,LIWA,INFOLI).NE. -200) .AND. 
     1      (PC2 .LT. IWA(14))) THEN
          PC2=PC2+1
          GO TO 2410
        ENDIF
        IF (PC2 .LT. IWA(14)) THEN
          CALL PUT1(PC2+1,0.0D0,VPF,PC1,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
          CALL PUT1(IWA(14)+1,0.0D0,VPF,PC2,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          PC1=PC1+1
          GO TO 2400
        ELSE IF (GETIWA(VPF,PC1-1,0,IWA,LIWA,INFOLI) .EQ. BEQ) THEN
          CALL PUT1(IWA(14)+1,0.0D0,VPF,PC1,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ENDIF
C
C stmt : MARKE CONTINUE '\n' {}
C
      ELSE IF (YYTMP .EQ. 91) THEN
        S1=YYV(YYPVT-2)
        IF (MC .GT. 0) THEN
          IHELP1=GETIWA(VIC,S1,0,IWA,LIWA,INFOLI)
          DO 2500 I=1,MC
            IF (MARKST(I,1) .EQ. IHELP1) THEN
              IERR=57
              LNUM=IHELP1
            ENDIF
 2500     CONTINUE
          IF (IERR .NE. 0) GO TO 9999
        ENDIF
        MC=MC+1
        IF (MC .GT. MAXMAR) THEN
          IERR=32
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(VIC,S1,0,IWA,LIWA,INFOLI)
        MARKST(MC,1)=IHELP1
        MARKST(MC,2)=IWA(14)+1
        CALL PUT1(CONTIN,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1            IERR)
        CALL PUT1(IHELP1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1            IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
C stmt : GOTO INUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 92) THEN
        S2=YYV(YYPVT-1)
        GC=GC+1
        IF (GC .GT. MAXMAR) THEN
          IERR=32
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(VIC,S2,0,IWA,LIWA,INFOLI)
        GOTOST(GC,1)=IHELP1
        GOTOST(GC,2)=IWA(14)+1
        CALL PUT1(GOTO,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1            IERR)
        CALL PUT1(IHELP1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1            IERR)
        CALL PUT1(0,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1            IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
C opt_else_if : opt_else_if ELSE IF '(' logic_expr ')' {}
C               THEN stmts '\n'
C
      ELSE IF (YYTMP .EQ. 93) THEN
        FUFLAG=0
        CALL PUT1(BEQ,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(-100,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C opt_else_if : opt_else_if ELSE IF '(' logic_expr ')'
C               THEN stmts '\n' {}
C
      ELSE IF (YYTMP .EQ. 94) THEN
        CALL PUT1(BRA,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(-200,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C opt_else : /* empty */ {}
C
      ELSE IF (YYTMP .EQ. 97) THEN
        CALL UNVPF(IWA(14)-2,LIWA,PLIWA,IWA,INFOLI)
C
C expr : expr '+' expr {}
C
      ELSE IF (YYTMP .EQ. 98) THEN
        CALL PUT1(ADD,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C expr : expr '-' expr {}
C
      ELSE IF (YYTMP .EQ. 99) THEN
        CALL PUT1(SUB,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C expr : expr '*' expr {}
C
      ELSE IF (YYTMP .EQ. 100) THEN
        CALL PUT1(MULT,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C expr : expr '/' expr {}
C
      ELSE IF (YYTMP .EQ. 101) THEN
        FUFLAG=1
        CALL PUT1(DIV,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C expr : expr '^' expr {}
C
      ELSE IF (YYTMP .EQ. 102) THEN
        FUFLAG=1
        IHELP1=GETIWA(VPF,IWA(14),0,IWA,LIWA,INFOLI)
        IF (IHELP1 .EQ. UMINUS) THEN
          IHELP1=GETIWA(VPF,IWA(14)-2,0,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. INUM) THEN
            CALL PUT1(POWER+128,0.0D0,VPF,IWA(14)-2,LIWA,PLIWA,IWA,LWA,
     1                PLWA,WA,INFOLI,IERR)
            IHELP2=GETIWA(VPF,IWA(14)-1,0,IWA,LIWA,INFOLI)
            IHELP3=-GETIWA(VIC,IHELP2,0,IWA,LIWA,INFOLI)
            CALL PUT1(IHELP3,0.0D0,VPF,IWA(14)-1,LIWA,PLIWA,IWA,LWA,
     1                PLWA,WA,INFOLI,IERR)
            CALL UNVPF(IWA(14)-1,LIWA,PLIWA,IWA,INFOLI)
          ELSE
            CALL PUT1(POWER,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ENDIF
        ELSE
          IHELP1=GETIWA(VPF,IWA(14)-1,0,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. INUM) THEN
            CALL PUT1(POWER+128,0.0D0,VPF,IWA(14)-1,LIWA,PLIWA,IWA,LWA,
     1                PLWA,WA,INFOLI,IERR)
            IHELP2=GETIWA(VPF,IWA(14),0,IWA,LIWA,INFOLI)
            IHELP3=GETIWA(VIC,IHELP2,0,IWA,LIWA,INFOLI)
            CALL PUT1(IHELP3,0.0D0,VPF,IWA(14),LIWA,PLIWA,IWA,LWA,PLWA,
     1                WA,INFOLI,IERR)
          ELSE
            CALL PUT1(POWER,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ENDIF
        ENDIF
C
C expr : '-' expr %prec UMINUS {}
C
      ELSE IF (YYTMP .EQ. 104) THEN
        CALL PUT1(UMINUS,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C expr : SUM {} '(' expr ',' ID IN ID ')'
C
      ELSE IF (YYTMP .EQ. 110) THEN
        FUFLAG=1
        CALL PUT1(SUM,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(0,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(0,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        YYV(YYPVT)=IWA(14)-2
        PSLVL=PSLVL+1
C
C expr : SUM '(' expr ',' ID IN ID ')' {}
C
      ELSE IF (YYTMP .EQ. 111) THEN
        FUFLAG=1
        S1=YYV(YYPVT-8)
        S6=YYV(YYPVT-3)
        S8=YYV(YYPVT-1)
        IF (SYMTYP(S6) .EQ. 0) THEN
          SYMTYP(S6)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMREF(S6)=IWA(15)
        ELSE IF (SYMTYP(S6) .NE. INDVAR) THEN
          IERR=8
          GO TO 9999
        ENDIF
        IF (SYMTYP(S8) .NE. INDEX) THEN
          IERR=8
          GO TO 9999
        ENDIF
        CALL PUT1(SYMREF(S6),0.0D0,VPF,S1+1,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(SYMREF(S8),0.0D0,VPF,S1+2,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(ENDSUM,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        PSLVL=PSLVL-1
C
C expr : PROD {} '(' expr ',' ID IN ID ')'
C
      ELSE IF (YYTMP .EQ. 112) THEN
        FUFLAG=1
        CALL PUT1(PROD,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(0,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(0,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        YYV(YYPVT)=IWA(14)-2
        PSLVL=PSLVL+1
C
C expr : PROD '(' expr ',' ID IN ID ')' {}
C
      ELSE IF (YYTMP .EQ. 113) THEN
        FUFLAG=1
        S1=YYV(YYPVT-8)
        S6=YYV(YYPVT-3)
        S8=YYV(YYPVT-1)
        IF (SYMTYP(S6) .EQ. 0) THEN
          SYMTYP(S6)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMREF(S6)=IWA(15)
        ELSE IF (SYMTYP(S6) .NE. INDVAR) THEN
          IERR=8
          GO TO 9999
        ENDIF
        IF (SYMTYP(S8) .NE. INDEX) THEN
          IERR=8
          GO TO 9999
        ENDIF
        CALL PUT1(SYMREF(S6),0.0D0,VPF,S1+1,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(SYMREF(S8),0.0D0,VPF,S1+2,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(ENDPRD,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
        PSLVL=PSLVL+1
C
C logic_expr : logic_expr AND logic_expr {}
C
      ELSE IF (YYTMP .EQ. 114) THEN
        CALL PUT1(AND,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C logic_expr : logic_expr OR logic_expr {}
C
      ELSE IF (YYTMP .EQ. 115) THEN
        CALL PUT1(OR,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C logic_expr : NOT logic_expr {}
C
      ELSE IF (YYTMP .EQ. 116) THEN
        CALL PUT1(NOT,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C logic_expr : expr RELOP expr {}
C
      ELSE IF (YYTMP .EQ. 118) THEN
        CALL PUT1(RELOP,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(YYV(YYPVT-1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C number : RNUM {}
C
      ELSE IF (YYTMP .EQ. 119) THEN
        FUFLAG=1
        CALL PUT1(RNUM,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(YYV(YYPVT),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C number : INUM {}
C
      ELSE IF (YYTMP .EQ. 120) THEN
        CALL PUT1(INUM,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(YYV(YYPVT),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C identifier : ID {}
C
      ELSE IF (YYTMP .EQ. 121) THEN
        S1=YYV(YYPVT)
        IF (SYMTYP(S1) .EQ. 0) THEN
          IF (PSLVL.EQ.0 .AND. PSCHK) THEN
            IERR=7
            GOTO 9999
          ENDIF
          SYMTYP(S1)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
          SYMREF(S1)=IWA(15)
        ENDIF
        IF (SYMTYP(S1) .EQ. INDVAR) THEN
          CALL PUT1(INDVAR,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
        ELSE IF (SYMTYP(S1) .EQ. REAL) THEN
          FUFLAG=1
          IHELP1=GETIWA(IRC,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 0) THEN
            CALL PUT1(REAL,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. INT) THEN
          IHELP1=GETIWA(IIC,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 0) THEN
            CALL PUT1(INT,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. VAR) THEN
          FUFLAG=1
          IHELP1=GETIWA(IVA,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 0) THEN
            CALL PUT1(VAR,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. FUNC) THEN
          IHELP1=GETIWA(IFN,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. -2) THEN
            CALL PUT1(-FUNC,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1)-IWA(10),0.0D0,VPF,0,LIWA,PLIWA,IWA,
     1                LWA,PLWA,WA,INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE IF (IHELP1 .EQ. -1) THEN
            FUFLAG=1
            CALL PUT1(-FUNC,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1)-IWA(10),0.0D0,VPF,0,LIWA,PLIWA,IWA,
     1                LWA,PLWA,WA,INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE IF (IHELP1 .EQ. 0) THEN
            FUFLAG=1
            CALL PUT1(FUNC,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,
     1                LWA,PLWA,WA,INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE
          IERR=8
          GO TO 9999
        ENDIF
C
C identifier : ID '(' ind_expr ')' {}
C
      ELSE IF (YYTMP .EQ. 122) THEN
        S1=YYV(YYPVT-3)
        IF (SYMTYP(S1) .EQ. REAL) THEN
          FUFLAG=1
          IHELP1=GETIWA(IRC,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 1) THEN
            CALL PUT1(REAL,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. INT) THEN
          IHELP1=GETIWA(IIC,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 1) THEN
            CALL PUT1(INT,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. VAR) THEN
          FUFLAG=1
          IHELP1=GETIWA(IVA,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 1) THEN
            CALL PUT1(VAR,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. FUNC) THEN
          FUFLAG=1
          IHELP1=GETIWA(IFN,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 1) THEN
            CALL PUT1(FUNC,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. 0) THEN
          IERR=7
          GO TO 9999
        ELSE
          IERR=8
          GO TO 9999
        ENDIF
C
C identifier : ID '(' ind_expr ',' ind_expr ')' {}
C
      ELSE IF (YYTMP .EQ. 123) THEN
        S1=YYV(YYPVT-5)
        IF (SYMTYP(S1) .EQ. REAL) THEN
          FUFLAG=1
          IHELP1=GETIWA(IRC,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 2) THEN
            CALL PUT1(REAL,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. INT) THEN
          IHELP1=GETIWA(IIC,SYMREF(S1),1,IWA,LIWA,INFOLI)
          IF (IHELP1 .EQ. 2) THEN
            CALL PUT1(INT,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,IERR)
            IF (IERR .NE. 0) GO TO 9999
          ELSE
            IERR=35
            GO TO 9999
          ENDIF
        ELSE IF (SYMTYP(S1) .EQ. 0) THEN
          IERR=7
          GO TO 9999
        ELSE
          IERR=8
          GO TO 9999
        ENDIF
C
C standard_function : STDRD {}
C
      ELSE IF (YYTMP .EQ. 124) THEN
        FUFLAG=1
        S1=YYV(YYPVT)
        IF (STDTYP(S1) .EQ. 0) THEN
          CALL PUT1(STDRD,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(S1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ELSE
          IERR=36
          GO TO 9999
        ENDIF
C
C standard_function : STDRD '(' expr ')' {}
C
      ELSE IF (YYTMP .EQ. 125) THEN
        FUFLAG=1
        S1=YYV(YYPVT-3)
        IF (STDTYP(S1) .EQ. 1) THEN
          CALL PUT1(STDRD,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(S1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ELSE
          IERR=36
          GO TO 9999
        ENDIF
C
C standard_function : STDRD '(' expr ',' expr ')' {}
C
      ELSE IF (YYTMP .EQ. 126) THEN
        FUFLAG=1
        S1=YYV(YYPVT-5)
        IF (STDTYP(S1) .EQ. 2) THEN
          CALL PUT1(STDRD,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(S1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ELSE
          IERR=36
          GO TO 9999
        ENDIF
C
C extern_function : EXTERN {}
C
      ELSE IF (YYTMP .EQ. 127) THEN
        FUFLAG=1
        S1=YYV(YYPVT)
        IF (EXTTYP(S1) .EQ. 0) THEN
          CALL PUT1(EXTERN,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(S1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ELSE
          IERR=36
          GO TO 9999
        ENDIF
C
C extern_function : EXTERN '(' ind_expr ')' {}
C
      ELSE IF (YYTMP .EQ. 128) THEN
        FUFLAG=1
        S1=YYV(YYPVT-3)
        IF (EXTTYP(S1) .EQ. 1) THEN
          CALL PUT1(EXTERN,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(S1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ELSE
          IERR=36
          GO TO 9999
        ENDIF
C
C extern_function : EXTERN '(' ind_expr ',' ind_expr ')' {}
C
      ELSE IF (YYTMP .EQ. 129) THEN
        FUFLAG=1
        S1=YYV(YYPVT-5)
        IF (EXTTYP(S1) .EQ. 2) THEN
          CALL PUT1(EXTERN,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(S1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          IF (IERR .NE. 0) GO TO 9999
        ELSE
          IERR=36
          GO TO 9999
        ENDIF
C
C interpolation_function : INTERP '(' expr ')' {}
C
      ELSE IF (YYTMP .EQ. 130) THEN
        FUFLAG=1
        S1=YYV(YYPVT-3)
        IF (SYMTYP(S1) .EQ. CONINT) THEN
          CALL PUT1(CONINT,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ELSE IF (SYMTYP(S1) .EQ. LININT) THEN
          CALL PUT1(LININT,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ELSE IF (SYMTYP(S1) .EQ. SPLINE) THEN
          CALL PUT1(SPLINE,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ELSE
          IERR=36
          GO TO 9999
        ENDIF
        CALL PUT1(SYMREF(S1),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        CALL PUT1(IWA(15),0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1            INFOLI,IERR)
        IF (IERR .NE. 0) GO TO 9999
C
C spl_interpolation_head : SPLINE ID '\n' {}
C
      ELSE IF (YYTMP .EQ. 131) THEN
        IF (SPFLAG .GT. 0) THEN
          IHELP1=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
          DIM=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)/5-1
          CALL SPLNES(WA(IHELP1),DIM,IERR,LNUM)
          IF (IERR .NE. 0) GO TO 9999
          SPFLAG=0
        ENDIF
        CALL CAS131(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1              LIWA,PLIWA,IWA,INFOLI,IPOL,SPFLAG)
        IF (IERR .NE. 0) GO TO 9999
C
C spline_declaration : RNUM RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 134) THEN
        CALL CAS134(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1              LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
        IF (IERR .NE. 0) GO TO 9999
C
C spline_declaration : RNUM '-' RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 135) THEN
        CALL CAS135(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1              LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
        IF (IERR .NE. 0) GO TO 9999
C
C spline_declaration : '-' RNUM RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 136) THEN
        CALL CAS136(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1              LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
        IF (IERR .NE. 0) GO TO 9999
C
C spline_declaration : '-' RNUM '-' RNUM '\n' {}
C
      ELSE IF (YYTMP .EQ. 137) THEN
        CALL CAS137(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1              LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
        IF (IERR .NE. 0) GO TO 9999
C
C
C
      ELSE
        IERR=26
        WRITE(*,*) 'YYPAR (2833) : unknown state ',YYTMP
        GO TO 9999
      ENDIF
C
C reset registers in driver code
C
      GO TO 10
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE INIRED (YYREDS)
      CHARACTER*100 YYREDS(0:138)
C
      YYREDS(0)='-no such reduction-'
      YYREDS(1)='module : declaration_blocks end_module'
      YYREDS(2)='declaration_blocks : declaration_blocks '//
     &  'declaration_block'
      YYREDS(3)='declaration_blocks : /* empty */'
      YYREDS(4)='declaration_block : param_head param_declarations'
      YYREDS(5)='declaration_block : index_head index_declarations'
      YYREDS(6)='declaration_block : real_head real_declarations'
      YYREDS(7)='declaration_block : integer_head integer_declarations'
      YYREDS(8)='declaration_block : table_head table_declarations'
      YYREDS(9)='declaration_block : con_interpolation_head '//
     &  'interpolation_declarations'
      YYREDS(10)='declaration_block : lin_interpolation_head '//
     &  'interpolation_declarations'
      YYREDS(11)='declaration_block : spl_interpolation_head '//
     &  'spline_declarations'
      YYREDS(12)='declaration_block : variable_head '//
     &  'variable_declarations'
      YYREDS(13)='declaration_block : infunc_head infunc_declarations'
      YYREDS(14)='declaration_block : function_head stmts'
      YYREDS(15)='param_head : PARAM ''\n'''
      YYREDS(16)='param_declarations : param_declarations '//
     &  'param_declaration'
      YYREDS(17)='param_declarations : /* empty */'
      YYREDS(18)='param_declaration : ID ''='' INUM ''\n'''
      YYREDS(19)='index_head : INDEX ''\n'''
      YYREDS(20)='index_declarations : index_declarations '//
     &  'index_declaration'
      YYREDS(21)='index_declarations : /* empty */'
      YYREDS(22)='index_declaration : ID ''='' index_delimiter '//
     &  'RANGE index_delimiter ''\n'''
      YYREDS(23)='index_declaration : ID ''='' INUM '','' INUM'
      YYREDS(24)='index_declaration : ID ''='' INUM '','' INUM '//
     &  'opt_inum ''\n'''
      YYREDS(25)='index_declaration : ID ''='' ind_expr '','' ID '//
     &  '''='' index_delimiter RANGE index_delimiter ''\n'''
      YYREDS(26)='index_delimiter : ID'
      YYREDS(27)='index_delimiter : INUM'
      YYREDS(28)='opt_inum : opt_inum '','' INUM'
      YYREDS(29)='opt_inum : /* empty */'
      YYREDS(30)='ind_expr : ind_expr ''+'' ind_expr'
      YYREDS(31)='ind_expr : ind_expr ''-'' ind_expr'
      YYREDS(32)='ind_expr : ind_expr ''*'' ind_expr'
      YYREDS(33)='ind_expr : ind_expr ''/'' ind_expr'
      YYREDS(34)='ind_expr : ''('' ind_expr '')'''
      YYREDS(35)='ind_expr : ''-'' ind_expr'
      YYREDS(36)='ind_expr : INUM'
      YYREDS(37)='ind_expr : ID'
      YYREDS(38)='real_head : REAL ''\n'''
      YYREDS(39)='real_declarations : real_declarations '//
     &  'real_declaration'
      YYREDS(40)='real_declarations : /* empty */'
      YYREDS(41)='real_declaration : ID ''='' expr ''\n'''
      YYREDS(42)='real_declaration : ID ''('' ID '')'' ''='' expr '//
     &  ''','' ID IN ID ''\n'''
      YYREDS(43)='real_declaration : ID ''('' INUM '')'' ''='' expr '//
     &  '''\n'''
      YYREDS(44)='real_declaration : ID ''('' ID '','' ID '')'' '//
     &  '''='' expr '','' IDIN ID '','' ID IN ID ''\n'''
      YYREDS(45)='real_declaration : ID ''('' INUM '','' INUM '')'' '//
     &  '''='' expr ''\n'''
      YYREDS(46)='integer_head : INT ''\n'''
      YYREDS(47)='integer_declarations : integer_declarations '//
     &  'integer_declaration'
      YYREDS(48)='integer_declarations : /* empty */'
      YYREDS(49)='integer_declaration : ID ''='' expr ''\n'''
      YYREDS(50)='integer_declaration : ID ''('' ID '')'' ''='' '//
     &  'expr '','' ID IN ID ''\n'''
      YYREDS(51)='integer_declaration : ID ''('' INUM '')'' ''='' '//
     &  'expr ''\n'''
      YYREDS(52)='integer_declaration : ID ''('' ID '','' ID '')'' '//
     &  '''='' expr '','' ID IN ID '','' ID IN ID ''\n'''
      YYREDS(53)='integer_declaration : ID ''('' INUM '','' INUM '//
     &  ''')'' ''='' expr ''\n'''
      YYREDS(54)='table_head : TABLE ID ''('' ID '')'' '','' ID IN '//
     &  'ID ''\n'''
      YYREDS(55)='table_head : TABLE ID ''('' ID '','' ID '')'' '//
     &  ''','' ID IN ID '','' ID IN ID ''\n'''
      YYREDS(56)='table_declarations : table_declarations '//
     &  'table_declaration'
      YYREDS(57)='table_declarations : /* empty */'
      YYREDS(58)='table_declaration : INUM RNUM ''\n'''
      YYREDS(59)='table_declaration : INUM ''-'' RNUM ''\n'''
      YYREDS(60)='table_declaration : INUM INUM RNUM ''\n'''
      YYREDS(61)='table_declaration : INUM INUM ''-'' RNUM ''\n'''
      YYREDS(62)='con_interpolation_head : CONINT ID ''\n'''
      YYREDS(63)='lin_interpolation_head : LININT ID ''\n'''
      YYREDS(64)='interpolation_declarations : '//
     &  'interpolation_declarations interpolation_declaration'
      YYREDS(65)='interpolation_declarations : /* empty */'
      YYREDS(66)='interpolation_declaration : RNUM RNUM ''\n'''
      YYREDS(67)='interpolation_declaration : RNUM ''-'' RNUM ''\n'''
      YYREDS(68)='interpolation_declaration : ''-'' RNUM RNUM ''\n'''
      YYREDS(69)='interpolation_declaration : ''-'' RNUM ''-'' RNUM '//
     &  '''\n'''
      YYREDS(70)='variable_head : VAR ''\n'''
      YYREDS(71)='variable_declarations : variable_declarations '//
     &  'variable_declaration'
      YYREDS(72)='variable_declarations : /* empty */'
      YYREDS(73)='variable_declaration : ID ''('' ID '')'' '','' ID '//
     &  'IN ID ''\n'''
      YYREDS(74)='variable_declaration : ID opt_id ''\n'''
      YYREDS(75)='opt_id : opt_id '','' ID'
      YYREDS(76)='opt_id : /* empty */'
      YYREDS(77)='infunc_head : INFUNC ''\n'''
      YYREDS(78)='infunc_declarations : infunc_declarations '//
     &  'infunc_declaration'
      YYREDS(79)='infunc_declarations : /* empty */'
      YYREDS(80)='infunc_declaration : ID opt_id ''\n'''
      YYREDS(81)='function_head : FUNC ID ''\n'''
      YYREDS(82)='function_head : FUNC ID ''('' ID '')'' '','' ID '//
     &  'IN ID ''\n'''
      YYREDS(83)='stmts : stmts stmt'
      YYREDS(84)='stmts : /* empty */'
      YYREDS(85)='stmt : ID ''='' expr ''\n'''
      YYREDS(86)='stmt : ID ''('' ID '')'' ''='' expr ''\n'''
      YYREDS(87)='stmt : IF'
      YYREDS(88)='stmt : IF ''('' logic_expr '')'''
      YYREDS(89)='stmt : IF ''('' logic_expr '')'' THEN ''\n'' stmts'
      YYREDS(90)='stmt : IF ''('' logic_expr '')'' THEN ''\n'' '//
     &  'stmts opt_else_if opt_else ENDIF ''\n'''
      YYREDS(91)='stmt : LABEL CONTINUE ''\n'''
      YYREDS(92)='stmt : GOTO INUM ''\n'''
      YYREDS(93)='opt_else_if : opt_else_if ELSE IF ''('' '//
     &  'logic_expr '')'''
      YYREDS(94)='opt_else_if : opt_else_if ELSE IF ''('' '//
     &  'logic_expr '')'' THEN ''\n'' stmts'
      YYREDS(95)='opt_else_if : /* empty */'
      YYREDS(96)='opt_else : ELSE ''\n'' stmts'
      YYREDS(97)='opt_else : /* empty */'
      YYREDS(98)='expr : expr ''+'' expr'
      YYREDS(99)='expr : expr ''-'' expr'
      YYREDS(100)='expr : expr ''*'' expr'
      YYREDS(101)='expr : expr ''/'' expr'
      YYREDS(102)='expr : expr ''^'' expr'
      YYREDS(103)='expr : ''('' expr '')'''
      YYREDS(104)='expr : ''-'' expr'
      YYREDS(105)='expr : number'
      YYREDS(106)='expr : identifier'
      YYREDS(107)='expr : standard_function'
      YYREDS(108)='expr : extern_function'
      YYREDS(109)='expr : interpolation_function'
      YYREDS(110)='expr : SUM'
      YYREDS(111)='expr : SUM ''('' expr '','' ID IN ID '')'''
      YYREDS(112)='expr : PROD'
      YYREDS(113)='expr : PROD ''('' expr '','' ID IN ID '')'''
      YYREDS(114)='logic_expr : logic_expr AND logic_expr'
      YYREDS(115)='logic_expr : logic_expr OR logic_expr'
      YYREDS(116)='logic_expr : NOT logic_expr'
      YYREDS(117)='logic_expr : ''('' logic_expr '')'''
      YYREDS(118)='logic_expr : expr RELOP expr'
      YYREDS(119)='number : RNUM'
      YYREDS(120)='number : INUM'
      YYREDS(121)='identifier : ID'
      YYREDS(122)='identifier : ID ''('' ind_expr '')'''
      YYREDS(123)='identifier : ID ''('' ind_expr '','' ind_expr '')'''
      YYREDS(124)='standard_function : STANDARD'
      YYREDS(125)='standard_function : STANDARD ''('' expr '')'''
      YYREDS(126)='standard_function : STANDARD ''('' expr '','' '//
     &  'expr '')'''
      YYREDS(127)='extern_function : EXTERN'
      YYREDS(128)='extern_function : EXTERN ''('' ind_expr '')'''
      YYREDS(129)='extern_function : EXTERN ''('' ind_expr '','' '//
     &  'ind_expr '')'''
      YYREDS(130)='interpolation_function : INTERP ''('' expr '')'''
      YYREDS(131)='spl_interpolation_head : SPLINE ID ''\n'''
      YYREDS(132)='spline_declarations : spline_declarations '//
     &  'spline_declaration'
      YYREDS(133)='spline_declarations : /* empty */'
      YYREDS(134)='spline_declaration : RNUM RNUM ''\n'''
      YYREDS(135)='spline_declaration : RNUM ''-'' RNUM ''\n'''
      YYREDS(136)='spline_declaration : ''-'' RNUM RNUM ''\n'''
      YYREDS(137)='spline_declaration : ''-'' RNUM ''-'' RNUM ''\n'''
      YYREDS(138)='end_module : END ''\n'''
      RETURN
      END
      SUBROUTINE CASE41(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER VEK4(4)
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1
      INTEGER IVAL
C
      INTEGER REAL
      PARAMETER (REAL=277)
C
        S1=YYV(YYPVT-3)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=REAL
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(0,FVAL,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
        VEK4(1)=0
        VEK4(2)=0
        VEK4(3)=0
        VEK4(4)=IWA(6)
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMREF(S1)=IWA(5)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE42(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,IHELP3,IHELP4,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3,S8,S10
      INTEGER IVAL,I
C
      INTEGER INDEX,REAL,INDVAR
      PARAMETER (INDEX=276,REAL=277,INDVAR=291)
C
        S1=YYV(YYPVT-10)
        S3=YYV(YYPVT-8)
        S8=YYV(YYPVT-3)
        S10=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=REAL
        SYMREF(S1)=IWA(5)+1
        IF (SYMTYP(S3) .EQ. 0) THEN
          SYMTYP(S3)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S3)=IWA(15)
        ELSE IF (SYMTYP(S3) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF (SYMTYP(S10) .NE. INDEX) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S10),3,IWA,LIWA,INFOLI)
        IF ((S3 .NE. S8) .OR. (IHELP1 .LE. 0)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        VEK4(1)=1
        VEK4(2)=GETIWA(IIS,SYMREF(S10),4,IWA,LIWA,INFOLI)
        VEK4(3)=0
        VEK4(4)=IWA(6)+1
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IRC,IWA(5),2,IWA,LIWA,INFOLI)
        DO 1300 I=1,IHELP1
          CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
 1300   CONTINUE
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S10),2,IWA,LIWA,INFOLI)
        DO 1310 I=1,IHELP1
          IHELP2=GETIWA(IIS,SYMREF(S10),5,IWA,LIWA,INFOLI)
          IHELP3=GETIWA(VIS,IHELP2+I-1,0,IWA,LIWA,INFOLI)
          CALL PUT1(IHELP3,0.0D0,IV,SYMREF(S3),LIWA,PLIWA,IWA,LWA,PLWA,
     1              WA,INFOLI,IERR)
          CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          IHELP1=GETIWA(IRC,IWA(5),2,IWA,LIWA,INFOLI)
          IHELP4=GETIWA(IV,SYMREF(S3),0,IWA,LIWA,INFOLI)
          CALL PUT1(0,FVAL,VRC,IWA(6)-IHELP1+IHELP4,LIWA,PLIWA,IWA,LWA,
     1              PLWA,WA,INFOLI,IERR)
 1310   CONTINUE
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE43(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER IHELP1,IHELP2,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3
      INTEGER IVAL
C
      INTEGER REAL
      PARAMETER (REAL=277)
C
        S1=YYV(YYPVT-6)
        S3=YYV(YYPVT-4)
        IF (SYMTYP(S1) .NE. REAL) THEN
          IERR=7
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(S1),2,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 1) .OR. 
     1      (IHELP1 .GT. IHELP2)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IRC,SYMREF(S1),4,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        CALL PUT1(0,FVAL,VRC,IHELP1+IHELP2-1,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1      INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE44(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,IHELP3,IHELP4,IHELP5,IHELP6,IHELP7
      INTEGER IHELP8,IHELP9,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3,S5,S10,S12,S14,S16
      INTEGER IVAL,I,J
C
      INTEGER INDEX,REAL,INDVAR
      PARAMETER (INDEX=276,REAL=277,INDVAR=291)
C
        S1=YYV(YYPVT-16)
        S3=YYV(YYPVT-14)
        S5=YYV(YYPVT-12)
        S10=YYV(YYPVT-7)
        S12=YYV(YYPVT-5)
        S14=YYV(YYPVT-3)
        S16=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=REAL
        SYMREF(S1)=IWA(5)+1
        IF (SYMTYP(S3) .EQ. 0) THEN
          SYMTYP(S3)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S3)=IWA(15)
        ELSE IF (SYMTYP(S3) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF (SYMTYP(S5) .EQ. 0) THEN
          SYMTYP(S5)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S5)=IWA(15)
        ELSE IF (SYMTYP(S5) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF ((SYMTYP(S12) .NE. INDEX) .OR. 
     1      (SYMTYP(S16) .NE. INDEX)) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S12),3,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIS,SYMREF(S16),3,IWA,LIWA,INFOLI)
        IF ((S3 .NE. S10) .OR. (S5 .NE. S14) .OR. 
     1      (IHELP1 .LE. 0) .OR. (IHELP2 .LE. 0)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        VEK4(1)=2
        VEK4(2)=GETIWA(IIS,SYMREF(S12),4,IWA,LIWA,INFOLI)
        VEK4(3)=GETIWA(IIS,SYMREF(S16),4,IWA,LIWA,INFOLI)
        VEK4(4)=IWA(6)+1
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IRC,IWA(5),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,IWA(5),3,IWA,LIWA,INFOLI)
        DO 1400 I=1,IHELP1*IHELP2
          CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
 1400   CONTINUE
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S12),2,IWA,LIWA,INFOLI)
        DO 1420 I=1,IHELP1
          IHELP2=GETIWA(IIS,SYMREF(S12),5,IWA,LIWA,INFOLI)
          IHELP3=GETIWA(VIS,IHELP2+I-1,0,IWA,LIWA,INFOLI)
          CALL PUT1(IHELP3,0.0D0,IV,SYMREF(S3),LIWA,PLIWA,IWA,LWA,PLWA,
     1              WA,INFOLI,IERR)
          IHELP4=GETIWA(IIS,SYMREF(S16),2,IWA,LIWA,INFOLI)
          DO 1410 J=1,IHELP4
            IHELP5=GETIWA(IIS,SYMREF(S16),5,IWA,LIWA,INFOLI)
            IHELP6=GETIWA(VIS,IHELP5+J-1,0,IWA,LIWA,INFOLI)
            CALL PUT1(IHELP6,0.0D0,IV,SYMREF(S5),LIWA,PLIWA,IWA,LWA,
     1                PLWA,WA,INFOLI,IERR)
            CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
            IF (IERR .NE. 0) THEN
              LNUM=LNUM-1
              GO TO 9999
            ENDIF
          IHELP1=GETIWA(IRC,IWA(5),2,IWA,LIWA,INFOLI)
          IHELP2=GETIWA(IRC,IWA(5),3,IWA,LIWA,INFOLI)
          IHELP7=GETIWA(IV,SYMREF(S3),0,IWA,LIWA,INFOLI)
          IHELP8=GETIWA(IRC,IWA(5),3,IWA,LIWA,INFOLI)
          IHELP9=GETIWA(IV,SYMREF(S5),0,IWA,LIWA,INFOLI)
          CALL PUT1(0,FVAL,VRC,
     1              IWA(6)-(IHELP1*IHELP2)+(IHELP7-1)*IHELP8+IHELP9,
     2              LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
 1410     CONTINUE
 1420   CONTINUE
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE45(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER IHELP1,IHELP2,IHELP3,IHELP4,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3,S5
      INTEGER IVAL
C
      INTEGER REAL
      PARAMETER (REAL=277)
C
        S1=YYV(YYPVT-8)
        S3=YYV(YYPVT-6)
        S5=YYV(YYPVT-4)
        IF (SYMTYP(S1) .NE. REAL) THEN
          IERR=7
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S5,0,IWA,LIWA,INFOLI)
        IHELP3=GETIWA(IRC,SYMREF(S1),2,IWA,LIWA,INFOLI)
        IHELP4=GETIWA(IRC,SYMREF(S1),3,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 1) .OR. (IHELP1 .GT. IHELP3) .OR.
     1      (IHELP2 .LT. 1) .OR. (IHELP2 .GT. IHELP4)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IRC,SYMREF(S1),4,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP3=GETIWA(IRC,SYMREF(S1),3,IWA,LIWA,INFOLI)
        IHELP4=GETIWA(VIC,S5,0,IWA,LIWA,INFOLI)
        CALL PUT1(0,FVAL,VRC,IHELP1+(IHELP2-1)*IHELP3+IHELP4-1,
     1            LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE49(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER VEK4(4)
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1
      INTEGER IVAL
C
      INTEGER INT
      PARAMETER (INT=278)
C
        S1=YYV(YYPVT-3)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=INT
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(IDNINT(FVAL),0.0D0,VIC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1      INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
        VEK4(1)=0
        VEK4(2)=0
        VEK4(3)=0
        VEK4(4)=IWA(4)
        CALL PUT4(VEK4,IIC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMREF(S1)=IWA(3)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE50(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,IHELP3,IHELP4,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3,S8,S10
      INTEGER IVAL
C
      INTEGER INDEX,INT,INDVAR,I
      PARAMETER (INDEX=276,INT=278,INDVAR=291)
C
        S1=YYV(YYPVT-10)
        S3=YYV(YYPVT-8)
        S8=YYV(YYPVT-3)
        S10=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=INT
        SYMREF(S1)=IWA(3)+1
        IF (SYMTYP(S3) .EQ. 0) THEN
          SYMTYP(S3)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S3)=IWA(15)
        ELSE IF (SYMTYP(S3) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF (SYMTYP(S10) .NE. INDEX) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S10),3,IWA,LIWA,INFOLI)
        IF ((S3 .NE. S8) .OR. (IHELP1 .LE. 0)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        VEK4(1)=1
        VEK4(2)=GETIWA(IIS,SYMREF(S10),4,IWA,LIWA,INFOLI)
        VEK4(3)=0
        VEK4(4)=IWA(4)+1
        CALL PUT4(VEK4,IIC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIC,IWA(3),2,IWA,LIWA,INFOLI)
        DO 1500 I=1,IHELP1
          CALL PUT1(0,0.0D0,VIC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
 1500   CONTINUE
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S10),2,IWA,LIWA,INFOLI)
        DO 1510 I=1,IHELP1
          IHELP2=GETIWA(IIS,SYMREF(S10),5,IWA,LIWA,INFOLI)
          IHELP3=GETIWA(VIS,IHELP2+I-1,0,IWA,LIWA,INFOLI)
          CALL PUT1(IHELP3,0.0D0,IV,SYMREF(S3),LIWA,PLIWA,IWA,LWA,PLWA,
     1              WA,INFOLI,IERR)
          CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          IHELP1=GETIWA(IIC,IWA(3),2,IWA,LIWA,INFOLI)
          IHELP4=GETIWA(IV,SYMREF(S3),0,IWA,LIWA,INFOLI)
          CALL PUT1(IDNINT(FVAL),0.0D0,VIC,IWA(4)-IHELP1+IHELP4,LIWA,
     1              PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
 1510   CONTINUE
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE51(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER IHELP1,IHELP2,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3
      INTEGER IVAL
C
      INTEGER INT
      PARAMETER (INT=278)
C
        S1=YYV(YYPVT-6)
        S3=YYV(YYPVT-4)
        IF (SYMTYP(S1) .NE. INT) THEN
          IERR=7
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIC,SYMREF(S1),2,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 0) .OR. (IHELP1 .GT. IHELP2)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIC,SYMREF(S1),4,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        CALL PUT1(IDNINT(FVAL),0.0D0,VIC,IHELP1+IHELP2-1,LIWA,PLIWA,
     1            IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE52(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,IHELP3,IHELP4,IHELP5,IHELP6,IHELP7
      INTEGER IHELP8,IHELP9,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3,S5,S10,S12,S14,S16
      INTEGER IVAL
C
      INTEGER INDEX,INT,INDVAR,I,J
      PARAMETER (INDEX=276,INT=278,INDVAR=291)
C
        S1=YYV(YYPVT-16)
        S3=YYV(YYPVT-14)
        S5=YYV(YYPVT-12)
        S10=YYV(YYPVT-7)
        S12=YYV(YYPVT-5)
        S14=YYV(YYPVT-3)
        S16=YYV(YYPVT-1)
        IF (SYMTYP(S1) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S1)=INT
        SYMREF(S1)=IWA(3)+1
        IF (SYMTYP(S3) .EQ. 0) THEN
          SYMTYP(S3)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S3)=IWA(15)
        ELSE IF (SYMTYP(S3) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF (SYMTYP(S5) .EQ. 0) THEN
          SYMTYP(S5)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S5)=IWA(15)
        ELSE IF (SYMTYP(S5) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF ((SYMTYP(S12) .NE. INDEX) .OR. 
     1      (SYMTYP(S16) .NE. INDEX)) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S12),3,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIS,SYMREF(S16),3,IWA,LIWA,INFOLI)
        IF ((S3 .NE. S10) .OR. (S5 .NE. S14) .OR. 
     1      (IHELP1 .LE. 0) .OR. (IHELP2 .LE. 0)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        VEK4(1)=2
        VEK4(2)=GETIWA(IIS,SYMREF(S12),4,IWA,LIWA,INFOLI)
        VEK4(3)=GETIWA(IIS,SYMREF(S16),4,IWA,LIWA,INFOLI)
        VEK4(4)=IWA(4)+1
        CALL PUT4(VEK4,IIC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIC,IWA(3),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IIC,IWA(3),3,IWA,LIWA,INFOLI)
        DO 1600 I=1,IHELP1*IHELP2
          CALL PUT1(0,0.0D0,VIC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
 1600   CONTINUE
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S12),2,IWA,LIWA,INFOLI)
        DO 1620 I=1,IHELP1
          IHELP2=GETIWA(IIS,SYMREF(S12),5,IWA,LIWA,INFOLI)
          IHELP3=GETIWA(VIS,IHELP2+I-1,0,IWA,LIWA,INFOLI)
          CALL PUT1(IHELP3,0.0D0,IV,SYMREF(S3),LIWA,PLIWA,IWA,LWA,PLWA,
     1              WA,INFOLI,IERR)
          IHELP4=GETIWA(IIS,SYMREF(S16),2,IWA,LIWA,INFOLI)
          DO 1610 J=1,IHELP4
            IHELP5=GETIWA(IIS,SYMREF(S16),5,IWA,LIWA,INFOLI)
            IHELP6=GETIWA(VIS,IHELP5+J-1,0,IWA,LIWA,INFOLI)
            CALL PUT1(IHELP6,0.0D0,IV,SYMREF(S5),LIWA,PLIWA,IWA,LWA,
     1                PLWA,WA,INFOLI,IERR)
            CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
            IF (IERR .NE. 0) THEN
              LNUM=LNUM-1
              GO TO 9999
            ENDIF
          IHELP1=GETIWA(IIC,IWA(3),2,IWA,LIWA,INFOLI)
          IHELP2=GETIWA(IIC,IWA(3),3,IWA,LIWA,INFOLI)
          IHELP7=GETIWA(IV,SYMREF(S3),0,IWA,LIWA,INFOLI)
          IHELP8=GETIWA(IIC,IWA(3),3,IWA,LIWA,INFOLI)
          IHELP9=GETIWA(IV,SYMREF(S5),0,IWA,LIWA,INFOLI)
          CALL PUT1(IDNINT(FVAL),0.0D0,VIC,
     1             IWA(4)-(IHELP1*IHELP2)+(IHELP7-1)*IHELP8+IHELP9,
     2             LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
 1610     CONTINUE
 1620   CONTINUE
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE53(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,PC)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
C
      INTEGER IHELP1,IHELP2,IHELP3,IHELP4,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      DOUBLE PRECISION FVAL
      INTEGER S1,S3,S5
      INTEGER IVAL
C
      INTEGER INT
      PARAMETER (INT=278)
C
        S1=YYV(YYPVT-8)        
        S3=YYV(YYPVT-6)        
        S5=YYV(YYPVT-4)        
        IF (SYMTYP(S1) .NE. INT) THEN
          IERR=7
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S5,0,IWA,LIWA,INFOLI)
        IHELP3=GETIWA(IIC,SYMREF(S1),2,IWA,LIWA,INFOLI)
        IHELP4=GETIWA(IIC,SYMREF(S1),3,IWA,LIWA,INFOLI)
        IF ((IHELP1 .LT. 1) .OR. (IHELP1 .GT. IHELP3) .OR.
     1      (IHELP2 .LT. 1) .OR. (IHELP2 .GT. IHELP4)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL PUT1(-1,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        CALL EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIC,SYMREF(S1),4,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(VIC,S3,0,IWA,LIWA,INFOLI)
        IHELP3=GETIWA(IIC,SYMREF(S1),3,IWA,LIWA,INFOLI)
        IHELP4=GETIWA(VIC,S5,0,IWA,LIWA,INFOLI)
        CALL PUT1(IDNINT(FVAL),0.0D0,VIC,
     1           IHELP1+(IHELP2-1)*IHELP3+IHELP4-1,
     2           LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CASE54(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,S9)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15),S9,I
C
      INTEGER VEK4(4)
      INTEGER IHELP1,GETIWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER S2,S4,S7
C
      INTEGER INDEX,REAL,INDVAR
      PARAMETER (INDEX=276,REAL=277,INDVAR=291)
C
        S2=YYV(YYPVT-8)
        S4=YYV(YYPVT-6)
        S7=YYV(YYPVT-3)
        S9=YYV(YYPVT-1)
        IF (SYMTYP(S2) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S2)=REAL
        SYMREF(S2)=IWA(5)+1
        IF (SYMTYP(S4) .EQ. 0) THEN
          SYMTYP(S4)=INDVAR
          CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
          IF (IERR .NE. 0) THEN
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
          SYMREF(S4)=IWA(15)
        ELSE IF (SYMTYP(S4) .NE. INDVAR) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IF (SYMTYP(S9) .NE. INDEX) THEN
          IERR=8
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IIS,SYMREF(S9),3,IWA,LIWA,INFOLI)
        IF ((S4 .NE. S7) .OR. (IHELP1 .LE. 0)) THEN
          IERR=33
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        VEK4(1)=1
        VEK4(2)=GETIWA(IIS,SYMREF(S9),4,IWA,LIWA,INFOLI)
        VEK4(3)=0
        VEK4(4)=IWA(6)+1
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        IHELP1=GETIWA(IRC,IWA(5),2,IWA,LIWA,INFOLI)
        DO 1700 I=1,IHELP1
          CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
 1700   CONTINUE
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CAS131(YYV,YYPVT,SYMTYP,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,INFOLI,IPOL,SPFLAG)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMTYP(MAXSYM),SYMREF(MAXSYM)
      INTEGER LIWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      INTEGER INFOLI(15)
      INTEGER IPOL,SPFLAG
C
      INTEGER VEK4(4)
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER S2
C
      INTEGER SPLINE
      PARAMETER (SPLINE=282)
C
        SPFLAG=1
        S2=YYV(YYPVT-1)
        IPOL=S2
        IF (SYMTYP(S2) .NE. 0) THEN
          IERR=4
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
        SYMTYP(S2)=SPLINE
        SYMREF(S2)=IWA(5)+1
        VEK4(1)=1
        VEK4(2)=0
        VEK4(3)=0
        VEK4(4)=IWA(6)+1
        CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
        IF (IERR .NE. 0) THEN
          LNUM=LNUM-1
          GO TO 9999
        ENDIF
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CAS134(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER IPOL
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,GETIWA
      DOUBLE PRECISION HELP3,HELP4,GETWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER S1,S2
C
        S1=YYV(YYPVT-2)
        S2=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S2,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S2,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
        ENDIF
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-9,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=5+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CAS135(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER IPOL
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,GETIWA
      DOUBLE PRECISION HELP3,HELP4,GETWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER S1,S3
C
        S1=YYV(YYPVT-3)
        S3=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=GETWA(VRC,S1,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP4,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ENDIF
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-9,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=5+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CAS136(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER IPOL
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,GETIWA
      DOUBLE PRECISION HELP3,HELP4,GETWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER S2,S3
C
        S2=YYV(YYPVT-2)
        S3=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=GETWA(VRC,S3,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,S2,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ENDIF
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-9,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=5+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CAS137(YYV,YYPVT,SYMREF,MAXSYM,IERR,LNUM,
     1                  LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IPOL)
C
      INTEGER YYPVT,YYV(0:149)
      INTEGER MAXSYM,SYMREF(MAXSYM)
      INTEGER LWA,LIWA,PLWA,PLIWA,IERR,LNUM
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER IPOL
C
      INTEGER VEK4(4)
      INTEGER IHELP1,IHELP2,GETIWA
      DOUBLE PRECISION HELP3,HELP4,GETWA
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER S2,S4
C
        S2=YYV(YYPVT-3)
        S4=YYV(YYPVT-1)
        IHELP1=GETIWA(IRC,SYMREF(IPOL),2,IWA,LIWA,INFOLI)
        IHELP2=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 2 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S4,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE IF ( ((IHELP1+2)-(IWA(6)-IHELP2)-1) .EQ. 1 ) THEN
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S4,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,IWA(6),LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     1              IERR)
        ELSE
          HELP3=-GETWA(VRC,S2,WA,LWA,INFOLI)
          HELP4=-GETWA(VRC,S4,WA,LWA,INFOLI)
          CALL PUT1(0,HELP3,VRC,S2,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
          CALL PUT1(0,HELP4,VRC,S4,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,IERR)
        ENDIF
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IHELP1 .GT. 0) THEN
          HELP4=GETWA(VRC,IWA(6)-9,WA,LWA,INFOLI)
          IF (HELP3 .LE. HELP4) THEN
            IERR= 63
            LNUM=LNUM-1
            GO TO 9999
          ENDIF
        ENDIF
        VEK4(1)=1
        VEK4(2)=5+IHELP1
        VEK4(3)=0
        VEK4(4)=GETIWA(IRC,SYMREF(IPOL),4,IWA,LIWA,INFOLI)
        CALL PUT4(VEK4,IRC,SYMREF(IPOL),LIWA,PLIWA,IWA,INFOLI,IERR)
C
 9999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE SYMINP (INPUT,SYMFIL,WA,LWA,IWA,LIWA,PLWA,PLIWA,
     1                   IERR,LNUM,MODE,NVAR,NFUNC,DBGFIL,DBGLVL)
      INTEGER INPUT,SYMFIL
      INTEGER LWA,LIWA
      DOUBLE PRECISION WA(LWA)
      INTEGER IWA(LIWA)
      INTEGER PLWA,PLIWA
      INTEGER IERR,LNUM,MODE,NVAR,NFUNC,DBGFIL,DBGLVL
C
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
      INTEGER INFOLI(15),VEK3(3),VEK4(4),VEK5(5),VEK7(7)
      INTEGER I,GETIWA
C
      INTEGER OUTPUT
C     PARAMETER (OUTPUT=7)
C
      INTEGER GSMDEP
      PARAMETER (GSMDEP=10)
C
      INTEGER HSMDEP
      PARAMETER (HSMDEP=10)
C
      INTEGER NOGRAD,GRAD,HESS
      PARAMETER (NOGRAD=0,GRAD=1,HESS=2)
C
C-----------------------------------------------------------------------
      CALL NEXTID(OUTPUT)
C-----------------------------------------------------------------------
      PLIWA=15
      PLWA=0
      NVAR=0
      NFUNC=0
      DO 1 I=1,3
        VEK3(I)=0
1     CONTINUE
      DO 2 I=1,4
        VEK4(I)=0
2     CONTINUE
      DO 3 I=1,5
        VEK5(I)=0
3     CONTINUE
      DO 4 I=1,7
        VEK7(I)=0
4     CONTINUE
      DO 8 I=1,LWA
        WA(I)=0.0D0
8     CONTINUE
      DO 9 I=1,LIWA
        IWA(I)=0
9     CONTINUE
C
      IERR=1
      LNUM=0
      OPEN(OUTPUT,STATUS='SCRATCH',ERR=99)
      CALL YYPRE(INPUT,OUTPUT,IERR,LNUM)
      REWIND(OUTPUT)
      IF (IERR .NE. 0) RETURN
      CALL YYPAR (OUTPUT,WA,LWA,IWA,LIWA,PLWA,PLIWA,MODE,IERR,LNUM,
     1            DBGFIL,DBGLVL)
      CLOSE(OUTPUT)
      IF (IERR .NE. 0) RETURN
C
      INFOLI(1)=15
      INFOLI(2)=IWA(1)*5+INFOLI(1)
      INFOLI(3)=IWA(2)+INFOLI(2)
      INFOLI(4)=IWA(3)*4+INFOLI(3)
      INFOLI(5)=IWA(4)+INFOLI(4)
      INFOLI(7)=IWA(5)*4+INFOLI(5)
      INFOLI(9)=IWA(7)*3+INFOLI(7)
      INFOLI(14)=IWA(9)*7+INFOLI(9)
      INFOLI(15)=IWA(14)+INFOLI(14)
C
      INFOLI(6)=0
      INFOLI(8)=IWA(6)+INFOLI(6)
      INFOLI(11)=IWA(8)+INFOLI(8)
      INFOLI(12)=IWA(11)+INFOLI(11)
      INFOLI(13)=IWA(12)*IWA(8)+INFOLI(12)
C
      IF (IWA(1) .EQ. 0) 
     1  CALL PUT5(VEK5,IIS,0,LIWA,PLIWA,IWA,INFOLI,IERR)
      IF (IWA(2) .EQ. 0) 
     1  CALL PUT1(0,0.0D0,VIS,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
      IF (IWA(3) .EQ. 0) 
     1  CALL PUT4(VEK4,IIC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
      IF (IWA(4) .EQ. 0) 
     1  CALL PUT1(0,0.0D0,VIC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
      IF (IWA(5) .EQ. 0) 
     1  CALL PUT4(VEK4,IRC,0,LIWA,PLIWA,IWA,INFOLI,IERR)
      IF (IWA(6) .EQ. 0) 
     1  CALL PUT1(0,0.0D0,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
      IF (IWA(7) .EQ. 0) 
     1  CALL PUT3(VEK3,IVA,0,LIWA,PLIWA,IWA,INFOLI,IERR)
      IF (IWA(8) .EQ. 0) 
     1  CALL PUT1(0,0.0D0,VVA,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
      IF (IWA(9) .EQ. 0) 
     1  CALL PUT7(VEK7,IFN,0,LIWA,PLIWA,IWA,INFOLI,IERR)
      IF (IWA(11) .EQ. 0)
     1  CALL PUT1(0,0.0D0,VFN,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
C     IF (IWA(12) .EQ. 0) 
C    1  CALL PUT1(0,0.0D0,VGR,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
C     IF (IWA(13) .EQ. 0) 
C    1  IWA(13)=IWA(13)+1
      IF (IWA(14) .EQ. 0) 
     1  CALL PUT1(0,0.0D0,VPF,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
      IF (IWA(15) .EQ. 0) 
     1  CALL PUT1(0,0.0D0,IV,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
C
      NVAR=IWA(8)
      NFUNC=IWA(10)
      DO 50 I=1,IWA(9)
      IF (GETIWA(IFN,I,1,IWA,LIWA,INFOLI).EQ.1) THEN
        NFUNC=NFUNC + GETIWA(IFN,I,3,IWA,LIWA,INFOLI)-1
      ENDIF
   50 CONTINUE
      IF (SYMFIL.GT.0) THEN
        DO 90 I=1,PLIWA-IWA(15)
          WRITE(SYMFIL,'(I6)',ERR=95) IWA(I)
 90     CONTINUE
        DO 91 I=1,IWA(6)
          WRITE(SYMFIL,'(D24.17)',ERR=95) WA(I)
 91     CONTINUE
      ENDIF
      IF (MODE .EQ. NOGRAD) THEN
        PLWA=IWA(6)+IWA(8)+IWA(11)+GSMDEP*IWA(8)+HSMDEP+2
      ELSE IF (MODE .EQ. GRAD) THEN
        PLWA=IWA(6)+IWA(8)+IWA(11)+IWA(11)*IWA(8)+GSMDEP*IWA(8)+
     1       HSMDEP+1
      ELSE IF (MODE .EQ. HESS) THEN
        PLWA=IWA(6)+IWA(8)+IWA(11)+IWA(11)*IWA(8)+GSMDEP*IWA(8)+
     1       IWA(11)*IWA(8)*IWA(8)+HSMDEP*IWA(8)*IWA(8)
      ELSE
        IERR=60
      ENDIF
      IF (PLWA .GT. LWA) THEN
        IERR=32
        RETURN
      ENDIF
      RETURN
 95   IERR=26
 99   RETURN
      END
C
C
C
      SUBROUTINE KLEIN (D,BUCHST)
      CHARACTER D 
      LOGICAL BUCHST
C 
      CHARACTER DAK,DZK
      DAK='a'
      DZK='z'
      BUCHST = ((ICHAR(DAK) .LE. ICHAR(D)) .AND.
     1          (ICHAR(D) .LE. ICHAR(DZK)))
      RETURN
      END
C
C
C
      SUBROUTINE UPCASE (D)
      CHARACTER D
C
      CHARACTER DAK,DAG
      DAK='a'
      DAG='A'
      D = CHAR(ICHAR(D)+ICHAR(DAG)-ICHAR(DAK))
      RETURN
      END
C
C
C
      CHARACTER FUNCTION GETCHR (FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
      INTEGER FD,LNUM,LPOS,SPOS
      CHARACTER STAR*5,CONT*1,LINE*66
      LOGICAL EOF
C
      CHARACTER D
      LOGICAL BUCHST
C
      IF (LPOS .GE. 67) THEN
 5      LNUM=LNUM+1
        READ(FD,'(A5,A1,A66)',END=10) STAR,CONT,LINE
        IF ((STAR(1:1) .EQ. 'C') .OR. (STAR(1:1) .EQ. 'c') .OR.
     1    (LINE .EQ. ' ')) GO TO 5
        IF (CONT .EQ. ' ' .OR. CONT .EQ. '0') THEN
          GETCHR='#'
          LPOS=1
          SPOS=2
          RETURN
        ELSE
          LPOS=1
          SPOS=2
        ENDIF
      ENDIF
      IF (STAR(1:1) .EQ. '*') THEN
        GETCHR='@'
        STAR(1:1)=' '
        RETURN
      ENDIF
      IF (SPOS .LE. 5) THEN
        D=STAR(SPOS:SPOS)
        CALL KLEIN(D,BUCHST)
        IF (BUCHST) CALL UPCASE(D)
        GETCHR=D
        SPOS=SPOS+1
        RETURN
      ENDIF
      D=LINE(LPOS:LPOS)
      LPOS=LPOS+1
      CALL KLEIN(D,BUCHST)
      IF (BUCHST) CALL UPCASE(D)
      GETCHR=D
      RETURN
 10   GETCHR='#'
      EOF=.TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE UNGTC (LPOS,IERR)
C
      INTEGER LPOS
      INTEGER IERR
C
      IF (LPOS .GT. 1) THEN
        LPOS=LPOS-1
      ELSE
        IERR=48
      ENDIF
      RETURN
      END
C
C
C
      INTEGER FUNCTION KEYWD (STRING,YYLVAL,IERR)
      CHARACTER*(*) STRING
      INTEGER YYLVAL,IERR
C
      INTEGER RELOP,AND,OR,NOT,ID,SUM,PROD,IN,IF,THEN,ELSE
      INTEGER ENDIF,STDRD,EXTERN,INTERP,GOTO,CONTIN
      PARAMETER (RELOP=258,AND=259,OR=260,NOT=261,ID=264)
      PARAMETER (SUM=265,PROD=266,IN=267,IF=268,THEN=269)
      PARAMETER (ELSE=270,ENDIF=271,STDRD=272,EXTERN=273)
      PARAMETER (INTERP=274,GOTO=287,CONTIN=289)
C
      INTEGER MAXSTD,MAXEXT
      PARAMETER (MAXSTD=17,MAXEXT=1)
C
      CHARACTER*6 STDNAM(MAXSTD),ALTNAM(MAXSTD),EXTNAM(MAXEXT)
      CHARACTER*4 RELNAM(6)
C
      INTEGER I
C
      DATA STDNAM /'DABS','DSQRT','DEXP','DLOG','DLOG10',
     1             'DSIN','DCOS','DTAN','DASIN','DACOS','DATAN',
     2             'DSINH','DCOSH','DTANH','DASINH','DACOSH','DATANH'/
      DATA ALTNAM /'ABS','SQRT','EXP','LOG','LOG10',
     1             'SIN','COS','TAN','ASIN','ACOS','ATAN',
     2             'SINH','COSH','TANH','ASINH','ACOSH','ATANH'/
      DATA EXTNAM /'XXXXXX'/
      DATA RELNAM /'.EQ.','.NE.','.LT.','.LE.','.GT.','.GE.'/
C
      IF (STRING(1:1) .NE. '.') THEN
        IF (STRING .EQ. 'SUM') THEN
          KEYWD=SUM
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'PROD') THEN
          KEYWD=PROD
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'IN') THEN
          KEYWD=IN
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'IF') THEN
          KEYWD=IF
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'THEN') THEN
          KEYWD=THEN
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'ELSE') THEN
          KEYWD=ELSE
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'ENDIF') THEN
          KEYWD=ENDIF
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'GOTO') THEN
          KEYWD=GOTO
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. 'CONTINUE') THEN
          KEYWD=CONTIN
          YYLVAL=0
          RETURN
        ENDIF
        DO 10 I=1,MAXSTD
          IF ((STRING .EQ. STDNAM(I)) .OR. (STRING .EQ. ALTNAM(I))) THEN
            KEYWD=STDRD
            YYLVAL=I
            RETURN
          ENDIF
 10     CONTINUE
        DO 20 I=1,MAXEXT
          IF (STRING .EQ. EXTNAM(I)) THEN
            KEYWD=EXTERN
            YYLVAL=I
            RETURN
          ENDIF
 20     CONTINUE
        KEYWD=ID
        YYLVAL=0
        RETURN
      ELSE
        DO 30 I=1,6
          IF (STRING .EQ. RELNAM(I)) THEN
            KEYWD=RELOP
            YYLVAL=I
            RETURN
          ENDIF
 30     CONTINUE
        IF (STRING .EQ. '.AND.') THEN
          KEYWD=AND
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. '.OR.') THEN
          KEYWD=OR
          YYLVAL=0
          RETURN
        ENDIF
        IF (STRING .EQ. '.NOT.') THEN
          KEYWD=NOT
          YYLVAL=0
          RETURN
        ENDIF
        KEYWD=0
        IERR=11
        RETURN
      ENDIF
      END
C
C
C
      INTEGER FUNCTION INSERT (STRING,STRLEN,TYPE,LIWA,PLIWA,IWA,LWA,
     1                         PLWA,WA,INFOLI,MAXSYM,SYMNAM,SYMTYP,
     2                         SYMREF,SYMEND,IERR)
      CHARACTER*(*) STRING
      INTEGER STRLEN,TYPE
      INTEGER LIWA,PLIWA,LWA,PLWA
      INTEGER IWA(LIWA),INFOLI(15)
      DOUBLE PRECISION WA(LWA)
      INTEGER MAXSYM
      CHARACTER*20 SYMNAM(MAXSYM)
      INTEGER SYMTYP(MAXSYM),SYMREF(MAXSYM),SYMEND
      INTEGER IERR
C
      INTEGER INUM,RNUM,ID,GETIWA
      PARAMETER (INUM=262,RNUM=263,ID=264)
C
      INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
      PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
      PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C
      INTEGER I,X
      INTEGER IHELP
      DOUBLE PRECISION Y,GETWA,RHELP
      CHARACTER*40 COPY
C
      INSERT=0
      IF ((TYPE .EQ. INUM) .OR. (TYPE .EQ. RNUM)) THEN
        COPY=' '
        COPY(40-STRLEN+1:40)=STRING(1:STRLEN)
      ENDIF
      IF (TYPE .EQ. INUM) THEN
        READ(COPY,'(I40)',ERR=11) X
        DO 10 I=1,IWA(4)
          IHELP=GETIWA(VIC,I,0,IWA,LIWA,INFOLI)
          IF (X .EQ. IHELP) THEN
            INSERT=I
            RETURN
          ENDIF
 10     CONTINUE
        CALL PUT1(X,0.0D0,VIC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) RETURN
        INSERT=IWA(4)
        RETURN
 11     IERR=22
        RETURN
C
      ELSE IF (TYPE .EQ. RNUM) THEN
        READ(COPY,'(D40.17)',ERR=21) Y
        DO 20 I=1,IWA(6)
          RHELP=GETWA(VRC,I,WA,LWA,INFOLI)
          IF (Y .EQ. RHELP) THEN
            INSERT=I
            RETURN
          ENDIF
 20     CONTINUE
        CALL PUT1(0,Y,VRC,0,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,IERR)
        IF (IERR .NE. 0) RETURN
        INSERT=IWA(6)
        RETURN
 21     IERR=23
        RETURN
C
      ELSE IF (TYPE .EQ. ID) THEN
        DO 30 I=1,SYMEND
          IF (STRING(1:STRLEN) .EQ. SYMNAM(I)) THEN
            INSERT=I
            RETURN
          ENDIF
 30     CONTINUE
        IF (SYMEND .GE. MAXSYM) THEN
          IERR=32
          RETURN
        ENDIF
        SYMEND=SYMEND+1
        SYMNAM(SYMEND)=STRING(1:STRLEN)
        SYMTYP(SYMEND)=0
        SYMREF(SYMEND)=0
        INSERT=SYMEND
        RETURN
C
      ELSE
        IERR=26
        WRITE(*,*) ' INSERT (5787) : unknown type ',TYPE
      ENDIF
      END
C
C
      INTEGER FUNCTION YYLEX (FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF,
     1                        YYLVAL,LIWA,PLIWA,IWA,LWA,PLWA,WA,INFOLI,
     2                        MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR,
     3                        DBGFIL,DBGLVL)
      INTEGER FD,LNUM,LPOS,SPOS
      CHARACTER STAR*5,CONT*1,LINE*66
      LOGICAL EOF
      INTEGER YYLVAL
      INTEGER LIWA,PLIWA,LWA,PLWA
      INTEGER IWA(LIWA),INFOLI(15)
      DOUBLE PRECISION WA(LWA)
      INTEGER MAXSYM
      CHARACTER*20 SYMNAM(MAXSYM)
      INTEGER SYMTYP(MAXSYM),SYMREF(MAXSYM),SYMEND
      INTEGER IERR,DBGFIL,DBGLVL
C
      INTEGER ADD,SUB,MULT,DIV,POWER,LEFT,RIGHT,COMMA,ASSIGN,NLINE
      INTEGER RANGE,RELOP,AND,OR,NOT,INUM,RNUM,ID,SUM,PROD,IN,IF,THEN
      INTEGER ELSE,ENDIF,STDRD,EXTERN,INTERP,PARAM,INDEX,REAL,INT,TABLE
      INTEGER CONINT,LININT,SPLINE,VAR,INFUNC,FUNC,END,GOTO,MARKE,CONTIN
      INTEGER UMINUS,INDVAR,ENDSUM,ENDPRD,BEQ,BRA,LABEL,VECTOR,ACTIVE
      PARAMETER (ADD=43,SUB=45,MULT=42,DIV=47,POWER=94,LEFT=40,RIGHT=41)
      PARAMETER (COMMA=44,ASSIGN=61,NLINE=10,RANGE=257,RELOP=258)
      PARAMETER (AND=259,OR=260,NOT=261,INUM=262,RNUM=263,ID=264)
      PARAMETER (SUM=265,PROD=266,IN=267,IF=268,THEN=269,ELSE=270)
      PARAMETER (ENDIF=271,STDRD=272,EXTERN=273,INTERP=274,PARAM=275)
      PARAMETER (INDEX=276,REAL=277,INT=278,TABLE=279,CONINT=280)
      PARAMETER (LININT=281,SPLINE=282,VAR=283,INFUNC=284,FUNC=285)
      PARAMETER (END=286,GOTO=287,MARKE=288,CONTIN=289,UMINUS=290)
      PARAMETER (INDVAR=291,ENDSUM=292,ENDPRD=293,BEQ=294,BRA=295)
      PARAMETER (LABEL=296,VECTOR=297,ACTIVE=298)
C
      INTEGER       TYPNUM
      PARAMETER (TYPNUM=12)
      CHARACTER*24  TYPNAM(TYPNUM)
      INTEGER       TYPLEN(TYPNUM)
C
      INTEGER      I,INSERT,J,KEYWD,STRLEN,K
      LOGICAL      BUCHST,ZIFFER,DEBUG
      CHARACTER    C,C2,GETCHR,STRING*40
C
      BUCHST(C) = (('A' .LE. C) .AND. (C .LE. 'Z') .OR. (C .EQ. '_'))
      ZIFFER(C) = (('0' .LE. C) .AND. (C .LE. '9'))
C
      DATA TYPNAM /'PARAMETER','SETOFINDICES','REALCONSTANT',
     1             'INTEGERCONSTANT','TABLE','CONINT','LININT',
     2             'SPLINE','VARIABLE','INDEX','FUNCTION','END'/
      DATA TYPLEN /9,12,12,15,5,6,6,6,8,5,8,3/
C
      DEBUG=((DBGLVL .EQ. 1) .OR. (DBGLVL .EQ. 3))
      STRLEN=0
      IF ((EOF) .OR. (IERR .NE. 0)) THEN
        YYLEX=-1
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: End of file'
        RETURN
      ENDIF
 10   C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
      IF (C .EQ. ' ') GO TO 10
 20   IF (C .EQ. '(') THEN
        YYLEX=LEFT
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: "("'
        RETURN
      ENDIF
 30   IF (C .EQ. ')') THEN
        YYLEX=RIGHT
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: ")"'
        RETURN
      ENDIF
 40   IF (C .EQ. '+') THEN
        YYLEX=ADD
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: "+"'
        RETURN
      ENDIF
 50   IF (C .EQ. '-') THEN
        YYLEX=SUB
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: "-"'
        RETURN
      ENDIF
 60   IF (C .EQ. '*') THEN
        C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
        IF (C .EQ. '*') THEN
          YYLEX=POWER
          YYLVAL=0
          IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: "**"'
          RETURN
        ELSE
          YYLEX=MULT
          YYLVAL=0
          IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: "*"'
          CALL UNGTC(LPOS,IERR)
          RETURN
        ENDIF
      ENDIF
 70   IF (C .EQ. '/') THEN
        YYLEX=DIV
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: "/"'
        RETURN
      ENDIF
 80   IF (C .EQ. '=') THEN
        YYLEX=ASSIGN
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: "="'
        RETURN
      ENDIF
 90   IF (C .EQ. ',') THEN
        YYLEX=COMMA
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: ","'
        RETURN
      ENDIF
 100  IF (C .NE. '.') GO TO 110
      C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
      IF (C .EQ. '.') THEN
        YYLEX=RANGE
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: ".."'
        YYLVAL=0
        RETURN
      ELSE IF (ZIFFER(C)) THEN
        STRING='.'
        STRLEN=1
        GO TO 113
      ELSE IF (.NOT. BUCHST(C)) THEN
        YYLEX=-1
        IERR=13
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: Range operator expected!'
        RETURN
      ENDIF
      STRING='.'//C
      STRLEN=2
 101  C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
      IF (BUCHST(C)) THEN
        STRLEN=STRLEN+1
        STRING(STRLEN:STRLEN)=C
        GO TO 101
      ELSE IF (C .EQ. '.') THEN
        STRLEN=STRLEN+1
        STRING(STRLEN:STRLEN)=C
      ENDIF
      YYLEX=KEYWD(STRING,YYLVAL,IERR)
      IF (IERR .NE. 0) THEN
        YYLEX=-1
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: Keyword expected!'
      ENDIF
      IF (DEBUG) THEN
      IF (YYLEX .EQ. RELOP) THEN
      WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Relational operator "',
     1                        STRING(1:STRLEN),'"'
      ELSE IF (YYLEX .EQ. AND) THEN
      WRITE(DBGFIL,'(A)') 'YYLEX: ".AND."'
      ELSE IF (YYLEX .EQ. OR) THEN
      WRITE(DBGFIL,'(A)') 'YYLEX: ".OR."'
      ELSE IF (YYLEX .EQ. NOT) THEN
      WRITE(DBGFIL,'(A)') 'YYLEX: ".NOT."'
      ELSE 
      WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Unknown keyword "',
     1                        STRING(1:STRLEN),'"'
      ENDIF
      ENDIF
      RETURN
 110  IF (.NOT. ZIFFER(C)) GO TO 120
      STRING=C
      STRLEN=1
      IF (LPOS .EQ. 1) THEN
 111    C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
        IF ((ZIFFER(C)) .AND. (SPOS .LE. 6)) THEN
          STRLEN=STRLEN+1
          STRING(STRLEN:STRLEN)=C
          GO TO 111
        ELSE
          YYLEX=MARKE
          YYLVAL=INSERT(STRING,STRLEN,INUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                  INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
          IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Label "',
     1                                       STRING(1:STRLEN),'"'
          IF ((SPOS .EQ. 6) .AND. (C .EQ. ' ')) LPOS=LPOS+1
          IF (SPOS .EQ. 6) THEN
            CALL UNGTC(LPOS,IERR)
          ELSE
            CALL UNGTC(SPOS,IERR)
          ENDIF
          IF (CONT .NE. ' ') IERR=31
          RETURN
        ENDIF
      ELSE
 112    C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
        IF (ZIFFER(C)) THEN
          STRLEN=STRLEN+1
          STRING(STRLEN:STRLEN)=C
          GO TO 112
        ELSE IF (C .EQ. 'D' .OR. C .EQ. 'E') THEN
          GO TO 114
        ELSE IF (C .NE. '.') THEN 
          YYLEX=INUM
          YYLVAL=INSERT(STRING,STRLEN,INUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                  INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
          IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Integer number "',
     1                                       STRING(1:STRLEN),'"'
          CALL UNGTC(LPOS,IERR)
          RETURN
        ENDIF
      ENDIF
      STRLEN=STRLEN+1
      STRING(STRLEN:STRLEN)='.'
      C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
      IF (ZIFFER(C)) THEN
        GO TO 113
      ELSE IF (BUCHST(C)) THEN
        IF (C .EQ. 'D' .OR. C .EQ. 'E') THEN
          GO TO 114
        ELSE
          STRING(STRLEN:STRLEN)=' '
          STRLEN=STRLEN-1
          YYLEX=INUM
          YYLVAL=INSERT(STRING,STRLEN,INUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                  INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
          IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Integer number "',
     1                                       STRING(1:STRLEN),'"'
          CALL UNGTC(LPOS,IERR)
          CALL UNGTC(LPOS,IERR)
          RETURN
        ENDIF
      ELSE IF (C .EQ. '.') THEN
        YYLEX=INUM
        STRING(STRLEN:STRLEN)=' '
        STRLEN=STRLEN-1
        YYLVAL=INSERT(STRING,STRLEN,INUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
        IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Integer number "',
     1                     STRING(1:STRLEN),'"'
        CALL UNGTC(LPOS,IERR)
        CALL UNGTC(LPOS,IERR)
        RETURN
      ELSE
        YYLEX=RNUM
        YYLVAL=INSERT(STRING,STRLEN,RNUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
        IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Real number "',
     1                     STRING(1:STRLEN),'"'
        CALL UNGTC(LPOS,IERR)
        RETURN
      ENDIF
 113  IF (ZIFFER(C)) THEN
        STRLEN=STRLEN+1
        STRING(STRLEN:STRLEN)=C
        C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
        GO TO 113
      ELSE IF (C .NE. 'D' .AND. C .NE. 'E') THEN
        STRLEN=STRLEN+2
        STRING(STRLEN-1:STRLEN)='D0'
        YYLEX=RNUM 
        YYLVAL=INSERT(STRING,STRLEN,RNUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
        IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Real number "',
     1                     STRING(1:STRLEN),'"'
        CALL UNGTC(LPOS,IERR)
        RETURN
      ENDIF
 114  STRLEN=STRLEN+1
      STRING(STRLEN:STRLEN)='D'
      C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
      IF (BUCHST(C)) THEN
        CALL UNGTC(LPOS,IERR)
        CALL UNGTC(LPOS,IERR)
        STRING(STRLEN:STRLEN)=' '
        STRLEN=STRLEN-1
        YYLEX=RNUM
        YYLVAL=INSERT(STRING,STRLEN,RNUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
        IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Real number "',
     1                     STRING(1:STRLEN),'"'
        RETURN
      ELSE IF (ZIFFER(C)) THEN
        GO TO 115
      ELSE IF (C .NE. '+' .AND. C .NE. '-') THEN
        YYLEX=-1
        IERR=23
        IF (DEBUG) WRITE(DBGFIL,'(A)') 
     1    'YYLEX: Wrong format for real number!'
        RETURN
      ENDIF
      STRLEN=STRLEN+1
      STRING(STRLEN:STRLEN)=C
      C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
 115  IF (ZIFFER(C)) THEN
        STRLEN=STRLEN+1
        STRING(STRLEN:STRLEN)=C
        C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
        GO TO 115
      ENDIF
      YYLEX=RNUM
      YYLVAL=INSERT(STRING,STRLEN,RNUM,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1              INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
      IF (DEBUG) WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Real number "',
     1                   STRING(1:STRLEN),'"'
      CALL UNGTC(LPOS,IERR)
      RETURN
 120  IF (BUCHST(C)) THEN
        STRING=C
        STRLEN=1
 121    C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
        IF (BUCHST(C) .OR. ZIFFER(C)) THEN
          STRLEN=STRLEN+1
          STRING(STRLEN:STRLEN)=C
          GO TO 121
        ENDIF
        CALL UNGTC(LPOS,IERR)
        YYLEX=KEYWD(STRING,YYLVAL,IERR)
        IF (DEBUG) THEN
        IF (YYLEX .EQ. SUM) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "SUM"'
        ELSE IF (YYLEX .EQ. PROD) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "PROD"'
        ELSE IF (YYLEX .EQ. IN) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "IN"'
        ELSE IF (YYLEX .EQ. IF) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "IF"'
        ELSE IF (YYLEX .EQ. THEN) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "THEN"'
        ELSE IF (YYLEX .EQ. ELSE) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "ELSE"'
        ELSE IF (YYLEX .EQ. ENDIF) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "ENDIF"'
        ELSE IF (YYLEX .EQ. GOTO) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "GOTO"'
        ELSE IF (YYLEX .EQ. CONTIN) THEN
        WRITE(DBGFIL,'(A)') 'YYLEX: "CONTINUE"'
        ELSE IF (YYLEX .EQ. STDRD) THEN
        WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Standard function "',
     1                      STRING(1:STRLEN),'"'
        ELSE IF (YYLEX .EQ. EXTERN) THEN
        WRITE(DBGFIL,'(A,A,A)') 'YYLEX: External function "',
     1                      STRING(1:STRLEN),'"'
        ELSE IF (YYLEX .EQ. ID) THEN
        WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Identifier "',
     1                      STRING(1:STRLEN),'"'
        ELSE 
        WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Unknown keyword "',
     1                      STRING(1:STRLEN),'"'
        ENDIF
        ENDIF
        IF (STRLEN .GT. 20) THEN
          YYLEX=-1
          IERR=27
          RETURN
        ENDIF
        IF (YYLEX .EQ. ID) THEN
          YYLVAL=INSERT(STRING,STRLEN,ID,LIWA,PLIWA,IWA,LWA,PLWA,WA,
     1                  INFOLI,MAXSYM,SYMNAM,SYMTYP,SYMREF,SYMEND,IERR)
          IF ( (SYMTYP(YYLVAL) .EQ. CONINT) .OR. 
     1         (SYMTYP(YYLVAL) .EQ. LININT) .OR.
     2         (SYMTYP(YYLVAL) .EQ. SPLINE) ) THEN
            YYLEX=INTERP
          ENDIF
        ENDIF
        RETURN
      ENDIF
 130  IF (C .EQ. '@') THEN
 131    C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
        IF (C .EQ. ' ') GO TO 131
        J=0
        DO 132 I=1,TYPNUM
          IF (C .EQ. TYPNAM(I)(1:1)) J=I
 132    CONTINUE
        K=0
        IF (C .EQ. 'I') THEN
 160      C2=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
          K=K+1
          IF (C2 .EQ. ' ') GO TO 160
 170      C2=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
          K=K+1
          IF (C2 .EQ. ' ') GO TO 170
          IF (C2 .EQ. 'T') J=4
          DO 180 I=1,K
            CALL UNGTC(LPOS,IERR)
 180      CONTINUE
        ENDIF
        IF (C .EQ. 'S') THEN
 190      C2=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
          K=K+1
          IF (C2 .EQ. ' ') GO TO 190
          IF (C2 .EQ. 'E') J=2
          DO 200 I=1,K
            CALL UNGTC(LPOS,IERR)
 200      CONTINUE
        ENDIF
        STRING=C
        DO 134 I=2,TYPLEN(J)
 133      C=GETCHR(FD,LNUM,STAR,CONT,LINE,LPOS,SPOS,EOF)
          IF (C .EQ. ' ') THEN
            IF ((STRING(1:3).EQ.'SET').OR.(STRING(1:3).EQ.'REA')
     /         .OR.(STRING(1:3).EQ.'INT')) LPOS=100
            GOTO 135
          ENDIF
          STRING(I:I)=C
 134    CONTINUE
 135    IF (STRING(1:3) .EQ. TYPNAM(J)(1:3)) THEN
          YYLEX=PARAM+(J-1)
          YYLVAL=0
          IF (DEBUG) THEN
          IF (j .EQ. 1) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "PARAMETER"'
          ELSE IF (j .EQ. 2) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "SET OF INDICES"'
          ELSE IF (j .EQ. 3) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "REAL CONSTANT"'
          ELSE IF (j .EQ. 4) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "INTEGER CONSTANT"'
          ELSE IF (j .EQ. 5) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "TABLE"'
          ELSE IF (j .EQ. 6) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "CONINT"'
          ELSE IF (j .EQ. 7) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "LININT"'
          ELSE IF (j .EQ. 8) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "SPLINE"'
          ELSE IF (j .EQ. 9) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "VARIABLE"'
          ELSE IF (j .EQ. 10) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "INDEX"'
          ELSE IF (j .EQ. 11) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "FUNCTION"'
          ELSE IF (j .EQ. 12) THEN 
          WRITE(DBGFIL,'(A)') 'YYLEX: Block "END"'
          ELSE 
          WRITE(DBGFIL,'(A,A,A)') 'YYLEX: Unknown block "',
     1                       STRING(1:STRLEN),'"'
          ENDIF
          ENDIF
          RETURN
        ELSE
          YYLEX=-1
          YYLVAL=0
          IERR=28
          RETURN
        ENDIF
      ENDIF
 140  IF (C .EQ. '#') THEN
        YYLEX=NLINE
        YYLVAL=0
        IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: Newline'
        RETURN
      ENDIF
 150  YYLEX=-1
      IF (DEBUG) WRITE(DBGFIL,'(A)') 'YYLEX: Unknown input character'
      IERR=29
      RETURN
      END
C
C
C
      SUBROUTINE PUT1(IVEK,RVEK,ART,OFFSET,LIWA,PLIWA,IWA,LWA,PLWA,
     1               WA,INFOLI,IERR)
C
      INTEGER ART,LIWA,PLIWA,LWA,PLWA,OFFSET
      INTEGER IVEK,IWA(LIWA),INFOLI(15),IERR
      DOUBLE PRECISION RVEK,WA(LWA)
C
      INTEGER VIS,VIC,VRC,VVA,VFN,VGR,VHE,VPF,IV
C     INTEGER LVIS,LIIC,LVIC,LIRC,LVRC,LVVA,LVFN,LVGR,LVHE,LVPF,LIV
      INTEGER LVIS,LIIC,LVIC,LIRC,LVRC,LVVA,LVFN,LVGR,LVPF,LIV
      PARAMETER (VIS=2,VIC=4,VRC=6,VVA=8,VFN=11,VGR=12,VHE=13)
      PARAMETER (VPF=14,IV=15)
C
      INTEGER I
C
      LVIS=INFOLI(2)
      LIIC=INFOLI(3)
      LVIC=INFOLI(4)
      LIRC=INFOLI(5)
      LVPF=INFOLI(14)
      LIV =INFOLI(15)
C
      LVRC=INFOLI(6)
      LVVA=INFOLI(8)
      LVFN=INFOLI(11)
      LVGR=INFOLI(12)
C     LVHE=INFOLI(13)
C
      IF (ART .EQ. VIS) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(2)=IWA(2)+1
          PLIWA=PLIWA+1
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-1 .GT. LIIC) THEN
            DO 100 I=PLIWA,LIIC+1+1,-1
              IWA(I)=IWA(I-1)
 100         CONTINUE
          ENDIF
          IWA(IWA(2)+LVIS)=IVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(2))) THEN
          IWA(LVIS+OFFSET)=IVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. VIC) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(4)=IWA(4)+1
          PLIWA=PLIWA+1
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-1 .GT. LIRC) THEN
            DO 300 I=PLIWA,LIRC+1+1,-1
              IWA(I)=IWA(I-1)
 300         CONTINUE
          ENDIF
          IWA(IWA(4)+LVIC)=IVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(4))) THEN
          IWA(LVIC+OFFSET)=IVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. VRC) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(6)=IWA(6)+1
          PLWA=PLWA+1
          IF (PLWA .GT. LWA) GO TO 9999
          IF (PLWA-1 .GT. LVVA) THEN
            DO 500 I=PLWA,LVVA+1+1,-1
              WA(I)=WA(I-1)
 500         CONTINUE
          ENDIF
          WA(IWA(6)+LVRC)=RVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(6))) THEN
          WA(LVRC+OFFSET)=RVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. VVA) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(8)=IWA(8)+1
          PLWA=PLWA+1
          IF (PLWA .GT. LWA) GO TO 9999
          IF (PLWA-1 .GT. LVFN) THEN
            DO 700 I=PLWA,LVFN+1+1,-1
              WA(I)=WA(I-1)
 700         CONTINUE
          ENDIF
          WA(IWA(8)+LVVA)=RVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(8))) THEN
          WA(LVVA+OFFSET)=RVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. VFN) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(11)=IWA(11)+1
          PLWA=PLWA+1
          IF (PLWA .GT. LWA) GO TO 9999
          IF (PLWA-1 .GT. LVGR) THEN
            DO 900 I=PLWA,LVGR+1+1,-1
              WA(I)=WA(I-1)
 900         CONTINUE
          ENDIF
          WA(IWA(11)+LVFN)=RVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(11))) THEN
          WA(LVFN+OFFSET)=RVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. VGR) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(12)=IWA(12)+1
          PLWA=PLWA+1
          IF (PLWA .GT. LWA) GO TO 9999
          IF (PLWA-1 .GT. IWA(12)-1+LVGR) THEN
            DO 1000 I=PLWA,IWA(12)-1+LVGR+1+1,-1
              WA(I)=WA(I-1)
 1000       CONTINUE
          ENDIF
          WA(IWA(12)+LVGR)=RVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(12))) THEN
          WA(LVGR+OFFSET)=RVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. VPF) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(14)=IWA(14)+1
          PLIWA=PLIWA+1
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-1 .GT. LIV) THEN
            DO 1100 I=PLIWA,LIV+1+1,-1
              IWA(I)=IWA(I-1)
 1100       CONTINUE
          ENDIF
          IWA(LVPF+IWA(14))=IVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(14))) THEN
          IWA(LVPF+OFFSET)=IVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. IV) THEN
        IF (OFFSET .EQ. 0) THEN
          IWA(15)=IWA(15)+1
          PLIWA=PLIWA+1
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-1 .GT. IWA(15)-1+LIV) THEN
            DO 1200 I=PLIWA,IWA(15)-1+LIV+1+1,-1
              IWA(I)=IWA(I-1)
 1200       CONTINUE
          ENDIF
          IWA(LIV+IWA(15))=IVEK
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(15))) THEN
          IWA(LIV+OFFSET)=IVEK
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE
        IERR=34
        RETURN
C
C
      ENDIF
C
      INFOLI(1)=15
      INFOLI(2)=IWA(1)*5+INFOLI(1)
      INFOLI(3)=IWA(2)+INFOLI(2)
      INFOLI(4)=IWA(3)*4+INFOLI(3)
      INFOLI(5)=IWA(4)+INFOLI(4)
      INFOLI(7)=IWA(5)*4+INFOLI(5)
      INFOLI(9)=IWA(7)*3+INFOLI(7)
      INFOLI(14)=IWA(9)*7+INFOLI(9)
      INFOLI(15)=IWA(14)+INFOLI(14)
C
      INFOLI(6)=0
      INFOLI(8)=IWA(6)+INFOLI(6)
      INFOLI(11)=IWA(8)+INFOLI(8)
      INFOLI(12)=IWA(11)+INFOLI(11)
      INFOLI(13)=IWA(12)*IWA(8)+INFOLI(12)
C 
      RETURN
C
 9999 CONTINUE
C     WRITE(*,*) 'OUT OF MEMORY'
      IERR=32
      RETURN
      END
C
C
C
      SUBROUTINE PUT3(IVEK,ART,OFFSET,LIWA,PLIWA,IWA,INFOLI,IERR)
C
      INTEGER ART,LIWA,PLIWA,OFFSET
      INTEGER IVEK(3),IWA(LIWA),INFOLI(15),IERR
C
      INTEGER IVA
      INTEGER LIVA,LIFN
      PARAMETER (IVA=7)
C
      INTEGER I,J
C
      LIVA=INFOLI(7)
      LIFN=INFOLI(9)
C
      IF (ART .EQ. IVA) THEN
C
        IF (OFFSET .EQ. 0) THEN
          IWA(7)=IWA(7)+1
          PLIWA=PLIWA+3
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-3 .GT. LIFN) THEN
            DO 100 I=PLIWA,LIFN+3+1,-1
              IWA(I)=IWA(I-3)
 100        CONTINUE
          ENDIF
          IF (IWA(7) .GT. 1) THEN
            DO 120 I=3,1,-1
              IWA(LIVA+IWA(7)*I)=IVEK(I)
              DO 110 J=1,IWA(7)-1
                IWA(LIVA+IWA(7)*I-J)=IWA(LIVA+IWA(7)*I-J-(I-1))
 110           CONTINUE
 120         CONTINUE
          ELSE
            DO 130 I=1,3
              IWA(LIVA+I) = IVEK(I)
 130         CONTINUE
          ENDIF
C
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(7))) THEN
          DO 140 I=1,3
            IWA(LIVA+IWA(7)*(I-1)+OFFSET)=IVEK(I)
 140      CONTINUE
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE
        IERR=34
        RETURN
C
C
      ENDIF
C
      INFOLI(1)=15
      INFOLI(2)=IWA(1)*5+INFOLI(1)
      INFOLI(3)=IWA(2)+INFOLI(2)
      INFOLI(4)=IWA(3)*4+INFOLI(3)
      INFOLI(5)=IWA(4)+INFOLI(4)
      INFOLI(7)=IWA(5)*4+INFOLI(5)
      INFOLI(9)=IWA(7)*3+INFOLI(7)
      INFOLI(14)=IWA(9)*7+INFOLI(9)
      INFOLI(15)=IWA(14)+INFOLI(14)
C
      INFOLI(6)=0
      INFOLI(8)=IWA(6)+INFOLI(6)
      INFOLI(11)=IWA(8)+INFOLI(8)
      INFOLI(12)=IWA(11)+INFOLI(11)
      INFOLI(13)=IWA(12)*IWA(8)+INFOLI(12)
C 
      RETURN
C
 9999 CONTINUE
C     WRITE(*,*) 'OUT OF MEMORY'
      IERR=32
      RETURN
      END
C
C
C
      SUBROUTINE PUT4(IVEK,ART,OFFSET,LIWA,PLIWA,IWA,INFOLI,IERR)
C
      INTEGER ART,LIWA,PLIWA,OFFSET
      INTEGER IVEK(4),IWA(LIWA),INFOLI(15),IERR
C
      INTEGER IIC,IRC
      INTEGER LIIC,LVIC,LIRC,LIVA
      PARAMETER (IIC=3,IRC=5)
C
      INTEGER I,J
C
      LIIC=INFOLI(3)
      LVIC=INFOLI(4)
      LIRC=INFOLI(5)
      LIVA=INFOLI(7)
C
      IF (ART .EQ. IIC) THEN
C
        IF (OFFSET .EQ. 0) THEN
          IWA(3)=IWA(3)+1
          PLIWA=PLIWA+4
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-4 .GT. LVIC) THEN
            DO 100 I=PLIWA,LVIC+1+4,-1
              IWA(I)=IWA(I-4)
 100        CONTINUE
          ENDIF
          IF (IWA(3) .GT. 1) THEN
            DO 120 I=4,1,-1
              IWA(LIIC+IWA(3)*I)=IVEK(I)
              DO 110 J=1,IWA(3)-1
                IWA(LIIC+IWA(3)*I-J)=IWA(LIIC+IWA(3)*I-J-(I-1))
 110          CONTINUE
 120        CONTINUE
          ELSE
            DO 130 I=1,4
              IWA(LIIC+I) = IVEK(I)
 130        CONTINUE
          ENDIF
C
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(3))) THEN
          DO 140 I=1,4
            IWA(LIIC+IWA(3)*(I-1)+OFFSET)=IVEK(I)
 140      CONTINUE
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE IF (ART .EQ. IRC) THEN
C
        IF (OFFSET .EQ. 0) THEN
          IWA(5)=IWA(5)+1
          PLIWA=PLIWA+4
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-4 .GT. LIVA) THEN
            DO 200 I=PLIWA,LIVA+4+1,-1
              IWA(I)=IWA(I-4)
 200        CONTINUE
          ENDIF
          IF (IWA(5) .GT. 1) THEN
            DO 220 I=4,1,-1
              IWA(LIRC+IWA(5)*I)=IVEK(I)
              DO 210 J=1,IWA(5)-1
                IWA(LIRC+IWA(5)*I-J)=IWA(LIRC+IWA(5)*I-J-(I-1))
 210           CONTINUE
 220         CONTINUE
          ELSE
            DO 230 I=1,4
              IWA(LIRC+I) = IVEK(I)
 230         CONTINUE
          ENDIF
C
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(5))) THEN
          DO 240 I=1,4
            IWA(LIRC+IWA(5)*(I-1)+OFFSET)=IVEK(I)
 240      CONTINUE
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE
        IERR=34
        RETURN
C
C
      ENDIF
C
      INFOLI(1)=15
      INFOLI(2)=IWA(1)*5+INFOLI(1)
      INFOLI(3)=IWA(2)+INFOLI(2)
      INFOLI(4)=IWA(3)*4+INFOLI(3)
      INFOLI(5)=IWA(4)+INFOLI(4)
      INFOLI(7)=IWA(5)*4+INFOLI(5)
      INFOLI(9)=IWA(7)*3+INFOLI(7)
      INFOLI(14)=IWA(9)*7+INFOLI(9)
      INFOLI(15)=IWA(14)+INFOLI(14)
C
      INFOLI(6)=0
      INFOLI(8)=IWA(6)+INFOLI(6)
      INFOLI(11)=IWA(8)+INFOLI(8)
      INFOLI(12)=IWA(11)+INFOLI(11)
      INFOLI(13)=IWA(12)*IWA(8)+INFOLI(12)
C 
      RETURN
C
 9999 CONTINUE
C     WRITE(*,*) 'OUT OF MEMORY'
      IERR=32
      RETURN
      END
C
C
C
      SUBROUTINE PUT5(IVEK,ART,OFFSET,LIWA,PLIWA,IWA,INFOLI,IERR)
C
      INTEGER ART,LIWA,PLIWA,OFFSET
      INTEGER IVEK(5),IWA(LIWA),INFOLI(15),IERR
C
      INTEGER IIS
      INTEGER LIIS,LVIS
      PARAMETER (IIS=1)
C
      INTEGER I,J
C
      LIIS=INFOLI(1)
      LVIS=INFOLI(2)
C
      IF (ART .EQ. IIS) THEN
C
        IF (OFFSET .EQ. 0) THEN
          IWA(1)=IWA(1)+1
          PLIWA=PLIWA+5
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-5 .GT. LVIS) THEN
            DO 10 I=PLIWA,LVIS+5+1,-1
              IWA(I)=IWA(I-5)
 10         CONTINUE
          ENDIF
          IF (IWA(1) .GT. 1) THEN
            DO 30 I=5,1,-1
              IWA(LIIS+IWA(1)*I)=IVEK(I)
              DO 20 J=1,IWA(1)-1
                IWA(LIIS+IWA(1)*I-J)=IWA(LIIS+IWA(1)*I-J-(I-1))
 20           CONTINUE
 30         CONTINUE
          ELSE
            DO 40 I=1,5
              IWA(LIIS+I) = IVEK(I)
 40         CONTINUE
          ENDIF
C
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(1))) THEN
          DO 140 I=1,5
            IWA(LIIS+IWA(1)*(I-1)+OFFSET)=IVEK(I)
 140      CONTINUE
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE
        IERR=34
        RETURN
C
C
      ENDIF
C
      INFOLI(1)=15
      INFOLI(2)=IWA(1)*5+INFOLI(1)
      INFOLI(3)=IWA(2)+INFOLI(2)
      INFOLI(4)=IWA(3)*4+INFOLI(3)
      INFOLI(5)=IWA(4)+INFOLI(4)
      INFOLI(7)=IWA(5)*4+INFOLI(5)
      INFOLI(9)=IWA(7)*3+INFOLI(7)
      INFOLI(14)=IWA(9)*7+INFOLI(9)
      INFOLI(15)=IWA(14)+INFOLI(14)
C
      INFOLI(6)=0
      INFOLI(8)=IWA(6)+INFOLI(6)
      INFOLI(11)=IWA(8)+INFOLI(8)
      INFOLI(12)=IWA(11)+INFOLI(11)
      INFOLI(13)=IWA(12)*IWA(8)+INFOLI(12)
C 
      RETURN
C
 9999 CONTINUE
C     WRITE(*,*) 'OUT OF MEMORY'
      IERR=32
      RETURN
      END
C
C
C
      SUBROUTINE PUT7(IVEK,ART,OFFSET,LIWA,PLIWA,IWA,INFOLI,IERR)
C
      INTEGER ART,LIWA,PLIWA,OFFSET
      INTEGER IVEK(7),IWA(LIWA),INFOLI(15),IERR
C
      INTEGER IFN
      INTEGER LIFN,LVPF
      PARAMETER (IFN=9)
C
      INTEGER I,J
C
      LIFN=INFOLI(9)
      LVPF=INFOLI(14)
C
      IF (ART .EQ. IFN) THEN
C
        IF (OFFSET .EQ. 0) THEN
          IWA(9)=IWA(9)+1
          PLIWA=PLIWA+7
          IF (PLIWA .GT. LIWA) GO TO 9999
          IF (PLIWA-7 .GT. LVPF) THEN
            DO 100 I=PLIWA,LVPF+7+1,-1
              IWA(I)=IWA(I-7)
 100        CONTINUE
          ENDIF
          IF (IWA(9) .GT. 1) THEN
            DO 120 I=7,1,-1
              IWA(LIFN+IWA(9)*I)=IVEK(I)
              DO 110 J=1,IWA(9)-1
                IWA(LIFN+IWA(9)*I-J)=IWA(LIFN+IWA(9)*I-J-(I-1))
 110           CONTINUE
 120         CONTINUE
          ELSE
            DO 130 I=1,7
              IWA(LIFN+I) = IVEK(I)
 130        CONTINUE
          ENDIF
C
        ELSE IF ((OFFSET .GE. 1) .AND. (OFFSET .LE. IWA(9))) THEN
          DO 140 I=1,7
            IWA(LIFN+IWA(9)*(I-1)+OFFSET)=IVEK(I)
 140      CONTINUE
        ELSE
          IERR=34
          RETURN
        ENDIF    
C
C
      ELSE
        IERR=34
        RETURN
C
C
      ENDIF
C
      INFOLI(1)=15
      INFOLI(2)=IWA(1)*5+INFOLI(1)
      INFOLI(3)=IWA(2)+INFOLI(2)
      INFOLI(4)=IWA(3)*4+INFOLI(3)
      INFOLI(5)=IWA(4)+INFOLI(4)
      INFOLI(7)=IWA(5)*4+INFOLI(5)
      INFOLI(9)=IWA(7)*3+INFOLI(7)
      INFOLI(14)=IWA(9)*7+INFOLI(9)
      INFOLI(15)=IWA(14)+INFOLI(14)
C
      INFOLI(6)=0
      INFOLI(8)=IWA(6)+INFOLI(6)
      INFOLI(11)=IWA(8)+INFOLI(8)
      INFOLI(12)=IWA(11)+INFOLI(11)
      INFOLI(13)=IWA(12)*IWA(8)+INFOLI(12)
C 
      RETURN
C
 9999 CONTINUE
C     WRITE(*,*) 'OUT OF MEMORY'
      IERR=32
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION GETWA(FELD,DIM,WA,LWA,INFOLI)
C
      INTEGER FELD,DIM
      INTEGER INFOLI(15),LWA
      DOUBLE PRECISION WA(LWA)
C
      INTEGER VRC,VVA,VFN,VGR,VHE
      INTEGER LVRC,LVVA,LVFN,LVGR,LVHE
C
      PARAMETER (VRC=6,VVA=8,VFN=11,VGR=12,VHE=13)
C
      LVRC=INFOLI(6)
      LVVA=INFOLI(8)
      LVFN=INFOLI(11)
      LVGR=INFOLI(12)
      LVHE=INFOLI(13)
C
      IF (FELD .EQ. VRC) THEN
        GETWA=WA(LVRC+DIM)
C
      ELSE IF (FELD .EQ. VVA) THEN
        GETWA=WA(LVVA+DIM)
C
      ELSE IF (FELD .EQ. VFN) THEN
        GETWA=WA(LVFN+DIM)
C
      ELSE IF (FELD .EQ. VGR) THEN
        GETWA=WA(LVGR+DIM)
C    
      ELSE IF (FELD .EQ. VHE) THEN
        GETWA=WA(LVHE+DIM)
C
      ELSE
        STOP 'GETWA: INTERNAL ERROR OF THE DYNAMIC PARSER.'
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE UNVPF(PC,LIWA,PLIWA,IWA,INFOLI)
C
      INTEGER PC,LIWA,PLIWA,IWA(LIWA),INFOLI(15)
C
      INTEGER I
C
      DO 100 I=1,IWA(15)
        IWA(INFOLI(14)+PC+I)=IWA(INFOLI(15)+I)
        IWA(INFOLI(15)+I)=0
 100  CONTINUE
      INFOLI(15)=INFOLI(14)+PC
      IWA(14)=PC
      PLIWA=INFOLI(15)+IWA(15)
      RETURN
      END 
C
C
C
      SUBROUTINE EVALCA(LIWA,IWA,LWA,WA,INFOLI,PC,IVAL,FVAL,IERR)
      INTEGER LWA,LIWA
      INTEGER IWA(LIWA)
      DOUBLE PRECISION WA(LWA)
      INTEGER INFOLI(15)
      INTEGER PC
      DOUBLE PRECISION FVAL
      INTEGER IVAL,IERR
C
      INTEGER MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN
      INTEGER MPVFN,MPVPF,MPIV,MPVHE,PVVAHE,MODE
      INTEGER DFVAR(1)
C
      INTEGER NOGRAD
      PARAMETER (NOGRAD=0)
C
        MPIIS=MAX(1,IWA(1))
        MPVIS=MAX(1,IWA(2))
        MPIIC=MAX(1,IWA(3))
        MPVIC=MAX(1,IWA(4))
        MPIRC=MAX(1,IWA(5))
        MPVRC=MAX(1,IWA(6))
        MPIVA=MAX(1,IWA(7))
        MPVVA=MAX(1,IWA(8))
        MPIFN=MAX(1,IWA(9))
        MPVFN=MAX(1,IWA(11))
        MPVPF=MAX(1,IWA(14))
        MPIV =MAX(1,IWA(15))
        MPVHE=1
        PVVAHE=1
        MODE=NOGRAD
        CALL EVAL(PC+1,IVAL,FVAL,MODE,IWA(8),1,MPIIS,MPVIS,MPIIC,
     1            MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN,MPVFN,1,MPVPF,
     2            MPIV,IWA(INFOLI(1)+1),IWA(INFOLI(2)+1),
     3            IWA(INFOLI(3)+1),IWA(INFOLI(4)+1),IWA(INFOLI(5)+1),
     4            WA(INFOLI(6)+1),IWA(INFOLI(7)+1),WA(INFOLI(8)+1),
     5            IWA(INFOLI(9)+1),WA(INFOLI(11)+1),WA(INFOLI(12)+1),
     6            WA(INFOLI(13)+1),IWA(INFOLI(14)+1),IWA(INFOLI(15)+1),
     7            WA(INFOLI(13)+MPVHE*PVVAHE*PVVAHE+1),
     8            WA(INFOLI(13)+MPVHE*PVVAHE*PVVAHE+MPVVA*10+1),
     9            DFVAR,0,IERR)
       RETURN
       END
C
C
C
      SUBROUTINE YYPRE (INPUT,OUTPUT,IERR,LNUM)
      INTEGER INPUT,OUTPUT,IERR,LNUM
C
      CHARACTER STAR*5,CONT*1,LINE*66
C
      INTEGER MAXMAC
      PARAMETER (MAXMAC=50)
      CHARACTER*20 MNAME(MAXMAC)
      INTEGER MSTART(MAXMAC),MLINES(MAXMAC),MFLAG(MAXMAC)
      INTEGER MCOUNT
C
      INTEGER MAXBUF
      PARAMETER (MAXBUF=1000)
      CHARACTER*72 EXPAND(MAXBUF),BUFFER(MAXBUF)
      INTEGER ECOUNT,BCOUNT
C
      CHARACTER*20 KEYWD
      INTEGER POS,LENGTH,I,J,K
C
C   INITIALIZE VARIABLES
C
      MCOUNT=0
      ECOUNT=0
      BCOUNT=0
      IERR=0
      LNUM=0
C
C   READ AND PROCESS ONE LINE FROM INPUT DEVICE
C
   10 READ(INPUT,'(A5,A1,A66)',END=90) STAR,CONT,LINE
      LNUM=LNUM+1
   12 IF (STAR(1:1).NE.'*') GOTO 80
      CALL GETID(LINE,1,POS,LENGTH)
      IF ((POS.EQ.0).OR.(LENGTH.GT.20)) GOTO 80
      KEYWD=LINE(POS:POS+LENGTH-1)
      K=POS+LENGTH
C
C   MACRO DEFINITION : READ THE NAME AND INSERT IT INTO THE TABLE
C
      IF ((KEYWD.EQ.'MACRO').OR.(KEYWD.EQ.'macro').OR.
     /       (KEYWD.EQ.'Macro')) THEN
        CALL GETID(LINE,K,POS,LENGTH)
        IF (POS.EQ.0) THEN
          IERR=68
          GOTO 90  
        ENDIF
        MCOUNT=MCOUNT+1
        IF (MCOUNT.GT.MAXMAC) THEN
          IERR=69
          GOTO 90
        ENDIF
        IF (LENGTH.GT.20) LENGTH=20
        MNAME(MCOUNT)=LINE(POS:POS+LENGTH-1)
        MSTART(MCOUNT)=ECOUNT+1
        MLINES(MCOUNT)=0
        MFLAG(MCOUNT)=0
C
C   ... READ THE MACRO EXPANSION (UP TO THE NEXT BLOCK DEFINITION)
C
   20   READ(INPUT,'(A5,A1,A66)',END=90) STAR,CONT,LINE
        LNUM=LNUM+1
        IF (STAR(1:1).NE.'*') THEN
          ECOUNT=ECOUNT+1
          IF (ECOUNT.GT.MAXBUF) THEN
            IERR=70
            GOTO 90
          ENDIF
          EXPAND(ECOUNT)=STAR//CONT//LINE        
          MLINES(MCOUNT)=MLINES(MCOUNT)+1
          GOTO 20
        ENDIF
        GOTO 12
C
C   FUNCTION : READ STATEMENTS INTO THE BUFFER
C
      ELSE IF (KEYWD.EQ.'FUNCTION') THEN
        WRITE(OUTPUT,'(3A)') STAR,CONT,LINE
        DO 42 I=1,MCOUNT
          MFLAG(I)=0
   42   CONTINUE
        BCOUNT=0
C
   40   READ(INPUT,'(A5,A1,A66)',END=90) STAR,CONT,LINE
        LNUM=LNUM+1
        IF (STAR(1:1).EQ.'C'.OR.LINE.EQ.' ') GOTO 44
C
C   NEXT BLOCK DEFINITION - WRITE STATEMENTS AND USED MACROS TO OUTPUT
C
        IF (STAR(1:1).EQ.'*') THEN
          DO 401 I=1,MCOUNT
          IF (MFLAG(I).NE.0) THEN
            DO 402 J=1,MLINES(I)
              WRITE(OUTPUT,'(A)') EXPAND(MSTART(I)+J-1)
  402       CONTINUE
          MFLAG(I)=0
          ENDIF
  401     CONTINUE
          DO 403 I=1,BCOUNT
            WRITE(OUTPUT,'(A)') BUFFER(I)
  403     CONTINUE
          BCOUNT=0
          GOTO 12
        ENDIF
C
C   SAVE THE CURRENT STATEMENT IN THE BUFFER
C
   44   BCOUNT=BCOUNT+1
        IF (BCOUNT.GT.MAXBUF) THEN
          IERR=71
          GOTO 90
        ENDIF
        BUFFER(BCOUNT)=STAR//CONT//LINE
C
C   LOOK FOR MACROS IN THE CURRENT STATEMENT (SET FLAGS IN THE TABLE)
C
        K=1
   45   CALL GETID(LINE,K,POS,LENGTH)
        IF (POS.GT.0) THEN
          K=POS+LENGTH
          IF (LENGTH.GT.20) LENGTH=20
          DO 47 I=1,MCOUNT
            IF (MNAME(I).EQ.LINE(POS:POS+LENGTH-1)) MFLAG(I)=1
   47     CONTINUE
          GOTO 45
        ENDIF 
        GOTO 40
      ENDIF
C
C   OTHER LINE - JUST PRINT IT
C
   80 WRITE(OUTPUT,'(3A)') STAR,CONT,LINE
      GOTO 10
   90 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE GETID (LINE,START,POS,LENGTH)
      CHARACTER*66 LINE
      INTEGER START,POS,LENGTH
C
      CHARACTER C
      INTEGER P,Q
C
      LOGICAL ISALPH,ISDIGI,ISALNU,ISLOWE
      CHARACTER TOUPPE
      ISALPH(C)=((C.GE.'A').AND.(C.LE.'Z').OR.(C.EQ.'_'))
      ISDIGI(C)=((C.GE.'0').AND.(C.LE.'9'))
      ISALNU(C)=(ISALPH(C).OR.ISDIGI(C))
      ISLOWE(C)=((C.GE.'a').AND.(C.LE.'z'))
      TOUPPE(C)=CHAR(ICHAR(C)+ICHAR('A')-ICHAR('a'))
C
      P=START
   10 C=LINE(P:P)
      IF (ISLOWE(C)) THEN
        C=TOUPPE(C)
        LINE(P:P)=C
      ENDIF
      IF (.NOT.ISALPH(C)) THEN
        P=P+1
        IF (P.LE.66) GOTO 10
      ENDIF
      IF (P.EQ.67) THEN
        POS=0
        LENGTH=0
        RETURN
      ENDIF
      POS=P
      Q=P
   20 C=LINE(Q:Q)
      IF (ISALNU(C)) THEN
        Q=Q+1
        IF (Q.LE.66) GOTO 20
      ENDIF
      LENGTH=Q-P
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NEXTID (NOUT)
      INTEGER NOUT
C
      LOGICAL OP
C
      NOUT=10
   10 INQUIRE(NOUT,OPENED=OP)
      IF (OP) THEN
	NOUT=NOUT+1
	GOTO 10
      ENDIF
      RETURN
      END
C*********************************************
C                                            *
C   PROGRAM   : PCOMP                        * 
C   MODULE    : S (RUNTIME EXECUTIVE)        *
C   ABSTRACT  : FORTRAN PRECOMPILER          *
C   KEY WORD  : AUTOMATIC DIFFERENTIATION    *
C   SOURCE    : PCOMP 2.3 by M.LIEPELT       *
C               PCOMP 3.0 by M.DOBMANN       *
C   COPYRIGHT : C.TRASSL, K.SCHITTKOWSKI     *
C               MATHEMATISCHES INSTITUT,     *
C               UNIVERSITAET BAYREUTH,       *
C               D-95440 BAYREUTH, GERMANY    *
C   DATE      : NOVEMBER 23, 1997            *
C   VERSION   : 5.3                          *
C                                            *
C*********************************************
C
C
C
C     SUBROUTINE SYMPRP (SYMFIL,WA,LWA,IWA,LIWA,UWA,UIWA,IERR,MODE,
C    /                   NVAR,NFUNC)
C     INTEGER SYMFIL
C     INTEGER LWA,LIWA
C     DOUBLE PRECISION WA(LWA)
C     INTEGER IWA(LIWA)
C     INTEGER UWA,UIWA
C     INTEGER IERR,MODE,NVAR,NFUNC
C
C**********************************************************************
C
C   S Y M P R P   -   LOAD INTERMEDIATE CODE GENERATED BY SYMINP FROM
C                     SYMFIL INTO WORKING ARRAYS.
C
C   PARAMETERS:
C      SYMFIL    - INPUT DEVICE; THE INTERMEDIATE CODE GENERATED BY
C                  SYMINP WAS WRITTEN TO THIS FILE AND IS NOW LOADED.
C      WA(LWA)   - REAL WORKING ARRAY, REQUIRED BY SYMPRP. ON RETURN,
C                  WA() CONTAINS THE INTERMEDIATE CODE.
C      IWA(LIWA) - INTEGER WORKING ARRAY, CF. WA().
C      UWA,UIWA  - INDICATE THE ACTUAL SPACE OF WA() AND IWA() THAT
C                  HAS BEEN USED BY THE SUBROUTINE.
C      IERR      - THE PARAMETER SHOWS THE REASON FOR TERMINATING THE
C                  SUBROUTINE. ON RETURN IERR COULD CONTAIN THE FOLLOW-
C                  ING VALUES:
C                  IERR = 0 : SUCCESSFUL TERMINATION.
C                  IERR > 0 : AN ERROR HAS BEEN DETECTED. FOR FURTHER
C                             INFORMATION CF. SUBROUTINE SYMERR.
C      MODE      - THE PARAMETER IS USED FOR THE RESERVATION OF SPACE
C                  FOR THE HESSIAN MATRIX
C      NVAR      - ON RETURN, NVAR CONTAINS THE NUMBER OF VARIABLES ON
C                  FUNCTION INPUT FILE
C      NFUNC     - ON RETURN, NFUNC CONTAINS THE NUMBER OF FUNCTIONS ON
C                  INPUT FILE
C
C**********************************************************************
C
C     INTEGER I,PWA,PIWA,PX,PIX
C
C     INTEGER GSMDEP,GETIWA,HILF
C     PARAMETER (GSMDEP=10)
C
C     INTEGER HSMDEP
C     PARAMETER (HSMDEP=10)
C
C     INTEGER NOGRAD,GRAD,HESS
C     PARAMETER (NOGRAD=0,GRAD=1,HESS=2)
C
C     INTEGER IIS,VIS,IIC,VIC,IRC,VRC,IVA,VVA,IFN,XFN,VFN,VGR,VHE,VPF,IV
C     PARAMETER (IIS=1,VIS=2,IIC=3,VIC=4,IRC=5,VRC=6,IVA=7,VVA=8)
C     PARAMETER (IFN=9,XFN=10,VFN=11,VGR=12,VHE=13,VPF=14,IV=15)
C     INTEGER INFOLI(15)
C
C     DO 10 I=1,15
C       READ(SYMFIL,'(I6)',ERR=100) IWA(I)
C10   CONTINUE
C     PIWA=IWA(1)*5+IWA(2)+IWA(3)*4+IWA(4)+IWA(5)*4+IWA(7)*3+
C    1     IWA(9)*7+IWA(14)+15
C     PWA=IWA(6)
C     IF (MODE .EQ. NOGRAD) THEN
C       PX=IWA(8)+IWA(11)+GSMDEP+HSMDEP+2
C     ELSE IF (MODE .EQ. GRAD) THEN
C       PX=IWA(8)+IWA(11)+IWA(12)*IWA(8)+GSMDEP*IWA(8)+HSMDEP+1
C     ELSE IF (MODE .EQ. HESS) THEN
C       PX=IWA(8)+IWA(11)+IWA(12)*IWA(8)+IWA(13)*IWA(8)*IWA(8)+
C    1     GSMDEP*IWA(8)+HSMDEP*IWA(8)*IWA(8)
C     ELSE
C       IERR=61
C     ENDIF
C     IF (MODE .EQ. NOGRAD) IWA(12)=0
C     IF (MODE .NE. HESS) IWA(13)=0
C     PIX=IWA(15)
C     IF ((PWA+PX .GT. LWA) .OR. (PIWA+PIX .GT. LIWA)) THEN
C       IERR=32
C       RETURN
C     ENDIF
C     DO 20 I=16,PIWA
C       READ(SYMFIL,'(I6)',ERR=100) IWA(I)
C20   CONTINUE
C
C     NVAR=IWA(8)
C     NFUNC=IWA(10)
C     INFOLI(1)=15
C     INFOLI(2)=IWA(1)*5+INFOLI(1)
C     INFOLI(3)=IWA(2)+INFOLI(2)
C     INFOLI(4)=IWA(3)*4+INFOLI(3)
C     INFOLI(5)=IWA(4)+INFOLI(4)
C     INFOLI(7)=IWA(5)*4+INFOLI(5)
C     INFOLI(9)=IWA(7)*3+INFOLI(7)
C     INFOLI(14)=IWA(9)*7+INFOLI(9)
C     INFOLI(15)=IWA(14)+INFOLI(14)
C
C     INFOLI(6)=0
C     INFOLI(8)=IWA(6)+INFOLI(6)
C     INFOLI(11)=IWA(8)+INFOLI(8)
C     INFOLI(12)=IWA(11)+INFOLI(11)
C     INFOLI(13)=IWA(12)+INFOLI(12)
C
C     DO 50 I=1,IWA(9)
C       HILF=GETIWA(IFN,I,1,IWA,LIWA,INFOLI)
C       IF (HILF.EQ.1) THEN
C         NFUNC=NFUNC + GETIWA(IFN,I,3,IWA,LIWA,INFOLI)-1
C       ENDIF
C  50 CONTINUE
C
C     DO 30 I=1,PWA
C       READ(SYMFIL,'(D24.17)',ERR=100) WA(I)
C30   CONTINUE
C     IERR=0
C     UWA=PWA+PX
C     UIWA=PIWA+PIX
C     RETURN
C100  IERR=26
C     RETURN
C     END
C
C
C
      SUBROUTINE SYMFUN (X,N,F,M,ACTIVE,WA,LWA,IWA,LIWA,DFX,DFXLEN,IERR)
      INTEGER N,M
      DOUBLE PRECISION X(N),F(M)
      LOGICAL ACTIVE(M)
      INTEGER LWA,LIWA
      DOUBLE PRECISION WA(LWA)
      INTEGER IWA(LIWA)
      INTEGER DFXLEN
      INTEGER DFX(*)
      INTEGER IERR
C
C**********************************************************************
C
C   S Y M F U N   -   EVALUATE SYMBOLICALLY DEFINED FUNCTIONS.
C
C   PARAMETERS:
C      X(N)        - ON INPUT, THE ONE-DIMENSIONAL ARRAY X HAS TO
C                    CONTAIN THE ARGUMENT THE FUNCTIONS ARE TO BE
C                    COMPUTED AT.
C      F(M)        - ON RETURN, F CONTAINS THE VALUES OF THE ACTIVE
C                    FUNCTIONS AT ARGUMENT X.
C      ACTIVE(M)   - THE LOGICAL ARRAY SPECIFIES WHICH OF THE M 
C                    FUNCTIONS ARE TO BE COMPUTED (ACTIVE(K) = .TRUE.).
C      WA(LWA)     - REAL WORKING ARRAY, CONTAINS THE INTERMEDIATE CODE 
C                    GENERATED BY SYMINP.
C      IWA(LIWA)   - INTEGER WORKING ARRAY, CF. WA(LWA). 
C      DFX(DFXLEN) - THE ARRAY SPECIFIES WHICH FIRST AND SECOND
C                    DERIVATIVES ARE TO BE COMPUTED BY CONTAINING THE
C                    NUMBER OF THE VARIABLES
C      IERR        - THE PARAMETER SHOWS THE REASON FOR TERMINATING THE 
C                    SUBROUTINE. ON RETURN IERR COULD CONTAIN THE
C                    FOLLOWING VALUES:
C                    IERR = 0 : SUCCESSFUL TERMINATION.
C                    IERR > 0 : AN ERROR HAS BEEN DETECTED. FOR FURTHER
C                             INFORMATION CF. SUBROUTINE SYMERR.
C
C**********************************************************************
C
      INTEGER PIIS,PVIS,PIIC,PVIC,PIRC,PVRC,PIVA,PVVA,PIFN,PXFN,PVFN
C     INTEGER PVGR,PVHE,PVPF,PIV,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA
      INTEGER PVGR,PVHE,PVPF,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA
      INTEGER MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN
      INTEGER MPVFN,MPVPF,MPIV,MODE
      INTEGER LVVA,LIFN,LVFN,LVGR,LVHE,LVPF,LIV,LGST,LHST
C
      INTEGER IVAL
      DOUBLE PRECISION FVAL
      INTEGER I,J,K,L,PC
C
      INTEGER GSMDEP
      PARAMETER (GSMDEP=10)
C
      INTEGER NOGRAD
      PARAMETER (NOGRAD=0)
C
      IERR=0
C
      PIIS=IWA(1)
      PVIS=IWA(2)
      PIIC=IWA(3)
      PVIC=IWA(4)
      PIRC=IWA(5)
      PVRC=IWA(6)
      PIVA=IWA(7)
      PVVA=IWA(8)
      PIFN=IWA(9)
      PXFN=IWA(10)
      PVFN=IWA(11)
      PVGR=IWA(12)
      PVHE=1
      PVPF=IWA(14)
C     PIV=IWA(15)
C
      LIIS=16
      LVIS=LIIS+PIIS*5
      LIIC=LVIS+PVIS
      LVIC=LIIC+PIIC*4
      LIRC=LVIC+PVIC
      LVRC=1
      LIVA=LIRC+PIRC*4
      LVVA=LVRC+PVRC
      LIFN=LIVA+PIVA*3
      LVFN=LVVA+PVVA
      LVGR=LVFN+PVFN
      LVHE=LVGR+PVGR*PVVA
      LVPF=LIFN+PIFN*7
      LIV=LVPF+PVPF
      LGST=LVHE+PVHE*1
      LHST=LGST+GSMDEP*PVVA
C
      IF (N .NE. PVVA) THEN
        IERR=43
        RETURN
      ENDIF
      IF (M .NE. PVFN-(PIFN-PXFN)) THEN
        IERR=44
        RETURN
      ENDIF
      DO 10 I=1,N
        WA(LVVA+I-1)=X(I)
 10   CONTINUE
      K=0
      DO 20 I=1,PXFN
        IF (IWA(LIFN+(I-1)+(1-1)*PIFN) .EQ. 0) THEN
          K=K+1
          IF (ACTIVE(K)) THEN
            PC=IWA(LIFN+(I-1)+(6-1)*PIFN)
            MPIIS=MAX(1,IWA(1))
            MPVIS=MAX(1,IWA(2))
            MPIIC=MAX(1,IWA(3))
            MPVIC=MAX(1,IWA(4))
            MPIRC=MAX(1,IWA(5))
            MPVRC=MAX(1,IWA(6))
            MPIVA=MAX(1,IWA(7))
            MPVVA=MAX(1,IWA(8))
            MPIFN=MAX(1,IWA(9))
            MPVFN=MAX(1,IWA(11))
            MPVPF=MAX(1,IWA(14))
            MPIV =MAX(1,IWA(15))
            MODE=NOGRAD
            DO 30 L=LIV,LIV+MPIV-1
              IWA(L)=0
 30         CONTINUE
            CALL EVAL(PC,IVAL,FVAL,MODE,PVVA,1,
     1                MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,
     2                MPIFN,MPVFN,PVHE,MPVPF,MPIV,IWA(LIIS),IWA(LVIS),
     3                IWA(LIIC),IWA(LVIC),IWA(LIRC),WA(LVRC),IWA(LIVA),
     4                WA(LVVA),IWA(LIFN),WA(LVFN),WA(LVGR),WA(LVHE),
     5                IWA(LVPF),IWA(LIV),WA(LGST),WA(LHST),DFX,DFXLEN,
     6                IERR)
            IF (IERR .NE. 0) RETURN
            F(K)=WA(LVFN+IWA(LIFN+(I-1)+(4-1)*PIFN)-1)
          ENDIF
        ELSE IF (IWA(LIFN+(I-1)+(1-1)*PIFN) .EQ. 1) THEN
          DO 15 J=1,IWA(LIFN+(I-1)+(3-1)*PIFN)
            K=K+1
            IF (ACTIVE(K)) THEN
              PC=IWA(LIFN+(I-1)+(6-1)*PIFN)
              MPIIS=MAX(1,IWA(1))
              MPVIS=MAX(1,IWA(2))
              MPIIC=MAX(1,IWA(3))
              MPVIC=MAX(1,IWA(4))
              MPIRC=MAX(1,IWA(5))
              MPVRC=MAX(1,IWA(6))
              MPIVA=MAX(1,IWA(7))
              MPVVA=MAX(1,IWA(8))
              MPIFN=MAX(1,IWA(9))
              MPVFN=MAX(1,IWA(11))
              MPVPF=MAX(1,IWA(14))
              MPIV =MAX(1,IWA(15))
              MODE=NOGRAD
              DO 40 L=LIV,LIV+MPIV-1
                IWA(L)=0
 40           CONTINUE
              IWA(LIV+IWA(LIFN+(I-1)+(2-1)*PIFN)-1)=J
              CALL EVAL(PC,IVAL,FVAL,MODE,PVVA,1,
     1                MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,
     2                MPIFN,MPVFN,PVHE,MPVPF,MPIV,IWA(LIIS),IWA(LVIS),
     3                IWA(LIIC),IWA(LVIC),IWA(LIRC),WA(LVRC),IWA(LIVA),
     4                WA(LVVA),IWA(LIFN),WA(LVFN),WA(LVGR),WA(LVHE),
     5                IWA(LVPF),IWA(LIV),WA(LGST),WA(LHST),DFX,DFXLEN,
     6                IERR)
              IF (IERR .NE. 0) RETURN
              F(K)=WA(LVFN+IWA(LIFN+(I-1)+(4-1)*PIFN)-1+(J-1))
            ENDIF
 15       CONTINUE
        ENDIF
 20   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE SYMGRA (X,N,F,M,DF,MMAX,ACTIVE,WA,LWA,IWA,LIWA,DFX,
     1                   DFXLEN,IERR)
      INTEGER N,M,MMAX
      DOUBLE PRECISION X(N),F(M),DF(MMAX,N)
      LOGICAL ACTIVE(M)
      INTEGER LWA,LIWA
      DOUBLE PRECISION WA(LWA)
      INTEGER IWA(LIWA)
      INTEGER DFXLEN
      INTEGER DFX(*)
      INTEGER IERR
C
C***********************************************************************
C
C   S Y M G R A   -   EVALUATE SYMBOLICALLY DEFINED FUNCTIONS AND
C                     CORRESPONDING GRADIENTS.
C
C   PARAMETERS:
C      X(N)       - ON INPUT, THE ONE-DIMENSIONAL ARRAY X HAS TO CONTAIN
C                   THE ARGUMENT THE FUNCTIONS ARE TO BE COMPUTED AT.
C      F(M)       - ON RETURN, F CONTAINS THE VALUES OF THE ACTIVE
C                   FUNCTIONS AT ARGUMENT X.
C      DF(MMAX,N) - ON RETURN, DF CONTAINS THE GRADIENTS OF THE ACTIVE
C                   FUNCTIONS AT ARGUMENT X. IN THE DRIVING PROGRAM, THE
C                   ROW DIMENSION OF DF HAS TO BE EQUAL TO MMAX.
C      ACTIVE(M)  - THE LOGICAL ARRAY SPECIFIES WHICH OF THE M FUNCTIONS
C                   ARE TO BE COMPUTED ( ACTIVE(K) = .TRUE. ).
C      WA(LWA)    - REAL WORKING ARRAY, CONTAINS THE INTERMEDIATE CODE
C                   GENERATED BY SYMINP.
C      IWA(LIWA)  - INTEGER WORKING ARRAY, CF. WA(LWA).
C      DFX(DFXLEN)- THE ARRAY SPECIFIES WHICH FIRST AND SECOND
C                   DERIVATIVES ARE TO BE COMPUTED BY CONTAINING THE
C                   NUMBER OF THE VARIABLES
C      IERR       - THE PARAMETER SHOWS THE REASON FOR TERMINATING THE
C                   SUBROUTINE. ON RETURN IERR COULD CONTAIN THE FOLLOW-
C                   ING VALUES:
C                   IERR = 0 : SUCCESSFUL TERMINATION.
C                   IERR > 0 : AN ERROR HAS BEEN DETECTED. FOR FURTHER
C                              INFORMATION CF. SUBROUTINE SYMERR.
C
C***********************************************************************
C
      INTEGER PIIS,PVIS,PIIC,PVIC,PIRC,PVRC,PIVA,PVVA,PIFN,PXFN,PVFN
C     INTEGER PVGR,PVHE,PVPF,PIV,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA
      INTEGER PVGR,PVHE,PVPF,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA
      INTEGER LVVA,LIFN,LVFN,LVGR,LVHE,LVPF,LIV,LGST,LHST
      INTEGER MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN
      INTEGER MPVFN,MPVPF,MPIV,MODE
C
      INTEGER IVAL
      DOUBLE PRECISION FVAL
      INTEGER I,J,K,L,PC
C
      INTEGER GSMDEP
      PARAMETER (GSMDEP=10)
C
      INTEGER GRAD
      PARAMETER (GRAD=1)
C
      IERR=0
C
      PIIS=IWA(1)
      PVIS=IWA(2)
      PIIC=IWA(3)
      PVIC=IWA(4)
      PIRC=IWA(5)
      PVRC=IWA(6)
      PIVA=IWA(7)
      PVVA=IWA(8)
      PIFN=IWA(9)
      PXFN=IWA(10)
      PVFN=IWA(11)
      PVGR=IWA(12)
      IF (PVGR .EQ. 0) THEN
        IERR=67
        RETURN
      ENDIF
      PVHE=1
      PVPF=IWA(14)
C     PIV=IWA(15)
      LIIS=16
      LVIS=LIIS+PIIS*5
      LIIC=LVIS+PVIS
      LVIC=LIIC+PIIC*4
      LIRC=LVIC+PVIC
      LVRC=1
      LIVA=LIRC+PIRC*4
      LVVA=LVRC+PVRC
      LIFN=LIVA+PIVA*3
      LVFN=LVVA+PVVA
      LVGR=LVFN+PVFN
      LVHE=LVGR+PVGR*PVVA
      LVPF=LIFN+PIFN*7
      LIV=LVPF+PVPF
      LGST=LVHE+PVHE*1
      LHST=LGST+GSMDEP*PVVA
      IF (N .NE. PVVA) THEN
        IERR=43
        RETURN
      ENDIF
      IF (M .NE. PVFN-(PIFN-PXFN)) THEN
        IERR=44
        RETURN
      ENDIF
      DO 10 I=1,N
        WA(LVVA+I-1)=X(I)
 10   CONTINUE
      K=0
      DO 20 I=1,PXFN
        IF (IWA(LIFN+(I-1)+(1-1)*PIFN) .EQ. 0) THEN
          K=K+1
          IF (ACTIVE(K)) THEN
            PC=IWA(LIFN+(I-1)+(6-1)*PIFN)
            MPIIS=MAX(1,IWA(1))
            MPVIS=MAX(1,IWA(2))
            MPIIC=MAX(1,IWA(3))
            MPVIC=MAX(1,IWA(4))
            MPIRC=MAX(1,IWA(5))
            MPVRC=MAX(1,IWA(6))
            MPIVA=MAX(1,IWA(7))
            MPVVA=MAX(1,IWA(8))
            MPIFN=MAX(1,IWA(9))
            MPVFN=MAX(1,IWA(11))
            MPVPF=MAX(1,IWA(14))
            MPIV =MAX(1,IWA(15))
            MODE=GRAD
            DO 30 L=LIV,LIV+MPIV-1
              IWA(L)=0
 30         CONTINUE
            CALL EVAL(PC,IVAL,FVAL,MODE,PVVA,1,
     1                MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,
     2                MPIFN,MPVFN,PVHE,MPVPF,MPIV,IWA(LIIS),IWA(LVIS),
     3                IWA(LIIC),IWA(LVIC),IWA(LIRC),WA(LVRC),IWA(LIVA),
     4                WA(LVVA),IWA(LIFN),WA(LVFN),WA(LVGR),WA(LVHE),
     5                IWA(LVPF),IWA(LIV),WA(LGST),WA(LHST),DFX,DFXLEN,
     6                IERR)
            IF (IERR .NE. 0) RETURN
            F(K)=WA(LVFN+IWA(LIFN+(I-1)+(4-1)*PIFN)-1)
            DO 11 L=1,PVVA
              DF(K,L)=WA(LVGR+(IWA(LIFN+(I-1)+(5-1)*PIFN)-1)+(L-1)*PVGR)
 11         CONTINUE
          ENDIF
        ELSE IF (IWA(LIFN+(I-1)+(1-1)*PIFN) .EQ. 1) THEN
          DO 15 J=1,IWA(LIFN+(I-1)+(3-1)*PIFN)
            K=K+1
            IF (ACTIVE(K)) THEN
              PC=IWA(LIFN+(I-1)+(6-1)*PIFN)
              MPIIS=MAX(1,IWA(1))
              MPVIS=MAX(1,IWA(2))
              MPIIC=MAX(1,IWA(3))
              MPVIC=MAX(1,IWA(4))
              MPIRC=MAX(1,IWA(5))
              MPVRC=MAX(1,IWA(6))
              MPIVA=MAX(1,IWA(7))
              MPVVA=MAX(1,IWA(8))
              MPIFN=MAX(1,IWA(9))
              MPVFN=MAX(1,IWA(11))
              MPVPF=MAX(1,IWA(14))
              MPIV =MAX(1,IWA(15))
              MODE=GRAD
              DO 40 L=LIV,LIV+MPIV-1
                IWA(L)=0
 40           CONTINUE
              IWA(LIV+IWA(LIFN+(I-1)+(2-1)*PIFN)-1)=J
              CALL EVAL(PC,IVAL,FVAL,MODE,PVVA,1,
     1                MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,
     2                MPIFN,MPVFN,PVHE,MPVPF,MPIV,IWA(LIIS),IWA(LVIS),
     3                IWA(LIIC),IWA(LVIC),IWA(LIRC),WA(LVRC),IWA(LIVA),
     4                WA(LVVA),IWA(LIFN),WA(LVFN),WA(LVGR),WA(LVHE),
     5                IWA(LVPF),IWA(LIV),WA(LGST),WA(LHST),DFX,DFXLEN,
     6                IERR)
              IF (IERR .NE. 0) RETURN
              F(K)=WA(LVFN+IWA(LIFN+(I-1)+(4-1)*PIFN)+(J-1)-1)
              DO 12 L=1,PVVA
                DF(K,L)=WA(LVGR+(IWA(LIFN+(I-1)+(5-1)*PIFN)-1+(J-1))+
     1              (L-1)*PVGR)
 12           CONTINUE
            ENDIF
 15       CONTINUE
        ENDIF
 20   CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE SYMHES (X,N,F,M,DF,DDF,MMAX,ACTIVE,WA,LWA,IWA,LIWA,
     1                   DFX,DFXLEN,IERR)
      INTEGER N,M,MMAX
      DOUBLE PRECISION X(N),F(M),DF(MMAX,N),DDF(MMAX,N*N)
      LOGICAL ACTIVE(M)
      INTEGER LWA,LIWA
      DOUBLE PRECISION WA(LWA)
      INTEGER IWA(LIWA)
      INTEGER DFXLEN
      INTEGER DFX(*)
      INTEGER IERR
C
C***********************************************************************
C
C   S Y M H E S   - EVALUATE SYMBOLICALLY DEFINED FUNCTIONS,
C                   CORRESPONDING GRADIENTS AND THE HESSIAN MATRIX
C
C   PARAMETERS:
C      X(N)          - ON INPUT, THE ONE-DIMENSIONAL ARRAY X HAS TO
C                      CONTAIN THE ARGUMENT THE FUNCTIONS ARE TO BE
C                      COMPUTED AT.
C      F(M)          - ON RETURN, F CONTAINS THE VALUES OF THE ACTIVE
C                      FUNCTIONS AT ARGUMENT X.
C      DF(MMAX,N)    - ON RETURN, DF CONTAINS THE GRADIENTS OF THE
C                      ACTIVE FUNCTIONS AT ARGUMENT X. IN THE DRIVING
C                      PROGRAM, THE ROW DIMENSION OF DF HAS TO BE EQUAL
C                      TO MMAX.
C      DDF(MMAX,N*N) - ON RETURN, DDF CONTAINS THE SECOND DERIVATIVES OF
C                      THE ACTIVE FUNCTIONS AT ARGUMENT X. IN THE
C                      DRIVING PROGRAM, THE ROW DIMENSION OF DDF HAS TO
C                      BE EQUAL TO MMAX
C      ACTIVE(M)     - THE LOGICAL ARRAY SPECIFIES WHICH OF THE M
C                      FUNCTIONS ARE TO BE COMPUTED (ACTIVE(K)=.TRUE.).
C      WA(LWA)       - REAL WORKING ARRAY, CONTAINS THE INTERMEDIATE
C                      CODE GENERATED BY SYMINP.
C      IWA(LIWA)     - INTEGER WORKING ARRAY, CF. WA(LWA).
C      DFX(DFXLEN)   - THE ARRAY SPECIFIES WHICH FIRST AND SECOND
C                      DERIVATIVES ARE TO BE COMPUTED BY CONTAINING THE
C                      NUMBER OF THE VARIABLES
C      IERR          - THE PARAMETER SHOWS THE REASON FOR TERMINATING
C                      THE SUBROUTINE. ON RETURN IERR COULD CONTAIN THE
C                      FOLLOWING VALUES:
C                      IERR = 0 : SUCCESSFUL TERMINATION.
C                      IERR > 0 : AN ERROR HAS BEEN DETECTED. FOR
C                                 FURTHER INFORMATION CF. SUBROUTINE
C                                 SYMERR.
C
C***********************************************************************
C
      INTEGER PIIS,PVIS,PIIC,PVIC,PIRC,PVRC,PIVA,PVVA,PIFN,PXFN,PVFN
C     INTEGER PVGR,PVHE,PVPF,PIV,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA
      INTEGER PVGR,PVHE,PVPF,LIIS,LVIS,LIIC,LVIC,LIRC,LVRC,LIVA
      INTEGER LVVA,LIFN,LVFN,LVGR,LVHE,LVPF,LIV,LGST,LHST
      INTEGER MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,MPIFN
      INTEGER MPVFN,MPVPF,MPIV,MODE
C
      INTEGER IVAL
      DOUBLE PRECISION FVAL
      INTEGER I,J,K,L,P,PC
C
      INTEGER GSMDEP
      PARAMETER (GSMDEP=10)
C
      INTEGER HESS
      PARAMETER (HESS=2)
C
      IERR=0
C
      PIIS=IWA(1)
      PVIS=IWA(2)
      PIIC=IWA(3)
      PVIC=IWA(4)
      PIRC=IWA(5)
      PVRC=IWA(6)
      PIVA=IWA(7)
      PVVA=IWA(8)
      PIFN=IWA(9)
      PXFN=IWA(10)
      PVFN=IWA(11)
      PVGR=IWA(12)
      PVHE=IWA(13)
      IF (PVHE .EQ. 0) THEN
        IERR=62
        RETURN
      ENDIF
      PVPF=IWA(14)
C     PIV=IWA(15)
      LIIS=16
      LVIS=LIIS+PIIS*5
      LIIC=LVIS+PVIS
      LVIC=LIIC+PIIC*4
      LIRC=LVIC+PVIC
      LVRC=1
      LIVA=LIRC+PIRC*4
      LVVA=LVRC+PVRC
      LIFN=LIVA+PIVA*3
      LVFN=LVVA+PVVA
      LVGR=LVFN+PVFN
      LVHE=LVGR+PVGR*PVVA
      LVPF=LIFN+PIFN*7
      LIV=LVPF+PVPF
      LGST=LVHE+PVHE*PVVA*PVVA
      LHST=LGST+GSMDEP*PVVA
      IF (N .NE. PVVA) THEN
        IERR=43
        RETURN
      ENDIF
      IF (M .NE. PVFN-(PIFN-PXFN)) THEN
        IERR=44
        RETURN
      ENDIF
      DO 10 I=1,N
        WA(LVVA+I-1)=X(I)
 10   CONTINUE
      K=0
      DO 30 I=1,PXFN
        IF (IWA(LIFN+(I-1)+(1-1)*PIFN) .EQ. 0) THEN
          K=K+1
          IF (ACTIVE(K)) THEN
            PC=IWA(LIFN+(I-1)+(6-1)*PIFN)
            MPIIS=MAX(1,IWA(1))
            MPVIS=MAX(1,IWA(2))
            MPIIC=MAX(1,IWA(3))
            MPVIC=MAX(1,IWA(4))
            MPIRC=MAX(1,IWA(5))
            MPVRC=MAX(1,IWA(6))
            MPIVA=MAX(1,IWA(7))
            MPVVA=MAX(1,IWA(8))
            MPIFN=MAX(1,IWA(9))
            MPVFN=MAX(1,IWA(11))
            MPVPF=MAX(1,IWA(14))
            MPIV =MAX(1,IWA(15))
            MODE=HESS
            DO 40 L=LIV,LIV+MPIV-1
              IWA(L)=0
 40         CONTINUE
            CALL EVAL(PC,IVAL,FVAL,MODE,PVVA,PVVA,
     1                MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,
     2                MPIFN,MPVFN,PVHE,MPVPF,MPIV,IWA(LIIS),IWA(LVIS),
     3                IWA(LIIC),IWA(LVIC),IWA(LIRC),WA(LVRC),IWA(LIVA),
     4                WA(LVVA),IWA(LIFN),WA(LVFN),WA(LVGR),WA(LVHE),
     5                IWA(LVPF),IWA(LIV),WA(LGST),WA(LHST),DFX,DFXLEN,
     6                IERR)
            IF (IERR .NE. 0) RETURN
            F(K)=WA(LVFN+IWA(LIFN+(I-1)+(4-1)*PIFN)-1)
            DO 11 L=1,PVVA
              DF(K,L)=WA(LVGR+(IWA(LIFN+(I-1)+(5-1)*PIFN)-1)+(L-1)*PVGR)
 11         CONTINUE
            DO 13 L=1,PVVA
              DO 12 P=1,PVVA
                DDF(K,(L-1)*PVVA+P)=WA(LVHE+(IWA(LIFN+(I-1)+(7-1)*PIFN)
     1              -1)+((L-1)*PVVA+P-1)*PVHE)
 12           CONTINUE
 13         CONTINUE
          ENDIF
        ELSE IF (IWA(LIFN+(I-1)+(1-1)*PIFN) .EQ. 1) THEN
          DO 20 J=1,IWA(LIFN+(I-1)+(3-1)*PIFN)
            K=K+1
            IF (ACTIVE(K)) THEN
              PC=IWA(LIFN+(I-1)+(6-1)*PIFN)
              MPIIS=MAX(1,IWA(1))
              MPVIS=MAX(1,IWA(2))
              MPIIC=MAX(1,IWA(3))
              MPVIC=MAX(1,IWA(4))
              MPIRC=MAX(1,IWA(5))
              MPVRC=MAX(1,IWA(6))
              MPIVA=MAX(1,IWA(7))
              MPVVA=MAX(1,IWA(8))
              MPIFN=MAX(1,IWA(9))
              MPVFN=MAX(1,IWA(11))
              MPVPF=MAX(1,IWA(14))
              MPIV =MAX(1,IWA(15))
              MODE=HESS
              DO 50 L=LIV,LIV+MPIV-1
                IWA(L)=0
 50           CONTINUE
              IWA(LIV+IWA(LIFN+(I-1)+(2-1)*PIFN)-1)=J
              CALL EVAL(PC,IVAL,FVAL,MODE,PVVA,PVVA,
     1                MPIIS,MPVIS,MPIIC,MPVIC,MPIRC,MPVRC,MPIVA,MPVVA,
     2                MPIFN,MPVFN,PVHE,MPVPF,MPIV,IWA(LIIS),IWA(LVIS),
     3                IWA(LIIC),IWA(LVIC),IWA(LIRC),WA(LVRC),IWA(LIVA),
     4                WA(LVVA),IWA(LIFN),WA(LVFN),WA(LVGR),WA(LVHE),
     5                IWA(LVPF),IWA(LIV),WA(LGST),WA(LHST),DFX,DFXLEN,
     6                IERR)
              IF (IERR .NE. 0) RETURN
              F(K)=WA(LVFN+IWA(LIFN+(I-1)+(4-1)*PIFN)+(J-1)-1)
              DO 14 L=1,PVVA
                DF(K,L)=WA(LVGR+(IWA(LIFN+(I-1)+(5-1)*PIFN)-1+(J-1))+
     1              (L-1)*PVGR)
 14           CONTINUE
              DO 16 L=1,PVVA
                DO 15 P=1,PVVA
                  DDF(K,(L-1)*PVVA+P)=WA(LVHE+(IWA(LIFN+(I-1)+
     1                  (7-1)*PIFN)-1+(J-1))+((L-1)*PVVA+P-1)*PVHE)
 15             CONTINUE
 16           CONTINUE
            ENDIF
 20       CONTINUE
        ENDIF
 30   CONTINUE
      RETURN
      END

