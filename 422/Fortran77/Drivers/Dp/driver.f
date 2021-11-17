      PROGRAM SPANTR
c
c  A driver program to demonstrate Algorithm 422.
c  This example shows how to find the minimum spanning
c  tree which connects 6 large cities in Michigan.
c  A graph which stores the distances in miles between
c  each city is read in as file spantr.in and is then
c  used to load array bgraph(6,6).  The output file
c  lists the edges which make up the minimal spanning
c  tree, these edges are stored in array mstree(2,6).
c  This driver program was compiled and tested using
c  Fortran 77.
c
c  Note: The input file spantr.in should contain only
c  eight lines; two header lines and then six lines
c  containing the six city names and the 6-by-6 matrix
c  of distances between these cities (in miles).
c
c     buf1           is an 80-char scratch buffer
c     sname(6)       holds the 14-char city names
c     bgraph(50,50)  holds the graph and the edges between the 6 nodes
c     mstree(2,50)   holds the nodes for the minimal spanning tree
c     numbls         is the number of edges in the minimal spanning tree
c     sumlen         is the sum of the edge lengths in the minimal
c                    spanning tree
c
c
C     .. Local Scalars ..
      DOUBLE PRECISION SUMLEN
      INTEGER I,J,NIN,NOUT,NUMBLS
      CHARACTER*80 BUF1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION BGRAPH(NMAX,NMAX)
      INTEGER MSTREE(2,NMAX)
      CHARACTER*14 SNAME(NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL DMTOMS
C     ..
C     OPEN (10,FILE='data',STATUS='old',ERR=30)
C     OPEN (20,FILE='res',STATUS='unknown',ERR=40)
c
C Use port library routine to set default input and
C output channels
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=50)
C     ..
C     .. External Functions ..
      INTEGER I1MACH
      EXTERNAL I1MACH
C     ..
      NIN = I1MACH(1)
      NOUT = I1MACH(2)
      READ (NIN,FMT='(a)') BUF1
      WRITE (NOUT,FMT='(//,1x,a)') BUF1
      READ (NIN,FMT='(a)') BUF1
      WRITE (NOUT,FMT='(1x,a)') BUF1
      DO 10 I = 1,6
          READ (NIN,FMT=9000) SNAME(I), (BGRAPH(I,J),J=1,6)
          WRITE (NOUT,FMT=9010) SNAME(I), (BGRAPH(I,J),J=1,6)
   10 CONTINUE
C     CLOSE (10)
c
      CALL DMTOMS(BGRAPH,6,MSTREE,NUMBLS,SUMLEN)
c
      WRITE (NOUT,FMT=
     +  '(//,'' number of edges in minimal spanning tree = '',    i8)')
     +  NUMBLS
      WRITE (NOUT,FMT='(/,'' sum of the edges (miles) = '',f8.4,//)')
     +  SUMLEN
      DO 20 I = 1,NUMBLS
          WRITE (NOUT,FMT='('' from,to: '',2i4,2x,2a18)') MSTREE(1,I),
     +      MSTREE(2,I),SNAME(MSTREE(1,I)),SNAME(MSTREE(2,I))
   20 CONTINUE
C     CLOSE (20)
      STOP
c
C  30 WRITE(NOUT,FMT='(//,'' Error opening spantr.in !! '',/)')
C     STOP

C  40 WRITE(NOUT,FMT='(//,'' Error opening spantr.out !! '',/)')
C     STOP

 9000 FORMAT (A14,6F7.1)
 9010 FORMAT (1X,A14,6F7.1)
      END
