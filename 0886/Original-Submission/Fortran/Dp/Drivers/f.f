      DOUBLE PRECISION FUNCTION F(FTYPE,X1,X2)
C
C     ******************************************************************
C
C                                                           From Padua2D
C
C                                    Marco Caliari and Stefano De Marchi
C                                         Department of Computer Science
C                                           University of Verona (Italy)
C                              {marco.caliari,stefano.demarchi}@univr.it
C
C                                                         Marco Vianello
C                             Department of Pure and Applied Mathematics
C                                            University of Padua (Italy)
C                                                   marcov@math.unipd.it
C                                                    November 14th, 2007
C
C       This function is an interface of Padua2D to the function TSTFN2
C     from CSHEP2D by Robert J. Renka.
C
C     ******************************************************************
C
C     PARAMS type                i=input, o=output, a=auxiliary           
C     
C     FTYPE  integer             i: function number
C     X1     double              i: first coordinate of the point
C     X2     double              i: second coordinate of the point
C
C     ******************************************************************
C
      INTEGER FTYPE
      DOUBLE PRECISION X1,X2
      IF ((FTYPE .LT. 1) .OR. (FTYPE .GT. 10)) THEN
         WRITE(6,*) 'IN _F: FTYPE ',FTYPE,' UNDEFINED'
         WRITE(6,*) 'STOP'
         STOP
      END IF
      CALL TSTFN2(FTYPE,X1,X2,0,F,F,F,F,F,F)
      RETURN    
      END
