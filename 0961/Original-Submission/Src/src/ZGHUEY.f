      SUBROUTINE ZGHUEY( S11, S12, H11, H12, CO, SI )
C
C     PURPOSE
C
C     To compute a unitary matrix Q for a complex regular 2-by-2 
C     skew-Hamiltonian/Hamiltonian pencil aS - bH with
C
C         (  S11  S12  )        (  H11  H12  )
C     S = (            ),   H = (            ),
C         (   0   S11' )        (   0  -H11' )
C
C     such that J Q' J' (aS - bH) Q is upper triangular but the
C     eigenvalues are in reversed order. The matrix Q is represented by
C
C         (  CO  SI  )
C     Q = (          ).
C         ( -SI' CO  )
C
C     The notation M' denotes the conjugate transpose of the matrix M.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     S11     (input) COMPLEX*16
C             Upper left element of the skew-Hamiltonian matrix S.
C
C     S12     (input) COMPLEX*16
C             Upper right element of the skew-Hamiltonian matrix S.
C
C     H11     (input) COMPLEX*16
C             Upper left element of the Hamiltonian matrix H.
C
C     H12     (input) COMPLEX*16
C             Upper right element of the Hamiltonian matrix H.
C
C     CO      (output) DOUBLE PRECISION
C             Upper left element of Q.
C
C     SI      (output) COMPLEX*16
C             Upper right element of Q.
C
C     METHOD
C
C     The algorithm uses unitary transformations as described on page 43
C     in [1].
C
C     REFERENCES
C
C     [1] Benner, P., Byers, R., Mehrmann, V. and Xu, H.
C         Numerical Computation of Deflating Subspaces of Embedded
C         Hamiltonian Pencils.
C         Tech. Rep. SFB393/99-15, Technical University Chemnitz,
C         Germany, June 1999.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically backward stable.
C
C     CONTRIBUTOR
C
C     M. Voigt, Technische Universitaet Chemnitz, Jan. 2009.
C
C     REVISIONS
C
C     V. Sima, Aug. 2009, Jan. 2010, June 2014.
C     M. Voigt, Jan. 2012, Jul. 2013.
C
C     KEYWORDS
C
C     Eigenvalue exchange, skew-Hamiltonian/Hamiltonian pencil, upper
C     triangular matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
      INTEGER            DKIND
      PARAMETER          ( DKIND = KIND( TWO ) )
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   CO
      COMPLEX(DKIND) ::  H11, H12, S11, S12, SI
C
C     .. Local Scalars ..
      COMPLEX(DKIND) ::  G, TMP
C
C     .. External Subroutines ..
      EXTERNAL           ZLARTG
C
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, DBLE
C
C     .. Executable Statements ..
C
      G = TWO*DBLE( H11*CONJG( S11 ) )
      CALL ZLARTG( CONJG( S11 )*H12 + S12*CONJG( H11 ), G, CO, SI, TMP )
C
      RETURN
C *** Last line of ZGHUEY ***
      END
