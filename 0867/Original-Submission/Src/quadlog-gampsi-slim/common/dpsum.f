      SUBROUTINE dpsum (sumneg, sumpos, term)
************************************************************************
*     (Double-precision partial sums)
*     Accumulate term in sumneg (if term < 0) or sumpos (if term >= 0).
*
*     If term is a NaN, it is accumulated in sumneg.
*     (27-Mar-2000)
************************************************************************
*
*     Parameter variables
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
*     Argument variables
*
      DOUBLE PRECISION    sumneg,      sumpos,      term
*
      IF (term .ge. zero) THEN
          sumpos = sumpos + term
      ELSE
          sumneg = sumneg + term
      END IF

      END
