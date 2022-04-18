*
*     Define the standard error codes returned by the quadrature
*     package.
*
*     Success:
*
      INTEGER             eokay
      PARAMETER           (eokay = 0)
*
*     Eigensolution could not be obtained:
*
      INTEGER             eeigen
      PARAMETER           (eeigen = 1)
*
*     Destructive overflow:
*
      INTEGER             eovflo
      PARAMETER           (eovflo = 2)
*
*     nquad out of range:
*
      INTEGER             enquad
      PARAMETER           (enquad = 3)
*
*     alpha out of range:
*
      INTEGER             ealpha
      PARAMETER           (ealpha = 4)
*
*     beta out of range:
*
      INTEGER             ebeta
      PARAMETER           (ebeta = 5)
