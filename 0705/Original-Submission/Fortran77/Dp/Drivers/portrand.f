      REAL FUNCTION RAND()
C..PORTABLE RANDOM NUMBER GENERATOR
C..USING THE RECURSION -
C..
C..       IX=IX*A MOD P
      INTEGER A,P,IX,B15,B16,XHI,XALO,LEFTLO,FHI,K
C..7**5, 2**15,  2**16, 2**31-1
      PARAMETER (A=16807, B15=32768, B16=65536,
     +           P=2147483647)
      COMMON /RANCOM/ IX
      SAVE /RANCOM/
C
C..GET 15 HIGH ORDER BITS OF IX
      XHI=IX/B16
C..GET 16 LO BITS OF IX AND FORM LO PRODUCT
      XALO=(IX-XHI*B16)*A
C..GET 15 HI ORDER BITS OF LO PRODUCT
      LEFTLO=XALO/B16
C..FORM THE 31 HIGHEST BITS OF FULL PRODUCT
      FHI=XHI*A + LEFTLO
C..GET OVERFLOW PAST 31ST BIT OF FULL PRODUCT
      K=FHI/B15
C..ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C..THE PARENTHESES ARE ESSENTIAL
      IX=(((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16)+K
C..ADD P BACK IN IF NECESSARY
      IF(IX.LT.0)IX=IX+P
C..MULTIPLY BY 1/(2**31-1)
      RAND=FLOAT(IX)*4.656612875E-10
      RETURN
      END

      subroutine seed(iseed)
c
c seed the portable random number generator
c
      integer iseed, ix
      common /rancom/ ix
      SAVE /RANCOM/
      ix = iseed
      end