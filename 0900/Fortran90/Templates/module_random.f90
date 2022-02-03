!*********************************************************************************************************************************
!*
!* MODULE: RANDOM
!*
!* PURPOSE: it generates random numbers and random vectors with prescribed distribution
!*
!* DEPENDENCIES:
!*               - PRECISION
!* 
!* GLOBALS:
!*          - LOW (integer): number of bytes for simple precision variables (from PRECISION)
!*          - HIGH (integer): number of bytes for double precision variables (from PRECISION)
!*          - PLO (integer): precision for logical variables (from PRECISION)
!*          - PCH (integer): precision for character variables (from PRECISION)
!*          - PIN (integer): precision for integer variables (from PRECISION)
!*          - PRE (integer): precision for real variables (from PRECISION)
!*
!*          - IDUM (integer): seed for random generators, it must be set to a negative value
!* 
!* AUTHOR: Germán Ariel Torres
!*         Facultad de Matemática, Astronomía y Física
!*         Universidad Nacional de Córdoba
!*         Argentina
!*
!* EMAIL: torres@famaf.unc.edu.ar
!*
!*********************************************************************************************************************************
module random
  !*
  use precision
  !*
  implicit none
  !*
  integer( kind = pin ) , public :: idum !* seed for random generators
  !*
contains
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_PARAMETERSUP
  !*
  !* PURPOSE: it sets the seed for the random generator with uniform distribution as a negative integer related to time
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !* 
  !* OUTPUTS: none
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine random_parametersup()
    implicit none
!!$    integer( kind = pin ) , dimension( 3_pin ) :: now
    !*----------------------------------------------------------------------------------------------------------------------------
    !* The "now" function is commented in order to get the same seed and reproduce the results plotted in the paper. 
    !* In a real experiment, you should not use the same seed every time, that is, you should remove the line "idum = -1_pin" and&
    !* & uncomment the other lines
    idum = - 1_pin
!!$    call itime( now )
!!$    idum = - now( 1_pin ) * now( 2_pin ) * now( 3_pin ) - 1_pin 
    !*----------------------------------------------------------------------------------------------------------------------------
  end subroutine random_parametersup
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_UNIFORM0D
  !*
  !* PURPOSE: long period (> 2 x 10^18) random number generator of L'Ecuyer with Bays-Durham shuffle and added safeguards. Return&
  !*          &s a uniform random deviate between 0.0 and 1.0 (exclusive of end-points values). Use global variable idum wich nee&
  !*          &ds to be initialized with a negative value. 
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - RANDOM_UNIFORM0D (real): random number with uniform distribution between 0.0 and 1.0
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: this is a translation of function ran2.f of Numerical Recipes adapted to generate random numbers in simple and dou&
  !*           &ble precision and translated to free format in Fortran 90
  !*
  !*******************************************************************************************************************************
  function random_uniform0d()
    implicit none    
    real( kind = pre ) :: random_uniform0d
    integer( kind = pin ) :: im1
    integer( kind = pin ) :: im2
    integer( kind = pin ) :: imm1
    integer( kind = pin ) :: ia1
    integer( kind = pin ) :: ia2
    integer( kind = pin ) :: iq1
    integer( kind = pin ) :: iq2
    integer( kind = pin ) :: ir1
    integer( kind = pin ) :: ir2
    integer( kind = pin ) :: ntab
    integer( kind = pin ) :: ndiv
    real( kind = pre ) :: am
    real( kind = pre ) :: eps
    real( kind = pre ) :: rnmx
    parameter ( im1  = 2147483563_pin )
    parameter ( im2  = 2147483399_pin )
    parameter ( am   = 1.0_pre / im1 )
    parameter ( imm1 = im1 - 1_pin )
    parameter ( ia1  = 40014_pin )
    parameter ( ia2  = 40692_pin )
    parameter ( iq1  = 53668_pin )
    parameter ( iq2  = 52774_pin )
    parameter ( ir1  = 12211_pin )
    parameter ( ir2  = 3791_pin )
    parameter ( ntab = 32_pin )
    parameter ( ndiv = 1_pin + imm1 / ntab )
    parameter ( eps  = 1.2e-7_pre )
    parameter ( rnmx = 1.0_pre - eps ) 
    integer( kind = pin ) :: idum2
    integer( kind = pin ) :: j
    integer( kind = pin ) :: k
    integer( kind = pin ) :: iv( ntab )
    integer( kind = pin ) :: iy
    save iv , iy , idum2
    data idum2/123456789_pin/ , iv/ntab*0_pin/ , iy/0_pin/
    if ( idum .le. 0_pin ) then
       idum = max( -idum , 1_pin )
       idum2 = idum
       do j = ntab + 8_pin , 1_pin , - 1_pin
          k = idum / iq1
          idum = ia1 * ( idum - k * iq1 ) - k * ir1
          if ( idum .lt. 0_pin ) idum = idum + im1
          if ( j .le. ntab ) iv( j ) = idum
       enddo
       iy = iv( 1_pin )
    endif
    k = idum / iq1
    idum = ia1 * ( idum - k * iq1 ) - k * ir1
    if ( idum .lt. 0_pin ) idum = idum + im1
    k = idum2 / iq2
    idum2 = ia2 * ( idum2 - k * iq2 ) - k * ir2
    if ( idum2 .lt. 0_pin ) idum2 = idum2 + im2
    j = 1_pin + iy / ndiv
    iy = iv( j ) - idum2
    iv( j ) = idum
    if( iy .lt. 1_pin ) iy = iy + imm1
    random_uniform0d = min( am * iy , rnmx )
    return
  end function random_uniform0d
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_NORMALSTD0D
  !*
  !* PURPOSE: it returns a normally (standard) distributed deviate with zero mean and unit variance, using RANDOM_UNIFORM0D as th&
  !*          &e source of uniform deviates.
  !*
  !* INPUTS: none
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - RANDOM_NORMALSTD0D (real): random number with normal (standard) distribution
  !*
  !* CALLS:
  !*        - RANDOM_UNIFORM0D
  !* 
  !* COMMENTS: this is a translation of function gasdev.f of the Numerical Recipes adapted to generate random numbers in simple a&
  !*           &nd double precision and translated to free format in Fortran 90
  !*
  !*******************************************************************************************************************************
  function random_normalstd0d()
    implicit none
    real( kind = pre ) :: random_normalstd0d
    integer( kind = pin ) :: iset
    real( kind = pre ) :: fac
    real( kind = pre ) :: gset
    real( kind = pre ) :: rsq
    real( kind = pre ) :: v1
    real( kind = pre ) :: v2
    save iset , gset
    data iset/0_pin/
    if ( iset .eq. 0_pin ) then
1      v1 = 2.0_pre * random_uniform0d() - 1.0_pre     
       v2 = 2.0_pre * random_uniform0d() - 1.0_pre
       rsq = v1**2_pin + v2**2_pin
       if ( rsq .ge. 1.0_pre .or. rsq .eq. 0.0_pre ) goto 1
       fac = sqrt( -2.0_pre * log( rsq ) / rsq )
       gset = v1 * fac
       random_normalstd0d = v2 * fac
       iset = 1_pin
    else
       random_normalstd0d = gset
       iset = 0_pin
    endif
    return
  end function random_normalstd0d
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_NORMAL0D
  !*
  !* PURPOSE: it returns a normally distributed deviate with MU mean and SIGMA**2 variance, using RANDOM_NORMALSTD0D as the sourc&
  !*          &e of normal deviates
  !*
  !* INPUTS:
  !*         - MU (real): mean
  !*         - SIGMA (real): standard deviation
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - RANDOM_NORMAL0D (real): random number with normal distribution with average MU and variance SIGMA**2
  !*
  !* CALLS:
  !*        - RANDOM_NORMALSTD0D
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  function random_normal0d( mu , sigma )
    implicit none
    real( kind = pre ) , intent( in ) :: mu
    real( kind = pre ) , intent( in ) :: sigma
    real( kind = pre ) :: random_normal0d
    real( kind = pre ) :: sigmaaux
    sigmaaux = sigma
    if ( sigmaaux < 0.0_pre ) then
       write( * , '("error: the variance is negative.")' )
       stop
    endif
    random_normal0d = mu + sigmaaux * random_normalstd0d()
  end function random_normal0d
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_NORMAL1DMATRIX
  !*
  !* PURPOSE: it returns a set of random vectors with normal distribution with prescribed mean and prescribed covariance matrix 
  !*
  !* INPUTS:
  !*         - MEAN (real array of size SIZEMEAN): mean
  !*         - COVARIANCE (real array of size SIZEMEAN x SIZEMEAN): covariance matrix
  !*         - SIZEMEAN (integer): size of the random vectors to be generated
  !*         - NUMBERSAMPLES (integer): number of samples to generate
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - SAMPLES (real array of size SIZEMEAN x NUMBERSAMPLES): each column of SAMPLES stores a random vector
  !*
  !* CALLS:
  !*        - DSYEV: (from LAPACK)
  !*        - SSYEV: (from LAPACK)
  !*        - DGEMV: (from BLAS)
  !*        - SGEMV: (from BLAS)
  !*        - RANDOM_NORMAL0D
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine random_normal1dmatrix( mean , covariance , sizemean , numbersamples , samples )
    implicit none
    integer( kind = pin ) , intent( in ) :: sizemean
    integer( kind = pin ) , intent( in ) :: numbersamples
    real( kind = pre ) , dimension( sizemean ) , intent( in ) :: mean
    real( kind = pre ) , dimension( sizemean , sizemean ) , intent( in ) :: covariance
    real( kind = pre ) , dimension( sizemean , numbersamples ) , intent( out ) :: samples
    integer( kind = pin ) :: lwork
    integer( kind = pin ) :: info
    integer( kind = pin ) :: i
    integer( kind = pin ) :: s
    real( kind = pre ) , allocatable , dimension( : ) :: work
    real( kind = pre ) , allocatable , dimension( : ) :: sigmavalues
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux4
    allocate( sigmavalues( sizemean ) )
    allocate( matrixaux1( sizemean , sizemean ) )
    matrixaux1 = covariance
    allocate( work( 1_pin ) )
    !* jobz  = 'V'
    !* uplo  = 'U'
    !* n     = sizemean
    !* a     = matrixaux1
    !* lda   = sizemean
    !* w     = sigmavalues
    !* work  = work
    !* lwork = - 1_pin
    !* info  = info
    if ( pre == high ) then
       call dsyev( 'V' , 'U' , sizemean , matrixaux1 , sizemean , sigmavalues , work , - 1_pin , info )
    else if ( pre == low ) then
       call ssyev( 'V' , 'U' , sizemean , matrixaux1 , sizemean , sigmavalues , work , - 1_pin , info )
    end if
    if ( info /= 0_pin ) then
       print*,'info /= 0'
       print*,'info = ' , info
       stop
    end if
    lwork = int( work( 1_pin ) )
    deallocate( work )
    allocate( work( lwork ) )
    !* jobz  = 'V'
    !* uplo  = 'U'
    !* n     = sizemean
    !* a     = matrixaux1
    !* lda   = sizemean
    !* w     = sigmavalues
    !* work  = work
    !* lwork = lwork
    !* info  = info
    if ( pre == high ) then
       call dsyev( 'V' , 'U' , sizemean , matrixaux1 , sizemean , sigmavalues , work , lwork , info )
    else if ( pre == low ) then
       call ssyev( 'V' , 'U' , sizemean , matrixaux1 , sizemean , sigmavalues , work , lwork , info )
    end if
    if ( info /= 0_pin ) then
       print*,'info /= 0'
       print*,'info = ' , info
       stop
    end if
    deallocate( work )
    allocate( matrixaux4( sizemean , numbersamples ) )
    do i = 1_pin , sizemean
       do s = 1_pin , numbersamples
          matrixaux4( i , s ) = random_normal0d( 0.0_pre , sqrt( sigmavalues( i ) ) )
       end do
    end do
    deallocate( sigmavalues )
    do s = 1_pin , numbersamples
       samples( 1_pin : sizemean , s ) = mean
       !* trans = 'N'
       !* m     = sizemean
       !* n     = sizemean
       !* alpha = 1.0_pre
       !* a     = matrixaux1
       !* lda   = sizemean
       !* x     = matrixaux4( 1_pin : sizemean , s )
       !* incx  = 1_pin
       !* beta  = 1.0_pre
       !* y     = samples( 1_pin : sizemean , s )
       !* incy  = 1_pin
       if ( pre == high ) then
          call dgemv( 'N' , sizemean , sizemean , 1.0_pre , matrixaux1 , sizemean , matrixaux4( 1_pin : sizemean , s ) , 1_pin , &
               &1.0_pre , samples( 1_pin : sizemean , s ) , 1_pin )
       else if ( pre == low ) then
          call sgemv( 'N' , sizemean , sizemean , 1.0_pre , matrixaux1 , sizemean , matrixaux4( 1_pin : sizemean , s ) , 1_pin , &
               &1.0_pre , samples( 1_pin : sizemean , s ) , 1_pin )
       end if
    end do
    deallocate( matrixaux4 )
    deallocate( matrixaux1 )
  end subroutine random_normal1dmatrix
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_NORMAL1DSQRTMATRIX
  !*
  !* PURPOSE: it returns a set of random vectors with normal distribution with prescribed mean and prescribed reduced rank square&
  !*          & root covariance matrix
  !*
  !* INPUTS:
  !*         - MEAN (real array of size SIZEMEAN): mean
  !*         - SQRTCOV (real array of size SIZEMEAN x NUMBERMODES): square root of the covariance matrix with reduced rank limite&
  !*                                                                &d to NUMBERMODES columns
  !*         - SIZEMEAN (integer): size of the random vectors to be generated
  !*         - NUMBERMODES (integer): number of columns of the reduced rank square root of the covariance matrix
  !*         - NUMBERSAMPLES (integer): number of samples to generate
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - SAMPLES (real array of size SIZEMEAN x NUMBERSAMPLES): each column of SAMPLES stores a random vector
  !*
  !* CALLS:
  !*        - DGEMV: (from BLAS)
  !*        - SGEMV: (from BLAS)
  !*        - RANDOM_NORMALSTD0D
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine random_normal1dsqrtmatrix( mean , sqrtcov , sizemean , numbermodes , numbersamples , samples )
    implicit none
    integer( kind = pin ) , intent( in ) :: sizemean
    integer( kind = pin ) , intent( in ) :: numbermodes
    integer( kind = pin ) , intent( in ) :: numbersamples
    real( kind = pre ) , dimension( sizemean ) , intent( in ) :: mean
    real( kind = pre ) , dimension( sizemean , numbermodes ) , intent( in ) :: sqrtcov
    real( kind = pre ) , dimension( sizemean , numbersamples ) , intent( out ) :: samples
    integer( kind = pin ) :: s
    integer( kind = pin ) :: i
    real( kind = pre ) , allocatable , dimension( : ) :: vectoraux1
    allocate( vectoraux1( sizemean ) )
    do s = 1_pin , numbersamples
       samples( 1_pin : sizemean , s ) = mean
       do i = 1_pin , numbermodes
          vectoraux1( i ) = random_normalstd0d()
       end do
       !* trans = 'N'
       !* m     = sizemean
       !* n     = numbermodes
       !* alpha = 1.0_pre
       !* a     = sqrtcov
       !* lda   = sizemean
       !* x     = vectoraux1
       !* incx  = 1_pin
       !* beta  = 1.0_pre
       !* y     = samples( 1_pin : sizemean , s )
       !* incy  = 1_pin
       if ( pre == high ) then
          call dgemv( 'N' , sizemean , numbermodes , 1.0_pre , sqrtcov , sizemean , vectoraux1 , 1_pin , 1.0_pre , samples( 1_pin&
               & : sizemean , s ) , 1_pin )
       else if ( pre == low ) then
          call sgemv( 'N' , sizemean , numbermodes , 1.0_pre , sqrtcov , sizemean , vectoraux1 , 1_pin , 1.0_pre , samples( 1_pin&
               & : sizemean , s ) , 1_pin )
       end if
    end do
    deallocate( vectoraux1 )
  end subroutine random_normal1dsqrtmatrix
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_MEANESTIMATOR1D
  !*
  !* PURPOSE: it returns the mean of a set of vectors stored by columns in a matrix
  !*
  !* INPUTS:
  !*         - SIZE (integer): size of the vectors to which we will compute the mean
  !*         - NUMBERSAMPLES (integer): number of vectors
  !*         - SAMPLES (real array of size SIZE x NUMBERSAMPLES): each column of SAMPLES is a vector
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - MEAN (real array of size SIZEMEAN): mean of the vectors stored in SAMPLES
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine random_meanestimator1d( size , numbersamples , samples , mean )
    implicit none
    integer( kind = pin ) , intent( in ) :: size
    integer( kind = pin ) , intent( in ) :: numbersamples
    real( kind = pre ) , dimension( size , numbersamples ) , intent( in ) :: samples
    real( kind = pre ) , dimension( size ) , intent( out ) :: mean
    integer( kind = pin ) :: i
    do i = 1_pin , size
       mean( i ) = sum( samples( i , 1_pin : numbersamples ) ) / numbersamples
    end do
  end subroutine random_meanestimator1d
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_COVARIANCEESTIMATOR1D
  !*
  !* PURPOSE: it returns the covariance matrix of a set of vectors stored by columns in a matrix
  !*
  !* INPUTS:
  !*         - SIZE (integer): size of the vectors to which we will compute the covariance matrix
  !*         - NUMBERSAMPLES (integer): number of random vectors
  !*         - SAMPLES (real array of size SIZE x NUMBERSAMPLES): each column of SAMPLES is a vector
  !*         - MEAN (real array of size SIZE): mean of the vectors stored by columns in SAMPLES
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - COVARIANCE (real array of size SIZE x SIZE): covariance matrix of the vectors stored by columns in the matrix SAM&
  !*                                                         &PLES
  !*
  !* CALLS:
  !*        - DGEMM: (from BLAS)
  !*        - SGEMM: (from BLAS)
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine random_covarianceestimator1d( size , numbersamples , samples , mean , covariance )
    implicit none
    integer( kind = pin ) , intent( in ) :: size
    integer( kind = pin ) , intent( in ) :: numbersamples
    real( kind = pre ) , dimension( size , numbersamples ) , intent( in ) :: samples
    real( kind = pre ) , dimension( size ) , intent( in ) :: mean
    real( kind = pre ) , dimension( size , size ) , intent( out ) :: covariance
    real( kind = pre ) , allocatable , dimension( : , : ) :: matrixaux1
    integer( kind = pin ) :: js
    allocate( matrixaux1( size , numbersamples ) )
    do js = 1_pin , numbersamples
       matrixaux1( 1_pin : size , js ) = samples( 1_pin : size , js ) - mean
    end do
    !* transa = 'N'
    !* transb = 'T'
    !* m      = size
    !* n      = size
    !* k      = numbersamples
    !* alpha  = 1.0_pre / ( numbersamples - 1.0_pre )
    !* a      = matrixaux1
    !* lda    = size
    !* b      = matrixaux1
    !* ldb    = size
    !* beta   = 0.0_pre
    !* c      = covariance
    !* ldc    = size
    if ( pre == high ) then
       call dgemm ( 'N' , 'T' , size , size , numbersamples , 1.0_pre / ( numbersamples - 1.0_pre ) , matrixaux1 , size , matrixa&
            &ux1 , size , 0.0_pre , covariance , size )
    else if ( pre == low ) then
       call sgemm ( 'N' , 'T' , size , size , numbersamples , 1.0_pre / ( numbersamples - 1.0_pre ) , matrixaux1 , size , matrixa&
            &ux1 , size , 0.0_pre , covariance , size )
    end if
    deallocate( matrixaux1 )
  end subroutine random_covarianceestimator1d
  !*******************************************************************************************************************************
  !* 
  !* FUNCTION/SUBROUTINE: RANDOM_SQRTCOVESTIMATOR1D
  !*
  !* PURPOSE: it returns the reduced rank square root of the covariance matrix of a set of vectors stored by columns in a matrix
  !*
  !* INPUTS:
  !*         - SIZE (integer): size of the vectors to which we will compute the reduced rank square root of the covariance matrix
  !*         - NUMBERSAMPLES (integer): number of vectors
  !*         - SAMPLES (real array of size SIZE x NUMBERSAMPLES): each column of SAMPLES is a vector
  !*         - MEAN (real array of size SIZE): mean of the vectors stored by columns in SAMPLES
  !*
  !* INPUTS/OUTPUTS: none
  !*
  !* OUTPUTS:
  !*          - SQRTCOVARIANCE (real array of size SIZE x NUMBERSAMPLES): reduced rank square root covariance matrix of the vecto&
  !*                                                                      &rs stored by columns in the matrix SAMPLES
  !*
  !* CALLS: none
  !* 
  !* COMMENTS: none
  !*
  !*******************************************************************************************************************************
  subroutine random_sqrtcovestimator1d( size , numbersamples , samples , mean , sqrtcovariance )
    implicit none
    integer( kind = pin ) , intent( in ) :: size
    integer( kind = pin ) , intent( in ) :: numbersamples
    real( kind = pre ) , dimension( size , numbersamples ) , intent( in ) :: samples
    real( kind = pre ) , dimension( size ) , intent( in ) :: mean
    real( kind = pre ) , dimension( size , numbersamples ) , intent( out ) :: sqrtcovariance
    integer( kind = pin ) :: s
    do s = 1_pin , numbersamples
       sqrtcovariance( 1_pin : size , s ) = ( samples( 1_pin : size , s ) - mean( 1_pin : size ) ) / sqrt( ( numbersamples - 1_pi&
            &n ) * 1.0_pre )
    end do
  end subroutine random_sqrtcovestimator1d
  !*
end module random
