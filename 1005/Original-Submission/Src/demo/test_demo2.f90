SUBROUTINE F(testcase, u, v)
  ! A test routine for the demo2_likeli function. Similar in structure to those
  ! used for testing the BLAS-RMD package.
  use rmd_dtesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  double precision, intent(in) :: u(*)
  double precision, intent(out) :: v(*)

  double precision, allocatable :: Sigma(:,:), A(:,:), X(:,:), &
      Sigmaa(:,:), Aa(:,:)
  integer n, r, i, dimu
  double precision :: demo2_likeli

  n = testcase % integers(1)
  r = testcase % integers(2)
  dimu = testcase % integers(3)
  
  allocate(Sigma(r,r), Sigmaa(r,r), Aa(r,r))

  A = reshape(u(1:r*r),[r,r])
  call rmd_dtrishape_v2m('l', r, u(r*r+1:dimu),Sigma, r)

  do i = 1,r
    Sigma(i,i+1:r) = Sigma(i+1:r,i)
  end do

  X = reshape(testcase % reals, [r,n])
  Sigmaa = 0
  Aa = 0

  v(1) = demo2_likeli(n,r, Sigma,A,X, Sigmaa,Aa)
END SUBROUTINE F

SUBROUTINE F_RMD(testcase, u, v, ua, va)
  use rmd_dtesttools
  implicit none
  type(fdata), intent(inout) :: testcase
  double precision, intent(in) :: u(*), v(*), va(*)
  double precision, intent(inout) :: ua(*)
  double precision, allocatable :: Sigma(:,:), A(:,:), X(:,:), Sigmaa(:,:), Aa(:,:)
  double precision usev, ll, demo2_likeli
  integer n, r, i, dimu
  usev = v(1); usev = usev
  usev = va(1); usev = usev

  n = testcase % integers(1)
  r = testcase % integers(2)
  dimu = testcase % integers(3)
  
  allocate(Sigma(r,r), Sigmaa(r,r), Aa(r,r))

  A = reshape(u(1:r*r),[r,r])

  call rmd_dtrishape_v2m('l', r, u(r*r+1:dimu),Sigma, r) 
  do i = 1,r
    Sigma(i,i+1:r) = Sigma(i+1:r,i)
  end do
  X = reshape(testcase % reals,[r,n])
  Sigmaa = 0
  Aa = 0
  ll = demo2_likeli(n,r, Sigma,A,X, Sigmaa,Aa)
  ll = ll
  ua(1:r*r) = reshape(Aa,[r*r])
  call rmd_dtrishape_m2v('l', r, ua(r*r+1:dimu), Sigmaa, r) 

end subroutine F_rmd

PROGRAM TEST_DEMO2
  use rmd_dtesttools
  implicit none
  type(fdata) :: testcase
  integer n, r, i, j, dimu
  double precision, allocatable :: &
    A(:,:), Sigma(:,:), X(:,:), Aa(:,:), Sigmaa(:,:), u(:)
  double precision tol
  character(30) :: filename
  procedure(F_interf) :: F
  procedure(F_rmd_interf) :: F_rmd
  call get_command_argument(1, filename)
  tol = 1.0d-6

  open(unit=10, file=filename, status='old')
  read(10,*) r
  read(10,*) n

  dimu = r*r + (r*(r+1))/2

  allocate(A(r,r), Sigma(r,r), X(r,n), Aa(r,r), Sigmaa(r,r), u(dimu))
  Sigmaa = 0d0
  Aa = 0d0
  
  read(10,*) ((X(i,j), j=1,n), i=1,r)
  read(10,*) Sigma
  read(10,*) ((A(i,j), j=1,r), i=1,r)

  u(1:r*r) = reshape(A, [r*r])
  u(r*r+1:dimu) = rmd_d_vech('l', Sigma)
  !call rmd_dtrishape_m2v('l', r, u(r*r+1),Sigma, r) 

  testcase = fdata(reshape(X,[n*r]),[n,r,dimu],[character::])

  call rmd_dtestf(F, F_rmd, u, dimu, 1, tol, testcase)

  close(10)
END PROGRAM TEST_DEMO2

