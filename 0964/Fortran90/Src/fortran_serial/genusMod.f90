!------------------------------------------!
!               Shared data                !
!------------------------------------------! 
Module genusmod

  Use iso_fortran_env, Only: error_unit, int8, int32, int64
  
  ! prevent impilicit typing
  Implicit None
  
  ! only compute_genus is expected to be called by the user
  Private 
  Public  :: compute_genus
  
Contains
!------------------------------------------!
! Subroutine: compute_genus                !
!             Computes number of faces,    ! 
!             edges,vertices, genus, and   ! 
!             Euler characteristic         !
!------------------------------------------!
  Subroutine compute_genus(a, g, f, e, v)
    
    Integer (int8),  Intent (In)  :: a(:, :, :)
    Integer (int64), Intent (Out) :: g, f, e, v
    
    Integer (int8), Dimension (2,2,2) :: cube1
    Integer (int8), Dimension (2,2)   :: cube2
    Integer (int32) :: i, j, k, nx, ny, nz
    Integer (int64) :: x
    
    nx = Size(a,1)
    ny = Size(a,2)
    nz = Size(a,3)
    
    f = 0 ! faces
    e = 0 ! edges
    v = 0 ! vertices
    
    Do k = 2, nz - 1
       Do j = 2, ny - 1
          Do i = 2, nx - 1 
             ! compute faces
             If ( a(i,j,k)==1 ) f = f + Count([a(i-1,j,k), a(i+1,j,k), a(i,j-1,k), a(i,j+1,k), a(i,j,k-1), a(i,j,k+1)]==0)
             
             ! compute vertices
             cube1 = a(i-1:i, j-1:j, k-1:k)
             If ( Any(cube1==0) .And. Any(cube1==1) ) v = v + getconnected1(cube1)
             
             ! compute edges
             cube2 = a(i-1:i, j-1:j, k)
             If ( Any(cube2==0) .And. Any(cube2==1) ) e = e + getconnected2(cube2)
             
             cube2 = a(i-1:i, j, k-1:k)
             If ( Any(cube2==0) .And. Any(cube2==1) ) e = e + getconnected2(cube2)
             
             cube2 = a(i, j-1:j, k-1:k)
             If ( Any(cube2==0) .And. Any(cube2==1) ) e = e + getconnected2(cube2)
             
          End Do
       End Do
    End Do
    
    x = f - e + v
    g = ( 2 - x )/2
    
    ! small sanity check
    If ( Mod(x,2_int64)/=0 ) Write (error_unit,'(A)') 'G not integer!'
    
  End Subroutine compute_genus
  
!--------------------------------------------!
! Function: getconnected1                    !
!           Computes increment of vertices   !
!--------------------------------------------!
  Pure Function getconnected1(cube1_) Result(dv)
    
    Integer (int8), Intent (In), Dimension(:,:,:) :: cube1_
    Integer (int64) :: dv
    
    Integer (int8), Dimension(2,2,2) :: cube1_aux
    Integer (int32) :: newgroup, i, j, k
    
    cube1_aux = cube1_
    dv        = 0
    
    ! degenerated case
    If (Sum(cube1_aux)==6 .And. ( (cube1_aux(1,1,1)==0 .And. cube1_aux(2,2,2)==0) &
                             .Or. (cube1_aux(1,1,2)==0 .And. cube1_aux(2,2,1)==0) &
                             .Or. (cube1_aux(1,2,1)==0 .And. cube1_aux(2,1,2)==0) &
                             .Or. (cube1_aux(1,2,2)==0 .And. cube1_aux(2,1,1)==0) ) ) Then
       dv = 2
    Else
       ! regular cases   
       Do k = 1, 2
          Do j = 1, 2
             Do i = 1, 2
                If (cube1_aux(i,j,k)/=1) Cycle
                newgroup = 1
                Call makecluster(i, j, k, newgroup, dv, cube1_aux)
             End Do
          End Do
       End Do
       
    End If
    
  Contains
!------------------------------------------!
! Subroutine: makecluster                  !
!             Recursive labeling           !
!------------------------------------------!
    Pure Recursive Subroutine makecluster(i, j, k, newgroup, dv, cube1_)
      
      Integer (int32), Intent (In)    :: i, j, k
      Integer (int32), Intent (Inout) :: newgroup
      Integer (int64), Intent (Inout) :: dv
      Integer (int8),  Intent (Inout), Dimension(:,:,:) :: cube1_
      
      ! if checked, skip  
      If ( cube1_(i,j,k)/=1 ) Return
      
      ! newgroup?                 
      If ( newgroup==1 ) Then
         dv        = dv + 1
         newgroup = 0
      End If
      
      ! add to group   
      cube1_(i, j, k) = -1
      
      ! compute connections   
      ! x
      If (i<=1) Call makecluster(i+1, j, k, newgroup, dv, cube1_)          
      If (i>=2) Call makecluster(i-1, j, k, newgroup, dv, cube1_)
      
      ! y
      If (j<=1) Call makecluster(i, j+1, k, newgroup, dv, cube1_)         
      If (j>=2) Call makecluster(i, j-1, k, newgroup, dv, cube1_)
      
      ! z
      If (k<=1) Call makecluster(i, j, k+1, newgroup, dv, cube1_)          
      If (k>=2) Call makecluster(i, j, k-1, newgroup, dv, cube1_)
      
    End Subroutine makecluster
    
  End Function getconnected1
      
!------------------------------------------!
! Function: getconnected2                  !
!           Computes increment of edges    !
!------------------------------------------!
  Pure Function getconnected2(cube2_) Result(de)
    
    Integer (int8), Dimension (2,2), Intent(in) :: cube2_
    Integer (int64) :: de
    
    de = Merge(2,1, &
         (     (cube2_(1,1)==1 .And. cube2_(2,2)==1 .And. cube2_(1,2)==0 .And. cube2_(2,1)==0) &
          .Or. (cube2_(1,1)==0 .And. cube2_(2,2)==0 .And. cube2_(1,2)==1 .And. cube2_(2,1)==1) &
         ) )        
  End Function getconnected2
  
End Module genusmod
