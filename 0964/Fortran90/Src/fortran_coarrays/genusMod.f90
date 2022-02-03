!------------------------------------------!
!               Shared data                !
!------------------------------------------! 
Module genusmod

  Use iso_fortran_env, Only: output_unit, error_unit, int8, int32, int64

  ! prevent impilicit typing
  Implicit None

  ! define variables
  Integer (int32)                 :: nx, ny, nz
  Integer (int32), Codimension[*] :: kfirst, klast
  Integer (int8),  Allocatable    :: a(:,:,:)
  
  ! only read_data and compute_genus are expected to be called by the user
  Private 
  Public  :: read_data, compute_genus
  
Contains

!------------------------------------------!
! Subroutine: read_data                    !
!             Reads and distributes data   !
!------------------------------------------!
  Subroutine read_data(filename)

    Character(*), Intent (In) :: filename

    Integer, Parameter :: file_unit=1
    Integer (int32)    :: dk, i
    Integer (int64)    :: position
      
    ! read data size
    critical
      Open (Unit=file_unit, File=filename, Action='read', Form='unformatted', &
            Access='stream', Status='old', Convert='big_endian')
      Read (file_unit) nx, ny, nz
    end critical

    ! print info
    If ( this_image()==1 ) Then       
       Write (output_unit, '(A,A)') '  reading: ', filename              
       Write (output_unit, '(A,I4,A,I4,A,I4)') '  nx =', nx, ' ny =', ny, ' nz =', nz
    End If

    ! domain decomposition
    dk        = nint( real(nz)/real(num_images()) )
    kfirst[1] = 1
    klast [1] = kfirst[1] + dk 
    Do i = 2, num_images()
       kfirst[i] = klast [i-1] - 1
       klast [i] = kfirst[i]   + dk 
    End Do
    klast[num_images()] = nz
    
    Allocate ( a(nx, ny, kfirst:klast) )

    ! each image reads its data
    position = int(3*4,int64) + int(kfirst-1,int64)*int(nx,int64)*int(ny,int64) + 1_int64
    Critical
      Read  (file_unit,Pos=position) a
      Close (file_unit)                        
    End Critical

    ! sync
    Sync All

  end Subroutine read_data
  
!------------------------------------------!
! Subroutine: compute_genus                !
!             Computes number of faces,    ! 
!             edges,vertices, genus, and   ! 
!             Euler characteristic         !
!------------------------------------------!
  Subroutine compute_genus(genus, faces, edges, vertices)

    Integer (int64), Intent (Out) :: genus, faces, edges, vertices
    
    Integer (int64), Codimension[*], save :: f, e, v
    Integer (int8),  Dimension (2,2,2)    :: cube1
    Integer (int8),  Dimension (2,2)      :: cube2
    Integer (int32) :: i, j, k
    Integer (int64) :: x
         
    f = 0 ! faces
    e = 0 ! edges
    v = 0 ! vertices

    Do k = kfirst+1, klast-1
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

    ! sum data from all images to image 1
    Sync All      
    If ( this_image()==1 ) Then
       Do i = 2, num_images()
          f = f + f[i]
          e = e + e[i]
          v = v + v[i]          
       End Do       
    End If
    Sync All
    
    ! distribute to all images
    faces    = f[1]
    edges    = e[1]
    vertices = v[1]        

    ! Euler charcateristic and genus
    x     = faces - edges + vertices
    genus = ( 2 - x )/2
    
    ! small sanity check
    If ( Mod(x,2_int64)/=0 ) Write (error_unit,'(A)') 'G not integer!'

    ! deallocate data
    Deallocate (a)
        
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
    
    Integer (int8), Dimension (2,2), Intent (In) :: cube2_
    Integer (int64) :: de
    
    de = Merge(2,1, &
         (     (cube2_(1,1)==1 .And. cube2_(2,2)==1 .And. cube2_(1,2)==0 .And. cube2_(2,1)==0) &
          .Or. (cube2_(1,1)==0 .And. cube2_(2,2)==0 .And. cube2_(1,2)==1 .And. cube2_(2,1)==1) &
         ) )        
  End Function getconnected2
  
End Module genusmod
