!-----------------------------------------------------------!
! Compute the genus of 3D objects. Serial version. VER 1.2  !
!                                                           !
! Reference:                                                !
! "An efficient algorithm to compute the genus of           !
! discrete 3D objects and applications to turbulent flows"  !
! Lozano-Duran et al.                                       !
! ACM TOMS                                                  !
!                                                           !
! Compile:                                                  !
! make genus                                                !
!                                                           !
! Tested with:                                              !
! Intel Fortran Studio XE 2016 16.0.0 20150815              !
!                                                           !
! Author:                                                   !
! Adrian Lozano-Duran October 2015                          !
!-----------------------------------------------------------!
Program main

  Use iso_fortran_env, Only : output_unit,error_unit,int8,int32,int64
  
  ! Prevent implicit typing
  Implicit None

  ! define variables
  Character (Len=:), Allocatable :: filename
  Integer   (int8) , Allocatable, Dimension (:, :, :) :: a
  Integer   (int32) :: nx, ny, nz
  Integer   (int64) :: genus, faces, edges, vertices
        
  ! header
  Write (output_unit, '(A)') '!----------------------------------!'
  Write (output_unit, '(A)') '!       WELCOME TO GENUSLAND!      !'
  Write (output_unit, '(A)') '!----------------------------------!'
  
  ! input
  filename = command_argument_string(argument_position=1)
  if ( Len(filename)<=0 ) Then 
     input : Block 
       Character(Len=:), Allocatable :: command_name
       command_name = command_argument_string(argument_position=0)
       Write(error_unit,'(A)')  "Usage:"
       Write(error_unit,'(3A)') "  ",command_name," <file-name> "
       Write(error_unit,'(3A)') "  Example: ",command_name," ../data/random.128.dat"
       error Stop 
     End Block input
  End If
  
  ! read data
  Write (output_unit, '(A,A)') '  reading: ', filename
  read_data : Block 
    Integer(int32), Parameter :: file_unit=1
    Open  (Unit=file_unit, File=filename, Action='read', Form='unformatted', &
         Access='stream', Status='old', Convert='big_endian')
    Read  (file_unit) nx, ny, nz
    Write (output_unit, '(A,I4,A,I4,A,I4)') '  nx =', nx, ' ny =', ny, ' nz =', nz
    Allocate (a(nx,ny,nz))
    Read  (file_unit) a
    Close (file_unit)
  end Block read_data
  
  ! compute results 
  compute : Block 
    Use genusmod, Only : compute_genus        
    Call compute_genus(a, genus, faces, edges, vertices)
  End Block compute

  ! output
  output : Block
    Write (output_unit, '(A12,I10)') ' Faces:    ', faces
    Write (output_unit, '(A12,I10)') ' Edges:    ', edges
    Write (output_unit, '(A12,I10)') ' Vertices: ', vertices
    Write (output_unit, '(A12,I10)') ' Genus:    ', genus
  End Block output
  
  Deallocate (a)
  
Contains
  Function command_argument_string(argument_position) Result(string)
    Integer (int32), Intent (In) :: argument_position
    Integer (int32) :: argument_length
    Character (Len=:), Allocatable :: string
    Call get_command_argument(number=argument_position,length=argument_length)
    string = Repeat(" ",Ncopies=argument_length)
    Call get_command_argument(number=argument_position,Value=string)
  End Function command_argument_string

End Program main

