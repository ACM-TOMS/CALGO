!-----------------------------------------------------------!
! Compute the genus of 3D objects. Coarrays version.        !
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

  Use iso_fortran_env, Only: output_unit,error_unit,int8,int32,int64
  
  ! prevent implicit typing
  Implicit None

  ! define variables
  Character (Len=:), Allocatable :: filename
  Integer   (int64) :: genus, faces, edges, vertices
  
  ! header
  If ( this_image()==1 ) Then
     Write (output_unit, '(A)') '!----------------------------------!'
     Write (output_unit, '(A)') '!       WELCOME TO GENUSLAND!      !'
     Write (output_unit, '(A)') '!----------------------------------!'
     Write (output_unit, '(A,I3)') 'Number of processors: ',num_images()
  End If
     
  ! input
  input : Block 
    Character(len=:), Allocatable :: command_name    
    filename = command_argument_string(argument_position=1)
    If ( Len(filename)<=0 ) Then 
       command_name = command_argument_string(argument_position=0)
       Write (error_unit,'(A)')  "Usage:"
       Write (error_unit,'(3A)') "  ",command_name," <file-name> "
       Write (error_unit,'(3A)') "  Example: ",command_name," ../data/random.128.dat"
       error Stop 
    End If
  End Block input
  
  ! read and distribute, data is store in 'a'
  read_dist : Block
    Use genusmod, Only: read_data
    call read_data(filename)
  End Block read_dist
    
  ! compute results 
  compute : Block 
    Use genusmod, Only: compute_genus    
    Call compute_genus(genus, faces, edges, vertices)
  End Block compute

  ! output results
  output : Block
    If ( this_image()==1 ) Then
       Write (output_unit, '(A12,I10)') ' Faces:    ', faces
       Write (output_unit, '(A12,I10)') ' Edges:    ', edges
       Write (output_unit, '(A12,I10)') ' Vertices: ', vertices
       Write (output_unit, '(A12,I10)') ' Genus:    ', genus
    End if
  End Block output
  
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

