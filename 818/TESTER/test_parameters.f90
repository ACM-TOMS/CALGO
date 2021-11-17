module mod_test_parameters
use blas_sparse_namedconstants
use properties
! *** The test matrices
      integer, parameter :: MATRIX_A = 1
      integer, parameter :: MATRIX_A_SU = 2
      integer, parameter :: MATRIX_A_SL = 3
      integer, parameter :: MATRIX_T = 4
      integer, parameter :: MATRIX_U = 5
      integer, parameter :: MATRIX_B = 6
      integer, parameter :: MATRIX_X = 7
      integer, parameter :: MATRIX_Y = 8
      integer, parameter :: MATRIX_Y_SU = 9
      integer, parameter :: MATRIX_Y_SL = 10
      integer, parameter :: MATRIX_M = 11
      integer, parameter :: MATRIX_N = 12
      integer, parameter :: MATRIX_Z = 13
      integer, parameter :: MATRIX_C = 14
      integer, parameter :: MATRIX_I = 15
      integer, parameter :: MATRIX_J = 16
! *** Which functionality of MV/MM mult./Tri solve should be called ?
      integer, parameter :: O_MATRIX = 0 !original
      integer, parameter :: T_MATRIX = 1 !transpose
      integer, parameter :: H_MATRIX = 2 !hermit


! *** Level 1 routines
real(kind=dp),dimension(3)::x_dot=(/1.,3.,5./)
real(kind=dp),dimension(6)::y_dot=(/2.,4.,6.,8.,10.,12./)
integer ,dimension(3),parameter ::indx_dot=(/2,3,4/)

real(kind=dp),dimension(3)::x_axpy=(/51.,53.,55./)
real(kind=dp),dimension(6)::y_axpy=(/52.,54.,56.,58.,60.,62./)
integer ,dimension(3)::indx_axpy=(/2,3,4/)
real(kind=dp),parameter::alpha=2.

real(kind=dp),dimension(3)::x_ga=(/31.,33.,35./)
real(kind=dp),dimension(6)::y_ga=(/32.,34.,36.,38.,40.,42./)
integer ,dimension(3)::indx_ga=(/2,3,4/)

real(kind=dp),dimension(3)::x_gz=(/11.,13.,15./)
real(kind=dp),dimension(6)::y_gz=(/12.,14.,16.,18.,20.,22./)
integer ,dimension(3) ::indx_gz=(/2,3,4/)

real(kind=dp),dimension(3)::x_sc=(/41.,43.,45./)
real(kind=dp),dimension(6)::y_sc=(/42.,44.,46.,48.,50.,52./)
integer ,dimension(3) ::indx_sc=(/2,3,4/)

complex(kind=dp),dimension(3)::xc_dot=(/(1.,1.),(3.,1.),(5.,1.)/)
complex(kind=dp),dimension(6)::yc_dot=(/(2.,1.),(4.,1.),(6.,1.),(8.,1.),(10.,1.),(12.,1.)/)
integer ,dimension(3),parameter ::indxc_dot=(/2,3,4/)

complex(kind=dp),dimension(3)::xc_axpy=(/(51.,1.),(53.,1.),(55.,1.)/)
complex(kind=dp),dimension(6)::yc_axpy=(/(52.,1.),(54.,1.),(56.,1.),(58.,1.),(60.,1.),(62.,1.)/)
integer ,dimension(3)::indxc_axpy=(/2,3,4/)
complex(kind=dp),parameter::alphac=(2.,1)

complex(kind=dp),dimension(3)::xc_ga=(/(31.,1.),(33.,1.),(35.,1.)/)
complex(kind=dp),dimension(6)::yc_ga=(/(32.,1.),(34.,1.),(36.,1.),(38.,1.),(40.,1.),(42.,1.)/)
integer ,dimension(3)::indxc_ga=(/2,3,4/)

complex(kind=dp),dimension(3)::xc_gz=(/(11.,1.),(13.,1.),(15.,1.)/)
complex(kind=dp),dimension(6)::yc_gz=(/(12.,1.),(14.,1.),(16.,1.),(18.,1.),(20.,1.),(22.,1.)/)
integer ,dimension(3) ::indxc_gz=(/2,3,4/)

complex(kind=dp),dimension(3)::xc_sc=(/(41.,1.),(43.,1.),(45.,1.)/)
complex(kind=dp),dimension(6)::yc_sc=(/(42.,1.),(44.,1.),(46.,1.),(48.,1.),(50.,1.),(52.,1.)/)
integer ,dimension(3) ::indxc_sc=(/2,3,4/)

! *** Level 2,3 and handle management routines
 !-----------------------------------------------------------------
      ! Test - Matrices :
      !
      !    /11  0 13 14  0\     /11  0 13 14  0\     /11  0 31 41 51\
      !    | 0  0 23 24  0|     | 0  0 23 24  0|     | 0  0 32 42 52|
      ! A= |31 32 33 34  0|A_SU=|13 23 33 34  0|A_SL=|31 32 33  0  0|
      !    | 0 42  0 44  0|     |14 24  0 44  0|     |41 42  0 44  0|
      !    \51 52  0  0 55/     \ 0  0  0  0 55/     \51 52  0  0 55/
      !
      !    / 1  4  7 10  13\  
      ! B= | 2  5  8 11  14| 
      !    \ 3  6  9 12  15/ 
      !
      !    / 1 i\
      ! C= |    |
      !    \-i 1/       
      !
      !    / 1 i\
      ! D= |    |
      !    \-i 1/       
      !
      ! I = A (integer storage)
      !
      ! J = T (integer storage)
      !
      !    /  1  2 |  4  7 \     /  1  0 |      \
      !    |  0  3 |  5  8 |     |  2  5 |      |
      ! M= |-------+-------|  N= |-------+------|
      !    |       |  6  9 |     |  3  6 | 8  0 |
      !    \       |  0 10 /     \  4  7 | 9 10 /
      !
      !    / 1  1  1  1  1\     / 1            \
      !    |    1  1  1  1|     | 1  1         |
      ! T= |       1  1  1|   U=| 1  1  1      |
      !    |          1   |     | 1  1  1  1   |
      !    \             1/     \ 1  1  1     1/
      !
      !    /11 12 |  0  0 | 15 16|  0  0\
      !    |21 22 |  0  0 | 25 26|  0  0| 
      !    |------+-------+------+------|
      ! X= | 0  0 | 33  0 | 35 36|  0  0|
      !    | 0  0 | 43 44 | 45 46|  0  0|
      !    |------+-------+------+------|
      !    |51 52 |  0  0 |  0  0|  0  0|
      !    \61 62 |  0  0 |  0  0|  0  0/  
      !
      !    / 1  3 |  5  7 |  0  0\
      !    | 3  4 |  6  8 |  0  0| 
      !    |------+-------+------|
      ! Y= | 9 11 | 13 15 |  0  0|
      !    |10 12 | 15 16 |  0  0|
      !    |------+-------+------|
      !    | 0  0 |  0  0 |  0  0|
      !    \ 0  0 |  0  0 |  0  0/  
      !
      !    / 4  2 |  0  0  0 |  1 |  0  0  0 | -1  1 \
      !    | 1  5 |  0  0  0 |  2 |  0  0  0 |  0 -1 | 
      !    |------+----------+----+----------+-------|
      !    | 0  0 |  6  1  2 |  2 |  0  0  0 |  0  0 |
      !    | 1  5 |  2  7  1 |  0 |  0  0  0 |  0  0 | 
      !    | 1  5 | -1  2  9 |  3 |  0  0  0 |  0  0 | 
      !    |------+----------+----+----------+-------|
      ! Z= | 2  1 |  3  4  5 | 10 |  4  3  2 |  0  0 |
      !    |------+----------+----+----------+-------|
      !    | 0  0 |  0  0  0 |  4 | 13  4  2 |  0  0 |
      !    | 0  0 |  0  0  0 |  3 |  3 11  3 |  0  0 |
      !    | 0  0 |  0  0  0 |  0 |  2  0  7 |  0  0 |
      !    |------+----------+----+----------+-------|
      !    | 8  4 |  0  0  0 |  0 |  0  0  0 | 25  3 |
      !    \-2  3 |  0  0  0 |  0 |  0  0  0 |  8 12 / 
      !
      !-----------------------------------------------------------------


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      !    /11  0 13 14  0\    
      !    | 0  0 23 24  0|    
      ! A= |31 32 33 34  0|
      !    | 0 42  0 44  0|     
      !    \51 52  0  0 55/     
      !
   real(kind=dp), dimension(14)  ::&
 A_val=(/ 11.,51.,31.,32.,34.,52.,13.,23.,33.,14.,24.,42.,55.,44. /)
      integer, dimension(14)  ::&
 A_indx= (/ 1,5,3,3,3,5,1,2,3,1,2,4,5,4 /)
      integer, dimension(14)  ::&
 A_jndx = (/ 1,1,1,2,4,2,3,3,3,4,4,2,5,4 /)
     integer ::A_m=5
     integer ::A_n=5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      !      /11  0 31  0 51\
      !      | 0  0 32 42 52|
      ! A_SL=|31 32 33  0  0|
      !      | 0 42  0 44  0|
      !      \51 52  0  0 55/
      !   

   real(kind=dp),dimension(5) ::A_SL_1=(/ 11., 0.,31., 0.,51./)
   integer,dimension(5) ::A_SL_1_INDX=(/1,2,3,4,5/)

   real(kind=dp),dimension(5) ::A_SL_2=(/0.,  0.,32., 42., 52./)
   integer,dimension(5) ::A_SL_2_INDX=(/1,2,3,4,5/)

   real(kind=dp),dimension(5) ::A_SL_3=(/31.,32.,33., 0.,0./)
   integer,dimension(5) ::A_SL_3_INDX=(/1,2,3,4,5/)

   real(kind=dp),dimension(5) ::A_SL_4=(/0.,42.,0., 44.,0./)
   integer,dimension(5) ::A_SL_4_INDX=(/1,2,3,4,5/)

   real(kind=dp),dimension(5) ::A_SL_5=(/51.,52.,0., 0.,55./)
   integer,dimension(5) ::A_SL_5_INDX=(/1,2,3,4,5/)

   integer ::A_SL_m=5
   integer,parameter::A_SL_n=5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   !     /11  0 13 14  0\  
   !     | 0  0 23 24  0| 
   !A_SU=|13 23 33 34  0|
   !     |14 24  0 44  0|    
   !     \ 0  0  0  0 55/    


      real(kind=dp),dimension(5,5),parameter::&
           A_SU=reshape((/11.,0.,0.,0.,0.,0.,0.,0.,0.,0.,13.,23.,33.,0.,0.,14.,24.,34.,44.,0.,0.,0.,0.,0.,55./),shape=(/5,5/))
      integer,dimension(5) ::A_SU_indx=(/1,2,3,4,5/)
      integer,dimension(5) ::A_SU_jndx=(/1,2,3,4,5/)
      integer,parameter::A_SU_m=5
      integer,parameter::A_SU_n=5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      !    / 1  4  7 10  13\  
      ! B= | 2  5  8 11  14| 
      !    \ 3  6  9 12  15/ 
      !

   real(kind=dp),dimension(5) ::B_1=(/ 1.,4.,7.,10.,13./)
   integer,dimension(5) ::B_1_JNDX=(/1,2,3,4,5/)

   real(kind=dp),dimension(5) ::B_2=(/2.,5.,8.,11.,14./)
   integer,dimension(5) ::B_2_JNDX=(/1,2,3,4,5/)

   real(kind=dp),dimension(5) ::B_3=(/3.,  6.,9.,12.,15./)
   integer,dimension(5) ::B_3_JNDX=(/1,2,3,4,5/)
   integer ::B_m=3
   integer ::B_n=5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    / 1 0 i\
      !C= |  0 0 0 |
      !    \-i 0 1/       

complex(kind=dp),dimension(4),parameter ::VAL_C=(/(1.,0.),(0.,1.),(0.,-1.),(1.,0.)/)
integer ,dimension(4)  ::INDX_C=(/1,1,3,3/)
integer ,dimension(4)  ::JNDX_C=(/1,3,1,3/)
integer,parameter::C_m=3
integer,parameter::C_n=3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    / 1 -i -i\
      !D= |  i  1 -i |
      !    \ i  i  1 /       

complex(kind=dp),dimension(9),parameter ::VAL_d=(/(1.,0.),(0.,-1.),(0.,-1.),(0.,1.),(1.,0.),(0.,-1.),(0.,1.),(0.,1.),(1.,0.)/)
integer ,dimension(9)  ::INDX_d=(/1,1,1,2,2,2,3,3,3/)
integer ,dimension(9)  ::JNDX_d=(/1,2,3,1,2,3,1,2,3/)
integer,parameter::d_m=3
integer,parameter::d_n=3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    / 1 -i -i \
     !DT= |     1  -i |
      !    \       1 /       

complex(kind=dp),dimension(6),parameter ::VAL_DT=(/(1.,0.),(0.,-1.),(0.,-1.),(1.,0.),(0.,-1.),(1.,0.)/)
integer ,dimension(6)  ::INDX_DT=(/1,1,1,2,2,3/)
integer ,dimension(6)  ::JNDX_DT=(/1,2,3,2,3,3/)
integer,parameter::DT_m=3
integer,parameter::DT_N=3
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    /11  0 13 14  0\    
      !    | 0  0 23 24  0|    
      ! I= |31 32 33 34  0|
      !    | 0 42  0 44  0|     
      !    \51 52  0  0 55/     
      !
   integer, dimension(14)  ::&
 I_val=(/ 11,51,31,32,34,52,13,23,33,14,24,42,55,44 /)
      integer, dimension(14)  ::&
 I_indx= (/ 1,5,3,3,3,5,1,2,3,1,2,4,5,4 /)
      integer, dimension(14)  ::&
 I_jndx = (/ 1,1,1,2,4,2,3,3,3,4,4,2,5,4 /)
     integer ::I_m=5
     integer ::I_n=5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    / 1  1  1  1  1\   
      !    |    1  1  1  1|   
      ! J= |       1  1  1|   
      !    |          1   |   
      !    \             1/  
      !

integer,dimension(14) ::J_VAL=1
integer,dimension(14) ::J_indx=(/1,1,2,1,2,3,1,2,3,4,1,2,3,5/)
integer,dimension(14) ::J_jndx=(/1,2,2,3,3,3,4,4,4,4,5,5,5,5/)
integer,parameter::J_m=5
integer,parameter::J_n=5
integer,parameter::J_nz=14

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    /  1  2 |  4  7 \  
      !    |  0  3 |  5  8 |    
      ! M= |-------+-------|  
      !    |       |  6  9 |   
      !    \       |  0 10 /    

real(kind=dp),dimension(2,2), parameter::&
 M11=reshape((/1,0,2,3/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter::&
 M12=reshape((/4,5,7,8/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
M22=reshape((/6,0,9,10/),shape=(/2,2/))


integer,parameter::M_Mb=2
integer,parameter::M_Nb=2
integer,parameter::M_k=2
integer,parameter::M_l=2
      !
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!    /  1  0 |      \
!    |  2  5 |      |
! N= |-------+------|
!    |  3  6 | 8  0 |
!    \  4  7 | 9 10 /

real(kind=dp),dimension(2,2), parameter::&
 N11=reshape((/1,2,0,5/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter::&
 N21=reshape((/3,4,6,7/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
N22=reshape((/8,9,0,10/),shape=(/2,2/))


integer,parameter::N_Mb=2
integer,parameter::N_Nb=2
integer,parameter::N_k=2
integer,parameter::N_l=2
      !
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    / 1  1  1  1  1\   
      !    |    1  1  1  1|   
      ! T= |       1  1  1|   
      !    |          1   |   
      !    \             1/  
      !

real(kind=dp),dimension(14) ::T_VAL=1.

integer,dimension(14) ::T_indx=(/1,1,2,1,2,3,1,2,3,4,1,2,3,5/)
integer,dimension(14) ::T_jndx=(/1,2,2,3,3,3,4,4,4,4,5,5,5,5/)
integer,parameter::T_m=5
integer,parameter::T_n=5
integer,parameter::T_nz=14

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!    / 1            \
!    | 1  1         |
!  U=| 1  1  1      |
!    | 1  1  1  1   |
!    \ 1  1  1     1/

real(kind=dp),dimension(14) ::U_VAL=1.

integer,dimension(14) ::U_jndx=(/1,1,2,1,2,3,1,2,3,4,1,2,3,5/)
integer,dimension(14) ::U_indx=(/1,2,2,3,3,3,4,4,4,4,5,5,5,5/)
integer,parameter::U_m=5
integer,parameter::U_n=5
integer,parameter::U_nz=14

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 

      !    /11 12 |  0  0 | 15 16|  0  0\
      !    |21 22 |  0  0 | 25 26|  0  0| 
      !    |------+-------+------+------|
      ! X= | 0  0 | 33  0 | 35 36|  0  0|
      !    | 0  0 | 43 44 | 45 46|  0  0|
      !    |------+-------+------+------|
      !    |51 52 |  0  0 |  0  0|  0  0|
      !    \61 62 |  0  0 |  0  0|  0  0/  

real(kind=dp),dimension(2,2), parameter::&
 X11=reshape((/11,21,12,22/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter::&
 X31=reshape((/51,61,52,62/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
X22=reshape((/33,43,0,44/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
X13=reshape((/15,25,16,26/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
X23=reshape((/35,45,36,46/),shape=(/2,2/))
integer,parameter::X_Mb=3
integer,parameter::X_Nb=4
integer,parameter::X_k=2
integer,parameter::X_l=2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


      !    / 1  3 |  5  7 |  0  0\
      !    | 3  4 |  6  8 |  0  0| 
      !    |------+-------+------|
      ! Y= | 9 11 | 13 15 |  0  0|
      !    |10 12 | 15 16 |  0  0|
      !    |------+-------+------|
      !    | 0  0 |  0  0 |  0  0|
      !    \ 0  0 |  0  0 |  0  0/  
      
real(kind=dp),dimension(2,2), parameter::&
 Y11=reshape((/1,3,3,4/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter::&
 Y21=reshape((/9,10,11,12/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
Y12=reshape((/5,6,7,8/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
Y22=reshape((/13,15,15,16/),shape=(/2,2/))

integer,parameter::Y_Mb=3
integer,parameter::Y_Nb=3
integer,parameter::Y_k=2
integer,parameter::Y_l=2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !       / 1  3 |  5  7 |  0  0\
      !       | 3  4 |  6  8 |  0  0| 
      !       |------+-------+------|
      ! Y_SU= | 5  6 | 13 15 |  0  0|
      !       | 7  8 | 15 16 |  0  0|
      !       |------+-------+------|
      !       | 0  0 |  0  0 |  0  0|
      !       \ 0  0 |  0  0 |  0  0/  
      
real(kind=dp),dimension(2,2), parameter::&
 Y_SU11=reshape((/1,3,3,4/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter::&
 Y_SU21=reshape((/5,7,6,8/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
Y_SU12=reshape((/5,6,7,8/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
Y_SU22=reshape((/13,15,15,16/),shape=(/2,2/))

integer,parameter::Y_SU_Mb=3
integer,parameter::Y_SU_Nb=3
integer,parameter::Y_SU_k=2
integer,parameter::Y_SU_l=2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !    / 1  3 |  9 10 |  0  0\
      !    | 3  4 | 11 12 |  0  0| 
      !    |------+-------+------|
    !Y_SL= | 9 11 | 13 15 |  0  0|
      !    |10 12 | 15 16 |  0  0|
      !    |------+-------+------|
      !    | 0  0 |  0  0 |  0  0|
      !    \ 0  0 |  0  0 |  0  0/  
      
real(kind=dp),dimension(2,2), parameter::&
 Y_SL11=reshape((/1,3,3,4/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter::&
 Y_SL21=reshape((/9,10,11,12/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
Y_SL12=reshape((/9,11,10,12/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
Y_SL22=reshape((/13,15,15,16/),shape=(/2,2/))

integer,parameter::Y_SL_Mb=3
integer,parameter::Y_SL_Nb=3
integer,parameter::Y_SL_k=2
integer,parameter::Y_SL_l=2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !  test for Variable blocks
      !
      !    / 4  2 |  0  0  0 |  1 |  0  0  0 | -1  1 \
      !    | 1  5 |  0  0  0 |  2 |  0  0  0 |  0 -1 | 
      !    |------+----------+----+----------+-------|
      !    | 0  0 |  6  1  2 |  2 |  0  0  0 |  0  0 |
      !    | 1  5 |  2  7  1 |  0 |  0  0  0 |  0  0 | 
      !    | 1  5 | -1  2  9 |  3 |  0  0  0 |  0  0 | 
      !    |------+----------+----+----------+-------|
      ! Z= | 2  1 |  3  4  5 | 10 |  4  3  2 |  0  0 |
      !    |------+----------+----+----------+-------|
      !    | 0  0 |  0  0  0 |  4 | 13  4  2 |  0  0 |
      !    | 0  0 |  0  0  0 |  3 |  3 11  3 |  0  0 |
      !    | 0  0 |  0  0  0 |  0 |  2  0  7 |  0  0 |
      !    |------+----------+----+----------+-------|
      !    | 8  4 |  0  0  0 |  0 |  0  0  0 | 25  3 |
      !    \-2  3 |  0  0  0 |  0 |  0  0  0 |  8 12 / 
      !
      !-----------------------------------------------------------------
real(kind=dp),dimension(2,2), parameter::&
 Z11=reshape((/4.0,1.0,2.0,5.0/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
 Z15=reshape((/-1.0,0.0,1.0,-1.0/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
 Z55=reshape((/25.0,8.0,3.0,12.0/),shape=(/2,2/))
real(kind=dp),dimension(2,2), parameter:: &
 Z51=reshape((/8.0,-2.0,4.0,3.0/),shape=(/2,2/))
real(kind=dp),dimension(3,3), parameter:: & 
Z22=reshape((/6.0,2.0,-1.0,1.0,7.0,2.0,2.0,1.0,9.0/),shape=(/3,3/))
real(kind=dp),dimension(3,3), parameter:: & 
Z44=reshape((/13.0,3.0,2.0,4.0,11.0,0.0,2.0,3.0,7.0/),shape=(/3,3/))
real(kind=dp),dimension(2,1), parameter:: & 
Z13=reshape((/1.0,2.0/),shape=(/2,1/))
real(kind=dp),dimension(1,3), parameter ::& 
Z32=reshape((/3.0,4.0,5.0/),shape=(/1,3/))
real(kind=dp),dimension(1,3), parameter ::& 
Z34=reshape((/4.0,3.0,2.0/),shape=(/1,3/))
real(kind=dp),dimension(3,1), parameter ::& 
Z43=reshape((/4.0,3.0,0.0/),shape=(/3,1/))
real(kind=dp),dimension(3,1), parameter ::& 
Z23=reshape((/2.0,0.0,3.0/),shape=(/3,1/))
real(kind=dp),dimension(1,2), parameter ::& 
Z31=reshape((/2.0,1.0/),shape=(/1,2/))
real(kind=dp),dimension(1,1), parameter ::&
 Z33=reshape((/10.0/),shape=(/1,1/))

integer ,dimension(5)  ::Z_kk=(/2,3,1,3,2/)
integer ,dimension(5)  ::Z_ll=(/2,3,1,3,2/)
integer,parameter::Z_Mb=5
integer,parameter::Z_Nb=5




     !---------------------------------------------------!
      ! Test - Results for multiplication with x(i)=dbl(i)! 
      !---------------------------------------------------!

      real(kind=dp) :: res_usdot = 62.
      real(kind=dp) :: res_usaxpy(6) =(/52.,156.,162.,168.,60.,62./)
      real(kind=dp) :: res_usga(3) =(/34.,36.,38./)
      real(kind=dp) :: res_usgz_x(3)=(/14.,16.,18./)
      real(kind=dp) :: res_usgz_y(6)=(/12.,0.,0.,0.,20.,22./)
      real(kind=dp) :: res_ussc(6)=(/42.,41.,43.,45.,50.,52./)
      complex(kind=dp)  :: resc_usdot=(65.,-9.)
      complex(kind=dp) :: resc_usaxpy(6) =(/(52.,1.),(155.,54.),(161.,56.),(167.,58.),(60.,1.),(62.,1.)/)
      complex(kind=dp) :: resc_usga(3) =(/(34.,1.),(36.,1.),(38.,1.)/)
      complex(kind=dp) :: resc_usgz_x(3)=(/(14.,1.),(16.,1.),(18.,1.)/)
      complex(kind=dp) :: resc_usgz_y(6)=(/(12.,1.),(0.,0.),(0.,0.),(0.,0.),(20.,1.),(22.,1.)/)
      complex(kind=dp) :: resc_ussc(6)=(/(42.,1.),(41.,1.),(43.,1.),(45.,1.),(50.,1.),(52.,1.)/)
      real(kind=dp) :: Ax(5) =(/106.,165.,330.,260.,430./)
      real(kind=dp) :: A_SUx(5)= (/106.,165.,294.,340.,275./)
      real(kind=dp) :: A_SLx(5) = (/359.,524.,194.,260.,430./)
      real(kind=dp) :: ATx(5) = (/359.,524.,158.,340.,275./)
      real(kind=dp) :: AT_SLx(5) = (/106.,165.,294.,340.,275./)
      real(kind=dp) :: AT_SUx(5) = (/359.,524.,194.,260.,430./)
      real(kind=dp) :: Bx(3) = (/ 135.,150.,165. /)
      real(kind=dp) :: BTx(5) = (/14.,32.,50.,68.,86./)
      real(kind=dp) :: Mx(4) = (/ 45.,53.,54.,40./)
      real(kind=dp) :: MTx(4) = (/1.,8.,32.,90./)
      real(kind=dp) :: Nx(4) = (/ 1.,12.,39.,85./)
      real(kind=dp) :: NTx(4) = (/30.,56.,60.,40./)
      real(kind=dp) :: Tx(5) =(/15.,14.,12.,4.,5./)
      real(kind=dp) :: TTx(5) = (/ 1.,3.,6.,10.,11. /)
      real(kind=dp) :: Ux(5) = (/ 1.,3.,6.,10.,11. /)
      real(kind=dp) :: UTx(5) = (/15.,14.,12.,4.,5./)
      real(kind=dp) :: Xx(6) =  (/206.,346.,490.,806.,155.,185./)
      real(kind=dp) :: XTx(8) = (/674.,688.,271.,176.,350.,360.,0.,0./)
      real(kind=dp) :: Yx(6) =  (/50.,61.,130.,143.,0.,0./)
      real(kind=dp) :: YTx(6) = (/74.,92.,116.,132.,0.,0./)
      real(kind=dp) :: Y_SUx(6) = (/50.,61.,116.,132.,0.,0./)
      real(kind=dp) :: Y_SLx(6) = (/74.,92.,130.,143.,0.,0./)
      real(kind=dp) :: Zx(11) = (/15.,12.,44.,39.,68.,184.,165.,154.,77.,299.,216./)
      real(kind=dp) :: ZTx(11) = (/76.,91.,39.,65.,85.,138.,157.,134.,113.,337.,161./)
      complex(kind=dp) :: Cx(3) =(/(0.,2.),(0.,0.),(2.,0.)/)
      complex(kind=dp) :: CTx(3)=(/(0.,2.),(0.,0.),(2.,0.)/) 
      complex(kind=dp) :: CHx(3)=(/(0.,2.),(0.,0.),(2.,0.)/) 
      complex(kind=dp) :: dx(3) =(/(3.,-1.),(1.,1.),(-1.,3.)/)
      complex(kind=dp) :: dTx(3)=(/(0.,2.),(0.,0.),(2.,0.)/) 
      complex(kind=dp) :: dHx(3)=(/(0.,2.),(0.,0.),(2.,0.)/) 

      complex(kind=dp) :: DTTx(3) =(/(3.,-1.),(2.,0.),(1.,1.)/)
      complex(kind=dp) :: DTHx(3) =(/(3.,-1.),(1.,1.),(-1.,3.)/)

      integer :: Ix(5) =(/106,165,330,260,430/)
      integer :: Jx(5) =(/15,14,12,4,5/)

end module mod_test_parameters
