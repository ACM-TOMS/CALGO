!
!    version 1.0 
!    built Fri Apr 16 18:03:05 CDT 2004
!
!  by B.R. Fabijonas
!     Department of Mathematics
!     Southern Methodist University
!     bfabi@smu.edu
!
!  see the file airy_README for an explanation
!***
!************************************************************************
!***
      module airy_functions_real_single
        implicit none
        private
        public :: airy_ai, airy_bi, airy_ai_zero, airy_bi_zero
        public :: airy_info, airy_aux, airy_aux_info
        integer, parameter :: prd = kind(0.0e0)
        include 'airy_head'
!*
!  interfaces 
!*
        interface airy_ai
          module procedure airy_air 
        end interface 
        interface airy_bi
          module procedure airy_bir 
        end interface 
        interface airy_ai_zero
  	  module procedure ai_zeror 
          module procedure ai_zerorv 
        end interface
        interface airy_bi_zero
  	  module procedure bi_zeror 
          module procedure bi_zerorv 
        end interface
        interface airy_info 
          module procedure airy_paramsr 
        end interface
        interface airy_aux                                          
  	  module procedure airy_auxr 
        end interface
        interface airy_aux_info 
          module procedure airy_paramsr_aux 
        end interface
        contains
         include 'airy_real'
         include 'airy_parameters'
      end module airy_functions_real_single 
!***
!************************************************************************
!***
      module airy_functions_complex_single
        implicit none
        private
        public :: airy_ai, airy_ai_zero, airy_bi_zero
        integer, parameter :: prd = kind(0.0e0)
        include 'airy_head'
!*
!  interfaces 
!*
        interface airy_ai
          module procedure airy_aic
          module procedure airy_ai_rayc 
        end interface 
        interface airy_ai_zero
   	  module procedure ai_zeroc 
          module procedure ai_zerocv 
        end interface
        interface airy_bi_zero
   	  module procedure bi_zeroc         
          module procedure bi_zerocv 
        end interface
        contains
         include 'airy_complex'
         include 'airy_parameters'
      end module airy_functions_complex_single 
!***
!************************************************************************
!***
      module airy_functions_real_double
        implicit none
        private
        public :: airy_ai, airy_bi, airy_ai_zero, airy_bi_zero
        public :: airy_info, airy_aux, airy_aux_info
        integer, parameter :: prd = kind(0.0d0)
        include 'airy_head'
!*
!  interfaces 
!*
        interface airy_ai
          module procedure airy_air 
        end interface 
        interface airy_bi
          module procedure airy_bir 
        end interface 
        interface airy_ai_zero
  	  module procedure ai_zeror 
          module procedure ai_zerorv 
        end interface
        interface airy_bi_zero
  	  module procedure bi_zeror 
          module procedure bi_zerorv 
        end interface
        interface airy_info 
          module procedure airy_paramsr 
        end interface
        interface airy_aux                                          
  	  module procedure airy_auxr 
        end interface
        interface airy_aux_info 
          module procedure airy_paramsr_aux 
        end interface
        contains
         include 'airy_real'
         include 'airy_parameters'
      end module airy_functions_real_double 
!***
!************************************************************************
!***
      module airy_functions_complex_double
        implicit none
        private
        public :: airy_ai, airy_ai_zero, airy_bi_zero
        integer, parameter :: prd = kind(0.0d0)
        include 'airy_head'
!*
!  interfaces 
!*
        interface airy_ai
          module procedure airy_aic
          module procedure airy_ai_rayc 
        end interface 
        interface airy_ai_zero
   	  module procedure ai_zeroc 
          module procedure ai_zerocv 
        end interface
        interface airy_bi_zero
   	  module procedure bi_zeroc         
          module procedure bi_zerocv 
        end interface
        contains
         include 'airy_complex'
         include 'airy_parameters'
      end module airy_functions_complex_double 
!***
!************************************************************************
!***
!      module airy_functions_real_quad
!        implicit none
!        private
!        public :: airy_ai, airy_bi, airy_ai_zero, airy_bi_zero
!        public :: airy_info, airy_aux, airy_aux_info
!        integer, parameter :: prd = kind(0.0q0)
!        include 'airy_head'
!!*
!!  interfaces 
!!*
!        interface airy_ai
!          module procedure airy_air 
!        end interface 
!        interface airy_bi
!          module procedure airy_bir 
!        end interface 
!        interface airy_ai_zero
!  	  module procedure ai_zeror 
!          module procedure ai_zerorv 
!        end interface
!        interface airy_bi_zero
!  	  module procedure bi_zeror 
!          module procedure bi_zerorv 
!        end interface
!        interface airy_info 
!          module procedure airy_paramsr 
!        end interface
!        interface airy_aux                                          
!  	  module procedure airy_auxr 
!        end interface
!        interface airy_aux_info 
!          module procedure airy_paramsr_aux 
!        end interface
!        contains
!         include 'airy_real'
!         include 'airy_parameters'
!      end module airy_functions_real_quad 
!!***
!!************************************************************************
!!***
!      module airy_functions_complex_quad
!        implicit none
!        private
!        public :: airy_ai, airy_ai_zero, airy_bi_zero
!        integer, parameter :: prd = kind(0.0q0)
!        include 'airy_head'
!!*
!!  interfaces 
!!*
!        interface airy_ai
!          module procedure airy_aic
!          module procedure airy_ai_rayc 
!        end interface 
!        interface airy_ai_zero
!   	  module procedure ai_zeroc 
!          module procedure ai_zerocv 
!        end interface
!        interface airy_bi_zero
!   	  module procedure bi_zeroc         
!          module procedure bi_zerocv 
!        end interface
!        contains
!         include 'airy_complex'
!         include 'airy_parameters'
!      end module airy_functions_complex_quad 
!!***
!!************************************************************************
!!***
