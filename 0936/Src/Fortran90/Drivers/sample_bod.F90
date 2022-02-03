! Time-stamp:  <2015-12-10 11:43:17 m>
!>> 2015-12-06 sample Krogh Set up for use as include for generic use
type, public :: sample_ty
 type(messy_ty) :: e
 integer ::  what ! Takes values from the first parameter statement below.
end type sample_ty

integer, parameter :: setup_sample=0, partial_message=1, finish_message=2
real(fp), parameter :: pi=3.1415926535897932384626433832795028841971693993_fp

contains
  subroutine sample(s)
    type(sample_ty) :: s
    select case (s%what)
    case(setup_sample) ! The setup call
      s%e%ename="sample"
      return
    case(partial_message) ! Partial print of error1
      call messy(s%e, "$E34Start of error message from sample what=$I.$N&
        &I want some $R!$C",&
        & idat=[ 1, s%what ], rdat=[pi])
    case(finish_message) ! Finish this error message
      call messy(s%e, "Finish with what=$I.", [ s%what ])
    case default
      call messy(s%e, "$E56Fatal error called sample with what = $I.",&
      & [99, s%what ])
    end select
    return
  end subroutine sample
