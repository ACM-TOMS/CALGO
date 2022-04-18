subroutine set_num_threads(n)
  integer, intent(in) :: n
  call openblas_set_num_threads(n)
end subroutine set_num_threads
