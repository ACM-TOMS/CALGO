subroutine set_num_threads(n)
  integer, intent(in) :: n
  call mkl_set_num_threads(n)
end subroutine set_num_threads
