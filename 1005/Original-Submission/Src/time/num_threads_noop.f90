subroutine set_num_threads(n)
  integer, intent(in) :: n
  if (n == 0) continue
  !no-op (refblas uses always just one thread)
end subroutine set_num_threads
