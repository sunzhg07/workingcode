program main
  integer :: a,b
  a=13
  write(*,'(b6.6)')a
  b=ibits(a,1,3)
  write(*,'(b6.6)')b
  end program main
