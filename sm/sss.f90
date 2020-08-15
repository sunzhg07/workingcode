program sm
  implicit none
  integer*16:: u,ur,nh,no
  integer*16:: i
  u=1
!  u=3
! do while((u>=0) .and.( u<=96))
!   ur=u .and. -u
!   nh=u+ur
!   no=u.neqv.nh
!   no=no/ur
!   no=ishft(no,-2)
!   u=nh .or. no
!   call decode(u)
!enddo



!  do i=1,128
!    u=ishft(u,-1)
!  write(*,*)u
!enddo
do while(u==1)
  read(*,*) i
  call decode(i)
enddo
end program sm

subroutine decode(a)
  implicit none
  integer*16,intent(in):: a
  integer*16:: b,c(35),i,d,e
  e=a
  b=1
  do i=1,35
    d= (e.and.b)
    if(d==1)then
    c(i)=1
   else
    c(i)=0
  endif
  e=ishft(e,-1)
  enddo
  write(*,'(32I1)')c
  end subroutine
