
program main
    implicit none
    integer:: k,i
    integer,allocatable:: a_min(:),a_max(:),t(:)
    logical :: more
    k=3
    allocate(a_min(k))
    allocate(a_max(k))
    allocate(t(k))
    a_min=(/2,1,2/)
    t=a_min
    a_max=(/4,3,4/)
    more=.true.
    write(*,*)t
    do while(more)
      do i=1,k
      t(i)=t(i)+1
      if(t(i)<=a_max(i))exit
      enddo
      if(i>1)t(1:i-1)=a_min(1:i-1)

      more=.false.
      do i=1,k
      if(t(i) < a_max(i))more=.true.
      enddo
      write(*,*)t
    enddo
    end program main


