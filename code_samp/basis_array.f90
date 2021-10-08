program main
    implicit none
    integer,parameter :: k =6
    integer*8:: n,i,il,ct
    integer:: a(k),ah(k),al(k)
    logical:: more
    ct=1
    n=40
    do i=1,k
    a(i)=i
    al(i)=i
    ah(i)=n-k+i
    enddo

    write(*,*)a
    more=.true.
    do while(more)
    ct=ct+1
    do i=1,k-1
    if(a(i)+1 < a(i+1))exit
    enddo
    il=i
    a(il)=a(il)+1
    if(il>1) a(1:il-1)=al(1:il-1)
    more=.false.
    do i=1,k
    if(a(i)/=ah(i))more=.true.
    enddo
    write(*,*)a
    enddo


    write(*,*)ct
    ct=1
    do i=1,k
    ct=ct*(n-i+1)
    enddo
    write(*,*)ct
    ct=ct/int(gamma(real(k+1)))
    write(*,*)ct


end program
