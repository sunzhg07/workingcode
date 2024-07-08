
program main
    implicit none
    integer:: n,k
    integer,allocatable:: a(:)
    logical :: more
    n=15
    k=10
    allocate(a(k))
    a(1)=n
    write(*,*)a
    more=.true.
    do while(more)
    call comb_next(n,k,a,more)
    enddo
    contains
    subroutine comb_next(n,k,a,more)
        implicit none
        integer,intent(in):: n,k
        integer,dimension(:):: a
        integer:: lnz,nnz,i
        logical:: more
        more=.false.
        do i=1,k
         if(a(i)==0)cycle
         lnz=i
         exit
        enddo

        if(a(lnz)==1)then
            a(lnz+1)=a(lnz+1)+1
            a(lnz)=0
        else 
            a(lnz+1)=a(lnz+1)+1
            nnz=a(lnz)
            if(lnz==1)then
                a(lnz)=nnz-1
            else
                a(lnz)=0
                a(1)=nnz-1
            endif
        endif
        if(a(k)/= n) more=.true.
        end subroutine comb_next
    end program main


