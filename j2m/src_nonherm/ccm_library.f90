!--------------------------------------------------------------------                                                                                                                
FUNCTION calc_objective(work,mean,work_add,istart,istop,where_add)
  IMPLICIT NONE
  REAL :: calc_objective
  REAL :: work(1:istop)
  REAL :: mean, work_add
  INTEGER :: istart,istop,where_add
  INTEGER :: i

  calc_objective = 0.0
  DO i=istart,istop
     IF(i.NE.where_add)THEN
       calc_objective=calc_objective + ABS(work(i)-mean)
     ELSE
       calc_objective=calc_objective+ ABS(work(i) + work_add - mean)
     ENDIF
  ENDDO
end FUNCTION calc_objective

! sort routine 

subroutine sort_channel( c, b, n )
  implicit none 
  integer, intent(in) :: n
  integer, intent(inout) :: c(n), b(n)
  integer :: i,j, tmp1, tmp2


  do i = 1, n 
     do j = i+1, n 
        if ( c(j) < c(i) ) then
           tmp1 = c(j)
           c(j) = c(i)
           c(i) = tmp1
           
           tmp2 = b(j) 
           b(j) = b(i)
           b(i) = tmp2
           

        end if
     end do
  end do
  
end subroutine sort_channel


!     swaps values of 2 integers
SUBROUTINE SWAP_AB(a, b)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: a, b
  INTEGER :: c
  c = a; a = b; b = c
END SUBROUTINE SWAP_AB

!
!     Function to check triangular relations     
!
LOGICAL FUNCTION triag(i,j,k)
  IMPLICIT NONE
  INTEGER :: i, j, k
  triag = ((i-j-k)*(i-ABS(j-k)) > 0)

END FUNCTION triag
!
!      Function to calculate norm of g-mat    
!
real*8 FUNCTION dij(ja,jb)
  IMPLICIT NONE
  INTEGER :: ja, jb
  IF(ja == jb ) THEN
     dij=sqrt(2.d0)
  ELSE
     dij=1.0
  ENDIF
  
END FUNCTION dij

!
!     Function to calculate phase factors (-1)**(n)
!
INTEGER FUNCTION iph(n)
  IMPLICIT NONE
  INTEGER :: n
  iph=(-1)**n
END FUNCTION iph


FUNCTION delta(a,b)
  IMPLICIT NONE
  INTEGER :: a,b
  REAL*8  :: delta
  delta=0.0
  if(a == b)delta=1.0
end FUNCTION delta

SUBROUTINE zero(n,a)
  IMPLICIT NONE
  INTEGER :: n,i
  REAL*8  :: a(n)
  do i=1,n
    a(i)=0.0
  end do
end subroutine zero


SUBROUTINE lapack_diag(h, cvec, cvecl,ceig, n )


  implicit none
  integer, intent(in) :: n
  complex*16, dimension(n,n), intent(in) :: h
  COMPLEX*16, DIMENSION(n,n), intent(out) :: cvec, cvecl
  COMPLEX*16, DIMENSION(n,n) ::  vl,vr
  COMPLEX*16, DIMENSION(n), intent(out) :: ceig
  DOUBLE PRECISION, DIMENSION(2*n) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  INTEGER :: lda, ldvl, ldvr, info, lwork
  CHARACTER*1 :: jobvl, jobvr
  complex*16 :: norm
  integer :: i,j

  jobvl = 'V' ;  jobvr = 'V';  lda = n
  ldvl = n;  ldvr = n;  lwork = 10000
  ceig = 0.; cvec = 0.; cvecl = 0.
!  write(*,*)h
  CALL zgeev( jobvl, jobvr, n, h, lda, ceig, cvecl, ldvl, cvec, ldvr, &
    work1, lwork, rwork, info )

  ! berggren normalization
  do i = 1, n
    norm = sum( cvec(:,i)*cvec(:,i) )
    cvec(:,i) = cvec(:,i)/sqrt(norm)
    norm = sum( cvecl(:,i)*cvec(:,i) )
    !write(*,*)norm
    if(real(norm) < 0.d0 )then
    cvecl(:,i)=-cvecl(:,i)/sqrt(-norm)
    else
    cvecl(:,i)=cvecl(:,i)/sqrt(norm)
    endif
  
     
  end do

do i=1,n
   do j=1,n
    norm=sum(cvecl(:,i)*cvecl(:,j))
   ! write(6,*)'ss',i,j, norm
enddo
enddo

!  call hf_sort_cmplx2(ceig,cvec,cvecl, n)
  call hf_sort_overlap(ceig,cvec,cvecl,n)


end SUBROUTINE lapack_diag


SUBROUTINE hf_sort_overlap(a,b,c, n)
  IMPLICIT NONE

  INTEGER :: i, j, n
  complex*16, DIMENSION(n), INTENT(INOUT) :: a
  complex*16, DIMENSION(n,n), INTENT(INOUT) :: b,c
  complex*16 :: temp1, temp2
  complex*16, DIMENSION(n) :: temp3

  do i=1,n
    temp2=b(i,i)
    do j=i,n
      if(abs(temp2) < abs(b(i,j)))then

           temp1 = a(i)
           a(i) = a(j)
           a(j) = temp1

           temp3(:) = b(:,i)
           b(:,i) = b(:,j)
           b(:,j) = temp3(:)


           temp3(:) = c(:,i)
           c(:,i) = c(:,j)
           c(:,j) = temp3(:)
      endif
    enddo

  enddo


END SUBROUTINE hf_sort_overlap



SUBROUTINE hf_sort_cmplx2(a,b,c, n)
  IMPLICIT NONE

  INTEGER :: i, j, n
  complex*16, DIMENSION(n), INTENT(INOUT) :: a
  complex*16, DIMENSION(n,n), INTENT(INOUT) :: b,c
  complex*16 :: temp1, temp2
  complex*16, DIMENSION(n) :: temp3

  DO i = 1, n
     DO j = 1, n
        IF ( real( a(i) )  < real( a(j) ) ) THEN

           temp1 = a(i)
           a(i) = a(j)
           a(j) = temp1

           temp3(:) = b(:,i)
           b(:,i) = b(:,j)
           b(:,j) = temp3(:)


           temp3(:) = c(:,i)
           c(:,i) = c(:,j)
           c(:,j) = temp3(:)

        END IF
     END DO
  END DO

END SUBROUTINE hf_sort_cmplx2



SUBROUTINE sqrtmat( aa, a, a_inv ,n ) 
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: aa
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: a, a_inv
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: x, x0, x1, x2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: i_mat, temp
  INTEGER :: i,j,k
  COMPLEX*16 :: d

  ALLOCATE( x(n+n,n+n), x0(n+n,n+n), x1(n+n,n+n),x2(n+n,n+n))
  ALLOCATE( i_mat(n,n),temp(n,n))
  ! setup real identity matrix only
  i_mat = (0.D0,0.D0)
  DO i = 1, n
     i_mat(i,i) = (1.d0, 0.d0)
  ENDDO
  DO i = 1, 2*n
     DO j = 1, 2*n
        x0(j,i) = (0.d0, 0.d0)
        x1(j,i) = (0.d0, 0.d0)
        x2(j,i) = (0.d0, 0.d0)
        x(j,i) = (0.d0, 0.d0)
     ENDDO
  ENDDO
  DO i = 1, n
     DO j = 1, n
        temp(j,i) = (0.d0, 0.d0)
     ENDDO
  ENDDO
  DO i = n+1, 2*n
     DO j = 1, n
        x0(j,i) = aa(j,i-n)
     ENDDO
  ENDDO
  DO i = 1, n
     DO j = n+1, 2*n
        x0(j,i) = i_mat(j-n,i)
     ENDDO
  ENDDO
  k = 0
  DO WHILE( MAXVAL(ABS(temp-aa)) > 1.D-14 .AND.  k < 1000 )
     x1 = x0 
     x2 = x0 
     CALL cmplxmatinv(x2,n+n,d)
     x = 0.5d0 * ( x1 + x2 ) 
     x0 = x 
     k = k + 1
     DO i = 1, n
        DO j =1, n
           a(i,j) = x(i,j+n)
           a_inv(i,j) = x(i+n,j) 
        ENDDO
     ENDDO
     temp = MATMUL( a,a )
  ENDDO
  DEALLOCATE(i_mat,temp); DEALLOCATE(x,x0,x1,x2)

END SUBROUTINE sqrtmat
!
!    F90 program library, adapted from Numerical Recipes
!    All functions have been translated to F90 from F77
!    This is the complex*16 version of the 
!    routines to do matrix inversion, from Numerical
!    Recipes, Teukolsky et al. Routines included
!    below are MATINV, LUDCMP and LUBKSB. See chap 2
!    of Numerical Recipes for further details
!
SUBROUTINE cmplxmatinv(a,n,d)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT)  :: a
  COMPLEX*16, ALLOCATABLE :: y(:,:)
  COMPLEX*16 :: d
  INTEGER, ALLOCATABLE :: indx(:)

  ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
  y=0.
  !     setup identity matrix
  DO i=1,n
     y(i,i)=(1.d0, 0.d0) 
  ENDDO
  !     LU decompose the matrix just once
  CALL  cmplxlu_decompose(a,n,indx,d)

  !     Find inverse by columns
  DO j=1,n
     CALL cmplxlu_linear_equation(a,n,indx,y(:,j))
  ENDDO
  !     The original matrix a was destroyed, now we equate it with the inverse y 
  a=y

  DEALLOCATE ( y ); DEALLOCATE ( indx )

END SUBROUTINE cmplxmatinv
!
!     Given an NxN matrix A(N,N), this routine replaces it by the LU 
!     decomposed one, where the matrix elements are stored in the same 
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
SUBROUTINE cmplxlu_decompose(a,n,indx,d)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  COMPLEX*16 :: sum, dum, tiny, d, aamax
  COMPLEX*16, DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  COMPLEX*16, ALLOCATABLE :: vv(:)

  tiny= ( 1.0D-20, 0.d0 ) 
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=(0.D0,0.D0)
     DO j=1,n
        IF (ABS(a(i,j)) > ABS(aamax) ) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (abs(aamax) < 1E-12) STOP 'Singular matrix.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=(0.D0,0.D0)
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (ABS( dum ) >= ABS( aamax) ) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(abs(a(j,j)) < 1E-13) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(abs(a(n,n)) < 1E-13 )  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE cmplxlu_decompose
!
!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion
!
SUBROUTINE cmplxlu_linear_equation(a,n,indx,b)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  COMPLEX*16 :: sum 
  COMPLEX*16, DIMENSION(n,n) :: a
  COMPLEX*16, DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE cmplxlu_linear_equation

