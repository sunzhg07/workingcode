program m_sm
  use m_const
  use m_sp
  use m_int
  use combine
  use mbasis
  implicit none
  call ini_files
  call basis
  call ini_combin
  call matrix
  write(*,*)'total dimension:', num_mbasis
  call lapack_diag(hm, vec, vecl,em, num_mbasis)
  call charge_radii(om, vec, vecl,oem, num_mbasis)

end program m_sm

subroutine ini_files
  use m_int
  use m_sp
  use m_const
  use m_file
  implicit none
  integer :: i,n
  open(file='sm.dat',unit=5)
  read(5,*);read(5,*)np,nn
  read(5,*);
  read(5,*)num_interaction
  read(5,*)
  read(5,*)jz
  read(5,*)
  read(5,*)onebd_int
  read(5,*)
  read(5,*)twobd_int
  read(5,*)
  read(5,*)threebd_int
  read(5,*)
  read(5,*)onebd_opr
  read(5,*)
  read(5,*)twobd_opr
  read(5,*)
  read(5,*)threebd_opr
  read(5,*)
  read(5,*)output
  read(5,*)
  read(5,*)sp_file

  open(unit=7,file=sp_file)
  read(7,*)tot_orbs
  call allocate_sp_array(all_orbits,tot_orbs)

  np_orb=0
  nn_orb=0
  all_orbits%total_orbits=tot_orbs
  do i=1,tot_orbs
    read(7,*)all_orbits%nn(i),&
      & all_orbits%nshell(i),&
      & all_orbits%ll(i),&
      & all_orbits%jj(i),&
      & all_orbits%mm(i),&
      & all_orbits%itzp(i)
    if(all_orbits%itzp(i)==-1)then
      np_orb=np_orb+1
    elseif(all_orbits%itzp(i)==1)then
      nn_orb=nn_orb+1
    else
      write(*,*)'error in single particle itzp'
    endif

  enddo
  close(7)
  allocate(obs(tot_orbs))
  obs=0
  n=0
  do i=1,all_orbits%total_orbits
    if(all_orbits%itzp(i)==1)then
      n=n+1
      obs(n)=i
    endif
  enddo
  do i=1,all_orbits%total_orbits
    if(all_orbits%itzp(i)==-1)then
      n=n+1
      obs(n)=i
    endif
  enddo
  write(*,*)obs


  open(unit=8,file=onebd_int)
  read(8,*)n
  call ini_mat(n,1)
  call sort_mat(1)
  open(unit=9,file=twobd_int)
  read(9,*)n
  call ini_mat(n,2)
  call sort_mat(2)
  open(unit=18,file=onebd_opr)
  read(18,*)n
  call ini_mat(n,4)
  call sort_mat(4)
  open(unit=19,file=twobd_opr)
  read(19,*)n
  call ini_mat(n,5)
  call sort_mat(5)

  if(num_interaction==3)then
    write(*,*)num_interaction
    open(unit=10,file=threebd_int)
    read(10,*)n
    call ini_mat(n,3)
    call sort_mat(3)
  !call very_mat

    open(unit=20,file=threebd_opr)
    read(20,*)n
    call ini_mat(n,6)
    call sort_mat(6)

  endif
end subroutine ini_files



subroutine basis
  use mbasis
  use m_sp
  use bit_op
  use m_const
  implicit none
  integer*16 :: low, high,one,tp
  integer*16:: u,ur,nh,no
  integer :: i,j,a,b,tmp
  one=1
  if(np==0)then 
    write(*,*) 'number protons 0'
    np_mbasis=1
    allocate(mbsp(1))
    allocate(mzp(1))
    mbsp(1)=0
    mzp(1)=0
  else
    tp=0
    low=0
    do i=1,np
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,np_orb-np)
    !write(*,*)low,high
    u=low
    np_mbasis=1
    do while((u>=low) .and.( u<high))
      ur=and(u , -u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh,no)
      np_mbasis=np_mbasis+1
      !   call decode(u)
    enddo
    !write(*,*)'proton basis',np_mbasis
    !call decode(low)
    !call decode(high)

    allocate(mbsp(np_mbasis))
    allocate(mzp(np_mbasis)); mzp=0
    np_mbasis=1
    tp=0
    low=0
    do i=1,np
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,np_orb-np)
    mbsp(1)=low
    u=low
    !write(*,*)low,high

    do while((u>=low) .and.( u<high))
      ur=and(u,-u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh, no)
      np_mbasis=np_mbasis+1
      mbsp(np_mbasis)=u
    enddo
    do i=1,np_mbasis
      tp=mbsp(i)
      tmp=0
      do j=1,np_orb
        if(and(tp,one) ==one)then
          tmp=tmp+all_orbits%mm(j+nn_orb)
        endif
        tp=ishft(tp,-1)
      enddo
      mzp(i)=tmp
    enddo


  endif

  !write(*,*)'proton basis', np_mbasis

  tp=0
  low=0
  if(nn==0)then 
    write(*,*) 'number protons 0'
    nn_mbasis=1
    allocate(mbsn(1))
    allocate(mzn(1)); mzn=0
    mbsn(1)=0
  else
    do i=1,nn
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,nn_orb-nn)

    !write(*,*)low,high,nn
    u=low
    nn_mbasis=1
    do while((u>=low) .and.( u<high))
      ur=and(u,-u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh, no)
      nn_mbasis=nn_mbasis+1
      ! call decode(u)
    enddo
    !call decode(low)
    !call decode(high)

    allocate(mbsn(nn_mbasis))
    allocate(mzn(nn_mbasis)); mzn=0

    tp=0
    low=0
    do i=1,nn
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,nn_orb-nn)

    !write(*,*)low,high,nn
    u=low
    nn_mbasis=1
    mbsn(1)=low
    do while((u>=low) .and.( u<high))
      ur=and(u,-u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh,no)
      nn_mbasis=nn_mbasis+1
      mbsn(nn_mbasis)=u
    enddo

    do i=1,nn_mbasis
      tp=mbsn(i)
      tmp=0
      do j=1,nn_orb
        if(and(tp,one) ==one)then
          tmp=tmp+all_orbits%mm(j)
        endif
        tp=ishft(tp,-1)
      enddo
      mzn(i)=tmp
    enddo


  endif


  !write(*,*)'number ',num_mbasis


  num_mbasis=0
  do a=1,np_mbasis
    do b=1,nn_mbasis
      if(mzp(a)+mzn(b)==jz)then
        num_mbasis=num_mbasis+1
      endif
    enddo
  enddo



  allocate(mbs(num_mbasis,2));mbs=0
  allocate(mz(num_mbasis));mz=0

  num_mbasis=0
  do a=1,np_mbasis
    do b=1,nn_mbasis
      if(mzp(a)+mzn(b)==jz)then
        num_mbasis=num_mbasis+1
        mbs(num_mbasis,1)=a
        mbs(num_mbasis,2)=b
      endif
    enddo
  enddo


  allocate (hm(num_mbasis,num_mbasis));hm=0.d0
  allocate (om(num_mbasis,num_mbasis));om=0.d0
  allocate (jm(num_mbasis,num_mbasis));jm=0.d0
  allocate (vec(num_mbasis,num_mbasis));vec=0.d0
  allocate (vecl(num_mbasis,num_mbasis));vecl=0.d0
  allocate (em(num_mbasis));em=0.d0
  allocate (oem(num_mbasis));em=0.d0

  allocate (idmaskn(max(nn_orb,np_orb))); idmaskn=0
  allocate (idmaskp(max(nn_orb,np_orb))); idmaskp=0
  idmaskn(1)=1
  idmaskp(1)=1
  do a=2,max(nn_orb,np_orb)
    idmaskn(a)=idmaskn(a-1)*2+1
    idmaskp(a)=idmaskp(a-1)*2+1
  enddo
  write(*,*)idmaskn




end subroutine basis

subroutine matrix
  use bit_op
  use mbasis
  use m_sp
  use m_int
  use m_const
  implicit none
  ! initial h
  integer*16:: bra,ket
  integer*16:: bra_p,bra_n,ket_p,ket_n,c1,c2,c3,c4,ucbit,nucbit
  integer*16:: nxor,one
  real*8::tmp,tmp1,val,val1
  integer*16:: ia,ib,ic,id,ie,ik,i,j,k,l,idx,ii,ij,il
  integer*16:: bt_maskp(20),bt_maskn(20),bt_mask(40),bit_change(6),ncbitn,ncbitp,bit_uc(40)
  integer :: sig

  do bra=1,num_mbasis
    bra_p=mbs(bra,1)
    bra_p=mbsp(bra_p)
    bra_n=mbs(bra,2)
    bra_n=mbsn(bra_n)
    do ket=1,num_mbasis
      ket_p=mbs(ket,1)
      ket_p=mbsp(ket_p)
      ket_n=mbs(ket,2)
      ket_n=mbsn(ket_n)
      c1=xor(bra_p,ket_p)
      c2=xor(bra_n,ket_n)
      nxor=popcnt(c1)+popcnt(c2)
     ! call decode(bra_n)
     ! call decode(ket_n)
     ! write(*,*)bra_n,ket_n,nxor
      one=01
      !write(*,*)'nxor',bra_p,ket_p,bra_n,ket_n,nxor
      bt_mask=0
      bt_maskn=0
      bt_maskp=0

      tmp=0.d0
      tmp1=0

      select case(nxor)
      case (0)
        ! diagonal
        ! add 1b 2b 3b
        c1=bra_p
        c2=bra_n
        bt_maskn=0
        bt_maskp=0
        bt_mask=0
        ib=0;ic=0
        !    write(*,*)c1,c2
        ib=0;ic=0
        do ia=1,max(nn_orb,np_orb)
          if(and(c1, one)== one)then
            ib=ib+1
            bt_mask(ib+nn)=ia+nn_orb
          endif
          if(and(c2, one)== one)then
            ic=ic+1
            bt_mask(ic)=ia
            !       write(*,*)'hi',ic,c2
          endif
          c1=ishft(c1,-1)
          c2=ishft(c2,-1)
        enddo

        if(ib/=np .or. ic/= nn)write(*,*)'fatal error in case 0'
        ic=ib+ic


        !write(*,*)'ibic',ib,ic,bra_p,bra_n


        ! add 1b interaction
        do id=1,ic
          idx=0
          i=bt_mask(id)
          idx=ishft(i,7)
          idx=idx+i
          val=0.0
          !call decode(idx)
          call fetch_mat(idx,1,val)
          tmp=tmp+val
          call fetch_mat(idx,4,val1)
          tmp1=tmp1+val1
        enddo


        ! add 2b interaction <ij|V|ij>
        do ib=1,ic-1
          do id=ib+1,ic
            idx=0
            i=bt_mask(ib)
            j=bt_mask(id)
            idx=idx+ishft(i,21)
            idx=idx+ishft(j,14)
            idx=idx+ishft(i,7)
            idx=idx+j
            val=0.0
            ! call decode(idx)
            call fetch_mat(idx,2,val)
            tmp=tmp+val
          call fetch_mat(idx,5,val1)
          tmp1=tmp1+val1
          enddo
        enddo

        ! add 3b interaction <ijk|V|ijk>
        if(num_interaction==3)then
          do ib=1,ic-2
            do id=ib+1,ic-1
              do ie=id+1,ic

                idx=0
                i=bt_mask(ib)
                j=bt_mask(id)
                k=bt_mask(ie)
                idx=idx+ishft(i,35)
                idx=idx+ishft(j,28)
                idx=idx+ishft(k,21)
                idx=idx+ishft(i,14)
                idx=idx+ishft(j,7)
                idx=idx+k
                val=0.0
                call fetch_mat(idx,3,val)
                tmp=tmp+val
          call fetch_mat(idx,6,val1)
          tmp1=tmp1+val1
              enddo
            enddo
          enddo
        endif

      case (2)
        bt_mask=0
        ! onebody+ twobody+ 3b
        c1=xor(bra_p , ket_p)
        c2=xor(bra_n , ket_n)

        if(c1/=0)then
          if(c2/=0)stop 'case 2'
          ucbit=and(not(c1) ,bra_p)
          c2=and(c1,ket_p)
          c1=and(c1,bra_p)
          i=0; j=0
          do ia=1,np_orb
            if(and(c1,one)==one)then
              if(i/=0)write(*,*)'fatal error'
              i=ia
            endif
            if(and(c2,one)==one)then
              if(j/=0)write(*,*)'fatal error'
              j=ia
            endif
            c1=ishft(c1,-1)
            c2=ishft(c2,-1)
          enddo

          i=i+nn_orb
          j=j+nn_orb
          ! i,j is the chaged bit in bra and ket

          ib=0
          bt_mask=0
          do ia=1,np_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib+nn)=ia+nn_orb
            endif
            ucbit=ishft(ucbit,-1)
          enddo
          nucbit=ib

          ucbit=bra_n
          ib=0
          do ia=1,nn_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib)=ia
            endif
            ucbit=ishft(ucbit,-1)
          enddo
          nucbit=nucbit+ib





        else
          if(c1/=0)stop 'case 2'
          ucbit=and(not(c2), bra_n)
          c1=and(c2,bra_n)
          c2=and(c2,ket_n)
          i=0; j=0
          do ia=1,max(nn_orb,np_orb)
            if(and(c1,one)==one)then
              if(i/=0)write(*,*)'fatal error'
              i=ia
            endif
            if(and(c2,one)==one)then
              if(j/=0)write(*,*)'fatal error'
              j=ia
            endif
            c1=ishft(c1,-1)
            c2=ishft(c2,-1)
          enddo


          ib=0
          bt_mask=0
          do ia=1,nn_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib)=ia
            endif
            ucbit=ishft(ucbit,-1)
          enddo
          ucbit=bra_p
          do ia=1,np_orb
            if(and(ucbit, one)== one)then
              ib=ib+1
              bt_mask(ib)=ia+nn_orb
            endif
            ucbit=ishft(ucbit,-1)
          enddo
          nucbit=ib


        endif


        ! 1b
        idx=0
        idx=idx+ishft(i,7)
        idx=idx+j
        val=0.0
        sig=0
        sig=ch_sign(bra,int(1,16),i)
        sig=ch_sign(ket,int(1,16),j)+sig
        call fetch_mat(idx,1,val)
        tmp=tmp+val*(-1)**sig
        call fetch_mat(idx,4,val1)
        tmp1=tmp1+val1*(-1)**sig

        ! 2b


        do ia=1,nucbit
          k=bt_mask(ia)
          sig=0
          idx=0
          ii=min(i,k)
          ij=max(i,k)
          idx=idx+ishft(ii,21)
          idx=idx+ishft(ij,14)
          sig=ch_sign(bra,ii,ij)
          ii=min(j,k)
          ij=max(j,k)
          sig=sig+ch_sign(ket,ii,ij)
          idx=idx+ishft(ii,7)
          idx=idx+ij
          val=0.0

          call fetch_mat(idx,2,val)
          tmp=tmp+val*(-1.0)**sig
        call fetch_mat(idx,5,val1)
        tmp1=tmp1+val1*(-1)**sig
        enddo

        ! 3b
        if(num_interaction==3)then
          do ia=1,nucbit-1
            k=bt_mask(ia)
            do ic=ia+1,nucbit
              l=bt_mask(ic)
              if(k>l)stop 'fatal error 3b'
              sig=0
              idx=0
              if(i<k)then
                ii=i
                ij=k
                ik=l
              else if(i>l)then
                ii=k
                ij=l
                ik=i
              else
                ii=k
                ij=i
                ik=l
              endif
              sig=ch_sign(bra,int(0,16),ii)
              sig=ch_sign(bra,ij,ik)+sig
              idx=0

              idx=idx+ishft(ii,35)
              idx=idx+ishft(ij,28)
              idx=idx+ishft(ik,21)
              if(j<k)then
                ii=j
                ij=k
                ik=l
              elseif(j>l)then
                ii=k
                ij=l
                ik=j
              else
                ii=k
                ij=j
                ik=l
              endif

              idx=idx+ishft(ii,14)
              idx=idx+ishft(ij,7)
              idx=idx+ik
              val=0.0
              sig=ch_sign(ket,int(0,16),ii)+sig
              sig=ch_sign(ket,ij,ik)+sig
              call fetch_mat(idx,3,val)
              tmp=tmp+val*(-1)**sig
        call fetch_mat(idx,6,val1)
        tmp1=tmp1+val1*(-1)**sig
            enddo
          enddo
        endif








      case (4)
        ! twobdy+3b
        bit_uc=0


        ! twobody
        c1=xor(bra_p , ket_p)
        c2=xor(bra_n , ket_n)
        ! get unchanged bit
        c3= and(not(c1) , bra_p)
        c4= and(not(c2) , bra_n)
        ncbitn=0
        do ia=1,np_orb
          if(and(c3,one) == one )then
            ncbitn=ncbitn+1
            bit_uc(ncbitn)=ia+nn_orb
          endif
          c3=ishft(c3,-1)
        enddo
        do ia=1,nn_orb
          if(and(c4,one) == one)then
            ncbitn=ncbitn+1
            bit_uc(ncbitn)=ia
          endif
          c4=ishft(c4,-1)
        enddo





        ! add twobody interaction
        c1=xor(bra_p , ket_p)
        c2=xor(bra_n , ket_n)

        c3=and(c1,bra_p)
        c4=and(c2,bra_n)
        ncbitp=0
        bit_change=0
        do ia=1,np_orb
          if(and(c3,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia+nn_orb
          endif
          c3=ishft(c3,-1)
        enddo

        do ia=1,nn_orb
          if(and(c4,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia
          endif
          c4=ishft(c4,-1)
        enddo

        if(ncbitp/=2)write(*,*)'error in case 4 <bra|ket>, a', ncbitp, and(c2,bra_n)

        c3=and(c1,ket_p)
        c4=and(c2,ket_n)
        do ia=1,np_orb
          if(and(c3,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia+nn_orb
          endif
          c3=ishft(c3,-1)
        enddo
        do ia=1,nn_orb
          if(and(c4,one)==one)then
            ncbitp=ncbitp+1
            bit_change(ncbitp)=ia
          endif
          c4=ishft(c4,-1)
        enddo
        if(ncbitp/=4)write(*,*)'error in case 4 <bra|ket> b'



        sig=0


        ii=min(bit_change(1),bit_change(2))
        ij=max(bit_change(1),bit_change(2))
        ik=min(bit_change(3),bit_change(4))
        il=max(bit_change(3),bit_change(4))
        sig=ch_sign(bra,ii,ij)+ch_sign(ket,ik,il)

        idx=0
        idx=idx+ishft(ii,21)
        idx=idx+ishft(ij,14)
        idx=idx+ishft(ik,7)
        idx=idx+il

        val=0.0
        call fetch_mat(idx,2,val)
        tmp=tmp+val*(-1.0)**sig
        call fetch_mat(idx,5,val1)
        tmp1=tmp1+val1*(-1)**sig
        ! add 3b
        if(num_interaction==3)then
          do ia=1,ncbitn
            idx=0
            val=0.d0
            sig=0
            j=bit_uc(ia)

            if(j<ii)then
              idx=idx+ishft(j,35)
              idx=idx+ishft(ii,28)
              idx=idx+ishft(ij,21)
            sig=sig+ch_sign(bra,int(0,16),j)
              sig=ch_sign(bra,ii,ij)+sig
            elseif(j>ij)then
              idx=idx+ishft(ii,35)
              idx=idx+ishft(ij,28)
              idx=idx+ishft(j,21)
            sig=sig+ch_sign(bra,int(0,16),ii)
              sig=ch_sign(bra,ij,j)+sig
            else
              idx=idx+ishft(ii,35)
              idx=idx+ishft( j,28)
              idx=idx+ishft(ij,21)
            sig=sig+ch_sign(bra,int(0,16),ii)
              sig=ch_sign(bra,j,ij)+sig
            endif

            if(j<ik)then
              idx=idx+ishft(j,14)
              idx=idx+ishft(ik,7)
              idx=idx+il
            sig=sig+ch_sign(ket,int(0,16),j)
            sig=sig+ch_sign(ket,ik,il)
            elseif(j>il)then
              idx=idx+ishft(ik,14)
              idx=idx+ishft(il,7)
              idx=idx+j
            sig=sig+ch_sign(ket,int(0,16),ik)
            sig=sig+ch_sign(ket,il,j)
            else
              idx=idx+ishft(ik,14)
              idx=idx+ishft(j,7)
              idx=idx+il
            sig=sig+ch_sign(ket,int(0,16),ik)
            sig=sig+ch_sign(ket,j,il)
            endif


            call fetch_mat(idx,3,val)
            tmp=tmp+val*(-1)**sig
        call fetch_mat(idx,6,val1)
        tmp1=tmp1+val1*(-1)**sig
          enddo
        endif

      case (6)

        if(num_interaction==3)then
          c1=xor(bra_p , ket_p)
          c2=xor(bra_n , ket_n)
          c1=and(c1,bra_p)
          c2=and(c2,bra_n)
          ie=0
          do ia=1,nn_orb
            if(and(c2,one)==one)then
              ie=ie+1
              bit_change(ie)=ia
            endif
            c2=ishft(c2,-1)
          enddo

          do ia=1,np_orb
            if(and(c1,one)==one)then
              ie=ie+1
              bit_change(ie)=ia+nn_orb
            endif
            c1=ishft(c1,-1)
          enddo
          if(ie/=3)stop 'erro in case 6'

          c1=xor(bra_p , ket_p)
          c2=xor(bra_n , ket_n)
          c1=and(c1 ,ket_p)
          c2=and(c2 ,ket_n)
          do ia=1,nn_orb
            if(and(c2,one)==one)then
              ie=ie+1
              bit_change(ie)=ia
            endif
            c2=ishft(c2,-1)
          enddo

          do ia=1,np_orb
            if(and(c1,one)==one)then
              ie=ie+1
              bit_change(ie)=ia+nn_orb
            endif
            c1=ishft(c1,-1)
          enddo


          if(ie/=6)write(*,*)'fatal error in <bra| ket> 3'

          idx=0
          sig=0
          do ib=1,3
            i=bit_change(ib)
            idx=idx+ishft(i,42-ib*7)
          enddo
        
            sig=sig+ch_sign(bra,int(0,16),bit_change(1))
            sig=sig+ch_sign(bra,bit_change(2),bit_change(3))

          do ib=4,6
            i=bit_change(ib)
            idx=idx+ishft(i,42-ib*7)
          enddo
            sig=sig+ch_sign(ket,int(0,16),bit_change(4))
            sig=sig+ch_sign(ket,bit_change(5),bit_change(6))

          val=0.0
          call fetch_mat(idx,3,val)
          tmp=tmp+val*(-1)**sig
        call fetch_mat(idx,6,val1)
        tmp1=tmp1+val1*(-1)**sig

        endif

        ! pure 3b
      end select


      hm(bra,ket)=tmp
      om(bra,ket)=tmp1


    enddo
  enddo

!write(*,*)om
!write(*,*)hm

end subroutine matrix


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

  call hf_sort_cmplx2(ceig,cvec,cvecl, n)
  do i=1,n
    write(*,*)i,'th state:', ceig(i)
  enddo


end SUBROUTINE lapack_diag


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

subroutine charge_radii(h, cvec,cvecl, ceig, n )
  implicit none
  integer, intent(in) :: n
 integer :: i
  complex*16, dimension(n,n), intent(in) :: h
  COMPLEX*16, DIMENSION(n,n), intent(in) :: cvec,cvecl
  COMPLEX*16, DIMENSION(n), intent(inout) :: ceig
complex*16,dimension(:,:),allocatable :: local_mat
allocate(local_mat(n,n));local_mat=0.d0

local_mat=matmul(transpose(cvecl),h)
local_mat=matmul(local_mat,cvec)
do i=1,n
ceig(i)=local_mat(i,i)
write(*,*)ceig(i)+3.6890138634223710
enddo
deallocate(local_mat)


end subroutine charge_radii
