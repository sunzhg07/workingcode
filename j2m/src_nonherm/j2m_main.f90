program main
  use j2m_input
  use single_particle_orbits
  use channel_info
  use ang_mom_functions
  implicit none
  call read_input
  call setup_sp_orbits
  call setup_channels
  call read_jcoupled_interactions
  call commons_to_angmom
  call hermitelize

  stop
  call setup_mscheme_interactions

  
end program main


subroutine setup_sp_orbits
  use single_particle_orbits
  use j2m_input
  implicit none
  integer:: n,i,nm,mz
  integer:: k,nn,ll,jj,tt
  character(len=20):: info
  open(unit=10, file=jjsp_file)
  read(10,*)n
  njobs_tot=n
  call allocate_sp_array(jsp,n)
  nm=0
  do i=1,n
    read(10,*)k, nn,ll,jj,tt,info
    jsp%nn(i)=nn
    jsp%ll(i)=ll
    jsp%jj(i)=jj
    jsp%itzp(i)=tt
    jsp%jorder(i)=i
    nm=nm+jj+1
  enddo

  call allocate_sp_array(msp,nm)

  nm=0
  do i=1,n
    nn=jsp%nn(i)
    ll=jsp%ll(i)
    jj=jsp%jj(i)
    tt=jsp%itzp(i)
    if(tt==-1)cycle
    do mz=-jj,jj,2
      nm=nm+1
    msp%nn(nm)=nn
    msp%ll(nm)=ll
    msp%jj(nm)=jj
    msp%itzp(nm)=tt
    msp%mm(nm)=mz
    msp%jorder(nm)=i
    enddo
  enddo

  do i=1,n
    nn=jsp%nn(i)
    ll=jsp%ll(i)
    jj=jsp%jj(i)
    tt=jsp%itzp(i)
    if(tt==1)cycle
    do mz=-jj,jj,2
      nm=nm+1
    msp%nn(nm)=nn
    msp%ll(nm)=ll
    msp%jj(nm)=jj
    msp%itzp(nm)=tt
    msp%mm(nm)=mz
    msp%jorder(nm)=i
    enddo
  enddo
  nmobs_tot=nm

end subroutine setup_sp_orbits


subroutine read_jcoupled_interactions
  use single_particle_orbits
  use operators
  use channel_info
  use j2m_input
  implicit none
  integer:: n,a,b,c,d,i,jt,ipar,itz,channel
  integer:: bra,ket,ja,jb,jc,jd
  real*8:: val,fact
  real*8:: dij
  integer:: iph

  open(unit=11,file=jj1b_file)
  open(unit=12,file=jj2b_file)

  read(11,*)n
  do i=1,n
    read(11,*)a,b,val
    vj1b(a,b)=val
  enddo
  close(11)

  read(12,*)n
  do i=1,n
    read(12,*)a,b,c,d,jt,val
    ipar=mod(jsp%ll(a)+jsp%ll(b),2)
    itz=(jsp%itzp(a)+jsp%itzp(b))/2
    channel=lookup_jchannel(jt,ipar,itz)
    if(channel==0)stop
    bra=jab_configs(channel)%ival(a,b)
    ket=jab_configs(channel)%ival(c,d)
  
    if(normalized)then
    vj2b%op(channel)%val(bra,ket)=val
  else
    vj2b%op(channel)%val(bra,ket)=val*dij(a,b)*dij(c,d)

  endif
    if(bra/=ket .and. hermit )then
    vj2b%op(channel)%val(ket,bra)=vj2b%op(channel)%val(bra,ket)
    endif

    if(itz==0)then
    bra=jab_configs(channel)%ival(b,a)
    ket=jab_configs(channel)%ival(d,c)
    ja=jsp%jj(a)
    jb=jsp%jj(b)
    jc=jsp%jj(c)
    jd=jsp%jj(d)
    vj2b%op(channel)%val(bra,ket)=val*iph((ja+jb+jc+jd)/2)
    if(hermit)then
    vj2b%op(channel)%val(ket,bra)=val*iph((ja+jb+jc+jd)/2)
  endif
    endif
  enddo

!@!  do channel=1,vj2b%len
!@!    n=size(lookup_jab_configs(channel)%ival,2)
!@!    do bra=1,n
!@!      a=lookup_jab_configs(channel)%ival(1,bra)
!@!      b=lookup_jab_configs(channel)%ival(2,bra)
!@!      do ket=1,n
!@!      c=lookup_jab_configs(channel)%ival(1,ket)
!@!      d=lookup_jab_configs(channel)%ival(2,ket)
!@!      val=vj2b%op(channel)%val(bra,ket)
!@!      if(abs(val)<1E-6)cycle
!@!      write(15,*)a,b,c,d,val
!@!      enddo
!@!    enddo
!@!  enddo

end subroutine read_jcoupled_interactions




subroutine setup_mscheme_interactions
  use single_particle_orbits
  use operators
  use channel_info
  use j2m_input
  use ang_mom_functions
  implicit none
  real*8:: angmon_fact,val
  integer:: a,b,c,d,ja,jb,jc,jd,la,lb,lc,ld,ita,itb,itc,itd,ma,mb,mc,md
  integer:: jmin,jmax,jt
  integer:: channel,bra,ket,bra1,ket1,channel1
  integer:: ia,ib,ic,id,ipar,itz,mt,n

  do channel=1,vm2b%len
    ipar=m_channels(channel)%ipar
    itz=m_channels(channel)%itz
    mt=m_channels(channel)%mm
    do bra=1,size(lookup_mab_configs(channel)%ival,2)
      ia=lookup_mab_configs(channel)%ival(1,bra)
      ib=lookup_mab_configs(channel)%ival(2,bra)
      ma=msp%mm(ia)
      mb=msp%mm(ib)
      ja=msp%jj(ia)
      jb=msp%jj(ib)
      a=msp%jorder(ia)
      b=msp%jorder(ib)
    do ket=1,size(lookup_mab_configs(channel)%ival,2)
      ic=lookup_mab_configs(channel)%ival(1,ket)
      id=lookup_mab_configs(channel)%ival(2,ket)
      mc=msp%mm(ic)
      md=msp%mm(id)
      jc=msp%jj(ic)
      jd=msp%jj(id)
      c=msp%jorder(ic)
      d=msp%jorder(id)

      jmin=max(abs(ja-jb)/2,abs(jc-jd)/2)
      jmax=min(abs(ja+jb)/2,abs(jc+jd)/2)
      val=0.d0
      do jt=jmin,jmax
        channel1=lookup_jchannel(jt,ipar,itz)
        if(channel1==0)cycle
        bra1=jab_configs(channel1)%ival(a,b)
        ket1=jab_configs(channel1)%ival(c,d)
        angmon_fact=cgc(ja,jb,2*jt,ma,mb,2*mt)*&
                cgc(jc,jd,2*jt,mc,md,2*mt)
        val=val+angmon_fact*vj2b%op(channel1)%val(bra1,ket1)
      enddo
      vm2b%op(channel)%val(bra,ket)=val
    enddo
  enddo

  enddo

  open(unit=13,file=m1b_file)
  open(unit=14,file=m2b_file)

  n=0
  do ia=1,nmobs_tot
     a=msp%jorder(ia)
     do ib=1,nmobs_tot
       b=msp%jorder(ib)
       if(msp%ll(ia)/=msp%ll(ib))cycle
       if(msp%mm(ia)/=msp%mm(ib))cycle
       if(msp%itzp(ia)/=msp%itzp(ib))cycle
       val=vj1b(a,b)
       if(abs(val)<1E-6)cycle
       n=n+1
     enddo
     enddo
  write(13,*)n
  do ia=1,nmobs_tot
     a=msp%jorder(ia)
     do ib=1,nmobs_tot
       b=msp%jorder(ib)
       if(msp%ll(ia)/=msp%ll(ib))cycle
       if(msp%mm(ia)/=msp%mm(ib))cycle
       if(msp%itzp(ia)/=msp%itzp(ib))cycle
       val=vj1b(a,b)
       if(abs(val)<1E-6)cycle
       write(13,'(2(I4,1X),F13.6)')ia,ib,val
     enddo
     enddo
  


  n=0
  do channel=1,vm2b%len
    do bra=1,size(lookup_mab_configs(channel)%ival,2)
    do ket=1,size(lookup_mab_configs(channel)%ival,2)
      val=vm2b%op(channel)%val(bra,ket)
      if(abs(val)<1E-6)cycle
      n=n+1
    enddo
  enddo
  enddo

  write(14,*)n

  do channel=1,vm2b%len
    do bra=1,size(lookup_mab_configs(channel)%ival,2)
      ia=lookup_mab_configs(channel)%ival(1,bra)
      ib=lookup_mab_configs(channel)%ival(2,bra)
    do ket=1,size(lookup_mab_configs(channel)%ival,2)
      ic=lookup_mab_configs(channel)%ival(1,ket)
      id=lookup_mab_configs(channel)%ival(2,ket)
      val=vm2b%op(channel)%val(bra,ket)
      if(abs(val)<1E-6)cycle
       write(14,'(4(I4,1X),F13.6)')ia,ib,ic,id,val
    enddo
  enddo
  enddo
  close(13)
  close(14)

  open(unit=15,file=msp_file)

  write(15,*)nmobs_tot
  do n=1,nmobs_tot
    write(15,'(7(I4,2X),A12)')n,msp%nn(n),&
      msp%ll(n),&
      msp%jj(n),&
      msp%mm(n),&
      msp%itzp(n),&
      msp%jorder(n),&
      '1  particle'
  enddo

end subroutine setup_mscheme_interactions


subroutine hermitelize
  use single_particle_orbits
  use operators
  use channel_info
  use j2m_input
  implicit none
  integer:: channel,bra,ket,ia,ib,ic,id,ja,jb,jc,jd
  integer:: bra1,ket1,jt
  real*8:: dij
  real*8:: tmp
  integer:: nconfs,ndim
  integer,allocatable:: valid_conf(:,:)
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: hh,v,ed,rvec,lvec,u,u_inv,hback
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: ee
  complex*16:: dd

  do channel=1, vj2b%len
    ndim=0
    nconfs=size(lookup_jab_configs(channel)%ival,2)
    do bra=1,nconfs
      ia=lookup_jab_configs(channel)%ival(1,bra)
      ib=lookup_jab_configs(channel)%ival(2,bra)
      if(ia>ib)cycle
      ndim=ndim+1
    enddo
    allocate(valid_conf(3,ndim))


    ndim=0
    do bra=1,nconfs
      ia=lookup_jab_configs(channel)%ival(1,bra)
      ib=lookup_jab_configs(channel)%ival(2,bra)
      if(ia>ib)cycle
      ndim=ndim+1
      valid_conf(1:2,ndim)=lookup_jab_configs(channel)%ival(:,bra)
      valid_conf(3,ndim)=bra
    enddo

    allocate(u(ndim,ndim));u=0.d0
    allocate(u_inv(ndim,ndim));u_inv=0.d0
    allocate(hh(ndim,ndim));hh=0.d0
    allocate(hback(ndim,ndim));hback=0.d0
    allocate(ed(ndim,ndim));ed=0.d0
    allocate(rvec(ndim,ndim));rvec=0.d0
    allocate(lvec(ndim,ndim));lvec=0.d0
    allocate(ee(ndim)); ee=0.d0
    do bra=1,ndim
      ia=valid_conf(1,bra)
      ib=valid_conf(2,bra)
      bra1=valid_conf(3,bra)
      do ket=1,ndim
      ic=valid_conf(1,ket)
      id=valid_conf(2,ket)
      ket1=valid_conf(3,ket)

      tmp=0.d0
      
      tmp=tmp+vj2b%op(channel)%val(bra1,ket1)/dij(ia,ib)/dij(ic,id)

                  if (bra1 == ket1) then
                    !! suppose 1bd off diagonal is 0
                     tmp = tmp + vj1b(ia, ia) + vj1b(ib, ib)
                     ed(bra,ket)=vj1b(ia, ia) + vj1b(ib, ib)
                   endif
          hh(bra,ket)=tmp

    enddo
  enddo
  hback=hh
  call lapack_diag(hh,rvec,lvec, ee, ndim)

  call cmplxmatinv(rvec,ndim,dd)
  
  lvec=0.
  lvec=matmul(transpose(rvec),rvec)
  
  call sqrtmat(lvec,u,u_inv ,ndim)

  lvec=matmul(u,hback)
  hh=matmul(lvec,u_inv)

!  call lapack_diag(hh,rvec,lvec, ee, ndim)
  hh=hh-ed

  jt=j_channels(channel)%jj
    do bra=1,ndim
      ia=valid_conf(1,bra)
      ib=valid_conf(2,bra)
      do ket=1,ndim
      ic=valid_conf(1,ket)
      id=valid_conf(2,ket)
      write(999,'(5(I4,2x),1F13.6)')ia,ib,ic,id,jt,real(hh(bra,ket))*dij(ia,ib)*dij(ic,id)
    enddo
  enddo
  

  DEALLOCATE(hh,rvec,lvec,ee,ed,u,u_inv,hback)
      

    deallocate(valid_conf)
  enddo
  end subroutine hermitelize
