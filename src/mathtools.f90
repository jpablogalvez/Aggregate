!======================================================================!
!
       module mathtools
!
       implicit none
!
       contains
!
!======================================================================!
!
!  NTUPLA2IDX - N-TUPLA to componentes InDeX
!
! This function 
!
       function ntupla2idx(ntype,idx) result(nidx)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(ntype),intent(in)  ::  idx     !
       integer,intent(in)                   ::  ntype   !
       integer                              ::  nidx    !
!
! Local variables
!
       integer                              ::  suma    !
       integer                              ::  prod    !
       integer                              ::  fact    !
       integer                              ::  i,j     !
!
! 
!
       nidx = 0
!
       do i = 1, ntype
!
         suma = 0
         do j = 1, i
           suma = suma + idx(j)
         end do
!
         prod = 1
         fact = 1
         do j = 1, i
           prod = prod*(suma+j-1)
           fact = fact*j
         end do
         prod = prod/fact
!
         nidx = nidx + prod
!
       end do
!    
       return
       end function ntupla2idx
!
!======================================================================!
!
! SETIDX - SET InDeX
!
! This subroutine 
!
       subroutine setidx(msize,mtype,nmax,mmon,nmon,imon,mgrpsmon,     & 
                         ngrpsmon,igrpsmon,mbodymon,nbodymon,ibodymon)
!
       use systeminf,  only:  rep
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nmax),intent(out)        ::  mmon      !
       integer,dimension(nmax),intent(out)        ::  mgrpsmon  !
       integer,dimension(nmax),intent(out)        ::  mbodymon  !
       integer,dimension(mtype,nmax),intent(out)  ::  nmon      !
       integer,dimension(mtype,nmax),intent(out)  ::  imon      !
       integer,dimension(mtype,nmax),intent(out)  ::  ngrpsmon  !
       integer,dimension(mtype,nmax),intent(out)  ::  igrpsmon  !
       integer,dimension(mtype,nmax),intent(out)  ::  nbodymon  !
       integer,dimension(mtype,nmax),intent(out)  ::  ibodymon  !
       integer,intent(in)                         ::  mtype     !
       integer,intent(in)                         ::  msize     !
       integer,intent(in)                         ::  nmax      !
!
! Local variables
!
       integer,dimension(mtype)                   ::  ntype     !
       integer                                    ::  num       !
       integer                                    ::  fact      !
       integer                                    ::  itype     !
       integer                                    ::  i,j       !
!
! Setting the number of monomer molecules for each aggregate identifier 
!
       mmon(:) = -1
!
       j = 0
       do i = 1, msize
!
         num  = 1
         fact = 1
         do itype = 1, mtype-1
           num  = num*(mtype+i-itype)
           fact = fact*itype
         end do
         num = num/fact
!
         mmon(j+1:j+num) = i
         j = j + num
!
       end do
!
! Generating all possible n-tuples from combinations with repetition 
!
       nmon(:,:) = -1
       ntype(:)  = 0
!
       call genntuples(1,mtype,ntype,msize,nmax,nmon)
!
! Setting the information of the aggregates size
!
       mgrpsmon(:) = -1
       mbodymon(:) = -1
!
       ngrpsmon(:,:) = -1
       nbodymon(:,:) = -1
!
       igrpsmon(:,:) = -1
       ibodymon(:,:) = -1
!
       do i = 1, nmax-1
!
         mgrpsmon(i) = 0
         mbodymon(i) = 0
!
         igrpsmon(1,i) = 0
         ibodymon(1,i) = 0 
!
         imon(1,i) = 0        
!
         do j = 1, mtype-1
!
           ngrpsmon(j,i) = nmon(j,i)*rep(j)%mgrps
           nbodymon(j,i) = nmon(j,i)*rep(j)%mbody 
!
           mgrpsmon(i) = mgrpsmon(i) + ngrpsmon(j,i)
           mbodymon(i) = mbodymon(i) + nbodymon(j,i)
!
           igrpsmon(j+1,i) = igrpsmon(j,i) + ngrpsmon(j,i)
           ibodymon(j+1,i) = ibodymon(j,i) + nbodymon(j,i) 
!
           imon(j+1,i) = imon(j,i) + nmon(j,i)
!
         end do
!
         ngrpsmon(mtype,i) = nmon(mtype,i)*rep(mtype)%mgrps
         nbodymon(mtype,i) = nmon(mtype,i)*rep(mtype)%mbody
!
         mgrpsmon(i) = mgrpsmon(i) + ngrpsmon(mtype,i)
         mbodymon(i) = mbodymon(i) + nbodymon(mtype,i)
! 
       end do
!    
       return
       end subroutine setidx
!
!======================================================================!
!
! NEXTCR - NEXT Combination with Repetition
!
       recursive subroutine genntuples(posi,mtype,ntype,msize,nmax,imon)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(mtype),intent(inout)        ::  ntype  !
       integer,dimension(mtype,nmax),intent(inout)   ::  imon   !
       integer,intent(in)                            ::  mtype  !
       integer,intent(in)                            ::  msize  !
       integer,intent(in)                            ::  nmax   !
       integer,intent(in)                            ::  posi   !
!
! Local variables
!
       integer                                       ::  i,j
!
!
!
       if ( posi .gt. mtype ) then
!
         if ( (sum(ntype).le.msize) .and. (sum(ntype).gt.0) ) then
           j = ntupla2idx(mtype,ntype)
           imon(:,j) = ntype(:)
         end if
!
       else
!
         do i = 0, msize
           ntype(posi) = i
           call genntuples(posi+1,mtype,ntype,msize,nmax,imon)
         end do
!
       end if
!    
       return
       end subroutine genntuples
!
!======================================================================!
!
! NEXTP - NEXT Permutation
!
! Adapted from: http://rosetta code.org/wiki/Permutations
!
       logical function nextp(n,a)
!
       implicit none
!
       integer,intent(in)                  ::  n
       integer,dimension(n),intent(inout)  ::  a
!
!     local variables:
!
       integer i,j,k,t
!
       i = n-1
10     if ( a(i) .lt. a(i+1) ) GOTO 20
       i = i-1
       if ( i .eq. 0 ) GOTO 20
       GOTO 10
20     j = i+1
       k = n
30     t = a(j)
       a(j) = a(k)
       a(k) = t
       j = j+1
       k = k-1
       if ( j .lt. k ) GOTO 30
       j = i
       if (j .ne. 0 ) GOTO 40
!      
       nextp = .false.
!      
       return
!
40     j = j+1
       if ( a(j) .lt. a(i) ) GOTO 40
       t = a(i)
       a(i) = a(j)
       a(j) = t
!      
       nextp = .true.
!      
       return
       end function nextp
!
!======================================================================!
!
       function factorial(n) result(fact)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)  ::  n
       integer             ::  fact
!
! Local variables
!
       integer             ::  i
!
! Computing the factorial of N
!
       fact = 1
!
       do i = 1, n
         fact = fact*i
       end do
!
       end function factorial
!
!======================================================================!
!
       end module mathtools
!
!======================================================================!
