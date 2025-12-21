!======================================================================!
!
       module lifetimes
!
       use omp_var,  only:  np,chunklife
!
       use omp_lib
!
       implicit none
!
       contains
!
!======================================================================!
!
! TRACKLIFE - TRACK LIFEtimes information
!
! This subroutine compares the blocks of two adjacency matrices given as
!  an input list (array representation)
!
       subroutine tracklife(nnode,rsize,rmol,tsize,tmol,rnagg,riagg,   &
                            rnmol,rimol,tnagg,tiagg,timol,life,death,  &
                            imap)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)   ::  rmol   !  Molecules identifier
       integer,dimension(nnode),intent(in)   ::  tmol   !  Molecules identifier
       integer,dimension(nnode),intent(out)  ::  imap   ! 
       integer,dimension(nnode),intent(in)   ::  rnagg  !  
       integer,dimension(nnode),intent(in)   ::  riagg  !
       integer,dimension(nnode),intent(in)   ::  tnagg  !  
       integer,dimension(nnode),intent(in)   ::  tiagg  !
       integer,dimension(nnode),intent(in)   ::  rnmol  !  
       integer,dimension(nnode),intent(in)   ::  rimol  !
       integer,dimension(nnode),intent(in)   ::  timol  !
       integer,intent(in)                    ::  nnode  !  Number of molecules
       integer,intent(in)                    ::  rsize  !  Maximum aggregate size
       integer,intent(in)                    ::  tsize  !  Maximum aggregate size
       logical,dimension(nnode),intent(out)  ::  life   !
       logical,dimension(nnode),intent(out)  ::  death  !
!
! Local variables
! 
       integer                               ::  iiagg  !
       integer                               ::  jiagg  !
       integer                               ::  qiagg  !
       integer                               ::  iimol  !
       integer                               ::  jimol  !
       integer                               ::  qimol  !
       integer                               ::  msize  !
       integer                               ::  isize  !
       integer                               ::  jsize  !
       integer                               ::  ni,nj  !
!
! Comparing the diagonal blocks of two adjacency matrices
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,death,imap)                              &
!$omp             private(ni)                                          &
!$omp             schedule(dynamic,chunklife)
!
       do ni = 1, nnode
         life(ni)  = .FALSE.
         death(ni) = .TRUE.
         imap(ni)  = 0
       end do
!
!$omp end parallel do                   
! 
       msize = tsize
       if ( rsize .lt. tsize ) msize = rsize
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,death,imap,tnagg,tmol,rmol,tiagg,riagg,  &
!$omp                    rnagg,rnmol,rimol,timol,msize)                &
!$omp             private(isize,qiagg,iimol,jimol,qimol,iiagg,ni,nj,   &
!$omp                     jiagg,jsize)                                 &
!$omp             schedule(dynamic,1)
!
       do isize = 2, msize
         if ( (tnagg(isize).ne.0) .and. (rnagg(isize).ne.0) ) then
! 
           iimol = timol(isize)         ! First target molecule of size isize
           jimol = rimol(isize)         ! First target molecule of size isize
!
           qimol = jimol + rnmol(isize) ! Maximum reference molecule of size isize
!
           qiagg = 1                    ! Index of reference aggregate to be compared
           iiagg = 1                    ! Index of target aggregate
!
! Compare aggregates while
!  i) The index of the target aggregate is lower than or equal to the number of target aggregates
!     Run over all the target aggregates
! ii) The index of the reference molecule is lower than the maximum index of the reference molecule
!     Run over all the reference molecules (i.e., reference aggregates)
!     There are no more reference molecules to compare
!    
           do while ( (iiagg.le.tnagg(isize)) .and. (jimol.lt.qimol) )
!
             ni = iimol + 1  ! Target molecule index
             nj = jimol + 1  ! Reference molecule index
!
! Starting from the first reference aggregate
! The first molecule of the first reference aggregate to be compared 
! must be equal or greater than the first target molecule 
!
             do jiagg = qiagg, rnagg(isize) 
!
               nj = jimol + 1
!
! Starting comparison of the indexes
!
               if ( tmol(ni) .lt. rmol(nj) ) then
!
! If first molecule of the target aggregate is lower than the reference molecule to be compared
! it will never appear in the reference aggregates 
!
                 qiagg = jiagg
                 exit 
!
               else if ( tmol(ni) .eq. rmol(nj) ) then
!
! If first molecule index matches then we have a candidate
!
                 life(tiagg(isize)+iiagg) = .TRUE.
                 jsize = 2
                 do while ( (jsize.le.isize) .and.                     &
                                              life(tiagg(isize)+iiagg) )
!
! If all indexes do not match we have a fake positive
!
                   if ( tmol(iimol+jsize) .ne. rmol(jimol+jsize) )  then
                     life(tiagg(isize)+iiagg) = .FALSE.
                   end if
!
                   jsize = jsize + 1
!
                 end do
!
! If all the indexes match we have a positive
!
                 if ( life(tiagg(isize)+iiagg) ) then
                   imap(tiagg(isize)+iiagg) = riagg(isize)+jiagg
                   death(tiagg(isize)+iiagg) = .FALSE.
                 end if
!
                 qiagg = jiagg + 1
                 jimol = jimol + isize
                 exit
!
               else
!
! Move to the first molecule of the next reference aggregate
!
                 jimol = jimol + isize
!
               end if 
!
             end do
!
! Update the target molecule and the target aggregate
!
             iimol = iimol + isize
             iiagg = iiagg + 1
!
           end do
!
         end if
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine tracklife
!
!======================================================================!
!
! CALCLIFE - CALCulate LIFEtimes
!
! This subroutine 
!
       subroutine calclife(nmax,avlife,nlife,nnode,life,nsize,nagg,iagg,    &
                           magg,iwont)
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),dimension(nmax),intent(inout)  ::  avlife   !  Average lifetimes
       integer,dimension(nmax),intent(inout)       ::  nlife    !  
       integer,dimension(nnode),intent(inout)       ::  life     !
       integer,dimension(nnode),intent(in)          ::  nagg     !
       integer,dimension(nnode),intent(in)          ::  iagg     !
       integer,intent(in)                           ::  nmax     ! 
       integer,intent(in)                           ::  nnode    !  Total number of molecules
       integer,intent(in)                           ::  nsize    ! 
       integer,intent(in)                           ::  magg     ! 
       logical,dimension(nnode),intent(in)          ::  iwont    !
!
! Local variables
!
       integer                                      ::  iiagg    !   
       integer                                      ::  inagg    !   
       integer                                      ::  isize    !   
!
! Averaging lifetimes
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,nagg)                                    &
!$omp             private(iiagg)                                       &
!$omp             schedule(dynamic,chunklife)
!
       do iiagg = nagg(1)+1, magg
         life(iiagg) = life(iiagg) + 1 ! FIXME: reduction (?)
       end do
!
!$omp end parallel do 
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(iwont,avlife,life,nlife,nagg,iagg)            &
!$omp             private(isize,inagg,iiagg)                           &
!$omp             schedule(dynamic,1)
!
       do isize = 2, nsize
         do inagg = 1, nagg(isize)
!
           iiagg = iagg(isize) + inagg
!
           if ( iwont(iiagg) ) then
             avlife(isize) = avlife(isize) + life(iiagg) ! FIXME: reduction (?)
             nlife(isize)  = nlife(isize)  + 1           ! FIXME: reduction (?)
             life(iiagg)   = 0
           end if
!
         end do
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine calclife
!
!======================================================================!
!
! NTRACKLIFE - N-components TRACK LIFEtimes information
!
! This subroutine compares the blocks of two adjacency matrices given as
!  an input list (array representation)
!
       subroutine ntracklife(mtype,nnode,nmax,nmon,rmidx,rsize,rmol,   &
                             tmidx,tsize,tmol,rnagg,riagg,rnmol,rimol, &
                             tnagg,tiagg,timol,life,death,imap)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nmax),intent(in)    ::  nmon   !  
       integer,dimension(nnode),intent(in)   ::  rmol   !  Molecules identifier
       integer,dimension(nnode),intent(in)   ::  tmol   !  Molecules identifier
       integer,dimension(nnode),intent(out)  ::  imap   ! 
       integer,dimension(nmax),intent(in)    ::  rnagg  !  
       integer,dimension(nmax),intent(in)    ::  riagg  !
       integer,dimension(nmax),intent(in)    ::  tnagg  !  
       integer,dimension(nmax),intent(in)    ::  tiagg  !
       integer,dimension(nmax),intent(in)    ::  rnmol  !  
       integer,dimension(nmax),intent(in)    ::  rimol  !
       integer,dimension(nmax),intent(in)    ::  timol  !
       integer,intent(in)                    ::  mtype  ! 
       integer,intent(in)                    ::  nnode  !  Number of molecules
       integer,intent(in)                    ::  nmax   !  Number of possible aggregates 
       integer,intent(in)                    ::  rmidx  !  Maximum aggregate identifier
       integer,intent(in)                    ::  tmidx  !  Maximum aggregate identifier
       integer,intent(in)                    ::  rsize  !  Maximum aggregate size
       integer,intent(in)                    ::  tsize  !  Maximum aggregate size
       logical,dimension(nnode),intent(out)  ::  life   !
       logical,dimension(nnode),intent(out)  ::  death  !
!
! Local variables
! 
       integer                               ::  iiagg  !
       integer                               ::  jiagg  !
       integer                               ::  qiagg  !
       integer                               ::  iimol  !
       integer                               ::  jimol  !
       integer                               ::  qimol  !
       integer                               ::  msize  !
       integer                               ::  isize  !
       integer                               ::  jsize  !
       integer                               ::  ni,nj  !
!
! Comparing the diagonal blocks of two adjacency matrices
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,death,imap)                              &
!$omp             private(ni)                                          &
!$omp             schedule(dynamic,chunklife)
!
       do ni = 1, nnode
         life(ni)  = .FALSE.
         death(ni) = .TRUE.
         imap(ni)  = 0
       end do
!
!$omp end parallel do                   
! 
       msize = tmidx
       if ( rmidx .lt. tmidx ) msize = rmidx 
       if ( nmax .lt. msize ) msize = nmax
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,death,imap,tnagg,tmol,rmol,tiagg,riagg,  &
!$omp                    rnagg,rnmol,rimol,timol,msize,nmon)           &
!$omp             private(isize,qiagg,iimol,jimol,qimol,iiagg,ni,nj,   &
!$omp                     jiagg,jsize)                                 &
!$omp             schedule(dynamic,1)
!
       do isize = mtype+1, msize
         if ( (tnagg(isize).ne.0) .and. (rnagg(isize).ne.0) ) then
! 
           iimol = timol(isize)         ! First target molecule of size isize
           jimol = rimol(isize)         ! First target molecule of size isize
!
           qimol = jimol + rnmol(isize) ! Maximum reference molecule of size isize
!
           qiagg = 1                    ! Index of reference aggregate to be compared
           iiagg = 1                    ! Index of target aggregate
!
! Compare aggregates while
!  i) The index of the target aggregate is lower than or equal to the number of target aggregates
!     Run over all the target aggregates
! ii) The index of the reference molecule is lower than the maximum index of the reference molecule
!     Run over all the reference molecules (i.e., reference aggregates)
!     There are no more reference molecules to compare
!    
           do while ( (iiagg.le.tnagg(isize)) .and. (jimol.lt.qimol) )
!
             ni = iimol + 1  ! Target molecule index
             nj = jimol + 1  ! Reference molecule index
!
! Starting from the first reference aggregate
! The first molecule of the first reference aggregate to be compared 
! must be equal or greater than the first target molecule 
!
             do jiagg = qiagg, rnagg(isize) 
!
               nj = jimol + 1
!
! Starting comparison of the indexes
!
               if ( tmol(ni) .lt. rmol(nj) ) then
!
! If first molecule of the target aggregate is lower than the reference molecule to be compared
! it will never appear in the reference aggregates 
!
                 qiagg = jiagg
                 exit 
!
               else if ( tmol(ni) .eq. rmol(nj) ) then
!
! If first molecule index matches then we have a candidate
!
                 life(tiagg(isize)+iiagg) = .TRUE.
                 jsize = 2
                 do while ( (jsize.le.nmon(isize)) .and.               &
                                              life(tiagg(isize)+iiagg) )
!
! If all indexes do not match we have a fake positive
!
                   if ( tmol(iimol+jsize) .ne. rmol(jimol+jsize) )  then
                     life(tiagg(isize)+iiagg) = .FALSE.
                   end if
!
                   jsize = jsize + 1
!
                 end do
!
! If all the indexes match we have a positive
!
                 if ( life(tiagg(isize)+iiagg) ) then
                   imap(tiagg(isize)+iiagg) = riagg(isize)+jiagg
                   death(tiagg(isize)+iiagg) = .FALSE.
                 end if
!
                 qiagg = jiagg + 1
                 jimol = jimol + nmon(isize)
                 exit
!
               else
!
! Move to the first molecule of the next reference aggregate
!
                 jimol = jimol + nmon(isize)
!
               end if 
!
             end do
!
! Update the target molecule and the target aggregate
!
             iimol = iimol + nmon(isize)
             iiagg = iiagg + 1
!
           end do
!
         end if
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine ntracklife
!
!======================================================================!
!
! NCALCLIFE - N-components CALCulate LIFEtimes
!
! This subroutine 
!
       subroutine ncalclife(mtype,nnode,nmax,avlife,nlife,life,nsize,  &
                            nagg,iagg,magg,iwont)
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),dimension(nmax),intent(inout)   ::  avlife   !  Average lifetimes
       integer,dimension(nmax),intent(inout)        ::  nlife    !  
       integer,dimension(nnode),intent(inout)       ::  life     !
       integer,dimension(nmax),intent(in)           ::  nagg     !
       integer,dimension(nmax),intent(in)           ::  iagg     !
       integer,intent(in)                           ::  mtype    ! 
       integer,intent(in)                           ::  nmax     !  
       integer,intent(in)                           ::  nnode    !  Total number of molecules
       integer,intent(in)                           ::  nsize    ! 
       integer,intent(in)                           ::  magg     ! 
       logical,dimension(nnode),intent(in)          ::  iwont    !
!
! Local variables
!
       integer                                      ::  iiagg    !   
       integer                                      ::  inagg    !   
       integer                                      ::  isize    !   
!
! Averaging lifetimes
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,nagg)                                    &
!$omp             private(iiagg)                                       &
!$omp             schedule(dynamic,chunklife)
!
       do iiagg = nagg(mtype+1)+1, magg
         life(iiagg) = life(iiagg) + 1 ! FIXME: reduction (?)
       end do
!
!$omp end parallel do 
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(iwont,avlife,life,nlife,nagg,iagg)            &
!$omp             private(isize,inagg,iiagg)                           &
!$omp             schedule(dynamic,1)
!
       do isize = mtype+1, nsize
         do inagg = 1, nagg(isize)
!
           iiagg = iagg(isize) + inagg
!
           if ( iwont(iiagg) ) then
             avlife(isize) = avlife(isize) + life(iiagg) ! FIXME: reduction (?)
             nlife(isize)  = nlife(isize)  + 1           ! FIXME: reduction (?)
             life(iiagg)   = 0
           end if
!
         end do
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine ncalclife
!
!======================================================================!
!
       end module lifetimes
!
!======================================================================!
