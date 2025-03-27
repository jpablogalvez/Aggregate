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
           qiagg = 1
! 
           iimol = timol(isize)
           jimol = rimol(isize)
!
           qimol = jimol + rnmol(isize)
!
           iiagg = 1
!
           do while ( (iiagg.le.tnagg(isize)) .and. (jimol.lt.qimol) )
!
             ni = iimol + 1
             nj = jimol + 1
!
             do jiagg = qiagg, rnagg(isize) 
!
               nj = jimol + 1
!
! Starting comparison of the indexes
!
               if ( tmol(ni) .lt. rmol(nj) ) then
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
                   if ( tmol(iimol+jsize) .ne. rmol(jimol+jsize) )  then   !&
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
                 jimol = jimol + isize
!
               end if 
!
             end do
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
       subroutine calclife(avlife,nlife,nnode,life,nsize,nagg,iagg,    &
                           magg,iwont)
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),dimension(nnode),intent(inout)  ::  avlife   !  Average lifetimes
       integer,dimension(nnode),intent(inout)       ::  nlife    !  
       integer,dimension(nnode),intent(inout)       ::  life     !
       integer,dimension(nnode),intent(in)          ::  nagg     !
       integer,dimension(nnode),intent(in)          ::  iagg     !
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
       end module lifetimes
!
!======================================================================!
