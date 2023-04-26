!======================================================================!
!
       module screening
       implicit none
!
       contains
!
!======================================================================!
!
! SCRNINT - SCReeNing INTeractions
!
! This subroutine determines if any interaction between two molecules is
!  present in the previous and current snapshots or in the current and
!  next snapshots
!
       subroutine scrnint(nnode,oldadj,adj,newadj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(inout)  ::  adj     !  Adjacency matrix of the current snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  oldadj  !  Adjacency matrix of the previous snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  newadj  !  Adjacency matrix of the next snapshot
       integer,intent(in)                            ::  nnode   !  Number of molecules
!
! Local variables
! 
       integer                                       ::  i,j     !
!
! Comparing the diagonal blocks of two adjacency matrices
!
       do i = 1, nnode-1
         do j = i+1, nnode
           if ( (oldadj(j,i).and.adj(j,i)) .or.                        &
                (adj(j,i).and.newadj(j,i)) .or.                        &
                                    (oldadj(j,i).and.newadj(j,i)) ) then
             adj(j,i) = .TRUE.
             adj(i,j) = .TRUE.
           else
             adj(j,i) = .FALSE.
             adj(i,j) = .FALSE.
           end if
         end do
       end do
!
       return
       end subroutine scrnint
!
!======================================================================!
!
! SCRNBLOCK - SCReeNing BLOCKs
!
! This subroutine compares the blocks of two adjacency matrices given as
!  an input list (array representation)
!
       subroutine scrnblock(nnode,rsize,rmol,tsize,tmol,rnagg,riagg,   &
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
       life(:)  = .FALSE.
       death(:) = .TRUE.
       imap(:)  = 0
! 
       msize = tsize
       if ( rsize .lt. tsize ) msize = rsize
!
       do isize = 2, msize
         if ( (tnagg(isize).ne.0) .and. (rnagg(isize).ne.0) ) then
!
!write(*,*) 'Starting aggregates of size',isize,':',tnagg(isize),rnagg(isize)
!write(*,*) '----------------------------------------'
!write(*,*)
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
!write(*,'(1X,A,1X,I6,1X,A,1X,I6,1X,A,1X,I6)') 'AGGREGATE  ',iiagg,'starts at',ni,':',tmol(ni)
!~ write(*,'(1X,A,7X,I6,1X,A,1X,4I6)') '  COMPARISON starts at',nj,':',rmol(nj),jimol,qimol,qiagg 
!
             do jiagg = qiagg, rnagg(isize) 
!
               nj = jimol + 1
!
!write(*,'(1X,A,3X,I6,1X,A,1X,I6,1X,A,I6,1X,I6)') &
!'  COMPARING',ni,'WITH',nj,':',tmol(ni),rmol(nj)
!
               if ( tmol(ni) .lt. rmol(nj) ) then
!
                 qiagg = jiagg
                 exit 
!
               else if ( tmol(ni) .eq. rmol(nj) ) then
!
!write(*,*) '    WE HAVE A CANDIDATE'
!
                 life(tiagg(isize)+iiagg) = .TRUE.
                 jsize = 2
                 do while ( (jsize.le.isize) .and.                     &
                                              life(tiagg(isize)+iiagg) )
!
!write(*,'(2(1X,A,1X,I6),1X,A,I6,1X,I6)')  &
!'    COMPARING',iimol+jsize,'WITH',jimol+jsize,':',tmol(iimol+jsize),rmol(jimol+jsize)
!
                   if ( tmol(iimol+jsize) .ne. rmol(jimol+jsize) )  then   !&
                     life(tiagg(isize)+iiagg) = .FALSE.
!
!write(*,*) '      FAKE POSITIVE'   ! FLAG: check for partial matching
!
                   end if
                   jsize = jsize + 1
                 end do
!
                 if ( life(tiagg(isize)+iiagg) ) then
                   imap(tiagg(isize)+iiagg) = riagg(isize)+jiagg
                   death(tiagg(isize)+iiagg) = .FALSE.
!
!write(*,*) '      WE HAVE A POSITIVE',tiagg(isize)+iiagg,'is',tiagg(isize)+jiagg
!
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
!write(*,*)
!
             iimol = iimol + isize
             iiagg = iiagg + 1
!
           end do
!
         end if
       end do
!
       return
       end subroutine scrnblock
!
!======================================================================!
!
       end module screening
!
!======================================================================!
