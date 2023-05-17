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
       use omp_var
!
       use omp_lib
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
!$omp parallel do shared(adj,oldadj,newadj)                            &
!$omp             private(i,j)                                         &
!$omp             schedule(dynamic,chunkscrn)
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
!$omp end parallel do                   
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
       use omp_lib
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
!$omp parallel do shared(life,death,imap,tnagg,tmol,rmol,tiagg,riagg,  &
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
       end subroutine scrnblock
!
!======================================================================!
!
       end module screening
!
!======================================================================!
