!======================================================================!
!
       module screening
       implicit none
!
       contains
!
!======================================================================!
!
! SCRNALRD - SCReeNing ALReaDy
!
! This subroutine 
!
!
       subroutine array2list(nnode,mol,agg,nsize,nagg,imol,            &
                             head,prev,list,next)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)   ::  mol      !  Molecules identifier
       integer,dimension(nnode),intent(in)   ::  agg      !  Aggregates size
       integer,dimension(nnode),intent(out)  ::  head     !
       integer,dimension(nnode),intent(out)  ::  list     !  
       integer,dimension(nnode),intent(out)  ::  prev     !
       integer,dimension(nnode),intent(out)  ::  next     !
       integer,dimension(nnode),intent(in)   ::  nagg     !
       integer,dimension(nnode),intent(in)   ::  imol     !    
       integer,intent(in)                    ::  nnode    !  Number of molecules
       integer,intent(in)                    ::  nsize    !  Maximum aggregate size
!
! Local variables
! 
       integer                               ::  oldnext  !
       integer                               ::  oldprev  !
       integer                               ::  isize    !  
       integer                               ::  iiagg    ! 
       integer                               ::  iimol    ! 
       integer                               ::  ni       !
!
! Changing from array representation to liked list representation
!
       head(:nsize) = 0
!
       do isize = 2, nsize
         if ( nagg(isize) .gt. 0 ) then
!
           ni = imol(isize)
!
           oldnext = 0
           oldprev = mol(ni+1)
!
           do iimol = 1, isize
!
             ni = ni + 1
!
             list(mol(ni)) = head(isize)
             head(isize)   = mol(ni)
!
           end do
!
           next(mol(ni)) = oldnext
           oldnext       = mol(ni)
!
           prev(oldprev) = mol(ni)
!
           if ( nagg(isize) .gt. 1 ) then
             do iiagg = 2, nagg(isize)
!
               oldprev = mol(ni+1)
!
               do iimol = 1, isize
!
                 ni = ni + 1
!
                 list(mol(ni)) = head(isize)
                 head(isize)   = mol(ni)
!
               end do
!
               prev(oldnext) = oldprev
!
               next(mol(ni)) = oldnext
               oldnext       = mol(ni)
! 
               prev(oldprev) = mol(ni)
!
             end do
           end if
!
           prev(mol(ni)) = 0
!
         end if
       end do
!
!~        head(:)  = 0
!~ !
!~        do isize = 2, nsize
!~          ni = imol(isize)
!~          do iiagg = 1, nagg(isize)
!~            do iimol = 1, isize
!~ !
!~              ni = ni + 1
!~ !
!~              list(mol(ni)) = head(isize)
!~              head(isize)   = mol(ni)          
!~ !
!~            end do
!~          end do
!~        end do
!
!!~        do isize = 2, msize
!!~          do iimol = 1, nmol(isize)
!!~ !
!!~            ni = imol(isize)+iimol
!!~ !
!!~            tail(ni)      = list(ni)
!!~            list(ni)      = head(agg(ni))
!!~            head(agg(ni)) = ni          
!!~ !
!!~          end do
!!~        end do
!
       return
       end subroutine array2list
!
!======================================================================!
!
! SCRNBLOCK - SCReeNing BLOCKs
!
! This subroutine compares the blocks of two adjacency matrices given as
!  an input list
!
!
       subroutine scrnblock(nnode,rsize,rmol,ragg,tsize,tmol,tagg,     &
                            rnagg,riagg,rnmol,rimol,tnagg,tiagg,tnmol, &
                            timol,life,death,imap)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)   ::  rmol   !  Molecules identifier
       integer,dimension(nnode),intent(in)   ::  ragg   !  Aggregates size
       integer,dimension(nnode),intent(in)   ::  tmol   !  Molecules identifier
       integer,dimension(nnode),intent(in)   ::  tagg   !  Aggregates size      
       integer,dimension(nnode),intent(out)  ::  imap   ! 
       integer,dimension(nnode),intent(in)   ::  rnagg  !  
       integer,dimension(nnode),intent(in)   ::  riagg  !
       integer,dimension(nnode),intent(in)   ::  tnagg  !  
       integer,dimension(nnode),intent(in)   ::  tiagg  !
       integer,dimension(nnode),intent(in)   ::  rnmol  !  
       integer,dimension(nnode),intent(in)   ::  rimol  !
       integer,dimension(nnode),intent(in)   ::  tnmol  !  
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
! SCRNALRD - SCReeNing ALReaDy
!
! This subroutine 
!
!
       subroutine scrnalrd(nnode,rsize,rmol,ragg,tsize,tmol,tagg,      &
                           rnagg,riagg,rnmol,rimol,tnagg,tiagg,tnmol,  &
                           timol,rlife,rdeath,tlife,tdeath)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)  ::  rmol    !  Molecules identifier
       integer,dimension(nnode),intent(in)  ::  ragg    !  Aggregates size
       integer,dimension(nnode),intent(in)  ::  tmol    !  
       integer,dimension(nnode),intent(in)  ::  tagg    ! 
       integer,dimension(nnode),intent(in)  ::  rnagg   !  
       integer,dimension(nnode),intent(in)  ::  riagg   !
       integer,dimension(nnode),intent(in)  ::  tnagg   !  
       integer,dimension(nnode),intent(in)  ::  tiagg   !
       integer,dimension(nnode),intent(in)  ::  rnmol   !  
       integer,dimension(nnode),intent(in)  ::  rimol   !
       integer,dimension(nnode),intent(in)  ::  tnmol   !  
       integer,dimension(nnode),intent(in)  ::  timol   !
       integer,intent(in)                   ::  nnode   !  Number of molecules
       integer,intent(in)                   ::  rsize   !  Maximum aggregate size
       integer,intent(in)                   ::  tsize   !  Maximum aggregate size
       logical,dimension(nnode),intent(in)  ::  rlife   !
       logical,dimension(nnode),intent(in)  ::  rdeath  !
       logical,dimension(nnode),intent(in)  ::  tlife   !
       logical,dimension(nnode),intent(in)  ::  tdeath  !
!
! Local variables
! 
       integer,dimension(nnode)             ::  qmol    !  Molecules identifier
       integer                              ::  iqueue  !  Empty position in the queue index
       integer                              ::  nqueue  !
       integer                              ::  iiagg   !  Target aggregate index
       integer                              ::  jiagg   !  Reference aggregate index
       integer                              ::  iimol   !  First target molecule index
       integer                              ::  jimol   !  First reference molecule index
       integer                              ::  isize   !  Size of the target aggregate
       integer                              ::  jsize   !  Size of the reference aggregate
       integer                              ::  qsize   !  Size of the aggregate in the queue
       integer                              ::  nsize   !
       integer                              ::  iisize  !  Target molecule index
       integer                              ::  jisize  !  Reference molecule index
       integer                              ::  qisize  !  
       integer                              ::  ni,nj   !  Molecules indexes
       integer                              ::  nq      !  Molecules in queue index
       logical                              ::  next    !  
!
! Breaking previous clusters in smaller motifs
!
!~ !       do isize = tsize, 3, -1
       do isize = 4, 3, -1
         if ( tnagg(isize) .ne. 0 ) then
!
write(*,'(1X,A,1X,I3)') 'Breaking aggregates of size',isize
write(*,'(1X,A)')       '**************************************'
write(*,*)
!
           iimol = timol(isize)
           do iiagg = 1, tnagg(isize)  
! 
!~ !             ni = iimol + 1
!
             if ( tdeath(tiagg(isize)+iiagg) ) then
! 
write(*,'(3X,A,X,I4,1X,A,1X,I6)') 'Starting with aggregate',iiagg,':',tiagg(isize)+iiagg
write(*,'(3X,A)')              '-------------------------------------------------'
write(*,*)
!
! Starts comparison from aggregates of lower size
               qsize  = isize
               iqueue = 1
               next   = .TRUE.
! Adding current aggregate to the queue
               do iisize = 1, qsize
                 qmol(iisize) = tmol(iimol+iisize)
               end do
!              
               do while ( next )
!
write(*,'(A,20I5)') '- Actual aggregate in the queue',qmol(:qsize)
write(*,*)
!
                 call breakagg(nnode,qsize,qmol,rsize,rmol,ragg,       &
                               rnagg,riagg,rnmol,rimol,                &
                               rlife,rdeath,next)
!
               end do
!
             end if
!
             iimol = iimol + isize
!
           end do                !  iiagg
!
         end if 
       end do                    !  isize
!
       return
       end subroutine scrnalrd
!
!======================================================================!
!
! BREAKAGG - BREAK AGGregates
!
! This subroutine 
!
!
       subroutine breakagg(nnode,qsize,qmol,rsize,rmol,ragg,rnagg,     &
                           riagg,rnmol,rimol,rlife,rdeath,next)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)     ::  rmol    !  Molecules identifier
       integer,dimension(nnode),intent(inout)  ::  qmol    !  Molecules identifier
       integer,dimension(nnode),intent(in)     ::  ragg    !  Aggregates size
       integer,dimension(nnode),intent(in)     ::  rnagg   !  
       integer,dimension(nnode),intent(in)     ::  riagg   !
       integer,dimension(nnode),intent(in)     ::  rnmol   !  
       integer,dimension(nnode),intent(in)     ::  rimol   !
       integer,intent(inout)                   ::  qsize   !
       integer,intent(in)                      ::  nnode   !  Number of molecules
       integer,intent(in)                      ::  rsize   !  Maximum aggregate size
       logical,dimension(nnode),intent(in)     ::  rlife   !
       logical,dimension(nnode),intent(in)     ::  rdeath  !
       logical,intent(out)                     ::  next    !                              
!
! Local variables
! 
       integer                                 ::  iqueue  !  Empty position in the queue index
       integer                                 ::  nqueue  !
       integer                                 ::  iiagg   !  Target aggregate index
       integer                                 ::  jiagg   !  Reference aggregate index
       integer                                 ::  iimol   !  First target molecule index
       integer                                 ::  jimol   !  First reference molecule index
       integer                                 ::  isize   !  Size of the target aggregate
       integer                                 ::  jsize   !  Size of the reference aggregate
       integer                                 ::  nsize   !
       integer                                 ::  iisize  !  Target molecule index
       integer                                 ::  jisize  !  Reference molecule index
       integer                                 ::  qisize  !  
       integer                                 ::  ni,nj   !  Molecules indexes
       integer                                 ::  nq      !  Molecules in queue index
!
! 
!
       next = .FALSE.
!
       nsize = qsize
!
       do jsize = qsize-1, 2, -1
         if ( rnagg(jsize) .ne. 0 ) then
!
           jimol = rimol(jsize)
           do jiagg = 1, rnagg(jsize)
!
!~ !                       nj    = jimol + 1
!
             if ( rdeath(riagg(jsize)+jiagg) ) then
! 
write(*,'(3X,A,X,I4,A,1X,I6)') '-> Comparing with aggregate',jiagg,':',riagg(jsize)+jiagg
write(*,'(8X,A,1X,20I5)') 'Queue     :',qmol(:nsize)
write(*,'(8X,A,1X,20I5)') 'Aggregate :',rmol(jimol+1:jimol+jsize)
write(*,*)
!
!                         
! Checking if the first element of the reference array is in the queue
!   If false then look for partial matching
!   If true  then move to the next element in the queue
!
! Loop over jisize elements
!   Loop over iisize elements
!   -> Break if there are not enough elements
!
               qisize = 1
               iisize = 1
!
               do while ( (iisize.le.nsize) .and. (qisize.le.jsize) )
!
! If we do not have enough elements in the current reference aggregate 
!  to do the comparison then exit the loop and go to the next aggregate
!
write(*,'(2(A,1X,I2,1X),A,I6)') 'Taking element',iisize,'of',nsize,'from the queue :',qmol(iisize)
!                         
                 do jisize = qisize, jsize
!
                   nj = jimol + jisize
!
write(*,'(2(A,1X,I2,1X),A,1X,I6,1X)') '  Comparing with molecule',jisize,'of',jsize,':',rmol(nj)
!~ !write(*,'(2(2(A,1X,I2,1X),A,1X,I6,1X))')                        &
!~ !'Comparing TAR',iisize,'of',qsize,':',qmol(iisize),'with REF',jisize,'of',jsize,':',rmol(jimol+jisize)
!
                   if ( qmol(iisize) .lt. rmol(nj) ) then
!
                     qisize = jisize
                     exit
!
                   else if ( qmol(iisize) .eq. rmol(nj) ) then
!
write(*,*) '  WE HAVE A MATCH'
!
!~                      next   = .TRUE.
                     qsize  = qsize - 1
                     qisize = jisize + 1
                     exit
!                               
                   end if
!
!~ !                             if ( tmol(ni) .eq. rmol(nj) ) then
!~ ! Updating the system information
!~ !                               iagg(nnmol)   = nagg
!~ !
!~ !                               imol(nnmol)   = jnode
!~ !                               nnmol         = nnmol + 1
!~ !
!~ !                               ntag          = ntag + 1
!~ ! Marking the node connected to node k as visited
!~ !                               notvis(jnode) = .FALSE.
!~ ! Adding to the queue the node connected to node k
!~ !                               queue(nqueue) = jnode
!~ ! Updating next position in the queue
!~ !                               nqueue        = nqueue + 1
!~ !                             end if
                 end do  !  iisize
!
!~                  nsize = nsize - 1  ! FLAG: check when we have to update nsize
!
write(*,*)
!
                 iisize = iisize + 1
!
               end do
!
!~                if ( next ) return
!
             end if
!
             jimol = jimol + jsize
!
           end do      !  jiagg
!
         end if
       end do          !  jsize
!
       return
       end subroutine breakagg
!
!======================================================================!
!
       end module screening
!
!======================================================================!
