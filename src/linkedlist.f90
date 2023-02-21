!======================================================================!
!
       module linkedlist
       implicit none
!
       contains
!
!======================================================================!
!
! SETLINKLIST - SET LINKed LIST
!
!~        subroutine setlinklist(natms,posi,drnei,nnode,list,nhead,head,  &
!~                               msubg,mgrps,ngrps,igrps,ncell,box)
!~ !
!~        use geometry, only ::  sminimgvec,                              &
!~                               scenvec
!~ !
!~        implicit none
!~ !
!~ ! Input/output variables
!~ !
!~        real(kind=4),dimension(3,natms),intent(in)     ::  posi    !
!~        real(kind=8),dimension(3,natms),intent(inout)  ::  drnei   !
!~        real(kind=4),dimension(3),intent(in)           ::  box     !
!~        integer,dimension(natms),intent(out)           ::  list    !
!~        integer,dimension(nhead),intent(out)           ::  head    !
!~        integer,dimension(nat),intent(in)              ::  ngrps   !
!~        integer,dimension(nat),intent(in)              ::  igrps   !
!~        integer,dimension(3)                           ::  ncell   !
!~        integer,intent(in)                             ::  msubg   !
!~        integer,intent(in)                             ::  mgrps   !
!~        integer,intent(in)                             ::  nnode   !
!~        integer,intent(in)                             ::  natms   !
!~        integer,intent(in)                             ::  nhead   ! 
!~ !
!~ ! Local variables
!~ !
!~        real(kind=8),dimension(3,msubg)                ::  smaux   !
!~        real(kind=8),dimension(3)                      ::  svaux   !
!~        integer                                        ::  icell   !
!~        integer                                        ::  iinode  !
!~        integer                                        ::  innode  !
!~        integer                                        ::  iigrps  !
!~        integer                                        ::  ingrps  !
!~        integer                                        ::  i,j     !
!~ !
!~ ! Updating the neighbour list
!~ ! ---------------------------
!~ !      
!~        head(:) = 0
!~ !
!~        do iinode = 1, nnode
!~ !
!~          innode = (iinode-1)*msubg
!~          svaux(:) = posi(:,innode+1)
!~ !
!~          do ingrps = 1, mgrps
!~            i = innode + igrps(ingrps)
!~            do iigrps = 1, ngrps(ingrps)
!~              j = igrps(ingrps)+iigrps
!~              smaux(:,j) = sminimgvec(svaux(:),posi(:,i+iigrps)),box)
!~              smaux(:,j) = svaux(:) + smaux(:,j)
!~            end do
!~          end do
!~ !
!~          svaux(:) = scenvec(3,msubg,smaux)
!~ !
!~          icell = 1 + int(svaux(1)/box(1)*ncell(1))                     &
!~                    + int(svaux(2)/box(2)*ncell(2))*ncell(1)            &
!~                    + int(svaux(3)/box(3)*ncell(3))*ncell(1)*ncell(2)
!~ !
!~          list(iinode) = head(icell)
!~          head(icell)  = iinode
!~ !
!~        end do
!~ !
!~        drnei(:,:) = 0.0d0  
!~ !
!~        return
!~        end subroutine setlinklist
!
!======================================================================!
!
! SETLINKLIST - SET LINKed LIST
!
       subroutine setlinklist(natms,posi,drnei,list,nhead,head,        &
                              ncell,box,step)
!
       implicit none
!
! Input/output variables
!
       real(kind=4),dimension(3,natms),intent(inout)  ::  posi    !
       real(kind=4),dimension(3,natms),intent(inout)  ::  drnei   !
       real(kind=4),dimension(3),intent(in)           ::  box     !
       integer,dimension(natms),intent(out)           ::  list    !
       integer,dimension(nhead),intent(out)           ::  head    !
       integer,dimension(3)                           ::  ncell   !
       integer,intent(in)                             ::  natms   !
       integer,intent(in)                             ::  nhead   ! 
       integer,intent(in)                             ::  step
!
! Local variables
!
       integer                                        ::  icell   !
       integer                                        ::  i       !
!
! Setting up the cellular linked list
! -----------------------------------
!      
       head(:) = 0
!
       do i = 1, natms
         icell = 1 + int(posi(1,i)/box(1)*ncell(1))                    &
                   + int(posi(2,i)/box(2)*ncell(2))*ncell(1)           &
                   + int(posi(3,i)/box(3)*ncell(3))*ncell(1)*ncell(2)
         list(i)     = head(icell)
         head(icell) = i
       end do
!
!~        drnei(:,:) = 0.0d0  
!
       return
       end subroutine setlinklist
!
!======================================================================!
!
! UPLINKLIST - UPdate LINKed LIST
!
       subroutine uplinklist(natms,posi,drnei,list,nhead,head,         &
                             ncell,maxdis,box,STEP)
!
       implicit none
!
! Input/output variables
!
       real(kind=4),dimension(3,natms),intent(inout)  ::  posi    !
       real(kind=4),dimension(3,natms),intent(inout)  ::  drnei   !
       real(kind=4),dimension(3),intent(in)           ::  box     !
       real(kind=8),intent(in)                        ::  maxdis  !
       integer,dimension(natms),intent(out)           ::  list    !
       integer,dimension(nhead),intent(out)           ::  head    !
       integer,dimension(3),intent(in)                ::  ncell   !
       integer,intent(in)                             ::  natms   !
       integer,intent(in)                             ::  nhead   !
       integer,intent(in)                :: step
!
! Local variables
!
       real(kind=4)                                   ::  dr      !
       real(kind=4)                                   ::  drmax1  !
       real(kind=4)                                   ::  drmax2  !
       integer                                        ::  i       !
!  
! Checking if the neighbour list has to be updated
! ------------------------------------------------
!
       drmax1 = 0.0d0
       drmax2 = 0.0d0
!
       do i = 1, natms
         dr = dot_product(drnei(:,i),drnei(:,i))
         if ( dr .gt. drmax1 ) then
           drmax2 = drmax1
           drmax1 = dr
         else
            if ( dr .gt. drmax2 ) drmax2 = dr
         end if
       end do
!
! Updating the neighbour list if needed
! -------------------------------------
!
       if ( (sqrt(drmax1)+sqrt(drmax2)) .gt. maxdis ) then
!~ write(*,*) 'UPDATING NEIGHBOUR LIST',step
         call setlinklist(natms,posi,drnei,list,nhead,head,ncell,box,step)   
       end if
!
       return
       end subroutine uplinklist
!
!======================================================================!
!
! GETDISP - GET DISPlacements
!
       subroutine getdisp(natms,drnei,posi,oldposi,box)
!
       use geometry,   only: sminimgvec
!
       implicit none
!
! Input/output variables
!
       real(kind=4),dimension(3,natms),intent(in)     ::  oldposi  !
       real(kind=4),dimension(3,natms),intent(in)     ::  posi     !
       real(kind=8),dimension(3,natms),intent(inout)  ::  drnei    !
       real(kind=4),dimension(3),intent(in)           ::  box      !
       integer,intent(in)                             ::  natms    !
!
! Local variables
!
       real(kind=4),dimension(3)                      ::  r        !
       integer                                        ::  i,j      !
!  
! Compute displacements between two snapshots
! -------------------------------------------
!
       do i = 1, natms
         r(:) = sminimgvec(oldposi(:,i),posi(:,i),box)
         do j = 1, 3
           drnei(j,i) = drnei(j,i) + r(j)
         end do
       end do
!
       return
       end subroutine getdisp
!
!======================================================================!
!
! SETNEIG - SET NEIGhbors
!
       subroutine setneig(nhead,cell,ncell)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nhead*13),intent(out)  ::  cell   !
       integer,dimension(3),intent(in)          ::  ncell  !
       integer,intent(in)                       ::  nhead  !
!
! Local variables
!
       integer                                  ::  i,j,k  !
       integer                                  ::  n      !
!  
! Storing half the nearest neighbours of each cell
! ------------------------------------------------
!
       do k = 1, ncell(3)
         do j = 1, ncell(2)
           do i = 1, ncell(1)
!
             n = (idxcell(i,j,k,ncell) - 1)*13
!
             cell(n+1)  = idxcell(i+1,j  ,k  ,ncell)
             cell(n+2)  = idxcell(i+1,j+1,k  ,ncell)
             cell(n+3)  = idxcell(i  ,j+1,k  ,ncell)
             cell(n+4)  = idxcell(i-1,j+1,k  ,ncell)
             cell(n+5)  = idxcell(i+1,j  ,k-1,ncell)
             cell(n+6)  = idxcell(i+1,j+1,k-1,ncell)
             cell(n+7)  = idxcell(i  ,j+1,k-1,ncell)
             cell(n+8)  = idxcell(i-1,j+1,k-1,ncell)
             cell(n+9)  = idxcell(i+1,j  ,k+1,ncell)
             cell(n+10) = idxcell(i+1,j+1,k+1,ncell)
             cell(n+11) = idxcell(i  ,j+1,k+1,ncell)
             cell(n+12) = idxcell(i-1,j+1,k+1,ncell)
             cell(n+13) = idxcell(i  ,j  ,k+1,ncell)  
!               
           end do
         end do
       end do
!
       return
       end subroutine setneig
!
!======================================================================!
!
! IDXCELL - InDeX CELL
!
       integer function idxcell(ix,iy,iz,ncell)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(3),intent(in)  ::  ncell  !
       integer,intent(in)               ::  ix     !
       integer,intent(in)               ::  iy     !
       integer,intent(in)               ::  iz     !
!   
       idxcell = 1 + mod(ix-1+ncell(1),ncell(1))                       &
                   + mod(iy-1+ncell(2),ncell(2))*ncell(1)              &
                   + mod(iz-1+ncell(3),ncell(3))*ncell(1)*ncell(2)
!
       return
       end function idxcell
!
!======================================================================!
!
       end module linkedlist
!
!======================================================================!
