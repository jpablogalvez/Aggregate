!======================================================================!
!
       module graphtools
       implicit none
!
       contains
!
!======================================================================!
!
! BUILDADJMOLLINK - BUILD ADJacency matrix MOLecule-based representation
!                   using a LINKed list
!
       subroutine buildadjmollink(nnode,adj,msubg,mgrps,thr,ngrps,     &
                                  igrps,grps,natms,posi,list,nhead,    &
                                  head,cell,ncell,box)
!
       use geometry,   only: sminimgvec
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(out)      ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(in)      ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)            ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=8),dimension(mgrps,mgrps),intent(in)  ::  thr     !
       integer,dimension(nhead*13),intent(in)          ::  cell    !
       integer,dimension(3),intent(in)                 ::  ncell   !
       integer,dimension(natms),intent(in)             ::  list    !
       integer,dimension(nhead),intent(in)             ::  head    !
       integer,dimension(mgrps),intent(in)             ::  ngrps   !  
       integer,dimension(mgrps),intent(in)             ::  igrps   ! 
       integer,dimension(msubg),intent(in)             ::  grps    ! 
       integer,intent(in)                              ::  msubg   !
       integer,intent(in)                              ::  mgrps   !  Number of subgroups
       integer,intent(in)                              ::  nnode   !  Number of residues
       integer,intent(in)                              ::  natms   ! 
       integer,intent(in)                              ::  nhead   !
!
! Local variables
!
       real(kind=4),dimension(3)                       ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                    ::  dis      !  Minimum image distance
       real(kind=4)                                    ::  mindis   !
       integer                                         ::  icell    !
       integer                                         ::  jcell    !
       integer                                         ::  iicell   !
       integer                                         ::  incell   !
       integer                                         ::  iinode   !
       integer                                         ::  innode   !
       integer                                         ::  jinode   !
       integer                                         ::  iigrps   !
       integer                                         ::  ingrps   !
       integer                                         ::  jnsubg   !
       integer                                         ::  ni,nj    !
       integer                                         ::  i        !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
       adj(:,:) = .FALSE. 
!
       do iinode = 1, nnode
         innode = (iinode-1)*msubg
         do ingrps = 1, mgrps           ! alternatively insubg = 1, msubg
           i = innode + igrps(ingrps)   ! alternatively ni = innode + insubg
           do iigrps = 1, ngrps(ingrps)
!
             ni = i + iigrps
!
             icell  = 1 + int(posi(1,ni)/box(1)*ncell(1))              &
                     + int(posi(2,ni)/box(2)*ncell(2))*ncell(1)        &
                     + int(posi(3,ni)/box(3)*ncell(3))*ncell(1)*ncell(2) 
!
! Loop over all atoms below atom NI in the current cell
!
             nj = list(ni)
             do while ( nj .gt. 0 )
               jinode = (nj - 1)/msubg + 1    ! FLAG: chek if this is correct
               if ( (jinode.ne.iinode) .and.                           &
                                        (.NOT.adj(jinode,iinode)) ) then
                 jnsubg = nj - (jinode - 1)*msubg
                 mindis = real(thr(grps(jnsubg),ingrps)) ! alternatively thr(grps(jngrps),grps(ingrps))
                 if ( mindis .gt. 1.0e-5 ) then
                   r   = sminimgvec(posi(:,ni),posi(:,nj),box)
                   dis = dot_product(r,r)
                   if ( dis .le. mindis ) then
                     adj(jinode,iinode) = .TRUE.
                     adj(iinode,jinode) = .TRUE.
                   end if
                 end if
               end if
               nj = list(nj)
             end do
!
! Loop over all atoms in the neighboring cells
!
             incell = 13*(icell - 1)
!
             do iicell = 1, 13
!
               jcell = cell(incell+iicell)
!
               nj = head(jcell)
               do while ( nj .gt. 0 )
                 jinode = (nj - 1)/msubg + 1    ! FLAG: chek if this is correct
                 if ( (jinode.ne.iinode) .and.                         &
                                        (.NOT.adj(jinode,iinode)) ) then
                   jnsubg = nj - (jinode - 1)*msubg
                   mindis = real(thr(grps(jnsubg),ingrps)) ! alternatively thr(grps(jngrps),grps(ingrps))
                   if ( mindis .gt. 1.0e-5 ) then
                     r   = sminimgvec(posi(:,ni),posi(:,nj),box)
                     dis = dot_product(r,r)
                     if ( dis .le. mindis ) then
                       adj(jinode,iinode) = .TRUE.
                       adj(iinode,jinode) = .TRUE.
                     end if
                   end if
                 end if
                 nj = list(nj)
               end do
             end do
!
           end do
         end do
       end do
!
       return
       end subroutine buildadjmollink
!
!======================================================================!
!
! BUILDADJMOLBUB - BUILD ADJacency matrix MOLecule-based representation
!                   using BUBles
!
!~        subroutine buildadjmolbub(nnode,adj,neidis,msubg,mgrps,thr,     &
!~                                  ngrps,igrps,grps,natms,posi,box)
!~ !
!~        use geometry,   only: sminimgvec
!~ !
!~        implicit none
!~ !
!~ ! Input/output variables
!~ !
!~        logical,dimension(nnode,nnode),intent(out)      ::  adj     !  Adjacency matrix
!~        real(kind=4),dimension(3,natms),intent(in)      ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
!~        real(kind=4),dimension(3),intent(in)            ::  box     !  Simulation box !FLAG: kind=8 to kind=4
!~        real(kind=8),dimension(mgrps,mgrps),intent(in)  ::  thr     !
!~        real(kind=8),intent(in)                         ::  neidis  !
!~        integer,dimension(mgrps),intent(in)             ::  ngrps   !  
!~        integer,dimension(mgrps),intent(in)             ::  igrps   ! 
!~        integer,dimension(msubg),intent(in)             ::  grps    ! 
!~        integer,intent(in)                              ::  msubg   !
!~        integer,intent(in)                              ::  mgrps   !  Number of subgroups
!~        integer,intent(in)                              ::  nnode   !  Number of residues
!~        integer,intent(in)                              ::  natms   ! 
!~ !
!~ ! Local variables
!~ !
!~        real(kind=4),dimension(3)                       ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
!~        real(kind=8)                                    ::  dis      !  Minimum image distance
!~        real(kind=8)                                    ::  mindis   !  Distance threshold between groups
!~        integer                                         ::  iinode   !
!~        integer                                         ::  innode   !
!~        integer                                         ::  jinode   !
!~        integer                                         ::  jnnode   !
!~        integer                                         ::  iigrps   !
!~        integer                                         ::  ingrps   !
!~        integer                                         ::  jigrps   !
!~        integer                                         ::  jngrps   !
!~        integer                                         ::  ni,nj    !
!~        integer                                         ::  i,j      !
!~ !
!~ ! Building the adjacency matrix in the molecule-based representation
!~ ! ------------------------------------------------------------------
!~ !
!~        adj(:,:) = .FALSE. 
!~ !
!~        do iinode = 1, nnode-1
!~          innode = (iinode-1)*msubg
!~          do jinode = iinode+1, nnode
!~            jnnode = (jinode-1)*msubg
!~ !
!~            ingrps = 0
!~            do while ( (ingrps.lt.mgrps) .or. (dis.gt.mindis) .or.      &
!~                                                        (dis.le.neidis) )
!~              ingrps = ingrps + 1
!~              i = innode + igrps(ingrps) 
!~ ! 
!~              iigrps = 0
!~              do while ( (iigrps.lt.ngrps(ingrps)) .or. (dis.gt.mindis) &
!~                                                   .or. (dis.le.neidis) )
!~                iigrps = iigrps + 1
!~                ni = i + iigrps
!~ !
!~                jngrps = 0
!~                do while ( (jngrps.lt.mgrps) .or. (dis.gt.mindis)       &
!~                                                   .or. (dis.le.neidis) )
!~                  jngrps = jngrps + 1
!~ !
!~                  mindis = thr(ingrps,jngrps)   
!~ !
!~                  if ( mindis .gt. 1.0e-6 ) then 
!~                    j = jnnode + igrps(jngrps)
!~ !
!~                    jigrps = 0
!~                    do while ( (jigrps.lt.ngrps(jngrps)) .or.           &
!~                                   (dis.gt.mindis) .or. (dis.le.neidis) )
!~                      jigrps = jigrps + 1
!~ !
!~                      nj = j + jigrps                    
!~ !
!~                      r   = sminimgvec(posi(:,ni),posi(:,nj),box)
!~                      dis = dot_product(r,r)
!~ !
!~                      if ( dis .le. mindis ) then
!~                        adj(iinode,jinode) = .TRUE.
!~                        adj(jinode,iinode) = .TRUE.
!~                      end if
!~ !
!~                    end do  !  jigrps
!~                  end if
!~               end do       !  jngrps
!~              end do        !  iigrps
!~            end do          !  ingrps
!~          end do            !  jinode
!~        end do              !  iinode
!~ !
!~        return
!~        end subroutine buildadjmolbub
!
!======================================================================!
!
! BUILDADJMOLBUB - BUILD ADJacency matrix in the MOLecule-based 
!                   representation using BUBles
!
       subroutine buildadjmolbub(nnode,adj,neidis,msubg,mgrps,nat,     &
                                 thr,ngrps,igrps,natms,posi,box)
!
       use geometry,   only: sminimgvec
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(out)  ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(in)  ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=8),dimension(nat,nat),intent(in)  ::  thr     !
       real(kind=8),intent(in)                     ::  neidis  !
       integer,dimension(nat),intent(in)           ::  ngrps   !  
       integer,dimension(nat),intent(in)           ::  igrps   ! 
       integer,intent(in)                          ::  msubg   !
       integer,intent(in)                          ::  mgrps   !  Number of subgroups
       integer,intent(in)                          ::  nnode   !  Number of residues
       integer,intent(in)                          ::  natms   ! 
       integer,intent(in)                          ::  nat     ! 
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=8)                                ::  dis      !  Minimum image distance
       real(kind=8)                                ::  mindis   !  Distance threshold between groups
       integer                                     ::  iinode   !
       integer                                     ::  innode   !
       integer                                     ::  jinode   !
       integer                                     ::  jnnode   !
       integer                                     ::  iigrps   !
       integer                                     ::  ingrps   !
       integer                                     ::  jigrps   !
       integer                                     ::  jngrps   !
       integer                                     ::  ni       !
       integer                                     ::  i,j      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
       adj(:,:) = .FALSE. 
!
       do iinode = 1, nnode-1
         innode = (iinode-1)*msubg
         do jinode = iinode+1, nnode
           jnnode = (jinode-1)*msubg
           do ingrps = 1, mgrps           
             i = innode + igrps(ingrps)  
             do iigrps = 1, ngrps(ingrps)
               ni = i + iigrps
               do jngrps = 1, mgrps
!
                 mindis = thr(jngrps,ingrps)   
!
                 if ( mindis .gt. 1.0e-6 ) then 
                   j = jnnode + igrps(jngrps)
!
                   do jigrps = 1, ngrps(jngrps)
!                     nj = j + jigrps                    
!
                     r   = sminimgvec(posi(:,ni),posi(:,j+jigrps),box)
                     dis = dot_product(r,r)
!
                     if ( dis .le. mindis ) then
                       adj(iinode,jinode) = .TRUE.
                       adj(jinode,iinode) = .TRUE.
                       GO TO 1000
                     end if
!
                     if ( dis .gt. neidis ) then
                       GO TO 1000
                     end if
!
                   end do  !  jigrps
                 end if
!
              end do       !  jngrps
             end do        !  iigrps
           end do          !  ingrps
1000       continue           
         end do            !  jinode
       end do              !  iinode
!
       return
       end subroutine buildadjmolbub
!
!======================================================================!
!
! BUILDADJMOL - BUILD ADJacency matrix MOLecule-based representation
!
       subroutine buildadjmol(ngrps,grps,subg,nnode,igrps,adj,thr,     &
                              maxdis,natconf,coord,natmol,box)
!
       use geometry
!
       implicit none
!
       include 'idxadj.h'
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(out)        ::  adj      !  Adjacency matrix
       real(kind=4),dimension(3,natconf),intent(inout)   ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=8),dimension(natmol,natmol),intent(in)  ::  thr      !  Distance threshold
       real(kind=4),dimension(3),intent(in)              ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=8),intent(in)                           ::  maxdis   !  Screening distance
       integer,dimension(natmol),intent(in)              ::  grps     !  Number of subgroups in each group
       integer,dimension(natmol),intent(in)              ::  subg     !  Number of atoms in each subgroup
       integer,dimension(natmol),intent(in)              ::  igrps    !  Atoms identifier
       integer,intent(in)                                ::  ngrps    !  Number of groups
       integer,intent(in)                                ::  nnode    !  Number of residues
       integer,intent(in)                                ::  natconf  !  Total number of atoms
       integer,intent(in)                                ::  natmol   !  Atoms per residue
!
! Local variables
!
       real(kind=4),dimension(3)                         ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=8)                                      ::  dist     !  Minimum image distance
       real(kind=8)                                      ::  mindis   !  Distance threshold between groups
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
       adj(:,:) = .FALSE. 
!
       do irenum = 1, nnode-1
         ingrps = (irenum-1)*natmol
         do jrenum = irenum+1, nnode
           jngrps = (jrenum-1)*natmol
!
           iisubg = 0
           ii     = 0
           do iigrps = 1, ngrps
             do insubg = 1, grps(iigrps)
               iisubg = iisubg + 1
               ii     = ii + subg(iisubg)
               i      = ingrps+igrps(ii)
!
               jisubg = 0
               jj     = 0
               do jigrps = 1, ngrps
                 mindis = thr(iigrps,jigrps)   
!
                 if ( mindis .gt. 1.0e-6 ) then 
                   do jnsubg = 1, grps(jigrps)
                     jisubg = jisubg + 1
                     jj     = jj + subg(jisubg)
                     j      = jngrps + igrps(jj)
!
                     r    = sminimgvec(coord(:,i),coord(:,j),box)
                     dist = dot_product(r,r)
!
                     if ( dist .le. mindis ) then
                       adj(jrenum,irenum) = .TRUE.
                       adj(irenum,jrenum) = .TRUE.
                       GO TO 1000
                     end if
!
                     if ( dist .gt. maxdis ) then
                       GO TO 1000
                     end if
                   end do
                 else
                   jisubg = jisubg + grps(jigrps)
                   jj     = jj + subg(jisubg)
                 end if
!
              end do
             end do
           end do
1000       continue           
         end do
       end do
!
       return
       end subroutine buildadjmol
!
!======================================================================!
!
! BUILDADJBODY - BUILD ADJacency matrix n-BODY simplified representation
!
!~        subroutine buildadjbody(natmol,nbody,body,grps,     &
!~                                igrps,adj,thr,nnode,imol,adjmol,   &
!~                                natconf,coord,box)
!~ !
!~        use geometry,   only: sminimgvec
!~ !
!~        implicit none
!~ !
!~        include 'idxadj.h'
!~ !
!~ ! Input/output variables
!~ !
!~        logical,dimension(nbody*nnode,nbody*nnode),intent(out) ::  adj      !  Adjacency matrix
!~        logical,dimension(nnode,nnode),intent(in)              ::  adjmol   !  Adjacency matrix
!~        real(kind=4),dimension(3,natconf),intent(in)           ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
!~        real(kind=8),dimension(natmol,natmol),intent(in)       ::  thr      !  Distance threshold
!~        real(kind=4),dimension(3),intent(in)                   ::  box      !  Simulation box !FLAG: kind=8 to kind=4
!~        integer,dimension(nnode),intent(in)                    ::  imol     !   
!~        integer,dimension(natmol),intent(in)                   ::  body     !  Number of groups in each body
!~        integer,dimension(natmol),intent(in)                   ::  grps     !  Number of subgroups in each group
!~        integer,dimension(natmol),intent(in)                   ::  igrps    !  Atoms identifier
!~        integer,intent(in)                                     ::  nbody    !  Number of bodies
!~        integer,intent(in)                                     ::  nnode    !  Number of residues
!~        integer,intent(in)                                     ::  natconf  !  Total number of atoms
!~        integer,intent(in)                                     ::  natmol   !  Atoms per residue
!~ !
!~ ! Local variables
!~ !
!~        real(kind=4),dimension(3)                              ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
!~        real(kind=8)                                           ::  dist     !  Minimum image distance
!~        real(kind=8)                                           ::  mindis   !  Distance threshold between groups
!~ !
!~ ! Building the adjacency matrix in the N-body simplified representation
!~ ! ---------------------------------------------------------------------
!~ !
!~        adj(:,:) = .FALSE. 
!~ !
!~        inbody = 0
!~        do irenum = 1, nnode-1
!write(*,'(1X,A,X,I4)') 'Starting IRESIDUE',irenum
!write(*,*) '----------------------------------'
!write(*,*)
!~          inat = (imol(irenum)-1)*natmol
!~ !
!~          jnbody = inbody + nbody           
!~          do jrenum = irenum+1, nnode
!~            if ( adjmol(jrenum,irenum) ) then
!write(*,'(1X,A,X,I4)') 'Comparing with JRESIDUE',jrenum
!write(*,*)
!~              jnat = (imol(jrenum)-1)*natmol
!~ !
!~              ingrps = 0
!~              ii     = 0
!~              do iibody = 1, nbody
!~                jngrps = 0
!~                jj     = 0
!~                do jibody = 1, nbody
!write(*,'(3X,6(X,A,X,I3),X,A,X,I6)') 'COMPARING inbody',inbody+iibody, &
!                                          'WITH jibody',jnbody+jibody
!write(*,*)
!~                  insubg = ii
!~                  do iigrps = 1, body(iibody)
!~                    do iisubg = 1, grps(ingrps+iigrps)
!~                      iat = inat + igrps(insubg+iisubg)    ! FLAG: now it takes the first atom in the group
!~ !
!write(*,'(5X,2(X,A,X,I3),X,A,X,I6,X,A,X,F4.2)') 'COMPARING iigrps',ingrps+iigrps,  &
!'iisubg',insubg+iisubg,'iat',iat
!~ !
!~                      jnsubg = jj
!~                      do jigrps = 1, body(jibody)
!~ !
!~                        mindis = thr(ingrps+iigrps,jngrps+jigrps)   
!~ !
!~                        if ( mindis .gt. 1.0e-6 ) then 
!~                          do jisubg = 1, grps(jngrps+jigrps)
!~                            r    = sminimgvec(coord(:,iat),coord(:,jnat+&
!~                                               igrps(jnsubg+jisubg)),box)
!~                            dist = dot_product(r,r)
!~ !
!write(*,'(5X,2(X,A,X,I3),X,A,X,I6,2(X,A,X,F4.2))') '     WITH jigrps',jngrps+jigrps,  &
!'jisubg',jnsubg+jisubg,'jat',jnat+igrps(jnsubg+jisubg),'MINDIS',sqrt(mindis),'DIST',sqrt(dist)
!~                            if ( dist .le. mindis ) then
!~                              adj(jnbody+jibody,inbody+iibody) = .TRUE.
!~                              adj(inbody+iibody,jnbody+jibody) = .TRUE.
!~ !
!~                              do j = jigrps, body(jibody)
!~                                jnsubg = jnsubg + grps(jngrps+j)
!~                              end do
!~ !
!~                              do i = iigrps, body(iibody)
!~                                insubg = insubg + grps(ingrps+i)
!~                              end do
!~ !
!~                              GO TO 1000
!~                            end if
!~                          end do
!~                        end if
!~                        jnsubg = jnsubg + grps(jngrps+jigrps)
!~                      end do
!write(*,*)
!~                    end do
!~                    insubg = insubg + grps(ingrps+iigrps)
!~                  end do
!~ 1000             continue
!~                  jj = jnsubg
!~                  jngrps = jngrps + body(jibody)
!~                end do
!~                ii = insubg 
!~                ingrps = ingrps + body(iibody)
!~              end do     
!~            end if 
!~            jnbody = jnbody + nbody 
!~          end do
!~            inbody = inbody + nbody 
!~        end do
!~ !
!~        return
!~        end subroutine buildadjbody
!
!======================================================================!
!
! FINDCOMPUNDIR - FIND COMPonents UNDIRected
!
! This subroutine finds the connected components in an undirected 
!  unweighted graph of NNODE vertices given as an adjacency matrix 
!  representaion ADJ(NODE,NODE) using Breadth First Search.  
!
       subroutine findcompundir(nnode,adj,imol,iagg,itag,maxagg,nmol,  &
                                nagg)
!
       implicit none
!
       include 'idxbfs.h'
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,dimension(nnode),intent(out)       ::  imol    !  Molecules identifier
       integer,dimension(nnode),intent(out)       ::  iagg    !  Aggregates identifier
       integer,dimension(nnode),intent(out)       ::  itag    !  Aggregates size
       integer,dimension(nnode),intent(out)       ::  nmol    !  Number of aggregates of each size
       integer,intent(out)                        ::  nagg    !  Number of aggregates
       integer,intent(in)                         ::  nnode   !  Number of molecules
       integer,intent(out)                        ::  maxagg  !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(nnode)                   ::  notvis  !  Nodes visited
       integer,dimension(nnode)                   ::  queue   !  Queue of connected nodes
       integer                                    ::  i       !  Indexes
!
! Performing Breadth First Search over the target molecules
! ---------------------------------------------------------
!
! Marking all the vertices as not visited
       notvis(:) = .TRUE.
! Initializing the molecules information
       imol(:) = 0   
       nmol(:) = 0
       nnmol   = 1
! Initializing the aggregates information
       nagg    = 0
       iagg(:) = 0
! Initializing the size information
       ntag    = 0
       nntag   = 0
       itag(:) = 0
! 
       maxagg  = 1
!
! Outer loop over each node
!
       do inode = 1, nnode
         if ( notvis(inode) ) then
! Marking head node as visited
           notvis(inode) = .FALSE.
! Updating the system information
           nagg        = nagg + 1
           iagg(nnmol) = nagg
!
           imol(nnmol) = inode
           nnmol       = nnmol + 1
!
           ntag        = 1
! Initializing queue
           queue(:) = 0
! Adding current node to the queue
           queue(1) = inode
! Initializing the queue counter
           iqueue = 1
! Setting the next position in the queue
           nqueue = 2
!
! Inner loop over the queue elements
!
           do while ( iqueue .lt. nqueue )
! Saving actual element in the queue
             knode = queue(iqueue)
! Checking the connection between actual queue element and the rest of nodes
             do jnode = inode + 1, nnode
! Checking if node j is connected to node k and has not been already visited
               if ( adj(jnode,knode) .and. notvis(jnode) ) then
! Updating the system information
                 iagg(nnmol)   = nagg
!
                 imol(nnmol)   = jnode
                 nnmol         = nnmol + 1
!
                 ntag          = ntag + 1
! Marking the node connected to node k as visited
                 notvis(jnode) = .FALSE.
! Adding to the queue the node connected to node k
                 queue(nqueue) = jnode
! Updating next position in the queue
                 nqueue        = nqueue + 1
               end if
             end do
! Updating the queue counter
             iqueue = iqueue + 1
           end do
! Saving the size of the aggregate found
           do i = nntag+1, nntag+ntag
             itag(i) = ntag
           end do
           nntag = nntag + ntag
! Update the number of aggregates of each size
           if ( ntag .gt. maxagg ) maxagg = ntag
           nmol(ntag)   = nmol(ntag)   + 1
         end if
       end do
!
       return
       end subroutine findcompundir
!
!======================================================================!
!
! CALCDEGUNDIR - CALCulate DEGrees UNDIRected
!
! This function calculates the degree of the NNODE vertices in an undi-
!  rected unweighted graph given as an adjacency matrix representaion 
!  ADJ(NODE,NODE) and stores them in the array DEGREE(NNODE)
!
       function calcdegundir(nnode,adj) result(degree)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !
       integer,intent(in)                         ::  nnode   !
       integer,dimension(nnode)                   ::  degree  !
!
! Local variables
!
       integer                                    ::  inode   !
       integer                                    ::  jnode   !
!
! Finding the degree of the vertices in a graph
! ---------------------------------------------
!
       do inode = 1, nnode
         degree(inode) = 0
! Traversing through row/column of each vertex 
         do jnode = 1, nnode
! If a path from this vertex to other exists then increment the degree
           if ( adj(jnode,inode) ) degree(inode) = degree(inode) + 1
         end do
       end do
!    
       return
       end function calcdegundir
!
!======================================================================!
!
! CHKTREE - CHecK TREE
!
! This function identifies wether an undirected unweighted graph of
!   NNODE vertices is a tree given the DEGREE(NNODE) of each vertex
!
! Seen on: https://www.geeksforgeeks.org/check-whether-given-degrees-vertices-represent-graph-tree/
!
       logical function chktree(nnode,degree)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)  ::  degree  !  Degree of each vertex
       integer,intent(in)                   ::  nnode   !  Number of vertices
!
! Local variables
!
       integer                              ::  degsum  !  Sum of all degrees
       integer                              ::  i       !
!
! Checking if the input graph is a tree 
! -------------------------------------
!
       chktree = .FALSE.
! Finding the sum of all degrees
       degsum = 0
       do i = 1, nnode
         degsum = degsum + degree(i)
       end do
! Graph is tree if the sum of all the degrees is equal to 2(n-1) (Handshaking Lemma)
       if ( degsum .eq. 2*(nnode-1) ) chktree = .TRUE.
!    
       return
       end function chktree
!
!======================================================================!
!
! CHKLTREE - CHecK Linear TREE
!
! This function identifies wether an undirected unweighted graph of 
!  NNODE vertices is a linear tree based on the DEGREE(NNODE) of each 
!  of its vertices
!
! Seen on: https://www.geeksforgeeks.org/check-if-a-given-tree-graph-is-linear-or-not/
!
       logical function chkltree(nnode,degree)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)  ::  degree  !  Degree of each vertex
       integer,intent(in)                   ::  nnode   !  Number of vertices
!
! Local variables
!
       integer                              ::  numdeg  !  Number of nodes with degree two
       integer                              ::  i       !  Indexes
!
! Checking if the input graph is a linear tree 
! --------------------------------------------
!
       chkltree = .FALSE.
! Counting the number of vertices with degree 2
       numdeg = 0
       do i = 1, nnode
         if ( degree(i) .gt. 2 ) return
         if ( degree(i) .eq. 2 ) numdeg = numdeg + 1
       end do
! The given tree would be linear only if NNODE-2 of its nodes have 
!  DEGREE=2 or the number of nodes is 1 
       if ( numdeg .eq. nnode-2 ) chkltree = .TRUE.
!    
       return
       end function chkltree
!
!======================================================================!
!
! FINDLONGT - FIND LONGest Tree
!
! This function finds the length of the longest path present in an undi-
!  rected unweighted tree of NNODE vertices given as an adjacency 
!  matrix representaion ADJ(NODE,NODE)
!
! Seen on:  https://www.geeksforgeeks.org/longest-path-undirected-tree/
!
       integer function findlongt(nnode,adj) result(maxdis)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,intent(in)                         ::  nnode   !  Number of vertices
!
! Local variables
!
       real(kind=4)                               ::  rndm    !  Random number
       integer                                    ::  src     !  Source node
       integer                                    ::  idx     !  Index of the farthest node
!
! Finding the length of the longest path of the input tree
! --------------------------------------------------------
!
       call random_number(rndm)
! Choice a random source node (not best performance, but avoids worst-case)
       src = 1 + floor(nnode*rndm)
! First BFS to find an endpoint of the longest path
       call findlongd(nnode,adj,src,idx,maxdis)
! second BFS from this endpoint to find the actual longest path
       call findlongd(nnode,adj,idx,src,maxdis)
!    
       return
       end function findlongt
!
!======================================================================!
!
! FINDLONGP - FIND LONGest Distance
!
! This function identifies the index IDX of an endpoint node with the 
!  longest distance MAXDIS from the source node SRC present in an undi-
!  rected unweighted tree of NNODE vertices given as an adjacency matrix
!  representaion ADJ(NODE,NODE) and returns both values
!
! Seen on: https://www.geeksforgeeks.org/longest-path-undirected-tree/
!
       subroutine findlongd(nnode,adj,src,idx,maxdis)
!
       implicit none
!
       include 'idxbfs.h'
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,intent(in)                         ::  nnode   !  Number of vertices
       integer,intent(in)                         ::  src     !  Source node
       integer,intent(out)                        ::  idx     !  Index of the farthest node
       integer,intent(out)                        ::  maxdis  !  Farthest node distance
!
! Local variables
!
       logical,dimension(nnode)                   ::  notvis  !  Nodes visited
       integer,dimension(nnode)                   ::  queue   !  Queue of connected nodes
       integer,dimension(nnode)                   ::  dist    !  Distances from the source node
!
! Finding farthest node and its distance from source node
! -------------------------------------------------------
!
       maxdis = 0
! Marking all distance with -1
       dist(:) = -1
! Distance of source node src from source node src will be 0
       dist(src) = 0
! Marking all the vertices as not visited
       notvis(:) = .TRUE.
! Marking source node as visited
       notvis(src) = .FALSE.
! Initializing queue
       queue(:) = 0
! Adding source node to the queue
       queue(1) = src
! Initializing the queue counter
       iqueue = 1
! Setting the next position in the queue
       nqueue = 2
! Continue until queue is not empty
       do while ( iqueue .lt. nqueue )
! Taking the next element in the queue
         knode = queue(iqueue)
! Checking the connection between actual queue element and the rest 
!  of nodes
         do jnode = 1, nnode
! Checking if node jnode is connected to node knode and has not been 
!  already visited
           if ( adj(jnode,knode) .and. notvis(jnode) ) then
! Marking the node connected to node knode as visited
             notvis(jnode) = .FALSE.
! Making distance of jnode one more than distance of knode
             dist(jnode) = dist(knode) + 1
!            
             if ( dist(jnode) .gt. maxdis ) then
               maxdis = dist(jnode)
               idx    = jnode
             end if
! Pushing node jnode into the queue only if it is not visited already
             queue(nqueue) = jnode
! Updating next position in the queue
             nqueue        = nqueue + 1
           end if
         end do
! Updating the queue counter
         iqueue = iqueue + 1
       end do
!    
       return
       end subroutine findlongd
!
!======================================================================!
!
! CHKSCYCLE - CHecK Simple CYCLE
!
! This function identifies wether an undirected unweighted graph of 
!  NNODE vertices is a simple cycle based on the DEGREE(NNODE) of each
!  of its vertices
!
       logical function chkscycle(nnode,degree)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nnode),intent(in)  ::  degree  !  Degree of each vertex
       integer,intent(in)                   ::  nnode   !  Number of vertices
!
! Local variables
!
       integer                              ::  i       !  Indexes
!
! Checking if the input graph is a single cycle 
! ---------------------------------------------
!
       chkscycle = .FALSE.
! Checking if any vertex in the graph has degree different from two
       do i = 1, nnode
         if ( degree(i) .ne. 2 ) return
       end do
! The given graph would be a simple cycle only if all the vertices 
!  are having DEGREE=2
       chkscycle = .TRUE.
!    
       return
       end function chkscycle
!
!======================================================================!
!
! FINDSHORTC - FIND SHORTest Cycle
!
! This function finds the length of the shortest cycle in an undirected
!  unweighted graph of NNODE vertices given as an adjacency matrix re-
!  presentaion ADJ(NODE,NODE). Screening of the leaf nodes based on the
!  DEGREE(NNODE) of its vertices is done such that they are not visited
!  during the Breadth First Search. If no cycle exists returns -1 
!
! Seen on: https://www.geeksforgeeks.org/shortest-cycle-in-an-undirected-unweighted-graph/
!
       integer function findshortc(nnode,adj,degree) result(ans)
!
       implicit none
!
       include 'idxbfs.h'
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,dimension(nnode),intent(in)        ::  degree  !  Degree of each vertex
       integer,intent(in)                         ::  nnode   !  Number of vertices
!
! Local variables
!
       logical,dimension(nnode)                   ::  vis     !  Nodes with degree greater than one
       integer,dimension(nnode)                   ::  queue   !  Queue of connected nodes
       integer,dimension(nnode)                   ::  dist    !  Distances from the source node
       integer,dimension(nnode)                   ::  parent  !  Parent nodes
!
! Parameters
!
       integer,parameter                          ::  nmax = int(1e9)  !  Maximum length
!
! Finding the length of the shortest cycle in the given graph
! -----------------------------------------------------------
!
       ans  = nmax
! Do not visit vertices with degree lower than two   ! FLAG: check also nodes of degree 2
       vis = .TRUE.
       do inode = 1, nnode
         if ( degree(inode) .lt. 2 ) vis(inode) = .FALSE.
       end do
! For every vertex, we check if it is possible to get the shortest cycle
!  involving this vertex
       do inode = 1, nnode
! Making distance maximum
         dist(:) = nmax
! Taking an imaginary parent
         parent(:) = -1
! Distance of source to source is 0
         dist(inode) = 0
! If the vertex has degree greater than or equal 2 start BFS
         if ( vis(inode) ) then
! Initializing queue
           queue(:) = 0
! Pushing source node to the queue
           queue(1) = inode
! Initializing the queue counter
           iqueue   = 1
! Setting the next position in the queue
           nqueue   = 2
! Continue until queue is not empty
           do while ( iqueue .lt. nqueue )
! Taking the next element in the queue
             knode = queue(iqueue)
! Checking the connection between actual queue element and the rest
!  of nodes
             do jnode = 1, nnode
! Checking if node j is connected to node k and has degree greater
!  than one
               if ( adj(jnode,knode) .and. vis(jnode) ) then
! If node j has not been already visited then update the information
                 if ( dist(jnode) .eq. nmax ) then
! Increasing distance by 1
                   dist(jnode) = 1 + dist(knode)
! Changing parent node
                   parent(jnode) = knode                               
! Pushing to the queue the node connected to node k
                   queue(nqueue) = jnode
! Updating next position in the queue
                   nqueue        = nqueue + 1
! If the vertex which is already visited comes again and jnode is not a 
!  parent of knode and viceversa then the cycle is present
                  else if ( (parent(knode).ne.jnode) .and.             &
                                         (parent(jnode).ne.knode) ) then
! Checking the length of the cycle 
                    ans = min(ans, dist(knode) + dist(jnode) + 1);
                 end if
               end if
             end do
! Updating the queue counter
             iqueue = iqueue + 1
           end do
         end if
       end do
!    
       return
       end function findshortc
!
!======================================================================!
!
       end module graphtools
!
!======================================================================!
