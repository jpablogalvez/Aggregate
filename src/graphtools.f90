!======================================================================!
!
       module graphtools
!
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
       use geometry,  only:  sminimgvec
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(out)      ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(in)      ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)            ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=4),dimension(mgrps,mgrps),intent(in)  ::  thr     !
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
! BUILDADJMOLBUB - BUILD ADJacency matrix
!                   in the MOLecule-based representation
!                   using BUBles
!
       subroutine buildadjmolbub(nnode,adj,neidis,msubg,mgrps,nat,     &
                                 ngrps,igrps,natms,posi,box)
!
       use omp_lib
!
       use thresholds, only:  thr
!
       use geometry,   only:  sminimgvec
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(out)  ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(in)  ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=4),intent(in)                     ::  neidis  !
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
       real(kind=4)                                ::  dis      !  Minimum image distance
       real(kind=4)                                ::  mindis   !  Distance threshold between groups
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
!$omp parallel do num_threads(np)                                      &
!$omp              shared(adj)                                         &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkadj)
!
       do i = 1, nnode
           adj(:,i) = .FALSE.
       end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,posi,box,thr,neidis,igrps,ngrps)          &
!$omp             private(r,dis,mindis,iinode,innode,jinode,jnnode,    &
!$omp                     iigrps,ingrps,jigrps,jngrps,ni,i,j)          &
!$omp             schedule(dynamic,chunkadj)
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
!$omp end parallel do
!
       return
       end subroutine buildadjmolbub
!
!======================================================================!
!
! NBUILDADJMOLBUB - N-components BUILD ADJacency matrix
!                    in the MOLecule-based representation
!                    using BUBles
!
       subroutine nbuildadjmolbub(mnode,adj,matms,posi,neidis,box)
!
       use omp_lib
!
       use systeminf,  only:  rep,mtype,nnode,inode,iatms,igrps
       use thresholds, only:  thr
!
       use geometry,   only:  sminimgvec
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(mnode,mnode),intent(out)  ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=4),intent(in)                     ::  neidis  !
       integer,intent(in)                          ::  mnode   !
       integer,intent(in)                          ::  matms   !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis      !  Minimum image distance
       real(kind=4)                                ::  mindis   !  Distance threshold between groups
       integer                                     ::  iinode   !
       integer                                     ::  innode   !
       integer                                     ::  jinode   !
       integer                                     ::  jnnode   !
       integer                                     ::  iigrps   !
       integer                                     ::  ingrps   !
       integer                                     ::  jigrps   !
       integer                                     ::  jngrps   !
       integer                                     ::  iiatms   !
       integer                                     ::  jiatms   !
       integer                                     ::  ithr     !
       integer                                     ::  jthr     !
       integer                                     ::  iadj     !
       integer                                     ::  jadj     !
       integer                                     ::  ni       !
       integer                                     ::  ii,jj    !
       integer                                     ::  i,j      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj)                                          &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkadj)
!
       do i = 1, mnode
         adj(:,i) = .FALSE.
       end do
!
!$omp end parallel do
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
         iiatms = iatms(ii)
         ithr   = igrps(ii)
         iadj   = inode(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,posi,box,thr,neidis,rep,nnode,ii,iiatms,  &
!$omp                    ithr,iadj)                                    &
!$omp             private(r,dis,mindis,iinode,innode,jinode,jnnode,    &
!$omp                     iigrps,ingrps,jigrps,jngrps,ni,i,j)          &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
           innode = iiatms + (iinode-1)*rep(ii)%msubg
           do jinode = iinode+1, nnode(ii)
             jnnode = iiatms + (jinode-1)*rep(ii)%msubg
!
             do ingrps = 1, rep(ii)%mgrps
               i = innode + rep(ii)%igrps(ingrps)
               do iigrps = 1, rep(ii)%ngrps(ingrps)
                 ni = i + iigrps
!
                 do jngrps = 1, rep(ii)%mgrps
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. 1.0e-6 ) then
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       r   = sminimgvec(posi(:,ni),posi(:,j+jigrps),box)
                       dis = dot_product(r,r)
!
                       if ( dis .le. mindis ) then
                         adj(iadj+iinode,iadj+jinode) = .TRUE.
                         adj(iadj+jinode,iadj+iinode) = .TRUE.
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
!
               end do        !  iigrps
             end do          !  ingrps
!
1000         continue
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
      end do                 !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
         iiatms = iatms(ii)
         ithr   = igrps(ii)
         iadj   = inode(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,posi,box,thr,neidis,rep,nnode,ii,         &
!$omp                    iiatms,ithr,iadj)                             &
!$omp             private(r,dis,mindis,iinode,innode,jinode,jnnode,    &
!$omp                     iigrps,ingrps,jigrps,jngrps,jiatms,jthr,     &
!$omp                     jadj,ni,jj,i,j)                              &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
           innode = iiatms + (iinode-1)*rep(ii)%msubg

           do jj = ii+1, mtype
             jiatms = iatms(jj)
             jthr   = igrps(jj)
             jadj   = inode(jj)
             do jinode = 1, nnode(jj)
               jnnode = jiatms + (jinode-1)*rep(jj)%msubg
!
               do ingrps = 1, rep(ii)%mgrps
                 i = innode + rep(ii)%igrps(ingrps)
                 do iigrps = 1, rep(ii)%ngrps(ingrps)
                   ni = i + iigrps
!
                   do jngrps = 1, rep(jj)%mgrps
!
                     mindis = thr(jthr+jngrps,ithr+ingrps)
!
                     if ( mindis .gt. 1.0e-6 ) then
                       j = jnnode + rep(jj)%igrps(jngrps)
!
                       do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                         r = sminimgvec(posi(:,ni),posi(:,j+jigrps),box)
                         dis = dot_product(r,r)
!
                         if ( dis .le. mindis ) then
                           adj(iadj+iinode,jadj+jinode) = .TRUE.
                           adj(jadj+jinode,iadj+iinode) = .TRUE.
                           GO TO 2000
                         end if
!
                         if ( dis .gt. neidis ) then
                           GO TO 2000
                         end if
!
                       end do  !  jigrps
                     end if
!
                  end do       !  jngrps
!
                 end do        !  iigrps
               end do          !  ingrps
!
2000           continue
!
             end do            !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
       return
       end subroutine nbuildadjmolbub
!
!======================================================================!
!
! NBUILDADJBODYBUB - N-components BUILD ADJacency matrix
!                     in the n-BODY simplified representation
!                     using bubles
!
       subroutine nbuildadjbodybub(mnode,node,madj,adj,matms,posi,     &
                                   box,mtype,nnode,inode,ibodymon)
!
       use omp_lib
!
       use systeminf,  only:  rep,iat,iatms,igrps
!
       use thresholds, only:  thr
       use parameters, only:  zero
!
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(madj,madj),intent(inout)  ::  adj       !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi      !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box       !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(mnode),intent(in)         ::  node      !
       integer,dimension(mtype),intent(in)         ::  nnode     !
       integer,dimension(mtype),intent(in)         ::  inode     !
       integer,dimension(mtype),intent(in)         ::  ibodymon  !
       integer,intent(in)                          ::  madj      !
       integer,intent(in)                          ::  mnode     !
       integer,intent(in)                          ::  mtype     !
       integer,intent(in)                          ::  matms     !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1       !  Minimum image distance
       real(kind=4)                                ::  dis2       !  Minimum image distance
       real(kind=4)                                ::  angle      !
       real(kind=4)                                ::  mindis     !  Distance threshold between groups
       real(kind=4)                                ::  minang     !  Angle threshold between groups
       integer                                     ::  iinode     !
       integer                                     ::  innode     !
       integer                                     ::  jinode     !
       integer                                     ::  jnnode     !
       integer                                     ::  ibody      !
       integer                                     ::  iibody     !
       integer                                     ::  inbody     !
       integer                                     ::  jbody      !
       integer                                     ::  jibody     !
       integer                                     ::  jnbody     !
       integer                                     ::  iigrps     !
       integer                                     ::  ingrps     !
       integer                                     ::  jigrps     !
       integer                                     ::  jngrps     !
       integer                                     ::  iiat       !
       integer                                     ::  jiat       !
       integer                                     ::  iiatms     !
       integer                                     ::  jiatms     !
       integer                                     ::  ithr       !
       integer                                     ::  jthr       !
       integer                                     ::  iadj       !
       integer                                     ::  jadj       !
       integer                                     ::  i,j        !
       integer                                     ::  ii,jj      !
       integer                                     ::  ni,nj      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,rep,nnode,inode,        &
!$omp                    ii,iiatms,iiat,ithr,ibody)                    &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inbody,iibody,jnbody,jibody,   &
!$omp                     iadj,jadj,ni,nj,i,j)                         &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jinode = iinode+1, nnode(ii)
!
             jnnode = iiatms + (node(inode(ii)+jinode)-1)*rep(ii)%msubg
             jadj   = ibody  + (jinode-1)*rep(ii)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
!
                 do jnbody = 1, rep(ii)%mbody
                 do jibody = 1, rep(ii)%nbody(jnbody)
!
                   jngrps = rep(ii)%ibody(jnbody) + jibody
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       nj   = j + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
                         adj(iadj+inbody,jadj+jnbody) = .TRUE.
                         adj(jadj+jnbody,iadj+inbody) = .TRUE.
                         GO TO 1000
                       end if
!
                     end do  !  jigrps
                   end if
!
                 end do      !  jibody
!
1000             continue
!
                 end do      !  jnbody
!
               end do        !  iigrps
             end do          !  iibody
             end do          !  inbody
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
       end do                !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,rep,nnode,inode,iatms,  &
!$omp                    iat,igrps,ibodymon,ii,iiatms,iiat,ithr,ibody) &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inbody,iibody,jnbody,jibody,   &
!$omp                     iadj,jadj,ni,nj,i,j,jj,jiatms,jiat,          &
!$omp                     jthr,jbody)                                  &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jj = ii+1, mtype
!
             jiatms = iatms(jj)
             jiat   = iat(jj)
             jthr   = igrps(jj)
             jbody  = ibodymon(jj)
!
           do jinode = 1, nnode(jj)
!
             jnnode = jiatms + (node(inode(jj)+jinode)-1)*rep(jj)%msubg
             jadj   = jbody  + (jinode-1)*rep(jj)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
!
                 do jnbody = 1, rep(jj)%mbody
                 do jibody = 1, rep(jj)%nbody(jnbody)
!
                   jngrps = rep(jj)%ibody(jnbody) + jibody
!
                   mindis = thr(jthr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     j = jnnode + rep(jj)%igrps(jngrps)
!
                     do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                       nj   = j + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
                         adj(iadj+inbody,jadj+jnbody) = .TRUE.
                         adj(jadj+jnbody,iadj+inbody) = .TRUE.
                         GO TO 2000
                       end if
!
                     end do    !  jigrps
                   end if
!
                 end do        !  jibody
!
2000             continue
!
                 end do        !  jnbody
!
               end do          !  iigrps
             end do            !  iibody
             end do            !  inbody
!
           end do              !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
!
       return
       end subroutine nbuildadjbodybub
!
!======================================================================!
!
! NBUILDADJBODYBUBDIR - N-components BUILD ADJacency matrix
!                        in the n-BODY simplified representation
!                        using bubles
!                        for a DIRected graph
!
       subroutine nbuildadjbodybubdir(mnode,node,madj,adj,matms,posi,  &
                                      box,mtype,nnode,inode,ibodymon)
!
       use omp_lib
!
       use systeminf,  only:  rep,iat,iatms,igrps
!
       use thresholds, only:  thr
       use parameters, only:  zero
!
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(madj,madj),intent(inout)  ::  adj       !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi      !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box       !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(mnode),intent(in)         ::  node      !
       integer,dimension(mtype),intent(in)         ::  nnode     !
       integer,dimension(mtype),intent(in)         ::  inode     !
       integer,dimension(mtype),intent(in)         ::  ibodymon  !
       integer,intent(in)                          ::  madj      !
       integer,intent(in)                          ::  mnode     !
       integer,intent(in)                          ::  mtype     !
       integer,intent(in)                          ::  matms     !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1       !  Minimum image distance
       real(kind=4)                                ::  dis2       !  Minimum image distance
       real(kind=4)                                ::  mindis     !  Distance threshold between groups
       integer                                     ::  iinode     !
       integer                                     ::  innode     !
       integer                                     ::  jinode     !
       integer                                     ::  jnnode     !
       integer                                     ::  ibody      !
       integer                                     ::  iibody     !
       integer                                     ::  inbody     !
       integer                                     ::  jbody      !
       integer                                     ::  jibody     !
       integer                                     ::  jnbody     !
       integer                                     ::  iigrps     !
       integer                                     ::  ingrps     !
       integer                                     ::  jigrps     !
       integer                                     ::  jngrps     !
       integer                                     ::  inei       !
       integer                                     ::  jnei       !
       integer                                     ::  innei      !
       integer                                     ::  jnnei      !
       integer                                     ::  iiat       !
       integer                                     ::  jiat       !
       integer                                     ::  iiatms     !
       integer                                     ::  jiatms     !
       integer                                     ::  ithr       !
       integer                                     ::  jthr       !
       integer                                     ::  iadj       !
       integer                                     ::  jadj       !
       integer                                     ::  i,j        !
       integer                                     ::  ii,jj      !
       integer                                     ::  ni,nj      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,rep,nnode,              &
!$omp                    inode,ii,iiatms,iiat,ithr,ibody)              &
!$omp             private(v21,v23,dis1,dis2,mindis,iinode,innode,      &
!$omp                     jinode,jnnode,iigrps,ingrps,jigrps,jngrps,   &
!$omp                     inei,jnei,innei,jnnei,inbody,iibody,jnbody,  &
!$omp                     jibody,iadj,jadj,ni,nj,i,j)                  &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jinode = iinode+1, nnode(ii)
!
             jnnode = iiatms + (node(inode(ii)+jinode)-1)*rep(ii)%msubg
             jnnei  = iiat   + (node(inode(ii)+jinode)-1)*rep(ii)%nat
             jadj   = ibody  + (jinode-1)*rep(ii)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jnbody = 1, rep(ii)%mbody
                 do jibody = 1, rep(ii)%nbody(jnbody)
!
                   jngrps = rep(ii)%ibody(jnbody) + jibody
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(ii)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
                         if ( rep(ii)%neiang(jnei) .ne. 0 ) then ! points towards j
!
                           adj(iadj+inbody,jadj+jnbody) = .TRUE.
                           GO TO 1000
!
                         else if ( rep(ii)%neiang(inei) .ne. 0 ) then ! points towards i
!
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 1000
!
                         else
!
                           adj(iadj+inbody,jadj+jnbody) = .TRUE.
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 1000
!
                         end if
                       end if
!
                     end do  !  jigrps
                   end if
!
                 end do      !  jibody
!
1000             continue
!
                 end do      !  jnbody
!
               end do        !  iigrps
             end do          !  iibody
             end do          !  inbody
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
       end do                !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,rep,nnode,inode,iatms,  &
!$omp                    iat,igrps,ibodymon,ii,iiatms,iiat,ithr,ibody) &
!$omp             private(v21,v23,dis1,dis2,mindis,iinode,innode,      &
!$omp                     jinode,jnnode,iigrps,ingrps,jigrps,jngrps,   &
!$omp                     inei,jnei,innei,jnnei,inbody,iibody,jnbody,  &
!$omp                     jibody,iadj,jadj,ni,nj,i,j,jj,jiatms,jiat,   &
!$omp                     jthr,jbody)                                  &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jj = ii+1, mtype
!
             jiatms = iatms(jj)
             jiat   = iat(jj)
             jthr   = igrps(jj)
             jbody  = ibodymon(jj)
!
           do jinode = 1, nnode(jj)
!
             jnnode = jiatms + (node(inode(jj)+jinode)-1)*rep(jj)%msubg
             jnnei  = jiat   + (node(inode(jj)+jinode)-1)*rep(jj)%nat
             jadj   = jbody  + (jinode-1)*rep(jj)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jnbody = 1, rep(jj)%mbody
                 do jibody = 1, rep(jj)%nbody(jnbody)
!
                   jngrps = rep(jj)%ibody(jnbody) + jibody
!
                   mindis = thr(jthr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     j = jnnode + rep(jj)%igrps(jngrps)
!
                     do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(jj)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
                         if ( rep(jj)%neiang(jnei) .ne. 0 ) then ! points towards j
!
                           adj(iadj+inbody,jadj+jnbody) = .TRUE. 
                           GO TO 2000
!
                         else if ( rep(ii)%neiang(inei) .ne. 0 ) then ! points towards i
!
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 2000
!
                         else
!
                           adj(iadj+inbody,jadj+jnbody) = .TRUE.
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 2000
!
                         end if
                       end if
!
                     end do    !  jigrps
                   end if
!
                 end do        !  jibody
!
2000             continue
!
                 end do        !  jnbody
!
               end do          !  iigrps
             end do            !  iibody
             end do            !  inbody
!
           end do              !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
!
       return
       end subroutine nbuildadjbodybubdir
!
!======================================================================!
!
! NBUILDADJGRPSBUB - N-components BUILD ADJacency matrix
!                     in the GRouPS representation
!                     using bubles
!
       subroutine nbuildadjgrpsbub(mnode,node,madj,adj,matms,posi,     &
                                   box,mtype,nnode,inode,igrpsmon)
!
       use omp_lib
!
       use systeminf,  only:  rep,iat,iatms,igrps
!
       use thresholds, only:  thr
       use parameters, only:  zero
!
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(madj,madj),intent(inout)  ::  adj       !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi      !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box       !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(mnode),intent(in)         ::  node      !
       integer,dimension(mtype),intent(in)         ::  nnode     !
       integer,dimension(mtype),intent(in)         ::  inode     !
       integer,dimension(mtype),intent(in)         ::  igrpsmon  !
       integer,intent(in)                          ::  madj      !
       integer,intent(in)                          ::  mnode     !
       integer,intent(in)                          ::  mtype     !
       integer,intent(in)                          ::  matms     !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1       !  Minimum image distance
       real(kind=4)                                ::  dis2       !  Minimum image distance
       real(kind=4)                                ::  angle      !
       real(kind=4)                                ::  mindis     !  Distance threshold between groups
       real(kind=4)                                ::  minang     !  Angle threshold between groups
       integer                                     ::  iinode     !
       integer                                     ::  innode     !
       integer                                     ::  jinode     !
       integer                                     ::  jnnode     !
       integer                                     ::  iigrps     !
       integer                                     ::  ingrps     !
       integer                                     ::  jigrps     !
       integer                                     ::  jngrps     !
       integer                                     ::  iiat       !
       integer                                     ::  jiat       !
       integer                                     ::  iiatms     !
       integer                                     ::  jiatms     !
       integer                                     ::  ithr       !
       integer                                     ::  jthr       !
       integer                                     ::  iiadj      !
       integer                                     ::  jiadj      !
       integer                                     ::  iadj       !
       integer                                     ::  jadj       !
       integer                                     ::  i,j        !
       integer                                     ::  ii,jj      !
       integer                                     ::  ni,nj      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         iiadj  = igrpsmon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,rep,nnode,inode,        &
!$omp                    ii,iiatms,iiat,ithr,iiadj)                    &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,iadj,jadj,ni,nj,i,j)           &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           iadj   = iiadj  + (iinode-1)*rep(ii)%mgrps
!
           do jinode = iinode+1, nnode(ii)
!
             jnnode = iiatms + (node(inode(ii)+jinode)-1)*rep(ii)%msubg
             jadj   = iiadj  + (jinode-1)*rep(ii)%mgrps
!
             do ingrps = 1, rep(ii)%mgrps
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
!
                 do jngrps = 1, rep(ii)%mgrps
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       nj   = j + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
                         adj(iadj+ingrps,jadj+jngrps) = .TRUE.
                         adj(jadj+jngrps,iadj+ingrps) = .TRUE.
                         GO TO 1000
                       end if
!
                     end do  !  jigrps
                   end if
!
1000               continue
!
                 end do      !  jngrps
!
               end do        !  iigrps
             end do          !  ingrps
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
       end do                !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         iiadj  = igrpsmon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,rep,nnode,inode,iatms,  &
!$omp                    iat,igrps,igrpsmon,ii,iiatms,iiat,ithr,iiadj) &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,iadj,jadj,ni,nj,i,j,jj,        &
!$omp                     jiatms,jiat,jthr,jiadj)                      &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           iadj   = iiadj  + (iinode-1)*rep(ii)%mgrps
!
           do jj = ii+1, mtype
!
             jiatms = iatms(jj)
             jiat   = iat(jj)
             jthr   = igrps(jj)
             jiadj  = igrpsmon(jj)
!
           do jinode = 1, nnode(jj)
!
             jnnode = jiatms + (node(inode(jj)+jinode)-1)*rep(jj)%msubg
             jadj   = jiadj  + (jinode-1)*rep(jj)%mgrps
!
             do ingrps = 1, rep(ii)%mgrps
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
!
                 do jngrps = 1, rep(jj)%mgrps
!
                   mindis = thr(jthr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     j = jnnode + rep(jj)%igrps(jngrps)
!
                     do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                       nj   = j + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
                         adj(iadj+ingrps,jadj+jngrps) = .TRUE.
                         adj(jadj+jngrps,iadj+ingrps) = .TRUE.
                         GO TO 2000
                       end if
!
                     end do    !  jigrps
                   end if
!
2000               continue
!
                 end do        !  jngrps
!
               end do          !  iigrps
             end do            !  ingrps
!
           end do              !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
!
       return
       end subroutine nbuildadjgrpsbub
!
!======================================================================!
!
! BUILDADJMOLANG - BUILD ADJacency matrix
!                   in the MOLecule-based representation
!                   using bubles and ANGle restraints
!
       subroutine buildadjmolang(nnode,adj,neidis,msubg,mgrps,nat,     &
                                 ngrps,igrps,natms,posi,box)
!
       use omp_lib
!
       use systeminf,  only:  xtcf
       use thresholds, only:  thr,thrang,neiang
!
       use parameters, only:  zero
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(out)  ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(in)  ::  posi    !  Subgroup coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=4),intent(in)                     ::  neidis  !
       integer,dimension(nat),intent(in)           ::  ngrps   !
       integer,dimension(nat),intent(in)           ::  igrps   !
       integer,intent(in)                          ::  msubg   !
       integer,intent(in)                          ::  mgrps   !  Number of subgroups
       integer,intent(in)                          ::  nnode   !  Number of residues
       integer,intent(in)                          ::  natms   !  Number of subgroups in the system
       integer,intent(in)                          ::  nat     !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21      !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23      !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1     !  Minimum image distance
       real(kind=4)                                ::  dis2     !  Minimum image distance
       real(kind=4)                                ::  angle    !  Minimum image distance
       real(kind=4)                                ::  mindis   !  Distance threshold between groups
       real(kind=4)                                ::  minang   !  Angle threshold between groups
       integer                                     ::  iinode   !
       integer                                     ::  innode   !
       integer                                     ::  jinode   !
       integer                                     ::  jnnode   !
       integer                                     ::  iigrps   !
       integer                                     ::  ingrps   !
       integer                                     ::  jigrps   !
       integer                                     ::  jngrps   !
       integer                                     ::  ineiang  !
       integer                                     ::  jneiang  !
       integer                                     ::  innei    !
       integer                                     ::  jnnei    !
       integer                                     ::  ni,nj    !
       integer                                     ::  i,j      !
       logical                                     ::  doang    !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj)                                          &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkadj)
!
       do i = 1, nnode
         adj(:,i) = .FALSE.
       end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,posi,box,thr,thrang,neiang,neidis,igrps,  &
!$omp                    ngrps)                                        &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,ineiang,jneiang,innei,jnnei,   &
!$omp                     doang,ni,nj,i,j)                             &
!$omp             schedule(dynamic,chunkadj)
!
       do iinode = 1, nnode-1
!
         innode = (iinode-1)*msubg
         innei  = (iinode-1)*nat
!
         do jinode = iinode+1, nnode
!
           jnnode = (jinode-1)*msubg
           jnnei  = (jinode-1)*nat
!
           do ingrps = 1, mgrps
!
             i = innode + igrps(ingrps)
!
             do iigrps = 1, ngrps(ingrps)
!
               ni = i + iigrps
               ineiang = igrps(ingrps) + iigrps
!
               do jngrps = 1, mgrps
!
                 mindis = thr(jngrps,ingrps)
!
                 if ( mindis .gt. zero ) then
!
                   minang = thrang(jngrps,ingrps)
!
                   j = jnnode + igrps(jngrps)
!
                   do jigrps = 1, ngrps(jngrps)
!
                     nj = j + jigrps
                     jneiang = igrps(jngrps) + jigrps
!
                     v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                     dis1 = dot_product(v21,v21)
!
                     if ( dis1 .le. mindis ) then
!
                       doang = .TRUE.
!
                       if ( minang .gt. zero ) then
!
                         if ( neiang(jneiang) .ne. 0 ) then
!
                           v21(:) = -v21(:)
!
                           doang = chkangle(v21,dis1,posi(:,nj),       &
                           xtcf%pos(:,jnnei+neiang(jneiang)),box,minang)
!
                         else if ( neiang(ineiang) .ne. 0 ) then
!
                           doang = chkangle(v21,dis1,posi(:,ni),       &
                           xtcf%pos(:,innei+neiang(ineiang)),box,minang)
!
                         end if
!
                       end if
!
                       if ( doang ) then
                         adj(iinode,jinode) = .TRUE.
                         adj(jinode,iinode) = .TRUE.
                         GO TO 1000
                       end if
!
                     end if
!
                     if ( dis1 .gt. neidis ) then
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
!$omp end parallel do
!
       return
       end subroutine buildadjmolang
!
!======================================================================!
!
! NBUILDADJMOLANG - N-components BUILD ADJacency matrix
!                    in the MOLecule-based representation
!                    using bubles and ANGle restraints
!
       subroutine nbuildadjmolang(mnode,adj,matms,posi,neidis,box)
!
       use omp_lib
!
       use systeminf,  only:  xtcf,rep,mtype,nnode,inode,iat,iatms,igrps
!
       use thresholds, only:  thr,thrang
       use parameters, only:  zero
!
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(mnode,mnode),intent(out)  ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=4),intent(in)                     ::  neidis  !
       integer,intent(in)                          ::  mnode   !
       integer,intent(in)                          ::  matms   !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21      !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23      !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1     !  Minimum image distance
       real(kind=4)                                ::  dis2     !  Minimum image distance
       real(kind=4)                                ::  angle    !
       real(kind=4)                                ::  mindis   !  Distance threshold between groups
       real(kind=4)                                ::  minang   !  Angle threshold between groups
       integer                                     ::  iinode   !
       integer                                     ::  innode   !
       integer                                     ::  jinode   !
       integer                                     ::  jnnode   !
       integer                                     ::  iigrps   !
       integer                                     ::  ingrps   !
       integer                                     ::  jigrps   !
       integer                                     ::  jngrps   !
       integer                                     ::  inei     !
       integer                                     ::  jnei     !
       integer                                     ::  innei    !
       integer                                     ::  jnnei    !
       integer                                     ::  iiat     !
       integer                                     ::  jiat     !
       integer                                     ::  iiatms   !
       integer                                     ::  jiatms   !
       integer                                     ::  ithr     !
       integer                                     ::  jthr     !
       integer                                     ::  iadj     !
       integer                                     ::  jadj     !
       integer                                     ::  ii,jj    !
       integer                                     ::  ni,nj    !
       integer                                     ::  i,j      !
       logical                                     ::  doang    !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj)                                          &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkadj)
!
       do i = 1, mnode
         adj(:,i) = .FALSE.
       end do
!
!$omp end parallel do
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         iadj   = inode(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,posi,box,thr,thrang,rep,neidis,nnode,ii,  &
!$omp                    iiat,iiatms,ithr,iadj)                        &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     ni,nj,i,j)                                   &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
!
           innode = iiatms + (iinode-1)*rep(ii)%msubg
           innei  = iiat   + (iinode-1)*rep(ii)%nat
!
           do jinode = iinode+1, nnode(ii)
!
             jnnode = iiatms + (jinode-1)*rep(ii)%msubg
             jnnei  = iiat   + (jinode-1)*rep(ii)%nat
!
             do ingrps = 1, rep(ii)%mgrps
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jngrps = 1, rep(ii)%mgrps
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(ithr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(ii)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(ii)%neiang(jnei) .ne. 0 ) then
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(ii)%neiang(jnei)),  &
                                                             box,minang)
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                           end if
!
                         end if
!
                         if ( doang ) then
                           adj(iadj+iinode,iadj+jinode) = .TRUE.
                           adj(iadj+jinode,iadj+iinode) = .TRUE.
                           GO TO 1000
                         end if
!
                       end if
!
                       if ( dis1 .gt. neidis ) then
                         GO TO 1000
                       end if
!
                     end do  !  jigrps
                   end if
!
                 end do      !  jngrps
               end do        !  iigrps
             end do          !  ingrps
!
1000         continue
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
       end do                !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         iadj   = inode(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,posi,box,thr,thrang,rep,neidis,nnode,ii,  &
!$omp                    iiat,iiatms,ithr,iadj)                        &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     ni,nj,i,j,jj,jiat,jiatms,jthr,jadj)          &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
!
           innode = iiatms + (iinode-1)*rep(ii)%msubg
           innei  = iiat   + (iinode-1)*rep(ii)%nat
!
           do jj = ii+1, mtype
!
             jiatms = iatms(jj)
             jiat   = iat(jj)
             jthr   = igrps(jj)
             jadj   = inode(jj)
!
           do jinode = 1, nnode(jj)
!
             jnnode = jiatms + (jinode-1)*rep(jj)%msubg
             jnnei  = jiat   + (jinode-1)*rep(jj)%nat
!
             do ingrps = 1, rep(ii)%mgrps
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jngrps = 1, rep(jj)%mgrps
!
                   mindis = thr(jthr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(jthr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(jj)%igrps(jngrps)
!
                     do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(jj)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(jj)%neiang(jnei) .ne. 0 ) then
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(jj)%neiang(jnei)),  &
                                                             box,minang)
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                           end if
!
                         end if
!
                         if ( doang ) then
                           adj(iadj+iinode,jadj+jinode) = .TRUE.
                           adj(jadj+jinode,iadj+iinode) = .TRUE.
                           GO TO 2000
                         end if
!
                       end if
!
                       if ( dis1 .gt. neidis ) then
                         GO TO 2000
                       end if
!
                     end do    !  jigrps
                   end if
!
                 end do        !  jngrps
               end do          !  iigrps
             end do            !  ingrps
!
2000         continue
!
           end do              !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
!
       return
       end subroutine nbuildadjmolang
!
!======================================================================!
!
! NBUILDADJBODYANG - N-components BUILD ADJacency matrix
!                     in the n-BODY simplified representation
!                     using bubles and ANGle restraints
!
       subroutine nbuildadjbodyang(mnode,node,madj,adj,matms,posi,     &
                                   box,mtype,nnode,inode,ibodymon)
!
       use omp_lib
!
       use systeminf,  only:  xtcf,rep,iat,iatms,igrps
!
       use thresholds, only:  thr,thrang
       use parameters, only:  zero
!
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(madj,madj),intent(inout)  ::  adj       !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi      !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box       !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(mnode),intent(in)         ::  node      !
       integer,dimension(mtype),intent(in)         ::  nnode     !
       integer,dimension(mtype),intent(in)         ::  inode     !
       integer,dimension(mtype),intent(in)         ::  ibodymon  !
       integer,intent(in)                          ::  madj      !
       integer,intent(in)                          ::  mnode     !
       integer,intent(in)                          ::  mtype     !
       integer,intent(in)                          ::  matms     !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1       !  Minimum image distance
       real(kind=4)                                ::  dis2       !  Minimum image distance
       real(kind=4)                                ::  angle      !
       real(kind=4)                                ::  mindis     !  Distance threshold between groups
       real(kind=4)                                ::  minang     !  Angle threshold between groups
       integer                                     ::  iinode     !
       integer                                     ::  innode     !
       integer                                     ::  jinode     !
       integer                                     ::  jnnode     !
       integer                                     ::  ibody      !
       integer                                     ::  iibody     !
       integer                                     ::  inbody     !
       integer                                     ::  jbody      !
       integer                                     ::  jibody     !
       integer                                     ::  jnbody     !
       integer                                     ::  iigrps     !
       integer                                     ::  ingrps     !
       integer                                     ::  jigrps     !
       integer                                     ::  jngrps     !
       integer                                     ::  inei       !
       integer                                     ::  jnei       !
       integer                                     ::  innei      !
       integer                                     ::  jnnei      !
       integer                                     ::  iiat       !
       integer                                     ::  jiat       !
       integer                                     ::  iiatms     !
       integer                                     ::  jiatms     !
       integer                                     ::  ithr       !
       integer                                     ::  jthr       !
       integer                                     ::  iadj       !
       integer                                     ::  jadj       !
       integer                                     ::  i,j        !
       integer                                     ::  ii,jj      !
       integer                                     ::  ni,nj      !
       logical                                     ::  doang      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,thrang,rep,nnode,       &
!$omp                    inode,ii,iiatms,iiat,ithr,ibody)              &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     inbody,iibody,jnbody,jibody,                 &
!$omp                     iadj,jadj,ni,nj,i,j)                         &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jinode = iinode+1, nnode(ii)
!
             jnnode = iiatms + (node(inode(ii)+jinode)-1)*rep(ii)%msubg
             jnnei  = iiat   + (node(inode(ii)+jinode)-1)*rep(ii)%nat
             jadj   = ibody  + (jinode-1)*rep(ii)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jnbody = 1, rep(ii)%mbody
                 do jibody = 1, rep(ii)%nbody(jnbody)
!
                   jngrps = rep(ii)%ibody(jnbody) + jibody
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(ithr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(ii)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(ii)%neiang(jnei) .ne. 0 ) then
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(ii)%neiang(jnei)),  &
                                                             box,minang)
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                           end if
!
                         end if
!
                         if ( doang ) then
                           adj(iadj+inbody,jadj+jnbody) = .TRUE.
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 1000
                         end if
!
                       end if
!
                     end do  !  jigrps
                   end if
!
                 end do      !  jibody
!
1000             continue
!
                 end do      !  jnbody
!
               end do        !  iigrps
             end do          !  iibody
             end do          !  inbody
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
       end do                !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,thrang,rep,nnode,       &
!$omp                    inode,iatms,iat,igrps,ibodymon,ii,iiatms,     &
!$omp                    iiat,ithr,ibody)                              &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     inbody,iibody,jnbody,jibody,                 &
!$omp                     iadj,jadj,ni,nj,i,j,jj,jiatms,jiat,          &
!$omp                     jthr,jbody)                                  &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jj = ii+1, mtype
!
             jiatms = iatms(jj)
             jiat   = iat(jj)
             jthr   = igrps(jj)
             jbody  = ibodymon(jj)
!
           do jinode = 1, nnode(jj)
!
             jnnode = jiatms + (node(inode(jj)+jinode)-1)*rep(jj)%msubg
             jnnei  = jiat   + (node(inode(jj)+jinode)-1)*rep(jj)%nat
             jadj   = jbody  + (jinode-1)*rep(jj)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jnbody = 1, rep(jj)%mbody
                 do jibody = 1, rep(jj)%nbody(jnbody)
!
                   jngrps = rep(jj)%ibody(jnbody) + jibody
!
                   mindis = thr(jthr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(jthr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(jj)%igrps(jngrps)
!
                     do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(jj)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(jj)%neiang(jnei) .ne. 0 ) then
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(jj)%neiang(jnei)),  &
                                                             box,minang)
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                           end if
!
                         end if
!
                         if ( doang ) then
                           adj(iadj+inbody,jadj+jnbody) = .TRUE.
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 2000
                         end if
!
                       end if
!
                     end do    !  jigrps
                   end if
!
                 end do        !  jibody
!
2000             continue
!
                 end do        !  jnbody
!
               end do          !  iigrps
             end do            !  iibody
             end do            !  inbody
!
           end do              !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
!
       return
       end subroutine nbuildadjbodyang
!
!======================================================================!
!
! NBUILDADJBODYANGDIR - N-components BUILD ADJacency matrix
!                        in the n-BODY simplified representation
!                        using bubles and ANGle restraints
!                        for a DIRected graph
!
       subroutine nbuildadjbodyangdir(mnode,node,madj,adj,matms,posi,  &
                                      box,mtype,nnode,inode,ibodymon)
!
       use omp_lib
!
       use systeminf,  only:  xtcf,rep,iat,iatms,igrps
!
       use thresholds, only:  thr,thrang
       use parameters, only:  zero
!
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(madj,madj),intent(inout)  ::  adj       !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi      !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box       !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(mnode),intent(in)         ::  node      !
       integer,dimension(mtype),intent(in)         ::  nnode     !
       integer,dimension(mtype),intent(in)         ::  inode     !
       integer,dimension(mtype),intent(in)         ::  ibodymon  !
       integer,intent(in)                          ::  madj      !
       integer,intent(in)                          ::  mnode     !
       integer,intent(in)                          ::  mtype     !
       integer,intent(in)                          ::  matms     !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1       !  Minimum image distance
       real(kind=4)                                ::  dis2       !  Minimum image distance
       real(kind=4)                                ::  angle      !
       real(kind=4)                                ::  mindis     !  Distance threshold between groups
       real(kind=4)                                ::  minang     !  Angle threshold between groups
       integer                                     ::  iinode     !
       integer                                     ::  innode     !
       integer                                     ::  jinode     !
       integer                                     ::  jnnode     !
       integer                                     ::  ibody      !
       integer                                     ::  iibody     !
       integer                                     ::  inbody     !
       integer                                     ::  jbody      !
       integer                                     ::  jibody     !
       integer                                     ::  jnbody     !
       integer                                     ::  iigrps     !
       integer                                     ::  ingrps     !
       integer                                     ::  jigrps     !
       integer                                     ::  jngrps     !
       integer                                     ::  inei       !
       integer                                     ::  jnei       !
       integer                                     ::  innei      !
       integer                                     ::  jnnei      !
       integer                                     ::  iiat       !
       integer                                     ::  jiat       !
       integer                                     ::  iiatms     !
       integer                                     ::  jiatms     !
       integer                                     ::  ithr       !
       integer                                     ::  jthr       !
       integer                                     ::  iadj       !
       integer                                     ::  jadj       !
       integer                                     ::  i,j        !
       integer                                     ::  ii,jj      !
       integer                                     ::  ni,nj      !
       logical                                     ::  doang      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,thrang,rep,nnode,       &
!$omp                    inode,ii,iiatms,iiat,ithr,ibody)              &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     inbody,iibody,jnbody,jibody,                 &
!$omp                     iadj,jadj,ni,nj,i,j)                         &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jinode = iinode+1, nnode(ii)
!
             jnnode = iiatms + (node(inode(ii)+jinode)-1)*rep(ii)%msubg
             jnnei  = iiat   + (node(inode(ii)+jinode)-1)*rep(ii)%nat
             jadj   = ibody  + (jinode-1)*rep(ii)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jnbody = 1, rep(ii)%mbody
                 do jibody = 1, rep(ii)%nbody(jnbody)
!
                   jngrps = rep(ii)%ibody(jnbody) + jibody
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(ithr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(ii)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(ii)%neiang(jnei) .ne. 0 ) then ! points towards j
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(ii)%neiang(jnei)),  &
                                                             box,minang)
!
                             adj(iadj+inbody,jadj+jnbody) = .TRUE.
                             GO TO 1000
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then ! points towards i
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                             adj(jadj+jnbody,iadj+inbody) = .TRUE.
                             GO TO 1000
!
                           end if
!
                         else
!
                           adj(iadj+inbody,jadj+jnbody) = .TRUE.
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 1000
!
                         end if
!
                       end if
!
                     end do  !  jigrps
                   end if
!
                 end do      !  jibody
!
1000             continue
!
                 end do      !  jnbody
!
               end do        !  iigrps
             end do          !  iibody
             end do          !  inbody
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
       end do                !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         ibody  = ibodymon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,thrang,rep,nnode,       &
!$omp                    inode,iatms,iat,igrps,ibodymon,ii,iiatms,     &
!$omp                    iiat,ithr,ibody)                              &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     inbody,iibody,jnbody,jibody,                 &
!$omp                     iadj,jadj,ni,nj,i,j,jj,jiatms,jiat,          &
!$omp                     jthr,jbody)                                  &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = ibody  + (iinode-1)*rep(ii)%mbody
!
           do jj = ii+1, mtype
!
             jiatms = iatms(jj)
             jiat   = iat(jj)
             jthr   = igrps(jj)
             jbody  = ibodymon(jj)
!
           do jinode = 1, nnode(jj)
!
             jnnode = jiatms + (node(inode(jj)+jinode)-1)*rep(jj)%msubg
             jnnei  = jiat   + (node(inode(jj)+jinode)-1)*rep(jj)%nat
             jadj   = jbody  + (jinode-1)*rep(jj)%mbody
!
             do inbody = 1, rep(ii)%mbody
             do iibody = 1, rep(ii)%nbody(inbody)
!
               ingrps = rep(ii)%ibody(inbody) + iibody
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jnbody = 1, rep(jj)%mbody
                 do jibody = 1, rep(jj)%nbody(jnbody)
!
                   jngrps = rep(jj)%ibody(jnbody) + jibody
!
                   mindis = thr(jthr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(jthr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(jj)%igrps(jngrps)
!
                     do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(jj)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(jj)%neiang(jnei) .ne. 0 ) then ! points towards j
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(jj)%neiang(jnei)),  &
                                                             box,minang)
!
                             if ( doang ) then
                               adj(iadj+inbody,jadj+jnbody) = .TRUE. 
                               GO TO 2000
                             end if
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then ! points towards i
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                             if ( doang ) then
                               adj(jadj+jnbody,iadj+inbody) = .TRUE.
                               GO TO 2000
                             end if
!
                           end if
!
                         else
!
                           adj(iadj+inbody,jadj+jnbody) = .TRUE.
                           adj(jadj+jnbody,iadj+inbody) = .TRUE.
                           GO TO 2000
!
                         end if
!
                       end if
!
                     end do    !  jigrps
                   end if
!
                 end do        !  jibody
!
2000             continue
!
                 end do        !  jnbody
!
               end do          !  iigrps
             end do            !  iibody
             end do            !  inbody
!
           end do              !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
!
       return
       end subroutine nbuildadjbodyangdir
!
!======================================================================!
!
! NBUILDADJGRPSANG - N-components BUILD ADJacency matrix
!                     in the groups representation
!                     using bubles and ANGle restraints
!
       subroutine nbuildadjgrpsang(mnode,node,madj,adj,matms,posi,     &
                                   box,mtype,nnode,inode,igrpsmon)
!
       use omp_lib
!
       use systeminf,  only:  xtcf,rep,iat,iatms,igrps
!
       use thresholds, only:  thr,thrang
       use parameters, only:  zero
!
       use geometry,   only:  sminimgvec
!
       use omp_var,    only:  np,chunkadj
!
       implicit none
!
! Input/output variables
!
       logical,dimension(madj,madj),intent(inout)  ::  adj       !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi      !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)        ::  box       !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(mnode),intent(in)         ::  node      !
       integer,dimension(mtype),intent(in)         ::  nnode     !
       integer,dimension(mtype),intent(in)         ::  inode     !
       integer,dimension(mtype),intent(in)         ::  igrpsmon  !
       integer,intent(in)                          ::  madj      !
       integer,intent(in)                          ::  mnode     !
       integer,intent(in)                          ::  mtype     !
       integer,intent(in)                          ::  matms     !
!
! Local variables
!
       real(kind=4),dimension(3)                   ::  v21        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                   ::  v23        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4)                                ::  dis1       !  Minimum image distance
       real(kind=4)                                ::  dis2       !  Minimum image distance
       real(kind=4)                                ::  angle      !
       real(kind=4)                                ::  mindis     !  Distance threshold between groups
       real(kind=4)                                ::  minang     !  Angle threshold between groups
       integer                                     ::  iinode     !
       integer                                     ::  innode     !
       integer                                     ::  jinode     !
       integer                                     ::  jnnode     !
       integer                                     ::  iigrps     !
       integer                                     ::  ingrps     !
       integer                                     ::  jigrps     !
       integer                                     ::  jngrps     !
       integer                                     ::  inei       !
       integer                                     ::  jnei       !
       integer                                     ::  innei      !
       integer                                     ::  jnnei      !
       integer                                     ::  iiat       !
       integer                                     ::  jiat       !
       integer                                     ::  iiatms     !
       integer                                     ::  jiatms     !
       integer                                     ::  ithr       !
       integer                                     ::  jthr       !
       integer                                     ::  iiadj      !
       integer                                     ::  jiadj      !
       integer                                     ::  iadj       !
       integer                                     ::  jadj       !
       integer                                     ::  i,j        !
       integer                                     ::  ii,jj      !
       integer                                     ::  ni,nj      !
       logical                                     ::  doang      !
!
! Building the adjacency matrix in the molecule-based representation
! ------------------------------------------------------------------
!
! Interactions within the same moleculetype
! .........................................
!
       do ii = 1, mtype
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         iiadj  = igrpsmon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,thrang,rep,nnode,       &
!$omp                    inode,ii,iiatms,iiat,ithr,iiadj)              &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     iadj,jadj,ni,nj,i,j)                         &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)-1
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = iiadj  + (iinode-1)*rep(ii)%mgrps
!
           do jinode = iinode+1, nnode(ii)
!
             jnnode = iiatms + (node(inode(ii)+jinode)-1)*rep(ii)%msubg
             jnnei  = iiat   + (node(inode(ii)+jinode)-1)*rep(ii)%nat
             jadj   = iiadj  + (jinode-1)*rep(ii)%mgrps
!
             do ingrps = 1, rep(ii)%mgrps
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jngrps = 1, rep(ii)%mgrps
!
                   mindis = thr(ithr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(ithr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(ii)%igrps(jngrps)
!
                     do jigrps = 1, rep(ii)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(ii)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(ii)%neiang(jnei) .ne. 0 ) then
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(ii)%neiang(jnei)),  &
                                                             box,minang)
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                           end if
!
                         end if
!
                         if ( doang ) then
                           adj(iadj+ingrps,jadj+jngrps) = .TRUE.
                           adj(jadj+jngrps,iadj+ingrps) = .TRUE.
                           GO TO 1000
                         end if
!
                       end if
!
                     end do  !  jigrps
                   end if
!
1000               continue
!
                 end do      !  jngrps
!
               end do        !  iigrps
             end do          !  ingrps
!
           end do            !  jinode
         end do              !  iinode
!
!$omp end parallel do
!
       end do                !  ii
!
! Interactions within different moleculetypes
! ...........................................
!
       do ii = 1, mtype-1
!
         iiatms = iatms(ii)
         iiat   = iat(ii)
         ithr   = igrps(ii)
         iiadj  = igrpsmon(ii)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(node,adj,posi,box,thr,thrang,rep,nnode,       &
!$omp                    inode,iatms,iat,iiadj,igrpsmon,ii,iiatms,     &
!$omp                    iiat,ithr)                                    &
!$omp             private(v21,v23,dis1,dis2,angle,minang,mindis,       &
!$omp                     iinode,innode,jinode,jnnode,iigrps,ingrps,   &
!$omp                     jigrps,jngrps,inei,jnei,innei,jnnei,doang,   &
!$omp                     iadj,jadj,ni,nj,i,j,jj,jiatms,jiat,          &          
!$omp                     jthr,jiadj)                                  &
!$omp             schedule(dynamic,chunkadj)
!
         do iinode = 1, nnode(ii)
!
           innode = iiatms + (node(inode(ii)+iinode)-1)*rep(ii)%msubg
           innei  = iiat   + (node(inode(ii)+iinode)-1)*rep(ii)%nat
           iadj   = iiadj  + (iinode-1)*rep(ii)%mgrps
!
           do jj = ii+1, mtype
!
             jiatms = iatms(jj)
             jiat   = iat(jj)
             jthr   = igrps(jj)
             jiadj  = igrpsmon(jj)
!
           do jinode = 1, nnode(jj)
!
             jnnode = jiatms + (node(inode(jj)+jinode)-1)*rep(jj)%msubg
             jnnei  = jiat   + (node(inode(jj)+jinode)-1)*rep(jj)%nat
             jadj   = jiadj  + (jinode-1)*rep(jj)%mgrps
!
             do ingrps = 1, rep(ii)%mgrps
!
               i = innode + rep(ii)%igrps(ingrps)
!
               do iigrps = 1, rep(ii)%ngrps(ingrps)
!
                 ni   = i + iigrps
                 inei = rep(ii)%igrps(ingrps) + iigrps
!
                 do jngrps = 1, rep(jj)%mgrps
!
                   mindis = thr(jthr+jngrps,ithr+ingrps)
!
                   if ( mindis .gt. zero ) then
!
                     minang = thrang(jthr+jngrps,ithr+ingrps)
!
                     j = jnnode + rep(jj)%igrps(jngrps)
!
                     do jigrps = 1, rep(jj)%ngrps(jngrps)
!
                       nj   = j + jigrps
                       jnei = rep(jj)%igrps(jngrps) + jigrps
!
                       v21  = sminimgvec(posi(:,ni),posi(:,nj),box)
                       dis1 = dot_product(v21,v21)
!
                       if ( dis1 .le. mindis ) then
!
                         doang = .TRUE.
!
                         if ( minang .gt. zero ) then
!
                           if ( rep(jj)%neiang(jnei) .ne. 0 ) then
!
                             v21(:) = -v21(:)
!
                             doang = chkangle(v21,dis1,posi(:,nj),     &
                              xtcf%pos(:,jnnei+rep(jj)%neiang(jnei)),  &
                                                             box,minang)
!
                           else if ( rep(ii)%neiang(inei) .ne. 0 ) then
!
                             doang = chkangle(v21,dis1,posi(:,ni),     &
                              xtcf%pos(:,innei+rep(ii)%neiang(inei)),  &
                                                             box,minang)
!
                           end if
!
                         end if
!
                         if ( doang ) then
                           adj(iadj+ingrps,jadj+jngrps) = .TRUE.
                           adj(jadj+jngrps,iadj+ingrps) = .TRUE.
                           GO TO 2000
                         end if
!
                       end if
!
                     end do    !  jigrps
                   end if
!
2000             continue
!
                 end do        !  jngrps
!
               end do          !  iigrps
             end do            !  ingrps
!
           end do              !  jinode
           end do              !  jj
!
         end do                !  iinode
!
!$omp end parallel do
!
       end do                  !  ii
!
!
       return
       end subroutine nbuildadjgrpsang
!
!======================================================================!
!
! FINDCOMPUNDIR - FIND COMPonents UNDIRected
!
! This subroutine finds the connected components in an undirected
!  unweighted graph of NNODE vertices given as an adjacency matrix
!  representaion ADJ(NODE,NODE) using Breadth First Search.
!
       subroutine findcompundir(nnode,adj,mol,tag,agg,maxagg,nagg,magg)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,dimension(nnode),intent(out)       ::  mol    !  Molecules identifier
       integer,dimension(nnode),intent(out)       ::  tag    !  Aggregates identifier
       integer,dimension(nnode),intent(out)       ::  agg    !  Aggregates size
       integer,dimension(nnode),intent(out)       ::  nagg    !  Number of aggregates of each size
       integer,intent(out)                        ::  magg    !  Number of aggregates
       integer,intent(in)                         ::  nnode   !  Number of molecules
       integer,intent(out)                        ::  maxagg  !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(nnode)                   ::  notvis  !  Nodes visited
       integer,dimension(nnode)                   ::  queue   !  Queue of connected nodes
       integer                                    ::  inode   !  Node index
       integer                                    ::  jnode   !  Node index
       integer                                    ::  knode   !  Node index
       integer                                    ::  iqueue  !  Queue index
       integer                                    ::  nqueue  !  Number of queue elements
       integer                                    ::  nnmol   !  Number of aggregates index
       integer                                    ::  ntag    !  Size of the aggregate
       integer                                    ::  nntag   !  Size of the aggregate index
       integer                                    ::  i       !  Indexes
!
! Performing Breadth First Search over the target molecules
! ---------------------------------------------------------
!
! Marking all the vertices as not visited
       notvis(:) = .TRUE.
! Initializing the molecules information
       mol(:)  = 0
       nnmol   = 1
! Initializing the aggregates information
       ntag    = 0
       nntag   = 0
       agg(:)  = 0
       nagg(:) = 0
! Initializing the size information
       magg    = 0
       tag(:)  = 0
!
       maxagg  = 1
!
! Outer loop over each node
! .........................
!
       do inode = 1, nnode
         if ( notvis(inode) ) then
! Marking head node as visited
           notvis(inode) = .FALSE.
! Updating the system information
           magg       = magg + 1
           tag(nnmol) = magg
!
           mol(nnmol) = inode
           nnmol      = nnmol + 1
!
           ntag       = 1
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
! ..................................
!
           do while ( iqueue .lt. nqueue )
! Saving actual element in the queue
             knode = queue(iqueue)
! Checking the connection between actual queue element and the rest of nodes
             do jnode = inode + 1, nnode
! Checking if node j is connected to node k and has not been already visited
               if ( notvis(jnode) .and. adj(jnode,knode) ) then
! Updating the system information
                 tag(nnmol)   = magg
!
                 mol(nnmol)    = jnode
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
             agg(i) = ntag
           end do
           nntag = nntag + ntag
! Update the number of aggregates of each size
           if ( ntag .gt. maxagg ) maxagg = ntag
           nagg(ntag) = nagg(ntag)   + 1
         end if
       end do
!
       return
       end subroutine findcompundir
!
!======================================================================!
!
! NFINDCOMPUNDIR - N-components FIND COMPonents UNDIRected
!
! This subroutine finds the connected components in an undirected
!  unweighted graph of NNODE vertices given as an adjacency matrix
!  representaion ADJ(NODE,NODE) using Breadth First Search.
!
       subroutine nfindcompundir(adj,mol,node,tag,agg,idx,itype,       &
                                 ntype,nagg,magg,maxagg,maxidx)
!
       use systeminf,   only:  mtype,nnode,inode,mnode
       use properties,  only:  msize,nmax
!
       use mathtools,   only:  ntupla2idx
!
       implicit none
!
! Input/output variables
!
       logical,dimension(mnode,mnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,dimension(mnode),intent(out)       ::  mol     !  Molecules identifier
       integer,dimension(mnode),intent(out)       ::  node    !  Molecules identifier
       integer,dimension(mnode),intent(out)       ::  tag     !  Aggregates identifier
       integer,dimension(mnode),intent(out)       ::  agg     !  Aggregates size
       integer,dimension(mnode),intent(out)       ::  idx     !
       integer,dimension(mnode),intent(out)       ::  itype   !
       integer,dimension(mnode),intent(out)       ::  ntype   !
       integer,dimension(nmax),intent(out)        ::  nagg    !  Number of aggregates of each size
       integer,intent(out)                        ::  magg    !  Number of aggregates
       integer,intent(out)                        ::  maxagg  !  Maximum aggregate size
       integer,intent(out)                        ::  maxidx  !  Maximum aggregate identifier
!
! Local variables
!
       logical,dimension(mnode)                   ::  notvis  !  Nodes visited
       integer,dimension(mnode,mtype)             ::  tmp     !
       integer,dimension(mnode,mtype)             ::  tmpnode !
       integer,dimension(mtype)                   ::  nntype  !
       integer,dimension(mnode)                   ::  queue   !  Queue of connected nodes
       integer,dimension(mnode)                   ::  qtype   !
       integer,dimension(mnode)                   ::  qnode   !
       integer                                    ::  iinode  !  Node index
       integer                                    ::  jinode  !  Node index
       integer                                    ::  kinode  !  Node index
       integer                                    ::  iqueue  !  Queue index
       integer                                    ::  nqueue  !  Number of queue elements
       integer                                    ::  nnmol   !  Current molecule being processed
       integer                                    ::  ntag    !  Current size of the aggregate
       integer                                    ::  nntag   !  Size of the aggregate index
       integer                                    ::  iidx    !  Aggregate identifier
       integer                                    ::  ni,nj   !  Node index
       integer                                    ::  nk      !  Node index
       integer                                    ::  ii,jj   !  Indexes
       integer                                    ::  kk      !  Indexes
       integer                                    ::  i       !  Indexes
!
! Performing Breadth First Search over the target molecules
! ---------------------------------------------------------
!
!~ write(*,*) 'STARTING BFS'
! Marking all the vertices as not visited
       notvis(:) = .TRUE.
! Initializing index information
       idx(:)   = 0
       ntype(:) = 0
       itype(:) = 0
! Initializing the molecules information
       mol(:)   = -1
       node(:)  = -1
       nnmol    = 1
! Initializing the aggregates information
       ntag    = 0
       nntag   = 0
       agg(:)  = 0
       nagg(:) = 0
! Initializing the size information
       magg    = 0
       tag(:)  = 0
!
       maxagg  = 1
       maxidx  = 1
!
! Outer loop over each node
! .........................
!
       do ii = 1, mtype
         do iinode = 1, nnode(ii)
!
           ni = inode(ii) + iinode
!
           if ( notvis(ni) ) then
! Marking head node as visited
             notvis(ni) = .FALSE.
! Updating the system information
             magg       = magg + 1
             tag(nnmol) = magg     ! TODO: check if needed
!
             mol(nnmol)  = ni      ! TODO: check if needed
             node(nnmol) = iinode  ! TODO: check if needed
             nnmol       = nnmol + 1
!
             ntag       = 1
! Initializing queue
             queue(:) = 0
             qtype(:) = 0
             qnode(:) = 0
! Adding current node to the queue
             queue(1) = ni
             qtype(1) = ii
             qnode(1) = iinode
! Initializing the queue counter
             iqueue = 1
! Setting the next position in the queue
             nqueue = 2
! Initializing temporary variables
             tmp(:,:)      = 0
             tmpnode(:,:)  = 0
             nntype(:)     = 0
!
             nntype(ii) = nntype(ii) + 1
             tmp(nntype(ii),ii) = ni
             tmpnode(nntype(ii),ii) = iinode
!
! Inner loop over the queue elements
! ..................................
!
             do while ( iqueue .lt. nqueue )
!
! Saving actual element in the queue
!
               nk     = queue(iqueue)
               kk     = qtype(iqueue)
               kinode = qnode(iqueue)
!~ write(*,*) 'taking from the queue node',nk,'type',kk
!
! Checking the connection between actual queue element and the rest of nodes
!
!   Connection within the same moleculetype as the outer loop node
!   ..............................................................
!
!               do nj = inode(ii)+1, inode(kk)+nnode(kk)
!                 jinode = nj - inode(kk)
               do nj = inode(ii)+iinode+1, inode(ii)+nnode(ii)
                 jinode = nj - inode(ii)
! Checking if node j is connected to node k and has not been already visited
!~ write(*,*) 'checking nk',nk,'nj',nj,'jtype',ii,':',adj(nj,nk)
                 if ( notvis(nj) .and. adj(nj,nk) ) then
!~ write(*,*) '  ACCEPTED',nj
!write(*,*) 'ADDING NODE',nk,nj
! Updating the system information
                   tag(nnmol) = magg
                   nnmol      = nnmol + 1
!
                   ntag       = ntag + 1
!
                   nntype(ii) = nntype(ii) + 1
                   tmp(nntype(ii),ii) = nj
                   tmpnode(nntype(ii),ii) = jinode
! Marking the node connected to node k as visited
                   notvis(nj) = .FALSE.
! Adding to the queue the node connected to node k
                   queue(nqueue) = nj
                   qtype(nqueue) = ii
                   qnode(nqueue) = jinode
! Updating next position in the queue
                   nqueue = nqueue + 1
                 end if
               end do
!
!   Connection within different moleculetypes to the outer loop node
!   ................................................................
!
               do jj = ii+1, mtype
                 do jinode = 1, nnode(jj)
                   nj = inode(jj) + jinode
! Checking if node j is connected to node k and has not been already visited
!~ write(*,*) 'checking nk',nk,'nj',nj,'jtype',jj,':',adj(nj,nk)
                   if ( notvis(nj) .and. adj(nj,nk) ) then
!write(*,*) 'ADDING NODE',nk,nj
!~ write(*,*) '  ACCEPTED'
! Updating the system information
                     tag(nnmol) = magg
                     nnmol      = nnmol + 1
!
                     ntag       = ntag + 1
!
                     nntype(jj) = nntype(jj) + 1
                     tmp(nntype(jj),jj) = nj
                     tmpnode(nntype(jj),jj) = jinode
! Marking the node connected to node k as visited
                     notvis(nj) = .FALSE.
! Adding to the queue the node connected to node k
                     queue(nqueue) = nj
                     qtype(nqueue) = jj
                     qnode(nqueue) = jinode
! Updating next position in the queue
                     nqueue = nqueue + 1
                   end if
                 end do
               end do
!
! Updating the queue counter
!
               iqueue = iqueue + 1
             end do
!~ write(*,*)
!
! Saving aggregate information
! ............................
!
! Computing the aggregate identifier
             iidx = ntupla2idx(mtype,nntype)
! Saving the information of the aggregate found
             agg(nntag+1:nntag+ntag) = ntag
             idx(nntag+1:nntag+ntag) = iidx
!~ if (ntag.gt.1) write(*,*) 'NEW AGG',ntag,iidx,':',nntype
!
             do i = 1, mtype
!
               mol(nntag+1:nntag+nntype(i))   = tmp(:nntype(i),i)
               node(nntag+1:nntag+nntype(i))  = tmpnode(:nntype(i),i)
               itype(nntag+1:nntag+nntype(i)) = i
               ntype(nntag+1:nntag+nntype(i)) = nntype(i)
!~ if (ntag.gt.1) write(*,*) 'mol',tmp(:nntype(i),i)
!
               nntag = nntag + nntype(i)
             end do
! Update the number of aggregates of each size
             if ( ntag .gt. maxagg ) maxagg = ntag
             if ( iidx .gt. maxidx ) maxidx = iidx
!
             if ( ntag .le. msize ) then
               nagg(iidx) = nagg(iidx) + 1
             else
               nagg(nmax) = nagg(nmax) + 1
             end if
!
           end if
!
         end do
       end do
!~ write(*,*)
!
       return
       end subroutine nfindcompundir
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
! BONDS2ADJ - BOND termS TO ADJacency matrix
!
! This subroutine
!
       subroutine bonds2adj(nbond,ibond,nat,adj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(out)  ::  adj    !  Adjacency matrix
       integer,dimension(2,nbond),intent(in)   ::  ibond  !  Atom indexes
       integer,intent(in)                      ::  nbond  !  Number of edges
       integer,intent(in)                      ::  nat    !  Number of nodes
!
! Local variables
!
       integer                                     ::  i  !  Index
!
! Building adjacency matrix according to Wiberg bond indexes
!
       adj(:,:) = .FALSE.
!
       do i = 1, nbond
         adj(ibond(1,i),ibond(2,i)) = .TRUE.
         adj(ibond(2,i),ibond(1,i)) = .TRUE.
       end do
!
       return
       end subroutine bonds2adj
!
!======================================================================!
!
! ADJ2ADJATMS - ADJacency matrix in canonical order
!                TO ADJacency matrix of reference AToMS
!
! This subroutine
!
       subroutine adj2adjatms(nat,inadj,matms,atms,outadj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)       ::  inadj   !  Adjacency matrix
       logical,dimension(matms,matms),intent(out)  ::  outadj  !  Adjacency matrix
       integer,dimension(matms),intent(in)         ::  atms    !  Indexes of reference nodes
       integer,intent(in)                          ::  nat     !  Number of nodes
       integer,intent(in)                          ::  matms   !  Number of reference nodes
!
! Local variables
!
       integer                                     ::  i,j     !  Index
!
! Building adjacency matrix
!
       outadj(:,:) = .FALSE.
!
       do i = 1, matms-1
         do j = i+1, matms
           if ( inadj(atms(i),atms(j))                                 &
                                      .or. (atms(i).eq.(atms(j))) ) then
             outadj(i,j) = .TRUE.
             outadj(j,i) = .TRUE.
           end if
         end do
       end do
!
       return
       end subroutine adj2adjatms
!
!======================================================================!
!
! REDUCEADJ - REDUCE ADJacency matrix to a lower representation
!
! This subroutine
!
       subroutine reduceadj(matms,adjatms,msubg,nsubg,isubg,adjsubg)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(matms,matms),intent(in)   ::  adjatms  !  Adjacency matrix
       logical,dimension(msubg,msubg),intent(out)  ::  adjsubg  !  Adjacency matrix
       integer,dimension(msubg),intent(in)         ::  nsubg    !
       integer,dimension(msubg),intent(in)         ::  isubg    !
       integer,intent(in)                          ::  matms    !  Number of reference nodes
       integer,intent(in)                          ::  msubg    !  Number of reference nodes
!
! Local variables
!
       integer                                     ::  i,j      !  Index
       integer                                     ::  ii,jj    !  Index
       integer                                     ::  iisubg   !  Index
       integer                                     ::  jjsubg   !  Index
!
! Building adjacency matrix
!
       adjsubg(:,:) = .FALSE.
!
       do i = 1, msubg-1
         do j = i+1, msubg
!
           do iisubg = 1, nsubg(i)
             ii = isubg(i) + iisubg
             do jjsubg = 1, nsubg(j)
               jj = isubg(j) + jjsubg
!
               if ( adjatms(ii,jj) ) then
                 adjsubg(i,j) = .TRUE.
                 adjsubg(j,i) = .TRUE.
                 GOTO 1000
               end if
!
             end do
           end do
!
1000       continue
!
         end do
       end do
!
       return
       end subroutine reduceadj
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
! CHKANGLE - CHecK ANGLE
!
! This subroutine
!
       logical function chkangle(v21,dis1,r2,r3,box,minang)
!
       use parameters
       use geometry,  only:  sminimgvec
!
       implicit none
!
! Input/output variables
!
       real(kind=4),dimension(3),intent(in)  ::  v21     !
       real(kind=4),dimension(3),intent(in)  ::  r2      !
       real(kind=4),dimension(3),intent(in)  ::  r3      !
       real(kind=4),dimension(3),intent(in)  ::  box     !
       real(kind=4),intent(in)               ::  dis1    !
       real(kind=4),intent(in)               ::  minang  !
!
! Local variables
!
       real(kind=4),dimension(3)             ::  v23     !
       real(kind=4)                          ::  dis2    !
       real(kind=4)                          ::  angle   !
!
!
!
       chkangle = .FALSE.
!
       v23(:) = sminimgvec(r2,r3,box)
       dis2   = dot_product(v23,v23)
!
       dis2 = sqrt(dis1*dis2)
!
       angle = (v21(1)*v23(1) + v21(2)*v23(2) + v21(3)*v23(3))/dis2
       angle = acos(angle)
!
       if ( angle .gt. minang ) chkangle = .TRUE.
!
       return
       end function chkangle
!
!======================================================================!
!
       end module graphtools
!
!======================================================================!
