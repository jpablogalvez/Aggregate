!======================================================================!
!
       module aggtools
!
       implicit none
!
       contains
!
!======================================================================!
!
! AGGSCRN - AGGregates SCReeNing algorithm
!
! This subroutine obtains the size and the number of aggregates present
!  in the system according to a given distance criteria.
! The interaction between two molecules is stored in the corresponding
!  adjacency matrix in the molecule-based representation. After block-
!  diagonalization of this matrix, the size aggregates is obtained from
!  the dimension of each block and the number of aggregates of each size
!  is obtained from the number of blocks of each size.
! The subroutine SCRNINT is employed to screen the interactions between
!  the molecules. To do so a three-steps procedure is carried out with
!  the following phylosophy: When an interaction exists, if I was or I 
!  will then I am, i.e., each interaction is anhilated if it is not pre-
!  sent in two consequtive configurations (removing collisions), and 
!  when an interaction does not exist, then if it was and it will then 
!  it is, i.e., one interaction can be created if it is present in the 
!  former and the following configurations (adding oscillations).
!
       subroutine aggscrn(sys,xtcf,nat,nnode,natms,thr,thr2,neidis,    &
                          pim,msize,pop,conc,frac,cin,volu,nsteps,     &
                          grptag,nbody,ngrps,nsubg,ibody,igrps,isubg,  &
                          body,grps,subg,atms,mbody,mgrps,msubg,matms, &
                          nprint,minstep,maxstep,nsolv,dopim,debug)
!
       use xdr, only: xtcfile
!
       use screening
       use datatypes
       use timings
       use printings
       use utils
!
       implicit none
!
! System information
!
       type(groinp),intent(in)                                  ::  sys      !  Monomer information
       type(xtcfile),intent(inout)                              ::  xtcf     !  Trajectory information
!
! Interaction criteria information
! 
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr      !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr2     !  Distance threshold
       real(kind=8),intent(in)                                  ::  neidis   !  Screening distance
!
! Average properties
!
       real(kind=8),dimension(mgrps,mgrps,msize-1),intent(out)  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(msize),intent(out)                ::  pop      !  Populations
       real(kind=8),dimension(msize),intent(out)                ::  conc     !  Concentrations
       real(kind=8),dimension(msize),intent(out)                ::  frac     !  Molar fractions
       real(kind=8),intent(out)                                 ::  cin      !  Stechiometric concentration
       real(kind=8),intent(out)                                 ::  volu     !  Simulation box volume
       integer,intent(out)                                      ::  nsteps   !  Number of snapshots analyzed
       integer,intent(in)                                       ::  msize    !  Maximum aggregate size
!
! Trajectory control variables
!     
       integer,intent(in)                                       ::  nprint   !  Populations printing interval
       integer,intent(in)                                       ::  minstep  !  First step for analysis
       integer,intent(in)                                       ::  maxstep  !  Last step for analysis
!
! Topological representations information
!
       character(len=lentag),dimension(nat),intent(in)          ::  grptag   !  Names of the groups
       integer,dimension(nat),intent(in)                        ::  body     !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  nbody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  ibody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  grps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  ngrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  igrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  subg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  nsubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  isubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  atms     !  Atoms identifier
       integer,intent(in)                                       ::  nat      !  Monomer atoms
       integer,intent(in)                                       ::  natms    !  Total number of subgroups in the system
       integer,intent(in)                                       ::  mbody    !  Number of bodies
       integer,intent(in)                                       ::  mgrps    !  Number of groups
       integer,intent(in)                                       ::  msubg    !  Number of subgroups
       integer,intent(in)                                       ::  matms    !  Number of interacting atoms in the monomer
!
! Program control flags
!
       logical,intent(in)                                       ::  dopim    !  PIM calculation flag
       logical,intent(in)                                       ::  debug    !  Debug mode
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable                       ::  adj      !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  oldadj   !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  newadj   !  Adjacency matrix
       integer,dimension(:),allocatable                         ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable                         ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  agg      !  Aggregates size 
       integer,dimension(:),allocatable                         ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  iagg     !  
       integer,dimension(:),allocatable                         ::  nmol     !  
       integer,dimension(:),allocatable                         ::  imol     !  
       integer                                                  ::  nnode    !  Total number of molecules
       integer                                                  ::  nsize    !  Actual maximum aggregate size
       integer                                                  ::  magg     !  Actual number of chemical species
       integer                                                  ::  actstep  !
       integer                                                  ::  newstep  !
!
! Declaration of time control variables
!
       integer                                                  ::  t1read   !  Initial reading time
       integer                                                  ::  t2read   !  Final reading time
       integer                                                  ::  t1adj    !  Initial building time
       integer                                                  ::  t2adj    !  Final building time
       integer                                                  ::  t1scrn   !  Initial screening time
       integer                                                  ::  t2scrn   !  Final screening time
       integer                                                  ::  t1pim    !  Initial PIM analysis time
       integer                                                  ::  t2pim    !  Final PIM analysis time
!
! AnalysisPhenolMD variables
!     
       integer,intent(in)                                       ::  nsolv    !
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable                  ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable                  ::  newposi  !  Auxiliary coordinates
       real(kind=4),dimension(3)                                ::  box      !
       real(kind=4),dimension(3)                                ::  newbox   !
       integer                                                  ::  i        !
!
! Initializing variables
!
       call system_clock(t1read)        
!
       nsteps  = 0
!
       pim(:,:,:) = 0.0d0
!
       pop(:)  = 0.0d0
       conc(:) = 0.0d0
       frac(:) = 0.0d0
!
       cin     = 0.0d0
       volu    = 0.0d0
!
! Allocating variables depending on system size
!
       allocate(adj(nnode,nnode))
       allocate(oldadj(nnode,nnode),newadj(nnode,nnode))
!
       allocate(mol(nnode),tag(nnode),agg(nnode))
       allocate(nmol(nnode),imol(nnode))
       allocate(nagg(nnode),iagg(nnode))
!
! Allocating variables depending on topological information 
!
       allocate(newposi(3,natms))
       allocate(posi(3,natms))
!
! Reading the first old-configuration
!
       do while ( xtcf%STEP .lt. minstep )
         call xtcf%read
         if ( xtcf%STAT .ne. 0 ) then
           write(*,*)
           write(*,*) 'ERROR: Not enough steps'
           write(*,*) 
           call print_end()
         end if
       end do
!
       box  = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the first old-configuration
!
       call system_clock(t1adj) 
!
       call buildadj(nnode,oldadj,natms,posi,xtcf%NATOMS,xtcf%pos,     &
                     sys%nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,       &
                     isubg,atms,box,neidis)
!
       call system_clock(t2adj) 
!
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)    
!
! Reading the first configuration
!
       call system_clock(t1read)
!
       call xtcf%read
!
       if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.ge.maxstep) ) then
         write(*,*)
         write(*,*) 'ERROR: Not enough steps'
         write(*,*) 
         call print_end()
       end if
!
       do while ( mod(xtcf%STEP-minstep,nprint) .ne. 0 ) 
         call xtcf%read
         if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.ge.maxstep) ) then
           write(*,*)
           write(*,*) 'ERROR: Not enough steps'
           write(*,*) 
           call print_end()
         end if
       end do
!
       box(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       actstep = xtcf%STEP
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the first configuration
!
       call system_clock(t1adj) 
!
       call buildadj(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,        &
                     sys%nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,       &
                     isubg,atms,box,neidis)
!
       call system_clock(t2adj) 
!
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Reading the first new-configuration
!
       call system_clock(t1read)
!
       call xtcf%read
!
       if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.gt.maxstep) ) then
         write(*,*)
         write(*,*) 'ERROR:: Not enough steps'
         write(*,*) 
         call print_end()
       end if
!
       do while ( mod(xtcf%STEP-minstep,nprint) .ne. 0 ) 
         call xtcf%read
         if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.gt.maxstep) ) then
           write(*,*)
           write(*,*) 'ERROR:: Not enough steps'
           write(*,*) 
           call print_end()
         end if
       end do
!
! Analyzing frames in the inverval [minstep,maxstep]
!
       do while ( (xtcf%STAT.eq.0) .and. (xtcf%STEP.le.maxstep) )
         if ( mod(xtcf%STEP-minstep,nprint) .eq. 0 ) then
!
           newbox(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
           newstep   = xtcf%STEP
! 
           nsteps = nsteps + 1
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the new-configuration
!
           call system_clock(t1adj) 
!
           call buildadj(nnode,newadj,natms,newposi,xtcf%NATOMS,       &
                         xtcf%pos,sys%nat,thr2,mgrps,ngrps,igrps,      &
                         msubg,nsubg,isubg,atms,newbox,neidis)
!
           call system_clock(t2adj) 
!
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Screening interactions between the molecules
!
           call system_clock(t1scrn)
!
           call scrnint(nnode,oldadj,adj,newadj)
!
           call system_clock(t2scrn)
!
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
! Block-diagonalizing the interaction-corrected adjacency matrix
!
           call blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg,       &
                          nmol,imol,magg)
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call printpop(actstep,nsize,msize,pop,conc,frac,cin,volu,   &
                         box,nnode,nagg,magg,nsolv,uniout)
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 

!~            if ( debug ) then
!~              call print_coord(xtcf,sys,outp,msize,nagg,nnode,mol,agg)
!~            end if
!
! Adding up the pairwise interaction matrix of the current snapshot
! 
           if ( dopim ) then
             call system_clock(t1pim)     
!
!~              call build_pim(ngrps,grps,nsubg,subg,sys%nat,atms,thr,   &
!~                             nnode,mol,agg,msize,nagg,xtcf%NATOMS,   &
!~                             xtcf%pos,(/xtcf%box(1,1),xtcf%box(2,2),    &
!~                             xtcf%box(3,3)/),pim,debug)
!
             call system_clock(t2pim)     
!
             tpim = tpim + dble(t2pim-t1pim)/dble(count_rate)   
           end if      
!
! Analyzing aggregates by their connectivity
!
!~            call analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,    &
!~                             sys%nat,atms,thr,nnode,adj,mol,agg,     & 
!~                             msize,nagg,xtcf%NATOMS,xtcf%pos,          &
!~                             (/xtcf%box(1,1),xtcf%box(2,2),             &
!~                               xtcf%box(3,3)/),posi,sys%mass,          &
!~                             sys%atname,xtcf%STEP,outp,debug)
!
           call system_clock(t1scrn)
!
           oldadj(:,:) = adj(:,:)
           adj(:,:)    = newadj(:,:)
!
           box(:)    = newbox(:)
!
           actstep = newstep
!
           posi(:,:) = newposi(:,:)
!
           call system_clock(t2scrn)
!
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
           call system_clock(t1read)
!
         end if
!
! Reading information of the next snapshot
!
         call xtcf%read
!
       end do
!
! Deallocate memory
!
       deallocate(newposi)
       deallocate(posi)
!
       deallocate(adj,oldadj,newadj)
!
       deallocate(mol,tag,agg)
       deallocate(nmol,imol,nagg,iagg)
!
! Saving timings
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       return
       end subroutine aggscrn
!
!======================================================================!
!
! AGGDIST - AGGregates DISTances algorithm
!
! This subroutine obtains the size and the number of aggregates present
!  in the system according to a given distance criteria.
! The interaction between two molecules is stored in the corresponding
!  adjacency matrix in the molecule-based representation. After block-
!  diagonalization of this matrix, the size aggregates is obtained from
!  the dimension of each block and the number of aggregates of each size
!  is obtained from the number of blocks of each size.
!
       subroutine aggdist(sys,xtcf,nat,nnode,natms,thr,thr2,neidis,    &
                          pim,msize,pop,conc,frac,cin,volu,nsteps,     &
                          grptag,nbody,ngrps,nsubg,ibody,igrps,isubg,  &
                          body,grps,subg,atms,mbody,mgrps,msubg,matms, &
                          nprint,minstep,maxstep,nsolv,dopim,debug)
!
       use xdr, only: xtcfile
!
       use screening
       use datatypes
       use timings
       use printings
       use utils
!
       implicit none
!
! System information
!
       type(groinp),intent(in)                                  ::  sys      !  Monomer information
       type(xtcfile),intent(inout)                              ::  xtcf     !  Trajectory information
!
! Interaction criteria information
! 
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr      !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr2     !  Distance threshold
       real(kind=8),intent(in)                                  ::  neidis   !  Screening distance
!
! Average properties
!
       real(kind=8),dimension(mgrps,mgrps,msize-1),intent(out)  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(msize),intent(out)                ::  pop      !  Populations
       real(kind=8),dimension(msize),intent(out)                ::  conc     !  Concentrations
       real(kind=8),dimension(msize),intent(out)                ::  frac     !  Molar fractions
       real(kind=8),intent(out)                                 ::  cin      !  Stechiometric concentration
       real(kind=8),intent(out)                                 ::  volu     !  Simulation box volume
       integer,intent(out)                                      ::  nsteps   !  Number of snapshots analyzed
       integer,intent(in)                                       ::  msize    !  Maximum aggregate size
!
! Trajectory control variables
!     
       integer,intent(in)                                       ::  nprint   !  Populations printing interval
       integer,intent(in)                                       ::  minstep  !  First step for analysis
       integer,intent(in)                                       ::  maxstep  !  Last step for analysis
!
! Topological representations information
!
       character(len=lentag),dimension(nat),intent(in)          ::  grptag   !  Names of the groups
       integer,dimension(nat),intent(in)                        ::  body     !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  nbody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  ibody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  grps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  ngrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  igrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  subg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  nsubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  isubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  atms     !  Atoms identifier
       integer,intent(in)                                       ::  nat      !  Monomer atoms
       integer,intent(in)                                       ::  natms    !  Total number of subgroups in the system
       integer,intent(in)                                       ::  mbody    !  Number of bodies
       integer,intent(in)                                       ::  mgrps    !  Number of groups
       integer,intent(in)                                       ::  msubg    !  Number of subgroups
       integer,intent(in)                                       ::  matms    !  Number of interacting atoms in the monomer
!
! Program control flags
!
       logical,intent(in)                                       ::  dopim    !  PIM calculation flag
       logical,intent(in)                                       ::  debug    !  Debug mode
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable                       ::  adj      !  Adjacency matrix
       integer,dimension(:),allocatable                         ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable                         ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  agg      !  Aggregates size 
       integer,dimension(:),allocatable                         ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  iagg     !  
       integer,dimension(:),allocatable                         ::  nmol     !  
       integer,dimension(:),allocatable                         ::  imol     !  
       integer                                                  ::  nnode    !  Total number of molecules
       integer                                                  ::  nsize    !  Actual maximum aggregate size
       integer                                                  ::  magg     !  Actual number of chemical species
!
! Declaration of time control variables
!
       integer                                                  ::  t1read   !  Initial reading time
       integer                                                  ::  t2read   !  Final reading time
       integer                                                  ::  t1adj    !  Initial building time
       integer                                                  ::  t2adj    !  Final building time
       integer                                                  ::  t1pim    !  Initial PIM analysis time
       integer                                                  ::  t2pim    !  Final PIM analysis time
!
! AnalysisPhenolMD variables
!     
       integer,intent(in)                                       ::  nsolv    !
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable                  ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(3)                                ::  box      !
       integer                                                  ::  i        !
!
! Initializing variables
!
       call system_clock(t1read)        
!
       nsteps  = 0
!
       pim(:,:,:) = 0.0d0
!
       pop(:)  = 0.0d0
       conc(:) = 0.0d0
       frac(:) = 0.0d0
!
       cin     = 0.0d0
       volu    = 0.0d0
!
! Allocating variables depending on system size
!
       allocate(adj(nnode,nnode))
!
       allocate(mol(nnode),tag(nnode),agg(nnode))
       allocate(nmol(nnode),imol(nnode))
       allocate(nagg(nnode),iagg(nnode))
!
! Allocating variables depending on topological information 
!
       allocate(posi(3,natms))
!
! Analyzing frames in the inverval [minstep,maxstep]
!
       do while ( (xtcf%STAT.eq.0) .and. (xtcf%STEP.le.maxstep) )
         if ( (mod(xtcf%STEP-minstep,nprint).eq.0) .and.               &
                                           (xtcf%STEP.ge.minstep) ) then
!
           box(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
! 
           nsteps = nsteps + 1
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the new-configuration
!
           call system_clock(t1adj) 
!
           call buildadj(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,    &
                         sys%nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,   &
                         isubg,atms,box,neidis)
!
           call system_clock(t2adj) 
!
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Block-diagonalizing the interaction-corrected adjacency matrix
!
           call blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg,       &
                          nmol,imol,magg)
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call printpop(xtcf%STEP,nsize,msize,pop,conc,frac,cin,volu, &
                         box,nnode,nagg,magg,nsolv,uniout)
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 

!~            if ( debug ) then
!~              call print_coord(xtcf,sys,outp,msize,nagg,nnode,mol,agg)
!~            end if
!
! Adding up the pairwise interaction matrix of the current snapshot
! 
           if ( dopim ) then
             call system_clock(t1pim)     
!
!~              call build_pim(ngrps,grps,nsubg,subg,sys%nat,atms,thr,   &
!~                             nnode,mol,agg,msize,nagg,xtcf%NATOMS,   &
!~                             xtcf%pos,(/xtcf%box(1,1),xtcf%box(2,2),    &
!~                             xtcf%box(3,3)/),pim,debug)
!
             call system_clock(t2pim)     
!
             tpim = tpim + dble(t2pim-t1pim)/dble(count_rate)   
           end if      
!
! Analyzing aggregates by their connectivity
!
!~            call analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,    &
!~                             sys%nat,atms,thr,nnode,adj,mol,agg,     & 
!~                             msize,nagg,xtcf%NATOMS,xtcf%pos,          &
!~                             (/xtcf%box(1,1),xtcf%box(2,2),             &
!~                               xtcf%box(3,3)/),posi,sys%mass,          &
!~                             sys%atname,xtcf%STEP,outp,debug)
!
           call system_clock(t1read)
!
         end if
!
! Reading information of the next snapshot
!
         call xtcf%read
!
       end do
!
! Deallocate memory
!
       deallocate(posi)
!
       deallocate(adj)
!
       deallocate(mol,tag,agg)
       deallocate(nmol,imol,nagg,iagg)
!
! Saving timings
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       return
       end subroutine aggdist
!
!======================================================================!
!
! AGGSCRNLIFE - AGGregates SCReeNing algorithm for LIFEtimes
!
! This subroutine obtains the size and the number of aggregates present
!  in the system according to a given distance criteria.
! The interaction between two molecules is stored in the corresponding
!  adjacency matrix in the molecule-based representation. After block-
!  diagonalization of this matrix, the size aggregates is obtained from
!  the dimension of each block and the number of aggregates of each size
!  is obtained from the number of blocks of each size.
! The subroutine SCRNINT is employed to screen the interactions between
!  the molecules. To do so a three-steps procedure is carried out with
!  the following phylosophy: When an interaction exists, if I was or I 
!  will then I am, i.e., each interaction is anhilated if it is not pre-
!  sent in two consequtive configurations (removing collisions), and 
!  when an interaction does not exist, then if it was and it will then 
!  it is, i.e., one interaction can be created if it is present in the 
!  former and the following configurations (adding oscillations).
! The lifetimes are calculated keeping track of the aggregates which
!  are present in the former and the previous configurations.
!
       subroutine aggscrnlife(sys,xtcf,nat,nnode,natms,thr,thr2,       &
                              neidis,pim,msize,pop,conc,frac,cin,volu, &
                              nsteps,grptag,nbody,ngrps,nsubg,ibody,   &
                              igrps,isubg,body,grps,subg,atms,mbody,   &
                              mgrps,msubg,matms,nprint,minstep,        &
                              maxstep,nsolv,dopim,debug)
!
       use xdr, only: xtcfile
!
       use screening
       use datatypes
       use timings
       use printings
       use utils
!
       implicit none
!
! System information
!
       type(groinp),intent(in)                                  ::  sys      !  Monomer information
       type(xtcfile),intent(inout)                              ::  xtcf     !  Trajectory information
!
! Interaction criteria information
! 
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr      !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr2     !  Distance threshold
       real(kind=8),intent(in)                                  ::  neidis   !  Screening distance
!
! Average properties
!
       real(kind=8),dimension(mgrps,mgrps,msize-1),intent(out)  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(msize),intent(out)                ::  pop      !  Populations
       real(kind=8),dimension(msize),intent(out)                ::  conc     !  Concentrations
       real(kind=8),dimension(msize),intent(out)                ::  frac     !  Molar fractions
       real(kind=8),intent(out)                                 ::  cin      !  Stechiometric concentration
       real(kind=8),intent(out)                                 ::  volu     !  Simulation box volume
       integer,intent(out)                                      ::  nsteps   !  Number of snapshots analyzed
       integer,intent(in)                                       ::  msize    !  Maximum aggregate size
!
! Trajectory control variables
!     
       integer,intent(in)                                       ::  nprint   !  Populations printing interval
       integer,intent(in)                                       ::  minstep  !  First step for analysis
       integer,intent(in)                                       ::  maxstep  !  Last step for analysis
!
! Topological representations information
!
       character(len=lentag),dimension(nat),intent(in)           ::  grptag   !  Names of the groups
       integer,dimension(nat),intent(in)                         ::  body     !  Number of groups in each body
       integer,dimension(nat),intent(in)                         ::  nbody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                         ::  ibody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                         ::  grps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                         ::  ngrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                         ::  igrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                         ::  subg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                         ::  nsubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                         ::  isubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                         ::  atms     !  Atoms identifier
       integer,intent(in)                                        ::  nat      !  Monomer atoms
       integer,intent(in)                                        ::  natms    !  Total number of subgroups in the system
       integer,intent(in)                                        ::  mbody    !  Number of bodies
       integer,intent(in)                                        ::  mgrps    !  Number of groups
       integer,intent(in)                                        ::  msubg    !  Number of subgroups
       integer,intent(in)                                        ::  matms    !  Number of interacting atoms in the monomer
!
! Program control flags
!
       logical,intent(in)                                        ::  dopim    !  PIM calculation flag
       logical,intent(in)                                        ::  debug    !  Debug mode
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable                       ::  adj      !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  oldadj   !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  newadj   !  Adjacency matrix
       integer,dimension(:),allocatable                         ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable                         ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  agg      !  Aggregates size 
       integer,dimension(:),allocatable                         ::  oldmol   !  Molecules identifier
       integer,dimension(:),allocatable                         ::  oldtag   !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  oldagg   !  Aggregates size 
       integer,dimension(:),allocatable                         ::  newmol   !  Molecules identifier
       integer,dimension(:),allocatable                         ::  newtag   !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  newagg   !  Aggregates size 
       integer,dimension(:),allocatable                         ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  iagg     !  
       integer,dimension(:),allocatable                         ::  nmol     !  
       integer,dimension(:),allocatable                         ::  imol     !  
       integer,dimension(:),allocatable                         ::  oldnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  oldiagg  !  Number of aggregates of lower size
       integer,dimension(:),allocatable                         ::  oldnmol  !  
       integer,dimension(:),allocatable                         ::  oldimol  ! 
       integer,dimension(:),allocatable                         ::  newnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  newiagg  !  
       integer,dimension(:),allocatable                         ::  newnmol  !  
       integer,dimension(:),allocatable                         ::  newimol  ! 
       integer,dimension(:),allocatable                         ::  wasmap   !
       integer,dimension(:),allocatable                         ::  willmap  !
       integer                                                  ::  nnode    !  Total number of molecules
       integer                                                  ::  nsize    !  Actual maximum aggregate size
       integer                                                  ::  oldsize  !  Old maximum aggregate size
       integer                                                  ::  newsize  !  New maximum aggregate size
       integer                                                  ::  magg     !  Actual number of chemical species
       integer                                                  ::  newmagg  !  New number of chemical species
       integer                                                  ::  oldmagg  !  Old number of chemical species
       integer                                                  ::  oldstep  !
       integer                                                  ::  actstep  !
       integer                                                  ::  newstep  !
       logical,dimension(:),allocatable                         ::  iam      !
       logical,dimension(:),allocatable                         ::  iwas     !
       logical,dimension(:),allocatable                         ::  iwill    !
       logical,dimension(:),allocatable                         ::  imnot    !
       logical,dimension(:),allocatable                         ::  iwasnt   !
       logical,dimension(:),allocatable                         ::  iwont    !
       logical,dimension(:),allocatable                         ::  oldiam   !
       logical,dimension(:),allocatable                         ::  oldnot   !
!
! Declaration of time control variables
!
       integer                                                  ::  t1read   !  Initial reading time
       integer                                                  ::  t2read   !  Final reading time
       integer                                                  ::  t1adj    !  Initial building time
       integer                                                  ::  t2adj    !  Final building time
       integer                                                  ::  t1scrn   !  Initial screening time
       integer                                                  ::  t2scrn   !  Final screening time
       integer                                                  ::  t1pim    !  Initial PIM analysis time
       integer                                                  ::  t2pim    !  Final PIM analysis time
!
! AnalysisPhenolMD variables
!     
       integer,intent(in)                                       ::  nsolv    !
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable                  ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable                  ::  newposi  !  Auxiliary coordinates
       real(kind=4),dimension(3)                                ::  box      !
       real(kind=4),dimension(3)                                ::  oldbox   !
       real(kind=4),dimension(3)                                ::  newbox   !
       integer                                                  ::  i        !
!
! Initializing variables
!
       call system_clock(t1read)        
!
       nsteps  = 0
!
       pim(:,:,:) = 0.0d0
!
       pop(:)  = 0.0d0
       conc(:) = 0.0d0
       frac(:) = 0.0d0
!
       cin     = 0.0d0
       volu    = 0.0d0
!
! Allocating variables depending on system size
!
       allocate(adj(nnode,nnode))
       allocate(oldadj(nnode,nnode),newadj(nnode,nnode))
!
       allocate(mol(nnode),tag(nnode),agg(nnode))
       allocate(nmol(nnode),imol(nnode))
       allocate(nagg(nnode),iagg(nnode))
!
       allocate(oldmol(nnode),oldtag(nnode),oldagg(nnode))
       allocate(oldnmol(nnode),oldimol(nnode))
       allocate(oldnagg(nnode),oldiagg(nnode))
!
       allocate(newmol(nnode),newtag(nnode),newagg(nnode))
       allocate(newnmol(nnode),newimol(nnode))
       allocate(newnagg(nnode),newiagg(nnode))
!
       allocate(iam(nnode),iwas(nnode),iwill(nnode))
       allocate(imnot(nnode),iwasnt(nnode),iwont(nnode))
       allocate(oldiam(nnode),oldnot(nnode))
!
       allocate(wasmap(nnode),willmap(nnode))
!
! Allocating variables depending on topological information 
!
       allocate(newposi(3,natms))
       allocate(posi(3,natms))
!
! Reading the first old-configuration
!
       do while ( xtcf%STEP .lt. minstep )
         call xtcf%read
         if ( xtcf%STAT .ne. 0 ) then
           write(*,*)
           write(*,*) 'ERROR: Not enough steps'
           write(*,*) 
           call print_end()
         end if
       end do
!
       oldbox  = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       oldstep = xtcf%STEP
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the first old-configuration
!
       call system_clock(t1adj) 
!
       call buildadj(nnode,oldadj,natms,posi,xtcf%NATOMS,xtcf%pos,     &
                     sys%nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,       &
                     isubg,atms,oldbox,neidis)
!
       call system_clock(t2adj) 
!
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)    
!
! Reading the first configuration
!
       call system_clock(t1read)
!
       call xtcf%read
!
       if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.ge.maxstep) ) then
         write(*,*)
         write(*,*) 'ERROR: Not enough steps'
         write(*,*) 
         call print_end()
       end if
!
       do while ( mod(xtcf%STEP-minstep,nprint) .ne. 0 ) 
         call xtcf%read
         if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.ge.maxstep) ) then
           write(*,*)
           write(*,*) 'ERROR: Not enough steps'
           write(*,*) 
           call print_end()
         end if
       end do
!
       box(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       actstep = xtcf%STEP
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the first configuration
!
       call system_clock(t1adj) 
!
       call buildadj(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,        &
                     sys%nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,       &
                     isubg,atms,box,neidis)
!
       call system_clock(t2adj) 
!
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Reading the first new-configuration
!
       call system_clock(t1read)
!
       call xtcf%read
!
       if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.gt.maxstep) ) then
         write(*,*)
         write(*,*) 'ERROR:: Not enough steps'
         write(*,*) 
         call print_end()
       end if
!
       do while ( mod(xtcf%STEP-minstep,nprint) .ne. 0 ) 
         call xtcf%read
         if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.gt.maxstep) ) then
           write(*,*)
           write(*,*) 'ERROR:: Not enough steps'
           write(*,*) 
           call print_end()
         end if
       end do
!
! Analyzing frames in the inverval [minstep,maxstep]
!
       do while ( (xtcf%STAT.eq.0) .and. (xtcf%STEP.le.maxstep) )
         if ( mod(xtcf%STEP-minstep,nprint) .eq. 0 ) then
!
           newbox(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
           newstep   = xtcf%STEP
! 
           nsteps = nsteps + 1
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the new-configuration
!
           call system_clock(t1adj) 
!
           call buildadj(nnode,newadj,natms,newposi,xtcf%NATOMS,       &
                         xtcf%pos,sys%nat,thr2,mgrps,ngrps,igrps,      &
                         msubg,nsubg,isubg,atms,newbox,neidis)
!
           call system_clock(t2adj) 
!
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Screening interactions between the molecules
!
           call system_clock(t1scrn)
!
           call scrnint(nnode,oldadj,adj,newadj)
!
           call system_clock(t2scrn)
!
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
! Block-diagonalizing the interaction-corrected adjacency matrix
!
           call blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg,       &
                          nmol,imol,magg)
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call printpop(actstep,nsize,msize,pop,conc,frac,cin,volu,   &
                         box,nnode,nagg,magg,nsolv,uniout)
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 

!~            if ( debug ) then
!~              call print_coord(xtcf,sys,outp,msize,nagg,nnode,mol,agg)
!~            end if
!
! Adding up the pairwise interaction matrix of the current snapshot
! 
           if ( dopim ) then
             call system_clock(t1pim)     
!
!~              call build_pim(ngrps,grps,nsubg,subg,sys%nat,atms,thr,   &
!~                             nnode,mol,agg,msize,nagg,xtcf%NATOMS,   &
!~                             xtcf%pos,(/xtcf%box(1,1),xtcf%box(2,2),    &
!~                             xtcf%box(3,3)/),pim,debug)
!
             call system_clock(t2pim)     
!
             tpim = tpim + dble(t2pim-t1pim)/dble(count_rate)   
           end if      
!
! Analyzing aggregates by their connectivity
!
!~            call analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,    &
!~                             sys%nat,atms,thr,nnode,adj,mol,agg,     & 
!~                             msize,nagg,xtcf%NATOMS,xtcf%pos,          &
!~                             (/xtcf%box(1,1),xtcf%box(2,2),             &
!~                               xtcf%box(3,3)/),posi,sys%mass,          &
!~                             sys%atname,xtcf%STEP,outp,debug)
!
           call system_clock(t1scrn)
!
           oldadj(:,:) = adj(:,:)
           adj(:,:)    = newadj(:,:)
!
           oldbox(:) = box(:)
           box(:)    = newbox(:)
!
           oldstep = actstep
           actstep = newstep
!
           posi(:,:) = newposi(:,:)
!
           call system_clock(t2scrn)
!
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
           call system_clock(t1read)
!
         end if
!
! Reading information of the next snapshot
!
         call xtcf%read
!
       end do
!
! Deallocate memory
!
       deallocate(newposi)
       deallocate(posi)
!
       deallocate(adj,oldadj,newadj)
!
       deallocate(mol,tag,agg)
       deallocate(nmol,imol,nagg,iagg)
!
       deallocate(oldmol,oldtag,oldagg)
       deallocate(oldnmol,oldimol,oldnagg,oldiagg)
!
       deallocate(newmol,newtag,newagg)
       deallocate(newnmol,newimol,newnagg,newiagg)
!
       deallocate(iam,iwas,iwill)
       deallocate(imnot,iwasnt,iwont)
       deallocate(oldiam,oldnot)
!
       deallocate(wasmap,willmap)
!
! Saving timings
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       return
       end subroutine aggscrnlife
!
!======================================================================!
!
! PRINTPOP - PRINT POPulations
!
! This subroutine 
!
       subroutine printpop(step,nsize,msize,pop,conc,frac,cin,volu,    &
                           box,nnode,nagg,magg,nsolv,iuni)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(3),intent(in)         ::  box    !  Simulation box dimensions
       real(kind=8),dimension(msize),intent(inout)  ::  pop    !  Populations
       real(kind=8),dimension(msize),intent(inout)  ::  conc   !  Concentrations
       real(kind=8),dimension(msize),intent(inout)  ::  frac   !  Molar fractions
       real(kind=8),intent(inout)                   ::  cin    !  Stechiometric concentration
       real(kind=8),intent(inout)                   ::  volu   !  Simulation box volume
       integer,dimension(nnode),intent(in)          ::  nagg   !
       integer,intent(in)                           ::  magg   !
       integer,intent(in)                           ::  nsolv  !
       integer,intent(in)                           ::  step   !
       integer,intent(in)                           ::  nnode  !  Total number of molecules
       integer,intent(in)                           ::  msize  !  Maximum aggregate size   
       integer,intent(in)                           ::  nsize  !  Actual maximum aggregate size
       integer,intent(in)                           ::  iuni   !
!
! Local variables
!
       real(kind=8),parameter                       ::  Na = 6.022140760E+23
       integer                                      ::  i      !  Index 
!
! Accumulating properties
!
       do i = 1, msize
         pop(i)  = pop(i)  + real(nagg(i))/magg*100
         frac(i) = frac(i) + real(nagg(i))/(magg+nsolv)
         conc(i) = conc(i) + real(nagg(i))/box(1)**3 
       end do     
!
       if ( nsize .gt. msize ) then
         do i = msize+1, nsize
           pop(msize)  = pop(msize)  + real(nagg(i))/magg*100
           frac(msize) = frac(msize) + real(nagg(i))/(magg+nsolv)
           conc(msize) = conc(msize) + real(nagg(i))/box(1)**3  
         end do
       end if                   
!
       cin = cin + real(nnode)/box(1)**3
!
       volu = volu + box(1)**3
!
! Printing populations of the current configuration
!
       write(iuni+1,'(I10,100(X,F12.8))') step,                        &
                                             real(nagg(:msize))/magg*100
       write(iuni+2,'(I10,100(X,F12.10))') step,                       &
                                         real(nagg(:msize))/(magg+nsolv)
       write(iuni+3,'(I10,100(X,F12.10))') step,                       &
                               real(nagg(:msize))/box(1)**3/(Na*1.0E-24)                            
!
       return
       end subroutine printpop
!
!======================================================================!
!
       end module aggtools
!
!======================================================================!
