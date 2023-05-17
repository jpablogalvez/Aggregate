!======================================================================!
!
       module aggtools
!
       implicit none
!
       private
       public   ::  aggdist,aggscrn,aggscrnlife
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
       subroutine aggscrn(xtcf,nat,nnode,natms,thr,thr2,neidis,pim,    &
                          msize,pop,conc,frac,cin,volu,nsteps,nbody,   &
                          ngrps,nsubg,ibody,igrps,isubg,body,grps,     &
                          subg,atms,mbody,mgrps,msubg,matms,nprint,    &
                          minstep,maxstep,nsolv,dopim,buildadjmol,debug)
!
       use xdr, only: xtcfile
!
       use omp_var
       use screening
       use datatypes
       use timings
       use printings
       use utils
!
       use omp_lib
!
       implicit none
!
! System information
!
       type(xtcfile),intent(inout)                              ::  xtcf      !  Trajectory information
!
! Interaction criteria information
! 
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr       !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr2      !  Distance threshold
       real(kind=8),intent(in)                                  ::  neidis    !  Screening distance
!
! Average properties
!
       real(kind=8),dimension(mgrps,mgrps,msize-1),intent(out)  ::  pim       !  Pairwise interaction matrix
       real(kind=8),dimension(msize),intent(out)                ::  pop       !  Populations
       real(kind=8),dimension(msize),intent(out)                ::  conc      !  Concentrations
       real(kind=8),dimension(msize),intent(out)                ::  frac      !  Molar fractions
       real(kind=8),intent(out)                                 ::  cin       !  Stechiometric concentration
       real(kind=8),intent(out)                                 ::  volu      !  Simulation box volume
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
       integer,intent(in)                                       ::  msize     !  Maximum aggregate size
!
! Trajectory control variables
!     
       integer,intent(in)                                       ::  nprint    !  Populations printing interval
       integer,intent(in)                                       ::  minstep   !  First step for analysis
       integer,intent(in)                                       ::  maxstep   !  Last step for analysis
!
! Topological representations information
!
       integer,dimension(nat),intent(in)                        ::  body      !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  nbody     !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  ibody     !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  grps      !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  ngrps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  igrps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  subg      !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  nsubg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  isubg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  atms      !  Atoms identifier
       integer,intent(in)                                       ::  nat       !  Monomer atoms
       integer,intent(in)                                       ::  natms     !  Total number of subgroups in the system
       integer,intent(in)                                       ::  mbody     !  Number of bodies
       integer,intent(in)                                       ::  mgrps     !  Number of groups
       integer,intent(in)                                       ::  msubg     !  Number of subgroups
       integer,intent(in)                                       ::  matms     !  Number of interacting atoms in the monomer
!
! Program control flags
!
       logical,intent(in)                                       ::  dopim     !  PIM calculation flag
       logical,intent(in)                                       ::  debug     !  Debug mode
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable                       ::  adj       !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  oldadj    !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  newadj    !  Adjacency matrix
       integer,dimension(:),allocatable                         ::  mol       !  Molecules identifier
       integer,dimension(:),allocatable                         ::  tag       !   Aggregates identifier
       integer,dimension(:),allocatable                         ::  agg       !  Aggregates size 
       integer,dimension(:),allocatable                         ::  nagg      !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  iagg      !  
       integer,dimension(:),allocatable                         ::  nmol      !  
       integer,dimension(:),allocatable                         ::  imol      !  
       integer,intent(in)                                       ::  nnode     !  Total number of molecules
       integer                                                  ::  nsize     !  Actual maximum aggregate size
       integer                                                  ::  magg      !  Actual number of chemical species
       integer                                                  ::  actstep   !
       integer                                                  ::  newstep   !
!
! External functions
!
       external                                                 ::  buildadjmol
!
! Declaration of time control variables
!
       real(kind=8)                                             ::  tinadj    !  Initial CPU building time
       real(kind=8)                                             ::  tfinadj   !  Final CPU building time
       real(kind=8)                                             ::  tinscrn   !  Initial CPU screening time
       real(kind=8)                                             ::  tfinscrn  !  Final CPU screening time
       integer                                                  ::  t1read    !  Initial reading time
       integer                                                  ::  t2read    !  Final reading time
       integer                                                  ::  t1adj     !  Initial building time
       integer                                                  ::  t2adj     !  Final building time
       integer                                                  ::  t1scrn    !  Initial screening time
       integer                                                  ::  t2scrn    !  Final screening time
       integer                                                  ::  t1pim     !  Initial PIM analysis time
       integer                                                  ::  t2pim     !  Final PIM analysis time
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
       call cpu_time(tinadj)
       call system_clock(t1adj) 
!
       call buildadj(nnode,oldadj,natms,posi,xtcf%NATOMS,xtcf%pos,nat, &
                     thr2,mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,    &
                     box,neidis,buildadjmol)
!
       call cpu_time(tfinadj)
       call system_clock(t2adj) 
!
       tcpuadj = tcpuadj + tfinadj - tinadj
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
       call cpu_time(tinadj)
       call system_clock(t1adj) 
!
       call buildadj(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,nat,    &
                     thr2,mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,    &
                     box,neidis,buildadjmol)
!
       call cpu_time(tfinadj)
       call system_clock(t2adj) 
!
       tcpuadj = tcpuadj + tfinadj - tinadj
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
           call cpu_time(tinadj)
           call system_clock(t1adj) 
!
           call buildadj(nnode,newadj,natms,newposi,xtcf%NATOMS,       &
                         xtcf%pos,nat,thr2,mgrps,ngrps,igrps,msubg,    &
                         nsubg,isubg,atms,newbox,neidis,buildadjmol)
!
           call cpu_time(tfinadj)
           call system_clock(t2adj) 
!
           tcpuadj = tcpuadj + tfinadj - tinadj
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Screening interactions between the molecules
!
           call cpu_time(tinscrn)
           call system_clock(t1scrn)
!
           call scrnint(nnode,oldadj,adj,newadj)
!
           call cpu_time(tfinscrn)
           call system_clock(t2scrn)
!
           tcpuscrn = tcpuscrn + tfinscrn - tinscrn
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
!~              call build_pim(ngrps,grps,nsubg,subg,nat,atms,thr,   &
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
!~                             nat,atms,thr,nnode,adj,mol,agg,     & 
!~                             msize,nagg,xtcf%NATOMS,xtcf%pos,          &
!~                             (/xtcf%box(1,1),xtcf%box(2,2),             &
!~                               xtcf%box(3,3)/),posi,sys%mass,          &
!~                             sys%atname,xtcf%STEP,outp,debug)
!
           call cpu_time(tinscrn)        ! TODO: parallelize this region
           call system_clock(t1scrn)
!
           box(:)    = newbox(:)
!
           actstep = newstep
!
!$omp parallel do shared(adj,oldadj,newadj)                            &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkscrn)
!
           do i = 1, nnode
             oldadj(:,i) = adj(:,i)
             adj(:,i)    = newadj(:,i)
           end do
!
!$omp end parallel do                   
!
!$omp parallel do shared(posi,newposi)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkscrn)
!
           do i = 1, natms
             posi(:,i) = newposi(:,i)
           end do
!
!$omp end parallel do                   
!
           call cpu_time(tfinscrn)
           call system_clock(t2scrn)
!
           tcpuscrn = tcpuscrn + tfinscrn - tinscrn
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
!  in the system according to a given interaction criteria depending on
!  the external subroutine BUILDADJMOL. This can be a distance criteria
!  or a distance criteria between atoms A-B combined with a threshold 
!  angle A-B-C where C is an atom covalently bonded to A or B
! The interaction between two molecules is stored in the corresponding
!  adjacency matrix in the molecule-based representation. After block-
!  diagonalization of this matrix, the size aggregates is obtained from
!  the dimension of each block and the number of aggregates of each size
!  is obtained from the number of blocks of each size.
!
       subroutine aggdist(xtcf,nat,nnode,natms,thr,thr2,neidis,pim,    &
                          msize,pop,conc,frac,cin,volu,nsteps,nbody,   &
                          ngrps,nsubg,ibody,igrps,isubg,body,grps,     &
                          subg,atms,mbody,mgrps,msubg,matms,nprint,    &
                          minstep,maxstep,nsolv,dopim,buildadjmol,debug)
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
       integer,intent(in)                                       ::  nnode    !  Total number of molecules
       integer                                                  ::  nsize    !  Actual maximum aggregate size
       integer                                                  ::  magg     !  Actual number of chemical species
!
! External functions
!
       external                                                 ::  buildadjmol
!
! Declaration of time control variables
!
       integer                                                  ::  t1read   !  Initial reading time
       integer                                                  ::  t2read   !  Final reading time
!
! AnalysisPhenolMD variables
!     
       integer,intent(in)                                       ::  nsolv    !
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable                  ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(3)                                ::  box      !
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
! Finding the aggregates present in the current configuration
!
           call loopdist(xtcf,nat,nnode,natms,thr,thr2,neidis,pim,     &
                         msize,pop,conc,frac,cin,volu,nbody,ngrps,     &
                         nsubg,ibody,igrps,isubg,body,grps,subg,atms,  &
                         mbody,mgrps,msubg,matms,nsolv,posi,box,adj,   &
                         mol,tag,agg,nagg,iagg,nmol,imol,dopim,        &
                         buildadjmol,debug)
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
! LOOPDIST - main statements in the LOOP of the DISTances algorithm
!
! This subroutine 
!
       subroutine loopdist(xtcf,nat,nnode,natms,thr,thr2,neidis,pim,   &
                           msize,pop,conc,frac,cin,volu,nbody,ngrps,   &
                           nsubg,ibody,igrps,isubg,body,grps,subg,     &
                           atms,mbody,mgrps,msubg,matms,nsolv,posi,    &
                           box,adj,mol,tag,agg,nagg,iagg,nmol,imol,    &
                           dopim,buildadjmol,debug)
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
       type(xtcfile),intent(inout)                                ::  xtcf     !  Trajectory information
       real(kind=4),dimension(3,natms),intent(out)                ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(3),intent(in)                       ::  box      !
!
! Interaction criteria information
! 
       real(kind=8),dimension(nat,nat),intent(in)                 ::  thr      !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(in)                 ::  thr2     !  Distance threshold
       real(kind=8),intent(in)                                    ::  neidis   !  Screening distance
!
! Average properties
!
       real(kind=8),dimension(mgrps,mgrps,msize-1),intent(inout)  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(msize),intent(inout)                ::  pop      !  Populations
       real(kind=8),dimension(msize),intent(inout)                ::  conc     !  Concentrations
       real(kind=8),dimension(msize),intent(inout)                ::  frac     !  Molar fractions
       real(kind=8),intent(inout)                                 ::  cin      !  Stechiometric concentration
       real(kind=8),intent(inout)                                 ::  volu     !  Simulation box volume
       integer,intent(in)                                         ::  msize    !  Maximum aggregate size
!
! Topological representations information
!
       integer,dimension(nat),intent(in)                          ::  body     !  Number of groups in each body
       integer,dimension(nat),intent(in)                          ::  nbody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                          ::  ibody    !  Number of groups in each body
       integer,dimension(nat),intent(in)                          ::  grps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                          ::  ngrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                          ::  igrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                          ::  subg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                          ::  nsubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                          ::  isubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                          ::  atms     !  Atoms identifier
       integer,intent(in)                                         ::  nat      !  Monomer atoms
       integer,intent(in)                                         ::  natms    !  Total number of subgroups in the system
       integer,intent(in)                                         ::  mbody    !  Number of bodies
       integer,intent(in)                                         ::  mgrps    !  Number of groups
       integer,intent(in)                                         ::  msubg    !  Number of subgroups
       integer,intent(in)                                         ::  matms    !  Number of interacting atoms in the monomer
!
! Program control flags
!
       logical,intent(in)                                         ::  dopim    !  PIM calculation flag
       logical,intent(in)                                         ::  debug    !  Debug mode
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(nnode,nnode),intent(out)                 ::  adj      !  Adjacency matrix
       integer,dimension(nnode),intent(out)                       ::  mol      !  Molecules identifier
       integer,dimension(nnode),intent(out)                       ::  tag      !  Aggregates identifier
       integer,dimension(nnode),intent(out)                       ::  agg      !  Aggregates size 
       integer,dimension(nnode),intent(out)                       ::  nagg     !  Number of aggregates of each size
       integer,dimension(nnode),intent(out)                       ::  iagg     !  
       integer,dimension(nnode),intent(out)                       ::  nmol     !  
       integer,dimension(nnode),intent(out)                       ::  imol     !  
       integer,intent(in)                                         ::  nnode    !  Total number of molecules
       integer                                                    ::  nsize    !  Actual maximum aggregate size
       integer                                                    ::  magg     !  Actual number of chemical species
!
! External functions
!
       external                                                   ::  buildadjmol  
!
! Declaration of time control variables
!
       real(kind=8)                                               ::  tinadj   !  Initial CPU building time
       real(kind=8)                                               ::  tfinadj  !  Final CPU building time
       integer                                                    ::  t1read   !  Initial reading time
       integer                                                    ::  t2read   !  Final reading time
       integer                                                    ::  t1adj    !  Initial building time
       integer                                                    ::  t2adj    !  Final building time
       integer                                                    ::  t1pim    !  Initial PIM analysis time
       integer                                                    ::  t2pim    !  Final PIM analysis time
!
! AnalysisPhenolMD variables
!     
       integer,intent(in)                                         ::  nsolv    !
!
! Building the adjacency matrix
!
       call cpu_time(tinadj)
       call system_clock(t1adj) 
!
       call buildadj(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,nat, &
                        thr2,mgrps,ngrps,igrps,msubg,nsubg,isubg,atms, &
                        box,neidis,buildadjmol)
!
       call cpu_time(tfinadj)
       call system_clock(t2adj) 
!
       tcpuadj = tcpuadj + tfinadj - tinadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Block-diagonalizing the adjacency matrix
!
       call blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg,nmol,      &
                      imol,magg)
!
! Printing the population of every aggregate
!
       call system_clock(t1read)
!         
       call printpop(xtcf%STEP,nsize,msize,pop,conc,frac,cin,volu,     &
                     box,nnode,nagg,magg,nsolv,uniout)
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 

!~        if ( debug ) then
!~          call print_coord(xtcf,sys,outp,msize,nagg,nnode,mol,agg)
!~        end if
!
! Adding up the pairwise interaction matrix of the current snapshot
! 
       if ( dopim ) then
         call system_clock(t1pim)     
!
!~          call build_pim(ngrps,grps,nsubg,subg,nat,atms,thr,   &
!~                         nnode,mol,agg,msize,nagg,xtcf%NATOMS,   &
!~                         xtcf%pos,box,pim,debug)
!
         call system_clock(t2pim)     
!
         tpim = tpim + dble(t2pim-t1pim)/dble(count_rate)   
       end if      
!
! Analyzing aggregates by their connectivity
!
!~        call analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,    &
!~                         nat,atms,thr,nnode,adj,mol,agg,     & 
!~                         msize,nagg,xtcf%NATOMS,xtcf%pos,          &
!~                         box,posi,sys%mass,          &
!~                         sys%atname,xtcf%STEP,outp,debug)
!
       return
       end subroutine loopdist
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
       subroutine aggscrnlife(xtcf,nat,nnode,natms,thr,thr2,neidis,    &
                              pim,msize,pop,conc,frac,cin,volu,nsteps, &
                              nbody,ngrps,nsubg,ibody,igrps,isubg,     &
                              body,grps,subg,atms,mbody,mgrps,msubg,   &
                              matms,nprint,minstep,maxstep,nsolv,      &
                              avlife,nlife,dopim,buildadjmol,debug)
!
       use xdr, only: xtcfile
!
       use omp_var
       use screening
       use datatypes
       use timings
       use printings
       use utils
!
       use omp_lib
!
       implicit none
!
! System information
!
       type(xtcfile),intent(inout)                              ::  xtcf      !  Trajectory information
!
! Interaction criteria information
! 
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr       !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(in)               ::  thr2      !  Distance threshold
       real(kind=8),intent(in)                                  ::  neidis    !  Screening distance
!
! Average properties
!
       real(kind=8),dimension(mgrps,mgrps,msize-1),intent(out)  ::  pim       !  Pairwise interaction matrix
       real(kind=8),dimension(msize),intent(out)                ::  pop       !  Populations
       real(kind=8),dimension(msize),intent(out)                ::  conc      !  Concentrations
       real(kind=8),dimension(msize),intent(out)                ::  frac      !  Molar fractions
       real(kind=8),intent(out)                                 ::  cin       !  Stechiometric concentration
       real(kind=8),intent(out)                                 ::  volu      !  Simulation box volume
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
       integer,intent(in)                                       ::  msize     !  Maximum aggregate size
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(nnode),intent(out)                ::  avlife    !
       integer,dimension(nnode),intent(out)                     ::  nlife     !
       integer,dimension(:),allocatable                         ::  life      !
       integer,dimension(:),allocatable                         ::  oldlife   !
       integer,dimension(:),allocatable                         ::  auxlife   !
!
! Trajectory control variables
!     
       integer,intent(in)                                       ::  nprint    !  Populations printing interval
       integer,intent(in)                                       ::  minstep   !  First step for analysis
       integer,intent(in)                                       ::  maxstep   !  Last step for analysis
!
! Topological representations information
!
       integer,dimension(nat),intent(in)                        ::  body      !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  nbody     !  Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  ibody     !   Number of groups in each body
       integer,dimension(nat),intent(in)                        ::  grps      !   Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  ngrps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  igrps     !  Number of subgroups in each group
       integer,dimension(nat),intent(in)                        ::  subg      !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  nsubg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  isubg     !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)                        ::  atms      !  Atoms identifier
       integer,intent(in)                                       ::  nat       !  Monomer atoms
       integer,intent(in)                                       ::  natms     !  Total number of subgroups in the system
       integer,intent(in)                                       ::  mbody     !  Number of bodies
       integer,intent(in)                                       ::  mgrps     !  Number of groups
       integer,intent(in)                                       ::  msubg     !  Number of subgroups
       integer,intent(in)                                       ::  matms     !  Number of interacting atoms in the monomer
!
! Program control flags
!
       logical,intent(in)                                       ::  dopim     !  PIM calculation flag
       logical,intent(in)                                       ::  debug     !  Debug mode
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable                       ::  adj       !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  newadj    !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  nextadj   !  Adjacency matrix
       integer,dimension(:),allocatable                         ::  mol       !  Molecules identifier
       integer,dimension(:),allocatable                         ::  tag       !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  agg       !  Aggregates size 
       integer,dimension(:),allocatable                         ::  oldmol    !  Molecules identifier
       integer,dimension(:),allocatable                         ::  oldtag    !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  oldagg    !  Aggregates size 
       integer,dimension(:),allocatable                         ::  newmol    !  Molecules identifier
       integer,dimension(:),allocatable                         ::  newtag    !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  newagg    !  Aggregates size 
       integer,dimension(:),allocatable                         ::  nagg      !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  iagg      !  
       integer,dimension(:),allocatable                         ::  nmol      !  
       integer,dimension(:),allocatable                         ::  imol      !  
       integer,dimension(:),allocatable                         ::  oldnagg   !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  oldiagg   !  Number of aggregates of lower size
       integer,dimension(:),allocatable                         ::  oldnmol   !  
       integer,dimension(:),allocatable                         ::  oldimol   ! 
       integer,dimension(:),allocatable                         ::  newnagg   !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  newiagg   !  
       integer,dimension(:),allocatable                         ::  newnmol   !  
       integer,dimension(:),allocatable                         ::  newimol   ! 
       integer,dimension(:),allocatable                         ::  wasmap    !
       integer,dimension(:),allocatable                         ::  willmap   !
       integer,intent(in)                                       ::  nnode     !  Total number of molecules
       integer                                                  ::  nsize     !  Actual maximum aggregate size
       integer                                                  ::  oldsize   !  Old maximum aggregate size
       integer                                                  ::  newsize   !  New maximum aggregate size
       integer                                                  ::  magg      !  Actual number of chemical species
       integer                                                  ::  newmagg   !  New number of chemical species
       integer                                                  ::  oldmagg   !  Old number of chemical species
       integer                                                  ::  oldstep   !
       integer                                                  ::  actstep   !
       integer                                                  ::  newstep   !
       integer                                                  ::  nextstep  !
       logical,dimension(:),allocatable                         ::  iwas      !
       logical,dimension(:),allocatable                         ::  iwill     !
       logical,dimension(:),allocatable                         ::  iwasnt    !
       logical,dimension(:),allocatable                         ::  iwont     !
       logical,dimension(:),allocatable                         ::  oldiwas   !
       logical,dimension(:),allocatable                         ::  oldiwill  !
!
! External functions
!
       external                                                 ::  buildadjmol  
!
! Declaration of time control variables
!
       real(kind=8)                                             ::  tinadj    !  Initial CPU building time
       real(kind=8)                                             ::  tfinadj   !  Final CPU building time
       real(kind=8)                                             ::  tinscrn   !  Initial CPU screening time
       real(kind=8)                                             ::  tfinscrn  !  Final CPU screening time
       real(kind=8)                                             ::  tinlife   !  Initial CPU lifetimes time
       real(kind=8)                                             ::  tfinlife  !  Final CPU lifetimes time
       integer                                                  ::  t1read    !  Initial reading time
       integer                                                  ::  t2read    !  Final reading time
       integer                                                  ::  t1adj     !  Initial building time
       integer                                                  ::  t2adj     !  Final building time
       integer                                                  ::  t1scrn    !  Initial screening time
       integer                                                  ::  t2scrn    !  Final screening time
       integer                                                  ::  t1life    !  Initial lifetimes time
       integer                                                  ::  t2life    !  Final lifetimes time
       integer                                                  ::  t1pim     !  Initial PIM analysis time
       integer                                                  ::  t2pim     !  Final PIM analysis time
!
! AnalysisPhenolMD variables
!     
       integer,intent(in)                                       ::  nsolv     !
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable                  ::  posi      !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable                  ::  newposi   !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable                  ::  nextposi  !  Auxiliary coordinates
       real(kind=4),dimension(3)                                ::  box       !
       real(kind=4),dimension(3)                                ::  newbox    !
       real(kind=4),dimension(3)                                ::  nextbox   !
       integer                                                  ::  i,j       !   
       integer                                                  ::  iiagg     !   
       integer                                                  ::  inagg     !   
       integer                                                  ::  isize     ! 
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
       allocate(newadj(nnode,nnode),nextadj(nnode,nnode))
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
       allocate(iwas(nnode),iwill(nnode))
       allocate(iwasnt(nnode),iwont(nnode))
       allocate(oldiwas(nnode),oldiwill(nnode))
!
       allocate(wasmap(nnode),willmap(nnode))
!
       allocate(life(nnode),auxlife(nnode),oldlife(nnode)) 
!
! Allocating variables depending on topological information 
!
       allocate(nextposi(3,natms))
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
       nextbox = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       oldstep = xtcf%STEP
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the first old-configuration
!
       call cpu_time(tinadj)
       call system_clock(t1adj) 
!
       call buildadj(nnode,nextadj,natms,posi,xtcf%NATOMS,xtcf%pos,    &
                     nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,isubg,     &
                     atms,nextbox,neidis,buildadjmol)
!
       call cpu_time(tfinadj)
       call system_clock(t2adj) 
!
       tcpuadj = tcpuadj + tfinadj - tinadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Block-diagonalizing the adjacency matrix of the old-configuration
!
       call blockdiag(nnode,nextadj,oldmol,oldtag,oldagg,oldsize,      &
                      oldnagg,oldiagg,oldnmol,oldimol,oldmagg)
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
       box(:)  = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       actstep = xtcf%STEP
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the first configuration
!
       call cpu_time(tinadj)
       call system_clock(t1adj) 
!
       call buildadj(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,        &
                     nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,       &
                     isubg,atms,box,neidis,buildadjmol)
!
       call cpu_time(tfinadj)
       call system_clock(t2adj) 
!
       tcpuadj = tcpuadj + tfinadj - tinadj
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
       newbox(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       newstep = xtcf%STEP
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the first new-configuration
! 
       call cpu_time(tinadj)
       call system_clock(t1adj) 
!
       call buildadj(nnode,newadj,natms,newposi,xtcf%NATOMS,xtcf%pos,  &
                     nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,isubg,     &
                     atms,newbox,neidis,buildadjmol)
!
       call cpu_time(tfinadj)
       call system_clock(t2adj) 
!
       tcpuadj = tcpuadj + tfinadj - tinadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Screening interactions between the molecules
!
       call cpu_time(tinscrn)
       call system_clock(t1scrn)
!
       call scrnint(nnode,nextadj,adj,newadj)
!
       call cpu_time(tfinscrn)
       call system_clock(t2scrn)
!
       tcpuscrn = tcpuscrn + tfinscrn - tinscrn
       tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
! Block-diagonalizing the first interaction-corrected adjacency matrix
! 
       call blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg,nmol,      &
                      imol,magg)
!
! Finding aggregates present in the old and the actual configurations
!
       call cpu_time(tinlife)
       call system_clock(t1life)
! 
       call scrnblock(nnode,oldsize,oldmol,nsize,mol,oldnagg,oldiagg,  &
                      oldnmol,oldimol,nagg,iagg,imol,iwas,iwasnt,wasmap)
!
       call cpu_time(tfinlife)
       call system_clock(t2life)
!
       tcpulife = tcpulife + tfinlife - tinlife
       tlife = tlife + dble(t2life-t1life)/dble(count_rate) 
!
! Reading the first next-configuration
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
       life(:)   = 0
       avlife(:) = 0.0d0
       nlife(:)  = 0
!
! Analyzing frames in the inverval [minstep,maxstep]
!
       do while ( (xtcf%STAT.eq.0) .and. (xtcf%STEP.le.maxstep) )
         if ( mod(xtcf%STEP-minstep,nprint) .eq. 0 ) then
!
           nextbox(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
           nextstep   = xtcf%STEP
! 
           nsteps = nsteps + 1
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Building adjacency matrix of the new-configuration
!
           call cpu_time(tinadj)
           call system_clock(t1adj) 
!
           call buildadj(nnode,nextadj,natms,nextposi,xtcf%NATOMS,     &
                         xtcf%pos,nat,thr2,mgrps,ngrps,igrps,msubg,    &
                         nsubg,isubg,atms,nextbox,neidis,buildadjmol)
!
           call cpu_time(tfinadj)
           call system_clock(t2adj) 
!
           tcpuadj = tcpuadj + tfinadj - tinadj
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Screening interactions between the molecules
!
           call cpu_time(tinscrn)
           call system_clock(t1scrn)
!
           call scrnint(nnode,adj,newadj,nextadj)
!
           call cpu_time(tfinscrn)
           call system_clock(t2scrn)
!
           tcpuscrn = tcpuscrn + tfinscrn - tinscrn
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
! Block-diagonalizing the interaction-corrected adjacency matrix
!
           call blockdiag(nnode,newadj,newmol,newtag,newagg,newsize,   &
                          newnagg,newiagg,newnmol,newimol,newmagg)
!
! Finding aggregates present in the new and the actual configurations
!
           call cpu_time(tinlife)
           call system_clock(t1life)
! 
           call scrnblock(nnode,newsize,newmol,nsize,mol,newnagg,      &
                          newiagg,newnmol,newimol,nagg,iagg,imol,      &
                          iwill,iwont,willmap)
! 
!$omp parallel do shared(life,iwas,nagg)                               &
!$omp             private(i)                                           &
!$omp             schedule(dynamic)
!
           do i = nagg(1)+1, magg
             if ( iwas(i) ) then
               life(i) = life(i) + 1
             else
               life(i) = 0
             end if
           end do
!
!$omp end parallel do                   
!
           call cpu_time(tfinlife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tfinlife - tinlife
           tlife = tlife + dble(t2life-t1life)/dble(count_rate) 
!
!
           if ( debug ) then
             write(*,'(1X,A)') 'FIRST SCREENING'
             write(*,'(1X,15("."))') 
             write(*,*)
             write(*,'(2X,A,X,I9)')      'Old step',oldstep
             write(*,'(2X,21("="))') 
             write(*,*)
             write(*,'(2X,A,X,I6)')     'Total number of entities : ', &
                                                                 oldmagg
             write(*,*)
             write(*,'(2X,A,20(X,I6))') 'Aggregates of each type  : ', &
                                                        oldnagg(oldsize)
             write(*,'(2X,A,20(X,I6))') '  Accumulation           : ', &
                                                        oldiagg(oldsize)
             write(*,*)
             write(*,'(2X,A,20(X,I6))') 'Molecules of each type   : ', &
                                                        oldnmol(oldsize)
             write(*,'(2X,A,20(X,I6))') '  Accumulation           : ', &
                                                        oldimol(oldsize)
             write(*,*)
!
             call print_info(oldnagg(1),nnode-oldnagg(1),              &
                             oldmol(oldnagg(1)+1:),                    &
                             oldtag(oldnagg(1)+1:),                    &    
                             oldagg(oldnagg(1)+1:),'mol','tag','agg')
!
             write(*,'(2X,A,X,I9)')      'Step number',actstep
             write(*,'(2X,21("="))') 
             write(*,*)
             write(*,'(2X,A,X,I6)')     'Total number of entities : ', &
                                                                    magg
             write(*,*)
             write(*,'(2X,A,20(X,I6))') 'Aggregates of each type  : ', & 
                                                             nagg(nsize)
             write(*,'(2X,A,20(X,I6))') '  Accumulation           : ', &
                                                             iagg(nsize)
             write(*,*)
             write(*,'(2X,A,20(X,I6))') 'Molecules of each type   : ', &
                                                             nmol(nsize)
             write(*,'(2X,A,20(X,I6))') '  Accumulation           : ', & 
                                                             imol(nsize)
             write(*,*)
!
             call print_info(nagg(1),nnode-nagg(1),                  &
                             mol(nagg(1)+1:),                        &
                             tag(nagg(1)+1:),                        &  
                             agg(nagg(1)+1:),'mol','tag','agg')
!
             write(*,'(2X,A,X,I9)')      'New step',newstep
             write(*,'(2X,21("="))') 
             write(*,*)
             write(*,'(2X,A,X,I6)')     'Total number of entities : ', &
                                                                 newmagg
             write(*,*)
             write(*,'(2X,A,20(X,I6))') 'Aggregates of each type  : ', &
                                                        newnagg(newsize)
             write(*,'(2X,A,20(X,I6))') '  Accumulation           : ', &
                                                        newiagg(newsize)
             write(*,*)
             write(*,'(2X,A,20(X,I6))') 'Molecules of each type   : ', &
                                                        newnmol(newsize)
             write(*,'(2X,A,20(X,I6))') '  Accumulation           : ', &
                                                        newimol(newsize)
             write(*,*)
!
             call print_info(newnagg(1),nnode-newnagg(1),              &
                             newmol(newnagg(1)+1:),                    &
                             newtag(newnagg(1)+1:),                    &  
                             newagg(newnagg(1)+1:),'mol','tag','agg')
!
             write(*,'(4X,A)') 'I was'
             do i = 2, nsize
               if ( (iagg(i)+nagg(i)) .lt. (iagg(i)+1) ) exit
               write(*,'(I3,2(1X,A,1X,I5),1X,A,20(5X,L1))') i,'from',iagg(i)+1,   &
               'to',iagg(i)+nagg(i),':',iwas(iagg(i)+1:iagg(i)+nagg(i))
               write(*,'(25X,20(1X,I5))') (wasmap(iagg(i)+1:iagg(i)+nagg(i))-oldnagg(1))
             end do
             write(*,*)
!
             write(*,'(4X,A)') 'I am'
             do i = 2, nsize
               if ( (iagg(i)+nagg(i)) .lt. (iagg(i)+1) ) exit
!~                write(*,'(I3,2(1X,A,1X,I5),1X,A,20(5X,L1))') i,'from',iagg(i)+1,   &
!                 'to',iagg(i)+nagg(i),':',iam(iagg(i)+1:iagg(i)+nagg(i))
               write(*,'(25X,20(1X,I5))') life(iagg(i)+1:iagg(i)+nagg(i))
             end do
             write(*,*)
!
             write(*,'(4X,A)') 'I will'
             do i = 2, nsize
               if ( (iagg(i)+nagg(i)) .lt. (iagg(i)+1) ) exit
               write(*,'(I3,2(1X,A,1X,I5),1X,A,20(5X,L1))') i,'from',iagg(i)+1,   &
              'to',iagg(i)+nagg(i),':',iwill(iagg(i)+1:iagg(i)+nagg(i))
               write(*,'(25X,20(1X,I5))') (willmap(iagg(i)+1:iagg(i)+nagg(i))-newnagg(1))
             end do
             write(*,*)
           end if
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
           call cpu_time(tinlife)
           call system_clock(t1life)
! 
           call lifetimes(avlife,nlife,nnode,life,nsize,nagg,iagg,     &
                          willmap,iwont)
!
           call cpu_time(tfinlife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tfinlife - tinlife
           tlife = tlife + dble(t2life-t1life)/dble(count_rate) 
!
           call cpu_time(tinscrn)
           call system_clock(t1scrn)     ! TODO: parallelize this region
!
           box(:)    = newbox(:)
           newbox(:) = nextbox(:)
!
           oldmagg = magg
           oldsize = nsize
!
!$omp parallel do default(shared)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkscrn)
!
           do i = 1, nnode
             adj(:,i)    = newadj(:,i)
             newadj(:,i) = nextadj(:,i)
!
             oldmol(i) = mol(i)
             oldagg(i) = agg(i)
             oldtag(i) = tag(i)
!
             oldnagg(i) = nagg(i)
             oldiagg(i) = iagg(i)
!
             oldnmol(i) = nmol(i)
             oldimol(i) = imol(i)
!
             oldiwas(i) = iwas(i)
             oldiwill(i) = iwill(i)
             oldlife(i) = life(i)
!
             wasmap(i) = 0
             iwas(i)   = .FALSE.
             iwasnt(i) = .TRUE.
!
!~            do i = nagg(1)+1, magg
!~              if ( willmap(i) .ne. 0 ) then
!~                wasmap(willmap(i)) = i
!~                iwas(willmap(i))   = .TRUE.
!~                iwasnt(willmap(i)) = .FALSE.
!~              end if
!~            end do
!~ !
             auxlife(i) = 0
           end do
!
!$omp end parallel do                   
!
!$omp parallel do default(shared)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkscrn)
!
           do i = nagg(1)+1, magg
             if ( willmap(i) .ne. 0 ) then
               wasmap(willmap(i)) = i
               auxlife(willmap(i)) = life(i)
               iwas(willmap(i))   = .TRUE.
               iwasnt(willmap(i)) = .FALSE.
             end if
           end do
!
!$omp end parallel do                   
!
             magg    = newmagg
             nsize   = newsize
!
             oldstep = actstep
             actstep = newstep
             newstep = nextstep
!
!$omp parallel do default(shared)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkscrn)
!
           do i = 1, nnode
             life(:) = auxlife(i)
!
             mol(i) = newmol(i)
             agg(i) = newagg(i)
             tag(i) = newtag(i)
!
             nagg(i) = newnagg(i)
             iagg(i) = newiagg(i)
!
             nmol(i) = newnmol(i)
             imol(i) = newimol(i)
           end do
!
!$omp end parallel do                   
!
!$omp parallel do shared(posi,newposi,nextposi)                        &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunkscrn)
!
           do i = 1, natms
             posi(:,i) = newposi(:,i)
             newposi(:,i) = nextposi(:,i)
           end do
!
!$omp end parallel do                   
!
           call cpu_time(tfinscrn)
           call system_clock(t2scrn)
!
           tcpuscrn = tcpuscrn + tfinscrn - tinscrn
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
! Averaging lifetimes of the last snapshot
!
!~        call system_clock(t2read)
!~ !
!~        tread = tread + dble(t2read-t1read)/dble(count_rate) 
!~ !
!~        call system_clock(t1life)
!
!~        do isize = 2, oldsize
!~          do inagg = 1, oldnagg(isize)
!~ !
!~            iiagg = oldiagg(isize) + inagg
!~ !
!~            if ( oldiwas(iiagg) .and. oldiwill(iiagg) ) then
!~              avlife(isize) = avlife(isize) + oldlife(iiagg)
!~              nlife(isize)  = nlife(isize) + 1
!~            end if
!~ !
!~          end do
!~        end do
!
!write(*,*) avlife(:3)
!write(*,*) nlife(:3)
!
!~        call system_clock(t2life)
!~ !
!~        tlife = tlife + dble(t2life-t1life)/dble(count_rate) 
!
!~        call system_clock(t1read)
!
! Deallocate memory
!
       deallocate(life,auxlife,oldlife)
!
       deallocate(newposi,nextposi)
       deallocate(posi)
!
       deallocate(adj,newadj,nextadj)
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
       deallocate(iwas,iwill)
       deallocate(iwasnt,iwont)
       deallocate(oldiwas,oldiwill)
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
       use parameters
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
       real(kind=8)                                 ::  dp1    !
       real(kind=8)                                 ::  dp2    !
       real(kind=8)                                 ::  dp3    !
       integer                                      ::  i      !  Index 
!
! Accumulating properties
!
       do i = 1, msize-1
         pop(i)  = pop(i)  + real(nagg(i))/magg*100
         frac(i) = frac(i) + real(nagg(i))/(magg+nsolv)
         conc(i) = conc(i) + real(nagg(i))/box(1)**3 
       end do     
!
       dp1 = 0.0d0
       dp2 = 0.0d0
       dp3 = 0.0d0
       if ( nsize .ge. msize ) then
         do i = msize, nsize
           dp1 = dp1 + real(nagg(i))
           dp2 = dp2 + real(nagg(i))
           dp3 = dp3 + real(nagg(i))
         end do
         pop(msize)  = pop(msize)  + dp1/magg*100
         frac(msize) = frac(msize) + dp2/(magg+nsolv)
         conc(msize) = conc(msize) + dp3/box(1)**3             
       end if                   
!
       cin = cin + real(nnode)/box(1)**3
!
       volu = volu + box(1)**3
!
! Printing populations of the current configuration
!
       write(iuni+1,'(I10,100(X,F10.6))') step,                        &
                              real(nagg(:msize-1))/magg*100,dp1/magg*100
       write(iuni+2,'(I10,100(X,F10.6))') step,                       &
                      real(nagg(:msize-1))/(magg+nsolv),dp2/(magg+nsolv)
       write(iuni+3,'(I10,100(X,F10.6))') step,                       &
                      real(nagg(:msize-1))/box(1)**3/(Na*1.0E-24),     &
                                              dp3/box(1)**3/(Na*1.0E-24)  
!
!~        do i = 1, msize
!~          pop(i)  = pop(i)  + real(nagg(i))/magg*100
!~          frac(i) = frac(i) + real(nagg(i))/(magg+nsolv)
!~          conc(i) = conc(i) + real(nagg(i))/box(1)**3 
!~        end do     
!~ !
!~        if ( nsize .gt. msize ) then ! FLAG: maybe here is a bug
!~          do i = msize+1, nsize
!~            pop(msize)  = pop(msize)  + real(nagg(i))/magg*100
!~            frac(msize) = frac(msize) + real(nagg(i))/(magg+nsolv)
!~            conc(msize) = conc(msize) + real(nagg(i))/box(1)**3  
!~          end do
!~        end if                   
!~ !
!~        cin = cin + real(nnode)/box(1)**3
!~ !
!~        volu = volu + box(1)**3
!~ !
!~ ! Printing populations of the current configuration
!~ !
!~        write(iuni+1,'(I10,100(X,F12.8))') step,                        &
!~                                              real(nagg(:msize))/magg*100
!~        write(iuni+2,'(I10,100(X,F12.10))') step,                       &
!~                                          real(nagg(:msize))/(magg+nsolv)
!~        write(iuni+3,'(I10,100(X,F12.10))') step,                       &
!~                                real(nagg(:msize))/box(1)**3/(Na*1.0E-24)                            
!~ !
       return
       end subroutine printpop
       !
!======================================================================!
!
! LIFETIMES - LIFETIMES calculation
!
! This subroutine 
!
       subroutine lifetimes(avlife,nlife,nnode,life,nsize,nagg,iagg,   &
                            willmap,iwont)
!
       use omp_lib
!
       implicit none
!
! Input/Output variables
!
       real(kind=8),dimension(nnode),intent(inout)  ::  avlife   !  Average lifetimes
       integer,dimension(nnode),intent(inout)       ::  nlife    !  
       integer,dimension(nnode),intent(inout)       ::  life     !
       integer,dimension(nnode),intent(in)          ::  willmap  !
       integer,dimension(nnode),intent(in)          ::  nagg     !
       integer,dimension(nnode),intent(in)          ::  iagg     !
       integer,intent(in)                           ::  nnode    !  Total number of molecules
       integer,intent(in)                           ::  nsize    ! 
       logical,dimension(nnode),intent(in)          ::  iwont    !
!
! Local variables
!
       integer,dimension(nnode)                     ::  auxlife
       integer                                      ::  iiagg    !   
       integer                                      ::  inagg    !   
       integer                                      ::  isize    !   
!
! Averaging lifetimes
!
!~        auxlife(:) = 0
!
!$omp parallel do shared(iwont,life,iagg,avlife,nlife,nagg)            &
!$omp             private(isize,inagg,iiagg)                           &
!$omp             schedule(dynamic,1)
!
       do isize = 2, nsize
         do inagg = 1, nagg(isize)
!
           iiagg = iagg(isize) + inagg
!
           if ( iwont(iiagg) ) then
             avlife(isize) = avlife(isize) + life(iiagg)
             nlife(isize)  = nlife(isize) + 1
!~            else
!~              auxlife(willmap(iiagg)) = life(iiagg)
           end if
!
         end do
       end do
!
!$omp end parallel do                   
!
!~ write(*,*) avlife(:3)
!~ write(*,*) nlife(:3)
!
!~        life(:) = auxlife(:)
!
       return
       end subroutine lifetimes
!
!======================================================================!
!
       end module aggtools
!
!======================================================================!
