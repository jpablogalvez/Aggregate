!======================================================================!
!
       module aggtools
!
       implicit none
!
       private
       public   ::  driver
!
       contains
!
!======================================================================!
!
! DRIVER - DRIVER selection algorithm
!
!  This subroutine defines the execution architecture of the algorithm
!
       subroutine driver(xtcf,ntype,rep,nnode,inode,nat,iat,natms,     &
                         iatms,ngrps,igrps,mnode,mat,matms,mgrps,thr,  &
                         thr2,neidis,msize,nmax,nsteps,nprint,minstep, &
                         maxstep,nsolv,avlife,nlife,dopim,schm,scrn,   &
                         doscrn,dolife,debug)
!
       use xdr,         only:  xtcfile
       use datatypes,   only:  repre
!
       use lengths,     only:  lenschm
!
       use properties
!
       use graphtools,  only:  buildadjmolbub,buildadjmolang
       use screening,   only:  scrnint,scrnosc,scrncol
!
       implicit none
!
! System information
!
       type(xtcfile),intent(inout)                              ::  xtcf      !  Trajectory information
       type(repre),dimension(ntype),intent(inout)               ::  rep       !
       integer,dimension(ntype),intent(in)                      ::  nnode     !  Total number of molecules of each size
       integer,dimension(ntype),intent(in)                      ::  inode     !  Initial molecule of each size
       integer,dimension(ntype),intent(in)                      ::  nat       !  Number of atoms in each molecule type
       integer,dimension(ntype),intent(in)                      ::  iat       !  Initial atoms of each molecule type
       integer,dimension(ntype),intent(in)                      ::  natms     !
       integer,dimension(ntype),intent(in)                      ::  iatms     !
       integer,dimension(ntype),intent(in)                      ::  ngrps     !  Number of active atoms in each representation type
       integer,dimension(ntype),intent(in)                      ::  igrps     !  Initial active atom of each representation type
       integer,intent(in)                                       ::  ntype     !  Number of molecules types
       integer,intent(in)                                       ::  mnode     !  Total number of molecules
       integer,intent(in)                                       ::  mat       !  Total number of atoms in the molecules
       integer,intent(in)                                       ::  matms     !  Total number of subgroups in the molecules
       integer,intent(in)                                       ::  mgrps     !  
       integer,intent(in)                                       ::  msize     !  Maximum aggregate size
       integer,intent(in)                                       ::  nmax      !
!
! Interaction criteria information
! 
       real(kind=8),dimension(mat,mat),intent(in)               ::  thr       !  Distance threshold
       real(kind=8),dimension(mat,mat),intent(in)               ::  thr2      !  Distance threshold
       real(kind=8),intent(in)                                  ::  neidis    !  Screening distance
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(nmax),intent(out)                ::  avlife    !
       integer,dimension(nmax),intent(out)                     ::  nlife     !
!
! Trajectory control variables
!     
       integer,intent(in)                                       ::  nprint    !  Populations printing interval
       integer,intent(in)                                       ::  minstep   !  First step for analysis
       integer,intent(in)                                       ::  maxstep   !  Last step for analysis
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
!
! Program control flags
!
       character(len=lenschm),intent(inout)                     ::  schm      !  Calculation scheme flag
       character(len=lenschm),intent(in)                        ::  scrn      !  Calculation screening flag
       logical,intent(in)                                       ::  doscrn    !  Screening calculation flag
       logical,intent(in)                                       ::  dolife    !  Lifetimes calculation flag
       logical,intent(in)                                       ::  dopim     !  PIM calculation flag
       logical,intent(in)                                       ::  debug     !  Debug mode
!
! AnalysisPhenolMD variables
!     
       integer,intent(in)                                       ::  nsolv     !
!
! Local variables
!
       character(len=lenschm)                                   ::  ch        !
!
! Declaration of procedure pointers
!
       abstract interface
!
       subroutine subadj(nnode,adj,neidis,msubg,mgrps,nat,thr,ngrps,   &
                         igrps,natms,posi,box)
!
       logical,dimension(nnode,nnode),intent(out)    ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(in)    ::  posi    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)          ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=8),dimension(nat,nat),intent(in)    ::  thr     !
       real(kind=8),intent(in)                       ::  neidis  !
       integer,dimension(nat),intent(in)             ::  ngrps   !  
       integer,dimension(nat),intent(in)             ::  igrps   ! 
       integer,intent(in)                            ::  msubg   !
       integer,intent(in)                            ::  mgrps   !  Number of subgroups
       integer,intent(in)                            ::  nnode   !  Number of residues
       integer,intent(in)                            ::  natms   ! 
       integer,intent(in)                            ::  nat     ! 
!
       end subroutine
!
       subroutine subscrn(nnode,oldadj,adj,newadj)
!
       logical,dimension(nnode,nnode),intent(inout)  ::  adj     !  Adjacency matrix of the current snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  oldadj  !  Adjacency matrix of the previous snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  newadj  !  Adjacency matrix of the next snapshot
       integer,intent(in)                            ::  nnode   !  Number of molecules
!
       end subroutine
!
       end interface
!
       procedure(subadj),pointer    ::  subbuildadj  => null()
       procedure(subscrn),pointer   ::  subscrnint   => null()
!~        procedure(subnadj),pointer   ::  subnbuildadj => null()
!~        procedure(subnscrn),pointer  ::  subnscrnint  => null()
!
! Selection of the aggregation algorithm ! FIXME: count reading time in this section
! --------------------------------------
! 
       ch = 'original'
       if ( doscrn ) ch = 'scrn'   
       if ( dolife ) ch = 'life'   
       if ( doscrn .and. dolife ) ch = 'scrnlife' 
!
       if ( ntype .eq. 1 ) then
!
! Homogeneous systems algorithm
!
         if ( trim(scrn) .eq. 'complete' ) then
           subscrnint => scrnint
         else if ( trim(scrn) .eq. 'collisions' ) then
           subscrnint => scrncol
         else if ( trim(scrn) .eq. 'oscillations' ) then
           subscrnint => scrnosc
         end if
!
         if ( trim(schm) .eq. 'distances' ) then
           subbuildadj => buildadjmolbub
         else if ( trim(schm) .eq. 'angles' ) then
           subbuildadj => buildadjmolang
         end if  
!
         select case (trim(ch))
           case ('original')
!
             call aggdist(xtcf,rep(1)%nat,nnode(1),natms(1),thr,thr2,  &
                          neidis,msize,nsteps,rep(1)%nbody,            &
                          rep(1)%ngrps,rep(1)%nsubg,rep(1)%ibody,      &
                          rep(1)%igrps,rep(1)%isubg,rep(1)%body,       &
                          rep(1)%grps,rep(1)%subg,rep(1)%atms,         &
                          rep(1)%mbody,rep(1)%mgrps,rep(1)%msubg,      &
                          rep(1)%matms,nprint,minstep,maxstep,nsolv,   &
                          dopim,subbuildadj,debug)
!
           case ('life')
!
             call agglife(xtcf,rep(1)%nat,nnode(1),natms(1),thr,thr2,  &
                          neidis,msize,nsteps,rep(1)%nbody,            &
                          rep(1)%ngrps,rep(1)%nsubg,rep(1)%ibody,      &
                          rep(1)%igrps,rep(1)%isubg,rep(1)%body,       &
                          rep(1)%grps,rep(1)%subg,rep(1)%atms,         &
                          rep(1)%mbody,rep(1)%mgrps,rep(1)%msubg,      &
                          rep(1)%matms,nprint,minstep,maxstep,nsolv,   &
                          avlife,nlife,dopim,subbuildadj,debug)
!
           case ('scrn')
!
             call aggscrn(xtcf,rep(1)%nat,nnode(1),natms(1),thr,thr2,  &
                          neidis,msize,nsteps,rep(1)%nbody,            &
                          rep(1)%ngrps,rep(1)%nsubg,rep(1)%ibody,      &
                          rep(1)%igrps,rep(1)%isubg,rep(1)%body,       &
                          rep(1)%grps,rep(1)%subg,rep(1)%atms,         &
                          rep(1)%mbody,rep(1)%mgrps,rep(1)%msubg,      &
                          rep(1)%matms,nprint,minstep,maxstep,nsolv,   &
                          dopim,subbuildadj,subscrnint,debug)
!
           case ('scrnlife')
!
             call aggscrnlife(xtcf,rep(1)%nat,nnode(1),natms(1),thr,   &
                              thr2,neidis,msize,nsteps,rep(1)%nbody,   &
                              rep(1)%ngrps,rep(1)%nsubg,rep(1)%ibody,  &
                              rep(1)%igrps,rep(1)%isubg,rep(1)%body,   &
                              rep(1)%grps,rep(1)%subg,rep(1)%atms,     &
                              rep(1)%mbody,rep(1)%mgrps,rep(1)%msubg,  &
                              rep(1)%matms,nprint,minstep,maxstep,     &
                              nsolv,avlife,nlife,dopim,subbuildadj,    &
                              subscrnint,debug)
!
         end select             
!
       else
!
! N-components systems algorithm
!
stop 'N-components algorithm not yet implemented!'
!~          if ( trim(scrn) .eq. 'complete' ) then
!~            subnscrnint => nscrnint
!~          else if ( trim(scrn) .eq. 'collisions' ) then
!~            subnscrnint => nscrncol
!~          else if ( trim(scrn) .eq. 'oscillations' ) then
!~            subnscrnint => nscrnosc
!~          end if
!~ !
!~          if ( trim(schm) .eq. 'distances' ) then
!~            subnbuildadj => nbuildadjmolbub
!~          else if ( trim(schm) .eq. 'angles' ) then
!~            subnbuildadj => nbuildadjmolang
!~          end if
!
!~          select case (trim(ch))
!~            case ('original')
!~ !
!~              call aggdist(xtcf,nat,nnode,natms,thr,thr2,neidis,msize,    &
!~                           nsteps,nbody,ngrps,nsubg,ibody,igrps,isubg,    &
!~                           body,grps,subg,atms,mbody,mgrps,msubg,matms,   &
!~                           nprint,minstep,maxstep,nsolv,dopim,            &
!~                           subbuildadj,debug)
!~ !
!~            case ('life')
!~ !
!~              call agglife(xtcf,nat,nnode,natms,thr,thr2,neidis,msize,    &
!~                           nsteps,nbody,ngrps,nsubg,ibody,igrps,isubg,    &
!~                           body,grps,subg,atms,mbody,mgrps,msubg,matms,   &
!~                           nprint,minstep,maxstep,nsolv,avlife,nlife,     &
!~                           dopim,subbuildadj,debug)
!~ !
!~            case ('scrn')
!~ !
!~              call aggscrn(xtcf,nat,nnode,natms,thr,thr2,neidis,msize,    &
!~                           nsteps,nbody,ngrps,nsubg,ibody,igrps,isubg,    &
!~                           body,grps,subg,atms,mbody,mgrps,msubg,matms,   &
!~                           nprint,minstep,maxstep,nsolv,dopim,            &
!~                           subbuildadj,subscrnint,debug)
!~ !
!~            case ('scrnlife')
!~ !
!~              call aggscrnlife(xtcf,nat,nnode,natms,thr,thr2,neidis,      &
!~                               msize,nsteps,nbody,ngrps,nsubg,ibody,      &
!~                               igrps,isubg,body,grps,subg,atms,mbody,     &
!~                               mgrps,msubg,matms,nprint,minstep,maxstep,  &
!~                               nsolv,avlife,nlife,dopim,subbuildadj,      &
!~                               subscrnint,debug)
!~ !
!~          end select     
!
       end if
!
       return
       end subroutine driver
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
! The subroutine SCREENING is employed to screen the interactions 
!  between the molecules. To do so a three-steps procedure is carried 
!  out with the following phylosophy: When an interaction exists, if I 
!  was or I will then I am, i.e., each interaction is anhilated if it is
!  not present in two consequtive configurations (removing collisions), 
!  and when an interaction does not exist, then if it was and it will 
!  then it is, i.e., one interaction can be created if it is present in 
!  the former and the following configurations (adding oscillations).
!
       subroutine aggscrn(xtcf,nat,nnode,natms,thr,thr2,neidis,msize,  &
                          nsteps,nbody,ngrps,nsubg,ibody,igrps,isubg,  &
                          body,grps,subg,atms,mbody,mgrps,msubg,matms, &
                          nprint,minstep,maxstep,nsolv,dopim,          &
                          buildadjmol,screen,debug)
!
       use xdr,        only:  xtcfile
!
       use properties
!
       use omp_var
       use datatypes
       use timings
!
       use printings
       use utils
!
       use screening
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
! External functions
!
       external                                                 ::  buildadjmol
       external                                                 ::  screen
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
       prob(:) = 0.0d0
       num(:)  = 0
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
           call screen(nnode,oldadj,adj,newadj)
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
           call printpop(actstep,nsize,msize,box,nnode,nagg,magg,      &
                         nsolv,uniout)
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
! Analyzing aggregates by their connectivity ! TODO: needs corrected adjacency matrix of aggs
!



!
! Keeping track the adjacency matrix information
!
           call cpu_time(tinscrn)
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
       subroutine aggdist(xtcf,nat,nnode,natms,thr,thr2,neidis,msize,  &
                          nsteps,nbody,ngrps,nsubg,ibody,igrps,isubg,  &
                          body,grps,subg,atms,mbody,mgrps,msubg,matms, &
                          nprint,minstep,maxstep,nsolv,dopim,          &
                          buildadjmol,debug)
!
       use xdr,       only:  xtcfile
!
       use properties
!
       use filenames, only:  outp,weight
       use datatypes
       use timings
!
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
! External functions
!
       external                                                 ::  buildadjmol
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
       prob(:) = 0.0d0
       num(:)  = 0.0d0
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
           call loopdist(xtcf,nat,nnode,natms,thr,thr2,neidis,msize,   &
                         nbody,ngrps,nsubg,ibody,igrps,isubg,body,     &
                         grps,subg,atms,mbody,mgrps,msubg,matms,nsolv, &
                         posi,box,adj,mol,tag,agg,nagg,iagg,nmol,imol, &
                         dopim,buildadjmol,debug)
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
       subroutine loopdist(xtcf,nat,nnode,natms,thr,thr2,neidis,msize, &
                           nbody,ngrps,nsubg,ibody,igrps,isubg,body,   &
                           grps,subg,atms,mbody,mgrps,msubg,matms,     &
                           nsolv,posi,box,adj,mol,tag,agg,nagg,iagg,   &
                           nmol,imol,dopim,buildadjmol,debug)
!
       use xdr,        only:  xtcfile
!
       use systeminf,  only:  sys
       use filenames,  only:  outp
!
       use datatypes
       use timings
!
       use printings
       use utils
!
       use thresholds,  only:  thrang,neiang
       use parameters,  only:  zero
!
       use geometry,    only:  sminimgvec
       use graphtools,  only:  chkangle
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
! External functions
!
       external                                                   ::  buildadjmol  
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
! Block-diagonalizing the adjacency matrix
!
       call blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg,nmol,      &
                      imol,magg)
!
! Printing the population of every aggregate
!
       call system_clock(t1read)
!         
       call printpop(xtcf%STEP,nsize,msize,box,nnode,nagg,magg,nsolv,  &
                     uniout)
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
! The subroutine SCREENING is employed to screen the interactions 
!  between the molecules. To do so a three-steps procedure is carried 
!  out with the following phylosophy: When an interaction exists, if I 
!  was or I will then I am, i.e., each interaction is anhilated if it is
!  not present in two consequtive configurations (removing collisions), 
!  and when an interaction does not exist, then if it was and it will 
!  then it is, i.e., one interaction can be created if it is present in 
!  the former and the following configurations (adding oscillations).
! The lifetimes are calculated keeping track of the aggregates which
!  are present in the former and the previous configurations.
!
       subroutine aggscrnlife(xtcf,nat,nnode,natms,thr,thr2,neidis,    &
                              msize,nsteps,nbody,ngrps,nsubg,ibody,    &
                              igrps,isubg,body,grps,subg,atms,mbody,   &
                              mgrps,msubg,matms,nprint,minstep,        &
                              maxstep,nsolv,avlife,nlife,dopim,        &
                              buildadjmol,screen,debug)
!
       use xdr,  only:  xtcfile
!
       use properties
!
       use omp_var
       use datatypes
       use timings
!
       use printings
       use utils
!
       use screening
       use lifetimes
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
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
       integer,intent(in)                                       ::  msize     !  Maximum aggregate size
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(nnode),intent(out)                ::  avlife    !
       integer,dimension(nnode),intent(out)                     ::  nlife     !
       integer,dimension(:),allocatable                         ::  life      !
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
! External functions
!
       external                                                 ::  buildadjmol  
       external                                                 ::  screen
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable                       ::  adj       !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  newadj    !  Adjacency matrix
       logical,dimension(:,:),allocatable                       ::  nextadj   !  Adjacency matrix
       integer,dimension(:),allocatable                         ::  mol       !  Molecules identifier
       integer,dimension(:),allocatable                         ::  tag       !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  agg       !  Aggregates size 
       integer,dimension(:),allocatable                         ::  newmol    !  Molecules identifier
       integer,dimension(:),allocatable                         ::  newtag    !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  newagg    !  Aggregates size 
       integer,dimension(:),allocatable                         ::  nagg      !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  iagg      !  
       integer,dimension(:),allocatable                         ::  nmol      !  
       integer,dimension(:),allocatable                         ::  imol      !  
       integer,dimension(:),allocatable                         ::  newnagg   !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  newiagg   !  
       integer,dimension(:),allocatable                         ::  newnmol   !  
       integer,dimension(:),allocatable                         ::  newimol   ! 
       integer,dimension(:),allocatable                         ::  willmap   !
       integer,intent(in)                                       ::  nnode     !  Total number of molecules
       integer                                                  ::  nsize     !  Actual maximum aggregate size
       integer                                                  ::  newsize   !  New maximum aggregate size
       integer                                                  ::  magg      !  Actual number of chemical species
       integer                                                  ::  newmagg   !  New number of chemical species
       integer                                                  ::  actstep   !
       integer                                                  ::  newstep   !
       integer                                                  ::  nextstep  !
       logical,dimension(:),allocatable                         ::  iwill     !
       logical,dimension(:),allocatable                         ::  iwont     !
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
       allocate(newmol(nnode),newtag(nnode),newagg(nnode))
       allocate(newnmol(nnode),newimol(nnode))
       allocate(newnagg(nnode),newiagg(nnode))
!
       allocate(iwill(nnode),iwont(nnode))
!
       allocate(willmap(nnode))
!
       allocate(life(nnode),auxlife(nnode)) 
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
       nextbox  = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       nextstep = xtcf%STEP
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
       call blockdiag(nnode,nextadj,newmol,newtag,newagg,newsize,      &
                      newnagg,newiagg,newnmol,newimol,newmagg)
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
       call screen(nnode,nextadj,adj,newadj)
!
       call cpu_time(tfinscrn)
       call system_clock(t2scrn)
!
       tcpuscrn = tcpuscrn + tfinscrn - tinscrn
       tscrn    = tscrn    + dble(t2scrn-t1scrn)/dble(count_rate) 
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
!$omp parallel do shared(life,avlife,nlife)                            &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
       do i = 1, nnode
         life(i)   = 0
         avlife(i) = 0.0d0
         nlife(i)  = 0
       end do
! 
!$omp end parallel do 
!
       call tracklife(nnode,newsize,newmol,nsize,mol,newnagg,newiagg,  &
                      newnmol,newimol,nagg,iagg,imol,iwill,iwont,      &
                      willmap)
!
!$omp parallel do shared(life,iwill,nagg)                              &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
       do i = nagg(1)+1, magg
         if ( iwill(i) ) life(i) = 1      
       end do
!
!$omp end parallel do 
!
       call cpu_time(tfinlife)
       call system_clock(t2life)
!
       tcpulife = tcpulife + tfinlife - tinlife
       tlife    = tlife    + dble(t2life-t1life)/dble(count_rate) 
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
           tadj    = tadj    + dble(t2adj-t1adj)/dble(count_rate) 
!
! Screening interactions between the molecules
!
           call cpu_time(tinscrn)
           call system_clock(t1scrn)
!
           call screen(nnode,adj,newadj,nextadj)
!
           call cpu_time(tfinscrn)
           call system_clock(t2scrn)
!
           tcpuscrn = tcpuscrn + tfinscrn - tinscrn
           tscrn    = tscrn    + dble(t2scrn-t1scrn)/dble(count_rate) 
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
           call tracklife(nnode,newsize,newmol,nsize,mol,newnagg,      &
                          newiagg,newnmol,newimol,nagg,iagg,imol,      &
                          iwill,iwont,willmap)
! 
           call cpu_time(tfinlife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tfinlife - tinlife
           tlife    = tlife    + dble(t2life-t1life)/dble(count_rate) 
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call printpop(actstep,nsize,msize,box,nnode,nagg,magg,      &
                         nsolv,uniout)
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
! Adding up lifetimes of the aggregates
!
           call cpu_time(tinlife)
           call system_clock(t1life)                  
!
           call calclife(avlife,nlife,nnode,life,nsize,nagg,iagg,      &
                         magg,iwont)
!
! Analyzing aggregates by their connectivity
!



!
! Keeping track the adjacency matrix information
!
           box(:)    = newbox(:)
           newbox(:) = nextbox(:)
!
!$omp parallel do shared(adj,newadj,auxlife)                           &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, nnode
             adj(:,i)    = newadj(:,i)
             newadj(:,i) = nextadj(:,i)
!
             auxlife(i) = 0
           end do
!
!$omp end parallel do 
!
!$omp parallel do shared(nagg,willmap,life,auxlife)                    &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = nagg(1)+1, magg
             if ( willmap(i) .ne. 0 ) auxlife(willmap(i)) = life(i)
           end do
!
!$omp end parallel do 
!
!$omp parallel do shared(life,auxlife,mol,agg,tag,nagg,iagg,nmol,imol, &
!$omp                    newmol,newagg,newtag,newnagg,newiagg,newnmol, &
!$omp                    newimol)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, nnode
             life(i) = auxlife(i)
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
           magg    = newmagg
           nsize   = newsize
!
           actstep = newstep
           newstep = nextstep
!
!$omp parallel do shared(posi,newposi,nextposi)                        &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, natms
             posi(:,i)    = newposi(:,i)
             newposi(:,i) = nextposi(:,i)
           end do
!
!$omp end parallel do       
!
           call cpu_time(tfinlife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tfinlife - tinlife
           tlife    = tlife    + dble(t2life-t1life)/dble(count_rate) 
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
! Deallocating memory
!
       deallocate(life,auxlife)
!
       deallocate(newposi,nextposi)
       deallocate(posi)
!
       deallocate(adj,newadj,nextadj)
!
       deallocate(mol,tag,agg)
       deallocate(nmol,imol,nagg,iagg)
!
       deallocate(newmol,newtag,newagg)
       deallocate(newnmol,newimol,newnagg,newiagg)
!
       deallocate(iwill,iwont)
!
       deallocate(willmap)
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
! AGGLIFE - AGGregates LIFEtimes
!
! This subroutine obtains the size and the number of aggregates present
!  in the system according to a given distance criteria.
! The interaction between two molecules is stored in the corresponding
!  adjacency matrix in the molecule-based representation. After block-
!  diagonalization of this matrix, the size aggregates is obtained from
!  the dimension of each block and the number of aggregates of each size
!  is obtained from the number of blocks of each size.
! The lifetimes are calculated keeping track of the aggregates which
!  are present in the former and the previous configurations.
!
       subroutine agglife(xtcf,nat,nnode,natms,thr,thr2,neidis,msize,  &
                          nsteps,nbody,ngrps,nsubg,ibody,igrps,isubg,  &
                          body,grps,subg,atms,mbody,mgrps,msubg,matms, &
                          nprint,minstep,maxstep,nsolv,avlife,nlife,   &
                          dopim,buildadjmol,debug)
!
       use xdr,  only:  xtcfile
!
       use properties
! 
       use omp_var
       use datatypes
       use timings
!
       use printings
       use utils
!
       use lifetimes
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
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
       integer,intent(in)                                       ::  msize     !  Maximum aggregate size
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(nnode),intent(out)                ::  avlife    !
       integer,dimension(nnode),intent(out)                     ::  nlife     !
       integer,dimension(:),allocatable                         ::  life      !
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
! External functions
!
       external                                                 ::  buildadjmol  
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable                       ::  adj       !  Adjacency matrix
       integer,dimension(:),allocatable                         ::  mol       !  Molecules identifier
       integer,dimension(:),allocatable                         ::  tag       !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  agg       !  Aggregates size 
       integer,dimension(:),allocatable                         ::  newmol    !  Molecules identifier
       integer,dimension(:),allocatable                         ::  newtag    !  Aggregates identifier
       integer,dimension(:),allocatable                         ::  newagg    !  Aggregates size 
       integer,dimension(:),allocatable                         ::  nagg      !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  iagg      !  
       integer,dimension(:),allocatable                         ::  nmol      !  
       integer,dimension(:),allocatable                         ::  imol      !  
       integer,dimension(:),allocatable                         ::  newnagg   !  Number of aggregates of each size
       integer,dimension(:),allocatable                         ::  newiagg   !  
       integer,dimension(:),allocatable                         ::  newnmol   !  
       integer,dimension(:),allocatable                         ::  newimol   ! 
       integer,dimension(:),allocatable                         ::  willmap   !
       integer,intent(in)                                       ::  nnode     !  Total number of molecules
       integer                                                  ::  nsize     !  Actual maximum aggregate size
       integer                                                  ::  newsize   !  New maximum aggregate size
       integer                                                  ::  magg      !  Actual number of chemical species
       integer                                                  ::  newmagg   !  New number of chemical species
       integer                                                  ::  actstep   !
       integer                                                  ::  newstep   !
       logical,dimension(:),allocatable                         ::  iwill     !
       logical,dimension(:),allocatable                         ::  iwont     !
!
! Declaration of time control variables
!
       real(kind=8)                                             ::  tinadj    !  Initial CPU building time
       real(kind=8)                                             ::  tfinadj   !  Final CPU building time
       real(kind=8)                                             ::  tinlife   !  Initial CPU lifetimes time
       real(kind=8)                                             ::  tfinlife  !  Final CPU lifetimes time
       integer                                                  ::  t1read    !  Initial reading time
       integer                                                  ::  t2read    !  Final reading time
       integer                                                  ::  t1adj     !  Initial building time
       integer                                                  ::  t2adj     !  Final building time
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
       real(kind=4),dimension(3)                                ::  box       !
       real(kind=4),dimension(3)                                ::  newbox    !
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
!
       allocate(mol(nnode),tag(nnode),agg(nnode))
       allocate(nmol(nnode),imol(nnode))
       allocate(nagg(nnode),iagg(nnode))
!
       allocate(newmol(nnode),newtag(nnode),newagg(nnode))
       allocate(newnmol(nnode),newimol(nnode))
       allocate(newnagg(nnode),newiagg(nnode))
!
       allocate(iwill(nnode),iwont(nnode))
!
       allocate(willmap(nnode))
!
       allocate(life(nnode),auxlife(nnode)) 
!
! Allocating variables depending on topological information 
!
       allocate(newposi(3,natms))
       allocate(posi(3,natms))
!
! Reading the first configuration
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
       box     = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       actstep = xtcf%STEP
!
!$omp parallel do shared(life,avlife,nlife)                            &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
       do i = 1, nnode
         life(i)   = 0
         avlife(i) = 0.0d0
         nlife(i)  = 0
       end do
! 
!$omp end parallel do 
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
! Block-diagonalizing the adjacency matrix of the first configuration
!
       call blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg,nmol,      &
                      imol,magg)
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
           call buildadj(nnode,adj,natms,newposi,xtcf%NATOMS,xtcf%pos, &
                         nat,thr2,mgrps,ngrps,igrps,msubg,nsubg,isubg, &
                         atms,newbox,neidis,buildadjmol)
!
           call cpu_time(tfinadj)
           call system_clock(t2adj) 
!
           tcpuadj = tcpuadj + tfinadj - tinadj
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate) 
!
! Block-diagonalizing the interaction-corrected adjacency matrix
!
           call blockdiag(nnode,adj,newmol,newtag,newagg,newsize,      &
                          newnagg,newiagg,newnmol,newimol,newmagg)
!
! Finding aggregates present in the new and the actual configurations
!
           call cpu_time(tinlife)
           call system_clock(t1life)
! 
           call tracklife(nnode,newsize,newmol,nsize,mol,newnagg,      &
                          newiagg,newnmol,newimol,nagg,iagg,imol,      &
                          iwill,iwont,willmap)
!
           call cpu_time(tfinlife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tfinlife - tinlife
           tlife = tlife + dble(t2life-t1life)/dble(count_rate) 
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call printpop(actstep,nsize,msize,box,nnode,nagg,magg,      &
                         nsolv,uniout)
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
! Adding up lifetimes of the aggregates
!
           call cpu_time(tinlife)
           call system_clock(t1life)
! 
           call calclife(avlife,nlife,nnode,life,nsize,nagg,iagg,      &
                         magg,iwont)
!
! Analyzing aggregates by their connectivity
!



!
! Keeping track the adjacency matrix information
!
           box(:)    = newbox(:)
!
!$omp parallel do shared(auxlife)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, nnode
             auxlife(i) = 0
           end do
!
!$omp end parallel do 
!
!$omp parallel do shared(nagg,willmap,life,auxlife)                    &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = nagg(1)+1, magg
             if ( willmap(i) .ne. 0 ) auxlife(willmap(i)) = life(i)
           end do
!
!$omp end parallel do 
!
!$omp parallel do shared(life,auxlife,mol,agg,tag,nagg,iagg,nmol,imol, &
!$omp                    newmol,newagg,newtag,newnagg,newiagg,newnmol, &
!$omp                    newimol)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, nnode
             life(i) = auxlife(i)
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
           magg    = newmagg
           nsize   = newsize
!
           actstep = newstep
!
!$omp parallel do shared(posi,newposi)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, natms
             posi(:,i) = newposi(:,i)
           end do
!
!$omp end parallel do                      
!
           call cpu_time(tfinlife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tfinlife - tinlife
           tlife    = tlife    + dble(t2life-t1life)/dble(count_rate) 
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
! Deallocating memory
!
       deallocate(life,auxlife)
!
       deallocate(newposi)
       deallocate(posi)
!
       deallocate(adj)
!
       deallocate(mol,tag,agg)
       deallocate(nmol,imol,nagg,iagg)
!
       deallocate(newmol,newtag,newagg)
       deallocate(newnmol,newimol,newnagg,newiagg)
!
       deallocate(iwill,iwont)
!
       deallocate(willmap)
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       return
       end subroutine agglife
!
!======================================================================!
!
! PRINTPOP - PRINT POPulations
!
! This subroutine 
!
       subroutine printpop(step,nsize,msize,box,nnode,nagg,magg,       &
                           nsolv,iuni)
!
       use parameters
       use properties
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(3),intent(in)         ::  box    !  Simulation box dimensions
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
       integer                                      ::  i      !  Index 
!
! Accumulating properties
!
       do i = 1, msize-1
         num(i)  = num(i)  + nagg(i)
         pop(i)  = pop(i)  + real(nagg(i))/magg*100
         frac(i) = frac(i) + real(nagg(i))/(magg+nsolv)
         conc(i) = conc(i) + real(nagg(i))/box(1)**3 
         prob(i) = prob(i) + real(i*nagg(i))/nnode*100 
       end do     
!
       dp1 = 0.0d0
       dp2 = 0.0d0
!
       if ( nsize .ge. msize ) then
         do i = msize, nsize
           dp1 = dp1 + real(nagg(i))
           dp2 = dp2 + real(i*nagg(i))
         end do
         num(msize)  = num(msize)  + int(dp1)
         pop(msize)  = pop(msize)  + dp1/magg*100
         frac(msize) = frac(msize) + dp1/(magg+nsolv)
         conc(msize) = conc(msize) + dp1/box(1)**3             
         prob(msize) = prob(msize) + dp2/nnode*100             
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
                      real(nagg(:msize-1))/(magg+nsolv),dp1/(magg+nsolv)
       write(iuni+3,'(I10,100(X,F10.6))') step,                        &
                      real(nagg(:msize-1))/box(1)**3/(Na*1.0E-24),     &
                                              dp1/box(1)**3/(Na*1.0E-24)  
       write(iuni+4,'(I10,100(X,F10.6))') step,                        & 
                            real(nagg(:msize-1))/nnode*100,dp2/nnode*100 ! TODO: print correct probability  
       write(iuni+5,'(I10,100(X,I10))') step,nagg(:msize-1),int(dp1)  
!
       return
       end subroutine printpop
!
!======================================================================!
!
       end module aggtools
!
!======================================================================!
