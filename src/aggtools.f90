!======================================================================!
!
       module aggtools
!
       implicit none
!
       private
       public   ::  driver
!
       logical  ::  doiso_global = .FALSE.
!
       contains
!
!======================================================================!
!
! DRIVER - DRIVER selection algorithm
!
!  This subroutine defines the execution architecture of the algorithm
!
       subroutine driver(neidis,nsteps,nprint,minstep,maxstep,nsolv,   &
                         avlife,nlife,dopim,doconf,domon,schm,scrn,    &
                         cconf,doscrn,dolife,docoord,doiso,debug)
!
       use systeminf,   only:  mnode,rep,xtcf,mtype,nnode,natms
       use properties,  only:  nmax,msize,cin
!
       use lengths,     only:  lenschm
!
       use graphtools,  only:  buildadjmolbub,buildadjmolang,          &
                               nbuildadjmolbub,nbuildadjmolang,        &
                               nbuildadjgrpsbub,nbuildadjgrpsang,      &
                               nbuildadjbodybub,nbuildadjbodyang,      &
                               nbuildadjbodybubdir,nbuildadjbodyangdir,&
                               nbuildadjgrpsmonbub,nbuildadjgrpsmonang,&
                               nbuildadjbodymonbub,nbuildadjbodymonang,&
                               nbuildadjbodymonbubdir,                 &
                               nbuildadjbodymonangdir
       use screening,   only:  scrnint,scrnosc,scrncol
!
       implicit none
!
! Interaction criteria information
!
       real(kind=4),intent(in)                     ::  neidis   !  Screening distance
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(nmax),intent(out)    ::  avlife   !
       integer,dimension(nmax),intent(out)         ::  nlife    !
!
! Trajectory control variables
!
       integer,intent(in)                          ::  nprint   !  Populations printing interval
       integer,intent(in)                          ::  minstep  !  First step for analysis
       integer,intent(in)                          ::  maxstep  !  Last step for analysis
       integer,intent(out)                         ::  nsteps   !  Number of snapshots analyzed
!
! Program control flags
!
       character(len=lenschm),intent(inout)        ::  schm     !  Calculation scheme flag
       character(len=lenschm),intent(in)           ::  scrn     !  Calculation screening flag
       character(len=lenschm),intent(in)           ::  cconf    !  Calculation conformations flag
       logical,intent(in)                          ::  doscrn   !  Screening calculation flag
       logical,intent(in)                          ::  dolife   !  Lifetimes calculation flag
       logical,intent(in)                          ::  dopim    !  PIM calculation flag
       logical,intent(in)                          ::  doconf   !  Conformational analysis flag
       logical,intent(in)                          ::  domon    !  Monomer intramolecular edges flag
       logical,intent(in)                          ::  docoord  !  Coordinate printing flag
       logical,intent(in)                          ::  doiso    !  On-the-fly isomorphism classification flag
       logical,intent(in)                          ::  debug    !  Debug mode
!
! AnalysisPhenolMD variables
!
       integer,intent(in)                          ::  nsolv    !
!
! Local variables
!
       character(len=lenschm)                      ::  ch       !
!
! Declaration of procedure pointers
!
       abstract interface
!
       subroutine subadj(nnode,adj,neidis,msubg,mgrps,nat,ngrps,igrps, &
                         natms,posi,box)
!
       logical,dimension(nnode,nnode),intent(out)  ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(in)  ::  posi    !  Atomic coordinates
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box
       real(kind=4),intent(in)                     ::  neidis  !
       integer,dimension(nat),intent(in)           ::  ngrps   !
       integer,dimension(nat),intent(in)           ::  igrps   !
       integer,intent(in)                          ::  msubg   !
       integer,intent(in)                          ::  mgrps   !  Number of subgroups
       integer,intent(in)                          ::  nnode   !  Number of residues
       integer,intent(in)                          ::  natms   !
       integer,intent(in)                          ::  nat     !
!
       end subroutine subadj
!
       subroutine nsubadj(mnode,adj,matms,posi,neidis,box)
!
       logical,dimension(mnode,mnode),intent(out)  ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,matms),intent(in)  ::  posi    !  Atomic coordinates
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box
       real(kind=4),intent(in)                     ::  neidis  !
       integer,intent(in)                          ::  mnode   !
       integer,intent(in)                          ::  matms   !
!
       end subroutine nsubadj
!
       subroutine nsubadjrep(mnode,node,madj,adj,matms,posi,           &
                             box,mtype,nnode,inode,ibodymon)
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
       end subroutine nsubadjrep
!
       subroutine nsubprint(nmax,nagg,imol,mnode,node,matms,posi,      &
                            box,buildadjrep,buildadjmon,domon)
!
! Input/output variables
!
       integer,dimension(nmax),intent(in)          ::  nagg    !  Number of aggregates of each size
       integer,dimension(nmax),intent(in)          ::  imol    !
       integer,dimension(mnode),intent(in)         ::  node    !  Molecules identifier
       real(kind=4),dimension(3,matms),intent(in)  ::  posi    !
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       integer,intent(in)                          ::  nmax    !
       integer,intent(in)                          ::  mnode   !
       integer,intent(in)                          ::  matms   !
       logical,intent(in)                          ::  domon   !
!
! External functions
!
       external                                    ::  buildadjrep
       external                                    ::  buildadjmon
!
       end subroutine nsubprint
!
       subroutine subscrn(nnode,oldadj,adj,newadj)
!
       logical,dimension(nnode,nnode),intent(inout)  ::  adj     !  Adjacency matrix of the current snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  oldadj  !  Adjacency matrix of the previous snapshot
       logical,dimension(nnode,nnode),intent(in)     ::  newadj  !  Adjacency matrix of the next snapshot
       integer,intent(in)                            ::  nnode   !  Number of molecules
!
       end subroutine subscrn
!
       end interface
!
       procedure(subadj),pointer      ::  subbuildadj     => null()
       procedure(nsubadj),pointer     ::  subnbuildadj    => null()
       procedure(nsubadjrep),pointer  ::  subnbuildadjrep => null()
       procedure(nsubadjrep),pointer  ::  subnbuildadjmon => null()
       procedure(nsubprint),pointer   ::  subnprintadjrep => null()
       procedure(subscrn),pointer     ::  subscrnint      => null()
!
       doiso_global = doiso
!
! Selection of the aggregation algorithm ! FIXME: count reading time in this section
! --------------------------------------
!
       ch = 'original'
       if ( doscrn ) ch = 'scrn'
       if ( dolife ) ch = 'life'
       if ( doscrn .and. dolife ) ch = 'scrnlife'
!
       if ( trim(scrn) .eq. 'complete' ) then
         subscrnint => scrnint
       else if ( trim(scrn) .eq. 'collisions' ) then
         subscrnint => scrncol
       else if ( trim(scrn) .eq. 'oscillations' ) then
         subscrnint => scrnosc
       end if
!
! Homogeneous systems algorithm
! .............................
!
!!!       if ( mtype .eq. 1 ) then
!!!!
!!!         if ( trim(schm) .eq. 'distances' ) then
!!!           subbuildadj => buildadjmolbub
!!!         else if ( trim(schm) .eq. 'angles' ) then
!!!           subbuildadj => buildadjmolang
!!!         end if
!!!!
!!!         select case (trim(ch))
!!!           case ('original')
!!!!
!!!             call aggdist(xtcf,rep(1)%nat,nnode(1),natms(1),neidis,    &
!!!                          msize,nsteps,rep(1)%ngrps,rep(1)%nsubg,      &
!!!                          rep(1)%igrps,rep(1)%isubg,rep(1)%atms,       &
!!!                          rep(1)%mgrps,rep(1)%msubg,nprint,minstep,    &
!!!                          maxstep,nsolv,dopim,doconf,cin(1),           &
!!!                          subbuildadj,debug)
!!!!
!!!           case ('life')
!!!!
!!!             call agglife(xtcf,rep(1)%nat,nnode(1),natms(1),neidis,    &
!!!                          msize,nsteps,rep(1)%nbody,rep(1)%ngrps,      &
!!!                          rep(1)%nsubg,rep(1)%ibody,rep(1)%igrps,      &
!!!                          rep(1)%isubg,rep(1)%body,rep(1)%grps,        &
!!!                          rep(1)%subg,rep(1)%atms,rep(1)%mbody,        &
!!!                          rep(1)%mgrps,rep(1)%msubg,rep(1)%matms,      &
!!!                          nprint,minstep,maxstep,nsolv,avlife,nlife,   &
!!!                          dopim,cin(1),subbuildadj,debug)
!!!!
!!!           case ('scrn')
!!!!
!!!             call aggscrn(xtcf,rep(1)%nat,nnode(1),natms(1),neidis,    &
!!!                          msize,nsteps,rep(1)%nbody,rep(1)%ngrps,      &
!!!                          rep(1)%nsubg,rep(1)%ibody,rep(1)%igrps,      &
!!!                          rep(1)%isubg,rep(1)%body,rep(1)%grps,        &
!!!                          rep(1)%subg,rep(1)%atms,rep(1)%mbody,        &
!!!                          rep(1)%mgrps,rep(1)%msubg,rep(1)%matms,      &
!!!                          nprint,minstep,maxstep,nsolv,dopim,cin(1),   &
!!!                          subbuildadj,subscrnint,debug)
!!!!
!!!           case ('scrnlife')
!!!!
!!!             call aggscrnlife(xtcf,rep(1)%nat,nnode(1),natms(1),       &
!!!                              neidis,msize,nsteps,rep(1)%nbody,        &
!!!                              rep(1)%ngrps,rep(1)%nsubg,rep(1)%ibody,  &
!!!                              rep(1)%igrps,rep(1)%isubg,rep(1)%body,   &
!!!                              rep(1)%grps,rep(1)%subg,rep(1)%atms,     &
!!!                              rep(1)%mbody,rep(1)%mgrps,rep(1)%msubg,  &
!!!                              rep(1)%matms,nprint,minstep,maxstep,     &
!!!                              nsolv,avlife,nlife,dopim,cin(1),         &
!!!                              subbuildadj,subscrnint,debug)
!!!!
!!!         end select
!!!!
!!!! N-components systems algorithm
!!!! ..............................
!!!!
!!!       else
!
         if ( trim(schm) .eq. 'distances' ) then
           subnbuildadj     => nbuildadjmolbub
         else if ( trim(schm) .eq. 'angles' ) then
           subnbuildadj     => nbuildadjmolang
         end if
!
         if ( trim(cconf) .eq. 'body' ) then
           subnprintadjrep => printadjbody
           if ( trim(schm) .eq. 'distances' ) then
             subnbuildadjrep => nbuildadjbodybub
             subnbuildadjmon => nbuildadjbodymonbub
           else if ( trim(schm) .eq. 'angles' ) then
             subnbuildadjrep => nbuildadjbodyang
             subnbuildadjmon => nbuildadjbodymonang
           end if
         else if ( trim(cconf) .eq. 'grps' ) then
           subnprintadjrep => printadjgrps
           if ( trim(schm) .eq. 'distances' ) then
             subnbuildadjrep => nbuildadjgrpsbub
             subnbuildadjmon => nbuildadjgrpsmonbub
           else if ( trim(schm) .eq. 'angles' ) then
             subnbuildadjrep => nbuildadjgrpsang
             subnbuildadjmon => nbuildadjgrpsmonang
           end if
         else if ( trim(cconf) .eq. 'bodydir' ) then
           subnprintadjrep => printadjbodydir
           if ( trim(schm) .eq. 'distances' ) then
             subnbuildadjrep => nbuildadjbodybubdir
             subnbuildadjmon => nbuildadjbodymonbubdir
           else if ( trim(schm) .eq. 'angles' ) then
             subnbuildadjrep => nbuildadjbodyangdir
             subnbuildadjmon => nbuildadjbodymonangdir
           end if
         end if
!
         select case (trim(ch))
           case ('original')
!
             call naggdist(neidis,nsteps,nprint,minstep,maxstep,       &
                           nsolv,dopim,doconf,domon,subnbuildadj,      &
                           subnbuildadjrep,subnbuildadjmon,            &
                           subnprintadjrep,docoord,debug)
!
           case ('life')
!
             call nagglife(neidis,nsteps,nprint,minstep,maxstep,       &
                           nsolv,avlife,nlife,dopim,doconf,domon,      &
                           subnbuildadj,subnbuildadjrep,               &
                           subnbuildadjmon,subnprintadjrep,docoord,debug)
!
           case ('scrn')
!
             call naggscrn(neidis,nsteps,nprint,minstep,maxstep,nsolv, &
                           dopim,doconf,domon,subnbuildadj,            &
                           subnbuildadjrep,subnbuildadjmon,            &
                           subscrnint,scrn,cconf,docoord,debug)
!
           case ('scrnlife')
!
             call naggscrnlife(neidis,nsteps,nprint,minstep,maxstep,  &
                               nsolv,avlife,nlife,dopim,doconf,domon, &
                               subnbuildadj,subnbuildadjrep,          &
                               subnbuildadjmon,subscrnint,scrn,cconf, &
                               docoord,debug)
!
         end select
!
!!!       end if
!
       return
       end subroutine driver
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
       subroutine aggdist(xtcf,nat,nnode,natms,neidis,msize,nsteps,    &
                          ngrps,nsubg,igrps,isubg,atms,mgrps,msubg,    &
                          nprint,minstep,maxstep,nsolv,dopim,doconf,   &
                          cin,buildadjmol,debug)
!
       use xdr,         only:  xtcfile
!
       use properties,  only:  pim,num,pop,frac,conc,volu
!
       use units,       only:  uniout
!
       use timings,     only:  count_rate,tread,tadj,tpim,tconf,       &
                               tcpuadj,tcpupim,tcpuconf
!
       implicit none
!
! System information
!
       type(xtcfile),intent(inout)                 ::  xtcf     !  Trajectory information
!
! Interaction criteria information
!
       real(kind=4),intent(in)                     ::  neidis   !  Screening distance
!
! Average properties
!
       integer,intent(out)                         ::  nsteps   !  Number of snapshots analyzed
       real(kind=8),intent(out)                    ::  cin      !
       integer,intent(in)                          ::  msize    !  Maximum aggregate size
!
! Trajectory control variables
!
       integer,intent(in)                          ::  nprint   !  Populations printing interval
       integer,intent(in)                          ::  minstep  !  First step for analysis
       integer,intent(in)                          ::  maxstep  !  Last step for analysis
!
! Topological representations information
!
       integer,dimension(nat),intent(in)           ::  ngrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)           ::  igrps    !  Number of subgroups in each group
       integer,dimension(nat),intent(in)           ::  nsubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)           ::  isubg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(in)           ::  atms     !  Atoms identifier
       integer,intent(in)                          ::  nat      !  Monomer atoms
       integer,intent(in)                          ::  natms    !  Total number of subgroups in the system
       integer,intent(in)                          ::  mgrps    !  Number of groups
       integer,intent(in)                          ::  msubg    !  Number of subgroups
!
! Program control flags
!
       logical,intent(in)                          ::  dopim    !  PIM calculation flag
       logical,intent(in)                          ::  doconf   !  Conformational analysis flag
       logical,intent(in)                          ::  debug    !  Debug mode
!
! External functions
!
       external                                    ::  buildadjmol
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable          ::  adj      !  Adjacency matrix
       integer,dimension(:),allocatable            ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable            ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable            ::  agg      !  Aggregates size
       integer,dimension(:),allocatable            ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable            ::  iagg     !
       integer,dimension(:),allocatable            ::  nmol     !
       integer,dimension(:),allocatable            ::  imol     !
       integer,intent(in)                          ::  nnode    !  Total number of molecules
       integer                                     ::  nsize    !  Actual maximum aggregate size
       integer                                     ::  magg     !  Actual number of chemical species

!
! Declaration of time control variables
!
       real(kind=8)                                ::  tinadj   !  Initial CPU building time
       real(kind=8)                                ::  tfinadj  !  Final CPU building time
       real(kind=8)                                ::  tinpim   !  Initial CPU PIM time
       real(kind=8)                                ::  tfinpim  !  Final CPU PIM time
       real(kind=8)                                ::  tinconf  !  Initial CPU conformational time
       real(kind=8)                                ::  tfinconf !  Final CPU conformational time
       integer                                     ::  t1read   !  Initial reading time
       integer                                     ::  t2read   !  Final reading time
       integer                                     ::  t1adj    !  Initial building time
       integer                                     ::  t2adj    !  Final building time
       integer                                     ::  t1pim    !  Initial PIM analysis time
       integer                                     ::  t2pim    !  Final PIM analysis time
       integer                                     ::  t1conf   !  Initial conformational analysis time
       integer                                     ::  t2conf   !  Final conformational analysis time
!
! AnalysisPhenolMD variables
!
       integer,intent(in)                          ::  nsolv    !
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
       num(:)  = 0
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
! -----------------------------------------------------------
!
! Building the adjacency matrix
!
           call cpu_time(tinadj)
           call system_clock(t1adj)
!
           call buildadj(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,    &
                         nat,mgrps,ngrps,igrps,msubg,nsubg,isubg,atms, &
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
                          imol,magg,debug)
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call printpop(xtcf%STEP,nsize,msize,box,nnode,nagg,magg,nsolv,  &
                         cin,uniout)
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
             call cpu_time(tinpim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfinpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfinpim - tinpim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
           end if
!
! Analyzing aggregates by their connectivity
!
           if ( doconf ) then
             call cpu_time(tinconf)
             call system_clock(t1conf)
!
!             call isomorphism()
!
             call cpu_time(tfinconf)
             call system_clock(t2conf)
!
             tcpuconf = tcpuconf + tfinconf - tinconf
             tconf    = tconf + dble(t2conf-t1conf)/dble(count_rate)
           end if
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
       subroutine agglife(xtcf,nat,nnode,natms,neidis,msize,nsteps,    &
                          nbody,ngrps,nsubg,ibody,igrps,isubg,body,    &
                          grps,subg,atms,mbody,mgrps,msubg,matms,      &
                          nprint,minstep,maxstep,nsolv,avlife,nlife,   &
                          dopim,cin,buildadjmol,debug)
!
       use xdr,         only:  xtcfile
       use omp_lib
!
       use properties,  only:  nmax,pim,num,pop,frac,conc,volu
!
       use omp_var
       use datatypes
       use timings
!
       use units,       only:  uniout
!
       use printings
       use utils
!
       use lifetimes
!
       implicit none
!
! System information
!
       type(xtcfile),intent(inout)                              ::  xtcf      !  Trajectory information
!
! Interaction criteria information
!
       real(kind=4),intent(in)                                  ::  neidis    !  Screening distance
!
! Average properties
!
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
       real(kind=8),intent(out)                                 ::  cin       !
       integer,intent(in)                                       ::  msize     !  Maximum aggregate size
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
! AnalysisPhenolMD variables
!
       integer,intent(in)                                       ::  nsolv     !
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
! Local lifetimes calculation variables
!
       integer,dimension(:),allocatable                         ::  life      !
       integer,dimension(:),allocatable                         ::  auxlife   !
!
! Declaration of time control variables
!
       real(kind=8)                                             ::  tinadj    !  Initial CPU building time
       real(kind=8)                                             ::  tfinadj   !  Final CPU building time
       real(kind=8)                                             ::  tinlife   !  Initial CPU lifetimes time
       real(kind=8)                                             ::  tfinlife  !  Final CPU lifetimes time
       real(kind=8)                                             ::  tinpim    !  Initial CPU PIM time
       real(kind=8)                                             ::  tfinpim   !  Final CPU PIM time
       integer                                                  ::  t1read    !  Initial reading time
       integer                                                  ::  t2read    !  Final reading time
       integer                                                  ::  t1adj     !  Initial building time
       integer                                                  ::  t2adj     !  Final building time
       integer                                                  ::  t1life    !  Initial lifetimes time
       integer                                                  ::  t2life    !  Final lifetimes time
       integer                                                  ::  t1pim     !  Initial PIM analysis time
       integer                                                  ::  t2pim     !  Final PIM analysis time
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable                  ::  posi      !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable                  ::  newposi   !  Auxiliary coordinates
       real(kind=4),dimension(3)                                ::  box       !
       real(kind=4),dimension(3)                                ::  newbox    !
       integer                                                  ::  i,j       !
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
       num(:)  = 0
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life)                                         &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
       do i = 1, nnode
         life(i)   = 0
       end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(avlife,nlife)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
       do i = 1, nmax
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
                     mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,         &
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
                      imol,magg,debug)
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
                         nat,mgrps,ngrps,igrps,msubg,nsubg,isubg,      &
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
                          newnagg,newiagg,newnmol,newimol,newmagg,debug)
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
                         nsolv,cin,uniout)
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
             call cpu_time(tinpim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfinpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfinpim - tinpim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
           end if
!
! Adding up lifetimes of the aggregates
!
           call cpu_time(tinlife)
           call system_clock(t1life)
!
           call calclife(nmax,avlife,nlife,nnode,life,nsize,nagg,iagg,      &
                         magg,iwont)
!
! Analyzing aggregates by their connectivity
!



!
! Keeping track the adjacency matrix information
!
           box(:)    = newbox(:)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(auxlife)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, nnode
             auxlife(i) = 0
           end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(nagg,willmap,life,auxlife)                    &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = nagg(1)+1, magg
             if ( willmap(i) .ne. 0 ) auxlife(willmap(i)) = life(i)
           end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,auxlife,mol,agg,tag,nagg,iagg,nmol,imol, &
!$omp                    newmol,newagg,newtag,newnagg,newiagg,newnmol, &
!$omp                    newimol)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, nnode
!
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(posi,newposi)                                 &
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
       subroutine aggscrn(xtcf,nat,nnode,natms,neidis,msize,nsteps,    &
                          nbody,ngrps,nsubg,ibody,igrps,isubg,body,    &
                          grps,subg,atms,mbody,mgrps,msubg,matms,      &
                          nprint,minstep,maxstep,nsolv,dopim,cin,      &
                          buildadjmol,screen,debug)
!
       use xdr,         only:  xtcfile
       use omp_lib
!
       use properties,  only:  pim,num,pop,frac,conc,volu
!
       use omp_var
       use datatypes
       use timings
       use units,       only:  uniout
!
       use printings
       use utils
!
       use screening
!
       implicit none
!
! System information
!
       type(xtcfile),intent(inout)                              ::  xtcf      !  Trajectory information
!
! Interaction criteria information
!
       real(kind=4),intent(in)                                  ::  neidis    !  Screening distance
!
! Average properties
!
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
       real(kind=8),intent(out)                                 ::  cin       !
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
       real(kind=8)                                             ::  tinpim    !  Initial CPU PIM time
       real(kind=8)                                             ::  tfinpim   !  Final CPU PIM time
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
      box(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
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
                     mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,box,     &
                     neidis,buildadjmol)
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
                     mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,box,     &
                     neidis,buildadjmol)
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
                         xtcf%pos,nat,mgrps,ngrps,igrps,msubg,nsubg,   &
                         isubg,atms,newbox,neidis,buildadjmol)
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
                          nmol,imol,magg,debug)
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call printpop(actstep,nsize,msize,box,nnode,nagg,magg,      &
                         nsolv,cin,uniout)
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
             call cpu_time(tinpim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfinpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfinpim - tinpim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,oldadj,newadj)                            &
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(posi,newposi)                                 &
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
       subroutine aggscrnlife(xtcf,nat,nnode,natms,neidis,msize,       &
                              nsteps,nbody,ngrps,nsubg,ibody,igrps,    &
                              isubg,body,grps,subg,atms,mbody,mgrps,   &
                              msubg,matms,nprint,minstep,maxstep,      &
                              nsolv,avlife,nlife,dopim,cin,            &
                              buildadjmol,screen,debug)
!
       use xdr,         only:  xtcfile
       use omp_lib
!
       use properties,  only:  nmax,pim,num,pop,frac,conc,volu
!
       use omp_var
       use datatypes
       use timings
!
       use units,       only:  uniout
!
       use printings
       use utils
!
       use screening
       use lifetimes
!
       implicit none
!
! System information
!
       type(xtcfile),intent(inout)                              ::  xtcf      !  Trajectory information
!
! Interaction criteria information
!
       real(kind=4),intent(in)                                  ::  neidis    !  Screening distance
!
! Average properties
!
       integer,intent(out)                                      ::  nsteps    !  Number of snapshots analyzed
       real(kind=8),intent(out)                                 ::  cin       !
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
       real(kind=8)                                             ::  tinpim    !  Initial CPU PIM time
       real(kind=8)                                             ::  tfinpim   !  Final CPU PIM time
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
       num(:)  = 0
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
                     nat,mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,     &
                     nextbox,neidis,buildadjmol)
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
                      newnagg,newiagg,newnmol,newimol,newmagg,debug)
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
                     mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,box,     &
                     neidis,buildadjmol)
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
                     nat,mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,     &
                     newbox,neidis,buildadjmol)
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
                      imol,magg,debug)
!
! Finding aggregates present in the old and the actual configurations
!
       call cpu_time(tinlife)
       call system_clock(t1life)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,avlife,nlife)                            &
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,iwill,nagg)                              &
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
                         xtcf%pos,nat,mgrps,ngrps,igrps,msubg,nsubg,   &
                         isubg,atms,nextbox,neidis,buildadjmol)
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
                          newnagg,newiagg,newnmol,newimol,newmagg,debug)
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
                         nsolv,cin,uniout)
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
             call cpu_time(tinpim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfinpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfinpim - tinpim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
           end if
!
! Adding up lifetimes of the aggregates
!
           call cpu_time(tinlife)
           call system_clock(t1life)
!
           call calclife(nmax,avlife,nlife,nnode,life,nsize,nagg,iagg,      &
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(adj,newadj,auxlife)                           &
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(nagg,willmap,life,auxlife)                    &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = nagg(1)+1, magg
             if ( willmap(i) .ne. 0 ) auxlife(willmap(i)) = life(i)
           end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,auxlife,mol,agg,tag,nagg,iagg,nmol,imol, &
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(posi,newposi,nextposi)                        &
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
! PRINTPOP - PRINT POPulations
!
! This subroutine
!
       subroutine printpop(step,nsize,msize,box,nnode,nagg,magg,       &
                           nsolv,cin,iuni)
!
       use properties,  only:  num,pop,frac,conc,volu
!
       use parameters,  only:  Na
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(3),intent(in)         ::  box    !  Simulation box dimensions
       real(kind=8),intent(inout)                   ::  cin    !
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
       write(iuni+5,'(I10,100(X,I10))') step,nagg(:msize-1),int(dp1)
!
       return
       end subroutine printpop
!
!======================================================================!
!
! NAGGDIST - N-components AGGregates DISTances algorithm
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
       subroutine naggdist(neidis,nsteps,nprint,minstep,maxstep,nsolv, &
                           dopim,doconf,domon,buildadjmol,             &
                           buildadjrep,buildadjmon,printadjrep,        &
                           docoord,debug)
!
       use systeminf,   only:  xtcf,mnode,matms,maxat,coord
       use properties,  only:  nmax,pim,num,pop,frac,conc,cin,volu
!
       use timings,     only:  count_rate,tread,tadj,tpim,tconf,       &
                               tcpuadj,tcpupim,tcpuconf
!
       use units,       only:  uniout
!
       implicit none
!
! Interaction criteria information
!
       real(kind=4),intent(in)                     ::  neidis   !  Screening distance
!
! Trajectory control variables
!
       integer,intent(in)                          ::  nprint   !  Populations printing interval
       integer,intent(in)                          ::  minstep  !  First step for analysis
       integer,intent(in)                          ::  maxstep  !  Last step for analysis
       integer,intent(out)                         ::  nsteps   !  Number of snapshots analyzed
!
! Program control flags
!
       logical,intent(in)                          ::  dopim    !  PIM calculation flag
       logical,intent(in)                          ::  doconf   !  Conformational analysis flag
       logical,intent(in)                          ::  domon    !  Monomer intramolecular edges flag
       logical,intent(in)                          ::  docoord  !  Coordinate printing flag
       logical,intent(in)                          ::  debug    !  Debug mode
!
! AnalysisPhenolMD variables
!
       integer,intent(in)                          ::  nsolv    !
!
! External functions
!
       external                                    ::  buildadjmol
       external                                    ::  buildadjrep
       external                                    ::  buildadjmon
       external                                    ::  printadjrep
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable          ::  adj      !  Adjacency matrix
       integer,dimension(:),allocatable            ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable            ::  node     !  Molecules identifier
       integer,dimension(:),allocatable            ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable            ::  agg      !  Aggregates size
       integer,dimension(:),allocatable            ::  idx      !  Aggregate identifier
       integer,dimension(:),allocatable            ::  ntype    !
       integer,dimension(:),allocatable            ::  itype    !
       integer,dimension(:),allocatable            ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable            ::  iagg     !
       integer,dimension(:),allocatable            ::  nmol     !
       integer,dimension(:),allocatable            ::  imol     !
       integer                                     ::  nsize    !  Actual maximum aggregate size
       integer                                     ::  midx     !  Actual maximum aggregate identifier
       integer                                     ::  magg     !  Actual number of chemical species
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable     ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(3)                   ::  box      !
!
! Declaration of time control variables
!
       real(kind=8)                                ::  tiadj    !  Initial CPU building time
       real(kind=8)                                ::  tfadj    !  Final CPU building time
       real(kind=8)                                ::  ticonf   !  Initial CPU conformational time
       real(kind=8)                                ::  tfconf   !  Final CPU conformational time
       real(kind=8)                                ::  tipim    !  Initial CPU PIM time
       real(kind=8)                                ::  tfpim    !  Final CPU PIM time
       integer                                     ::  t1read   !  Initial reading time
       integer                                     ::  t2read   !  Final reading time
       integer                                     ::  t1adj    !  Initial building time
       integer                                     ::  t2adj    !  Final building time
       integer                                     ::  t1pim    !  Initial PIM analysis time
       integer                                     ::  t2pim    !  Final PIM analysis time
       integer                                     ::  t1conf   !  Initial conformational analysis time
       integer                                     ::  t2conf   !  Final conformational analysis time
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
       num(:)  = 0
!
       cin(:)  = 0.0d0
       volu    = 0.0d0
!
! Allocating variables depending on system size
!
       allocate(adj(mnode,mnode))
!
       allocate(mol(mnode),node(mnode),tag(mnode),agg(mnode))
       allocate(idx(mnode),ntype(mnode),itype(mnode))
!
       allocate(nmol(nmax),imol(nmax))
       allocate(nagg(nmax),iagg(nmax))
!
! Allocating variables depending on topological information
!
       allocate(posi(3,matms))
!
       coord => xtcf%pos(:,:)
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
! ...........................................................
!
! Building the adjacency matrix
!
           call cpu_time(tiadj)
           call system_clock(t1adj)
!
           call nbuildadj(adj,posi,xtcf%pos,box,neidis,buildadjmol)
!
           call cpu_time(tfadj)
           call system_clock(t2adj)
!
           tcpuadj = tcpuadj + tfadj - tiadj
           tadj    = tadj    + dble(t2adj-t1adj)/dble(count_rate)
!
! Block-diagonalizing the adjacency matrix
!
           call nblockdiag(adj,mol,node,tag,agg,idx,ntype,itype,       &
                           nsize,nagg,iagg,nmol,imol,magg,midx,debug)
!
! Analyzing aggregates by their connectivity
!
           if ( doconf ) then
             call cpu_time(ticonf)
             call system_clock(t1conf)
!
             call printadjrep(nmax,nagg,imol,mnode,node,matms,posi,    &
                              box,buildadjrep,buildadjmon,domon)
!
             call cpu_time(tfconf)
             call system_clock(t2conf)
!
             tcpuconf = tcpuconf + tfconf - ticonf
             tconf    = tconf    + dble(t2conf-t1conf)/dble(count_rate)
           end if
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call nprintpop(xtcf%STEP,box,nagg,imol,mnode,agg,itype,    &
                          magg,nsolv,uniout)
!
          if ( docoord ) then
            call nprint_coord(xtcf%STEP,xtcf%pos,box,maxat,nmax,nagg, &
                              imol,mnode,mol,node,agg,itype)
          end if
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Adding up the pairwise interaction matrix of the current snapshot
!
           if ( dopim ) then
             call cpu_time(tipim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfpim - tipim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
           end if
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
       deallocate(mol,tag,agg,idx,ntype,itype)
       deallocate(nmol,imol,nagg,iagg)
!
! Saving timings
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate)
!
       return
       end subroutine naggdist
!
!======================================================================!
!
! NAGGLIFE - N-components AGGregates LIFEtimes
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
       subroutine nagglife(neidis,nsteps,nprint,minstep,maxstep,nsolv, &
                           avlife,nlife,dopim,doconf,domon,            &
                           buildadjmol,buildadjrep,buildadjmon,        &
                           printadjrep,docoord,debug)
!
       use omp_lib
!
       use systeminf,   only:  xtcf,rep,mtype,mnode,matms,mmon,maxat,  &
                               coord
       use properties,  only:  nmax,pim,num,pop,frac,conc,cin,volu
!
       use timings,     only:  count_rate,tread,tadj,tlife,tpim,tconf, &
                               tcpuadj,tcpupim,tcpulife,tcpuconf
!
       use units,       only:  uniout
!
       use printings,   only:  print_end
!
       use lifetimes,   only:  ntracklife,ncalclife
!
       use omp_var,     only:  np,chunklife
!
       implicit none
!
! Interaction criteria information
!
       real(kind=4),intent(in)                         ::  neidis    !  Screening distance
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(nmax),intent(out)        ::  avlife    !
       integer,dimension(nmax),intent(out)             ::  nlife     !
!
! Trajectory control variables
!
       integer,intent(in)                              ::  nprint    !  Populations printing interval
       integer,intent(in)                              ::  minstep   !  First step for analysis
       integer,intent(in)                              ::  maxstep   !  Last step for analysis
       integer,intent(out)                             ::  nsteps    !  Number of snapshots analyzed
!
! Program control flags
!
       logical,intent(in)                              ::  dopim     !  PIM calculation flag
       logical,intent(in)                              ::  doconf    !  Conformational analysis flag
       logical,intent(in)                              ::  domon     !  Monomer intramolecular edges flag
       logical,intent(in)                              ::  docoord   !  Coordinate printing flag
       logical,intent(in)                              ::  debug     !  Debug mode
!
! AnalysisPhenolMD variables
!
       integer,intent(in)                              ::  nsolv     !
!
! External functions
!
       external                                        ::  buildadjmol
       external                                        ::  buildadjrep
       external                                        ::  buildadjmon
       external                                        ::  printadjrep
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable              ::  adj       !  Adjacency matrix
       logical,dimension(:),allocatable                ::  iwill     !
       logical,dimension(:),allocatable                ::  iwont     !
       integer,dimension(:),allocatable                ::  mol       !  Molecules identifier
       integer,dimension(:),allocatable                ::  node      !  Molecules identifier
       integer,dimension(:),allocatable                ::  tag       !  Aggregates identifier
       integer,dimension(:),allocatable                ::  agg       !  Aggregates size
       integer,dimension(:),allocatable                ::  idx       !  Aggregate identifier
       integer,dimension(:),allocatable                ::  ntype     !
       integer,dimension(:),allocatable                ::  itype     !
       integer,dimension(:),allocatable                ::  newmol    !  Molecules identifier
       integer,dimension(:),allocatable                ::  newnode   !  Molecules identifier
       integer,dimension(:),allocatable                ::  newtag    !  Aggregates identifier
       integer,dimension(:),allocatable                ::  newagg    !  Aggregates size
       integer,dimension(:),allocatable                ::  newidx    !  Aggregate identifier
       integer,dimension(:),allocatable                ::  newntype  !
       integer,dimension(:),allocatable                ::  newitype  !
       integer,dimension(:),allocatable                ::  nagg      !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  iagg      !
       integer,dimension(:),allocatable                ::  nmol      !
       integer,dimension(:),allocatable                ::  imol      !
       integer,dimension(:),allocatable                ::  newnagg   !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  newiagg   !
       integer,dimension(:),allocatable                ::  newnmol   !
       integer,dimension(:),allocatable                ::  newimol   !
       integer,dimension(:),allocatable                ::  willmap   !
       integer                                         ::  nsize     !  Actual maximum aggregate size
       integer                                         ::  newsize   !  New maximum aggregate size
       integer                                         ::  midx      !  Actual maximum aggregate identifier
       integer                                         ::  newmidx   !  New maximum aggregate identifier
       integer                                         ::  magg      !  Actual number of chemical species
       integer                                         ::  newmagg   !  New number of chemical species
       integer                                         ::  actstep   !
       integer                                         ::  newstep   !
!
! Local lifetimes calculation variables
!
       integer,dimension(:),allocatable                ::  life      !
       integer,dimension(:),allocatable                ::  auxlife   !
!
! Local variables
!
       real(kind=4),dimension(:,:),allocatable         ::  posi      !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable         ::  newposi   !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable,target  ::  tmpposi   !  Auxiliary coordinates
       real(kind=4),dimension(3)                       ::  box       !
       real(kind=4),dimension(3)                       ::  newbox    !
       integer                                         ::  i,j       !
!
! Declaration of time control variables
!
       real(kind=8)                                    ::  tiadj     !  Initial CPU building time
       real(kind=8)                                    ::  tfadj     !  Final CPU building time
       real(kind=8)                                    ::  tipim     !  Initial CPU PIM time
       real(kind=8)                                    ::  tfpim     !  Final CPU PIM time
       real(kind=8)                                    ::  ticonf    !  Initial CPU conformational time
       real(kind=8)                                    ::  tfconf    !  Final CPU conformational time
       real(kind=8)                                    ::  tilife    !  Initial CPU lifetimes time
       real(kind=8)                                    ::  tflife    !  Final CPU lifetimes time
       integer                                         ::  t1read    !  Initial reading time
       integer                                         ::  t2read    !  Final reading time
       integer                                         ::  t1adj     !  Initial building time
       integer                                         ::  t2adj     !  Final building time
       integer                                         ::  t1pim     !  Initial PIM analysis time
       integer                                         ::  t2pim     !  Final PIM analysis time
       integer                                         ::  t1conf    !  Initial conformational analysis time
       integer                                         ::  t2conf    !  Final conformational analysis time
       integer                                         ::  t1life    !  Initial lifetimes time
       integer                                         ::  t2life    !  Final lifetimes time
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
       num(:)  = 0
!
       cin     = 0.0d0
       volu    = 0.0d0
!
! Allocating variables depending on system size
!
       allocate(adj(mnode,mnode))
!
       allocate(mol(mnode),node(mnode),tag(mnode),agg(mnode))
       allocate(idx(mnode),ntype(mnode),itype(mnode))
!
       allocate(newmol(mnode),newnode(mnode))
       allocate(newtag(mnode),newagg(mnode))
       allocate(newidx(mnode),newntype(mnode),newitype(mnode))
!
       allocate(nmol(nmax),imol(nmax))
       allocate(nagg(nmax),iagg(nmax))
!
       allocate(newnmol(nmax),newimol(nmax))
       allocate(newnagg(nmax),newiagg(nmax))
!
       allocate(iwill(mnode),iwont(mnode))
!
       allocate(willmap(mnode))
!
       allocate(life(mnode),auxlife(mnode))
!
! Allocating variables depending on topological information
!
       allocate(newposi(3,matms))
       allocate(posi(3,matms))
!
       allocate(tmpposi(3,maxat))
!
       coord => tmpposi(:,:)
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,tmpposi)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
       do i = 1, maxat
         tmpposi(:,i) = xtcf%pos(:,i)
       end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life)                                         &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
       do i = 1, mnode
         life(i)   = 0
       end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(avlife,nlife)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
       do i = 1, nmax
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
       call cpu_time(tiadj)
       call system_clock(t1adj)
!
       call nbuildadj(adj,posi,xtcf%pos,box,neidis,buildadjmol)
!
       call cpu_time(tfadj)
       call system_clock(t2adj)
!
       tcpuadj = tcpuadj + tfadj - tiadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
! Block-diagonalizing the adjacency matrix of the first configuration
!
        call nblockdiag(adj,mol,node,tag,agg,idx,ntype,itype,nsize,    &
                        nagg,iagg,nmol,imol,magg,midx,debug)
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
           call cpu_time(tiadj)
           call system_clock(t1adj)
!
           call nbuildadj(adj,newposi,xtcf%pos,newbox,neidis,buildadjmol)
!
           call cpu_time(tfadj)
           call system_clock(t2adj)
!
           tcpuadj = tcpuadj + tfadj - tiadj
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
! Block-diagonalizing the adjacency matrix of the new-configuration
!
           call nblockdiag(adj,newmol,newnode,newtag,newagg,newidx,    &
                           newntype,newitype,newsize,newnagg,newiagg,  &
                           newnmol,newimol,newmagg,newmidx,debug)
!
! Finding aggregates present in the new and the actual configurations
!
           call cpu_time(tilife)
           call system_clock(t1life)
!
           call ntracklife(mtype,mnode,nmax,mmon,newmidx,newsize,      &
                           newmol,midx,nsize,mol,newnagg,newiagg,      &
                           newnmol,newimol,nagg,iagg,imol,iwill,       &
                           iwont,willmap)
!
           call cpu_time(tflife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tflife - tilife
           tlife = tlife + dble(t2life-t1life)/dble(count_rate)
!
! Analyzing aggregates by their connectivity
!
           if ( doconf ) then
             call cpu_time(ticonf)
             call system_clock(t1conf)
!
             call printadjrep(nmax,nagg,imol,mnode,node,matms,posi,    &
                              box,buildadjrep,buildadjmon,domon)
!
             call cpu_time(tfconf)
             call system_clock(t2conf)
!
             tcpuconf = tcpuconf + tfconf - ticonf
             tconf    = tconf    + dble(t2conf-t1conf)/dble(count_rate)
           end if
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call nprintpop(actstep,box,nagg,imol,mnode,agg,itype,      &
                          magg,nsolv,uniout)
!
           if ( docoord ) then
             call nprint_coord(actstep,tmpposi,box,maxat,nmax,nagg,    &
                               imol,mnode,mol,node,agg,itype)
           end if
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Adding up the pairwise interaction matrix of the current snapshot
!
           if ( dopim ) then
             call cpu_time(tipim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfpim - tipim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
           end if
!
! Adding up lifetimes of the aggregates
!
           call cpu_time(tilife)
           call system_clock(t1life)
!
           call ncalclife(mtype,mnode,nmax,avlife,nlife,life,nmax,     &
                          nagg,iagg,magg,iwont)
!
! Keeping track the adjacency matrix information
!
           box(:) = newbox(:)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,tmpposi)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
           do i = 1, maxat
             tmpposi(:,i) = xtcf%pos(:,i)
           end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(auxlife)                                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, mnode
             auxlife(i) = 0
           end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(nagg,willmap,life,auxlife)                    &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = nagg(1)+1, magg
             if ( willmap(i) .ne. 0 ) auxlife(willmap(i)) = life(i)
           end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(life,auxlife,mol,agg,tag,newmol,              &
!$omp                    newagg,newtag)                                &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, mnode
             life(i) = auxlife(i)
             mol(i)  = newmol(i)
             agg(i)  = newagg(i)
             tag(i)  = newtag(i)
           end do
!
!$omp end parallel do
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(nagg,iagg,nmol,imol,newnagg,                  &
!$omp                    newiagg,newnmol,newimol)                      &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
           do i = 1, nmax
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
           midx    = newmidx
           nsize   = newsize
!
           actstep = newstep
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(posi,newposi)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
!
           do i = 1, matms
             posi(:,i) = newposi(:,i)
           end do
!
!$omp end parallel do
!
           call cpu_time(tflife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tflife - tilife
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
       deallocate(posi,newposi)
!
       deallocate(adj)
!
       deallocate(mol,tag,agg,idx,ntype,itype)
       deallocate(nmol,imol,nagg,iagg)
!
       deallocate(newmol,newtag,newagg,newidx,newntype,newitype)
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
       end subroutine nagglife
!
!======================================================================!
!
! NAGGSCRN - N-components AGGregates SCReeNing algorithm
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
       subroutine naggscrn(neidis,nsteps,nprint,minstep,maxstep,nsolv, &
                           dopim,doconf,domon,buildadjmol,buildadjrep, &
                           buildadjmon,screen,scrn,cconf,docoord,debug)
!
       use omp_lib
!
       use systeminf,   only:  xtcf,mnode,matms,maxat,                 &
                               coord
       use properties,  only:  nmax,pim,num,pop,frac,conc,cin,volu
       use lengths,     only:  lenschm
!
       use timings,     only:  count_rate,tread,tadj,tscrn,tshift,    &
                               tpim,tconf,trepscrn,tcpuadj,           &
                               tcpuscrn,tcpushift,tcpupim,tcpuconf,   &
                               tcpurepscrn
!
       use units,       only:  uniout
!
       use printings,   only:  print_end
!
       use omp_var,     only:  np,chunklife
!
       implicit none
!
! Interaction criteria information
!
       real(kind=4),intent(in)                         ::  neidis   !  Screening distance
!
! Trajectory control variables

!
       integer,intent(in)                              ::  nprint   !  Populations printing interval
       integer,intent(in)                              ::  minstep  !  First step for analysis
       integer,intent(in)                              ::  maxstep  !  Last step for analysis
       integer,intent(out)                             ::  nsteps   !  Number of snapshots analyzed
!
! Program control flags
!
       logical,intent(in)                              ::  dopim    !  PIM calculation flag
       logical,intent(in)                              ::  doconf   !  Conformational analysis flag
       logical,intent(in)                              ::  domon    !  Monomer intramolecular edges flag
       logical,intent(in)                              ::  docoord  !  Coordinate printing flag
       logical,intent(in)                              ::  debug    !  Debug mode
       character(len=lenschm),intent(in)               ::  scrn     !  Calculation screening flag
       character(len=lenschm),intent(in)               ::  cconf    !  Calculation conformations flag
!
! AnalysisPhenolMD variables
!
       integer,intent(in)                              ::  nsolv    !
!
! External functions
!
       external                                        ::  buildadjmol
       external                                        ::  buildadjrep
       external                                        ::  buildadjmon
       external                                        ::  screen
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),pointer                  ::  adj      !  Adjacency matrix
       logical,dimension(:,:),pointer                  ::  oldadj   !  Adjacency matrix
       logical,dimension(:,:),pointer                  ::  newadj   !  Adjacency matrix
       logical,dimension(:,:),pointer                  ::  tmpadj   !  Auxiliary adjacency pointer
       logical,dimension(:,:),allocatable,target       ::  adjbuf1  !  Adjacency matrix buffer
       logical,dimension(:,:),allocatable,target       ::  adjbuf2  !  Adjacency matrix buffer
       logical,dimension(:,:),allocatable,target       ::  adjbuf3  !  Adjacency matrix buffer
       integer,dimension(:),allocatable                ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable                ::  node     !  Molecules identifier
       integer,dimension(:),allocatable                ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable                ::  agg      !  Aggregates size
       integer,dimension(:),allocatable                ::  idx      !  Aggregate identifier
       integer,dimension(:),allocatable                ::  ntype    !
       integer,dimension(:),allocatable                ::  itype    !
       integer,dimension(:),allocatable                ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  iagg     !
       integer,dimension(:),allocatable                ::  nmol     !
       integer,dimension(:),allocatable                ::  imol     !
       integer,dimension(:),allocatable                ::  oldmol   !  Molecules identifier
       integer,dimension(:),allocatable                ::  oldnode  !  Molecules identifier
       integer,dimension(:),allocatable                ::  oldtag   !  Aggregates identifier
       integer,dimension(:),allocatable                ::  oldagg   !  Aggregates size
       integer,dimension(:),allocatable                ::  oldidx   !  Aggregate identifier
       integer,dimension(:),allocatable                ::  oldntype !
       integer,dimension(:),allocatable                ::  olditype !
       integer,dimension(:),allocatable                ::  oldnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  oldiagg  !
       integer,dimension(:),allocatable                ::  oldnmol  !
       integer,dimension(:),allocatable                ::  oldimol  !
       integer,dimension(:),allocatable                ::  newmol   !  Molecules identifier
       integer,dimension(:),allocatable                ::  newnode  !  Molecules identifier
       integer,dimension(:),allocatable                ::  newtag   !  Aggregates identifier
       integer,dimension(:),allocatable                ::  newagg   !  Aggregates size
       integer,dimension(:),allocatable                ::  newidx   !  Aggregate identifier
       integer,dimension(:),allocatable                ::  newntype !
       integer,dimension(:),allocatable                ::  newitype !
       integer,dimension(:),allocatable                ::  newnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  newiagg  !
       integer,dimension(:),allocatable                ::  newnmol  !
       integer,dimension(:),allocatable                ::  newimol  !
       integer                                         ::  oldsize  !  Old maximum aggregate size
       integer                                         ::  oldmidx  !  Old maximum aggregate index
       integer                                         ::  oldmagg  !  Old number of chemical species
       integer                                         ::  newsize  !  New maximum aggregate size
       integer                                         ::  newmidx  !  New maximum aggregate index
       integer                                         ::  newmagg  !  New number of chemical species
       integer                                         ::  nsize    !  Actual maximum aggregate size
       integer                                         ::  midx     !  Actual maximum aggregate index
       integer                                         ::  magg     !  Actual number of chemical species
       integer                                         ::  actstep  !
       integer                                         ::  newstep  !
!
! Aggregates information in the selected virtual representation
!
       logical,dimension(:,:),allocatable              ::  sysbase  !  Global representation template matrix
       logical,dimension(:,:),pointer                  ::  oldrep   !  Previous screened representation matrix
       logical,dimension(:,:),pointer                  ::  sysrep   !  Current representation matrix
       logical,dimension(:,:),pointer                  ::  newrep   !  Next representation matrix
       logical,dimension(:,:),pointer                  ::  tmprep   !  Auxiliary representation pointer
       logical,dimension(:,:),allocatable,target       ::  repbuf1  !  Global representation matrix buffer
       logical,dimension(:,:),allocatable,target       ::  repbuf2  !  Global representation matrix buffer
       logical,dimension(:,:),allocatable,target       ::  repbuf3  !  Global representation matrix buffer
       integer,dimension(:),allocatable                ::  irepnode !  First representation index of each molecule
       integer,dimension(:),allocatable                ::  nrepnode !  Number of representation sites per molecule
       integer                                         ::  msysrep  ! Total size of the global representation matrix
!
! Local variables
!
       real(kind=4),dimension(:,:),pointer             ::  oldposi  !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  newposi  !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  tmpposi1 !  Auxiliary coordinates pointer
       real(kind=4),dimension(:,:),allocatable,target  ::  posibuf1 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  posibuf2 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  posibuf3 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),pointer             ::  oldtmpposi !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  tmpposi   !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  newtmpposi !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  tmpcoord !  Auxiliary coordinates pointer
       real(kind=4),dimension(:,:),allocatable,target  ::  coordbuf1 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  coordbuf2 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  coordbuf3 !  Auxiliary coordinates buffer
       real(kind=4),dimension(3)                       ::  oldbox   !
       real(kind=4),dimension(3)                       ::  box      !
       real(kind=4),dimension(3)                       ::  newbox   !
       logical                                         ::  dobdir   !
       logical                                         ::  dobody   !
       integer                                         ::  i        !
!
! Declaration of time control variables
!
       real(kind=8)                                    ::  tiadj    !  Initial CPU building time
       real(kind=8)                                    ::  tfadj    !  Final CPU building time
       real(kind=8)                                    ::  tiscrn   !  Initial CPU screening time
       real(kind=8)                                    ::  tfscrn   !  Final CPU screening time
       real(kind=8)                                    ::  tishift  !  Initial CPU state-shift time
       real(kind=8)                                    ::  tfshift  !  Final CPU state-shift time
       real(kind=8)                                    ::  ticonf   !  Initial CPU conformational time
       real(kind=8)                                    ::  tfconf   !  Final CPU conformational time
       real(kind=8)                                    ::  tipim    !  Initial CPU PIM time
       real(kind=8)                                    ::  tfpim    !  Final CPU PIM time
       integer                                         ::  t1read   !  Initial reading time
       integer                                         ::  t2read   !  Final reading time
       integer                                         ::  t1adj    !  Initial building time
       integer                                         ::  t2adj    !  Final building time
       integer                                         ::  t1scrn   !  Initial screening time
       integer                                         ::  t2scrn   !  Final screening time
       integer                                         ::  t1shift  !  Initial state-shift time
       integer                                         ::  t2shift  !  Final state-shift time
       integer                                         ::  t1pim    !  Initial PIM analysis time
       integer                                         ::  t2pim    !  Final PIM analysis time
       integer                                         ::  t1conf   !  Initial conformational analysis time
       integer                                         ::  t2conf   !  Final conformational analysis time
!
! Initializing variables
!
       call system_clock(t1read)
!
       nsteps  = 0
!
       dobdir  = trim(cconf) .eq. 'bodydir'
       dobody  = dobdir .or. (trim(cconf) .eq. 'body')
!
       pim(:,:,:) = 0.0d0
!
       pop(:)  = 0.0d0
       conc(:) = 0.0d0
       frac(:) = 0.0d0
       num(:)  = 0
!
       cin     = 0.0d0
       volu    = 0.0d0
!
! Allocating variables depending on system size
!
       allocate(adjbuf1(mnode,mnode))
       allocate(adjbuf2(mnode,mnode))
       allocate(adjbuf3(mnode,mnode))
       oldadj => adjbuf1
       adj    => adjbuf2
       newadj => adjbuf3
!
       allocate(mol(mnode),node(mnode),tag(mnode),agg(mnode))
       allocate(idx(mnode),ntype(mnode),itype(mnode))
       allocate(oldmol(mnode),oldnode(mnode),oldtag(mnode),            &
                oldagg(mnode))
       allocate(oldidx(mnode),oldntype(mnode),olditype(mnode))
       allocate(newmol(mnode),newnode(mnode),newtag(mnode),            &
                newagg(mnode))
       allocate(newidx(mnode),newntype(mnode),newitype(mnode))
!
       allocate(nmol(nmax),imol(nmax))
       allocate(nagg(nmax),iagg(nmax))
       allocate(oldnmol(nmax),oldimol(nmax))
       allocate(oldnagg(nmax),oldiagg(nmax))
       allocate(newnmol(nmax),newimol(nmax))
       allocate(newnagg(nmax),newiagg(nmax))
!
       if ( doconf ) then
         allocate(irepnode(mnode),nrepnode(mnode))
         call setsysrepidx(dobody,msysrep,irepnode,nrepnode)
         allocate(sysbase(msysrep,msysrep))
         allocate(repbuf1(msysrep,msysrep))
         allocate(repbuf2(msysrep,msysrep))
         allocate(repbuf3(msysrep,msysrep))
         oldrep => repbuf1
         sysrep => repbuf2
         newrep => repbuf3
         call buildsysbaseadjrep(dobody,msysrep,irepnode,sysbase)
       end if
!
! Allocating variables depending on topological information
!
       allocate(posibuf1(3,matms))
       allocate(posibuf2(3,matms))
       allocate(posibuf3(3,matms))
       oldposi => posibuf1
       posi    => posibuf2
       newposi => posibuf3
!
       allocate(coordbuf1(3,maxat))
       allocate(coordbuf2(3,maxat))
       allocate(coordbuf3(3,maxat))
       oldtmpposi => coordbuf1
       tmpposi    => coordbuf2
       newtmpposi => coordbuf3
!
       coord => tmpposi(:,:)
!
! Initializing the old configuration
! ----------------------------------
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
      oldbox(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,oldtmpposi)                              &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
      do i = 1, maxat
        oldtmpposi(:,i) = xtcf%pos(:,i)
      end do
!
!$omp end parallel do
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Building adjacency matrix of the first old-configuration
!
       call cpu_time(tiadj)
       call system_clock(t1adj)
!
       call nbuildadj(oldadj,oldposi,xtcf%pos,oldbox,neidis,           &
                      buildadjmol)
!
       call cpu_time(tfadj)
       call system_clock(t2adj)
!
       tcpuadj = tcpuadj + tfadj - tiadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
       if ( doconf ) then
         call nblockdiag(oldadj,oldmol,oldnode,oldtag,oldagg,oldidx,   &
                         oldntype,olditype,oldsize,oldnagg,oldiagg,    &
                         oldnmol,oldimol,oldmagg,oldmidx,debug)
         call buildsysadjrep(nmax,oldnagg,oldimol,mnode,oldmol,        &
                             oldnode,msysrep,irepnode,                 &
                             nrepnode,oldrep,matms,oldposi,            &
                             oldtmpposi,oldbox,buildadjrep,            &
                             buildadjmon,domon,dobody)
       end if
!
! Initializing the actual configuration
! -------------------------------------
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,tmpposi)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
       do i = 1, maxat
         tmpposi(:,i) = xtcf%pos(:,i)
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
       call cpu_time(tiadj)
       call system_clock(t1adj)
!
       call nbuildadj(adj,posi,xtcf%pos,box,neidis,buildadjmol)
!
       call cpu_time(tfadj)
       call system_clock(t2adj)
!
       tcpuadj = tcpuadj + tfadj - tiadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
       if ( doconf ) then
         call nblockdiag(adj,mol,node,tag,agg,idx,ntype,itype,         &
                         nsize,nagg,iagg,nmol,imol,magg,midx,debug)
         call buildsysadjrep(nmax,nagg,imol,mnode,mol,node,            &
                             msysrep,irepnode,nrepnode,sysrep,matms,   &
                             posi,tmpposi,box,buildadjrep,             &
                             buildadjmon,domon,dobody)
       end if
!
! Processing the remaining trajectory
! -----------------------------------
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,newtmpposi)                              &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
          do i = 1, maxat
            newtmpposi(:,i) = xtcf%pos(:,i)
          end do
!
!$omp end parallel do
!
          call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Building adjacency matrix of the new-configuration
!
           call cpu_time(tiadj)
           call system_clock(t1adj)
!
           call nbuildadj(newadj,newposi,xtcf%pos,newbox,neidis,       &
                          buildadjmol)
!
           call cpu_time(tfadj)
           call system_clock(t2adj)
!
           tcpuadj = tcpuadj + tfadj - tiadj
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
           if ( doconf ) then
             call nblockdiag(newadj,newmol,newnode,newtag,newagg,       &
                             newidx,newntype,newitype,newsize,newnagg,  &
                             newiagg,newnmol,newimol,newmagg,newmidx,   &
                             debug)
             call buildsysadjrep(nmax,newnagg,newimol,mnode,newmol,    &
                                 newnode,msysrep,irepnode,nrepnode,    &
                                 newrep,matms,newposi,newtmpposi,      &
                                 newbox,buildadjrep,buildadjmon,       &
                                 domon,dobody)
           end if
!
! Screening interactions between the molecules
!
           call cpu_time(tiscrn)
           call system_clock(t1scrn)
!
           call screen(mnode,oldadj,adj,newadj)
!
           call cpu_time(tfscrn)
           call system_clock(t2scrn)
!
           tcpuscrn = tcpuscrn + tfscrn - tiscrn
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate)
!
! Block-diagonalizing the interaction-corrected adjacency matrix
!
           call nblockdiag(adj,mol,node,tag,agg,idx,ntype,itype,       &
                           nsize,nagg,iagg,nmol,imol,magg,midx,debug)
!
! Analyzing aggregates by their screened connectivity
!
           if ( doconf ) then
             call cpu_time(ticonf)
             call system_clock(t1conf)
!
             if ( .not. dobody ) then
!
! Screening actual adyacency matrix in the topological representation
!
               call scrnadjrep(scrn,msysrep,oldrep,sysrep,newrep,      &
                    sysbase,.FALSE.,screen)
             else
               call scrnadjrep(scrn,msysrep,oldrep,sysrep,newrep,      &
                    sysbase,dobdir,screen)
             end if
!
             call cpu_time(tfconf)
             call system_clock(t2conf)
!
             tcpurepscrn = tcpurepscrn + tfconf - ticonf
             trepscrn = trepscrn + dble(t2conf-t1conf)/dble(count_rate)
!
! Printing corrected adjacency matrix in the topological representation
!
             call cpu_time(ticonf)
             call system_clock(t1conf)
!
             if ( .not. dobody ) then
               call printsysadjrep(nmax,nagg,imol,mnode,mol,           &
                    msysrep,irepnode,nrepnode,sysrep,.FALSE.,          &
                    .FALSE.,domon)
             else
               call printsysadjrep(nmax,nagg,imol,mnode,mol,           &
                    msysrep,irepnode,nrepnode,sysrep,.TRUE.,           &
                    dobdir,domon)
             end if
!
             call cpu_time(tfconf)
             call system_clock(t2conf)
!
             tcpuconf = tcpuconf + tfconf - ticonf
             tconf = tconf + dble(t2conf-t1conf)/dble(count_rate)
           end if
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call nprintpop(actstep,box,nagg,imol,mnode,agg,itype,       &
                          magg,nsolv,uniout)
!
           if ( docoord ) then
             call nprint_coord(actstep,tmpposi,box,maxat,nmax,nagg,    &
                               imol,mnode,mol,node,agg,itype)
           end if
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Adding up the pairwise interaction matrix of the current snapshot
!
           if ( dopim ) then
             call cpu_time(tipim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfpim - tipim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
           end if
!
! Keeping track the adjacency matrix information
!
           call cpu_time(tishift)
           call system_clock(t1shift)
!
           oldbox(:) = box(:)
           box(:)    = newbox(:)
           actstep   = newstep
!
           tmpadj => oldadj
           oldadj => adj
           adj    => newadj
           newadj => tmpadj
!
           tmpcoord   => oldtmpposi
           oldtmpposi => tmpposi
           tmpposi    => newtmpposi
           newtmpposi => tmpcoord
           coord      => tmpposi(:,:)
!
           tmpposi1 => oldposi
           oldposi  => posi
           posi     => newposi
           newposi  => tmpposi1
!
           oldmol(:)  = mol(:)
           oldnode(:) = node(:)
           oldtag(:)  = tag(:)
           oldagg(:)  = agg(:)
           oldidx(:)  = idx(:)
!
           oldntype(:) = ntype(:)
           olditype(:) = itype(:)
           oldnagg(:)  = nagg(:)
           oldiagg(:)  = iagg(:)
           oldnmol(:)  = nmol(:)
           oldimol(:)  = imol(:)
!
           oldsize = nsize
           oldmidx = midx
           oldmagg = magg
!
           if ( doconf ) then
             tmprep => oldrep
             oldrep => sysrep
             sysrep => newrep
             newrep => tmprep
           end if
!
           call cpu_time(tfshift)
           call system_clock(t2shift)
!
           tcpushift = tcpushift + tfshift - tishift
           tshift = tshift + dble(t2shift-t1shift)/dble(count_rate)
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
       deallocate(posibuf1,posibuf2,posibuf3)
       deallocate(coordbuf1,coordbuf2,coordbuf3)
!
       deallocate(adjbuf1,adjbuf2,adjbuf3)
!
       if ( doconf ) then
         deallocate(sysbase,repbuf1,repbuf2,repbuf3)
         deallocate(irepnode,nrepnode)
       end if
!
       deallocate(mol,tag,agg,idx,ntype,itype)
       deallocate(oldmol,oldnode,oldtag,oldagg,oldidx,oldntype,olditype)
       deallocate(newmol,newnode,newtag,newagg,newidx,newntype,newitype)
       deallocate(nmol,imol,nagg,iagg)
       deallocate(oldnmol,oldimol,oldnagg,oldiagg)
       deallocate(newnmol,newimol,newnagg,newiagg)
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate)
!
       return
       end subroutine naggscrn
!
!======================================================================!
!
! NAGGSCRNLIFE - N-components AGGregates SCReeNing LIFEtime algorithm
!
       subroutine naggscrnlife(neidis,nsteps,nprint,minstep,maxstep,   &
                               nsolv,avlife,nlife,dopim,doconf,domon,  &
                               buildadjmol,buildadjrep,buildadjmon,    &
                               screen,scrn,cconf,docoord,debug)
!
       use omp_lib
!
       use systeminf,   only:  xtcf,mtype,mnode,matms,mmon,maxat,      &
                               coord
       use properties,  only:  nmax,pim,num,pop,frac,conc,cin,volu
       use lengths,     only:  lenschm
!
       use timings,     only:  count_rate,tread,tadj,tscrn,tshift,    &
                               tpim,tconf,trepscrn,tlife,             &
                               tcpuadj,tcpuscrn,tcpushift,tcpupim,    &
                               tcpuconf,tcpurepscrn,tcpulife
!
       use units,       only:  uniout
!
       use printings,   only:  print_end
!
       use lifetimes,   only:  ntracklife,ncalclife
!
       use omp_var,     only:  np,chunklife
!
       implicit none
!
! Interaction criteria information
!
       real(kind=4),intent(in)                         ::  neidis   !  Screening distance
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(nmax),intent(out)        ::  avlife   !
       integer,dimension(nmax),intent(out)             ::  nlife    !
!
! Trajectory control variables

!
       integer,intent(in)                              ::  nprint   !  Populations printing interval
       integer,intent(in)                              ::  minstep  !  First step for analysis
       integer,intent(in)                              ::  maxstep  !  Last step for analysis
       integer,intent(out)                             ::  nsteps   !  Number of snapshots analyzed
!
! Program control flags
!
       logical,intent(in)                              ::  dopim    !  PIM calculation flag
       logical,intent(in)                              ::  doconf   !  Conformational analysis flag
       logical,intent(in)                              ::  domon    !  Monomer intramolecular edges flag
       logical,intent(in)                              ::  docoord  !  Coordinate printing flag
       logical,intent(in)                              ::  debug    !  Debug mode
       character(len=lenschm),intent(in)               ::  scrn     !  Calculation screening flag
       character(len=lenschm),intent(in)               ::  cconf    !  Calculation conformations flag
!
! AnalysisPhenolMD variables
!
       integer,intent(in)                              ::  nsolv    !
!
! External functions
!
       external                                        ::  buildadjmol
       external                                        ::  buildadjrep
       external                                        ::  buildadjmon
       external                                        ::  screen
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),pointer                  ::  adj      !  Adjacency matrix
       logical,dimension(:,:),pointer                  ::  oldadj   !  Adjacency matrix
       logical,dimension(:,:),pointer                  ::  newadj   !  Adjacency matrix
       logical,dimension(:,:),pointer                  ::  tmpadj   !  Auxiliary adjacency pointer
       logical,dimension(:,:),allocatable,target       ::  adjbuf1  !  Adjacency matrix buffer
       logical,dimension(:,:),allocatable,target       ::  adjbuf2  !  Adjacency matrix buffer
       logical,dimension(:,:),allocatable,target       ::  adjbuf3  !  Adjacency matrix buffer
       integer,dimension(:),allocatable                ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable                ::  node     !  Molecules identifier
       integer,dimension(:),allocatable                ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable                ::  agg      !  Aggregates size
       integer,dimension(:),allocatable                ::  idx      !  Aggregate identifier
       integer,dimension(:),allocatable                ::  ntype    !
       integer,dimension(:),allocatable                ::  itype    !
       integer,dimension(:),allocatable                ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  iagg     !
       integer,dimension(:),allocatable                ::  nmol     !
       integer,dimension(:),allocatable                ::  imol     !
       integer,dimension(:),allocatable                ::  oldmol   !  Molecules identifier
       integer,dimension(:),allocatable                ::  oldnode  !  Molecules identifier
       integer,dimension(:),allocatable                ::  oldtag   !  Aggregates identifier
       integer,dimension(:),allocatable                ::  oldagg   !  Aggregates size
       integer,dimension(:),allocatable                ::  oldidx   !  Aggregate identifier
       integer,dimension(:),allocatable                ::  oldntype !
       integer,dimension(:),allocatable                ::  olditype !
       integer,dimension(:),allocatable                ::  oldnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  oldiagg  !
       integer,dimension(:),allocatable                ::  oldnmol  !
       integer,dimension(:),allocatable                ::  oldimol  !
       integer,dimension(:),allocatable                ::  newmol   !  Molecules identifier
       integer,dimension(:),allocatable                ::  newnode  !  Molecules identifier
       integer,dimension(:),allocatable                ::  newtag   !  Aggregates identifier
       integer,dimension(:),allocatable                ::  newagg   !  Aggregates size
       integer,dimension(:),allocatable                ::  newidx   !  Aggregate identifier
       integer,dimension(:),allocatable                ::  newntype !
       integer,dimension(:),allocatable                ::  newitype !
       integer,dimension(:),allocatable                ::  newnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable                ::  newiagg  !
       integer,dimension(:),allocatable                ::  newnmol  !
       integer,dimension(:),allocatable                ::  newimol  !
       logical,dimension(:),allocatable                ::  iwill    !
       logical,dimension(:),allocatable                ::  iwont    !
       integer,dimension(:),allocatable                ::  willmap  !
       integer,dimension(:),allocatable                ::  life     !
       integer,dimension(:),allocatable                ::  auxlife  !
       integer                                         ::  oldsize  !  Old maximum aggregate size
       integer                                         ::  oldmidx  !  Old maximum aggregate index
       integer                                         ::  oldmagg  !  Old number of chemical species
       integer                                         ::  newsize  !  New maximum aggregate size
       integer                                         ::  newmidx  !  New maximum aggregate index
       integer                                         ::  newmagg  !  New number of chemical species
       integer                                         ::  nsize    !  Actual maximum aggregate size
       integer                                         ::  midx     !  Actual maximum aggregate index
       integer                                         ::  magg     !  Actual number of chemical species
       integer                                         ::  actstep  !
       integer                                         ::  newstep  !
!
! Aggregates information in the selected virtual representation
!
       logical,dimension(:,:),allocatable              ::  sysbase  !  Global representation template matrix
       logical,dimension(:,:),pointer                  ::  oldrep   !  Previous screened representation matrix
       logical,dimension(:,:),pointer                  ::  sysrep   !  Current representation matrix
       logical,dimension(:,:),pointer                  ::  newrep   !  Next representation matrix
       logical,dimension(:,:),pointer                  ::  tmprep   !  Auxiliary representation pointer
       logical,dimension(:,:),allocatable,target       ::  repbuf1  !  Global representation matrix buffer
       logical,dimension(:,:),allocatable,target       ::  repbuf2  !  Global representation matrix buffer
       logical,dimension(:,:),allocatable,target       ::  repbuf3  !  Global representation matrix buffer
       integer,dimension(:),allocatable                ::  irepnode !  First representation index of each molecule
       integer,dimension(:),allocatable                ::  nrepnode !  Number of representation sites per molecule
       integer                                         ::  msysrep  ! Total size of the global representation matrix
!
! Local variables
!
       real(kind=4),dimension(:,:),pointer             ::  oldposi  !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  newposi  !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  tmpposi1 !  Auxiliary coordinates pointer
       real(kind=4),dimension(:,:),allocatable,target  ::  posibuf1 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  posibuf2 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  posibuf3 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),pointer             ::  oldtmpposi !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  tmpposi   !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  newtmpposi !  Auxiliary coordinates
       real(kind=4),dimension(:,:),pointer             ::  tmpcoord !  Auxiliary coordinates pointer
       real(kind=4),dimension(:,:),allocatable,target  ::  coordbuf1 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  coordbuf2 !  Auxiliary coordinates buffer
       real(kind=4),dimension(:,:),allocatable,target  ::  coordbuf3 !  Auxiliary coordinates buffer
       real(kind=4),dimension(3)                       ::  oldbox   !
       real(kind=4),dimension(3)                       ::  box      !
       real(kind=4),dimension(3)                       ::  newbox   !
       logical                                         ::  dobdir   !
       logical                                         ::  dobody   !
       logical                                         ::  haslife  !
       integer                                         ::  i        !
!
! Declaration of time control variables
!
       real(kind=8)                                    ::  tiadj    !  Initial CPU building time
       real(kind=8)                                    ::  tfadj    !  Final CPU building time
       real(kind=8)                                    ::  tiscrn   !  Initial CPU screening time
       real(kind=8)                                    ::  tfscrn   !  Final CPU screening time
       real(kind=8)                                    ::  tishift  !  Initial CPU state-shift time
       real(kind=8)                                    ::  tfshift  !  Final CPU state-shift time
       real(kind=8)                                    ::  ticonf   !  Initial CPU conformational time
       real(kind=8)                                    ::  tfconf   !  Final CPU conformational time
       real(kind=8)                                    ::  tipim    !  Initial CPU PIM time
       real(kind=8)                                    ::  tfpim    !  Final CPU PIM time
       real(kind=8)                                    ::  tilife   !  Initial CPU lifetimes time
       real(kind=8)                                    ::  tflife   !  Final CPU lifetimes time
       integer                                         ::  t1read   !  Initial reading time
       integer                                         ::  t2read   !  Final reading time
       integer                                         ::  t1adj    !  Initial building time
       integer                                         ::  t2adj    !  Final building time
       integer                                         ::  t1scrn   !  Initial screening time
       integer                                         ::  t2scrn   !  Final screening time
       integer                                         ::  t1shift  !  Initial state-shift time
       integer                                         ::  t2shift  !  Final state-shift time
       integer                                         ::  t1pim    !  Initial PIM analysis time
       integer                                         ::  t2pim    !  Final PIM analysis time
       integer                                         ::  t1life   !  Initial lifetimes time
       integer                                         ::  t2life   !  Final lifetimes time
       integer                                         ::  t1conf   !  Initial conformational analysis time
       integer                                         ::  t2conf   !  Final conformational analysis time
!
! Initializing variables
!
       call system_clock(t1read)
!
       nsteps  = 0
!
       dobdir  = trim(cconf) .eq. 'bodydir'
       dobody  = dobdir .or. (trim(cconf) .eq. 'body')
       haslife = .FALSE.
!
       pim(:,:,:) = 0.0d0
!
       pop(:)  = 0.0d0
       conc(:) = 0.0d0
       frac(:) = 0.0d0
       num(:)  = 0
!
       avlife(:) = 0.0d0
       nlife(:)  = 0
!
       cin     = 0.0d0
       volu    = 0.0d0
!
! Allocating variables depending on system size
!
       allocate(adjbuf1(mnode,mnode))
       allocate(adjbuf2(mnode,mnode))
       allocate(adjbuf3(mnode,mnode))
       oldadj => adjbuf1
       adj    => adjbuf2
       newadj => adjbuf3
!
       allocate(mol(mnode),node(mnode),tag(mnode),agg(mnode))
       allocate(idx(mnode),ntype(mnode),itype(mnode))
       allocate(oldmol(mnode),oldnode(mnode),oldtag(mnode),            &
                oldagg(mnode))
       allocate(oldidx(mnode),oldntype(mnode),olditype(mnode))
       allocate(newmol(mnode),newnode(mnode),newtag(mnode),            &
                newagg(mnode))
       allocate(newidx(mnode),newntype(mnode),newitype(mnode))
!
       allocate(nmol(nmax),imol(nmax))
       allocate(nagg(nmax),iagg(nmax))
       allocate(oldnmol(nmax),oldimol(nmax))
       allocate(oldnagg(nmax),oldiagg(nmax))
       allocate(newnmol(nmax),newimol(nmax))
       allocate(newnagg(nmax),newiagg(nmax))
!
       allocate(iwill(mnode),iwont(mnode),willmap(mnode))
       allocate(life(mnode),auxlife(mnode))
       life(:) = 0
       auxlife(:) = 0
!
       if ( doconf ) then
         allocate(irepnode(mnode),nrepnode(mnode))
         call setsysrepidx(dobody,msysrep,irepnode,nrepnode)
         allocate(sysbase(msysrep,msysrep))
         allocate(repbuf1(msysrep,msysrep))
         allocate(repbuf2(msysrep,msysrep))
         allocate(repbuf3(msysrep,msysrep))
         oldrep => repbuf1
         sysrep => repbuf2
         newrep => repbuf3
         call buildsysbaseadjrep(dobody,msysrep,irepnode,sysbase)
       end if
!
! Allocating variables depending on topological information
!
       allocate(posibuf1(3,matms))
       allocate(posibuf2(3,matms))
       allocate(posibuf3(3,matms))
       oldposi => posibuf1
       posi    => posibuf2
       newposi => posibuf3
!
       allocate(coordbuf1(3,maxat))
       allocate(coordbuf2(3,maxat))
       allocate(coordbuf3(3,maxat))
       oldtmpposi => coordbuf1
       tmpposi    => coordbuf2
       newtmpposi => coordbuf3
!
       coord => tmpposi(:,:)
!
! Initializing the old configuration
! ----------------------------------
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
      oldbox(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
!
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,oldtmpposi)                              &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
      do i = 1, maxat
        oldtmpposi(:,i) = xtcf%pos(:,i)
      end do
!
!$omp end parallel do
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Building adjacency matrix of the first old-configuration
!
       call cpu_time(tiadj)
       call system_clock(t1adj)
!
       call nbuildadj(oldadj,oldposi,xtcf%pos,oldbox,neidis,           &
                      buildadjmol)
!
       call cpu_time(tfadj)
       call system_clock(t2adj)
!
       tcpuadj = tcpuadj + tfadj - tiadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
       if ( doconf ) then
         call nblockdiag(oldadj,oldmol,oldnode,oldtag,oldagg,oldidx,   &
                         oldntype,olditype,oldsize,oldnagg,oldiagg,    &
                         oldnmol,oldimol,oldmagg,oldmidx,debug)
         call buildsysadjrep(nmax,oldnagg,oldimol,mnode,oldmol,        &
                             oldnode,msysrep,irepnode,                 &
                             nrepnode,oldrep,matms,oldposi,            &
                             oldtmpposi,oldbox,buildadjrep,            &
                             buildadjmon,domon,dobody)
       end if
!
! Initializing the actual configuration
! -------------------------------------
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,tmpposi)                                 &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
       do i = 1, maxat
         tmpposi(:,i) = xtcf%pos(:,i)
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
       call cpu_time(tiadj)
       call system_clock(t1adj)
!
       call nbuildadj(adj,posi,xtcf%pos,box,neidis,buildadjmol)
!
       call cpu_time(tfadj)
       call system_clock(t2adj)
!
       tcpuadj = tcpuadj + tfadj - tiadj
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
       if ( doconf ) then
         call nblockdiag(adj,mol,node,tag,agg,idx,ntype,itype,         &
                         nsize,nagg,iagg,nmol,imol,magg,midx,debug)
         call buildsysadjrep(nmax,nagg,imol,mnode,mol,node,            &
                             msysrep,irepnode,nrepnode,sysrep,matms,   &
                             posi,tmpposi,box,buildadjrep,             &
                             buildadjmon,domon,dobody)
       end if
!
! Processing the remaining trajectory
! -----------------------------------
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
!$omp parallel do num_threads(np)                                      &
!$omp             shared(xtcf,newtmpposi)                              &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,chunklife)
          do i = 1, maxat
            newtmpposi(:,i) = xtcf%pos(:,i)
          end do
!
!$omp end parallel do
!
          call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Building adjacency matrix of the new-configuration
!
           call cpu_time(tiadj)
           call system_clock(t1adj)
!
           call nbuildadj(newadj,newposi,xtcf%pos,newbox,neidis,       &
                          buildadjmol)
!
           call cpu_time(tfadj)
           call system_clock(t2adj)
!
           tcpuadj = tcpuadj + tfadj - tiadj
           tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)
!
           if ( doconf ) then
             call nblockdiag(newadj,newmol,newnode,newtag,newagg,       &
                             newidx,newntype,newitype,newsize,newnagg,  &
                             newiagg,newnmol,newimol,newmagg,newmidx,   &
                             debug)
             call buildsysadjrep(nmax,newnagg,newimol,mnode,newmol,    &
                                 newnode,msysrep,irepnode,nrepnode,    &
                                 newrep,matms,newposi,newtmpposi,      &
                                 newbox,buildadjrep,buildadjmon,       &
                                 domon,dobody)
           end if
!
! Screening interactions between the molecules
!
           call cpu_time(tiscrn)
           call system_clock(t1scrn)
!
           call screen(mnode,oldadj,adj,newadj)
!
           call cpu_time(tfscrn)
           call system_clock(t2scrn)
!
           tcpuscrn = tcpuscrn + tfscrn - tiscrn
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate)
!
! Block-diagonalizing the interaction-corrected adjacency matrix
!
           call nblockdiag(adj,mol,node,tag,agg,idx,ntype,itype,       &
                           nsize,nagg,iagg,nmol,imol,magg,midx,debug)
!
! Finding aggregates present in consecutive screened configurations
!
           call cpu_time(tilife)
           call system_clock(t1life)
!
           if ( haslife ) then
             call ntracklife(mtype,mnode,nmax,mmon,midx,nsize,mol,    &
                             oldmidx,oldsize,oldmol,nagg,iagg,nmol,   &
                             imol,oldnagg,oldiagg,oldimol,            &
                             iwill,iwont,willmap)
           else
             iwill(:) = .FALSE.
             iwont(:) = .FALSE.
             willmap(:) = 0
           end if
!
           call cpu_time(tflife)
           call system_clock(t2life)
!
           tcpulife = tcpulife + tflife - tilife
           tlife = tlife + dble(t2life-t1life)/dble(count_rate)
!
! Analyzing aggregates by their screened connectivity
!
           if ( doconf ) then
!
! Screening actual adjacency matrix in the topological representation
!
             call cpu_time(ticonf)
             call system_clock(t1conf)
!
             if ( .not. dobody ) then
               call scrnadjrep(scrn,msysrep,oldrep,sysrep,newrep,      &
                    sysbase,.FALSE.,screen)
             else
               call scrnadjrep(scrn,msysrep,oldrep,sysrep,newrep,      &
                    sysbase,dobdir,screen)
             end if
!
             call cpu_time(tfconf)
             call system_clock(t2conf)
!
             tcpurepscrn = tcpurepscrn + tfconf - ticonf
             trepscrn = trepscrn + dble(t2conf-t1conf)/dble(count_rate)
!
! Printing actual adjacency matrix in the topological representation
!
             call cpu_time(ticonf)
             call system_clock(t1conf)
!
             if ( .not. dobody ) then
               call printsysadjrep(nmax,nagg,imol,mnode,mol,           &
                    msysrep,irepnode,nrepnode,sysrep,.FALSE.,          &
                    .FALSE.,domon)
             else
               call printsysadjrep(nmax,nagg,imol,mnode,mol,           &
                    msysrep,irepnode,nrepnode,sysrep,.TRUE.,           &
                    dobdir,domon)
             end if
!
             call cpu_time(tfconf)
             call system_clock(t2conf)
!
             tcpuconf = tcpuconf + tfconf - ticonf
             tconf = tconf + dble(t2conf-t1conf)/dble(count_rate)
           end if
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           call nprintpop(actstep,box,nagg,imol,mnode,agg,itype,       &
                          magg,nsolv,uniout)
!
           if ( docoord ) then
             call nprint_coord(actstep,tmpposi,box,maxat,nmax,nagg,    &
                               imol,mnode,mol,node,agg,itype)
           end if
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate)
!
! Adding up the pairwise interaction matrix of the current snapshot
!
           if ( dopim ) then
             call cpu_time(tipim)
             call system_clock(t1pim)
!
!             call build_pim()
!
             call cpu_time(tfpim)
             call system_clock(t2pim)
!
             tcpupim = tcpupim + tfpim - tipim
             tpim    = tpim    + dble(t2pim-t1pim)/dble(count_rate)
           end if
!
! Adding up lifetimes of the aggregates
!
           call cpu_time(tilife)
           call system_clock(t1life)
!
           if ( haslife ) then
             call ncalclife(mtype,mnode,nmax,avlife,nlife,life,nmax,  &
                            oldnagg,oldiagg,oldmagg,iwont)
           end if
!
! Keeping track the adjacency matrix information
!
           call cpu_time(tishift)
           call system_clock(t1shift)
!
           oldbox(:) = box(:)
           box(:)    = newbox(:)
           actstep   = newstep
!
           tmpadj => oldadj
           oldadj => adj
           adj    => newadj
           newadj => tmpadj
!
           tmpcoord   => oldtmpposi
           oldtmpposi => tmpposi
           tmpposi    => newtmpposi
           newtmpposi => tmpcoord
           coord      => tmpposi(:,:)
!
           tmpposi1 => oldposi
           oldposi  => posi
           posi     => newposi
           newposi  => tmpposi1
!
           auxlife(:) = 0
           if ( haslife ) then
             do i = oldnagg(mtype+1)+1, oldmagg
               if ( willmap(i) .ne. 0 ) auxlife(willmap(i)) = life(i)
             end do
           end if
           life(:) = auxlife(:)
!
           oldmol(:)  = mol(:)
           oldnode(:) = node(:)
           oldtag(:)  = tag(:)
           oldagg(:)  = agg(:)
           oldidx(:)  = idx(:)
!
           oldntype(:) = ntype(:)
           olditype(:) = itype(:)
           oldnagg(:)  = nagg(:)
           oldiagg(:)  = iagg(:)
           oldnmol(:)  = nmol(:)
           oldimol(:)  = imol(:)
!
           oldsize = nsize
           oldmidx = midx
           oldmagg = magg
           haslife = .TRUE.
!
           if ( doconf ) then                                         
             tmprep => oldrep
             oldrep => sysrep
             sysrep => newrep
             newrep => tmprep
           end if                                                     
!
           call cpu_time(tfshift)
           call system_clock(t2shift)
!
           tcpushift = tcpushift + tfshift - tishift
           tshift = tshift + dble(t2shift-t1shift)/dble(count_rate)
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
       deallocate(posibuf1,posibuf2,posibuf3)
       deallocate(coordbuf1,coordbuf2,coordbuf3)
!
       deallocate(adjbuf1,adjbuf2,adjbuf3)
!
       if ( doconf ) then                                            
         deallocate(sysbase,repbuf1,repbuf2,repbuf3)
         deallocate(irepnode,nrepnode)                               
       end if                                                        
!
       deallocate(mol,tag,agg,idx,ntype,itype)
       deallocate(oldmol,oldnode,oldtag,oldagg,oldidx,oldntype,olditype)
       deallocate(newmol,newnode,newtag,newagg,newidx,newntype,newitype)
       deallocate(nmol,imol,nagg,iagg)
       deallocate(oldnmol,oldimol,oldnagg,oldiagg)
       deallocate(newnmol,newimol,newnagg,newiagg)
       deallocate(iwill,iwont,willmap)
       deallocate(life,auxlife)
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate)
!
       return
       end subroutine naggscrnlife
!
!======================================================================!
!
! NPRINTPOP - N-components PRINT POPulations
!
! This subroutine
!
       subroutine nprintpop(step,box,nagg,imol,mnode,agg,itype,magg,  &
                            nsolv,iuni)
!
       use systeminf,   only:  nnode,mtype,nmon
       use properties,  only:  nmax,num,pop,frac,conc,cin,volu,       &
                                xsize,psize,pqsize,homxsize,mixxsize, &
                                hompsize,mixpsize,                    &
                                sumagg,summol,sumqmol,hompop,homconc, &
                                homnum,mixpop,mixconc,mixnum,mixapop, &
                                mixaconc,mixanum
!
       use parameters,  only:  Na
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(3),intent(in)  ::  box    !  Simulation box dimensions
       integer,dimension(nmax),intent(in)    ::  nagg   !
       integer,dimension(nmax),intent(in)    ::  imol   !
       integer,intent(in)                    ::  mnode  !
       integer,dimension(mnode),intent(in)   ::  agg    !
       integer,dimension(mnode),intent(in)   ::  itype  !
       integer,intent(in)                    ::  magg   !
       integer,intent(in)                    ::  nsolv  !
       integer,intent(in)                    ::  step   !
       integer,intent(in)                    ::  iuni   !
!
! Local variables
!
       integer,dimension(mtype)              ::  comp   !
       integer,dimension(mtype,mnode)        ::  ihom   !
       integer,dimension(mnode)              ::  imix   !
       integer                               ::  i      !  Index
       integer                               ::  j      !  Index
       integer                               ::  k      !  Index
       integer                               ::  nsp    !  Number of species in aggregate
       integer                               ::  qhom   !  Homogeneous aggregate species
       integer                               ::  isize  !  Aggregate size
!
       ihom(:,:) = 0
       imix(:)   = 0
!
! Accumulating properties
!
       do i = 1, nmax
         num(i)  = num(i)  + nagg(i)
         pop(i)  = pop(i)  + real(nagg(i))/magg*100
         frac(i) = frac(i) + real(nagg(i))/(magg+nsolv)
         conc(i) = conc(i) + real(nagg(i))/box(1)**3
       end do
!
! Accumulating denominators and size distributions
!
       sumagg = sumagg + dble(magg)
       summol = summol + dble(sum(nnode(:)))
       sumqmol(:) = sumqmol(:) + dble(nnode(:))
!
       do i = 1, nmax-1
         if ( nagg(i) .eq. 0 ) cycle
         isize = sum(nmon(:,i))
         xsize(isize) = xsize(isize) + dble(nagg(i))
         psize(isize) = psize(isize) + dble(isize*nagg(i))
         nsp  = 0
         qhom = 0
         do j = 1, mtype
           pqsize(j,isize) = pqsize(j,isize) + dble(nmon(j,i)*nagg(i))
           if ( nmon(j,i) .gt. 0 ) then
             nsp  = nsp + 1
             qhom = j
           end if
         end do
         if ( nsp .eq. 1 ) then
           ihom(qhom,isize) = ihom(qhom,isize) + nagg(i)
           homxsize(qhom,isize) = homxsize(qhom,isize) +             &
                                  dble(nagg(i))
           hompsize(qhom,isize) = hompsize(qhom,isize) +             &
                                   dble(isize*nagg(i))
         else
           imix(isize) = imix(isize) + nagg(i)
           mixxsize(isize) = mixxsize(isize) + dble(nagg(i))
           mixpsize(isize) = mixpsize(isize) + dble(isize*nagg(i))
           mixapop(i)  = mixapop(i)  + dble(nagg(i))/dble(magg)*100.0d0
           mixaconc(i) = mixaconc(i) + dble(nagg(i))/dble(box(1)**3)
           mixanum(i)  = mixanum(i)  + dble(nagg(i))
         end if
       end do
!
       if ( nagg(nmax) .gt. 0 ) then
         k = imol(nmax)
         do i = 1, nagg(nmax)
           isize = agg(k+1)
           comp(:) = 0
           do j = k+1, k+isize
             comp(itype(j)) = comp(itype(j)) + 1
           end do
!
           xsize(isize) = xsize(isize) + 1.0d0
           psize(isize) = psize(isize) + dble(isize)
           nsp  = 0
           qhom = 0
           do j = 1, mtype
             pqsize(j,isize) = pqsize(j,isize) + dble(comp(j))
             if ( comp(j) .gt. 0 ) then
               nsp  = nsp + 1
               qhom = j
             end if
           end do
           if ( nsp .eq. 1 ) then
             ihom(qhom,isize) = ihom(qhom,isize) + 1
             homxsize(qhom,isize) = homxsize(qhom,isize) + 1.0d0
             hompsize(qhom,isize) = hompsize(qhom,isize) + dble(isize)
           else
             imix(isize) = imix(isize) + 1
             mixxsize(isize) = mixxsize(isize) + 1.0d0
             mixpsize(isize) = mixpsize(isize) + dble(isize)
           end if
!
           k = k + isize
         end do
       end if
!
       do i = 1, mtype
         hompop(i,:)  = hompop(i,:)  + dble(ihom(i,:))/dble(magg)*100.0d0
         homconc(i,:) = homconc(i,:) + dble(ihom(i,:))/dble(box(1)**3)
         homnum(i,:)  = homnum(i,:)  + dble(ihom(i,:))
       end do
!
       mixpop(:)  = mixpop(:)  + dble(imix(:))/dble(magg)*100.0d0
       mixconc(:) = mixconc(:) + dble(imix(:))/dble(box(1)**3)
       mixnum(:)  = mixnum(:)  + dble(imix(:))
!
       cin(:) = cin(:) + real(nnode(:))/box(1)**3
!
       volu = volu + box(1)**3
!
! Printing populations of the current configuration
!
       write(iuni+1,'(I10,100(X,F10.6))') step,real(nagg(:))/magg*100
       write(iuni+2,'(I10,100(X,F10.6))') step,                        &
                                              real(nagg(:))/(magg+nsolv)
       write(iuni+3,'(I10,100(X,F10.6))') step,                        &
                                    real(nagg(:))/box(1)**3/(Na*1.0E-24)
       write(iuni+5,'(I10,100(X,I7))') step,nagg(:)
!
       return
       end subroutine nprintpop
!
!======================================================================!
!
! NPRINT_COORD - N-components PRINT aggregate COORDinates
!
! This subroutine prints the atomic coordinates of each aggregate in XYZ
! format for the current snapshot.
!
       subroutine nprint_coord(step,rcoord,box,natot,nmax,nagg,imol,  &
                               mnode,mol,node,agg,itype)
!
       use systeminf,   only:  sys,iat,mmon,mtype,agglabel
       use filenames,   only:  outp
       use lengths,     only:  lenout
       use geometry,    only:  sminimgvec
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(3,natot),intent(in)  ::  rcoord  !
       real(kind=4),dimension(3),intent(in)        ::  box     !
       integer,dimension(nmax),intent(in)          ::  nagg    !
       integer,dimension(nmax),intent(in)          ::  imol    !
       integer,dimension(mnode),intent(in)         ::  mol     !
       integer,dimension(mnode),intent(in)         ::  node    !
       integer,dimension(mnode),intent(in)         ::  agg     !
       integer,dimension(mnode),intent(in)         ::  itype   !
       integer,intent(in)                          ::  step    !
       integer,intent(in)                          ::  natot   !
       integer,intent(in)                          ::  nmax    !
       integer,intent(in)                          ::  mnode   !
!
! Local variables
!
       character(len=lenout+128),allocatable,save, &
                              dimension(:)         ::  used    !
       character(len=lenout+128)                   ::  fname   !
       character(len=lenout+64)                    ::  label   !
       character(len=16)                           ::  aux     !
       logical                                     ::  newfile !
       integer,dimension(mtype)                    ::  comp    !
       real(kind=4),dimension(3)                   ::  bar     !  Geometrical center
       real(kind=4),dimension(3)                   ::  refpos  !  Aggregate reference position
       real(kind=4),dimension(3)                   ::  molpos  !  Rebuilt molecule reference
       real(kind=4),dimension(3)                   ::  rpos    !  Rebuilt atom position
       integer                                     ::  iagg    !  Aggregate index
       integer                                     ::  iocc    !  Aggregate occurrence
       integer                                     ::  im      !  Molecule index
       integer                                     ::  ia      !  Atom index
       integer                                     ::  i       !  Index
       integer                                     ::  imini   !  Initial molecule index
       integer                                     ::  imax    !  Final molecule index
       integer                                     ::  isize   !  Aggregate size
       integer                                     ::  q       !  Molecule type
       integer                                     ::  iatom   !  Initial atom index
       integer                                     ::  natagg  !  Number of atoms
       integer                                     ::  iuni    !  Output unit
       integer,save                                ::  nused=0 !
!
       if ( .not. allocated(used) ) then
         allocate(used(mnode))
         used(:) = ''
       end if
!
       do iagg = 1, nmax
         if ( nagg(iagg) .eq. 0 ) cycle
!
         imini = imol(iagg)
!
         do iocc = 1, nagg(iagg)
           if ( iagg .eq. nmax ) then
             isize = agg(imini+1)
           else
             isize = mmon(iagg)
           end if
!
           imax   = imini + isize
           natagg = 0
           comp(:) = 0
!
           do im = imini+1, imax
             q      = itype(im)
             comp(q) = comp(q) + 1
             natagg = natagg + sys(q)%nat
           end do
!
           label = ''
           if ( iagg .eq. nmax ) then
             do q = 1, mtype
               if ( comp(q) .eq. 0 ) cycle
               write(aux,'(I0)') q
               label = trim(label)//'_q'//trim(adjustl(aux))
               write(aux,'(I0)') comp(q)
               label = trim(label)//'n'//trim(adjustl(aux))
             end do
           else
             label = agglabel(iagg)
           end if
!
           write(aux,'(I0)') iagg
           fname = trim(outp)//'_id_'//                                &
                   trim(adjustl(aux))//trim(label)//'.xyz'
!
           newfile = .TRUE.
           do i = 1, nused
             if ( trim(used(i)) .eq. trim(fname) ) then
               newfile = .FALSE.
               exit
             end if
           end do
!
           if ( newfile ) then
             nused = nused + 1
             used(nused) = fname
             open(newunit=iuni,file=trim(fname),status='replace',     &
                  action='write')
           else
             open(newunit=iuni,file=trim(fname),status='old',         &
                  position='append',action='write')
           end if
!
           write(iuni,'(I10)') natagg
           write(iuni,'(A,I0,A)',advance='no') 'snapshot ',step,      &
                                               ' molecule_ids'
           do im = imini+1, imax
             write(iuni,'(1X,I0)',advance='no') mol(im)
           end do
           write(iuni,*)
!
           q      = itype(imini+1)
           iatom  = iat(q) + (node(imini+1)-1)*sys(q)%nat
           refpos = rcoord(:,iatom+1)
           bar(:) = 0.0
!
           do im = imini+1, imax
             q     = itype(im)
             iatom = iat(q) + (node(im)-1)*sys(q)%nat
             molpos(:) = refpos(:) + sminimgvec(refpos,               &
                                  rcoord(:,iatom+1),box)
             do ia = 1, sys(q)%nat
               rpos(:) = molpos(:) + sminimgvec(rcoord(:,iatom+1),     &
                                                 rcoord(:,iatom+ia),box)
               bar(:) = bar(:) + rpos(:)
             end do
           end do
!
           bar(:) = bar(:)/real(natagg)
!
           do im = imini+1, imax
             q     = itype(im)
             iatom = iat(q) + (node(im)-1)*sys(q)%nat
             molpos(:) = refpos(:) + sminimgvec(refpos,                &
                                  rcoord(:,iatom+1),box)
             do ia = 1, sys(q)%nat
               rpos(:) = molpos(:) + sminimgvec(rcoord(:,iatom+1),     &
                                                 rcoord(:,iatom+ia),box)
               write(iuni,'(A2,3(1X,F14.6))') sys(q)%atsymb(ia),       &
                                              (rpos(1)-bar(1))*10,     &
                                              (rpos(2)-bar(2))*10,     &
                                              (rpos(3)-bar(3))*10
             end do
           end do
!
           close(iuni)
!
           imini = imax
         end do
       end do
!
       return
       end subroutine nprint_coord
!
!======================================================================!
!
! SETSYSREPIDX - SET SYStem REPresentation InDeXes                   
!                                                                     
       subroutine setsysrepidx(dobody,msysrep,irepnode,nrepnode)      
!                                                                     
       use systeminf,   only:  rep,mtype,nnode,inode                  
!                                                                     
       implicit none                                                  
!                                                                     
! Input/output variables                                              
!                                                                     
       integer,dimension(:),intent(out)  ::  irepnode                 ! First representation index of each molecule
       integer,dimension(:),intent(out)  ::  nrepnode                 ! Number of representation sites per molecule
       integer,intent(out)               ::  msysrep                  ! Total size of global representation matrix 
       logical,intent(in)                ::  dobody                   ! Body representation flag 
!                                                                     
! Local variables                                                     
!                                                                     
       integer                           ::  q                        ! Molecule type index 
       integer                           ::  imol                     ! Local molecule index within a type
       integer                           ::  molid                    ! Global molecule identifier
       integer                           ::  nrep                     ! Representation sites in one molecule
!                                                                      
       msysrep = 0                                                     
       do q = 1, mtype                                                 
         if ( dobody ) then                                            
           nrep = rep(q)%mbody                                         
         else                                                          
           nrep = rep(q)%mgrps                                         
         end if                                                        
         do imol = 1, nnode(q)                                         
           molid = inode(q) + imol                                     
           irepnode(molid) = msysrep + 1                               
           nrepnode(molid) = nrep                                      
           msysrep = msysrep + nrep                                    
         end do                                                        
       end do                                                          
!                                                                      
       return                                                          
       end subroutine setsysrepidx                                     
!                                                                       
!======================================================================!
!                                                                       
! BUILDSYSBASEADJREP - BUILD SYStem BASE ADJacency in REPresentation    
!                                                                       
       subroutine buildsysbaseadjrep(dobody,msysrep,irepnode,base)    
!                                                                    
       use systeminf,   only:  rep,mtype,nnode,inode                
!                                                                       
       implicit none                                               
!                                                                       
! Input/output variables                                                
!                                                                       
       logical,dimension(msysrep,msysrep),intent(out)  ::  base        ! Global template adjacency matrix 
       integer,dimension(:),intent(in)                 ::  irepnode    ! First representation index of each molecule
       integer,intent(in)                              ::  msysrep     ! Total size of global representation matrix
       logical,intent(in)                              ::  dobody      ! Body representation flag
!                                                                        
! Local variables                                                        
!                                                                        
       integer                                         ::  q           ! Molecule type index 
       integer                                         ::  imol        ! Local molecule index within a type
       integer                                         ::  molid       ! Global molecule identifier 
       integer                                         ::  i0          ! Initial global representation index
       integer                                         ::  nrep        ! Representation sites in one molecule
!                                                                      
       base(:,:) = .FALSE.                                             
       do q = 1, mtype                                                 
         if ( dobody ) then                                            
           nrep = rep(q)%mbody                                         
         else                                                          
           nrep = rep(q)%mgrps                                         
         end if                                                        
         do imol = 1, nnode(q)                                         
           molid = inode(q) + imol                                     
           i0 = irepnode(molid)                                        
           if ( dobody ) then                                          
             base(i0:i0+nrep-1,i0:i0+nrep-1) = rep(q)%adjbody(:,:)     
           else                                                        
             base(i0:i0+nrep-1,i0:i0+nrep-1) = rep(q)%adjgrps(:,:)     
           end if                                                      
         end do                                                        
       end do                                                          
!                                                                      
       return                                                          
       end subroutine buildsysbaseadjrep                               
!                                                                       
!======================================================================!
!                                                                       
! BUILDSYSADJREP - BUILD SYStem ADJacency in aggregate REPresentation   
!                                                                       
       subroutine buildsysadjrep(nmax,nagg,imol,mnode,mol,node,        &
                                 msysrep,irepnode,nrepnode,sysadj,     &
                                 matms,posi,coordmat,box,buildadjrep,  &
                                 buildadjmon,domon,dobody)              
!                                                                       
       use systeminf,   only:  mtype,mmon,nmon,imon,mgrpsmon,          &
                               mbodymon,igrpsmon,ibodymon,adjgrps,     &
                               adjbody                                  
!                                                                       
       implicit none                                                   
!                                                                       
! Input/output variables                                                
!                                                                       
       logical,dimension(msysrep,msysrep),intent(out) ::  sysadj       ! Global raw representation adjacency matrix 
       real(kind=4),dimension(3,matms),intent(in)     ::  posi         ! Reduced coordinates for this frame 
       real(kind=4),dimension(:,:),target,intent(in)  ::  coordmat     ! Atomistic coordinates for this frame
       real(kind=4),dimension(3),intent(in)           ::  box          ! Simulation box for this frame 
       integer,dimension(nmax),intent(in)             ::  nagg         ! Number of aggregates per stoichiometry 
       integer,dimension(nmax),intent(in)             ::  imol         ! Initial molecule-list index per stoichiometry
       integer,dimension(mnode),intent(in)            ::  mol          ! Global molecule identifiers in aggregate lists
       integer,dimension(mnode),intent(in)            ::  node         ! Local molecule identifiers in aggregate lists
       integer,dimension(:),intent(in)                ::  irepnode     ! First representation index of each molecule 
       integer,dimension(:),intent(in)                ::  nrepnode     ! Number of representation sites per molecule
       integer,intent(in)                             ::  nmax         ! Maximum aggregate identifier plus overflow
       integer,intent(in)                             ::  mnode        ! Total number of molecules 
       integer,intent(in)                             ::  msysrep      ! Total size of global representation matrix
       integer,intent(in)                             ::  matms        ! Total number of reduced coordinate sites 
       logical,intent(in)                             ::  domon        ! Monomer intramolecular edge flag 
       logical,intent(in)                             ::  dobody       ! Body representation flag 
!                                                                       
! External functions                                                    
!                                                                       
       external                                       ::  buildadjrep  ! Aggregate representation builder
       external                                       ::  buildadjmon  ! Monomer intramolecular edge builder
!                                                                      
! Local variables                                                      
!                                                                      
       logical,dimension(:,:),allocatable             ::  locadj       ! Local aggregate representation adjacency matrix
       integer                                        ::  firstagg     ! First aggregate identifier to process
       integer                                        ::  iagg         ! Stoichiometric aggregate identifier
       integer                                        ::  iocc         ! Aggregate occurrence index
       integer                                        ::  istart       ! Initial molecule-list index for current aggregate
       integer                                        ::  madj         ! Local aggregate matrix dimension
       integer                                        ::  im           ! Local source molecule position
       integer                                        ::  jm           ! Local target molecule position
       integer                                        ::  iloc         ! Initial source index in local matrix 
       integer                                        ::  jloc         ! Initial target index in local matrix 
       integer                                        ::  isys         ! Initial source index in global matrix
       integer                                        ::  jsys         ! Initial target index in global matrix
       integer                                        ::  ni           ! Source molecule representation size 
       integer                                        ::  nj           ! Target molecule representation size
!                                                                      
       sysadj(:,:) = .FALSE.                                           
       firstagg = mtype + 1                                            
!
       if ( domon ) firstagg = 1                                       
!                                                                      
       do iagg = firstagg, nmax-1                                      
!
         if ( nagg(iagg) .eq. 0 ) cycle                                
!
         if ( dobody ) then                                            
           madj = mbodymon(iagg)                                       
         else                                                          
           madj = mgrpsmon(iagg)                                       
         end if   
!         
         allocate(locadj(madj,madj))                                   
!
         istart = imol(iagg)                                           
         do iocc = 1, nagg(iagg)                                       
           if ( dobody ) then                                          
             call buildconfadj(mmon(iagg),node(istart+1:istart+        & 
                  mmon(iagg)),madj,adjbody(iagg)%adj,locadj,matms,     & 
                  posi,coordmat,box,mtype,nmon(:,iagg),imon(:,iagg),   & 
                  ibodymon(:,iagg),buildadjrep,buildadjmon,domon)        
           else                                                        
             call buildconfadj(mmon(iagg),node(istart+1:istart+        & 
                  mmon(iagg)),madj,adjgrps(iagg)%adj,locadj,matms,     & 
                  posi,coordmat,box,mtype,nmon(:,iagg),imon(:,iagg),   & 
                  igrpsmon(:,iagg),buildadjrep,buildadjmon,domon)        
           end if                                                      
           iloc = 1                                                    
           do im = 1, mmon(iagg)                                       
             isys = irepnode(mol(istart+im))                           
             ni = nrepnode(mol(istart+im))                             
             jloc = 1                                                  
             do jm = 1, mmon(iagg)                                     
               jsys = irepnode(mol(istart+jm))                         
               nj = nrepnode(mol(istart+jm))                           
               sysadj(isys:isys+ni-1,jsys:jsys+nj-1) =                 &
                    locadj(iloc:iloc+ni-1,jloc:jloc+nj-1)              
               jloc = jloc + nj                                        
             end do                                                    
             iloc = iloc + ni                                          
           end do                                                      
           istart = istart + mmon(iagg)                                
         end do   
!         
         deallocate(locadj)                                            
!
       end do                                                          
!                                                                      
       return                                                          
       end subroutine buildsysadjrep                                   
!                                                                       
!======================================================================!
!                                                                       
! PRINTSYSADJREP - PRINT SYStem-screened aggregate REPresentation       
!                                                                       
       subroutine printsysadjrep(nmax,nagg,imol,mnode,mol,msysrep,     &
                                 irepnode,nrepnode,sysadj,dobody,      &
                                 directed,domon)
!                                                                       
       use systeminf,   only:  mtype,mmon,mgrpsmon,mbodymon,adjgrps,   &
                               adjbody                                  
       use properties,  only:  num
       use units,       only:  uniadj
       use isotools,    only:  classify_iso_graph
!                                                                       
       implicit none                                                   
!                                                                       
! Input/output variables                                                
!                                                                       
       logical,dimension(msysrep,msysrep),intent(in) ::  sysadj        ! Global screened representation adjacency matrix 
       integer,dimension(nmax),intent(in)            ::  nagg          ! Number of aggregates per stoichiometry 
       integer,dimension(nmax),intent(in)            ::  imol          ! Initial molecule-list index per stoichiometry 
       integer,dimension(mnode),intent(in)           ::  mol           ! Global molecule identifiers in aggregate lists 
       integer,dimension(:),intent(in)               ::  irepnode      ! First representation index of each molecule 
       integer,dimension(:),intent(in)               ::  nrepnode      ! Number of representation sites per molecule 
       integer,intent(in)                            ::  nmax          ! Maximum aggregate identifier plus overflow 
       integer,intent(in)                            ::  mnode         ! Total number of molecules 
       integer,intent(in)                            ::  msysrep       ! Total size of global representation matrix 
       logical,intent(in)                            ::  dobody        ! Body representation flag 
       logical,intent(in)                            ::  directed      ! Directed graph flag
       logical,intent(in)                            ::  domon         ! Monomer intramolecular edge flag
!                                                                     
! Local variables                                                     
!                                                                     
       logical,dimension(:,:),allocatable            ::  locadj        ! Local matrix extracted for one aggregate 
       integer                                       ::  firstagg      ! First aggregate identifier to print 
       integer                                       ::  iagg          ! Stoichiometric aggregate identifier 
       integer                                       ::  iocc          ! Aggregate occurrence index 
       integer                                       ::  istart        ! Initial molecule-list index for current aggregate 
       integer                                       ::  madj          ! Local aggregate matrix dimension 
       integer                                       ::  im            ! Local source molecule position 
       integer                                       ::  jm            ! Local target molecule position 
       integer                                       ::  iloc          ! Initial source index in local matrix 
       integer                                       ::  jloc          ! Initial target index in local matrix 
       integer                                       ::  isys          ! Initial source index in global matrix 
       integer                                       ::  jsys          ! Initial target index in global matrix 
       integer                                       ::  ni            ! Source molecule representation size 
       integer                                       ::  nj            ! Target molecule representation size 
       integer                                       ::  ii            ! Matrix row index
       integer                                       ::  jj            ! Matrix column index
!                                                                        
       firstagg = mtype + 1                                            
       if ( domon ) firstagg = 1                                       
!                                                                        
       do iagg = firstagg, nmax-1                                      
!
         if ( nagg(iagg) .eq. 0 ) cycle                                
         if ( dobody ) then                                            
           madj = mbodymon(iagg)                                       
         else                                                          
           madj = mgrpsmon(iagg)                                       
         end if                                                        
!
         allocate(locadj(madj,madj))                                   
!
         if ( .not. doiso_global ) then
           if ( num(iagg) .eq. 0 ) then
             if ( dobody ) then
               open(unit=uniadj,file=trim(adjbody(iagg)%outp),         &
                    action='write')
               write(uniadj,*) trim(adjbody(iagg)%lab)
               do ii = 1, madj
                 write(uniadj,*) (adjbody(iagg)%adj(ii,jj),jj=1,madj)
               end do
             else
               open(unit=uniadj,file=trim(adjgrps(iagg)%outp),         &
                    action='write')
               write(uniadj,*) trim(adjgrps(iagg)%lab)
               do ii = 1, madj
                 write(uniadj,*) (adjgrps(iagg)%adj(ii,jj),jj=1,madj)
               end do
             end if
           else
             if ( dobody ) then
               open(unit=uniadj,file=trim(adjbody(iagg)%outp),         &
                    position='append',action='write')
             else
               open(unit=uniadj,file=trim(adjgrps(iagg)%outp),         &
                    position='append',action='write')
             end if
           end if
         end if
!
         istart = imol(iagg)                                           
         do iocc = 1, nagg(iagg)                                       
!
           iloc = 1                                                    
!
           do im = 1, mmon(iagg)                                       
             isys = irepnode(mol(istart+im))                           
             ni = nrepnode(mol(istart+im))                             
             jloc = 1                                                  
             do jm = 1, mmon(iagg)                                     
               jsys = irepnode(mol(istart+jm))                         
               nj = nrepnode(mol(istart+jm))                           
               locadj(iloc:iloc+ni-1,jloc:jloc+nj-1) =                 & 
                    sysadj(isys:isys+ni-1,jsys:jsys+nj-1)              
               jloc = jloc + nj                                        
             end do                                                    
             iloc = iloc + ni                                          
           end do 
!           
           if ( doiso_global ) then
             if ( dobody ) then
               call classify_iso_graph(iagg,adjbody(iagg)%lab,        &
                    adjbody(iagg)%outp,locadj,directed)
             else
               call classify_iso_graph(iagg,adjgrps(iagg)%lab,        &
                    adjgrps(iagg)%outp,locadj,.FALSE.)
             end if
           else
             do ii = 1, madj
               write(uniadj,*) (locadj(ii,jj),jj=1,madj)
             end do
           end if
!
           istart = istart + mmon(iagg)                                
!
         end do
!
         if ( .not. doiso_global ) close(uniadj)
!
         deallocate(locadj)   
!         
       end do                                                          
!                                                                        
       return                                                          
       end subroutine printsysadjrep                                   
!                                                                        
!======================================================================! 
!                                                                        
! BUILDCONFADJ - BUILD CONFormational ADJacency matrix
!
       subroutine buildconfadj(isize,node,madj,base,adj,matms,posi,   &
                               coordmat,box,mtype,nnode,inode,        &
                               iadjmon,buildadjrep,buildadjmon,domon)
!
       use systeminf,   only:  coord
!
       implicit none
!
! Input/output variables
!
       logical,dimension(madj,madj),intent(in)       ::  base     !
       logical,dimension(madj,madj),intent(out)      ::  adj      !
       real(kind=4),dimension(3,matms),intent(in)    ::  posi     !
       real(kind=4),dimension(:,:),target,intent(in) ::  coordmat !
       real(kind=4),dimension(3),intent(in)          ::  box      !
       integer,dimension(isize),intent(in)           ::  node     !
       integer,dimension(mtype),intent(in)           ::  nnode    !
       integer,dimension(mtype),intent(in)           ::  inode    !
       integer,dimension(mtype),intent(in)           ::  iadjmon  !
       integer,intent(in)                            ::  isize    !
       integer,intent(in)                            ::  madj     !
       integer,intent(in)                            ::  matms    !
       integer,intent(in)                            ::  mtype    !
       logical,intent(in)                            ::  domon    !
!
! External functions
!
       external                                      ::  buildadjrep
       external                                      ::  buildadjmon
!
       coord => coordmat(:,:)
       adj(:,:) = base(:,:)
!
       if ( domon ) then
         call buildadjmon(isize,node,madj,adj,matms,posi,box,mtype,    &
                          nnode,inode,iadjmon)
       end if
!
       call buildadjrep(isize,node,madj,adj,matms,posi,box,mtype,      &
                        nnode,inode,iadjmon)
!
       return
       end subroutine buildconfadj
!
!======================================================================!
!
! SCRNADJREP - SCReeN ADJacency matrix in aggregate REPresentation
!
       subroutine scrnadjrep(scrn,madj,oldadj,adj,newadj,base,        &
                             directed,screen)
!
       use lengths,     only:  lenschm
!
       implicit none
!
! Input/output variables
!
       character(len=lenschm),intent(in)           ::  scrn      !
       logical,dimension(madj,madj),intent(in)     ::  oldadj    !
       logical,dimension(madj,madj),intent(inout)  ::  adj       !
       logical,dimension(madj,madj),intent(in)     ::  newadj    !
       logical,dimension(madj,madj),intent(in)     ::  base      !
       logical,intent(in)                          ::  directed  !
       integer,intent(in)                          ::  madj      !
!
! External functions
!
       external                                    ::  screen
!
! Local variables
!
       integer                                     ::  i,j       !
!
       if ( directed ) then
         do i = 1, madj
           do j = 1, madj
             select case (trim(scrn))
               case ('collisions')
                 adj(i,j) = adj(i,j) .and. (oldadj(i,j) .or.          &
                                            newadj(i,j))
               case ('oscillations')
                 adj(i,j) = adj(i,j) .or. (oldadj(i,j) .and.          &
                                           newadj(i,j))
               case default
                 adj(i,j) = (oldadj(i,j) .and. adj(i,j)) .or.         &
                            (adj(i,j) .and. newadj(i,j)) .or.         &
                            (oldadj(i,j) .and. newadj(i,j))
             end select
           end do
         end do
       else
         call screen(madj,oldadj,adj,newadj)
       end if
!
       adj(:,:) = adj(:,:) .or. base(:,:)
!
       return
       end subroutine scrnadjrep
!
!======================================================================!
!
! PRINTADJBODY - PRINT ADJacency matrix
!                 in the n-BODY simplified representation
!
       subroutine printadjbody(nmax,nagg,imol,mnode,node,matms,posi,   &
                               box,buildadjbody,buildadjmon,domon)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nmax),intent(in)          ::  nagg     !  Number of aggregates of each size
       integer,dimension(nmax),intent(in)          ::  imol     !
       integer,dimension(mnode),intent(in)         ::  node     !  Molecules identifier
       real(kind=4),dimension(3,matms),intent(in)  ::  posi     !
       real(kind=4),dimension(3),intent(in)        ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       integer,intent(in)                          ::  nmax     !
       integer,intent(in)                          ::  mnode    !
       integer,intent(in)                          ::  matms    !
       logical,intent(in)                          ::  domon    !
!
! External functions
!
       external                                    ::  buildadjbody
       external                                    ::  buildadjmon
!
       call classadjbody(nmax,nagg,imol,mnode,node,matms,posi,box,     &
                         buildadjbody,buildadjmon,domon,.FALSE.)
!
       return
       end subroutine printadjbody
!
!======================================================================!
!
! PRINTADJBODYDIR - PRINT ADJacency matrix
!                    in the directed n-BODY simplified representation
!
       subroutine printadjbodydir(nmax,nagg,imol,mnode,node,matms,     &
                                  posi,box,buildadjbody,buildadjmon,   &
                                  domon)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nmax),intent(in)          ::  nagg     !  Number of aggregates of each size
       integer,dimension(nmax),intent(in)          ::  imol     !
       integer,dimension(mnode),intent(in)         ::  node     !  Molecules identifier
       real(kind=4),dimension(3,matms),intent(in)  ::  posi     !
       real(kind=4),dimension(3),intent(in)        ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       integer,intent(in)                          ::  nmax     !
       integer,intent(in)                          ::  mnode    !
       integer,intent(in)                          ::  matms    !
       logical,intent(in)                          ::  domon    !
!
! External functions
!
       external                                    ::  buildadjbody
       external                                    ::  buildadjmon
!
       call classadjbody(nmax,nagg,imol,mnode,node,matms,posi,box,     &
                         buildadjbody,buildadjmon,domon,.TRUE.)
!
       return
       end subroutine printadjbodydir
!
!======================================================================!
!
! CLASSADJBODY - CLASSify ADJacency matrix
!                in the n-BODY simplified representation
!
       subroutine classadjbody(nmax,nagg,imol,mnode,node,matms,posi,   &
                               box,buildadjbody,buildadjmon,domon,     &
                               directed)
!
       use systeminf,   only:  mtype,mmon,nmon,imon,mbodymon,ibodymon, &
                               adjbody,tmpbody
       use properties,  only:  num
       use units,       only:  uniadj
       use isotools,    only:  classify_iso_graph
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nmax),intent(in)          ::  nagg     !  Number of aggregates of each size
       integer,dimension(nmax),intent(in)          ::  imol     !
       integer,dimension(mnode),intent(in)         ::  node     !  Molecules identifier
       real(kind=4),dimension(3,matms),intent(in)  ::  posi     !
       real(kind=4),dimension(3),intent(in)        ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       integer,intent(in)                          ::  nmax     !
       integer,intent(in)                          ::  mnode    !
       integer,intent(in)                          ::  matms    !
       logical,intent(in)                          ::  domon    !
       logical,intent(in)                          ::  directed !
!
! External functions
!
       external                                    ::  buildadjbody
       external                                    ::  buildadjmon
!
! Local variables
!
       integer                                     ::  firstagg !  Indexes
       integer                                     ::  iagg     !  Indexes
       integer                                     ::  madj     !  Indexes
       integer                                     ::  i,j      !  Indexes
       integer                                     ::  ii,jj    !  Indexes
!
! Classifying adj matrix of the aggregates in the N-body simplified representation
! --------------------------------------------------------------------------------
!
       firstagg = mtype + 1
       if ( domon ) firstagg = 1
!
       do iagg = firstagg, nmax-1
         if ( nagg(iagg) .eq. 0 ) cycle
!
         madj = mbodymon(iagg)
!
         if ( .not. doiso_global ) then
           if ( num(iagg) .eq. 0 ) then
             open(unit=uniadj,file=trim(adjbody(iagg)%outp),           &
                  action='write')
             write(uniadj,*) trim(adjbody(iagg)%lab)
             do ii = 1, madj
               write(uniadj,*) (adjbody(iagg)%adj(ii,jj),jj=1,madj)
             end do
           else
             open(unit=uniadj,file=trim(adjbody(iagg)%outp),           &
                  position='append',action='write')
           end if
         end if
!
         j = imol(iagg)
!
         do i = 1, nagg(iagg)
!
           tmpbody(iagg)%adj(:,:) = adjbody(iagg)%adj(:,:)
!
           if ( domon ) then
             call buildadjmon(mmon(iagg),node(j+1:j+mmon(iagg)),madj,  &
                              tmpbody(iagg)%adj,matms,posi,box,mtype,  &
                              nmon(:,iagg),imon(:,iagg),               &
                              ibodymon(:,iagg))
           end if
!
           call buildadjbody(mmon(iagg),node(j+1:j+mmon(iagg)),madj,   &
                             tmpbody(iagg)%adj,matms,posi,box,mtype,   &
                             nmon(:,iagg),imon(:,iagg),ibodymon(:,iagg))
!
           j = j + mmon(iagg)
!
           if ( doiso_global ) then
             call classify_iso_graph(iagg,adjbody(iagg)%lab,          &
                  adjbody(iagg)%outp,tmpbody(iagg)%adj,directed)
           else
             do ii = 1, madj
               write(uniadj,*) (tmpbody(iagg)%adj(ii,jj),jj=1,madj)
             end do
           end if
!
         end do
!
         if ( .not. doiso_global ) close(uniadj)
!
       end do
!
       return
       end subroutine classadjbody
!
!======================================================================!
!
! PRINTADJGRPS - PRINT ADJacency matrix
!                 in the GRouPS representation
!
       subroutine printadjgrps(nmax,nagg,imol,mnode,node,matms,posi,   &
                               box,buildadjgrps,buildadjmon,domon)
!
       use systeminf,   only:  mtype,mmon,nmon,imon,mgrpsmon,igrpsmon, &
                               adjgrps,tmpgrps
       use properties,  only:  num
       use units,       only:  uniadj
       use isotools,    only:  classify_iso_graph
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nmax),intent(in)          ::  nagg    !  Number of aggregates of each size
       integer,dimension(nmax),intent(in)          ::  imol    !
       integer,dimension(mnode),intent(in)         ::  node    !  Molecules identifier
       real(kind=4),dimension(3,matms),intent(in)  ::  posi    !
       real(kind=4),dimension(3),intent(in)        ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       integer,intent(in)                          ::  nmax    !
       integer,intent(in)                          ::  mnode   !
       integer,intent(in)                          ::  matms   !
       logical,intent(in)                          ::  domon   !
!
! External functions
!
       external                                    ::  buildadjgrps
       external                                    ::  buildadjmon
!
! Local variables
!
       integer                                     ::  firstagg !  Indexes
       integer                                     ::  iagg    !  Indexes
       integer                                     ::  madj    !  Indexes
       integer                                     ::  i,j     !  Indexes
       integer                                     ::  ii,jj   !  Indexes
!
! Classifying adj matrix of the aggregates in the groups representation
! --------------------------------------------------------------------
!
       firstagg = mtype + 1
       if ( domon ) firstagg = 1
!
       do iagg = firstagg, nmax-1
         if ( nagg(iagg) .eq. 0 ) cycle
!
         madj = mgrpsmon(iagg)
!
         if ( .not. doiso_global ) then
           if ( num(iagg) .eq. 0 ) then
             open(unit=uniadj,file=trim(adjgrps(iagg)%outp),           &
                  action='write')
             write(uniadj,*) trim(adjgrps(iagg)%lab)
             do ii = 1, madj
               write(uniadj,*) (adjgrps(iagg)%adj(ii,jj),jj=1,madj)
             end do
           else
             open(unit=uniadj,file=trim(adjgrps(iagg)%outp),           &
                  position='append',action='write')
           end if
         end if
!
         j = imol(iagg)
!
         do i = 1, nagg(iagg)
!
           tmpgrps(iagg)%adj(:,:) = adjgrps(iagg)%adj(:,:)
!
           if ( domon ) then
             call buildadjmon(mmon(iagg),node(j+1:j+mmon(iagg)),madj,  &
                              tmpgrps(iagg)%adj,matms,posi,box,mtype,  &
                              nmon(:,iagg),imon(:,iagg),igrpsmon(:,iagg))
           end if
!
           call buildadjgrps(mmon(iagg),node(j+1:j+mmon(iagg)),madj,   &
                             tmpgrps(iagg)%adj,matms,posi,box,mtype,   &
                             nmon(:,iagg),imon(:,iagg),igrpsmon(:,iagg))
!
           j = j + mmon(iagg)
!
           if ( doiso_global ) then
             call classify_iso_graph(iagg,adjgrps(iagg)%lab,          &
                  adjgrps(iagg)%outp,tmpgrps(iagg)%adj,.FALSE.)
           else
             do ii = 1, madj
               write(uniadj,*) (tmpgrps(iagg)%adj(ii,jj),jj=1,madj)
             end do
           end if
!
         end do
!
         if ( .not. doiso_global ) close(uniadj)
!
       end do
!
       return
       end subroutine printadjgrps
!
!======================================================================!
!
       end module aggtools
!
!======================================================================!
