!======================================================================!
!
       program aggregate
!
       use xdr, only: xtcfile
!
       use input_section
       use linkedlist
       use screening
       use datatypes
       use timings
       use utils
!
       implicit none
!
! System information
!
       type(groinp)                               ::  sys      !  Monomer information
       type(xtcfile)                              ::  xtcf     !  Trajectory information
!
! Input/output files
!
       character(len=leninp)                      ::  traj     !  Trajectory file name
       character(len=leninp)                      ::  conf     !  Configuration file name
       character(len=leninp)                      ::  inp      !  General input file name
       character(len=lenout)                      ::  outp     !  Populations file name
!
! Interaction criteria information
! 
       real(kind=8),dimension(:,:),allocatable    ::  thr      !  Distance threshold
       real(kind=8),dimension(:,:),allocatable    ::  thr2     !  Distance threshold
!
! Neighbour list information
!
       real(kind=4),dimension(:,:),allocatable    ::  drnei    !  Atomic displacements
       real(kind=8)                               ::  neidis   !  Screening distance
       real(kind=8)                               ::  cutdis   !  Cutoff distance
       real(kind=8)                               ::  maxdis   !  Maximum displacement
       integer,dimension(:),allocatable           ::  head     !  Starting neighbours in each subcell
       integer,dimension(:),allocatable           ::  list     !  Index of atoms
       integer,dimension(:),allocatable           ::  cell     !  Neighbouring cells
       integer,dimension(3)                       ::  ncell    !  Unit cells per direction
       integer                                    ::  nhead    !
       logical                                    ::  donei    !
!
! Average properties
!
       real(kind=8),dimension(:,:,:),allocatable  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(:),allocatable      ::  pop      !  Populations
       real(kind=8),dimension(:),allocatable      ::  conc     !  Concentrations
       real(kind=8),dimension(:),allocatable      ::  frac     !  Molar fractions
       real(kind=4),dimension(3)                  ::  box      !
       real(kind=4),dimension(3)                  ::  oldbox   !
       real(kind=8)                               ::  cin      !  Stechiometric concentration
       real(kind=8)                               ::  volu     !  Simulation box volume
!
! Aggregates information in the molecule-based representation
!
       logical,dimension(:,:),allocatable         ::  adj      !  Adjacency matrix
       integer,dimension(:),allocatable           ::  mol      !  Molecules identifier
       integer,dimension(:),allocatable           ::  tag      !  Aggregates identifier
       integer,dimension(:),allocatable           ::  agg      !  Aggregates size 
       integer,dimension(:),allocatable           ::  oldmol   !  Molecules identifier
       integer,dimension(:),allocatable           ::  oldtag   !  Aggregates identifier
       integer,dimension(:),allocatable           ::  oldagg   !  Aggregates size 
       integer,dimension(:),allocatable           ::  newmol   !  Molecules identifier
       integer,dimension(:),allocatable           ::  newtag   !  Aggregates identifier
       integer,dimension(:),allocatable           ::  newagg   !  Aggregates size 
       integer,dimension(:),allocatable           ::  nagg     !  Number of aggregates of each size
       integer,dimension(:),allocatable           ::  iagg     !  
       integer,dimension(:),allocatable           ::  nmol     !  
       integer,dimension(:),allocatable           ::  imol     !  
       integer,dimension(:),allocatable           ::  oldnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable           ::  oldiagg  !  
       integer,dimension(:),allocatable           ::  oldnmol  !  
       integer,dimension(:),allocatable           ::  oldimol  ! 
       integer,dimension(:),allocatable           ::  newnagg  !  Number of aggregates of each size
       integer,dimension(:),allocatable           ::  newiagg  !  
       integer,dimension(:),allocatable           ::  newnmol  !  
       integer,dimension(:),allocatable           ::  newimol  ! 
       integer,dimension(:),allocatable           ::  wasmap   !
       integer,dimension(:),allocatable           ::  willmap  !
       integer,dimension(:),allocatable           ::  molhead  !
       integer,dimension(:),allocatable           ::  mollist  !
       integer,dimension(:),allocatable           ::  molprev  !
       integer,dimension(:),allocatable           ::  molnext  !
       integer,dimension(:),allocatable           ::  oldhead  !
       integer,dimension(:),allocatable           ::  oldlist  !
       integer,dimension(:),allocatable           ::  oldprev  !
       integer,dimension(:),allocatable           ::  oldnext  !
       integer                                    ::  nnode    !  Total number of molecules
       integer                                    ::  msize    !  Maximum aggregate size
       integer                                    ::  nsize    !  
       integer                                    ::  newsize  !  
       integer                                    ::  oldsize  !  
       integer                                    ::  magg     !  Number of chemical species
       integer                                    ::  newmagg  !  Number of chemical species
       integer                                    ::  oldmagg  !  Number of chemical species
       integer                                    ::  oldstep  !
       integer                                    ::  actstep  !
       integer                                    ::  newstep  !
       logical,dimension(:),allocatable           ::  iam      !
       logical,dimension(:),allocatable           ::  iwas     !
       logical,dimension(:),allocatable           ::  iwill    !
       logical,dimension(:),allocatable           ::  imnot    !
       logical,dimension(:),allocatable           ::  iwasnt   !
       logical,dimension(:),allocatable           ::  iwont    !
       logical,dimension(:),allocatable           ::  oldiam   !
       logical,dimension(:),allocatable           ::  oldnot   !
!
! Topological representations information
!
       real(kind=4),dimension(:,:),allocatable    ::  posi     !  Auxiliary coordinates
       real(kind=4),dimension(:,:),allocatable    ::  oldposi  !  Auxiliary coordinates
       character(len=leninp)                      ::  tgrp     !  Groups file title
       character(len=8),dimension(:),allocatable  ::  grptag   !  Names of the groups
       integer,dimension(:),allocatable           ::  body     !  Number of groups in each body
       integer,dimension(:),allocatable           ::  grps     !  Number of subgroups in each group
       integer,dimension(:),allocatable           ::  subg     !  Number of atoms in each subgroup
       integer,dimension(:),allocatable           ::  atms     !  Atoms identifier
       integer,dimension(:),allocatable           ::  nbody    !  Number of groups in each body
       integer,dimension(:),allocatable           ::  ngrps    !  Number of subgroups in each group
       integer,dimension(:),allocatable           ::  nsubg    !  Number of atoms in each subgroup
       integer,dimension(:),allocatable           ::  ibody    !  Number of groups in each body
       integer,dimension(:),allocatable           ::  igrps    !  Number of subgroups in each group
       integer,dimension(:),allocatable           ::  isubg    !  Number of atoms in each subgroup
       integer                                    ::  natms    !
       integer                                    ::  mbody    !  Number of bodies
       integer                                    ::  mgrps    !  Number of groups
       integer                                    ::  msubg    !  Number of subgroups
       integer                                    ::  matms    !
!
! Trajectory control variables
!
       integer                                    ::  nprint   !  Populations printing interval
       integer                                    ::  minstep  !  First step for analysis
       integer                                    ::  maxstep  !  Last step for analysis
       integer                                    ::  nsteps   !  Number of snapshots analyzed
 !
! AnalysisPhenolMD variables
!
       real(kind=8),dimension(9,9,3)              ::  table    !  Number of aggregates of each type
       integer                                    ::  nsolv    !
!
! Declaration of time control variables
!
       integer                                    ::  t1,t2    !  CPU times
       integer                                    ::  t1read   !  Initial reading time
       integer                                    ::  t2read   !  Final reading time
       integer                                    ::  t1scrn   !  Initial screening time
       integer                                    ::  t2scrn   !  Final screening time
!~        integer                                    ::  t1pim    !  Initial PIM analysis time
!~        integer                                    ::  t2pim    !  Final PIM analysis time
!~        integer                                    ::  t1conf   !  Initial conformational analysis time
!~        integer                                    ::  t2conf   !  Final conformational analysis time
!
! Program control flags
!
       logical                                    ::  seed     !  Random seed flag
       logical                                    ::  dopim    !  PIM calculation flag
       logical                                    ::  debug    !  Debug mode
!
! Auxiliary variables
!
       integer                                    ::  lin      !  
       integer                                    ::  lfin     ! 
       integer                                    ::  io       !  Status
       integer                                    ::  i,j,k    !  Indexes
       real(kind=8),parameter                     ::  Na = 6.022140760E+23
!
! Printing header
!
       write(*,*)
!~        write(*,'(5X,76("#"))')
!~        write(*,'(5X,14("#"),3X,A,3X,13("#"))') 'RandConf - Random'//  &
!~                               ' Configurations Generator '
!~        write(*,'(5X,76("#"))')   
!~        write(*,*)
!~        write(*,'(7X,A)') 'Welcome to RandConf, a very simple progr'//  &
!~                          'am for the random generation of'
!~        write(*,'(7X,A)') ' starting configurations of molecular sy'//  &
!~                          'stems for Atomistic Simulations'
!~        write(*,*)
!~        write(*,'(1X,84("-"))')
!~        write(*,*)    
!
       call system_clock(count_max=count_max,count_rate=count_rate)
       call system_clock(t1)  
       call system_clock(t1read)        
!
       call print_start()
!
       lin  = 45
       lfin = 90
!
! Initializing timings
!
       tcpu  = 0.0d0
       tread = 0.0d0
       tadj  = 0.0d0
       tbfs  = 0.0d0
       tsort = 0.0d0
       tscrn = 0.0d0
       tpim  = 0.0d0
       tconf = 0.0d0
!
! Initializing random number generator
!
       call rndmseed(seed)
!
! Reading command line options
!
       call command_line(traj,conf,inp,outp,nprint,minstep,maxstep,    & 
                         msize,neidis,cutdis,maxdis,dopim,seed,debug)
!
!  General settings and fatal errors check
!
       io = index(traj,'.')
       if ( io .eq. 0 )  traj = trim(traj)//'.xtc'
!
       io = index(conf,'.')
       if ( io .eq. 0 ) conf = trim(conf)//'.gro'
!
       if ( outp(len_trim(outp)-3:) .eq. '.dat' )                      &
                                          outp = outp(:len_trim(outp)-4)
!
! Printing summary of the input information
!
       call print_title(6,1,'Input information','.')
       write(*,*)
       call line_str(6,2,'General input file name',lin,':',            &
                     trim(inp),lfin)
       call line_str(6,2,'Trajectory file name',lin,':',               &
                     trim(traj),lfin)
       call line_str(6,2,'Configuration file name',lin,':',            &
                     trim(conf),lfin)
       call line_str(6,2,'Output file name',lin,':',                   &
                     trim(outp),lfin)
       write(*,*)
       call line_dp(6,3,'Screening distance',lin,':','F5.2',           &
                    neidis,lfin)
       call line_int(6,2,'Maximum aggregate size',lin,':','I4',        &
                     msize,lfin)
       call line_int(6,2,'First step for analysis',lin,':','I12',      &
                     minstep,lfin)
       call line_int(6,2,'Last step for analysis',lin,':','I12',       &
                     maxstep,lfin)
       call line_int(6,2,'Printing interval',lin,':','I12',            &
                     nprint,lfin)
       if ( debug ) then
         call line_str(6,2,'Debug mode',lin,':','ON',lfin)
       else
         call line_str(6,2,'Debug mode',lin,':','OFF',lfin)
       end if
       write(*,*)
!
! Processing Gromacs input file
!
       call read_gro(conf,sys)    ! FLAG: Input .gro .top .xyz ...
!
! Opening output files
!
       open(unit=uniout+1,file=trim(outp)//'_pop.dat',action='write')
       open(unit=uniout+2,file=trim(outp)//'_frac.dat',action='write')
       open(unit=uniout+3,file=trim(outp)//'_conc.dat',action='write')
!
! Opening trajectory file
!
       if ( traj(len_trim(traj)-3:) .eq. '.xtc' ) then 
         call xtcf%init(trim(traj))
         call xtcf%read
       else
         write(*,*) 'Incorrect extension!'
         write(*,*)
         call print_end()
       end if
!
! Allocating memory depending on the system size
!
       nnode  = xtcf%NATOMS/sys%nat
!
       box(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
!
       allocate(adj(nnode,nnode))
!
       allocate(thr(sys%nat,sys%nat))
       allocate(thr2(sys%nat,sys%nat))
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
       allocate(molhead(nnode),mollist(nnode))
       allocate(molprev(nnode),molnext(nnode))
!
       allocate(oldhead(nnode),oldlist(nnode))
       allocate(oldprev(nnode),oldnext(nnode))
!
       allocate(nbody(sys%nat),ngrps(sys%nat),nsubg(sys%nat))
       allocate(ibody(sys%nat),igrps(sys%nat),isubg(sys%nat))
       allocate(body(sys%nat),grps(sys%nat),subg(sys%nat))
       allocate(grptag(sys%nat),atms(sys%nat))
!
! Processing general input file
!
       call read_inp(inp,sys%nat,tgrp,grptag,nbody,ngrps,nsubg,atms,   &
                     mbody,mgrps,msubg,matms,thr,thr2)
!
       natms = msubg*nnode
!         
       ibody(:) = 0
       igrps(:) = 0
       isubg(:) = 0
!
       body(:) = 0
       grps(:) = 0
       subg(:) = 0
! Setting up igrps and ngrps arrays
       k = 1
       do i = 1, mgrps
         if ( i .lt. mgrps ) then
           igrps(i+1) = igrps(i) + ngrps(i)
         end if
!
         do j = 1, ngrps(i)
           grps(k) = i
           k = k + 1
         end do
       end do
! Setting up isubg and nsubg arrays
       do i = 1, msubg-1
         isubg(i+1) = isubg(i) + nsubg(i)
       end do
!
       do i = 1, mgrps
         do j = 1, ngrps(i)
           io = igrps(i) + j
           do k = 1, nsubg(io)
             subg(isubg(io)+k) = i 
           end do
         end do
       end do
! Setting up ibody and nbody arrays
       do i = 1, mbody-1
         ibody(i+1) = ibody(i) + nbody(i)
       end do
!
       do i = 1, mbody
         do j = 1, nbody(i)
           io = ibody(i) + j
           do k = 1, ngrps(io)
             body(igrps(io)+k) = i
            end do
         end do
       end do
!
write(*,'(A,20I3)') 'mbody',mbody
write(*,'(A,20I3)') 'nbody',nbody
write(*,'(A,20I3)') 'ibody',ibody
write(*,'(A,20I3)') 'body ',body
write(*,*) 
write(*,'(A,20I3)') 'mgrps',mgrps
write(*,'(A,20I3)') 'ngrps',ngrps
write(*,'(A,20I3)') 'igrps',igrps
write(*,'(A,20I3)') 'grps ',grps
write(*,*) 
write(*,'(A,20I3)') 'msubg',msubg
write(*,'(A,20I3)') 'nsubg',nsubg
write(*,'(A,20I3)') 'isubg',isubg
write(*,'(A,20I3)') 'subg ',subg
write(*,*)       
write(*,'(A,20I3)') 'atms ',atms
write(*,*)       
!
! Allocating memory depending on the topological definition 
!
       allocate(oldposi(3,natms))
       allocate(posi(3,natms))
!
       allocate(pim(mgrps,mgrps,msize-1))
       allocate(pop(msize),conc(msize),frac(msize))
!
! Initializing variables
!
       do i = 1, 3
         ncell(i) = int(xtcf%box(i,i)/neidis)
       end do
!
       nhead = ncell(1)*ncell(2)*ncell(3)
!
       allocate(head(nhead))
       allocate(list(natms))
       allocate(cell(nhead*13))
       allocate(drnei(3,natms))
!
       drnei(:,:) = 100.0d0
!
       neidis = neidis**2
       cutdis = cutdis**2
       maxdis = maxdis**2
!
!~        neidis = maxval(thr2)
!~        neidis = neidis**2
!
       thr(:,:)  = thr(:,:)**2
       thr2(:,:) = thr2(:,:)**2        
!
! Computing the populations of the aggregates
!
       write(*,'(4X,A)') 'Computing the populations of the aggrega'//  & 
                                              'tes along the trajectory'
       write(*,'(4X,A)') 'Please wait, this may take a while...'
       write(*,*)
!
!~!        call setcoord(nnode,msubg,nsubg,isubg,atms,natms,posi,          &
!~!                      xtcf%NATOMS,xtcf%pos,sys%nat,box)
!
!~        call setneig(nhead,cell,ncell)
!
!~!        call setlinklist(natms,posi,drnei,list,nhead,head,ncell,box,xtcf%STEP) 
!
!~!        call buildadjmollink(nnode,adj,neidis,msubg,mgrps,              &
!~!                             thr(:mgrps,:mgrps),ngrps(:mgrps),          &
!~!                             igrps(:mgrps),grps(:msubg),natms,          &
!~!                             posi,list,nhead,head,cell,ncell,box)  
!
!~!        call xtcf%read
!
       donei = .TRUE.
!
       nsteps  = 0
!
       pim(:,:,:)   = 0.0d0
       table(:,:,:) = 0.0d0
!
       pop(:)  = 0.0d0
       conc(:) = 0.0d0
       frac(:) = 0.0d0
!
       cin     = 0.0d0
       volu    = 0.0d0
!
       nsolv = 10000     ! FLAG: change this value
!
! Reading the first old configuration
!
       do while ( xtcf%STEP .lt. minstep )
         call xtcf%read
       end do
!
       oldbox  = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
       oldstep = xtcf%STEP
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       call blockdiag(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,       &
                      sys%nat,thr2,atms,mgrps,grps,ngrps,igrps,        &
                      msubg,subg,nsubg,isubg,oldmol,oldtag,oldagg,     &
                      oldsize,oldnagg,oldiagg,oldnmol,oldimol,oldmagg, &
                      oldbox,neidis,debug)
!
!~        call array2list(nnode,oldmol,oldagg,oldsize,oldnagg,oldimol,    &
!~                        molhead,molprev,mollist,molnext)
!
!~        write(*,'(2X,A,X,I6)')     'Total number of entities : ',oldmagg
!~        write(*,*)
!~        write(*,'(2X,A,20(X,I6))') 'Aggregates of each type  : ',oldnagg(oldsize)
!~        write(*,'(2X,A,20(X,I6))') '  Accumulation           : ',oldiagg(oldsize)
!~        write(*,*)
!~        write(*,'(2X,A,20(X,I6))') 'Molecules of each type   : ',oldnmol(oldsize)
!~        write(*,'(2X,A,20(X,I6))') '  Accumulation           : ',oldimol(oldsize)
!~        write(*,*)
!~ !
!~        call print_info(oldnagg(1),nnode-oldnagg(1),              &
!~                        oldmol(oldnagg(1)+1:),                    &
!~                        oldtag(oldnagg(1)+1:),                    &  
!~                        oldagg(oldnagg(1)+1:),'mol','tag','agg')
!~ !
!~        call print_test(nnode,0,nnode,                                  &
!~                        molhead,                                        &
!~                        mollist,                                        &    
!~                        molprev,                                        &
!~                        molnext,                                        &
!~                        'head','list','prev','next')
!~ !
!~ stop 'Debugging'
!
! Reading the first configuration
!
       call system_clock(t1read)
!
       call xtcf%read
!
       do while ( mod(xtcf%STEP-minstep,nprint) .ne. 0 ) 
         call xtcf%read
         if ( (xtcf%STAT.ne.0) .or. (xtcf%STEP.gt.maxstep) ) then
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
       call blockdiag(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,       &
                      sys%nat,thr2,atms,mgrps,grps,ngrps,igrps,        &
                      msubg,subg,nsubg,isubg,mol,tag,agg,nsize,        &
                      nagg,iagg,nmol,imol,magg,box,neidis,debug)
!
       call system_clock(t1scrn)
!
       call scrnblock(nnode,oldsize,oldmol,oldagg,nsize,mol,agg,       &
                      oldnagg,oldiagg,oldnmol,oldimol,nagg,iagg,       &
                      nmol,imol,iwas,iwasnt,wasmap) 
!
       oldiam(:) = .FALSE.
       oldnot(:) = .TRUE.
!
       do i = nagg(1)+1, magg
         if ( wasmap(i) .ne. 0 ) then
           oldiam(wasmap(i)) = .TRUE.
           oldnot(wasmap(i)) = .FALSE.
         end if
       end do
!
       call system_clock(t2scrn)
!
       tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
! Reading the first new configuration
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
! Analyzing frames in the inverval (minstep,maxstep)
!
       do while ( (xtcf%STAT.eq.0) .and. (xtcf%STEP.le.maxstep) )
         if ( mod(xtcf%STEP-minstep,nprint) .eq. 0 ) then
!
           box(:) = (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/)
           newstep = xtcf%STEP
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 
! 
           nsteps = nsteps + 1
!
! Block diagonalization of the adjacency matrix of the new configuration
!
           call blockdiag(nnode,adj,natms,posi,xtcf%NATOMS,xtcf%pos,   &
                          sys%nat,thr2,atms,mgrps,grps,ngrps,igrps,    &
                          msubg,subg,nsubg,isubg,newmol,newtag,newagg, &
                          newsize,newnagg,newiagg,newnmol,newimol,     &
                          newmagg,box,neidis,debug)
!
! First screening: if I was and I will then I am
!
           call system_clock(t1scrn)
!
           call scrnblock(nnode,newsize,newmol,newagg,nsize,mol,agg,   &
                          newnagg,newiagg,newnmol,newimol,nagg,iagg,   &
                          nmol,imol,iwill,iwont,willmap)
!
           iam(:) = .FALSE.
           imnot(:) = .TRUE.
!
           do i = nagg(1)+1, magg
             if ( iwas(i) .or. iwill(i) ) then
               iam(i)   = .TRUE.
               imnot(i) = .FALSE.
             end if
           end do
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
             write(*,'(4X,A)') 'OLD I am'
             do i = 2, oldsize
               if ( (oldiagg(i)+oldnagg(i)) .lt. (oldiagg(i)+1) ) exit
               write(*,'(I3,2(1X,A,1X,I5),1X,A,20(5X,L1))') i,'from',oldiagg(i)+1,   &
               'to',oldiagg(i)+oldnagg(i),':',oldiam(oldiagg(i)+1:oldiagg(i)+oldnagg(i))
             end do
             write(*,*)
!
             write(*,'(4X,A)') 'I was'
             do i = 2, nsize
               if ( (iagg(i)+nagg(i)) .lt. (iagg(i)+1) ) exit
               write(*,'(I3,2(1X,A,1X,I5),1X,A,20(5X,L1))') i,'from',iagg(i)+1,   &
               'to',iagg(i)+nagg(i),':',iwas(iagg(i)+1:iagg(i)+nagg(i))
               write(*,'(25X,20(1X,I5))') wasmap(iagg(i)+1:iagg(i)+nagg(i))
             end do
             write(*,*)
!
             write(*,'(4X,A)') 'I am'
             do i = 2, nsize
               if ( (iagg(i)+nagg(i)) .lt. (iagg(i)+1) ) exit
               write(*,'(I3,2(1X,A,1X,I5),1X,A,20(5X,L1))') i,'from',iagg(i)+1,   &
                 'to',iagg(i)+nagg(i),':',iam(iagg(i)+1:iagg(i)+nagg(i))
             end do
             write(*,*)
!
             write(*,'(4X,A)') 'I will'
             do i = 2, nsize
               if ( (iagg(i)+nagg(i)) .lt. (iagg(i)+1) ) exit
               write(*,'(I3,2(1X,A,1X,I5),1X,A,20(5X,L1))') i,'from',iagg(i)+1,   &
              'to',iagg(i)+nagg(i),':',iwill(iagg(i)+1:iagg(i)+nagg(i))
               write(*,'(25X,20(1X,I5))') willmap(iagg(i)+1:iagg(i)+nagg(i))
             end do
             write(*,*)
!
!~              call print_coord(xtcf,sys,outp,msize,                     &
!~                               newnagg,nnode,newmol,newagg)
           end if
!
! Changing from the array representation to the linked-list representa-
!   tion (their array representation can change after the screening)
!
           call array2list(nnode,oldmol,oldagg,oldsize,oldnagg,oldimol,&
                           oldhead,oldprev,oldlist,oldnext)
!
           call array2list(nnode,mol,agg,nsize,nagg,imol,              &
                           molhead,molprev,mollist,molnext)
!
! Second screening: although i was or i wasnt i can still be if i wont
!
           call scrnalrd(nnode,nsize,mol,agg,oldsize,oldmol,oldagg,    &
                         nagg,iagg,nmol,imol,oldnagg,oldiagg,oldnmol,  &
                         oldimol,iwas,iwasnt,oldiam,oldnot)
!
           if ( debug ) then
             write(*,'(1X,A)') 'SECOND SCREENING'
             write(*,'(1X,15("."))') 
             write(*,*)
           end if
!
! Third screening: if im not maybe i was already
!

!
           if ( debug ) then
             write(*,'(1X,A)') 'THIRD SCREENING'
             write(*,'(1X,15("."))') 
             write(*,*)
           end if
!
           call system_clock(t2scrn)
!
           tscrn = tscrn + dble(t2scrn-t1scrn)/dble(count_rate) 
!
! Printing the population of every aggregate
!
           call system_clock(t1read)
!
           do i = 1, msize
             pop(i)  = pop(i)  + real(oldnagg(i))/oldmagg*100
             frac(i) = frac(i) + real(oldnagg(i))/(oldmagg+nsolv)
             conc(i) = conc(i) + real(oldnagg(i))/oldbox(1)**3 
           end do     
!
           if ( oldsize .gt. msize ) then
             do i = msize+1, oldsize
               pop(msize)  = pop(msize)  + real(oldnagg(i))/oldmagg*100
               frac(msize) = frac(msize)                               &
                                      + real(oldnagg(i))/(oldmagg+nsolv)
               conc(msize) = conc(msize) + real(oldnagg(i))/oldbox(1)**3  
             end do
           end if                   
!
           write(uniout+1,'(I10,100(X,F12.8))') oldstep,               &
                                       real(oldnagg(:msize))/oldmagg*100
           write(uniout+2,'(I10,100(X,F12.10))') oldstep,              &
                                   real(oldnagg(:msize))/(oldmagg+nsolv)
           write(uniout+3,'(I10,100(X,F12.10))') oldstep,              &
                         real(oldnagg(:msize))/oldbox(1)**3/(Na*1.0E-24)                            
!
           cin = cin + real(nnode)/oldbox(1)**3
!
           volu = volu + oldbox(1)**3
!
! Printing summary of the results
! 

!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Adding up the pairwise interaction matrix of the current snapshot
! 
!~            if ( dopim ) then
!~              call system_clock(t1pim)     
!~ !
!~              call build_pim(ngrps,grps,nsubg,subg,sys%nat,atms,thr,   &
!~                             nnode,mol,agg,msize,nagg,xtcf%NATOMS,   &
!~                             xtcf%pos,(/xtcf%box(1,1),xtcf%box(2,2),    &
!~                             xtcf%box(3,3)/),pim,debug)
!~ !
!~              call system_clock(t2pim)     
!~ !
!~              tpim = tpim + dble(t2pim-t1pim)/dble(count_rate)   
!~            end if      
!
! Analyzing aggregates by their connectivity
!~ !
!~            call analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,    &
!~                             sys%nat,atms,thr,nnode,adj,mol,agg,     & 
!~                             msize,nagg,xtcf%NATOMS,xtcf%pos,          &
!~                             (/xtcf%box(1,1),xtcf%box(2,2),             &
!~                               xtcf%box(3,3)/),posi,sys%mass,          &
!~                             sys%atname,xtcf%STEP,outp,debug)
!
           call system_clock(t1scrn)
!
           oldbox(:) = box(:)
!
           oldmol(:) = mol(:)
           oldagg(:) = agg(:)
           oldtag(:) = tag(:)
!
           oldnagg(:) = nagg(:)
           oldiagg(:) = iagg(:)
           oldnmol(:) = nmol(:)
           oldimol(:) = imol(:)
!
           oldmagg = magg
!
           oldiam(:) = iam(:)
           oldnot(:) = imnot(:)
!
           wasmap(:) = 0
           iwas(:)   = .FALSE.
           iwasnt(:) = .TRUE.
!
           do i = nagg(1)+1, magg
             if ( willmap(i) .ne. 0 ) then
               wasmap(willmap(i)) = i
               iwas(willmap(i))   = .TRUE.
               iwasnt(willmap(i)) = .FALSE.
             end if
           end do
!
           mol(:) = newmol(:)
           agg(:) = newagg(:)
           tag(:) = newtag(:)
!
           nagg(:) = newnagg(:)
           iagg(:) = newiagg(:)
           nmol(:) = newnmol(:)
           imol(:) = newimol(:)
!
           magg = newmagg
!
           oldstep = actstep
           actstep = newstep
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
! Close the file
       call xtcf%close
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
! Averaging pairwise interaction matrix
!  
!~        if ( dopim ) then
!~          pim(:,:,:) = pim(:,:,:)/nsteps ! FLAG: not need this operation
!~ !
!~          do i = 1, msize-1
!~            dpaux = 0.0d0
!~            do j = 1, ngrps
!~              do k = j, ngrps
!~                dpaux = dpaux + pim(k,j,i)
!~              end do
!~            end do
!~            pim(:,:,i) = pim(:,:,i)/dpaux*100
!~          end do
!~        end if
!
! Averaging populations
!
       pop(:) = pop(:)/nsteps
!
       frac(:) = frac(:)/nsteps
!
       conc(:) = conc(:)/nsteps/(Na*1.0E-24)
!
       cin = cin/nsteps/(Na*1.0E-24)
!
       volu = volu/nsteps
!~ !
!~        table(:,:,:) = table(:,:,:)/nsteps
!~ !
!~        do i = 1, 3
!~          dpaux = 0.0d0
!~          do j = 1, 9
!~            do k = 1, 9
!~              dpaux = dpaux + table(k,j,i)
!~            end do
!~          end do
!~          table(:,:,i) = table(:,:,i)/dpaux*100
!~        end do
!
! Printing summary of the results
!
       write(*,'(1X,A)') 'Output information'
       write(*,'(1X,18("-"))')
       write(*,*)
       write(*,'(1X,A,3X,I11)')      'Number of frames analyzed : ',   &
                                                                  nsteps
       write(*,*)
       write(*,'(1X,A,20(X,F12.8))') 'Initial concentration     : ',cin
       write(*,'(1X,A,20(X,F12.6))') 'Average volume            : ',volu
       write(*,*)
       write(*,'(1X,A,100(X,F6.2))')  'Global populations        : ',  &
                                                                  pop(:)
       write(*,'(1X,A,100(X,D12.6))') 'Global fractions          : ',  &
                                                                 frac(:)
       write(*,'(1X,A,100(X,D12.6))') 'Global concentrations     : ',  &
                                                                 conc(:)
       write(*,*) 
!
!~        do i = 1, 3
!~          write(*,'(3X,A,X,I3)') 'Populations of the aggregates belonging to type',i+1
!~          write(*,'(3X,51("-"))') 
!~          write(*,'(8X,3(1X,A))') 'IDXOH','IDXPh','Population'
!~          do j = 1, 9
!~            do k = 1, 9
!~              if ( table(k,j,i) .gt. 1e-3 ) then
!~                write(*,'(8X,2(1X,I4),1X,F6.2)') k,j,table(k,j,i)
!~              end if
!~            end do
!~          end do
!~          write(*,*)
!~        end do
!~ !
!~        if ( dopim ) then
!~          do i = 1, msize-1
!~            write(*,'(3X,A,X,I3)') 'Printing PIM for aggregates bel'//  &
!~                                                     'onging to type',i+1
!~            write(*,'(3X,50("-"))') 
!~            write(*,'(8X,20(X,A8))') (adjustr(grptag(j)),j=1,ngrps)
!~            do j = 1, ngrps
!~              write(*,'(A8,20(X,F8.2))') adjustr(grptag(j)),            &
!~                                                       (pim(k,j,i),k=1,j)
!~            end do
!~            write(*,*)
!~          end do
!~        end if
!
! Deallocate memory
!
       deallocate(oldposi)
       deallocate(posi)
!
       deallocate(drnei)
       deallocate(head)
       deallocate(list)
       deallocate(cell)
!
       deallocate(thr,thr2)
!
       deallocate(adj,mol,tag,agg)
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
       deallocate(molhead,mollist,molprev,molnext)
!
       deallocate(oldhead,oldlist,oldprev,oldnext)
!
       deallocate(nbody,ngrps,nsubg)
       deallocate(ibody,igrps,isubg)
       deallocate(body,grps,subg)
       deallocate(grptag,atms)
!
       deallocate(pop,conc,frac)
       deallocate(pim)
!
       close(uniout+1)
       close(uniout+2)
       close(uniout+3)
!
! Printing timings
!
       call system_clock(t2)            ! FLAG: 3 seconds lost
!
       tcpu = dble(t2-t1)/dble(count_rate)
!
       call print_time(6,1,'Total CPU time',31,tcpu) 
       write(*,'(1X,56("-"))')
!
       call print_time(6,1,'Total reading time',31,tread)
       call print_time(6,1,'Total building time',31,tadj)
       call print_time(6,1,'Total BFS time',31,tbfs)
       call print_time(6,1,'Total sorting time',31,tsort)
       call print_time(6,1,'Total screening time',31,tscrn)
!
       if ( dopim ) call print_time(6,1,'Total PIM time',31,tpim)
!
       write(*,*)
!
! Printing finishing date 
!    
       call print_end()
!
       end program aggregate
!
!======================================================================!
!
       subroutine command_line(traj,conf,inp,outp,nprint,minstep,      &
                               maxstep,msize,neidis,cutdis,maxdis,    &
                               dopim,seed,debug)
!
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(out)         ::  traj     !  Trajectory file name
       character(len=lenout),intent(out)         ::  outp     !  Populations file name
       character(len=leninp),intent(out)         ::  inp      !  Groups file name
       character(len=leninp),intent(out)         ::  conf     !  Structure file name
       logical,intent(out)                       ::  seed     !  Random seed flag
       logical,intent(out)                       ::  dopim    !  PIM calculation flag
       real(kind=8),intent(out)                  ::  neidis   !  Screening distance
       real(kind=8),intent(out)                  ::  cutdis   !  Cutoff distance
       real(kind=8),intent(out)                  ::  maxdis   !  Maximum displacement
       integer,intent(out)                       ::  msize    !  Maximum aggregate size
       integer,intent(out)                       ::  nprint   !  Populations printing steps interval
       integer,intent(out)                       ::  minstep  !  First step for analysis
       integer,intent(out)                       ::  maxstep  !  Last step for analysis
       logical,intent(out)                       ::  debug    !  Debug mode
!
! Local variables
!
       character(len=lencmd)                     ::  cmd     !  Command executed
       character(len=lenarg)                     ::  code    !  Executable name
       character(len=lenarg)                     ::  arg     !  Argument read
       character(len=lenarg)                     ::  next    !  Next argument to be read
       integer                                   ::  io      !  Status
       integer                                   ::  i       !  Index
!
! Setting defaults
!
       inp     = 'aggregate.inp'
       traj    = 'md.xtc'
       conf    = 'conf.gro'
       outp    = 'md.dat'
       nprint  = 1
       minstep = 0
       maxstep = 999999999
       msize  = 10
       neidis  = 1.5d0
       cutdis  = 0.5d0
       seed    = .FALSE.
       dopim   = .TRUE.
       debug   = .FALSE.
!
! Reading command line options
! ----------------------------
!
       call get_command_argument(0,code)
       call get_command(cmd)
! Checking if any argument has been introduced
       if ( command_argument_count().eq. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  No argument introduced on'//    &
                              ' command-line'
         write(*,*)
         write(*,'(3X,A)')    'Please, to request help execute'
             write(*,*)
         write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
! Reading command line options
       i = 1
       do
         call get_command_argument(i,arg)
         if ( len_trim(arg) == 0 ) exit
         i = i+1
         select case ( arg )
           case ('-f','-file','--file','-i','-inp','--inp','--input') 
             call get_command_argument(i,inp,status=io)
             call check_arg(inp,io,arg,cmd)
             i = i + 1
           case ('-t','-traj','--traj','--trajectory') 
             call get_command_argument(i,traj,status=io)
             call check_arg(traj,io,arg,cmd)
             io = index(traj,'.')
             if ( io .eq. 0 ) then
               outp = trim(traj)//'.dat'
             else 
               outp = traj(:len_trim(traj)-4)//'.dat'
             end if             
             i = i + 1
           case ('-c','-conf','--conf','--configuration')
             call get_command_argument(i,conf,status=io)
             call check_arg(conf,io,arg,cmd)
             i = i + 1
           case ('-o','-p','-outp','-pop','--pop','--outp',            &
                 '--population','--populations','--output')
             call get_command_argument(i,outp,status=io)
             call check_arg(outp,io,arg,cmd)
             i = i + 1
           case ('-d','-nd','-neidis','--neidis','--neighbour-distance')
             call get_command_argument(i,next,status=io)
             read(next,*) neidis
             maxdis = neidis - cutdis
             i = i + 1
           case ('-cd','-cutdis','--cutdis','--cutoff-distance')
             call get_command_argument(i,next,status=io)
             read(next,*) cutdis
             maxdis = neidis - cutdis

             i = i + 1
           case ('-md','-maxdis','--maxdis','--maximum-displacement')
             call get_command_argument(i,next,status=io)
             read(next,*) maxdis
             i = i + 1
           case ('-m','-msize','--msize','--maximum-size')
             call get_command_argument(i,next,status=io)
             read(next,*) msize
             i = i + 1
           case ('-n','-nprint','--nprint','--print-steps')
             call get_command_argument(i,next,status=io)
             read(next,*) nprint
             i = i + 1
           case ('-min','-minstep','--minstep','--minimum-step')
             call get_command_argument(i,next,status=io)
             read(next,*) minstep
             i = i + 1
           case ('-max','-maxstep','--maxstep','--maximum-step')
             call get_command_argument(i,next,status=io)
             read(next,*) maxstep     ! FLAG: if value not introduced print error
             i = i + 1
           case ('-pim','--pim','--do-pim')
             dopim = .TRUE.
           case ('-nopim','--nopim','--no-pim')
             dopim = .FALSE.
           case ('-seed','--random-seed')
             seed = .TRUE.
           case ('-v','--debug','--verbose')
             debug = .TRUE.
           case ('-h','-help','--help')
             call print_help()
             call print_end()
           case default
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)')    'ERROR:  Unknown statements from'//  &
                                  ' command line'
             write(*,*)
             write(*,'(4X,A)')     trim(cmd)
             write(*,*)
             write(*,'(3X,2(A))') 'Unrecognised command-line option'// &
                                  '  :  ', arg
             write(*,'(3X,A)')    'Please, to request help execute'
             write(*,*)
             write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
             write(*,'(2X,68("="))')
             write(*,*)
             call print_end()
         end select
       end do
!
       return
       end subroutine command_line
!
!======================================================================!
!
       subroutine print_help()
!
       implicit none
!
       write(*,'(1X,A)') 'Command-line options'
       write(*,'(1X,20("-"))')
       write(*,*)
       write(*,'(2X,A)') '-h,--help             Print usage inform'//  &
                                                         'tion and exit'
       write(*,'(2X,A)') '-f,--file             Input file name'
       write(*,'(2X,A)') '-t,--trajectory       Trajectory file name'
       write(*,'(2X,A)') '-c,--configuration    Configuration file name'
       write(*,'(2X,A)') '-p,--populations      Populations file name'
       write(*,'(2X,A)') '-n,--nprint           Printing steps interval'
       write(*,'(2X,A)') '-min,--minimum-step   First step to be a'//  &
                                                               'nalysed'
       write(*,'(2X,A)') '-max,--maximum-step   Last step to be an'//  &
                                                                'alysed'
       write(*,'(2X,A)') '-m,--msize            Maximum aggregate size'
       write(*,'(2X,A)') '-d,--neidis           Neighbour list cutoff'
       write(*,'(2X,A)') '-[no]pim              Compute pairwise i'//  &
                                                     'nteraction matrix'
       write(*,'(2X,A)') '-v,--verbose          Debug mode'
       write(*,*)
!
       return
       end subroutine print_help
!
!======================================================================!
!
! SETCOORD - SET COORDinates based on the subgroup-based representation
!
       subroutine setcoord(nnode,msubg,nsubg,isubg,atms,natms,fcoord,  &
                           natsys,rcoord,natmol,box)
!
       use geometry,   only: sminimgvec,                               &
                             scenvec
!
       implicit none
!
! Input/output variables
!
       real(kind=4),dimension(3,natsys),intent(in)  ::  rcoord   !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3,natms),intent(out)  ::  fcoord   !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)         ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(natmol),intent(in)         ::  nsubg    !  
       integer,dimension(natmol),intent(in)         ::  isubg    !   
       integer,dimension(natmol),intent(in)         ::  atms     !  Atoms identifier
       integer,intent(in)                           ::  msubg    !  Number of subgroups
       integer,intent(in)                           ::  nnode    !  Number of residues
       integer,intent(in)                           ::  natsys   !  Total number of atoms
       integer,intent(in)                           ::  natms    ! 
       integer,intent(in)                           ::  natmol   !  Atoms per residue
!
! Local variables
!
       real(kind=4),dimension(3,natmol)             ::  atcoord  !  Subgroup coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                    ::  cofm     !  Center of mass coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                    ::  svaux    !  Auxiliary single precision vector
       integer                                      ::  iinode   !
       integer                                      ::  innode   !
       integer                                      ::  jnnode   !
       integer                                      ::  iisubg   !
       integer                                      ::  insubg   !
       integer                                      ::  i,j      !
!
! Saving coordinates based on the subgroup-based representation
! -------------------------------------------------------------
!
       do iinode = 1, nnode
         innode = (iinode-1)*natmol
         jnnode = (iinode-1)*msubg
         do insubg = 1, msubg
!
           j = jnnode + insubg
!
           if ( nsubg(insubg) .gt. 1 ) then
!
             svaux(:) = rcoord(:,innode+atms(isubg(insubg)+1))
             do iisubg = 1, nsubg(insubg)
               atcoord(:,iisubg) = sminimgvec(svaux(:),                &
                         rcoord(:,innode+atms(isubg(insubg)+iisubg)),  &
                                                                    box)
               atcoord(:,iisubg) = svaux(:) + atcoord(:,iisubg)
!
!~ write(*,*) iinode,innode,innode+atms(isubg(insubg)+iisubg)
             end do
!
             cofm(:) = scenvec(3,nsubg(insubg),                        &
                               atcoord(:,:nsubg(insubg)))
!
             do i = 1, 3
               if ( cofm(i) .ge. box(i) ) then
                 cofm(i) = cofm(i) - box(i)			
               else if ( cofm(i) .lt. -1.0E-6 ) then
                 cofm(i) = cofm(i) + box(i)
               end if 
             end do 
!
             fcoord(:,j) = cofm(:)            
!
           else
!
             fcoord(:,j) = rcoord(:,innode+atms(isubg(insubg)+1))
!
             do i = 1, 3
               if ( fcoord(i,j) .ge. box(i) ) then
                 fcoord(i,j) = fcoord(i,j) - box(i)			
               else if ( fcoord(i,j) .lt. -1.0E-6 ) then
                 fcoord(i,j) = fcoord(i,j) + box(i)
               end if 
             end do 
!
!~ write(*,*) iinode,innode,innode+atms(isubg(insubg)+1)
           end if    
         end do
!~ write(*,*)
       end do
!
       return
       end subroutine setcoord
!
!======================================================================!
!
       subroutine print_coord(xtcf,sys,outp,msize,nagg,nnode,mol,agg)
!
       use xdr,       only: xtcfile
       use geometry,  only: sminimgvec
       use datatypes
!
       implicit none
!
       include 'inout.h'
!
! Input/output variables
!
       type(xtcfile),intent(inout)           ::  xtcf    !  xtc file informacion   
       type(groinp),intent(in)               ::  sys     !  System information
       character(len=lenout),intent(in)      ::  outp    !  Output file name
       integer,dimension(msize),intent(in)  ::  nagg    !  Number of aggregates of each size
       integer,dimension(nnode),intent(in)   ::  mol    !  Molecule identifier
       integer,dimension(nnode),intent(in)   ::  agg    !  Aggregates size
       integer,intent(in)                    ::  msize  !  Maximum aggregate size
       integer,intent(in)                    ::  nnode   !  Number of residues
!
! Local variables
!
       type(xtcfile)                         ::  xtco    !  xtc file informacion
       character(len=lenout)                 ::  straux  !  Auxiliary string
       character(len=lenout)                 ::  aux     !  Auxiliary string
       real(kind=4),dimension(3)             ::  svaux   !  Auxiliary single precision vector
       real(kind=8),dimension(3)             ::  cofm    !  Center of mass vector
       real(kind=8)                          ::  mass    !  Total mass of the aggregate
       integer                               ::  nsize   !  Size of the previous printed aggregate
       integer                               ::  i,j     !  Indexes
       integer                               ::  m,n     !  Indexes
       integer                               ::  p,q,r,s !  Indexes
!
! Priting coordinates of the aggregates
! -------------------------------------
!
       write(aux,*) xtcf%STEP
       aux = adjustl(aux)
!
       straux = trim(outp)//'_'//trim(aux)//'.xtc'
       aux    = trim(outp)//'_'//trim(aux)//'.xyz'
! Printing global coordinates in xtc format
       call xtco%init(straux,'w')
       call xtco%write(xtcf%natoms,xtcf%step,xtcf%time,xtcf%box,       &
                                                     xtcf%pos,xtcf%prec)
       call xtco%close
! Printing global coordinates in xyz format
       open(unit=uniinp,file=trim(aux),action='write')
!
       write(uniinp,*) xtcf%NATOMS
       write(uniinp,*) sys%title
!
       do i = 1, xtcf%NATOMS
         write(uniinp,*) sys%atname(modulo(i-1,sys%nat)+1),            &
                                                        xtcf%pos(:,i)*10
       end do
!
       close(uniinp)
! Printing coordinates of the complexes by size
       n     = nagg(1)
       nsize = 0
!
       do i = 2, msize
         do j = 1, nagg(i)
           n = n + agg(n)
           if ( nsize .ne. agg(n) ) then
             write(straux,*) agg(n)
             straux = adjustl(straux)
! 
             write(aux,*) xtcf%STEP
             aux = adjustl(aux)
!
             straux = trim(outp)//'_'//trim(aux)//'_'//trim(straux)//  &
                                                                  '.xyz'
             open(unit=uniinp,file=trim(straux),action='write')
!
             nsize = agg(n)
           end if
!
           cofm(:)  = 0.0d0
           mass     = 0.0d0
           svaux(:) = xtcf%pos(:,(mol(n)-1)*sys%nat+1) 
!
           m = n - 1
           do p = 1, agg(n)
             mass = mass + sys%totm
             m = m + 1
             r = (mol(m)-1)*sys%nat
             do q = 1, sys%nat
               r = r + 1 
               xtcf%pos(:,r) = sminimgvec(svaux(:),xtcf%pos(:,r),      &
                                          (/ xtcf%box(1,1),            &
                                             xtcf%box(2,2),            &
                                             xtcf%box(3,3) /) )
               xtcf%pos(:,r) = svaux(:) + xtcf%pos(:,r)
               do s = 1, 3
                 cofm(s) = cofm(s) + sys%mass(q)*xtcf%pos(s,r)
               end do
             end do
           end do
           cofm(:) = cofm(:)/mass
!
           write(uniinp,*) sys%nat*agg(n)
           write(uniinp,'(20(X,I5))') mol(n:n+agg(n)-1)
!
           m = n - 1
           do p = 1, agg(n)
             m = m + 1
             r = (mol(m)-1)*sys%nat
             do q = 1, sys%nat
               r = r + 1
               write(uniinp,'(A5,3(1X,F12.8))') sys%atname(q),         &
                                            (xtcf%pos(:,r) - cofm(:))*10
             end do
           end do
!
           if ( nsize .ne. agg(n+1) ) close(uniinp)
!
         end do
       end do
!
       close(uniinp)
!
! r  -> i
! q  -> iat
! p  -> irenum
! n  -> iitag
! m  -> iimol
! i  -> itype
! j  -> inagg
!
       return
       end subroutine print_coord
!
!======================================================================!
!
! BLOCKDIAG - BLOCK DIAGonalization
!
! This subroutine block diagonalizes the adjacency matrix 
!  ADJ(NNODE,NNODE) of an undirected unweighted graph of NNODE vertices.
!  The subroutine FINDCOMPUNDIR is employed to find the connected compo-
!  nents in the graph using BFS and the QUICKSHORT algorithm is used to
!  find the basis of nodes that block diagonalizes the adjacency matrix
!  by sorting the blocks by size
!
!
       subroutine blockdiag(nnode,adj,natms,posi,matms,inposi,nat,     &
                            thr,atms,mgrps,grps,ngrps,igrps,msubg,     &
                            subg,nsubg,isubg,mol,tag,agg,nsize,nagg,   &
                            iagg,nmol,imol,magg,box,neidis,debug)
!
       use timings,    only: tadj,tbfs,tsort,count_rate
       use graphtools, only: buildadjmolbub,                           &
                             findcompundir
       use sorting,    only: ivvqsort,                                 &
                             ivqsort,                                  &
                             iqsort
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(out)   ::  adj     !  Adjacency matrix
       real(kind=4),dimension(3,natms),intent(out)  ::  posi    !  
       real(kind=4),dimension(3,matms),intent(in)   ::  inposi  !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=8),dimension(nat,nat),intent(in)   ::  thr     !
       real(kind=4),dimension(3),intent(in)         ::  box     !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=8),intent(in)                      ::  neidis  !
       integer,dimension(nnode),intent(out)         ::  mol     !  Molecules identifier
       integer,dimension(nnode),intent(out)         ::  tag     !  Aggregates identifier
       integer,dimension(nnode),intent(out)         ::  agg     !  Aggregates size
       integer,dimension(nat),intent(in)            ::  ngrps   !  
       integer,dimension(nat),intent(in)            ::  igrps   ! 
       integer,dimension(nat),intent(in)            ::  grps    ! 
       integer,dimension(nat),intent(in)            ::  nsubg   !  
       integer,dimension(nat),intent(in)            ::  isubg   ! 
       integer,dimension(nat),intent(in)            ::  subg    ! 
       integer,dimension(nat),intent(in)            ::  atms    ! 
       integer,dimension(nnode),intent(out)         ::  nagg    !  Number of aggregates of each size
       integer,dimension(nnode),intent(out)         ::  iagg    !
       integer,dimension(nnode),intent(out)         ::  nmol    !
       integer,dimension(nnode),intent(out)         ::  imol    !
       integer,intent(in)                           ::  nnode   !  Number of molecules
       integer,intent(in)                           ::  matms   !
       integer,intent(in)                           ::  natms   ! 
       integer,intent(in)                           ::  nat     ! 
       integer,intent(out)                          ::  nsize   !  Maximum aggregate size
       integer,intent(in)                           ::  msubg   !
       integer,intent(in)                           ::  mgrps   !  Number of subgroups
       integer,intent(out)                          ::  magg    !  Number of aggregates
       logical,intent(in)                           ::  debug   !  Debug mode
!
! Local variables
! 
       integer                                      ::  i,j,k   !  Indexes
       integer                                      ::  ni,nj   !
!
! Declaration of time control variables
!
       integer                                      ::  t1adj   !  Initial building time
       integer                                      ::  t2adj   !  Final building time
       integer                                      ::  t1bfs   !  Initial BFS time
       integer                                      ::  t2bfs   !  Final BFS time
       integer                                      ::  t1sort  !  Initial sorting time
       integer                                      ::  t2sort  !  Final sorting time
!
! Building the adjacency matrix for the current snapshot
!
       call system_clock(t1adj) 
!
!~            oldposi(:,:) = posi(:,:)
!
       call setcoord(nnode,msubg,nsubg,isubg,atms,natms,posi,      &
                     matms,inposi,nat,box)
!
!~!        call getdisp(natms,drnei,posi,oldposi,box)
!
!~!            call uplinklist(natms,posi,drnei,list,nhead,head,ncell,     &
!~!                            maxdis,box,xtcf%STEP) 
!
!~            call setlinklist(natms,posi,drnei,list,nhead,head,ncell,    &
!~                             box,xtcf%STEP)   
!
!~            call buildadjmollink(nnode,adj,neidis,msubg,mgrps,          &
!~                                 thr2(:mgrps,:mgrps),ngrps(:mgrps),     &
!~                                 igrps(:mgrps),grps(:msubg),natms,      &
!~                                 posi,list,nhead,head,cell,ncell,box)
!
       call buildadjmolbub(nnode,adj,neidis,msubg,mgrps,nat,thr,       &
                           ngrps,igrps,grps,natms,posi,box)
!
!~            call buildadjmol(mgrps,ngrps,msubg,nsubg,nnode,atms,adj,    &
!~                             thr2,neidis,matms,posi,nat,                &
!~                             sys%mass,box,debug)
!
       call system_clock(t2adj) 
!
       tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)    
!
! Block-diagonalizing the adjacency matrix 
!
       call system_clock(t1bfs)     
!
       call findcompundir(nnode,adj,mol,tag,agg,nsize,nagg,magg,debug)
!
       call system_clock(t2bfs)     
!
       tbfs = tbfs + dble(t2bfs-t1bfs)/dble(count_rate)
! 
! Sorting molecules and aggregate identifiers based on the size of 
!  the aggregates
!
       call system_clock(t1sort)     
!
       iagg(1) = 0
       imol(1) = 0
       do i = 1, nsize-1
         iagg(i+1) = iagg(i) + nagg(i)
         nmol(i)   = i*nagg(i) 
         imol(i+1) = imol(i) + nmol(i)
       end do
       nmol(nsize) = nnode - imol(nsize)
!
!~ write(*,*) 'nagg',nagg
!~ write(*,*) 'iagg',iagg
!~ write(*,*) 'nmol',nmol
!~ write(*,*) 'imol',imol
!~ write(*,*)
!
       call ivvqsort(nnode,agg,tag,mol,1,nnode)
!
! Sorting molecules based on their aggregate identifier
!
       do i = 2, nsize-1
         if ( nagg(i) .gt. 1 ) then
!~ write(*,*) 'sorting from',imol(i)+1,'to',imol(i+1)
           call ivqsort(nnode,tag,mol,imol(i)+1,imol(i+1))
         end if
       end do
!
       if ( nagg(nsize) .gt. 1 )                                       &
                         call ivqsort(nnode,tag,mol,imol(nsize)+1,nnode)
!
! Starts old......
!
!~        if ( nagg(msize) .gt. 1 ) then
!~          ni = imol(msize) + 1
!~          i  = agg(ni)
!~          nj = ni + i
!~          k  = 1
!~ !!~ write(*,*) 'initial ni',ni,'i',i
!~          do while ( nj .le. nnode )
!~            j  = agg(nj)
!~ !!~ write(*,*) '  new nj',nj,'j',j
!~            if ( j .ne. i ) then
!~ !!~ write(*,*) 'sorting from',ni,'to',nj-1,'k',k
!~              if ( k .gt. 1 ) call ivqsort(nnode,tag,mol,ni,nj-1)
!~              k  = 0
!~              i  = j            
!~              ni = nj
!~            end if
!~            k  = k + 1
!~            nj = nj + j
!~          end do
!~ !!~ write(*,*) 'sorting from',ni,'to',nnode,'k',k      
!~          if ( k .gt. 1 ) call ivqsort(nnode,tag,mol,ni,nnode)
!~        end if
!
! ..... end old
!
!~        j = nagg(1)
!~        do i = 2, msize-1
!~          if ( nagg(i) .eq. 0 ) cycle
!~          if ( j .ge. nnode ) exit
!~          call ivqsort(nnode,tag,mol,j+1,j+i*nagg(i))
!~          j = j + i*nagg(i)
!~        end do 
!
!~        i = j ! j = imol(msize)
!~        do while ( i .lt. nnode )
!~          k = agg(i+1)
!~          if ( (i+k+1) .le. nnode ) then
!~            if ( agg(i+k+1) .ne. k ) then
!~              call ivqsort(nnode,tag,mol,j+1,i+k)
!~              j = i + k
!~            end if
!~          else
!~            call ivqsort(nnode,tag,mol,j+1,nnode)
!~          end if
!~          i = i + k
!~        end do
!
! Sorting molecules based on their canonical order
!
       do i = 2, nsize
!~          k = imol(i) + 1  ! FLAG: old
         k = imol(i)
         do j = 1, nagg(i)
!~ write(*,*) 'sorting from',k,'to',k+agg(k)-1
!~            call iqsort(nnode,mol,k,k+agg(k)-1)  ! FLAG: old
           call iqsort(nnode,mol,k+1,k+i)
           k = k + i
         end do
       end do
!
       call system_clock(t2sort)     
!
       tsort = tsort + dble(t2sort-t1sort)/dble(count_rate)
!
       return
       end subroutine blockdiag
!
!======================================================================!
!
! Seen on: Gaines & Di Tommaso. Pharmaceutics, 2018, 10. 10.3390/pharmaceutics10010012
!
!~        subroutine build_pim(ngrps,grps,nsubg,subg,natmol,igrps,thr,    &
!~                             nnode,mol,agg,msize,nagg,natconf,       &
!~                             coord,box,pim,debug)
!~ !
!~        use geometry
!~ !
!~        implicit none
!~ !
!~        include 'idxadj.h'
!~        include 'idxpim.h'
!~ !
!~ ! Input/output variables
!~ !
!~        real(kind=8),dimension(ngrps,ngrps,msize-1),intent(inout)  ::  pim      !  Pairwise interaction matrix
!~        real(kind=8),dimension(natmol,natmol),intent(in)            ::  thr      !  Distance threshold
!~        real(kind=4),dimension(3,natconf),intent(in)                ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
!~        real(kind=4),dimension(3),intent(in)                        ::  box      !  Simulation box !FLAG: kind=8 to kind=4
!~        integer,dimension(msize),intent(in)                        ::  nagg     !  Number of aggregates of each size
!~        integer,dimension(nnode),intent(in)                         ::  mol     !  Molecules identifier
!~        integer,dimension(nnode),intent(in)                         ::  agg     !  Aggregates size
!~        integer,dimension(natmol),intent(in)                        ::  grps     !  Number of subgroups in each group
!~        integer,dimension(natmol),intent(in)                        ::  subg     !  Number of atoms in each subgroup
!~        integer,dimension(natmol),intent(in)                        ::  igrps    !  
!~        integer,intent(in)                                          ::  ngrps    !  Number of groups
!~        integer,intent(in)                                          ::  nsubg    !  Number of subgroups
!~        integer,intent(in)                                          ::  nnode    !  Number of residues
!~        integer,intent(in)                                          ::  msize   !  Maximum aggregate size
!~        integer,intent(in)                                          ::  natconf  !  Total number of atoms
!~        integer,intent(in)                                          ::  natmol   !  Atoms per residue
!~        logical,intent(in)                                          ::  debug    !  Debug mode
!~ !
!~ ! Local variables
!~ !
!~        character(len=54)                                           ::  fmt1     !  Format string
!~        real(kind=4),dimension(3)                                   ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
!~        real(kind=8)                                                ::  dist     !  Minimum image distance
!~        real(kind=8)                                                ::  mindis   !  Distance threshold between groups
!~ !
!~ ! Building pairwise interaction matrix in the subgroup-based representation
!~ ! -------------------------------------------------------------------------
!~ !
!~        fmt1 = '(8X,X,A,X,I2,2(X,A,X,I4),3(X,A,X,I2),X,A,X,I6)'
!~ !
!~        iiagg = nagg(1)
!~ !
!~        do itype = 1, msize-1
!~          if ( nagg(itype+1) .eq. 0 ) cycle
!~          do inagg = 1, nagg(itype+1)
!~            iiagg = iiagg + agg(iiagg)
!~            iimol = iiagg - 1
!~            do irenum = 1, agg(iiagg)-1
!~              iimol  = iimol + 1
!~              ingrps = (mol(iimol)-1)*natmol
!~ !
!~              jimol  = iimol 
!~              do jrenum = irenum+1, agg(iiagg)
!~                jimol  = jimol + 1
!~                jngrps = (mol(jimol)-1)*natmol
!~ !
!~                iisubg = 0
!~                ii     = 0
!~                do iigrps = 1, ngrps
!~                  do insubg = 1, grps(iigrps)
!~                    iisubg = iisubg + 1 
!~                    ii     = ii + subg(iisubg)
!~                    i      = ingrps + igrps(ii)
!~ !
!~                    jisubg = 0
!~                    jj     = 0
!~                    do jigrps = 1, ngrps
!~                      mindis = thr(iigrps,jigrps)
!~                      if ( mindis .gt. 1.0e-6 ) then 
!~                        do jnsubg = 1, grps(jigrps)
!~                          jisubg = jisubg + 1
!~                          jj     = jj + subg(jisubg)
!~                          j      = jngrps + igrps(jj)
!~ !
!~                          r    = sminimgvec(coord(:,i),coord(:,j),box)
!~                          dist = dot_product(r,r)
!~ !
!~                          if ( dist .le. mindis ) then
!~                            pim(iigrps,jigrps,itype) =                  &  ! FLAG: alternatively, to save pim in triangular form
!~                                         pim(iigrps,jigrps,itype) + 1.0d0  !       1) find max(jigrps,iigrps) and min(jigrps,iigrps)
!~                            pim(jigrps,iigrps,itype) =                  &  !       2) then assign pim(max,min) 
!~                                         pim(jigrps,iigrps,itype) + 1.0d0       
!~                          end if
!~                        end do
!~                      else
!~                        jisubg = jisubg + grps(jigrps)
!~                        jj     = jj + subg(jisubg)
!~                      end if
!~                    end do            
!~                  end do
!~                end do
!~              end do
!~            end do
!~          end do
!~        end do
!~ !
!~        return
!~        end subroutine build_pim
!
!======================================================================!
!
!~        subroutine analyze_agg(ngrps,grps,nsubg,subg,natmol,igrps,thr,  &
!~                               nnode,adj,imol,itag,msize,nagg,natconf, &
!~                               coord,box,debug)
!~ !
!~        use graphtools, only: calcdegundir,                             &
!~                              chktree,                                  &
!~                              chkltree,                                 &             
!~                              chkscycle,                                &
!~                              findshortc,                               &
!~                              findlongt
!~        use geometry,   only: minimgvec
!~ !
!~        implicit none
!~ !
!~        include 'idxadj.h'
!~        include 'idxpim.h'
!~        include 'idxbfs.h'
!~ !
!~ ! Input/output variables
!~ !
!~        logical,dimension(nnode,nnode),intent(in)         ::  adj      !  Adjacency matrix in the molecule representation
!~        real(kind=8),dimension(natmol,natmol),intent(in)  ::  thr      !  Distance threshold
!~        real(kind=4),dimension(3,natconf),intent(in)      ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
!~        real(kind=4),dimension(3),intent(in)              ::  box      !  Simulation box !FLAG: kind=8 to kind=4
!~        integer,dimension(msize),intent(in)              ::  nagg     !  Number of aggregates of each size
!~        integer,dimension(nnode),intent(in)               ::  imol     !  Molecules identifier
!~        integer,dimension(nnode),intent(in)               ::  itag     !  Aggregates size
!~        integer,dimension(natmol),intent(in)              ::  grps     !  Number of subgroups in each group
!~        integer,dimension(natmol),intent(in)              ::  subg     !  Number of atoms in each subgroup
!~        integer,dimension(natmol),intent(in)              ::  igrps    !  
!~        integer,intent(in)                                ::  ngrps    !  Number of groups
!~        integer,intent(in)                                ::  nsubg    !  Number of subgroups
!~        integer,intent(in)                                ::  nnode    !  Number of residues
!~        integer,intent(in)                                ::  msize   !  Maximum aggregate size
!~        integer,intent(in)                                ::  natconf  !  Total number of atoms
!~        integer,intent(in)                                ::  natmol   !  Atoms per residue
!~        logical,intent(in)                                ::  debug    !  Debug mode
!~ !
!~ ! Local variables
!~ !
!~        logical,dimension(:,:),allocatable                ::  auxadj   !  Auxiliary adjacency matrix
!~        integer,dimension(:),allocatable                  ::  degree   !  Degree of each vertex
!~        character(len=54)                                 ::  fmt1     !  Format string
!~        real(kind=4),dimension(3)                         ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
!~        real(kind=8)                                      ::  dist     !  Minimum image distance
!~        real(kind=8)                                      ::  mindis   !  Distance threshold between groups
!~        integer,allocatable,dimension(:)                  ::  icmol    !  Molecules-in-cycle identifier
!~        integer,allocatable,dimension(:)                  ::  icycle   !  Cycles identifier
!~        integer,allocatable,dimension(:)                  ::  ictag    !  Cycles size
!~        integer,allocatable,dimension(:)                  ::  ncmol    !  Number of cycles of each size
!~        integer                                           ::  ncycle   !  Number of simple cycles       logical                                           ::  chk      !  Checking variable
!~        logical                                           ::  chklt    !  Linear tree checking variable
!~        logical                                           ::  chkbt    !  Binary tree checking variable
!~ !
!~ ! Analyzing aggregates of size greater than 2
!~ ! -------------------------------------------
!~ !
!~        iitag = nagg(1) + nagg(2)*2 - 1
!~ !
!~        do itype = 2, msize-1
!~          if ( nagg(itype+1) .eq. 0 ) cycle
!~          do inagg = 1, nagg(itype+1)
!~            iitag = iitag + itag(iitag)
!~            iimol = iitag - 1
!~ !
!~ ! Saving the adjacency matrix of the current aggregate in the molecule-
!~ !  based representation
!~ !
!~            allocate(auxadj(itag(iitag),itag(iitag)),degree(itag(iitag)))
!~            allocate(icmol(itag(iitag)),icycle(itag(iitag)),            &
!~                     ictag(itag(iitag)),ncmol(itag(iitag)) )
!~ !
!~            auxadj(:,:) = .FALSE.
!~ !
!~            write(*,'(2X,A)')          '-------------------------'
!~            write(*,'(2X,A,20(X,I5))') 'Saving block of molecules',(imol(i),i=iimol+1,iimol+itag(iitag))
!~            write(*,'(2X,A)')          '-------------------------'
!~            write(*,*)
!~ !
!~            do irenum = 1, itag(iitag)-1
!~              iimol = iimol + 1
!~ !
!~              jimol = iimol 
!~              do jrenum = irenum+1, itag(iitag)
!~                jimol = jimol + 1
!~ !
!~                auxadj(irenum,jrenum) = adj(imol(iimol),imol(jimol))
!~                auxadj(jrenum,irenum) = auxadj(irenum,jrenum)
!~              end do
!~            end do
!~ !
!~            do i = 1, itag(iitag)
!~              write(*,'(5X,20L2)') (auxadj(i,j),j=1,itag(iitag))
!~            end do
!~            write(*,*)
!~ !
!~ ! Analyzing adjacency matrix of the current aggregate
!~ !
!~ ! Calculating the degree of each vertex
!~            degree  = calcdegundir(itag(iitag),auxadj)
!~ !
!~            write(*,'(6X,A,20(X,I2))') 'Degrees',degree
!~            write(*,*)
!~ ! Checking if the aggregate is a tree
!~            if ( chktree(itag(iitag),degree) ) then 
!~              write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inagg,'of type',itype+1,'is a tree'
!~ ! Checking if the aggregate forms a linear tree
!~              if ( chkltree(itag(iitag),degree) ) then 
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inagg,'of type',itype+1,'is a linear tree'
!~                write(*,'(4X,A,X,I4)') 'The length of the longest chain is',itag(iitag)
!~              else 
!~ ! If the tree is not linear then find the longest molecular chain
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inagg,'of type',itype+1,'is a n-ary tree'
!~ !
!~                write(*,'(4X,A,X,I4)') 'The length of the longest chain is',findlongt(itag(iitag),auxadj)
!~              end if
!~ ! If the graph is not a tree then it is a cyclic graph
!~            else
!~              write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inagg,'of type',itype+1,'is a cycle'
!~ ! Checking if the aggregate forms a single cycle through all nodes of the component
!~              if ( chkscycle(itag(iitag),degree) ) then 
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inagg,'of type',itype+1,'is a simple cycle'
!~                write(*,'(4X,A,X,I4)') 'The length of the shortest simple cycle is',itag(iitag)
!~              else 
!~ ! If the aggregate does not form a single cycle then find the shortest simple cycle
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inagg,'of type',itype+1,'is a n-cycle'
!~ !                 
!~                write(*,'(4X,A,X,I4)') 'The length of the shortest simple cycle is',findshortc(itag(iitag),auxadj,degree)
!~              end if
!~ ! 
!~            end if
!~            write(*,*)
!~ !
!~            deallocate(auxadj,degree)
!~            deallocate(icmol,icycle,ictag,ncmol)
!~ !
!~          end do
!~        end do
!~ !
!~        return
!~        end subroutine analyze_agg
!
!======================================================================!
!
!~        subroutine analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,  &
!~                               natmol,igrps,thr,nnode,adj,imol,itag,    &
!~                               msize,nagg,natconf,coord,box,posi,      &
!~                               atmass,atname,step,outp,debug)
!~ !
!~        use geometry,   only: sminimgvec
!~        use graphtools, only: buildadjbody
!~ !
!~        implicit none
!~ !
!~        include 'inout.h'
!~        include 'idxpim.h'
!~        include 'idxadj.h'
!~ !
!~ ! Input/output variables
!~ !
!~        real(kind=8),dimension(9,9,3),intent(inout)       ::  table    !  Number of aggregates of each type
!~        logical,dimension(nnode,nnode),intent(in)         ::  adj      !  Adjacency matrix in the molecule representation
!~        real(kind=8),dimension(natmol,natmol),intent(in)  ::  thr      !  Distance threshold
!~        real(kind=4),dimension(3,natconf),intent(inout)   ::  posi     !  Auxiliary coordinates
!~        real(kind=4),dimension(3,natconf),intent(in)      ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
!~        real(kind=4),dimension(3),intent(in)              ::  box      !  Simulation box !FLAG: kind=8 to kind=4
!~        real(kind=8),dimension(natmol),intent(in)         ::  atmass   !
!~        character(len=lenout),intent(in)                  ::  outp     !  Output file name
!~        character(len=5),dimension(natmol),intent(in)     ::  atname   !      
!~        integer,dimension(msize),intent(in)              ::  nagg     !  Number of aggregates of each size
!~        integer,dimension(nnode),intent(in)               ::  imol     !  Molecules identifier
!~        integer,dimension(nnode),intent(in)               ::  itag     !  Aggregates size
!~        integer,dimension(natmol),intent(in)              ::  body     !  Number of groups in each body
!~        integer,dimension(natmol),intent(in)              ::  grps     !  Number of subgroups in each group
!~        integer,dimension(natmol),intent(in)              ::  subg     !  Number of atoms in each subgroup
!~        integer,dimension(natmol),intent(in)              ::  igrps    !  Atoms identifier
!~        integer,intent(in)                                ::  nbody    !  Number of bodies
!~        integer,intent(in)                                ::  ngrps    !  Number of groups
!~        integer,intent(in)                                ::  nsubg    !  Number of subgroups
!~        integer,intent(in)                                ::  nnode    !  Number of residues
!~        integer,intent(in)                                ::  msize   !  Maximum aggregate size
!~        integer,intent(in)                                ::  natconf  !  Total number of atoms
!~        integer,intent(in)                                ::  natmol   !  Atoms per residue
!~        integer,intent(in)                                ::  step     !
!~        logical,intent(in)                                ::  debug    !  Debug mode
!~ !
!~ ! Local variables
!~ !
!~        real(kind=8),dimension(3)                         ::  cofm     !  Center of mass vector
!~        real(kind=4),dimension(3)                         ::  svaux    !  Auxiliary single precision vector
!~        real(kind=8)                                      ::  mass     !  Total mass of the aggregate
!~        logical,dimension(:,:),allocatable                ::  adjmol   !  Adjacency matrix in the molecule-based rep.  ! FLAG: rename as adjaux
!~        logical,dimension(:,:),allocatable                ::  adjaux   !  Auxiliary adjacency matrix                   ! FLAG: rename as adjbody
!~        integer,dimension(:),allocatable                  ::  degree   !  Degree of each vertex
!~        integer                                           ::  ijmol
!~        integer                                           ::  idxbody  !  Index for the body-body interactions
!~        integer                                           ::  idxoh    !  Index for the OH-OH interactions
!~        integer                                           ::  idxph    !  Index for the Ph-Ph interactions
!~        integer                                           ::  uniconf  !
!~        character(len=lenout)                             ::  aux      !  Auxiliary string
!~ !
!~ ! Analyzing aggregates of size lower or equal than 4
!~ ! --------------------------------------------------
!~ !
!~        iitag = nagg(1)
!~ !
!~        do itype = 2, 4
!~          if ( nagg(itype) .eq. 0 ) cycle
!~          do inagg = 1, nagg(itype)
!~            iitag = iitag + itag(iitag)
!~            iimol = iitag - 1
!~ !
!~            allocate(adjmol(itag(iitag),itag(iitag)),degree(itag(iitag)))
!~            allocate(adjaux(itag(iitag)*nbody,itag(iitag)*nbody))
!~ !
!~            adjmol(:,:) = .FALSE.
!~ !
!~ ! Saving the adjacency matrix in the molecule-based representation
!~ !
!~            i = iimol
!~            do irenum = 1, itag(iitag)-1
!~              i = i + 1
!~ !
!~              jimol = i 
!~              do jrenum = irenum+1, itag(iitag)
!~                jimol = jimol + 1
!~ !
!~                adjmol(irenum,jrenum) = adj(imol(i),imol(jimol))
!~                adjmol(jrenum,irenum) = adjmol(irenum,jrenum)
!~              end do
!~            end do
!~ !
!~ !~!            write(*,'(2X,A)')          '-------------------------'
!~ !~!            write(*,'(2X,A,20(X,I5))') 'Saving block of molecules',(imol(i),i=iimol+1,iimol+itag(iitag))
!~ !~!            write(*,'(2X,A)')          '-------------------------'
!~ !~!            write(*,*)
!~ !~!            write(*,'(2X,2(A,X,I4,X))') 'Analyzing aggregate',inagg,'of type',itype
!~ !~!            write(*,*)
!~ !~!
!~ !~!            write(*,'(5X,A)') 'Adjacency matrix in the molecule-based representation'
!~ !~!            write(*,'(5X,A)') '-----------------------------------------------------'
!~ !~!            do i = 1, itag(iitag)
!~ !~!              write(*,'(5X,20L2)') (adjmol(i,j),j=1,itag(iitag))
!~ !~!            end do
!~ !~!            write(*,*)
!~ !
!~ ! Building the adjacency matrix in the N-body simplified representation
!~ !
!~            call buildadjbody(natmol,nbody,body,ngrps,grps,nsubg,       &
!~                              subg,igrps,adjaux,thr,itag(iitag),        &
!~                              imol(iimol+1:iimol+itag(iitag)),          &
!~                              adjmol,natconf,coord,box,debug)
!~ !
!~ !~!            write(*,'(5X,A)') 'Adjacency matrix in the N-body simplified representation'
!~ !~!            write(*,'(5X,A)') '--------------------------------------------------------'
!~ !~!            do i = 1, itag(iitag)*nbody
!~ !~!              write(*,'(5X,20L2)') (adjaux(i,j),j=1,itag(iitag)*nbody)
!~ !~!            end do
!~ !~!            write(*,*)
!~ !
!~ ! Analyzing the adjacency matrix for each body-body interaction
!~ !
!~            do iibody = 1, nbody
!~ !
!~              adjmol(:,:) = .FALSE.
!~ !
!~              inbody = 0
!~              do irenum = 1, itag(iitag)-1
!~                jnbody = inbody + nbody
!~                do jrenum = irenum+1, itag(iitag)
!~                  adjmol(jrenum,irenum) = adjaux(jnbody+iibody,inbody+iibody)
!~                  adjmol(irenum,jrenum) = adjmol(jrenum,irenum)
!~ !~!                  write(*,*) irenum,jrenum,inbody+iibody,jnbody+iibody,adjaux(inbody+iibody,jnbody+iibody)
!~ !
!~                  jnbody = jnbody + nbody
!~                end do
!~                inbody = inbody + nbody
!~              end do
!~ !
!~ !~!              write(*,'(7X,A,X,I2)') 'Adjacency matrix for bodies of type',iibody
!~ !~!              write(*,'(7X,A)')      '---------------------------------------'
!~ !~!              do i = 1, itag(iitag)
!~ !~!                write(*,'(7X,20L2)') (adjmol(i,j),j=1,itag(iitag))
!~ !~!              end do
!~ !~!              write(*,*)
!~ !
!~ ! Removing rows and columns that contain only zeros (disconnected bodies)
!~ !
!~              i  = 1
!~              j  = itag(iitag)
!~ !
!~              do irenum = 1, itag(iitag)
!~                degree(i) = 0
!~ ! Computing the degree of vertex i
!~                do jrenum = 1, j
!~                  if ( adjmol(jrenum,i) ) degree(i) = degree(i) + 1
!~                end do
!~ ! Pushing disconnected bodies to the end of the arrays
!~                if ( degree(i) .eq. 0 ) then
!~ ! Permuting columns of the adjacency matrix
!~                  adjmol(:,i) = adjmol(:,j)
!~                  adjmol(:,j) = .FALSE.
!~ ! Permuting rows of the adjacency matrix
!~                  adjmol(i,:) = adjmol(j,:)
!~                  adjmol(j,:) = .FALSE.
!~ ! Permuting degrees of the current node i and the last node j
!~                  degree(j) = 0
!~ !
!~                  j = j - 1
!~                else
!~                  i = i + 1
!~                end if
!~              end do
!~ !
!~ !~!              write(*,'(7X,A,20(X,I2))') 'Size   ',j
!~ !~!              write(*,'(7X,A,20(X,I2))') 'Degrees',degree
!~ !~!              write(*,*)
!~ !~!              write(*,'(9X,A,X,I2)') 'New adjacency matrix for bodies of type',iibody
!~ !~!              write(*,'(9X,A)')      '-------------------------------------------'
!~ !~!              do irenum = 1, itag(iitag)
!~ !~!                write(*,'(9X,20L2)') (adjmol(irenum,jrenum),            &
!~ !~!                                                    jrenum=1,itag(iitag))
!~ !~!              end do
!~ !~!              write(*,*)
!~ !
!~ ! Classifying aggregate according to its interactions
!~ !
!~              if ( iibody .eq. 1 ) then
!~                idxoh = idxbody(j,adjmol(:j,:j),degree(:j))
!~              else
!~                idxph = idxbody(j,adjmol(:j,:j),degree(:j))
!~              end if
!~ !
!~            end do  
!~ !
!~            table(idxoh,idxph,itype-1) = table(idxoh,idxph,itype-1)     &
!~                                                                  + 1.0d0
!~ !
!~            uniconf = itype*100 + idxoh*10 + idxph
!~ !
!~            write(aux,*) uniconf
!~            aux = adjustl(aux)
!~ !
!~            aux = trim(outp)//'_'//trim(aux)//'.xyz'
!~ !
!~            cofm(:)  = 0.0d0
!~            mass     = 0.0d0
!~            svaux(:) = posi(:,(imol(iitag)-1)*natmol+1)
!~ !
!~            ijmol = iimol
!~            do irenum = 1, itag(iitag)
!~              ijmol = ijmol + 1
!~              i     = (imol(ijmol)-1)*natmol
!~              do iat = 1, natmol
!~                i = i + 1
!~ ! 
!~                mass = mass + atmass(iat)
!~ !
!~                posi(:,i) = sminimgvec(svaux(:),posi(:,i),box)
!~                posi(:,i) = svaux(:) + posi(:,i)
!~ !
!~                do j = 1, 3
!~                  cofm(j) = cofm(j) + atmass(iat)*posi(j,i)
!~                end do
!~              end do
!~            end do
!~ !
!~            cofm(:) = cofm(:)/mass           
!~ !
!~            open(unit=uniconf,file=trim(aux),position='append',         &
!~                 action='write')
!~ !
!~            write(uniconf,*) natmol*itag(iitag)
!~            write(uniconf,'(1X,A,1X,I12,20(X,I5))')'STEP=',step,        &
!~                                          imol(iitag:iitag+itag(iitag)-1)
!~ !
!~            ijmol = iimol
!~            do irenum = 1, itag(iitag)
!~              ijmol = ijmol + 1
!~              i     = (imol(ijmol)-1)*natmol
!~              do iat = 1, natmol
!~                i = i + 1
!~ ! 
!~                write(uniconf,'(A5,3(1X,F12.8))') atname(iat),          &
!~                                                 (posi(:,i) - cofm(:))*10
!~ !
!~              end do
!~            end do
!~ !
!~            close(uniconf)
!~ ! 
!~            deallocate(adjmol,degree)
!~            deallocate(adjaux)
!~ !
!~          end do
!~        end do
!~ !
!~        return
!~        end subroutine analyze_agg
!
!======================================================================!
!
!~        integer function idxbody(nnode,adj,degree)
!~ !
!~        use graphtools, only: chktree,                                  &
!~                              chkltree,                                 &             
!~                              chkscycle
!~ !
!~        implicit none
!~ !
!~ ! Input/output variables
!~ !
!~        logical,dimension(nnode,nnode),intent(in)  ::  adj      !  Adjacency matrix in the N-body simplified rep.
!~        integer,dimension(nnode),intent(in)        ::  degree   !  Degree of each vertex
!~        integer,intent(in)                         ::  nnode    !  Number of bodies
!~ !
!~ ! Local variables
!~ !
!~        real(kind=8)                               ::  dpaux1   !
!~        real(kind=8)                               ::  dpaux2   !
!~        integer,dimension(nnode)                   ::  degaux   !  Degree of each vertex of the largest block
!~        integer,dimension(nnode)                   ::  ibody    !  Bodies identifier (mol)
!~        integer,dimension(nnode)                   ::  iblck    !  Blocks identifier (agg)
!~        integer,dimension(nnode)                   ::  isize    !  Blocks size       (tag)
!~        integer,dimension(nnode)                   ::  nsize    !  Number of blocks of each size (nagg)
!~        integer                                    ::  nblck    !  Number of blocks (magg)
!~        integer                                    ::  msize    !  Size of the largest block
!~        integer                                    ::  i        !  Index
!~        logical                                    ::  false    !
!~ !
!~ ! Classifying aggregate according to its interactions
!~ ! ---------------------------------------------------
!~ !
!~        false = .FALSE.
!~ !
!~        if ( nnode .eq. 0 ) then
!~          idxbody = 9
!~        else if ( nnode .gt. 2 ) then
!~ ! Block-diagonalizing the adjacency matrix for bodies of type iibody
!~          call blockdiag(nnode,adj,ibody,iblck,isize,nnode,nsize,       &
!~                         nblck,dpaux1,dpaux2,1,false)
!~ !
!~          msize = isize(nnode)
!~ !
!~ !~!          write(*,'(9X,A,X,20(X,I3))') 'ibody : ',ibody(:)
!~ !~!          write(*,'(9X,A,X,20(X,I3))') 'iblck : ',iblck(:)
!~ !~!          write(*,'(9X,A,X,20(X,I3))') 'isize : ',isize(:)
!~ !~!          write(*,'(9X,A,X,20(X,I3))') 'nsize : ',nsize(:)
!~ !~!          write(*,'(9X,A,X,20(X,I3))') 'msize : ',msize
!~ !~!          write(*,'(9X,A,X,20(X,I3))') 'nblck : ',nblck
!~ !~!          write(*,*)
!~ ! Marking the aggregate according to its topology
!~          select case ( msize )
!~            case (2)
!~ !
!~              idxbody = 2
!~ !
!~            case (3)
!~ ! Saving the degrees of the largest block vertices
!~              do i = 1, msize
!~                degaux(i) = degree(ibody(nnode-msize+i))
!~              end do
!~ !
!~ !~!              write(*,'(9X,A,20(X,I2))') 'Degrees',degaux(:msize)
!~ !~!              write(*,*)
!~ ! Checking if the aggregate is linear or cyclic
!~              if ( chktree(msize,degaux(:msize)) ) then 
!~ !~!                write(*,'(9X,A)') 'The bodies form a linear tree'
!~ !~!                write(*,*)
!~                idxbody = 3
!~              else
!~ !~!                write(*,'(9X,A)') 'The bodies form a simple cycle'
!~ !~!                write(*,*)
!~                idxbody = 4
!~              end if
!~ !
!~            case (4)
!~ ! Saving the degrees of the largest block vertices
!~              do i = 1, msize
!~                degaux(i) = degree(ibody(nnode-msize+i))
!~              end do
!~ !
!~ !~!              write(*,'(9X,A,20(X,I2))') 'Degrees',degaux(:msize)
!~ !~!              write(*,*)
!~ ! Checking if the aggregate is linear or cyclic
!~              if ( chktree(msize,degaux(:msize)) ) then 
!~                if ( chkltree(msize,degaux(:msize)) ) then 
!~ !~!                  write(*,'(9X,A)') 'The bodies form a linear tree'
!~ !~!                  write(*,*)
!~                  idxbody = 6
!~                else
!~ !~!                  write(*,'(9X,A)') 'The bodies form a tree'
!~ !~!                  write(*,*)
!~                  idxbody = 5
!~                end if
!~              else
!~                if ( chkscycle(msize,degaux(:msize)) ) then 
!~ !~!                  write(*,'(9X,A)') 'The bodies form a simple cycle'
!~ !~!                  write(*,*)
!~                  idxbody = 8
!~                else
!~ !~!                  write(*,'(9X,A)') 'The bodies form a cycle'
!~ !~!                  write(*,*)
!~                  idxbody = 7
!~                end if
!~              end if
!~ !
!~          end select
!~ !
!~        else
!~ !
!~          idxbody = 1
!~ !
!~        end if
!~ !
!~        return
!~        end function idxbody
!
!======================================================================!
