!======================================================================!
!
       program aggregate
!
       use xdr,           only:  xtcfile
!
       use systeminf
       use properties
       use filenames
       use thresholds
!
       use parameters
       use timings
       use lengths,       only:  leninp,lentag,lenschm
       use omp_var
!
       use aggtools,      only:  driver
       use graphtools,    only:  bonds2adj
       use printings
       use utils,         only:  rndmseed
       use input_section, only:  read_inp
       use gmx_files,     only:  read_gro,top_parser,count_moltype
!
       use omp_lib
!
       implicit none
!
       include 'inout.h'
!
! Aggregates information in the molecule-based representation
!
       integer,dimension(:),allocatable                ::  nnode    !  Total number of molecules of each size
       integer,dimension(:),allocatable                ::  inode    !  Initial molecule of each size
       integer,dimension(:),allocatable                ::  nat      !  Number of atoms in each molecule type
       integer,dimension(:),allocatable                ::  iat      !  Initial atoms of each molecule type
       integer,dimension(:),allocatable                ::  natms    !
       integer,dimension(:),allocatable                ::  iatms    !
       integer,dimension(:),allocatable                ::  ngrps    !  Number of active atoms in each representation type
       integer,dimension(:),allocatable                ::  igrps    !  Initial active atom of each representation type
       integer                                         ::  ntype    !  Number of molecules types
       integer                                         ::  mnode    !  Total number of molecules
       integer                                         ::  maxat    !  Total number of atoms in the system
       integer                                         ::  mat      !  Total number of atoms in the molecules
       integer                                         ::  matms    !  Total number of subgroups in the molecules
       integer                                         ::  mgrps    !  
       integer                                         ::  msize    !  Maximum aggregate size
       integer                                         ::  nmax     !  
!
! Trajectory control variables
!     
       integer                                         ::  nprint   !  Populations printing interval
       integer                                         ::  minstep  !  First step for analysis
       integer                                         ::  maxstep  !  Last step for analysis
       integer                                         ::  nsteps   !  Number of snapshots analyzed
!
! Lifetimes calculation variables
!
       real(kind=8),dimension(:),allocatable           ::  avlife   !
       integer,dimension(:),allocatable                ::  nlife    !
!
! AnalysisPhenolMD variables
!     
       integer                                         ::  nsolv    !
!
! Declaration of time control variables
!
       real(kind=8)                                    ::  tin      !  Initial CPU time
       real(kind=8)                                    ::  tfin     !  Final CPU time
       integer                                         ::  t1,t2    !  Wall times
       integer                                         ::  t1read   !  Initial reading time
       integer                                         ::  t2read   !  Final reading time
!
! Program control flags
!
       character(len=lenschm)                          ::  schm     !  Calculation scheme flag
       character(len=lenschm)                          ::  scrn     !  Calculation screening flag
       logical                                         ::  doscrn   !  Screening calculation flag
       logical                                         ::  dolife   !  Lifetimes calculation flag
       logical                                         ::  dopim    !  PIM calculation flag
       logical                                         ::  seed     !  Random seed flag
       logical                                         ::  debug    !  Debug mode
!
! OpenMP variables
!
       integer                                         ::  myid     !
       integer                                         ::  nproc    !
!
! Auxiliary variables
!
       integer                                         ::  lin      !  
       integer                                         ::  lfin     ! 
       integer                                         ::  io       !  Status
       integer                                         ::  i,j,k    !  Indexes
       integer                                         ::  itype    !  Indexes
!
! Printing header
!
       write(*,*)
       call print_start()
!
! Initializing line lengths
!
       lin  = 45
       lfin = 90
!
! Initializing timings
!
       call cpu_time(tin)
!
       call system_clock(count_max=count_max,count_rate=count_rate)
       call system_clock(t1)  
       call system_clock(t1read)        
!
       tcpu  = 0.0d0
       twall = 0.0d0
       tread = 0.0d0
       tadj  = 0.0d0
       tbfs  = 0.0d0
       tsort = 0.0d0
       tscrn = 0.0d0
       tpim  = 0.0d0
       tconf = 0.0d0
!
       tcpuadj  = 0.0d0
       tcpusort = 0.0d0
       tcpubfs  = 0.0d0
       tcpuscrn = 0.0d0
       tcpulife = 0.0d0
!
! Reading command line options
!
       call command_line(traj,conf,inp,outp,nprint,minstep,maxstep,    & 
                         msize,neidis,schm,scrn,dopim,dolife,doscrn,   &
                         seed,chunkadj,chunkscrn,chunklife,weight,     &
                         nsolv,debug)
!
!
! 
       if ( trim(outp) .eq. '[none]' ) outp = inp(:len_trim(inp)-4)//'.dat'
!
! Printing the number of threads available
!
!$omp parallel
!
      myid  = OMP_GET_THREAD_NUM()
      nproc = OMP_GET_NUM_THREADS()
!
      if ( myid .eq. 0 ) then
        write(*,'(2X,A,1X,I2,1X,A)') 'Running over',nproc,             &
                                                        'OpenMP threads'
        write(*,*)
      end if
!
!$omp end parallel
!
! Initializing random number generator
!
       call rndmseed(seed)
!
! Printing input data file
!
       call print_inp(uniinp,inp,lfin)
!
! Printing summary of the input information
!
       call print_title(6,1,'Input information','-')
       write(*,*)
       call line_str(6,2,'General input file name',lin,':',            &
                     trim(inp),lfin)
       call line_str(6,2,'Trajectory file name',lin,':',trim(traj),lfin)
       call line_str(6,2,'Configuration file name',lin,':',            &
                     trim(conf),lfin)
       call line_str(6,2,'Output file name',lin,':',trim(outp),lfin)
       write(*,*)
!
       call line_str(6,2,'Restraints algorithm',lin,':',trim(schm),lfin)
       if ( doscrn ) call line_str(6,2,'Screening algorithm',lin,':',  &
                                   trim(scrn),lfin)
       call line_log(6,2,'Screening calculation',lin,':',doscrn,lfin)
       call line_log(6,2,'Lifetimes calculation',lin,':',dolife,lfin)
       write(*,*)
!
       call line_dp(6,2,'Screening distance',lin,':','F5.2',           &
                    neidis,lfin)
       call line_int(6,2,'Number of additional molecules',lin,':',     &
                     'I6',nsolv,lfin)
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
! Opening trajectory file
!
       if ( traj(len_trim(traj)-3:) .eq. '.xtc' ) then 
         call xtcf%init(trim(traj))
         call xtcf%read
       else
         write(*,*) 'Incorrect extension!'   !  TODO: Handle .xtc .trr .gro .xyz ...
         write(*,*)
         call print_end()
       end if
!
! Processing monomer configuration input file
!
       if ( conf(len_trim(conf)-3:) .eq. '.gro' ) then  
!
         allocate(sys(1),rep(1))
         allocate(nnode(1))        ! TODO: only reads 1 gro file, only one molecule type
!                                  !        input # gro files and # number of molecules
         call read_gro(uniinp,conf,sys(1))  
!
         nnode(1) = xtcf%natoms/sys(1)%nat
!
       else if ( conf(len_trim(conf)-3:) .eq. '.top' ) then
!
         call count_moltype(uniinp,conf,ntype)
!
         allocate(sys(ntype),rep(ntype))
!
         allocate(nnode(ntype),inode(ntype))
         allocate(nat(ntype),iat(ntype))
         allocate(natms(ntype),iatms(ntype))
         allocate(ngrps(ntype),igrps(ntype))
!
         call top_parser(uniinp,conf,ntype,nnode)  
!
       else
         write(*,*) 'Incorrect extension!' !  TODO: Input .gro .top .xyz ...
         write(*,*)
         call print_end()
       end if
!
! Building adjacency matrix of the monomer
!
       if ( conf(len_trim(conf)-3:) .eq. '.gro' ) then  ! TODO: guess adj from atomic radii, etc.
!
         stop 'Adjacency matrix from gro file is not yet implemented!'
!
!~          do i = 1, ntype
!~            allocate(sys(i)%adj(sys(i)%nat,sys(i)%nat))
!~            call radii2adj(sys(i)%nbond,sys(i)%ibond,sys(i)%nat,        &
!~                           sys(i)%adj)
!~          end do
!
       else if ( conf(len_trim(conf)-3:) .eq. '.top' ) then
!
         do i = 1, ntype
!
           allocate(sys(i)%adj(sys(i)%nat,sys(i)%nat))
!
           call bonds2adj(sys(i)%nbond,sys(i)%ibond,sys(i)%nat,        &
                          sys(i)%adj)
!
         end do
!
       end if
!
! Opening output files
!
       if ( outp(len_trim(outp)-3:) .eq. '.dat' ) then 
         outp = outp(:len_trim(outp)-4)
       end if
!
       open(unit=uniout+1,file=trim(outp)//'_pop.dat',action='write')
       open(unit=uniout+2,file=trim(outp)//'_frac.dat',action='write')
       open(unit=uniout+3,file=trim(outp)//'_conc.dat',action='write')
       open(unit=uniout+4,file=trim(outp)//'_prob.dat',action='write')
       open(unit=uniout+5,file=trim(outp)//'_num.dat',action='write')
!
! Setting up system size information
!
       mnode = 0
       mat   = 0
       maxat = 0
!
       do itype = 1, ntype
!
         mnode = mnode + nnode(itype)
         mat   = mat   + sys(itype)%nat
!
         nat(itype) = sys(itype)%nat*nnode(itype)
         maxat      = maxat + sys(itype)%nat*nnode(itype)
!
         rep(itype)%nat = sys(itype)%nat
         rep(itype)%iat = 0 
!
       end do
!
       inode(:) = 0
       iat(:)   = 0
!
       do itype = 2, ntype
!
         inode(itype)   = inode(itype-1)   + nnode(itype-1)
         iat(itype)     = iat(itype-1)     + nat(itype-1)
         rep(itype)%iat = rep(itype-1)%iat + rep(itype-1)%nat
!
       end do
!
! Allocating variables depending on system size
!
       allocate(thr(mat,mat))
       allocate(thr2(mat,mat))
!
       allocate(thrang(mat,mat))
!
       do i = 1, ntype
!
         allocate(rep(i)%nbody(sys(i)%nat))
         allocate(rep(i)%ngrps(sys(i)%nat))
         allocate(rep(i)%nsubg(sys(i)%nat))
!
         allocate(rep(i)%ibody(sys(i)%nat))
         allocate(rep(i)%igrps(sys(i)%nat))
         allocate(rep(i)%isubg(sys(i)%nat))
!
         allocate(rep(i)%body(sys(i)%nat))
         allocate(rep(i)%grps(sys(i)%nat))
         allocate(rep(i)%subg(sys(i)%nat))
!
         allocate(rep(i)%grptag(sys(i)%nat))
         allocate(rep(i)%atms(sys(i)%nat))
! 
         allocate(rep(i)%neiang(sys(i)%nat))
!
       end do
!
! Setting up neiang array (first neighbour index)   TODO: find nei of H atoms with degree 1
!
!~        call setneiang
!
! Processing general input file
!
       call read_inp(inp,ntype,rep,mat,thr,thr2,thrang)
!
! Initializing topological information
!         
       ngrps(:) = 0
!
       matms = 0
       mgrps = 0
!
       do itype = 1, ntype
!
         natms(itype) = rep(itype)%msubg*nnode(itype)
         ngrps(itype) = rep(itype)%mgrps
!
         matms = matms + natms(itype)
         mgrps = mgrps + ngrps(itype)
!
         rep(itype)%ibody(:) = 0
         rep(itype)%igrps(:) = 0
         rep(itype)%isubg(:) = 0
!
         rep(itype)%body(:) = 0
         rep(itype)%grps(:) = 0
         rep(itype)%subg(:) = 0
!
! Setting up igrps and ngrps arrays
!
         k = 1
         do i = 1, rep(itype)%mgrps
           if ( i .lt. rep(itype)%mgrps ) then
             rep(itype)%igrps(i+1) = rep(itype)%igrps(i) + rep(itype)%ngrps(i)
           end if
!
           do j = 1, rep(itype)%ngrps(i)
             rep(itype)%grps(k) = i
             k = k + 1
           end do
         end do
!
! Setting up isubg and nsubg arrays
!
         do i = 1, rep(itype)%msubg-1
           rep(itype)%isubg(i+1) = rep(itype)%isubg(i) + rep(itype)%nsubg(i)
         end do
!
         do i = 1, rep(itype)%mgrps
           do j = 1, rep(itype)%ngrps(i)
             io = rep(itype)%igrps(i) + j
             do k = 1, rep(itype)%nsubg(io)
               rep(itype)%subg(rep(itype)%isubg(io)+k) = i 
             end do
           end do
         end do
!
! Setting up ibody and nbody arrays
!
         do i = 1, rep(itype)%mbody-1
           rep(itype)%ibody(i+1) = rep(itype)%ibody(i) + rep(itype)%nbody(i)
         end do
!
         do i = 1, rep(itype)%mbody
           do j = 1, rep(itype)%nbody(i)
             io = rep(itype)%ibody(i) + j
             do k = 1, rep(itype)%ngrps(io)
               rep(itype)%body(rep(itype)%igrps(io)+k) = i
              end do
           end do
         end do
!
       end do
!
! Setting up global iatms and igrps 
!
       igrps(:) = 0
       iatms(:) = 0
       do itype = 2, ntype
         iatms(itype) = iatms(itype-1) + natms(itype-1)
         igrps(itype) = igrps(itype-1) + ngrps(itype-1)
       end do
!
! Fatal error checking
!
       if ( maxat .ne. xtcf%natoms ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Total number of atoms does not match'
         write(*,*) 
         write(*,'(3X,A,I6)') 'Number of atoms in the topology file:   ',maxat
         write(*,'(3X,A,I6)') 'Number of atoms in the trajectory file: ',xtcf%natoms
         write(*,'(3X,A)')    'Please, check your topology file, '//trim(conf)
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       end if
!
! Printing information
!
       if ( debug ) then
!
         write(*,'(A)') 'System size information'
         write(*,'(A)') '-----------------------'
         write(*,*)
         write(*,'(A,20I5)') 'mnode  ',mnode
         write(*,'(A,20I5)') 'mat    ',mat
         write(*,'(A,20I5)') 'nnode  ',nnode(:)
         write(*,'(A,20I5)') 'inode  ',inode(:)
         write(*,*)
         write(*,'(A,20I5)') 'mgrps  ',mgrps
         write(*,'(A,20I5)') 'ngrps  ',ngrps(:)
         write(*,'(A,20I5)') 'igrps  ',igrps(:)
         write(*,*)
         write(*,'(A,20I5)') 'nat    ',nat(:)
         write(*,'(A,20I5)') 'iat    ',iat(:)
         write(*,*)
         write(*,'(A,20I5)') 'matms  ',matms
         write(*,'(A,20I5)') 'natms  ',natms(:)
         write(*,'(A,20I5)') 'iatms  ',iatms(:)
         write(*,*)
!
         do itype = 1, ntype
           write(*,'(A,I4)') 'Molecular system:',itype
           write(*,'(A)')    '.....................'
           write(*,*)
           write(*,'(A,20I4)') 'mbody  ',rep(itype)%mbody
           write(*,'(A,20I4)') 'nbody  ',rep(itype)%nbody
           write(*,'(A,20I4)') 'ibody  ',rep(itype)%ibody
           write(*,'(A,20I4)') 'body   ',rep(itype)%body
           write(*,*) 
           write(*,'(A,20I4)') 'mgrps  ',rep(itype)%mgrps
           write(*,'(A,20I4)') 'ngrps  ',rep(itype)%ngrps
           write(*,'(A,20I4)') 'igrps  ',rep(itype)%igrps
           write(*,'(A,20I4)') 'grps   ',rep(itype)%grps
           write(*,*) 
           write(*,'(A,20I4)') 'msubg  ',rep(itype)%msubg
           write(*,'(A,20I4)') 'nsubg  ',rep(itype)%nsubg
           write(*,'(A,20I4)') 'isubg  ',rep(itype)%isubg
           write(*,'(A,20I4)') 'subg   ',rep(itype)%subg
           write(*,*)       
           write(*,'(A,20I4)') 'matms  ',rep(itype)%matms
           write(*,'(A,20I4)') 'atms   ',rep(itype)%atms
           write(*,'(A,20I4)') 'neiang ',rep(itype)%neiang
           write(*,*)  
           write(*,'(A,20I4)') 'nat    ',rep(itype)%nat
           write(*,'(A,20I4)') 'iat    ',rep(itype)%iat
           write(*,*)
         end do
!
         write(*,'(A)') 'Distance thresholds'
         write(*,'(A)') '-------------------'
         write(*,*)
         write(*,'(25(2X,A5))') (rep(i)%grptag(:rep(i)%mgrps),i=1,ntype)
         do i = 1, mgrps
           write(*,'(25(F6.2,1X))') (thr(i,j),j=1,mgrps)
         end do
         write(*,*)
!
         write(*,'(A)') 'Angles thresholds'
         write(*,'(A)') '-----------------'
         write(*,*)
         write(*,'(25(2X,A5))') (rep(i)%grptag(:rep(i)%mgrps),i=1,ntype)
         do i = 1, mgrps
           write(*,'(25(F6.2,1X))') (thrang(i,j)*180.0/pi,j=1,mgrps)
         end do
         write(*,*)
!
       end if 
!
! Allocating variables depending on topological information 
!
       j = 1
       do i = 1, ntype
         j = j*i
       end do
!
       nmax = 1
       do i = 1, ntype
         nmax = nmax*(msize + i)
       end do
       nmax = nmax/j - 1
!
       allocate(pim(mgrps,mgrps,nmax-1))
       allocate(pop(nmax),conc(nmax),frac(nmax))
       allocate(prob(nmax),num(nmax))
!
       if ( dolife ) allocate(avlife(nmax),nlife(nmax))
!
! Setting distance variables
!
       neidis = neidis**2
!
       thr(:,:)  = thr(:,:)**2
       thr2(:,:) = thr2(:,:)**2   
!
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate)      
!
! Computing the populations of the aggregates
!
       write(*,'(1X,A)') 'Computing the populations of the aggrega'//  & 
                                              'tes along the trajectory'
       write(*,'(1X,A)') 'Please wait, this may take a while...'
       write(*,*)
       CALL FLUSH()
!
       call driver(xtcf,ntype,rep,nnode,inode,nat,iat,natms,iatms,     &
                   ngrps,igrps,mnode,mat,matms,mgrps,thr,thr2,neidis,  &
                   msize,nmax,nsteps,nprint,minstep,maxstep,nsolv,     &
                   avlife,nlife,dopim,schm,scrn,doscrn,dolife,debug)
!
       call system_clock(t1read)        
!
! Closing the trajectory file
!
       call xtcf%close
!
! Averaging pairwise interaction matrix
!  
!       if ( dopim ) then
!         pim(:,:,:) = pim(:,:,:)/nsteps ! FLAG: not need this operation
!
!         do i = 1, msize-1
!           dpaux = 0.0d0
!           do j = 1, mgrps
!             do k = j, mgrps
!               dpaux = dpaux + pim(k,j,i)
!             end do
!           end do
!           pim(:,:,i) = pim(:,:,i)/dpaux*100
!         end do
!       end if
!
! Averaging populations
!
!~        pop(:) = pop(:)/nsteps
!~ !
!~        prob(:) = prob(:)/nsteps
!~ !
!~        frac(:) = frac(:)/nsteps
!~ !
!~        conc(:) = conc(:)/nsteps/(Na*1.0E-24)
!~ !
!~        cin = cin/nsteps/(Na*1.0E-24)
!~ !
!~        volu = volu/nsteps
!
! Averaging lifetimes
!
!~        if ( dolife ) then
!~          do i = 1, nnode
!~            if ( nlife(i) .ne. 0 ) avlife(i) = avlife(i)/nlife(i)
!~          end do 
!~        end if
!
! Printing summary of the results
!
!~        call print_title(6,1,'Output information','-')
!~        write(*,*)
!~        call line_int(6,2,'Number of frames analyzed',lin,':','I12',    &
!~                      nsteps,lfin)
!~        write(*,*)
!~        call line_dp(6,2,'Initial concentration',lin,':','F12.8',       &
!~                     cin,lfin)
!~        call line_dp(6,2,'Average volume',lin,':','F12.6',volu,lfin)
!~        write(*,*)
!~        write(*,'(1X,A,100(X,F6.2))')  'Global populations        : ',  &
!~                                                                   pop(:)
!~        write(*,'(1X,A,100(X,F6.2))')  'Global probabilities      : ',  &
!~                                                                  prob(:)
!~        write(*,'(1X,A,100(X,D12.6))') 'Global fractions          : ',  &
!~                                                                  frac(:)
!~        write(*,'(1X,A,100(X,D12.6))') 'Global concentrations     : ',  &
!~                                                                  conc(:)
!~        write(*,'(1X,A,100(X,D12.6))') 'Number of samples         : ',  &
!~                                                                   num(:)
!~        write(*,*) 
!~        if ( dolife ) then
!~          write(*,'(1X,A,100(X,D12.6))') 'Average lifetimes        '//  &
!~                                                     ' : ',avlife(:msize)
!~          write(*,*) 
!~        end if
!
!       if ( dopim ) then
!         do i = 1, msize-1
!           write(*,'(3X,A,X,I3)') 'Printing PIM for aggregates bel'//  &
!                                                    'onging to type',i+1
!           write(*,'(3X,50("-"))') 
!           write(*,'(8X,20(X,A8))') (adjustr(grptag(j)),j=1,mgrps)
!           do j = 1, mgrps
!             write(*,'(A8,20(X,F8.2))') adjustr(grptag(j)),            &
!                                                      (pim(k,j,i),k=1,j)
!           end do
!           write(*,*)
!         end do
!       end if
!
! Deallocate memory
!
       deallocate(sys,rep)
       deallocate(nnode,inode,nat,iat,natms,iatms,ngrps,igrps)
!
       deallocate(thr,thr2)
       deallocate(thrang)
!
!~        deallocate(pop,conc,frac)
!~        deallocate(prob,num)
!~        deallocate(pim)
!
!~        if ( dolife ) deallocate(avlife,nlife)
!
       close(uniout+1)
       close(uniout+2)
       close(uniout+3)
       close(uniout+4)
       close(uniout+5)
!
! Printing timings
!
       call cpu_time(tfin)
       call system_clock(t2)
!
       call system_clock(t2read)
!
       tcpu  = tfin - tin
       twall = dble(t2-t1)/dble(count_rate)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       call print_time(6,1,'Total CPU time',35,tcpu) 
       call print_speed(6,1,'Total wall time',35,twall,tcpu) 
       write(*,'(1X,68("-"))')
!
       call print_time(6,1,'Total reading time',35,tread)
       call print_speed(6,1,'Total building time',35,tadj,tcpuadj)
!       
       if ( doscrn ) call print_speed(6,1,'Total screening time',35,   &
                                      tscrn,tcpuscrn)
!
       call print_time(6,1,'Total BFS time',35,tbfs)
       call print_time(6,1,'Total sorting time',35,tsort)
!       
       if ( dolife ) call print_speed(6,1,'Total lifetimes calcula'//  &
                                     'tion time',35,tlife,tcpulife)
!
       if ( dopim ) call print_speed(6,1,'Total PIM time',35,tpim,     &
                                     tcpupim)
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
                               maxstep,msize,neidis,schm,scrn,dopim,   &
                               dolife,doscrn,seed,chunkadj,chunkscrn,  &
                               chunklife,weight,nsolv,debug)
!
       use lengths, only: leninp,lenout,lenschm,lencmd,lenarg
!
       use printings
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(out)         ::  traj       !  Trajectory file name
       character(len=lenout),intent(out)         ::  outp       !  Populations file name
       character(len=leninp),intent(out)         ::  inp        !  Groups file name
       character(len=leninp),intent(out)         ::  conf       !  Structure file name
       character(len=lenschm),intent(out)        ::  schm       !  Calculation scheme flag
       character(len=lenschm),intent(out)        ::  scrn       !  Calculation screening flag
       character(len=leninp),intent(out)         ::  weight     !  Weights file name
       logical,intent(out)                       ::  seed       !  Random seed flag
       logical,intent(out)                       ::  dopim      !  PIM calculation flag
       logical,intent(out)                       ::  dolife     !  Lifetimes calculation flag
       logical,intent(out)                       ::  doscrn     !  Screening calculation flag
       real(kind=8),intent(out)                  ::  neidis     !  Screening distance
       integer,intent(out)                       ::  chunkadj   !
       integer,intent(out)                       ::  chunkscrn  !
       integer,intent(out)                       ::  chunklife  !
       integer,intent(out)                       ::  nsolv      !  
       integer,intent(out)                       ::  msize      !  Maximum aggregate size
       integer,intent(out)                       ::  nprint     !  Populations printing steps interval
       integer,intent(out)                       ::  minstep    !  First step for analysis
       integer,intent(out)                       ::  maxstep    !  Last step for analysis
       logical,intent(out)                       ::  debug      !  Debug mode
!
! Local variables
!
       character(len=lencmd)                     ::  cmd        !  Command executed
       character(len=lenarg)                     ::  code       !  Executable name
       character(len=lenarg)                     ::  arg        !  Argument read
       character(len=lenarg)                     ::  next       !  Next argument to be read
       integer                                   ::  io         !  Status
       integer                                   ::  i          !  Index
!
! Setting defaults
!
       inp     = 'aggregate.inp'
       traj    = 'md.xtc'
       conf    = 'conf.gro'
       weight  = '[none]'
       outp    = '[none]'
!
       schm    = 'angles'
       scrn    = 'complete'
!
       doscrn  = .FALSE.
       dolife  = .FALSE.
!
       chunkadj  = 5
       chunkscrn = 1
       chunklife = 1
!
       nprint  = 1
       minstep = 0
       maxstep = 999999999
!
       msize = 10
       nsolv = 10000
!
       neidis  = 1.5d0
!
       dopim   = .TRUE.
!
       seed    = .FALSE.
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
           case ('-schm','--schm','--scheme') 
             call get_command_argument(i,next,status=io)
             call check_arg(next,io,arg,cmd)
             next = lowercase(next)
!
             select case (trim(next))
               case ('dis','dist','distance','distances','original')
                 schm = 'distances'
                 doscrn = .FALSE.
                 dolife = .FALSE.
               case ('ang','angle','angles','restraint','restraints')
                 schm = 'angles'
                 doscrn = .FALSE.
                 dolife = .FALSE.
               case ('dis-scrn','dis-screen','dis-screening',          &
                     'dist-scrn','dist-screen','dist-screening',       &
                     'distance-scrn','distance-screen',                &
                     'distance-screening','distances-scrn',            &
                     'distances-screen','distances-screening')
                 schm   = 'distances'
                 doscrn = .TRUE.
                 dolife = .FALSE.
               case ('ang-scrn','ang-screen','ang-screening',          &
                     'angle-scrn','angle-screen','angle-screening',    &
                     'angles-scrn','angles-screen','angles-screening')
                 schm   = 'angles'
                 doscrn = .TRUE.
                 dolife = .FALSE.
               case ('dis-life','dis-lifetime','dis-lifetimes',        &
                     'dist-life','dist-lifetime','dist-lifetimes',     &
                     'distance-life','distance-lifetime',              &
                     'distance-lifetimes','distances-life',            &
                     'distances-lifetime','distances-lifetimes')
                 schm   = 'distances'
                 doscrn = .FALSE.
                 dolife = .TRUE.
               case ('ang-life','ang-lifetime','ang-lifetimes',        &
                     'angle-life','angle-lifetime','angle-lifetimes',  &
                     'angles-life','angles-lifetime','angles-lifetimes')
                 schm   = 'angles'
                 doscrn = .FALSE.
                 dolife = .TRUE.
               case ('screening-lifetimes','scrnlife','dis-scrnlife',  &
                     'dist-scrnlife','distance-scrnlife',              &
                     'distances-scrnlife','dis-scrn-life',             &
                     'dist-scrn-life','distance-scrn-life',            &
                     'distances-scrn-life')
                 schm   = 'distances'
                 doscrn = .TRUE.
                 dolife = .TRUE.
               case ('ang-scrn-life','angle-scrn-life',                &
                     'angles-scrn-life','ang-scrnlife',                &
                     'angle-scrnlife','angles-scrnlife')
                 schm   = 'angles'
                 doscrn = .TRUE.
                 dolife = .TRUE.
               case default
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)') 'ERROR:  Invalid value introduc'//  &
                                                'ed for --scheme option'
                 write(*,*)
                 write(*,'(3X,A)') 'Unrecognised value     : '//       &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, check the available op'//  &
                                                   'tions in the manual'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-r','--restraint','--restraints') 
             call get_command_argument(i,next,status=io)
             call check_arg(next,io,arg,cmd)
             next = lowercase(next)
!
             select case (trim(next))
               case ('dis','dist','distance','distances')
                 schm = 'distances'
               case ('ang','angle','angles','restraint','restraints')
                 schm = 'angles'
               case default
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)') 'ERROR:  Invalid value introduc'//  &
                                            'ed for --restraints option'
                 write(*,*)
                 write(*,'(3X,A)') 'Unrecognised value     : '//       &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between:  "dist'//  &
                                                    'ances" or "angles"'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-s','--screen-scheme','--screening-scheme') 
             call get_command_argument(i,next,status=io)
             call check_arg(next,io,arg,cmd)
             next = lowercase(next)
!
             select case (trim(next))
               case ('scrnint','int','complete')
                 scrn = 'complete'
               case ('scrncol','col','collision','collisions')
                 scrn = 'collisions'
               case ('scrnosc','osc','oscillation','oscillations')
                 scrn = 'oscillations'
               case default
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)') 'ERROR:  Invalid value introduc'//  &
                                            'ed for --restraints option'
                 write(*,*)
                 write(*,'(3X,A)') 'Unrecognised value     : '//       &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between:  "comp'//  &
                                 'lete", "collisions" or "oscillations"'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-f','-file','--file','-i','-inp','--inp','--input') 
             call get_command_argument(i,inp,status=io)
             call check_arg(inp,io,arg,cmd)
             i = i + 1
           case ('-w','-lw','--log-weights','--log-weight','--weight','--weights') 
             call get_command_argument(i,weight,status=io)
             call check_arg(weight,io,arg,cmd)
             i = i + 1
           case ('-t','-traj','--traj','--trajectory') 
             call get_command_argument(i,traj,status=io)
             call check_arg(traj,io,arg,cmd)
             if ( len_trim(outp) .eq. 0 ) then
               io = index(traj,'.')
               if ( io .eq. 0 ) then
                 outp = trim(traj)
               else 
                 outp = traj(:len_trim(traj)-4)
               end if             
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
             i = i + 1
           case ('-nsolv','--nsolv','--nsolvent')
             call get_command_argument(i,next,status=io)
             read(next,*) nsolv
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
           case ('-cadj','-chunkadj','--chunkadj')
             call get_command_argument(i,next,status=io)
             read(next,*) chunkadj
             i = i + 1
           case ('-cscrn','-chunkscrn','--chunkscrn','--chunkscreening')
             call get_command_argument(i,next,status=io)
             read(next,*) chunkscrn
             i = i + 1
           case ('-clife','-chunklife','--chunklife','--chunklifetimes')
             call get_command_argument(i,next,status=io)
             read(next,*) chunklife
             i = i + 1
           case ('-pim','--pim','--do-pim')
             dopim = .TRUE.
           case ('-nopim','--nopim','--no-pim')
             dopim = .FALSE.
           case ('-life','-lifetime','-lifetimes','--lifetime',        &
                                                          '--lifetimes')
             dolife = .TRUE.
           case ('-nolife','-nolifetime','-nolifetimes',               &
                 '--nolifetime','--nolifetimes','--no-lifetime',       &
                                                       '--no-lifetimes')
             dolife = .FALSE.
           case ('-scrn','-screen','--screen','--screening')
             scrn   = 'complete'
             doscrn = .TRUE.
           case ('-noscrn','-noscreen','--noscreen','--noscreening')
             doscrn = .FALSE.
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
       write(*,'(2X,A)') '-lw,--log-weights     Weights file name'
       write(*,*)
       write(*,'(2X,A)') '-n,--nprint           Printing steps interval'
       write(*,'(2X,A)') '-min,--minimum-step   First step to be a'//  &
                                                               'nalysed'
       write(*,'(2X,A)') '-max,--maximum-step   Last step to be an'//  &
                                                                'alysed'
       write(*,*)
       write(*,'(2X,A)') '-nsolv,--nsolvent     Number of solvent '//  &
                                                             'molecules'
       write(*,'(2X,A)') '-m,--msize            Maximum aggregate size'
       write(*,'(2X,A)') '-d,--neidis           Neighbour list cutoff'
       write(*,*)
       write(*,'(2X,A)') '-schm,--scheme        Aggregates identif'//  &
                                                         'ier algorithm'
       write(*,*)
       write(*,'(2X,A)') '-r,--restraints       Interaction criter'//  &
                                       'ia algorithm [distances,angles]'
       write(*,'(2X,A)') '-s,--screen-scheme    Screening algorith'//  &
                                   'm [complete,collisions,oscillations'
       write(*,*)
       write(*,'(2X,A)') '-[no]life             Compute lifetimes'
       write(*,'(2X,A)') '-[no]scrn             Screen interactions'
       write(*,'(2X,A)') '-[no]pim              Compute pairwise i'//  &
                                                     'nteraction matrix'
       write(*,*)
       write(*,'(2X,A)') '-cadj,--chunkadj      Chunk size for the'//  &
                                         ' adjacency matrix calculation'
       write(*,'(2X,A)') '-cscrn,--chunkscrn    Chunk size for the'//  &
                                                            ' screening'
       write(*,'(2X,A)') '-clife,--chunklife    Chunk size for the'//  &
                                                ' lifetimes calculation'
       write(*,*)
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
       use geometry,   only: sminimgvec,scenvec
       use omp_var
!
       use omp_lib
!
       implicit none
!
! Input/output variables
!
       real(kind=4),dimension(3,natsys),intent(in)  ::  rcoord   !  Input coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3,natms),intent(out)  ::  fcoord   !  Output coordinates !FLAG: kind=8 to kind=4
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
       real(kind=4),dimension(3)                    ::  svaux    !  Auxiliary single precision vector
       integer                                      ::  iinode   !
       integer                                      ::  innode   !
       integer                                      ::  jnnode   !
       integer                                      ::  iisubg   !
       integer                                      ::  insubg   !
       integer                                      ::  j        !
!
! Saving coordinates based on the subgroup-based representation
! -------------------------------------------------------------
!
!$omp parallel do shared(fcoord,rcoord,nsubg,isubg,atms)               &
!$omp             private(atcoord,svaux,iinode,innode,jnnode,iisubg,   &
!$omp                     insubg,j)                                    &
!$omp             schedule(dynamic,chunkadj)
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
             end do
!
             fcoord(:,j) = scenvec(3,nsubg(insubg),                    &
                                   atcoord(:,:nsubg(insubg)))
!           
           else
!
             fcoord(:,j) = rcoord(:,innode+atms(isubg(insubg)+1))
!
           end if    
         end do
!
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine setcoord
!
!======================================================================!
!
! BUILDADJ - BUILD ADJacency matrix
!
! This subroutine builds the adjacency matrix in the molecule-based
!  representation ADJ(NNODE,NNODE) for a given snapshot.
! First, positions of the particles in the physical group-based repre-
!  sentation are obtained, and then the adjacency matrix is computed
!  according to the interaction criteria specified through the subrou-
!  tine BUILDADJMOL
!
       subroutine buildadj(nnode,adj,natms,posi,matms,inposi,nat,thr,  &
                           mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,   &
                           box,neidis,buildadjmol)
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
       integer,dimension(nat),intent(in)            ::  ngrps   !  
       integer,dimension(nat),intent(in)            ::  igrps   ! 
       integer,dimension(nat),intent(in)            ::  nsubg   !  
       integer,dimension(nat),intent(in)            ::  isubg   ! 
       integer,dimension(nat),intent(in)            ::  atms    ! 
       integer,intent(in)                           ::  nnode   !  Number of molecules
       integer,intent(in)                           ::  matms   !
       integer,intent(in)                           ::  natms   ! 
       integer,intent(in)                           ::  nat     ! 
       integer,intent(in)                           ::  msubg   !
       integer,intent(in)                           ::  mgrps   !  Number of subgroups
!
! External functions
!
       external                                     ::  buildadjmol  
!
! Building the adjacency matrix for the current snapshot
!
       call setcoord(nnode,msubg,nsubg,isubg,atms,natms,posi,matms,    &
                     inposi,nat,box)
!
       call buildadjmol(nnode,adj,neidis,msubg,mgrps,nat,thr,ngrps,    &
                        igrps,natms,posi,box)
!
       return
       end subroutine buildadj
!
!======================================================================!
!
! BLOCKDIAG - BLOCK DIAGonalization
!
! This subroutine block diagonalizes an input adjacency matrix 
!  ADJ(NNODE,NNODE) of an undirected unweighted graph of NNODE vertices.
! The subroutine FINDCOMPUNDIR is employed to find the connected compo-
!  nents of the graph using BFS and the QUICKSHORT algorithm is used to
!  find the basis of nodes that block diagonalizes the adjacency matrix
!  sorting the blocks by size and identifier, and then the molecules
!  forming the aggregate according to their canonical order.
!
       subroutine blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg, &
                            nmol,imol,magg)
!
       use timings,     only: tbfs,tcpubfs,tsort,tcpusort,count_rate
       use graphtools,  only: findcompundir
       use sorting,     only: ivvqsort,ivqsort,iqsort
!
       use utils,       only: print_info
!
       use omp_lib
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)   ::  adj       !  Adjacency matrix
       integer,dimension(nnode),intent(out)        ::  mol       !  Molecules identifier
       integer,dimension(nnode),intent(out)        ::  tag       !  Aggregates identifier
       integer,dimension(nnode),intent(out)        ::  agg       !  Aggregates size
       integer,dimension(nnode),intent(out)        ::  nagg      !  Number of aggregates of each size
       integer,dimension(nnode),intent(out)        ::  iagg      !
       integer,dimension(nnode),intent(out)        ::  nmol      !
       integer,dimension(nnode),intent(out)        ::  imol      !
       integer,intent(in)                          ::  nnode     !  Number of molecules
       integer,intent(out)                         ::  nsize     !  Maximum aggregate size
       integer,intent(out)                         ::  magg      !  Number of aggregates
!
! Local variables
! 
       integer                                     ::  i,j,k     !  Indexes
!
! Declaration of time control variables
!
       real(kind=8)                                ::  tinbfs    !  Initial CPU BFS time
       real(kind=8)                                ::  tfinbfs   !  Final CPU BFS time
       real(kind=8)                                ::  tinsort   !  Initial CPU sorting time
       real(kind=8)                                ::  tfinsort  !  Final CPU sorting time
       integer                                     ::  t1bfs     !  Initial BFS time
       integer                                     ::  t2bfs     !  Final BFS time
       integer                                     ::  t1sort    !  Initial sorting time
       integer                                     ::  t2sort    !  Final sorting time  
!
! Block-diagonalizing the adjacency matrix 
!
       call cpu_time(tinbfs)
       call system_clock(t1bfs)     
!
       call findcompundir(nnode,adj,mol,tag,agg,nsize,nagg,magg)
!
       call cpu_time(tfinbfs)
       call system_clock(t2bfs)     
!
       tcpubfs = tcpubfs + tfinbfs - tinbfs
       tbfs    = tbfs    + dble(t2bfs-t1bfs)/dble(count_rate)
!
! Sorting the blocks of the adjacency matrix according to their size,
!  identifier, and canonical order of the constituent molecules
!
       call cpu_time(tinsort)
       call system_clock(t1sort)
! 
! Setting up iagg, imol and nmol arrays
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
! Sorting molecules and aggregate identifiers based on the size of 
!  the aggregates
!
       call ivvqsort(nnode,agg,tag,mol,1,nnode)
!
! Sorting molecules based on their aggregate identifier
!
!~ !$omp parallel do shared(nagg,tag,mol,imol)                            &
!~ !$omp             private(i)                                           &
!~ !$omp             schedule(dynamic,1)
!
       do i = 2, nsize-1
         if ( nagg(i) .gt. 1 ) then
           call ivqsort(nnode,tag,mol,imol(i)+1,imol(i+1))
         end if
       end do
!
!~ !$omp end parallel do                   
!
       if ( nagg(nsize) .gt. 1 )                                       &
                         call ivqsort(nnode,tag,mol,imol(nsize)+1,nnode)
!
! Sorting molecules based on their canonical order
!
!~ !$omp parallel do shared(nagg,mol,imol)                                &
!~ !$omp             private(i,j,k)                                       &
!~ !$omp             schedule(dynamic,1)
!
       do i = 2, nsize
         k = imol(i)
         do j = 1, nagg(i)
           call iqsort(nnode,mol,k+1,k+i)
           k = k + i
         end do
       end do
!
!~ !$omp end parallel do                   
!
       call cpu_time(tfinsort)
       call system_clock(t2sort)     
!
       tcpusort = tcpusort + tfinsort - tinsort
       tsort    = tsort    + dble(t2sort-t1sort)/dble(count_rate)
!
       return
       end subroutine blockdiag
!
!======================================================================!
!
! Seen on: Gaines & Di Tommaso. Pharmaceutics, 2018, 10. 10.3390/pharmaceutics10010012
!~ !
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
