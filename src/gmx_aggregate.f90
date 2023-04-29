!======================================================================!
!
       program aggregate
!
       use xdr, only: xtcfile
!
       use utils, only:  rndmseed
       use input_section
       use datatypes
       use timings
       use printings
       use lengths
       use aggtools
!
       use omp_lib
!
       implicit none
!
       include 'inout.h'
!
! System information
!
       type(groinp)                                    ::  sys      !  Monomer information
       type(xtcfile)                                   ::  xtcf     !  Trajectory information
!
! Input/output files
!
       character(len=leninp)                           ::  traj     !  Trajectory file name
       character(len=leninp)                           ::  conf     !  Configuration file name
       character(len=leninp)                           ::  inp      !  General input file name
       character(len=lenout)                           ::  outp     !  Output file name
!
! Interaction criteria information
! 
       real(kind=8),dimension(:,:),allocatable         ::  thr      !  Distance threshold
       real(kind=8),dimension(:,:),allocatable         ::  thr2     !  Distance threshold
       real(kind=8)                                    ::  neidis   !  Screening distance
!
! Average properties
!
       real(kind=8),dimension(:,:,:),allocatable       ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(:),allocatable           ::  pop      !  Populations
       real(kind=8),dimension(:),allocatable           ::  conc     !  Concentrations
       real(kind=8),dimension(:),allocatable           ::  frac     !  Molar fractions
       real(kind=8)                                    ::  cin      !  Stechiometric concentration
       real(kind=8)                                    ::  volu     !  Simulation box volume
!
! Aggregates information in the molecule-based representation
!
       integer                                         ::  nnode    !  Total number of molecules
       integer                                         ::  msize    !  Maximum aggregate size
!
! Topological representations information
!
       character(len=leninp)                           ::  tgrp     !  Groups file title
       character(len=lentag),dimension(:),allocatable  ::  grptag   !  Names of the groups
       integer,dimension(:),allocatable                ::  body     !  Number of groups in each body
       integer,dimension(:),allocatable                ::  nbody    !  Number of groups in each body
       integer,dimension(:),allocatable                ::  ibody    !  Number of groups in each body
       integer,dimension(:),allocatable                ::  grps     !  Number of subgroups in each group
       integer,dimension(:),allocatable                ::  ngrps    !  Number of subgroups in each group
       integer,dimension(:),allocatable                ::  igrps    !  Number of subgroups in each group
       integer,dimension(:),allocatable                ::  subg     !  Number of atoms in each subgroup
       integer,dimension(:),allocatable                ::  nsubg    !  Number of atoms in each subgroup
       integer,dimension(:),allocatable                ::  isubg    !  Number of atoms in each subgroup
       integer,dimension(:),allocatable                ::  atms     !  Atoms identifier
       integer                                         ::  natms    !  Total number of subgroups in the system
       integer                                         ::  mbody    !  Number of bodies
       integer                                         ::  mgrps    !  Number of groups
       integer                                         ::  msubg    !  Number of subgroups
       integer                                         ::  matms    !  Number of interacting atoms in the monomer
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
       real(kind=8),dimension(:),allocatable           ::  avlife    !
       integer,dimension(:),allocatable                ::  nlife     !
!
! AnalysisPhenolMD variables
!     
       real(kind=8),dimension(9,9,3)                   ::  table    !  Number of aggregates of each type
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
       character(len=lentag)                           ::  schm     !  Calculation scheme flag
       logical                                         ::  dopim    !  PIM calculation flag
       logical                                         ::  dolife   !  Lifetimes calculation flag
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
!
       real(kind=8),parameter                          ::  Na = 6.022140760E+23
!
! Printing header
!
       write(*,*)
       call print_start()
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
! Initializing line lengths
!
       lin  = 45
       lfin = 90
!
       nsolv = 10000
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
                         msize,neidis,schm,dopim,dolife,seed,debug)
!
! Initializing random number generator
!
       call rndmseed(seed)
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
       call line_str(6,2,'Calculation scheme',lin,':',trim(schm),lfin)
       call line_log(6,2,'Lifetimes calculation',lin,':',dolife,lfin)
       write(*,*)
!
       call line_dp(6,2,'Screening distance',lin,':','F5.2',           &
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
       call read_gro(conf,sys)               !  TODO: Input .gro .top .xyz ...
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
         write(*,*) 'Incorrect extension!'   !  TODO: Handle .xtc .trr .gro .xyz ...
         write(*,*)
         call print_end()
       end if
!
! Allocating variables depending on system size
!
       nnode  = xtcf%NATOMS/sys%nat
!
       allocate(thr(sys%nat,sys%nat))
       allocate(thr2(sys%nat,sys%nat))
!
       allocate(nbody(sys%nat),ngrps(sys%nat),nsubg(sys%nat))
       allocate(ibody(sys%nat),igrps(sys%nat),isubg(sys%nat))
       allocate(body(sys%nat),grps(sys%nat),subg(sys%nat))
       allocate(grptag(sys%nat),atms(sys%nat))
!
       if ( dolife ) allocate(avlife(nnode),nlife(nnode))
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
!
! Setting up igrps and ngrps arrays
!
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
!
! Setting up isubg and nsubg arrays
!
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
!
! Setting up ibody and nbody arrays
!
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
       if ( debug ) then
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
         write(*,'(A,20I3)') 'atms ',matms
         write(*,'(A,20I3)') 'atms ',atms
         write(*,*)  
       end if     
!
! Allocating variables depending on topological information 
!
       allocate(pim(mgrps,mgrps,msize-1))
       allocate(pop(msize),conc(msize),frac(msize))
!
! Setting distance variables
!
       neidis = neidis**2
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
       call system_clock(t2read)
!
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       select case (trim(schm))
         case ('original')
!
           call aggdist(xtcf,sys%nat,nnode,natms,thr,thr2,neidis,pim,  &
                        msize,pop,conc,frac,cin,volu,nsteps,nbody,     &
                        ngrps,nsubg,ibody,igrps,isubg,body,grps,subg,  &
                        atms,mbody,mgrps,msubg,matms,nprint,minstep,   &
                        maxstep,nsolv,dopim,debug)
!
         case ('screening')
!
           if ( dolife ) then
             call aggscrnlife(xtcf,sys%nat,nnode,natms,thr,thr2,       &
                              neidis,pim,msize,pop,conc,frac,cin,      &
                              volu,nsteps,nbody,ngrps,nsubg,ibody,     &
                              igrps,isubg,body,grps,subg,atms,mbody,   &
                              mgrps,msubg,matms,nprint,minstep,        &
                              maxstep,nsolv,avlife,nlife,dopim,debug)
           else
             call aggscrn(xtcf,sys%nat,nnode,natms,thr,thr2,neidis,    &
                          pim,msize,pop,conc,frac,cin,volu,nsteps,     &
                          nbody,ngrps,nsubg,ibody,igrps,isubg,body,    &
                          grps,subg,atms,mbody,mgrps,msubg,matms,      &
                          nprint,minstep,maxstep,nsolv,dopim,debug)
           end if
!
       end select
!
       call system_clock(t1read)        
!
! Closing the trajectory file
!
       call xtcf%close
!
! Averaging pairwise interaction matrix
!  
!~        if ( dopim ) then
!~          pim(:,:,:) = pim(:,:,:)/nsteps ! FLAG: not need this operation
!~ !
!~          do i = 1, msize-1
!~            dpaux = 0.0d0
!~            do j = 1, mgrps
!~              do k = j, mgrps
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
!
! Averaging lifetimes
!
       if ( dolife ) then
         do i = 1, nnode
           if ( nlife(i) .ne. 0 ) avlife(i) = avlife(i)/nlife(i)
         end do 
      end if
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
       call print_title(6,1,'Output information','-')
       write(*,*)
       call line_int(6,2,'Number of frames analyzed',lin,':','I12',    &
                     nsteps,lfin)
       write(*,*)
       call line_dp(6,2,'Initial concentration',lin,':','F12.8',       &
                    cin,lfin)
       call line_dp(6,2,'Average volume',lin,':','F12.6',volu,lfin)
       write(*,*)
       write(*,'(1X,A,100(X,F6.2))')  'Global populations        : ',  &
                                                                  pop(:)
       write(*,'(1X,A,100(X,D12.6))') 'Global fractions          : ',  &
                                                                 frac(:)
       write(*,'(1X,A,100(X,D12.6))') 'Global concentrations     : ',  &
                                                                 conc(:)
       write(*,*) 
       if ( dolife ) then
         write(*,'(1X,A,100(X,D12.6))') 'Average lifetimes        '//  &
                                                    ' : ',avlife(:msize)
         write(*,*) 
       end if
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
!
!~        if ( dopim ) then
!~          do i = 1, msize-1
!~            write(*,'(3X,A,X,I3)') 'Printing PIM for aggregates bel'//  &
!~                                                     'onging to type',i+1
!~            write(*,'(3X,50("-"))') 
!~            write(*,'(8X,20(X,A8))') (adjustr(grptag(j)),j=1,mgrps)
!~            do j = 1, mgrps
!~              write(*,'(A8,20(X,F8.2))') adjustr(grptag(j)),            &
!~                                                       (pim(k,j,i),k=1,j)
!~            end do
!~            write(*,*)
!~          end do
!~        end if
!
! Deallocate memory
!
       deallocate(nbody,ngrps,nsubg)
       deallocate(ibody,igrps,isubg)
       deallocate(body,grps,subg)
       deallocate(grptag,atms)
!
       deallocate(pop,conc,frac)
       deallocate(pim)
!
       if ( dolife ) deallocate(avlife,nlife)
!
       close(uniout+1)
       close(uniout+2)
       close(uniout+3)
!
! Printing timings
!
       call cpu_time(tfin)
!
       call system_clock(t2)
       call system_clock(t2read)
!
       tcpu = tfin - tin
!
       twall = dble(t2-t1)/dble(count_rate)
       tread = tread + dble(t2read-t1read)/dble(count_rate) 
!
       call print_time(6,1,'Total CPU time',35,tcpu) 
       call print_speed(6,1,'Total wall time',35,twall,tcpu) 
       write(*,'(1X,68("-"))')
!
       call print_time(6,1,'Total reading time',35,tread)
       call print_speed(6,1,'Total building time',35,tadj,tcpuadj)
!       
       if ( (trim(schm).eq.'screening') .or.                           &
                                       (trim(schm).eq.'scrnlife') ) then
         call print_speed(6,1,'Total screening time',35,tscrn,tcpuscrn)
       end if
!
       call print_speed(6,1,'Total BFS time',35,tbfs,tcpubfs)
       call print_speed(6,1,'Total sorting time',35,tsort,tcpusort)
!       
       if ( dolife ) call print_speed(6,1,'Total lifetimes calcula'//  &
                                     'tion time',35,tlife,tcpulife)
!
       if ( dopim ) call print_time(6,1,'Total PIM time',35,tpim)
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
                               maxstep,msize,neidis,schm,dopim,dolife, &
                               seed,debug)
!
       use printings
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
       character(len=lentag),intent(out)         ::  schm     !  Calculation scheme flag
       logical,intent(out)                       ::  seed     !  Random seed flag
       logical,intent(out)                       ::  dopim    !  PIM calculation flag
       logical,intent(out)                       ::  dolife   !  Lifetimes calculation flag
       real(kind=8),intent(out)                  ::  neidis   !  Screening distance
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
       outp    = ''
!
       schm    = 'screening'
!
       nprint  = 1
       minstep = 0
       maxstep = 999999999
!
       msize  = 10
!
       neidis  = 1.5d0
!
       dopim   = .TRUE.
       dolife  = .FALSE.
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
           case ('-s','--scheme') 
             call get_command_argument(i,next,status=io)
             call check_arg(next,io,arg,cmd)
             next = lowercase(next)
!
             select case (trim(next))
               case ('dist','distance','distances','original')
                 schm = 'original'
               case ('scrn','screen','screening')
                 schm   = 'screening'
                 dolife = .FALSE.
               case ('screening-lifetimes','scrnlife')
                 schm   = 'screening'
                 dolife = .TRUE.
               case default
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)') 'ERROR:  Invalid value introduc'//  &
                                                 'ed for --scheme option'
                 write(*,*)
                 write(*,'(3X,A)') 'Unrecognised value     : '//       &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between:  "orig'//  &
                                      'inal", "screening" or "scrnlife"'
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
           case ('-life','-lifetime','-lifetimes','--lifetime',        &
                                                          '--lifetimes')
             dolife = .TRUE.
           case ('-nolife','-nolifetime','-nolifetimes',               &
                 '--nolifetime','--nolifetimes','--no-lifetime',       &
                                                       '--no-lifetimes')
             dolife = .FALSE.
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
       write(*,*)
       write(*,'(2X,A)') '-n,--nprint           Printing steps interval'
       write(*,'(2X,A)') '-min,--minimum-step   First step to be a'//  &
                                                               'nalysed'
       write(*,'(2X,A)') '-max,--maximum-step   Last step to be an'//  &
                                                                'alysed'
       write(*,*)
       write(*,'(2X,A)') '-m,--msize            Maximum aggregate size'
       write(*,'(2X,A)') '-d,--neidis           Neighbour list cutoff'
       write(*,*)
       write(*,'(2X,A)') '-s,--scheme           Aggregates identif'//  &
                                                         'ier algorithm'
       write(*,'(2X,A)') '-[no]pim              Compute pairwise i'//  &
                                                     'nteraction matrix'
       write(*,'(2X,A)') '-[no]life             Compute lifetimes'
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
       use geometry,   only: sminimgvec,                               &
                             scenvec
!
       use omp_lib
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
!$omp             schedule(dynamic,5)
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
             fcoord(:,j) = scenvec(3,nsubg(insubg),                    &
                                   atcoord(:,:nsubg(insubg)))
!           
!
           else
!
             fcoord(:,j) = rcoord(:,innode+atms(isubg(insubg)+1))
!
!~ write(*,*) iinode,innode,innode+atms(isubg(insubg)+1)
           end if    
         end do
!~ write(*,*)
       end do
!
!$omp end parallel do                   
!
       return
       end subroutine setcoord
!
!======================================================================!
!
!~        subroutine print_coord(xtcf,sys,outp,msize,nagg,nnode,mol,agg)
!~ !
!~        use xdr,       only: xtcfile
!~        use geometry,  only: sminimgvec
!~        use datatypes
!~ !
!~        implicit none
!~ !
!~        include 'inout.h'
!~ !
!~ ! Input/output variables
!~ !
!~        type(xtcfile),intent(inout)           ::  xtcf    !  xtc file informacion   
!~        type(groinp),intent(in)               ::  sys     !  System information
!~        character(len=lenout),intent(in)      ::  outp    !  Output file name
!~        integer,dimension(msize),intent(in)  ::  nagg    !  Number of aggregates of each size
!~        integer,dimension(nnode),intent(in)   ::  mol    !  Molecule identifier
!~        integer,dimension(nnode),intent(in)   ::  agg    !  Aggregates size
!~        integer,intent(in)                    ::  msize  !  Maximum aggregate size
!~        integer,intent(in)                    ::  nnode   !  Number of residues
!~ !
!~ ! Local variables
!~ !
!~        type(xtcfile)                         ::  xtco    !  xtc file informacion
!~        character(len=lenout)                 ::  straux  !  Auxiliary string
!~        character(len=lenout)                 ::  aux     !  Auxiliary string
!~        real(kind=4),dimension(3)             ::  svaux   !  Auxiliary single precision vector
!~        real(kind=8),dimension(3)             ::  cofm    !  Center of mass vector
!~        real(kind=8)                          ::  mass    !  Total mass of the aggregate
!~        integer                               ::  nsize   !  Size of the previous printed aggregate
!~        integer                               ::  i,j     !  Indexes
!~        integer                               ::  m,n     !  Indexes
!~        integer                               ::  p,q,r,s !  Indexes
!~ !
!~ ! Priting coordinates of the aggregates
!~ ! -------------------------------------
!~ !
!~        write(aux,*) xtcf%STEP
!~        aux = adjustl(aux)
!~ !
!~        straux = trim(outp)//'_'//trim(aux)//'.xtc'
!~        aux    = trim(outp)//'_'//trim(aux)//'.xyz'
!~ ! Printing global coordinates in xtc format
!~        call xtco%init(straux,'w')
!~        call xtco%write(xtcf%natoms,xtcf%step,xtcf%time,xtcf%box,       &
!~                                                      xtcf%pos,xtcf%prec)
!~        call xtco%close
!~ ! Printing global coordinates in xyz format
!~        open(unit=uniinp,file=trim(aux),action='write')
!~ !
!~        write(uniinp,*) xtcf%NATOMS
!~        write(uniinp,*) sys%title
!~ !
!~        do i = 1, xtcf%NATOMS
!~          write(uniinp,*) sys%atname(modulo(i-1,sys%nat)+1),            &
!~                                                         xtcf%pos(:,i)*10
!~        end do
!~ !
!~        close(uniinp)
!~ ! Printing coordinates of the complexes by size
!~        n     = nagg(1)
!~        nsize = 0
!~ !
!~        do i = 2, msize
!~          do j = 1, nagg(i)
!~            n = n + agg(n)
!~            if ( nsize .ne. agg(n) ) then
!~              write(straux,*) agg(n)
!~              straux = adjustl(straux)
!~ ! 
!~              write(aux,*) xtcf%STEP
!~              aux = adjustl(aux)
!~ !
!~              straux = trim(outp)//'_'//trim(aux)//'_'//trim(straux)//  &
!~                                                                   '.xyz'
!~              open(unit=uniinp,file=trim(straux),action='write')
!~ !
!~              nsize = agg(n)
!~            end if
!~ !
!~            cofm(:)  = 0.0d0
!~            mass     = 0.0d0
!~            svaux(:) = xtcf%pos(:,(mol(n)-1)*sys%nat+1) 
!~ !
!~            m = n - 1
!~            do p = 1, agg(n)
!~              mass = mass + sys%totm
!~              m = m + 1
!~              r = (mol(m)-1)*sys%nat
!~              do q = 1, sys%nat
!~                r = r + 1 
!~                xtcf%pos(:,r) = sminimgvec(svaux(:),xtcf%pos(:,r),      &
!~                                           (/ xtcf%box(1,1),            &
!~                                              xtcf%box(2,2),            &
!~                                              xtcf%box(3,3) /) )
!~                xtcf%pos(:,r) = svaux(:) + xtcf%pos(:,r)
!~                do s = 1, 3
!~                  cofm(s) = cofm(s) + sys%mass(q)*xtcf%pos(s,r)
!~                end do
!~              end do
!~            end do
!~            cofm(:) = cofm(:)/mass
!~ !
!~            write(uniinp,*) sys%nat*agg(n)
!~            write(uniinp,'(20(X,I5))') mol(n:n+agg(n)-1)
!~ !
!~            m = n - 1
!~            do p = 1, agg(n)
!~              m = m + 1
!~              r = (mol(m)-1)*sys%nat
!~              do q = 1, sys%nat
!~                r = r + 1
!~                write(uniinp,'(A5,3(1X,F12.8))') sys%atname(q),         &
!~                                             (xtcf%pos(:,r) - cofm(:))*10
!~              end do
!~            end do
!~ !
!~            if ( nsize .ne. agg(n+1) ) close(uniinp)
!~ !
!~          end do
!~        end do
!~ !
!~        close(uniinp)
!~ !
!~ ! r  -> i
!~ ! q  -> iat
!~ ! p  -> irenum
!~ ! n  -> iitag
!~ ! m  -> iimol
!~ ! i  -> itype
!~ ! j  -> inagg
!~ !
!~        return
!~        end subroutine print_coord
!
!======================================================================!
!
! BUILDADJ - BUILD ADJacency matrix
!
! This subroutine builds the adjacency matrix in the molecule-based
!  representation ADJ(NNODE,NNODE) for a given snapshot.
! First, positions of the particles in the physical group-based repre-
!  sentation are obtained, and then the adjacency matrix is computed.
!
       subroutine buildadj(nnode,adj,natms,posi,matms,inposi,nat,thr,  &
                           mgrps,ngrps,igrps,msubg,nsubg,isubg,atms,   &
                           box,neidis)
!
       use graphtools, only: buildadjmolbub
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
! Building the adjacency matrix for the current snapshot
!
       call setcoord(nnode,msubg,nsubg,isubg,atms,natms,posi,          & 
                     matms,inposi,nat,box)
!
       call buildadjmolbub(nnode,adj,neidis,msubg,mgrps,nat,thr,       &
                           ngrps,igrps,natms,posi,box)
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
       use timings,    only: tbfs,tcpubfs,tsort,count_rate
       use graphtools, only: findcompundir
       use sorting,    only: ivvqsort,                                 &
                             ivqsort,                                  &
                             iqsort
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
       tbfs = tbfs + dble(t2bfs-t1bfs)/dble(count_rate)
!
! Sorting the blocks of the adjacency matrix according to their size,
!  identifier, and canonical order of the constituent molecules
!
       call system_clock(t1sort)         ! TODO: parallelize this region
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
!$omp parallel do shared(nagg,tag,mol,imol)                            &
!$omp             private(i)                                           &
!$omp             schedule(dynamic,1)
!
       do i = 2, nsize-1
         if ( nagg(i) .gt. 1 ) then
           call ivqsort(nnode,tag,mol,imol(i)+1,imol(i+1))
         end if
       end do
!
!$omp end parallel do                   
!
       if ( nagg(nsize) .gt. 1 )                                       &
                         call ivqsort(nnode,tag,mol,imol(nsize)+1,nnode)
!
! Sorting molecules based on their canonical order
!
!$omp parallel do shared(nagg,mol,imol)                                &
!$omp             private(i,j,k)                                       &
!$omp             schedule(dynamic,1)
!
       do i = 2, nsize
         k = imol(i)
         do j = 1, nagg(i)
           call iqsort(nnode,mol,k+1,k+i)
           k = k + i
         end do
       end do
!
!$omp end parallel do                   
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
