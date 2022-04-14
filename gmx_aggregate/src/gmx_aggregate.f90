!======================================================================!
!
       program aggregate

! Use the xdr interface
       use xdr, only: xtcfile!,trrfile
!
       use input_section
       use graphtools
       use datatypes
       use utils
!
       implicit none
!
       include 'parameters.h'
       include 'timings.h'
       include 'info.h'
!
       type(groinp)                               ::  sys      !  Monomer information
       character(len=leninp)                      ::  traj     !  Trajectory file name
       character(len=leninp)                      ::  conf     !  Configuration file name
       character(len=leninp)                      ::  inp      !  General input file name
       character(len=lenout)                      ::  outp     !  Populations file name
       character(len=leninp)                      ::  tgrp     !  Groups file title
       logical,dimension(:,:),allocatable         ::  adj      !  Adjacency matrix
       logical                                    ::  dopim    !  PIM calculation flag
       real(kind=8),dimension(:,:,:),allocatable  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(:,:),allocatable    ::  thr      !  Distance threshold
       real(kind=8),dimension(:),allocatable      ::  pop      !  Populations
       real(kind=8),dimension(9,9,3)              ::  table    !  Number of aggregates of each type
       real(kind=8)                               ::  maxdis   !  Screening distance
       real(kind=8)                               ::  dpaux    !  Auxiliary double precision number
       integer,dimension(:),allocatable           ::  nmol     !  Number of aggregates of each size
       integer,dimension(:),allocatable           ::  imol     !  Molecules identifier
       integer,dimension(:),allocatable           ::  iagg     !  Aggregates identifier
       integer,dimension(:),allocatable           ::  itag     !  Aggregates size 
       character(len=8),dimension(:),allocatable  ::  grptag   !  Names of the groups
       integer,dimension(:),allocatable           ::  body     !  Number of groups in each body
       integer,dimension(:),allocatable           ::  grps     !  Number of subgroups in each group
       integer,dimension(:),allocatable           ::  subg     !  Number of atoms in each subgroup
       integer,dimension(:),allocatable           ::  igrps    !  Atoms identifier
       integer                                    ::  nbody    !  Number of bodies
       integer                                    ::  ngrps    !  Number of groups
       integer                                    ::  nsubg    !  Number of subgroups
       integer                                    ::  nagg     !  Number of chemical species
       integer                                    ::  nnode    !  Total number of molecules
       integer                                    ::  nprint   !  Populations printing interval
       integer                                    ::  minstep  !  First step for analysis
       integer                                    ::  maxstep  !  Last step for analysis
       integer                                    ::  nsteps   !  Number of snapshots analyzed
       integer                                    ::  maxagg   !  Maximum aggregate size
       integer                                    ::  io       !  Status
       integer                                    ::  i,j,k    !  Indexes
       logical                                    ::  seed     !  Random seed flag
       logical                                    ::  debug    !  Debug mode
! Declaration of a variable of type xtcfile
       type(xtcfile)                              ::  xtcf     !  xtc file informacion
! Declaration of a variable of type trrfile
!~        type(trrfile)                              ::  trr      !  trr file information
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
! Initializing timings
!
       tcpu  = 0.0d0
       tread = 0.0d0
       tadj  = 0.0d0
       tbfs  = 0.0d0
       tsort = 0.0d0
       tpim  = 0.0d0
!
! Initializing random number generator
!
       call rndmseed(seed)
!
! Reading command line options
!
       call command_line(traj,conf,inp,outp,nprint,minstep,maxstep,    & 
                         maxagg,maxdis,dopim,seed,debug)
!
!  General settings and fatal errors check
!
       io = index(traj,'.')
       if ( io .eq. 0 )  traj = trim(traj)//'.xtc'
!
       io = index(conf,'.')
       if ( io .eq. 0 ) conf = trim(conf)//'.gro'
!
       io = index(outp,'.')
       if ( io .eq. 0 ) outp = trim(outp)//'.dat'
!
       maxdis = maxdis**2
!
! Printing summary of the input information
!
       write(*,'(1X,A)') 'Input information'
       write(*,'(1X,17("-"))')
       write(*,*)
       write(*,'(2X,A,2X,A)')    'General input file name    :',       &
                                                               trim(inp)
       write(*,'(2X,A,2X,A)')    'Trajectory file name       :',       &
                                                              trim(traj)
       write(*,'(2X,A,2X,A)')    'Structure file name        :',       & 
                                                              trim(conf)
       write(*,'(2X,A,2X,A)')    'Populations file name      :',       &
                                                              trim(outp)
       write(*,*)
       write(*,'(2X,A,9X,F4.2)') 'Screening distance         :',       &
                                                            sqrt(maxdis)  
       write(*,'(2X,A,9X,I4)')   'Maximum aggregate size     :',maxagg
       write(*,'(2X,A,2X,I11)')  'First step for analysis    :',minstep
       write(*,'(2X,A,2X,I11)')  'Last step for analysis     :',maxstep
       write(*,'(2X,A,2X,I11)')  'Printing interval          :',nprint
       if ( debug ) then
         write(*,'(2X,A,11X,A)') 'Debug mode                 :','ON'
       else
         write(*,'(2X,A,10X,A)') 'Debug mode                 :','OFF'
       end if
       write(*,*)
!
! Processing Gromacs input file
!
       call read_gro(conf,sys)
! Opening populations output file
       open(unit=uniout,file=trim(outp),action='write')
!
       if ( traj(len_trim(traj)-3:) .eq. '.xtc' ) then 
! Initialize it with the names of xtc files you want to read in and write out
         call xtcf%init(trim(traj))
! Read in each configuration. Everything is stored in the xtcfile type
         call xtcf%read
!
         nnode  = xtcf%NATOMS/sys%nat
!
         allocate(thr(sys%nat,sys%nat))
         allocate(adj(nnode,nnode),imol(nnode),iagg(nnode),            &
                  itag(nnode),nmol(maxagg),pop(maxagg))
         allocate(grptag(sys%nat),body(sys%nat),grps(sys%nat),         &
                  subg(sys%nat),igrps(sys%nat))
!
! Processing general input file
!
         call read_inp(inp,sys%nat,tgrp,grptag,body,grps,subg,igrps,   &
                       nbody,ngrps,nsubg,thr)
!
         allocate(pim(ngrps,ngrps,maxagg-1))
!
         pim(:,:,:) = 0
         thr(:,:)   = thr(:,:)**2
!
         call system_clock(t2read)
!
         tread = tread + dble(t2read-t1read)/dble(count_rate)         
!
! Computing the populations of the aggregates
!
         write(*,'(4X,A)') 'Computing the populations of the aggre'//  & 
                                            'gates along the trajectory'
         write(*,'(4X,A)') 'Please wait, this may take a while...'
         write(*,*)
!
         nsteps = 0
         pop(:) = 0.0d0
!
         table(:,:,:) = 0.0d0
!
! Analyzing frames in the inverval [minstep,maxstep]
!
         do while ( (xtcf%STAT.eq.0) .and. (xtcf%STEP.le.maxstep)      &
                                          .and. (xtcf%STEP.ge.minstep) )
           if ( mod(xtcf%STEP,nprint) .eq. 0 ) then
! Update configurations counter
             nsteps = nsteps + 1
!
! Building the adjacency matrix for the current snapshot
!
             call system_clock(t1adj)     
!
             call buildadjmol(ngrps,grps,nsubg,subg,nnode,igrps,adj,   &
                              thr,maxdis,xtcf%NATOMS,xtcf%pos,sys%nat, &
                              sys%mass,(/xtcf%box(1,1),xtcf%box(2,2),  &
                              xtcf%box(3,3)/),debug)
!
             call system_clock(t2adj) 
!
             tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)    
!
! Block-diagonalizing the adjacency matrix 
!
             call blockdiag(nnode,adj,imol,iagg,itag,maxagg,nmol,      &
                            nagg,tbfs,tsort,count_rate,debug)
!
! Printing the population of every aggregate
!
             write(uniout,'(I10,20(X,F6.2))') xtcf%STEP,               &
                                                  real(nmol(:))/nagg*100
!
             pop(:) = pop(:) + real(nmol(:))/nagg*100
!
! Printing summary of the results
! 
             if ( debug ) then
               write(*,'(X,A,X,I9)')      'Step number',xtcf%STEP
               write(*,'(X,21("="))') 
               write(*,'(2X,A,X,I6)')     'Total number of entitie'//  &
                                                             's : ',nagg
               write(*,'(2X,A,20(X,I6))') 'Aggregates of each type'//  &
                                                             '  : ',nmol
               write(*,*)
!
               call print_info(nmol(1),nnode-nmol(1),                  &
                               imol(nmol(1)+1:),iagg(nmol(1)+1:),      &
                               itag(nmol(1)+1:),'imol','iagg','itag')
!
               call print_coord(xtcf,sys,outp,maxagg,                  &
                                nmol,nnode,imol,itag)
             end if
!
! Adding up the pairwise interaction matrix of the current snapshot
! 
             if ( dopim ) then
               call system_clock(t1pim)     
!
               call build_pim(ngrps,grps,nsubg,subg,sys%nat,igrps,thr, &
                              nnode,imol,itag,maxagg,nmol,xtcf%NATOMS, &
                              xtcf%pos,(/xtcf%box(1,1),xtcf%box(2,2),  &
                              xtcf%box(3,3)/),pim,debug)
!
               call system_clock(t2pim)     
!
               tpim = tpim + dble(t2pim-t1pim)/dble(count_rate)   
             end if      
!
!
! Analyzing aggregates by their connectivity
!
             call analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,  &
                              sys%nat,igrps,thr,nnode,adj,imol,itag,   & 
                              maxagg,nmol,xtcf%NATOMS,xtcf%pos,        &
                              (/xtcf%box(1,1),xtcf%box(2,2),           &
                                xtcf%box(3,3)/),debug)
!
           end if
!
! Reading information of the next snapshot
!
           call system_clock(t1read)
!
           call xtcf%read
!
           call system_clock(t2read)
!
           tread = tread + dble(t2read-t1read)/dble(count_rate)         
!
         end do
! Close the file
         call xtcf%close
!
       else
         write(*,*) 'Incorrect extension!'
         write(*,*)
         call print_end()
       end if
!
! Averaging pairwise interaction matrix
!  
       if ( dopim ) then
         pim(:,:,:) = pim(:,:,:)/nsteps
!
         do i = 1, maxagg-1
           dpaux = 0.0d0
           do j = 1, ngrps
             do k = j, ngrps
               dpaux = dpaux + pim(k,j,i)
             end do
           end do
           pim(:,:,i) = pim(:,:,i)/dpaux*100
         end do
       end if
!
! Averaging populations
!
       pop(:) = pop(:)/nsteps
!
       table(:,:,:) = table(:,:,:)/nsteps
!
       do i = 1, 3
         dpaux = 0.0d0
         do j = 1, 9
           do k = 1, 9
             dpaux = dpaux + table(k,j,i)
           end do
         end do
         table(:,:,i) = table(:,:,i)/dpaux*100
       end do
!
! Printing summary of the results
!
       write(*,'(1X,A)') 'Output information'
       write(*,'(1X,18("-"))')
       write(*,*)
       write(*,'(1X,A,3X,I11)')     'Number of frames analyzed : ',    &
                                                                  nsteps
       write(*,'(1X,A,20(X,F6.2))') 'Global populations        : ',    &
                                                                  pop(:)
       write(*,*)
!
       do i = 1, 3
         write(*,'(3X,A,X,I3)') 'Populations of the aggregates belonging to type',i+1
         write(*,'(3X,51("-"))') 
         write(*,'(8X,A5)')          'IDXOH'
         write(*,'(1X,A5,20(X,I6))') 'IDXPh',(k,k=1,9)
         do j = 1, 9
           write(*,'(I6,20(X,F6.2))') j,(table(k,j,i),k=1,9)
         end do
         write(*,*)
       end do
!
       if ( dopim ) then
         do i = 1, maxagg-1
           write(*,'(3X,A,X,I3)') 'Printing PIM for aggregates bel'//  &
                                                    'onging to type',i+1
           write(*,'(3X,50("-"))') 
           write(*,'(8X,20(X,A8))') (adjustr(grptag(j)),j=1,ngrps)
           do j = 1, ngrps
             write(*,'(A8,20(X,F8.2))') adjustr(grptag(j)),            &
                                                      (pim(k,j,i),k=1,j)
           end do
           write(*,*)
         end do
       end if
!
! Deallocate memory
!
       deallocate(thr)
       deallocate(adj,imol,iagg,itag,nmol,pop)
       deallocate(grptag,body,grps,subg,igrps)
       deallocate(pim)
!
       close(uniout)
!
! Printing timings
!
       call system_clock(t2)            ! FLAG: 3 seconds lost
!
       tcpu = dble(t2-t1)/dble(count_rate)
!
       write(*,'(1X,A,3(X,I2,X,A))') 'Total CPU time           ',      &
                                      int(tcpu/(60*60)),'hours',       &
                                      mod(int(tcpu/60),60),'minutes',  &
                                      mod(int(tcpu),60),'seconds'  
       write(*,'(1X,56("-"))')
!
       write(*,'(1X,A,3(X,I2,X,A))') 'Total reading time       ',      &
                                      int(tread/(60*60)),'hours',      &
                                      mod(int(tread/60),60),'minutes', &
                                      mod(int(tread),60),'seconds'  
       write(*,'(1X,A,3(X,I2,X,A))') 'Total building time      ',      &
                                      int(tadj/(60*60)),'hours',       &
                                      mod(int(tadj/60),60),'minutes',  &
                                      mod(int(tadj),60),'seconds'  
       write(*,'(1X,A,3(X,I2,X,A))') 'Total BFS time           ',      &
                                      int(tbfs/(60*60)),'hours',       &
                                      mod(int(tbfs/60),60),'minutes',  &
                                      mod(int(tbfs),60),'seconds'  
       write(*,'(1X,A,3(X,I2,X,A))') 'Total sorting time       ',      &
                                      int(tsort/(60*60)),'hours',      &
                                      mod(int(tsort/60),60),'minutes', &
                                      mod(int(tsort),60),'seconds' 
       if ( dopim ) then
         write(*,'(1X,A,3(X,I2,X,A))') 'Total PIM analysis time  ',    &
                                        int(tpim/(60*60)),'hours',     &
                                        mod(int(tpim/60),60),'minutes',&
                                        mod(int(tpim),60),'seconds'  
       end if
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
                               maxstep,maxagg,maxdis,dopim,seed,debug)
!
       use utils
!
       implicit none
!
       include 'parameters.h'
       include 'info.h'
!
! Input/output variables
!
       character(len=leninp),intent(out)         ::  traj     !  Trajectory file name
       character(len=lenout),intent(out)         ::  outp     !  Populations file name
       character(len=leninp),intent(out)         ::  inp      !  Groups file name
       character(len=leninp),intent(out)         ::  conf     !  Structure file name
       logical                                   ::  dopim    !  PIM calculation flag
       real(kind=8),intent(out)                  ::  maxdis   !  Screening distance
       integer,intent(out)                       ::  maxagg   !  Maximum aggregate size
       integer,intent(out)                       ::  nprint   !  Populations printing steps interval
       integer,intent(out)                       ::  minstep  !  First step for analysis
       integer,intent(out)                       ::  maxstep  !  Last step for analysis
       logical,intent(out)                       ::  seed     !  Random seed flag
       logical,intent(out)                       ::  debug    !  Debug mode
!
! Local variables
!
       character(len=lencmd)                     ::  cmd     !  Command executed
       character(len=lenarg)                     ::  code    !  Executable name
       character(len=lenarg)                     ::  arg     !  Argument read
       character(len=lenarg)                     ::  next    !  Next argument to be read
       integer                                   ::  naux    !  Auxliary variable
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
       maxagg  = 10
       maxdis  = 1.5d0
       dopim   = .TRUE.
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
           case ('-d','-maxdis','--maxdis','--maximum-distance')
             call get_command_argument(i,next,status=io)
             read(next,*) maxdis
             i = i + 1
           case ('-a','-maxagg','--maxagg','--maximum-size')
             call get_command_argument(i,next,status=io)
             read(next,*) maxagg
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
           case ('-s','-seed','--seed')
             seed  = .TRUE.
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
       write(*,'(2X,A)') '-a,--maxagg           Maximum aggregate size'
       write(*,'(2X,A)') '-d,--maxdis           Cutoff distance'
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
       subroutine print_coord(xtcf,sys,outp,maxagg,nmol,nnode,imol,itag)
!
       use xdr,       only: xtcfile
       use geometry,  only: minimgvec
       use datatypes
!
       implicit none
!
       include 'info.h'
!~        include 'idxadj.h'
!
! Input/output variables
!
       type(xtcfile),intent(inout)           ::  xtcf    !  xtc file informacion   
       type(groinp),intent(in)               ::  sys     !  System information
       character(len=lenout),intent(in)      ::  outp    !  Output file name
       integer,dimension(maxagg),intent(in)  ::  nmol    !  Number of aggregates of each size
       integer,dimension(nnode),intent(in)   ::  imol    !  Molecule identifier
       integer,dimension(nnode),intent(in)   ::  itag    !  Aggregates size
       integer,intent(in)                    ::  maxagg  !  Maximum aggregate size
       integer,intent(in)                    ::  nnode   !  Number of residues
!
! Local variables
!
       type(xtcfile)                         ::  xtco    !  xtc file informacion
       character(len=lenout)                 ::  straux  !  Auxiliary string
       character(len=lenout)                 ::  aux     !  Auxiliary string
       real(kind=4),dimension(3)             ::  svaux   !  Auxiliary single precision vector       real(kind=8)                               ::  mass     !  Total mass of the aggregate
       real(kind=8),dimension(3)             ::  cofm    !  Center of mass vector
       real(kind=8)                          ::  mass    !  Total mass of the aggregate
       integer                               ::  nsize   !  Size of the previous printed aggregate
       integer                               ::  i,j,k    !  Indexes
       integer                               ::  l,m,n    !  Indexes
       integer                               ::  p,q,r,s  !  Indexes
!
! Priting coordinates of the aggregates
! -------------------------------------
!
       write(aux,*) xtcf%STEP
       aux = adjustl(aux)
!
       straux = outp(:len_trim(outp)-4)//'_'//trim(aux)//'.xtc'
       aux    = outp(:len_trim(outp)-4)//'_'//trim(aux)//'.xyz'
! Printing global coordinates in xtc format
       call xtco%init(straux,'w')
       call xtco%write(xtcf%natoms,xtcf%step,xtcf%time,xtcf%box,       &
                                                     xtcf%pos,xtcf%prec)
       call xtco%close
! Printing global coordinates in xyz format
       open(unit=uniinp,file=aux,action='write')
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
       n     = nmol(1)
       nsize = 0
!
       do i = 2, maxagg
         do j = 1, nmol(i)
           n = n + itag(n)
           if ( nsize .ne. itag(n) ) then
             write(straux,*) itag(n)
             straux = adjustl(straux)
! 
             write(aux,*) xtcf%STEP
             aux = adjustl(aux)
!
             straux = outp(:len_trim(outp)-4)//'_'//trim(aux)//        &
                                               '_'//trim(straux)//'.xyz'
             open(unit=uniinp,file=straux,action='write')
!
             nsize = itag(n)
           end if
!
           cofm(:)  = 0.0d0
           mass     = 0.0d0
           svaux(:) = xtcf%pos(:,(imol(n)-1)*sys%nat+1)
!
           m = n - 1
           do p = 1, itag(n)
             mass = mass + sys%totm
             m = m + 1
             r = (imol(m)-1)*sys%nat
             do q = 1, sys%nat
               r = r + 1 
               xtcf%pos(:,r) = minimgvec(svaux(:),xtcf%pos(:,r),       &
                                         (/ xtcf%box(1,1),             &
                                            xtcf%box(2,2),             &
                                            xtcf%box(3,3) /) )
               xtcf%pos(:,r) = svaux(:) + xtcf%pos(:,r)
               do s = 1, 3
                 cofm(s) = cofm(s) + sys%mass(q)*xtcf%pos(s,r)
               end do
             end do
           end do
           cofm(:) = cofm(:)/mass
!
           write(uniinp,*) sys%nat*itag(n)
           write(uniinp,'(20(X,I5))') imol(n:n+itag(n)-1)
!
           m = n - 1
           do p = 1, itag(n)
             m = m + 1
             r = (imol(m)-1)*sys%nat
             do q = 1, sys%nat
               r = r + 1
               write(uniinp,*) sys%atname(q),(xtcf%pos(:,r) - cofm(:))*10
             end do
           end do
!
           if ( nsize .ne. itag(n+1) ) close(uniinp)
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
! j  -> inmol
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
!  sorting the blocks by size
!
!
       subroutine blockdiag(nnode,adj,imol,iagg,itag,maxagg,nmol,nagg, &
                            tbfs,tsort,count_rate,debug)
!
       use graphtools, only: findcompundir
       use sorting,    only: ivvqsort,                                 &
                             ivqsort
!
       implicit none
!
       include 'timings.h'
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,dimension(nnode),intent(out)       ::  imol    !  Molecules identifier
       integer,dimension(nnode),intent(out)       ::  iagg    !  Aggregates identifier
       integer,dimension(nnode),intent(out)       ::  itag    !  Aggregates size
       integer,dimension(maxagg),intent(out)      ::  nmol    !  Number of aggregates of each size
       integer,intent(in)                         ::  nnode   !  Number of molecules
       integer,intent(in)                         ::  maxagg  !  Maximum aggregate size
       integer,intent(out)                        ::  nagg    !  Number of aggregates
       logical,intent(in)                         ::  debug   !  Debug mode
!
! Local variables
! 
       integer                                    ::  i,j,k   !  Indexes
!
! Block-diagonalizing the adjacency matrix
!
       call system_clock(t1bfs)     
!
       call findcompundir(nnode,adj,imol,iagg,itag,maxagg,nmol,nagg,   &
                          debug)
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
       call ivvqsort(nnode,itag,iagg,imol,1,nnode)
!
! Sorting molecules based on their aggregate identifier
!
       j = nmol(1)
       do i = 2, maxagg-1
         if ( nmol(i) .eq. 0 ) cycle
         if ( j .ge. nnode ) exit
         call ivqsort(nnode,iagg,imol,j+1,j+i*nmol(i))
         j = j + i*nmol(i)
       end do 
!
       i = j
       do while ( i .lt. nnode )
         k = itag(i+1)
         if ( (i+k+1) .le. nnode ) then
           if ( itag(i+k+1) .ne. k ) then
             call ivqsort(nnode,iagg,imol,j+1,i+k)
             j = i + k
           end if
         else
           call ivqsort(nnode,iagg,imol,j+1,nnode)
         end if
         i = i + k
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
       subroutine build_pim(ngrps,grps,nsubg,subg,natmol,igrps,thr,    &
                            nnode,imol,itag,maxagg,nmol,natconf,       &
                            coord,box,pim,debug)
!
       use geometry
!
       implicit none
!
       include 'idxadj.h'
       include 'idxpim.h'
!
! Input/output variables
!
       real(kind=8),dimension(ngrps,ngrps,maxagg-1),intent(inout)  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(natmol,natmol),intent(in)            ::  thr      !  Distance threshold
       real(kind=4),dimension(3,natconf),intent(in)                ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)                        ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(maxagg),intent(in)                        ::  nmol     !  Number of aggregates of each size
       integer,dimension(nnode),intent(in)                         ::  imol     !  Molecules identifier
       integer,dimension(nnode),intent(in)                         ::  itag     !  Aggregates size
       integer,dimension(natmol),intent(in)                        ::  grps     !  Number of subgroups in each group
       integer,dimension(natmol),intent(in)                        ::  subg     !  Number of atoms in each subgroup
       integer,dimension(natmol),intent(in)                        ::  igrps    !  
       integer,intent(in)                                          ::  ngrps    !  Number of groups
       integer,intent(in)                                          ::  nsubg    !  Number of subgroups
       integer,intent(in)                                          ::  nnode    !  Number of residues
       integer,intent(in)                                          ::  maxagg   !  Maximum aggregate size
       integer,intent(in)                                          ::  natconf  !  Total number of atoms
       integer,intent(in)                                          ::  natmol   !  Atoms per residue
       logical,intent(in)                                          ::  debug    !  Debug mode
!
! Local variables
!
       character(len=54)                                           ::  fmt1     !  Format string
       real(kind=4),dimension(3)                                   ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=8)                                                ::  dist     !  Minimum image distance
       real(kind=8)                                                ::  mindis   !  Distance threshold between groups
!
! Building pairwise interaction matrix 
! ------------------------------------
!
       fmt1 = '(8X,X,A,X,I2,2(X,A,X,I4),3(X,A,X,I2),X,A,X,I6)'
!
       iitag = nmol(1)
!
       do itype = 1, maxagg-1
         if ( nmol(itype+1) .eq. 0 ) cycle
         do inmol = 1, nmol(itype+1)
           iitag = iitag + itag(iitag)
           iimol = iitag - 1
           do irenum = 1, itag(iitag)-1
             iimol  = iimol + 1
             ingrps = (imol(iimol)-1)*natmol
!
             jimol  = iimol 
             do jrenum = irenum+1, itag(iitag)
               jimol  = jimol + 1
               jngrps = (imol(jimol)-1)*natmol
!
               iisubg = 0
               ii     = 0
               do iigrps = 1, ngrps
                 do insubg = 1, grps(iigrps)
                   iisubg = iisubg + 1 
                   ii     = ii + subg(iisubg)
                   i      = ingrps + igrps(ii)
!
                   jisubg = 0
                   jj     = 0
                   do jigrps = 1, ngrps
                     mindis = thr(iigrps,jigrps)
                     if ( mindis .gt. 1.0e-6 ) then 
                       do jnsubg = 1, grps(jigrps)
                         jisubg = jisubg + 1
                         jj     = jj + subg(jisubg)
                         j      = jngrps + igrps(jj)
!
                         r    = minimgvec(coord(:,i),coord(:,j),box)
                         dist = dot_product(r,r)
!
                         if ( dist .le. mindis ) then
                           pim(iigrps,jigrps,itype) =                  &  ! FLAG: alternatively, to save pim in triangular form
                                        pim(iigrps,jigrps,itype) + 1.0d0  !       1) find max(jigrps,iigrps) and min(jigrps,iigrps)
                           pim(jigrps,iigrps,itype) =                  &  !       2) then assign pim(max,min) 
                                        pim(jigrps,iigrps,itype) + 1.0d0       
                         end if
                       end do
                     else
                       jisubg = jisubg + grps(jigrps)
                       jj     = jj + subg(jisubg)
                     end if
                   end do            
                 end do
               end do
             end do
           end do
         end do
       end do
!
       return
       end subroutine build_pim
!
!======================================================================!
!
!~        subroutine analyze_agg(ngrps,grps,nsubg,subg,natmol,igrps,thr,  &
!~                               nnode,adj,imol,itag,maxagg,nmol,natconf, &
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
!~        integer,dimension(maxagg),intent(in)              ::  nmol     !  Number of aggregates of each size
!~        integer,dimension(nnode),intent(in)               ::  imol     !  Molecules identifier
!~        integer,dimension(nnode),intent(in)               ::  itag     !  Aggregates size
!~        integer,dimension(natmol),intent(in)              ::  grps     !  Number of subgroups in each group
!~        integer,dimension(natmol),intent(in)              ::  subg     !  Number of atoms in each subgroup
!~        integer,dimension(natmol),intent(in)              ::  igrps    !  
!~        integer,intent(in)                                ::  ngrps    !  Number of groups
!~        integer,intent(in)                                ::  nsubg    !  Number of subgroups
!~        integer,intent(in)                                ::  nnode    !  Number of residues
!~        integer,intent(in)                                ::  maxagg   !  Maximum aggregate size
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
!~        iitag = nmol(1) + nmol(2)*2 - 1
!~ !
!~        do itype = 2, maxagg-1
!~          if ( nmol(itype+1) .eq. 0 ) cycle
!~          do inmol = 1, nmol(itype+1)
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
!~              write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inmol,'of type',itype+1,'is a tree'
!~ ! Checking if the aggregate forms a linear tree
!~              if ( chkltree(itag(iitag),degree) ) then 
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inmol,'of type',itype+1,'is a linear tree'
!~                write(*,'(4X,A,X,I4)') 'The length of the longest chain is',itag(iitag)
!~              else 
!~ ! If the tree is not linear then find the longest molecular chain
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inmol,'of type',itype+1,'is a n-ary tree'
!~ !
!~                write(*,'(4X,A,X,I4)') 'The length of the longest chain is',findlongt(itag(iitag),auxadj)
!~              end if
!~ ! If the graph is not a tree then it is a cyclic graph
!~            else
!~              write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inmol,'of type',itype+1,'is a cycle'
!~ ! Checking if the aggregate forms a single cycle through all nodes of the component
!~              if ( chkscycle(itag(iitag),degree) ) then 
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inmol,'of type',itype+1,'is a simple cycle'
!~                write(*,'(4X,A,X,I4)') 'The length of the shortest simple cycle is',itag(iitag)
!~              else 
!~ ! If the aggregate does not form a single cycle then find the shortest simple cycle
!~                write(*,'(4X,2(A,X,I4,X),A)') 'Aggregate',inmol,'of type',itype+1,'is a n-cycle'
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
       subroutine analyze_agg(table,nbody,body,ngrps,grps,nsubg,subg,  &
                              natmol,igrps,thr,nnode,adj,imol,itag,    &
                              maxagg,nmol,natconf,coord,box,debug)
!
       use graphtools, only: buildadjbody
!
       implicit none
!
       include 'idxpim.h'
       include 'idxadj.h'
!
! Input/output variables
!
       real(kind=8),dimension(9,9,3),intent(inout)       ::  table    !  Number of aggregates of each type
       logical,dimension(nnode,nnode),intent(in)         ::  adj      !  Adjacency matrix in the molecule representation
       real(kind=8),dimension(natmol,natmol),intent(in)  ::  thr      !  Distance threshold
       real(kind=4),dimension(3,natconf),intent(in)      ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3),intent(in)              ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       integer,dimension(maxagg),intent(in)              ::  nmol     !  Number of aggregates of each size
       integer,dimension(nnode),intent(in)               ::  imol     !  Molecules identifier
       integer,dimension(nnode),intent(in)               ::  itag     !  Aggregates size
       integer,dimension(natmol),intent(in)              ::  body     !  Number of groups in each body
       integer,dimension(natmol),intent(in)              ::  grps     !  Number of subgroups in each group
       integer,dimension(natmol),intent(in)              ::  subg     !  Number of atoms in each subgroup
       integer,dimension(natmol),intent(in)              ::  igrps    !  Atoms identifier
       integer,intent(in)                                ::  nbody    !  Number of bodies
       integer,intent(in)                                ::  ngrps    !  Number of groups
       integer,intent(in)                                ::  nsubg    !  Number of subgroups
       integer,intent(in)                                ::  nnode    !  Number of residues
       integer,intent(in)                                ::  maxagg   !  Maximum aggregate size
       integer,intent(in)                                ::  natconf  !  Total number of atoms
       integer,intent(in)                                ::  natmol   !  Atoms per residue
       logical,intent(in)                                ::  debug    !  Debug mode
!
! Local variables
!
       logical,dimension(:,:),allocatable                ::  adjmol   !  Adjacency matrix in the molecule-based rep.  ! FLAG: rename as adjaux
       logical,dimension(:,:),allocatable                ::  adjaux   !  Auxiliary adjacency matrix                   ! FLAG: rename as adjbody
       integer,dimension(:),allocatable                  ::  degree   !  Degree of each vertex
       integer                                           ::  idxbody  !  Index for the body-body interactions
       integer                                           ::  idxoh    !  Index for the OH-OH interactions
       integer                                           ::  idxph    !  Index for the Ph-Ph interactions
!
! Analyzing aggregates of size lower or equal than 4
! --------------------------------------------------
!
       iitag = nmol(1)
!
       do itype = 2, 4
         if ( nmol(itype) .eq. 0 ) cycle
         do inmol = 1, nmol(itype)
           iitag = iitag + itag(iitag)
           iimol = iitag - 1
!
           allocate(adjmol(itag(iitag),itag(iitag)),degree(itag(iitag)))
           allocate(adjaux(itag(iitag)*nbody,itag(iitag)*nbody))
!
           adjmol(:,:) = .FALSE.
!
! Saving the adjacency matrix in the molecule-based representation
!
           i = iimol
           do irenum = 1, itag(iitag)-1
             i = i + 1
!
             jimol = i 
             do jrenum = irenum+1, itag(iitag)
               jimol = jimol + 1
!
               adjmol(irenum,jrenum) = adj(imol(i),imol(jimol))
               adjmol(jrenum,irenum) = adjmol(irenum,jrenum)
             end do
           end do
!
!~            write(*,'(2X,A)')          '-------------------------'
!~            write(*,'(2X,A,20(X,I5))') 'Saving block of molecules',(imol(i),i=iimol+1,iimol+itag(iitag))
!~            write(*,'(2X,A)')          '-------------------------'
!~            write(*,*)
!~            write(*,'(2X,2(A,X,I4,X))') 'Analyzing aggregate',inmol,'of type',itype
!~            write(*,*)
!~ !
!~            write(*,'(5X,A)') 'Adjacency matrix in the molecule-based representation'
!~            write(*,'(5X,A)') '-----------------------------------------------------'
!~            do i = 1, itag(iitag)
!~              write(*,'(5X,20L2)') (adjmol(i,j),j=1,itag(iitag))
!~            end do
!~            write(*,*)
!
! Building the adjacency matrix in the N-body simplified representation
!
           call buildadjbody(natmol,nbody,body,ngrps,grps,nsubg,       &
                             subg,igrps,adjaux,thr,itag(iitag),        &
                             imol(iimol+1:iimol+itag(iitag)),          &
                             adjmol,natconf,coord,box,debug)
!
!~            write(*,'(5X,A)') 'Adjacency matrix in the N-body simplified representation'
!~            write(*,'(5X,A)') '--------------------------------------------------------'
!~            do i = 1, itag(iitag)*nbody
!~              write(*,'(5X,20L2)') (adjaux(i,j),j=1,itag(iitag)*nbody)
!~            end do
!~            write(*,*)
!
! Analyzing the adjacency matrix for each body-body interaction
!
           do iibody = 1, nbody
!
             adjmol(:,:) = .FALSE.
!
             inbody = 0
             do irenum = 1, itag(iitag)-1
               jnbody = inbody + nbody
               do jrenum = irenum+1, itag(iitag)
                 adjmol(jrenum,irenum) = adjaux(jnbody+iibody,inbody+iibody)
                 adjmol(irenum,jrenum) = adjmol(jrenum,irenum)
!~                  write(*,*) irenum,jrenum,inbody+iibody,jnbody+iibody,adjaux(inbody+iibody,jnbody+iibody)
!
                 jnbody = jnbody + nbody
               end do
               inbody = inbody + nbody
             end do
!
!~              write(*,'(7X,A,X,I2)') 'Adjacency matrix for bodies of type',iibody
!~              write(*,'(7X,A)')      '---------------------------------------'
!~              do i = 1, itag(iitag)
!~                write(*,'(7X,20L2)') (adjmol(i,j),j=1,itag(iitag))
!~              end do
!~              write(*,*)
!
! Removing rows and columns that contain only zeros (disconnected bodies)
!
             i  = 1
             j  = itag(iitag)
!
             do irenum = 1, itag(iitag)
               degree(i) = 0
! Computing the degree of vertex i
               do jrenum = 1, j
                 if ( adjmol(jrenum,i) ) degree(i) = degree(i) + 1
               end do
! Pushing disconnected bodies to the end of the arrays
               if ( degree(i) .eq. 0 ) then
! Permuting columns of the adjacency matrix
                 adjmol(:,i) = adjmol(:,j)
                 adjmol(:,j) = .FALSE.
! Permuting rows of the adjacency matrix
                 adjmol(i,:) = adjmol(j,:)
                 adjmol(j,:) = .FALSE.
! Permuting degrees of the current node i and the last node j
                 degree(j) = 0
!
                 j = j - 1
               else
                 i = i + 1
               end if
             end do
!
!~              write(*,'(7X,A,20(X,I2))') 'Size   ',j
!~              write(*,'(7X,A,20(X,I2))') 'Degrees',degree
!~              write(*,*)
!~              write(*,'(9X,A,X,I2)') 'New adjacency matrix for bodies of type',iibody
!~              write(*,'(9X,A)')      '-------------------------------------------'
!~              do irenum = 1, itag(iitag)
!~                write(*,'(9X,20L2)') (adjmol(irenum,jrenum),            &
!~                                                    jrenum=1,itag(iitag))
!~              end do
!~              write(*,*)
!
! Classifying aggregate according to its interactions
!
             if ( iibody .eq. 1 ) then
               idxoh = idxbody(j,adjmol(:j,:j),degree(:j))
             else
               idxph = idxbody(j,adjmol(:j,:j),degree(:j))
             end if
!
           end do  
!
           table(idxoh,idxph,itype-1) = table(idxoh,idxph,itype-1)     &
                                                                 + 1.0d0
!
           deallocate(adjmol,degree)
           deallocate(adjaux)
!
         end do
       end do
!
       return
       end subroutine analyze_agg
!
!======================================================================!
!
       integer function idxbody(nnode,adj,degree)
!
       use graphtools, only: chktree,                                  &
                             chkltree,                                 &             
                             chkscycle
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj      !  Adjacency matrix in the N-body simplified rep.
       integer,dimension(nnode),intent(in)        ::  degree   !  Degree of each vertex
       integer,intent(in)                         ::  nnode    !  Number of bodies
!
! Local variables
!
       real(kind=8)                               ::  dpaux1   !
       real(kind=8)                               ::  dpaux2   !
       integer,dimension(nnode)                   ::  degaux   !  Degree of each vertex of the largest block
       integer,dimension(nnode)                   ::  ibody    !  Bodies identifier (imol)
       integer,dimension(nnode)                   ::  iblck    !  Blocks identifier (iagg)
       integer,dimension(nnode)                   ::  isize    !  Blocks size       (itag)
       integer,dimension(nnode)                   ::  nsize    !  Number of blocks of each size (nmol)
       integer                                    ::  nblck    !  Number of blocks (nagg)
       integer                                    ::  msize    !  Size of the largest block
       integer                                    ::  i        !  Index
       logical                                    ::  false    !
!
! Classifying aggregate according to its interactions
! ---------------------------------------------------
!
       false = .FALSE.
!
       if ( nnode .eq. 0 ) then
         idxbody = 9
       else if ( nnode .gt. 2 ) then
! Block-diagonalizing the adjacency matrix for bodies of type iibody
         call blockdiag(nnode,adj,ibody,iblck,isize,nnode,nsize,       &
                        nblck,dpaux1,dpaux2,1,false)
!
         msize = isize(nnode)
!
!~          write(*,'(9X,A,X,20(X,I3))') 'ibody : ',ibody(:)
!~          write(*,'(9X,A,X,20(X,I3))') 'iblck : ',iblck(:)
!~          write(*,'(9X,A,X,20(X,I3))') 'isize : ',isize(:)
!~          write(*,'(9X,A,X,20(X,I3))') 'nsize : ',nsize(:)
!~          write(*,'(9X,A,X,20(X,I3))') 'msize : ',msize
!~          write(*,'(9X,A,X,20(X,I3))') 'nblck : ',nblck
!~          write(*,*)
! Marking the aggregate according to its topology
         select case ( msize )
           case (2)
!
             idxbody = 2
!
           case (3)
! Saving the degrees of the largest block vertices
             do i = 1, msize
               degaux(i) = degree(ibody(nnode-msize+i))
             end do
!
!~              write(*,'(9X,A,20(X,I2))') 'Degrees',degaux(:msize)
!~              write(*,*)
! Checking if the aggregate is linear or cyclic
             if ( chktree(msize,degaux(:msize)) ) then 
!~                write(*,'(9X,A)') 'The bodies form a linear tree'
!~                write(*,*)
               idxbody = 3
             else
!~                write(*,'(9X,A)') 'The bodies form a simple cycle'
!~                write(*,*)
               idxbody = 4
             end if
!
           case (4)
! Saving the degrees of the largest block vertices
             do i = 1, msize
               degaux(i) = degree(ibody(nnode-msize+i))
             end do
!
!~              write(*,'(9X,A,20(X,I2))') 'Degrees',degaux(:msize)
!~              write(*,*)
! Checking if the aggregate is linear or cyclic
             if ( chktree(msize,degaux(:msize)) ) then 
               if ( chkltree(msize,degaux(:msize)) ) then 
!~                  write(*,'(9X,A)') 'The bodies form a linear tree'
!~                  write(*,*)
                 idxbody = 6
               else
!~                  write(*,'(9X,A)') 'The bodies form a tree'
!~                  write(*,*)
                 idxbody = 5
               end if
             else
               if ( chkscycle(msize,degaux(:msize)) ) then 
!~                  write(*,'(9X,A)') 'The bodies form a simple cycle'
!~                  write(*,*)
                 idxbody = 8
               else
!~                  write(*,'(9X,A)') 'The bodies form a cycle'
!~                  write(*,*)
                 idxbody = 7
               end if
             end if
!
         end select
!
       else
!
         idxbody = 1
!
       end if
!
       return
       end function idxbody
!
!======================================================================!
