!======================================================================!
!
       program aggregate

! Use the xdr interface
       use xdr, only: xtcfile,trrfile
!
       use sorting
       use datatypes
       use utils
!
       implicit none
!
       include 'parameters.h'
       include 'timings.h'
       include 'info.h'
!
       character(len=leninp)                      ::  inp      !  Trajectory file name
       character(len=leninp)                      ::  conf     !  Configuration file name
       character(len=leninp)                      ::  fgrp     !  Groups file name
       character(len=leninp)                      ::  tgrp     !  Groups file tile
       character(len=lenout)                      ::  outp     !  Populations file name
       logical,dimension(:,:),allocatable         ::  adj      !  Adjacency matrix
       logical                                    ::  dopim    !  PIM calculation flag
       real(kind=8),dimension(:,:,:),allocatable  ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(:,:),allocatable    ::  thr      !  Distance threshold
       real(kind=8),dimension(:),allocatable      ::  pop      !  Populations
       real(kind=8)                               ::  maxdis   !  Screening distance
       real(kind=8)                               ::  sumpim   !  Sum of the elements of a matrix
       integer,dimension(:),allocatable           ::  nmol     !  Number of aggregates of each size
       integer,dimension(:),allocatable           ::  imol     !  Molecules identifier
       integer,dimension(:),allocatable           ::  iagg     !  Aggregates identifier
       integer,dimension(:),allocatable           ::  itag     !  Aggregates size 
       character(len=8),dimension(:),allocatable  ::  grptag   !  Names of the groups
       integer,dimension(:),allocatable           ::  grps     !  Number of subgroups in each group
       integer,dimension(:),allocatable           ::  subg     !  Number of atoms in each subgroup
       integer,dimension(:),allocatable           ::  igrps    !  Atoms identifier
       integer                                    ::  ngrps    !  Number of groups
       integer                                    ::  nsubg    !  Number of subgroups
       integer                                    ::  nagg     !  Number of chemical species
       integer                                    ::  nnode    !  Total number of molecules
       integer                                    ::  nprint   !  Populations printing interval
       integer                                    ::  maxstep  !  Maximum step for analysis
       integer                                    ::  nsteps   !  Number of snapshots analyzed
       integer                                    ::  maxagg   !  Maximum aggregate size
       integer                                    ::  io       !  Status
       integer                                    ::  i,j,k    !  Indexes
       logical                                    ::  debug    !  Debug mode
! Declaration of a variable of type xtcfile
       type(xtcfile)                              ::  xtcf     !  xtc file informacion
! Declaration of a variable of type trrfile
       type(trrfile)                              ::  trr      !  trr file information
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
       tcpu  = 0.0d0
       tread = 0.0d0
       tadj  = 0.0d0
       tbfs  = 0.0d0
       tsort = 0.0d0
       tpim  = 0.0d0
!
! Reading command line options
!
       call command_line(inp,conf,fgrp,outp,nprint,maxstep,maxagg,     & 
                         maxdis,dopim,debug)
!
!  General settings and fatal errors check
!
       io = index(inp,'.')
       if ( io .eq. 0 )  inp = trim(inp)//'.xtc'
!
       io = index(conf,'.')
       if ( io .eq. 0 ) conf = trim(conf)//'.gro'
!
       io = index(outp,'.')
       if ( io .eq. 0 ) outp = trim(outp)//'.dat'
!
       maxdis = maxdis**2
!
! Processing Gromacs input file
!
       call read_gro(conf)
!
! Printing summary of the input information
!
       write(*,'(1X,A)') 'Input information'
       write(*,'(1X,17("-"))')
       write(*,*)
       write(*,'(2X,A,2X,A)')    'General input file name    :',       &
                                                              trim(fgrp)
       write(*,'(2X,A,2X,A)')    'Trajectory file name       :',       &
                                                               trim(inp)
       write(*,'(2X,A,2X,A)')    'Structure file name        :',       & 
                                                              trim(conf)
       write(*,'(2X,A,2X,A)')    'Populations file name      :',       &
                                                              trim(outp)
       write(*,*)
       write(*,'(2X,A,9X,F4.2)') 'Screening distance         :',maxdis  
       write(*,'(2X,A,9X,I4)')   'Maximum aggregate size     :',maxagg
       write(*,'(2X,A,2X,I11)')  'Printing interval          :',nprint
       write(*,'(2X,A,2X,I11)')  'Maximum step for analysis  :',maxstep
       if ( debug ) then
         write(*,'(2X,A,11X,A)') 'Debug mode                 :','ON'
       else
         write(*,'(2X,A,10X,A)') 'Debug mode                 :','OFF'
       end if
       write(*,*)
!
! Computing the populations of the aggregates
!
       write(*,'(4X,A)') 'Computing the populations of the aggrega'//  & 
                                              'tes along the trajectory'
       write(*,'(4X,A)') 'Please wait, this may take a while...'
       write(*,*)
! Opening populations output file
       open(unit=uniout,file=trim(outp),action='write')
!
       if ( inp(len_trim(inp)-3:) .eq. '.xtc' ) then 
! Initialize it with the names of xtc files you want to read in and write out
         call xtcf%init(trim(inp))
! Read in each configuration. Everything is stored in the xtcfile type
         call xtcf%read
!
         nnode  = xtcf%NATOMS/sys%nat
!
         allocate(thr(sys%nat,sys%nat))
         allocate(adj(nnode,nnode),imol(nnode),iagg(nnode),            &
                  itag(nnode),nmol(maxagg),pop(maxagg))
         allocate(grptag(sys%nat),grps(sys%nat),subg(sys%nat),         &
                  igrps(sys%nat))
!
! Processing general input file
!
         call read_inp(fgrp,sys%nat,tgrp,grptag,grps,subg,igrps,       &
                       ngrps,nsubg,thr)
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
         nsteps = 0
         pop(:) = 0.0d0
!
         do while ( (xtcf%STAT.eq.0) .and. (xtcf%STEP.le.maxstep) )
           if ( mod(xtcf%STEP,nprint) .eq. 0 ) then
! Update configurations counter
             nsteps = nsteps + 1
!
! Building the adjacency matrix for the current snapshot
!
             call system_clock(t1adj)     
!
             call build_adj(ngrps,grps,nsubg,subg,nnode,igrps,adj,thr, &
                            maxdis,xtcf%NATOMS,xtcf%pos,sys%nat,       &
                            sys%mass,(/xtcf%box(1,1),xtcf%box(2,2),    &
                            xtcf%box(3,3)/),debug)
!
!~              if ( debug ) then
!~                write(*,*)
!~                write(*,*) 'Molecule-based Adjacency matrix'
!~                write(*,*) '-------------------------------'
!~                write(*,*)
!~                do j = 1, nnode
!~                  write(*,*) (adj(i,j),i=1,nnode)
!~                end do
!~                write(*,*)
!~              end if
!
             call system_clock(t2adj) 
!
             tadj = tadj + dble(t2adj-t1adj)/dble(count_rate)    
!
! Block-diagonalizing the adjacency matrix
!
             call system_clock(t1bfs)     
!
             call blockdiag(nnode,adj,imol,iagg,itag,maxagg,nmol,nagg, &
                            debug)
! 
!~              if ( debug ) then
!~                write(*,*)
!~                write(*,*) 'BFS completed'
!~                write(*,*) '-------------'
!~                write(*,*)
!~                call print_info(nmol(1),nnode-nmol(1),                  &
!~                                imol(nmol(1)+1:),iagg(nmol(1)+1:),      &
!~                                itag(nmol(1)+1:),'imol-raw','iagg-raw', &
!~                                'itag-raw')
!~                write(*,*)
!~              end if
!
             call system_clock(t2bfs)     
!
             tbfs = tbfs + dble(t2bfs-t1bfs)/dble(count_rate)
! 
! Sorting molecules and aggregate identifiers based on the size of the aggregates
!
             call system_clock(t1sort)     
!
             call ivvqsort(nnode,itag,iagg,imol,1,nnode)
!
!~              if ( debug ) then
!~                write(*,*) 'Sorting based on the aggregate size com'//  &
!~                                                                 'pleted'
!~                write(*,'(45("="))') 
!~                write(*,*)
!~                call print_info(nmol(1),nnode-nmol(1),                  &
!~                                imol(nmol(1)+1:),iagg(nmol(1)+1:),      &
!~                                itag(nmol(1)+1:),'imol-size',           &
!~                                'iagg-size','itag-size')
!~                write(*,*)
!~              end if
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
!~              if ( debug ) then
!~                write(*,*) 'Sorting based on the aggregate numer co'//  &
!~                                                                'mpleted'
!~                write(*,'(46("="))')
!~                write(*,*)
!~                call print_info(nmol(1),nnode-nmol(1),                  &
!~                                imol(nmol(1)+1:),iagg(nmol(1)+1:),      &
!~                                itag(nmol(1)+1:),'imol-num','iagg-num', &
!~                                'itag-num')
!~                write(*,*)
!~              end if
!
             call system_clock(t2sort)     
!
             tsort = tsort + dble(t2sort-t1sort)/dble(count_rate)
!
! Printing population of every aggregate
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
               call print_coord(xtcf,outp,maxagg,nmol,nnode,imol,itag)
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
       else if ( inp(len_trim(inp)-3:) .eq. '.trr' ) then 
! Initialize it with the names of trr files you want to read in
         call trr%init(trim(inp))
! Read in each configuration. Everything is stored in the trrfile type 
         call trr%read
!
         nnode = xtcf%NATOMS/sys%nat
!
         allocate(adj(nnode,nnode))
!
         do while ( trr % STAT == 0 )
!
! Building the adjacency matrix for the current snapshot
!
!~            call build_adj(nnode,adj,thr,maxdis,trr%pos,trr%NATOMS,     &
!~                           sys%nat,(/trr%box(1,1),trr%box(2,2),         &
!~                           trr%box(3,3)/),debug)
!
!~            if ( debug ) then 
!~              write(*,*)
!~              write(*,*) 'Molecule-based Adjacency matrix'
!~              write(*,*) '-------------------------------'
!~              write(*,*)
!~              do j = 1, trr%NATOMS/sys%nat
!~                write(*,'(3X,13I2)') (adj(i,j),i=1,trr%NATOMS/sys%nat)
!~              end do
!~              write(*,*)
!~            end if
!
! Block-diagonalizing the adjacency matriz
!
           write(*,*) 'Not yet! Try with the .xtc file...'
           write(*,*)
           call print_end()
! Reading information of the next snapshot
           call trr%read
         end do
! Close the file
         call trr%close
       else
         write(*,*) 'Incorrect extension'
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
           if ( nmol(i+1) .eq. 0 ) cycle
           sumpim = 0.0d0
           do j = 1, ngrps
             do k = j, ngrps
               sumpim = sumpim + pim(k,j,i)
             end do
           end do
           pim(:,:,i) = pim(:,:,i)/sumpim*100
         end do
       end if
!
! Printing summary of the results
!
       pop(:) = pop(:)/nsteps
       write(*,'(1X,A)') 'Output information'
       write(*,'(1X,18("-"))')
       write(*,*)
       write(*,'(1X,A,3X,I11)')     'Number of frames analyzed : ',    &
                                                                  nsteps
       write(*,'(1X,A,20(X,F6.2))') 'Global populations        : ',    &
                                                                  pop(:)
       write(*,*)
!
       if ( dopim ) then
         do i = 1, maxagg-1
           if ( nmol(i+1) .eq. 0 ) cycle
           write(*,'(3X,A,X,I3)') 'Printing PIM for aggregates bel'//  &
                                                   'ongging to type',i+1
           write(*,'(3X,50("-"))') 
           write(*,'(8X,20(X,A8))') (adjustr(grptag(j)),j=1,ngrps)
           do j = 1, ngrps
             write(*,'(A8,20(X,F8.2))') adjustr(grptag(j)),                  &
                                                      (pim(k,j,i),k=1,j)
           end do
           write(*,*)
         end do
       end if
!
!
! Deallocate memory
!
       deallocate(thr)
       deallocate(adj,imol,iagg,itag,nmol,pop)
       deallocate(grptag,grps,subg,igrps)
       deallocate(pim)
!
       close(uniout)
!
! Printing timings
!
       call system_clock(t2)	        ! FLAG: 3 seconds lost
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
! Printing finishing date     
       call print_end()
!
       end program aggregate
!
!======================================================================!
!
       subroutine command_line(inp,conf,fgrp,outp,nprint,maxstep,      &
                               maxagg,maxdis,dopim,debug)
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
       character(len=leninp),intent(out)         ::  inp      !  Trajectory file name
       character(len=lenout),intent(out)         ::  outp     !  Populations file name
       character(len=leninp),intent(out)         ::  fgrp     !  Groups file name
       character(len=leninp),intent(out)         ::  conf     !  Structure file name
       logical                                   ::  dopim    !  PIM calculation flag
       real(kind=8),intent(out)                  ::  maxdis   !  Screening distance
       integer,intent(out)                       ::  maxagg   !  Maximum aggregate size
       integer,intent(out)                       ::  nprint   !  Populations printing steps interval
       integer,intent(out)                       ::  maxstep  !  Maximum step for analysis
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
       fgrp    = 'aggregate.inp'
       inp     = 'md.xtc'
       conf    = 'conf.gro'
       outp    = 'md.dat'
       nprint  = 1
       maxstep = 999999999
       maxagg  = 10
       maxdis  = 1.5d0
       dopim   = .TRUE.
       debug   = .FALSE.
!
! Reading command line
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
             call get_command_argument(i,fgrp,status=io)
             call check_arg(fgrp,io,arg,cmd)
             i = i + 1
           case ('-t','-traj','--traj','--trajectory') 
             call get_command_argument(i,inp,status=io)
             call check_arg(inp,io,arg,cmd)
             io = index(inp,'.')
             if ( io .eq. 0 ) then
               outp = trim(inp)//'.dat'
             else 
               outp = inp(:len_trim(inp)-4)//'.dat'
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
           case ('-m','-maxstep','--maxstep','--maximum-step')
             call get_command_argument(i,next,status=io)
             read(next,*) maxstep
             i = i + 1
           case ('-pim','--pim','--do-pim')
             dopim = .TRUE.
           case ('-nopim','--nopim','--no-pim')
             dopim = .FALSE.
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
       write(*,'(2X,A)') '-t,--trajectory       Trajectory file name'
       write(*,'(2X,A)') '-c,--configuration    Configuration file name'
       write(*,'(2X,A)') '-g,--groups-file      Groups file name'
       write(*,'(2X,A)') '-p,--populations      Populations file name'
       write(*,'(2X,A)') '-n,--nprint           Populations printi'//  &
                                                     'ng steps interval'
       write(*,'(2X,A)') '-m,--maxstep          Maximum step to be'//  &
                                                             ' analysed'
       write(*,'(2X,A)') '-a,--maxagg           Maximum aggregate size'
       write(*,'(2X,A)') '-d,--maxdis           Screening distance'
       write(*,'(2X,A)') '-v,--verbose          Debug mode'
       write(*,*)
!
       return
       end subroutine print_help
!
!======================================================================!
!
       subroutine read_gro(conf)
!
       use datatypes
       use utils
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
!
       character(len=leninp),intent(in)         ::  conf     !  Structure file name

!
! Local variables
!
       character(len=lenarg)                    ::  straux  !  Auxiliary string
       character(len=5)                         ::  aux     !
       integer                                  ::  io      !  Input/Output status
       integer                                  ::  i,j,k   !  Indexes
!
! Reading Gromacs configuration file
!
       open(unit=uniinp,file=trim(conf),action='read',                 &
            status='old',iostat=io)
       if ( io .ne. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Missing input file'
         write(*,*)
         write(*,'(3X,A)')   'Input file '//trim(conf)//               &
                                   ' not found in the current directory'
         write(*,'(2X,68("="))')
         call print_end()
       end if
!
         read(uniinp,'(A)') sys%title
         read(uniinp,*)     sys%nat
       allocate(sys%renum(sys%nat)   ,  &
                sys%rename(sys%nat)  ,  &
                sys%atname(sys%nat)  ,  &
                sys%atnum(sys%nat)   ,  &
                sys%mass(sys%nat))
!
       do k = 1, sys%nat
         read(uniinp,'(I5,2A5,I5,3F8.3)') sys%renum(k),   &
                                          sys%rename(k),  &
                                          sys%atname(k),  &
                                          sys%atnum(k)
       end do
         read(uniinp,*) sys%latvec
       close(uniinp)
!
       do j = 1, sys%nat
         aux    = adjustl(sys%atname(j))
         straux = ''
         do
           select case ( aux(1:1) )
             case ( 'a':'z','A':'Z')
               straux = trim(straux)//aux(1:1)
               aux    = aux(2:)
             case default
               exit
           end select
         end do
!
         sys%atname(j) = straux
!
         select case ( straux )
           case ( 'H' )
             sys%mass(j) = 1.007825d0
           case ( 'HE' )
             sys%mass(j) = 4.002602d0  ! Not exact
           case ( 'LI' )
             sys%mass(j) = 6.941d0     ! Not exact
           case ( 'BE' )
             sys%mass(j) = 9.012182d0  ! Not exact
           case ( 'B' )
             sys%mass(j) = 10.811d0    ! Not exact
           case ( 'C' )
             sys%mass(j) = 12.0d0
           case ( 'N' )
             sys%mass(j) = 14.003074d0
           case ( 'O' )
             sys%mass(j) = 15.994915d0
           case ( 'F' )
             sys%mass(j) = 18.998403d0 ! Not exact
           case ( 'NE' )
             sys%mass(j) = 20.1797d0   ! Not exact
           case ( 'CL' )
             sys%mass(j) = 35.453d0    ! Not exact
           case ( 'AR' )
             sys%mass(j) = 39.948d0    ! Not exact
           case ( 'KR' )
             sys%mass(j) = 83.798d0    ! Not exact
           case default
             write(*,*) straux, 'Not yet!'
             call exit(0)
         end select
       end do
!
! Computing total mass of the molecule
       sys%totm = 0.0d0
       do i = 1, sys%nat
         sys%totm = sys%totm + sys%mass(i)
       end do
!
       return
       end subroutine read_gro
!
!======================================================================!
!
       subroutine read_inp(fgrp,nat,tgrp,grptag,grps,subg,igrps,       &
                           ngrps,nsubg,thr)
!
       use utils
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
!
       character(len=leninp),intent(in)             ::  fgrp    !  Groups file name
       character(len=leninp),intent(out)            ::  tgrp    !  Groups file title
       integer,intent(in)                           ::  nat     !  Number of atoms in the molecule
       character(len=8),dimension(nat),intent(out)  ::  grptag  !  Names of the groups
       real(kind=8),dimension(nat,nat),intent(out)  ::  thr     !  Distance threshold
       integer,dimension(nat),intent(out)           ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(out)           ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(out)           ::  igrps   !  Atoms identifier
       integer,intent(out)                          ::  ngrps   !  Number of groups
       integer,intent(out)                          ::  nsubg   !  Number of subgroups
!
! Local variables
!
       character(len=leninp)                        ::  line    !
       character(len=leninp)                        ::  key     !
       integer                                      ::  io      !  Input/Output status
       integer                                      ::  i,j,k   !  Indexes
!
! Setting defaults
!
       thr(:,:) = 0.0d0
!
! Reading general input file
!
       open(unit=uniinp,file=trim(fgrp),action='read',                 &
            status='old',iostat=io)
       if ( io .ne. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Missing input file'
         write(*,*)
         write(*,'(3X,A)')   'Input file '//trim(fgrp)//               &
                                   ' not found in the current directory'
         write(*,'(2X,68("="))')
         call print_end()
       end if
! Reading input file block    
       do
! Reading input file line
         read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
         if ( io /= 0 ) exit
!~          write(*,'(A)') trim(line)  ! FLAG: dump of input data file
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Reading the different input file blocks       
         select case (key)
           case ('**GRPSFILE')
!~              write(*,*) 
!~              write(*,*) 'Reading **GRPSFILE block'
!~              write(*,*) 
             call read_grpsfile('**GRPSFILE',nat,tgrp,grptag,grps,     &
                                subg,igrps,ngrps,nsubg)      
!      
           case ('**THRESHOLD')
!~              write(*,*) 
!~              write(*,*) 'Reading **THRESHOLD block'
!~              write(*,*) 
             call read_threshold('**THRESHOLD',nat,thr,ngrps,grptag)
!
           case default
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)') 'ERROR:  Unknown block from input file'
             write(*,*) 
             write(*,'(3X,A)') 'Block '//trim(key)//' not known'
             write(*,'(2X,68("="))')
             write(*,*) 
             call print_end()
         end select  
       end do
! Closing the input file     
       close(uniinp)
!
       return
       end subroutine read_inp
!
!======================================================================!
!
       subroutine read_grpsfile(blck,nat,tgrp,grptag,grps,subg,igrps,  &
                                ngrps,nsubg)
!
       use utils
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
!
       character(len=*),intent(in)                  ::  blck    !  Block name
       character(len=leninp),intent(out)            ::  tgrp    !  Groups file title
       integer,intent(in)                           ::  nat     !  Number of atoms in the molecule
       character(len=8),dimension(nat),intent(out)  ::  grptag  !  Names of the groups
       integer,dimension(nat),intent(out)           ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(out)           ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(out)           ::  igrps   !  Atoms identifier
       integer,intent(out)                          ::  ngrps   !  Number of groups
       integer,intent(out)                          ::  nsubg   !  Number of subgroups
!
! Local variables
!
       character(len=leninp)                        ::  line    !
       character(len=leninp)                        ::  key     !
       integer                                      ::  posi    !
       integer                                      ::  aux     !
       integer                                      ::  io      !  Input/Output status
       integer                                      ::  i,j,k   !  Indexes
!
! Reading GRPSFILE block options 
!   
       tgrp = 'GRPSFILE Title'
!
       grps(:)  = 0
       subg(:)  = 0
       igrps(:) = 0
!
       ngrps = 0
       nsubg = 0
!
       aux = 0
!
       do
         read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
         if ( io /= 0 ) call endblck(blck)
!~        write(*,'(A)') trim(line) ! FLAG: dump of input data file
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.TITLE')
!~              write(*,*) 
!~              write(*,*) 'Reading .TITLE option'
!~              write(*,*)
             read(uniinp,*) tgrp
             tgrp = adjustl(tgrp)
!
           case ('*GRP')
!~              write(*,*) 
!~              write(*,*) 'Reading *GRP section'
!~              write(*,*)
             call read_grp('*GRP',nat,grptag,grps,subg,igrps,ngrps,    &
                           nsubg,aux)
!
           case ('*SGRP')
!~              write(*,*) 
!~              write(*,*) 'Reading *SGRP section'
!~              write(*,*)
             call read_grp('*SGRP',nat,grptag,grps,subg,igrps,ngrps,   &
                           nsubg,aux)
!
           case ('**END')
!~              write(*,*) 
!~              write(*,*) 'Exiting from **GRPSFILE block'
!~              write(*,*)
             return
           case default
             call unksect(key,blck)
         end select  
       end do
!
       return
       end subroutine read_grpsfile
!
!======================================================================!
!
       subroutine read_grp(sect,nat,grptag,grps,subg,igrps,ngrps,nsubg,&
                           aux)
!
       use utils
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
!
       character(len=*),intent(in)                    ::  sect    !  Section name
       integer,intent(in)                             ::  nat     !  Number of atoms in the molecule
       character(len=8),dimension(nat),intent(inout)  ::  grptag  !  Names of the groups
       integer,dimension(nat),intent(inout)           ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(inout)           ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(inout)           ::  igrps   !  Atoms identifier
       integer,intent(inout)                          ::  ngrps   !  Number of groups
       integer,intent(inout)                          ::  nsubg   !  NUmber of subgroups
       integer,intent(inout)                          ::  aux     !
!
! Local variables
!
       character(len=leninp)                          ::  line    ! 
       character(len=leninp)                          ::  key     ! 
       character(len=leninp)                          ::  arg     !  
       integer                                        ::  natgrp  !
       integer                                        ::  posi    ! 
       integer                                        ::  io      !  Input/Output status
       integer                                        ::  old     !
       integer                                        ::  i,j,k   !  Indexes
!
! Reading GRP section keywords 
!
       do
         read(uniinp,'(A)',iostat=io) line  
!
         if ( io /= 0 ) call endsect(sect)
!~          write(*,'(A)') trim(line)   ! FLAG: dump of input data file
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
!
! Procesing the keywords line
!
         ngrps = ngrps + 1  
!
         do
           if ( len_trim(line) == 0 ) return
! Saving the keyword 
           posi = scan(key,'=') 
           if ( posi .ne. 0 ) then 
             line = key(posi+1:)  
             line = adjustl(line)
!
             key  = key(:posi-1)
             key  = lowercase(key)
           else
             call errkey('section',sect)
           end if
! Saving the arguments
           select case (key)
             case ('atoms')
               call chkkeyarg(key,line,arg)
!
               read(arg,*) natgrp
!
               if ( trim(sect) .eq. '*GRP' ) then
                 nsubg       = nsubg + 1
                 grps(ngrps) = 1
                 subg(nsubg) = natgrp
               else if ( trim(sect) .eq. '*SGRP' ) then
                 old   = nsubg
                 nsubg = nsubg + natgrp
                 grps(ngrps) = natgrp
                 subg(old+1:nsubg) = 1
               end if    
!
               read(uniinp,*) igrps(aux+1:aux+natgrp)
!
               aux = aux + natgrp
!        
             case ('name')
               call chkkeyarg(key,line,arg)
               grptag(ngrps) = arg
!
             case default
               call unkkey(key,sect)
           end select  
!
           key = line
!
         end do
       end do
!
       return
       end subroutine read_grp
!
!======================================================================!
!
       subroutine read_threshold(blck,nat,thr,ngrps,grptag)
!
       use utils
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
!
       real(kind=8),dimension(nat,nat),intent(inout)  ::  thr     !
       character(len=8),dimension(nat),intent(in)     ::  grptag  !
       character(len=*),intent(in)                    ::  blck    !
       integer,intent(in)                             ::  nat     !
       integer,intent(in)                             ::  ngrps   !  Number of groups

!
! Local variables
!
       character(len=leninp)                          ::  line    !
       character(len=leninp)                          ::  key     !
       character(len=8)                               ::  caux1   !
       character(len=8)                               ::  caux2   !
       real(kind=8)                                   ::  daux    !
       integer                                        ::  iaux1   !
       integer                                        ::  iaux2   !
       integer                                        ::  posi    !
       integer                                        ::  aux     !
       integer                                        ::  io      !  Input/Output status
       integer                                        ::  i,j,k   !  Indexes
!
! Reading THRESHOLD block options 
!
       do
         read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
         if ( io /= 0 ) call endblck(blck)
!~        write(*,'(A)') trim(line) ! FLAG: dump of input data file
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.THR')
!~              write(*,*) 
!~              write(*,*) 'Reading .THR option'
!~              write(*,*)
             read(uniinp,*) daux
             thr(:,:) = daux
!
           case ('.VALUES')
!~              write(*,*) 
!~              write(*,*) 'Reading .VALUES option'
!~              write(*,*)
             do
               read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
               if ( io /= 0 ) call endkey('.VALUES')
! If reads a white line then reads the next line
               if ( len_trim(line) == 0 ) cycle
! Processing the line read 
               call chkcomment(line,key)
! If the line just contains a comment reads the next line
               if ( len_trim(key) == 0 ) cycle
! If the block has finished exit
               if ( uppercase(key(:5)) .eq. '**END' ) return  ! FLAG: if .THR specified after .VALUES code breaks 
!
               posi  = scan(key,' ') 
!
               caux1 = key(:posi-1)
               iaux1 = findcv(ngrps,grptag,caux1)  
               if ( iaux1 .eq. 0 ) call errkey('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
!
               caux2 = key(:posi-1)
               iaux2 = findcv(ngrps,grptag,caux2)  
               if ( iaux1 .eq. 0 ) call errkey('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
!
               key = key(:posi-1)
               read(key,*) daux
!
               thr(iaux1,iaux2) = daux
               thr(iaux2,iaux1) = daux
             end do
!
           case ('**END')
!~              write(*,*) 
!~              write(*,*) 'Exiting from **THRESHOLD block'
!~              write(*,*)
             return
!
           case default
             call unksect(key,blck)
         end select  
       end do
!
       return
       end subroutine read_threshold
!
!======================================================================!
!
       subroutine build_adj(ngrps,grps,nsubg,subg,nnode,igrps,adj,thr, &
                            maxdis,natconf,coord,natmol,mass,box,debug)
!
       use geometry
!
       implicit none
!
       include 'idxadj.h'
!
! Input/output variables
!
       real(kind=4),dimension(3,natconf),intent(inout)   ::  coord    !  Atomic coordinates !FLAG: kind=8 to kind=4
       real(kind=8),dimension(natmol,natmol),intent(in)  ::  thr      !  Distance threshold
       real(kind=8),dimension(natmol),intent(in)         ::  mass     !  Atomic masses
       real(kind=4),dimension(3),intent(in)              ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=8),intent(in)                           ::  maxdis   !  Screening distance
       logical,dimension(nnode,nnode),intent(out)        ::  adj      !  Adjacency matrix
       integer,dimension(natmol),intent(in)              ::  grps     !  Number of subgroups in each group
       integer,dimension(natmol),intent(in)              ::  subg     !  Number of atoms in each subgroup
       integer,dimension(natmol),intent(in)              ::  igrps    !
       integer,intent(in)                                ::  ngrps    !  Number of groups
       integer,intent(in)                                ::  nsubg    !  Number of subgroups
       integer,intent(in)                                ::  nnode    !  Number of residues
       integer,intent(in)                                ::  natconf  !  Total number of atoms
       integer,intent(in)                                ::  natmol   !  Atoms per residue
       logical,intent(in)                                ::  debug    !  Debug mode
!
! Local variables
!
       character(len=54)                                 ::  fmt1     !  Format string
!~        logical,dimension(nnode)                          ::  check    !
       real(kind=4),dimension(3,natmol)                  ::  atcoord  !  Subgroup coordinates !FLAG: kind=8 to kind=4
       real(kind=8),dimension(natmol)                    ::  atmass   !  Subgroup masses
       real(kind=8),dimension(3)                         ::  cofm     !  Center of mass coordinates !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                         ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3)                         ::  svaux    !  Auxiliary single precision vector
       real(kind=8)                                      ::  mindis   !  Distance threshold between groups
       real(kind=8)                                      ::  dist     !  Minimum image distance
!
! Saving coordinates based on the n-body simplified representation
!
!~ write(*,'(X,A,20(X,I4))') 'grps   :',grps
!~ write(*,'(X,A,20(X,I4))') 'subg   :',subg
!~ write(*,'(X,A,20(X,I4))') 'isubg  :',igrps
!
       do irenum = 1, nnode
!~ write(*,*)
!~ write(*,'(X,A,X,I4)') 'Starting residue',irenum
!~ write(*,'(X,21("-"))')
         ingrps = (irenum-1)*natmol
         iisubg = 0
         ii     = 0
         do iigrps = 1, ngrps
!~ write(*,*)
!~ write(*,'(2X,2(X,A,X,I4))') 'Printing atoms belonging to group',iigrps,'of residue',irenum
!~ write(*,'(3X,A,X,I4)') 'Number of subgroups',grps(iigrps)
           svaux(:) = coord(:,ingrps+igrps(ii+1))
           do insubg = 1, grps(iigrps)
             iisubg = iisubg + 1
!~ write(*,*)
!~ write(*,'(4X,2(X,A,X,I4))') 'Number of atoms in subgroup',iisubg,'is',subg(iisubg)
!~ write(*,*)
             if ( subg(iisubg) .ne. 1 ) then
               i = ii
               do iat = 1, subg(iisubg)
                 i = i + 1
!
                 atcoord(:,iat) = minimgvec(svaux(:),                  &
                                           coord(:,ingrps+igrps(i)),   &
                                           box)
                 atcoord(:,iat) = svaux(:) + atcoord(:,iat)
!
                 atmass(iat)    = mass(igrps(i))
               end do
!
               cofm(:) = cofm_vector(subg(iisubg),                     &
                                     atcoord(:,:subg(iisubg)),         &
                                     atmass(:subg(iisubg)) )
! 
               i = ii
               do iat = 1, subg(iisubg)
                 i = i + 1
                 coord(:,ingrps+igrps(i)) = cofm(:)
               end do            
             end if
!
             ii = ii + subg(iisubg)
!
!~ write(*,'(6X,3(X,A,X,I4),X,A,X,I8)') 'group',iigrps,'subgroup',iisubg, &
!~                                          'ii',ii,'atom',ingrps+igrps(ii)
           end do
         end do
       end do
!~ write(*,*)
!
! Building adjacency matrix
!
       adj(:,:) = .FALSE. 
!
       do irenum = 1, nnode-1
         ingrps = (irenum-1)*natmol
!
!~          if ( debug ) then
!~            write(*,*)
!~            write(*,'(3X,3(A,X,I5,X))') 'RESIDUE',irenum,'atom',        &
!~                                                      (irenum-1)*natmol+1
!~            write(*,'(3X,22("-"))')
!~          end if
!
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
                 mindis = thr(iigrps,jigrps)   ! FLAG: if mindis eq 0 cycle
                 do jnsubg = 1, grps(jigrps)
                   jisubg = jisubg + 1
                   jj     = jj + subg(jisubg)
                   j      = jngrps+igrps(jj)
!
                   r    = minimgvec(coord(:,i),coord(:,j),box)
                   dist = dot_product(r,r)
!
!~                    if ( debug ) then
!~                      fmt1 = '(2(X,A,X,I4,3(X,A,X,I2),X,A,X,I5),X,A'//  &
!~                                                                 ',F6.3)'
!~                      write(*,fmt1) 'COMPARING irenum',irenum,          &
!~                                              'iigrps',iigrps,          &
!~                                              'iisubg',iisubg,          &
!~                                              'iat',ii,'i',i,           &
!~                                         'WITH jrenum',jrenum,          &
!~                                              'jigrps',jigrps,          &
!~                                              'jisubg',jisubg,          &
!~                                              'jat',jj,'j',j,           &
!~                                         'DIST',dist
!~                    end if
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
!
                 end do
               end do
             end do
           end do
1000       continue           
         end do
       end do
!
       if ( debug ) write(*,*)
!
! OLD ALGORITHM v2
! ----------------
!
!~        adj(:,:) = .FALSE. 
!~ !
!~        do irenum = 1, nnode-1
!~ !!~            if ( debug ) then
!~ !!~              write(*,*)
!~ !!~              write(*,'(3X,3(A,X,I4,X))') 'RESIDUE',irenum,'atom',      &
!~ !!~                                                      (irenum-1)*natmol+1
!~ !!~              write(*,'(3X,22("-"))')
!~ !!~            end if
!~          do jrenum = irenum+1, nnode
!~            iat  = 0
!~            i    = (irenum-1)*natmol
!~            dist = thr + 1.0e-10
!~            do while ( (iat.lt.natmol) .and. (dist.ge.thr) )
!~              iat = iat + 1
!~              i   = i + 1
!~              jat = 0
!~              j   = (jrenum-1)*natmol 
!~              do while ( (jat.lt.natmol) .and. (dist.le.maxdis) .and.   &
!~                                                          (dist.ge.thr) ) 
!~                jat = jat + 1  
!~                j   = j + 1
!~ !
!~                r    = minimgvec(coord(:,i),coord(:,j),box)
!~                dist = sqrt(dot_product(r,r))
!~ !
!~ !!~                if ( debug ) then
!~ !!~                  write(*,'(6(X,A,X,I5),X,A,F6.3)') 'Comparing irenum', &
!~ !!~                                                    irenum,'iat',iat,   &
!~ !!~                                                    'i',i,'with jrenum',&
!~ !!~                                                    jrenum,'jat',jat,   &
!~ !!~                                                    'j',j,'dist',dist
!~ !!~                end if
!~ !
!~                if ( dist .le. thr ) then
!~                  adj(jrenum,irenum) = .TRUE.
!~                  adj(irenum,jrenum) = .TRUE.
!~                end if
!~              end do
!~            end do
!~          end do
!~        end do
!
!~        if ( debug ) write(*,*)
!
! OLD ALGORITHM v1
! ----------------
!
!~        adj(:,:) = .FALSE.
!~        dist     = thr + 0.1
!~        i = 1
!~        do irenum = 1, nnode-1
!~          check(:) = .TRUE.
!~          do iat = 1, natmol
!~            if ( debug ) then
!~              write(*,*)
!~              write(*,'(3X,3(A,X,I4,X))') 'RESIDUE',irenum,'atom',iat,  &
!~                                                                    'i',i
!~              write(*,'(3X,29("-"))')
!~            end if
!~            jj = irenum*natmol + 1
!~            do jrenum = irenum + 1, nnode
!~              if ( (.not. adj(jrenum,irenum)) .and. check(jrenum) ) then
!~                j    = jj
!~                jat  = 1
!~                dist = thr + 1.0e-10
!~                do while ( (jat.le.natmol) .and. (dist.le.maxdis) .and. &
!~                           (dist.ge.thr) ) 
!~                  r    = minimgvec(coord(:,i),coord(:,j),box)
!~                  dist = sqrt(dot_product(r,r))
!~ !
!~                  if ( debug ) then
!~                    write(*,'(1X,4(A,X,I4,X),A,F6.3)') 'Comparing w'//  &
!~                          'ith residue',jrenum,'atom',jat,'j',j,'jj',   &
!~                                                           jj,'dist',dist
!~                  end if
!~ !
!~                  if ( dist .le. thr ) then
!~                    adj(jrenum,irenum) = .TRUE.
!~                    adj(irenum,jrenum) = .TRUE.
!~                  end if
!~                  j   = j + 1
!~                  jat = jat + 1
!~                end do
!~                if ( dist .gt. maxdis ) check(jrenum) = .FALSE.
!~              end if
!~              jj = jj + natmol
!~            end do
!~            i = i + 1
!~          end do
!~        end do
!
       return
       end subroutine build_adj
!
!======================================================================!
!
       subroutine blockdiag(nnode,adj,imol,iagg,itag,maxagg,nmol,nagg, &
                            debug)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(inout)  ::  adj     !  Adjacency matrix
       integer,dimension(nnode),intent(out)          ::  imol    !  Molecules identifier
       integer,dimension(nnode),intent(out)          ::  iagg    !  Aggregates identifier
       integer,dimension(nnode),intent(out)          ::  itag    !  Aggregates size
       integer,dimension(maxagg),intent(out)         ::  nmol    !  Number of aggregates of each size
       integer,intent(out)                           ::  nagg    !  Number of chemical species
       integer,intent(in)                            ::  nnode   !  Number of residues
       integer,intent(in)                            ::  maxagg  !  Maximum aggregate size
       logical,intent(in)                            ::  debug   !  Debug mode
!
! Local variables
!
       logical,allocatable,dimension(:)              ::  notvis  !  Nodes visited
       integer,allocatable,dimension(:)              ::  queue   !  Queue of connected nodes
       integer                                       ::  iqueue  !  Queue index
       integer                                       ::  jqueue  !  Queue index
       integer                                       ::  inode   !  Node index
       integer                                       ::  jnode   !  Node index
       integer                                       ::  knode   !  Node index
       integer                                       ::  nqueue  !  Queue elements
       integer                                       ::  inmol   !  Number of aggregates index
       integer                                       ::  ntag    !  Size of the aggregate
       integer                                       ::  intag   !  Size of the aggregate index
       integer                                       ::  i,j,k   !  Indexes
!
! Performing Breadth First Search over the target molecules
!
       allocate(notvis(nnode),queue(nnode))
! Mark all the vertices as not visited
       notvis  = .TRUE.
! Initializing the molecules information
       imol(:) = 0   
       nmol(:) = 0
       inmol   = 1
! Initializing the aggregates information
       nagg    = 0
       iagg(:) = 0
! Initializing the size information
       ntag    = 0
       intag   = 0
       itag(:) = 0
!
! Outer loop over each node
!
       do inode = 1, nnode
         if ( notvis(inode) ) then
! Marking head node as visited
           notvis(inode) = .FALSE.
! Updating the system information
           nagg        = nagg + 1
           iagg(inmol) = nagg
!
           imol(inmol) = inode
           inmol       = inmol + 1
!
           ntag        = 1
! Initializing queue
           queue = 0
! Adding current node to the queue
           queue(1) = inode
! Initializing the queue counter
           iqueue = 1
! Setting the next element in the queue
           nqueue = 2
!~            if ( debug ) then
!~              write(*,*)
!~              write(*,'(X,A,I4)') 'Starting outer loop cycle in node',  &
!~                                                                    inode
!~              write(*,*)          '====================================='
!~              write(*,*)
!~              write(*,'(X,A,20L3)') 'Visited               :',notvis
!~              write(*,'(X,A,20I3)') 'Queue                 :',queue
!~              write(*,*) 'New aggregate found'
!~              write(*,'(X,A,I3)')   'Number of aggregates  :',nagg
!~              write(*,'(X,A,20I3)') 'Identifier            :',iagg
!~              write(*,*)
!~            end if
!
! Inner loop over the queue elements
!
           do while ( iqueue .lt. nqueue )
! Saving actual element in the queue
             knode = queue(iqueue)
! Checking if node k was already visited
!~                if ( debug ) then
!~                  write(*,*)
!~                  write(*,'(X,A,I4)') 'Taking from the queue node',knode
!~                  write(*,*)          '-------------------------------'
!~                  write(*,*)
!~                  write(*,'(3X,A,I4)') 'Starting inner loop'
!~                  write(*,'(3X,A)')    '-------------------'
!~                  write(*,*)
!~                end if
! Check the connection between actual queue element and the rest of nodes
               do jnode = inode + 1, nnode
!~                  if ( debug ) write(*,'(4X,A,I4,A)') 'Checking if '//  &
!~                                              'node',jnode,' was visited'
! Checking if node j is connected to node k or has been already visited
                 if ( adj(jnode,knode) .and. notvis(jnode) ) then
! Updating the system information
                   iagg(inmol)   = nagg
!
                   imol(inmol)   = jnode
                   inmol         = inmol + 1
!
                   ntag          = ntag + 1
! Marking the node connected to node k as visited
                   notvis(jnode) = .FALSE.
! Adding to the queue the node connected to node k
                   queue(nqueue) = jnode
! Updating next element in the queue
                   nqueue        = nqueue + 1
!~                    if ( debug ) then
!~                      write(*,*)
!~                      write(*,'(4X,A)')      'New connection found'
!~                      write(*,'(4X,A,I4)')   'Adding to the queue e'//  &
!~                                                           'lement',jnode
!~                      write(*,*)
!~                      write(*,'(4X,A,20I3)') 'Queue                '//  &
!~                                                               ' :',queue
!~                      write(*,'(4X,A,I3)')   'Number of aggregates '//  & 
!~                                                                ' :',nagg
!~                      write(*,'(4X,A,20I3)') 'Identifier           '//  &
!~                                                                ' :',iagg
!~                      write(*,*)
!~                    end if
                 end if
               end do
! Updating the queue counter
             iqueue = iqueue + 1
           end do
! Saving the size of the aggregate found
           do i = intag+1, intag+ntag
             itag(i) = ntag
           end do
           intag = intag + ntag
! Update the number of aggregates of each size
           if ( ntag .lt. maxagg ) then
             nmol(ntag)   = nmol(ntag)   + 1
           else
             nmol(maxagg) = nmol(maxagg) + 1
           end if
         end if
       end do
!
       deallocate(notvis,queue)
!
       return
       end subroutine blockdiag
!
!======================================================================!
!
       subroutine print_coord(xtcf,outp,maxagg,nmol,nnode,imol,itag)
!
       use xdr,       only: xtcfile
       use geometry,  only: minimgvec
       use datatypes
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
!
       type(xtcfile),intent(inout)           ::  xtcf    !  xtc file informacion   
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
       k = 0                  !  FLAG: remove the line
       n     = nmol(1)
       nsize = 0
!
       do i = 2, maxagg
!~ write(*,*)
!~ write(*,'(X,A,I4)') 'Molecules belongging to type',i
!~ write(*,'(X,A)')    '--------------------------------'
!
         k = k + nmol(i-1)    !  FLAG: remove the line
         do j = 1, nmol(i)
           l = k + j          !  FLAG: remove the line
           n = n + itag(n)
!  
!~ write(*,'(3X,2(A,I6))') 'Aggregate number',l, ' starts in molecule',n
!~ write(*,'(3x,2(A,I6))') 'Coordinates of molecule',imol(n),' start in atom',(imol(n)-1)*sys%nat+1
!
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
!~ write(*,*) 'Opening file',straux
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
!~ write(*,'(3X,2(A,I6))') 'Molecule index',imol(m),' for molecule',m
             do q = 1, sys%nat
               r = r + 1
!~ write(*,'(3X,A,I6)') 'Atom index    ',r
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
! Seen on: Gaines & Di Tommaso. Pharmaceutics, 2018, 10. 
!          10.3390/pharmaceutics10010012
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
       real(kind=8)                                                ::  mindis   !  Distance threshold between groups
       real(kind=8)                                                ::  dist     !  Minimum image distance
       integer                                                     ::  itype    !  Aggregates type
!
! Building pairwise interaction matrix 
!
       iitag = nmol(1)
!
       do itype = 1, maxagg-1
         if ( nmol(itype+1) .eq. 0 ) cycle
!
!~          write(*,'(3X,A,X,I4)') 'Building PIM for aggregates belon'//  &
!~                                                  'gging to type',itype+1
!~          write(*,'(3X,51("-"))') 
!~          write(*,*)
!
         do inmol = 1, nmol(itype+1)
           iitag = iitag + itag(iitag)
!  
!~ write(*,'(5X,A,X,I6)')    'NEW AGGREGATE starts in molecule',iitag
!~ write(*,*)
!
           iimol = iitag - 1
           do irenum = 1, itag(iitag)-1
             iimol  = iimol + 1
             ingrps = (imol(iimol)-1)*natmol
!
!~ write(*,'(6X,5(X,A,X,I4),X,A,X,I6)') 'TYPE',itype+1,'AGGREGATE',inmol, &
!~      'irenum',irenum,'iimol',iimol,'imol',imol(iimol),'ingrps',ingrps+1
!~ write(*,*)
!
             jimol = iimol 
             do jrenum = irenum+1, itag(iitag)
               jimol  = jimol + 1
               jngrps = (imol(jimol)-1)*natmol
!
               iisubg = 0
               ii     = 0
               do iigrps = 1, ngrps
                 do insubg = 1, grps(iigrps)
                   iisubg = iisubg + 1 
                   i      = ii
                   do iat = 1, subg(iisubg)
                     i = i + 1
!
                     fmt1 = '(8X,X,A,X,I2,2(X,A,X,I4),3(X,A,X,I2),'//  &
                                                             'X,A,X,I6)'
!~                      write(*,fmt1) 'COMPARING '//                      &
!~                                    'irenum',irenum,                    &
!~                                     'iimol',iimol,                     &
!~                                      'imol',imol(iimol),               &
!~                                     'igrps',iigrps,                    &
!~                                     'isubg',iisubg,                    &
!~                                       'iat',i,                         &
!~                                         'i',ingrps+igrps(i)
!~                      write(*,*)
!
                     jisubg = 0
                     jj     = 0
                     do jigrps = 1, ngrps
                       mindis = thr(iigrps,jigrps)   ! FLAG: if mindis eq 0 cycle
                       do jnsubg = 1, grps(jigrps)
                         jisubg = jisubg + 1
                         j      = jj
                         do jat = 1, subg(jisubg)
                           j = j + 1
!
!~                            write(*,fmt1) '     WITH '//                &
!~                                          'jrenum',jrenum,              &
!~                                           'jimol',jimol,               &
!~                                            'jmol',imol(jimol),         &
!~                                           'jgrps',jigrps,              &
!~                                           'jsubg',jisubg,              &
!~                                             'jat',j,                   &
!~                                               'j',jngrps+igrps(j)
!
                           r    = minimgvec(coord(:,ingrps+igrps(i)),  &
                                            coord(:,jngrps+igrps(j)),  &
                                            box)
                           dist = dot_product(r,r)
!
                           if ( dist .le. mindis ) then
                             pim(iigrps,jigrps,itype) =                &  ! FLAG: alternatively, to save pim in triangular form
                                            pim(iigrps,jigrps,itype) + 1  !       1) find max(jigrps,iigrps) and min(jigrps,iigrps)
                             pim(jigrps,iigrps,itype) =                &  !       2) then assign pim(max,min) 
                                            pim(jigrps,iigrps,itype) + 1       
                           end if
                         end do
                         jj = jj + subg(jisubg)
!~ !                     j  = jngrps+igrps(jj)   ! FLAG: value of j if we remove the loop
                       end do
                     end do
!~ write(*,*)
                   end do            
                   ii = ii + subg(iisubg)
!~ !               i  = ingrps+igrps(ii)   ! FLAG: value of i if we remove the loop
                 end do
               end do
!~ write(*,*)
             end do
           end do
         end do
       end do
!


!~        do irenum = 1, nnode-1
!~          ingrps = (irenum-1)*natmol
!~ !
!         if ( debug ) then
!           write(*,*)
!           write(*,'(3X,3(A,X,I4,X))') 'RESIDUE',irenum,'atom',        &
!                                                     (irenum-1)*natmol+1
!           write(*,'(3X,22("-"))')
!         end if
!~ !
!~          do jrenum = irenum+1, nnode
!~            jngrps = (jrenum-1)*natmol
!~ !
!~            iisubg = 0
!~            ii     = 0
!~            do iigrps = 1, ngrps
!~              do insubg = 1, grps(iigrps)
!~                iisubg = iisubg + 1
!~                ii     = ii + subg(iisubg)
!~                i      = ingrps+igrps(ii)
!~ !
!~                jisubg = 0
!~                jj     = 0
!~                do jigrps = 1, ngrps
!~                  mindis = thr(iigrps,jigrps)
!~                  do jnsubg = 1, grps(jigrps)
!~                    jisubg = jisubg + 1
!~                    jj     = jj + subg(jisubg)
!~                    j      = jngrps+igrps(jj)
!~ !
!~                    r    = minimgvec(coord(:,i),coord(:,j),box)
!~                    dist = sqrt(dot_product(r,r))
!~ !
!                   if ( debug ) then
!                     fmt1 = '(2(X,A,X,I4,3(X,A,X,I2),X,A,X,I4),X,A'//  &
!                                                                ',F6.3)'
!                     write(*,fmt1) 'COMPARING irenum',irenum,          &
!                                             'iigrps',iigrps,          &
!                                             'iisubg',iisubg,          &
!                                             'iat',ii,'i',i,           &
!                                        'WITH jrenum',jrenum,          &
!                                             'jigrps',jigrps,          &
!                                             'jisubg',jisubg,          &
!                                             'jat',jj,'j',j,           &
!                                        'DIST',dist
!                   end if
!~ !
!~                    if ( dist .le. mindis ) then
!~ write(*,*) 'Connection found'
!~                    end if
!~ !
!~                  end do
!~                end do
!~              end do
!~            end do
!~ 1000       continue           
!~          end do
!~        end do
!~ !
!~        if ( debug ) write(*,*)
     
!
!       iitag = 0
!       do i = 1, itype-1
!         iitag = iitag + nmol(i)
!write(*,*) iitag
!       end do
!
!       do inmol = 1, nmol(itype)
!         iitag = iitag + itag(iitag)
!  
!write(*,*)
!write(*,'(5X,A,X,I6)') 'New aggregate starts in molecule',imol(iitag)
!
!         iimol = iitag - 1
!         do irenum = 1, itag(iitag)
!           iimol = iimol + 1
!
!~          if ( debug ) then
!~            write(*,*)
!~            write(*,'(3X,3(A,X,I4,X))') 'RESIDUE',irenum,'atom',        &
!~                                                      (irenum-1)*natmol+1
!~            write(*,'(3X,22("-"))')
!~          end if
!
!
!write(*,*)
!write(*,'(6X,4(X,A,X,I6))') 'RESIDUE',irenum,'molecule',iimol,'index',imol(iimol),'atom',(imol(iimol)-1)*natmol+1
!write(*,*)
!~            ingrps = (imol(iimol)-1)*natmol
!~            iisubg = 0
!~            ii     = 0
!~            do iigrps = 1, ngrps
!~ write(*,*)
!~ write(*,'(8X,2(X,A,X,I4))') 'Printing atoms belonging to group',iigrps,'of residue',irenum
!~ write(*,'(9X,A,X,I4)') 'Number of subgroups',grps(iigrps)
!~              do insubg = 1, grps(iigrps)
!~                iisubg = iisubg + 1
!~ !
!~ write(*,*)
!~ write(*,'(11X,2(X,A,X,I4))') 'Number of atoms in subgroup',iisubg,'is',subg(iisubg)
!~ write(*,*)
!~ ! 
!~                i = ii
!~                do iat = 1, subg(iisubg)
!~                  i = i + 1
!~ write(*,'(13X,A,I6)') 'Atom index    ',ingrps+igrps(i)
!~ !
!write(*,'(6X,3(X,A,X,I4),X,A,X,I8)') 'group',iigrps,'subgroup',iisubg, &
!                                         'ii',ii,'atom',ingrps+igrps(i)
!~ !
!                   if ( debug ) then
!                     fmt1 = '(2(X,A,X,I4,3(X,A,X,I2),X,A,X,I4),X,A'//  &
!                                                                ',F6.3)'
!                     write(*,fmt1) 'COMPARING irenum',irenum,          &
!                                             'iigrps',iigrps,          &
!                                             'iisubg',iisubg,          &
!                                             'iat',ii,'i',i,           &
!                                        'WITH jrenum',jrenum,          &
!                                             'jigrps',jigrps,          &
!                                             'jisubg',jisubg,          &
!                                             'jat',jj,'j',j,           &
!                                        'DIST',dist
!                   end if
!~                end do            
!~ !
!~                ii = ii + subg(iisubg)
!               i      = ingrps+igrps(ii)   ! FLAG: value of i if we remove the loop
!~ !
!~              end do
!~            end do
!         end do
!       end do
!
       return
       end subroutine build_pim
!
!======================================================================!
