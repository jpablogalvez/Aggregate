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
       include 'info.h'
!
       character(len=leninp)                  ::  inp     !  Input file name
       character(len=leninp)                  ::  conf    !  Configuration file name
       character(len=lenout)                  ::  outp    !  Output file name
       logical,dimension(:,:),allocatable     ::  adj     !  Adjacency matrix
       real(kind=8),dimension(:),allocatable  ::  pop     !  Populations
       real(kind=8)                           ::  thr     !  Threshold distance
       real(kind=8)                           ::  maxdis  !  Screening distance
       integer,dimension(:),allocatable       ::  nmol    !  Number of aggregates of each size
       integer,dimension(:),allocatable       ::  imol    !  Molecule identifier
       integer,dimension(:),allocatable       ::  iagg    !  Aggregates identifier
       integer,dimension(:),allocatable       ::  itag    !  Aggregates size       
       integer                                ::  nagg    !  Number of chemical species
       integer                                ::  nnode   !  Total number of molecules
       integer                                ::  nprint  !  Populations printing interval
       integer                                ::  nsteps  !  Number of snapshots analyzed
       integer                                ::  maxagg  !  Maximum aggregate size
       integer                                ::  io      !  Status
       integer                                ::  i,j,k   !  Indexes
! Declaration of the time control variables
       real(kind=8)                           ::  tcpu    !  Total CPU time
       real(kind=8)                           ::  tread   !  Total reading time
       real(kind=8)                           ::  tadj    !  Total adjacency matrix building time
       real(kind=8)                           ::  tbfs    !  Total BFS time
       real(kind=8)                           ::  tsort   !  Total sorting time
       integer                                ::  t1,t2   !  CPU times
       integer                                ::  t1read  !  Initial reading time
       integer                                ::  t2read  !  Final reading time
       integer                                ::  t1adj   !  Initial building time
       integer                                ::  t2adj   !  Final building time
       integer                                ::  t1bfs   !  Initial BFS time
       integer                                ::  t2bfs   !  Final BFS time
       integer                                ::  t1sort  !  Initial sorting time
       integer                                ::  t2sort  !  Final sorting time
       logical                                ::  debug   !  Debug mode
! Declaration of a variable of type xtcfile
       type(xtcfile)                          ::  xtcf    !  xtc file informacion
! Declaration of a variable of type trrfile
       type(trrfile)                          ::  trr     !  trr file information
! Declaration of system_clock variables 
       integer                                ::  count_rate
       integer                                ::  count_max 
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
       tcpu  = 0.0d0
       tread = 0.0d0
       tadj  = 0.0d0
       tbfs  = 0.0d0
       tsort = 0.0d0
!
       call system_clock(count_max=count_max,count_rate=count_rate)
       call system_clock(t1)  
       call system_clock(t1read)        
!
       call print_start()
!
! Reading command line options
!
       call command_line(inp,conf,outp,nprint,maxagg,thr,maxdis,debug)
!
!  General settings and fatal errors check
!
       io = index(inp,'.')
       if ( io .eq. 0 ) inp = trim(inp)//'.xtc'
!
       io = index(conf,'.')
       if ( io .eq. 0 ) conf = trim(conf)//'.gro'
!
! Processing Gromacs input files
!
       call read_gro(conf)
!
! Printing summary of the input information
!
       write(*,'(1X,A)') 'Input information'
       write(*,'(1X,17("-"))')
       write(*,*)
       write(*,'(2X,A,2X,A)')    'Input file name         :',trim(inp)
       write(*,'(2X,A,2X,A)')    'Structure file name     :',trim(conf)
       write(*,'(2X,A,2X,A)')    'Populations file name   :',trim(outp)
       write(*,*)
       write(*,'(2X,A,4X,F4.2)') 'Threshold distance      :',thr
       write(*,'(2X,A,4X,I4)')   'Maximum aggregate size  :',maxagg
       write(*,'(2X,A,2X,I6)')   'Printing interval       :',nprint
       if ( debug ) then
         write(*,'(2X,A)')       'Debug mode              :      ON'
       else
         write(*,'(2X,A)')       'Debug mode              :     OFF'
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
         allocate(adj(nnode,nnode),imol(nnode),iagg(nnode),            &
                  itag(nnode),nmol(maxagg),pop(maxagg))
!
         nsteps = 0
         pop(:) = 0.0d0
!
         call system_clock(t2read)
!
         tread = tread + dble(t2read-t1read)/dble(count_rate)         
!
         do while ( xtcf%STAT == 0 )
           if ( mod(xtcf%STEP,nprint) .eq. 0 ) then
! Update configurations counter
             nsteps = nsteps + 1
!
! Building the adjacency matrix for the current snapshot
!
             call system_clock(t1adj)     
!
             call build_adj(debug,nnode,adj,thr,maxdis,xtcf%pos,       &
                            xtcf%NATOMS,sys%nat,                       &
                          (/xtcf%box(1,1),xtcf%box(2,2),xtcf%box(3,3)/))
!
             if ( debug ) then
               write(*,*)
               write(*,*) 'Molecule-based Adjacency matrix'
               write(*,*) '-------------------------------'
               write(*,*)
               do j = 1, nnode
                 write(*,*) (adj(i,j),i=1,nnode)
               end do
               write(*,*)
             end if
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
             if ( debug ) then
               write(*,*) 'Sorting based on the aggregate size completed'
               write(*,*) '---------------------------------------------'
               write(*,*)
               write(*,'(X,A,20I3)') 'imol-size : ',imol
               write(*,'(X,A,20I3)') 'iagg-size : ',iagg
               write(*,'(X,A,20I3)') 'itag-size : ',itag
               write(*,*)
             end if
!
! Sorting molecules based on their aggregate identifier
!
             j = 0
             do i = 1, maxagg
               if ( j .ge. nnode ) exit
               call ivqsort(nnode,iagg,imol,j+1,j+i*nmol(i))
               j = j + i*nmol(i)
             end do 
! 
             if ( debug ) then
               write(*,*) 'Sorting based on the aggregate numer completed'
               write(*,*) '----------------------------------------------'
               write(*,*)
               write(*,'(X,A,20I3)') 'imol-num  : ',imol
               write(*,'(X,A,20I3)') 'iagg-num  : ',iagg
               write(*,'(X,A,20I3)') 'itag-num  : ',itag
               write(*,*)
             end if
!
             call system_clock(t2sort)     
!
             tsort = tsort + dble(t2sort-t1sort)/dble(count_rate)
!
! Printing population of every aggregate
!
             write(uniout,'(I10,20(X,F6.2))') xtcf%STEP,real(nmol(:))/   &
                                                                nagg*100
!
             pop(:) = pop(:) + real(nmol(:))/nagg*100
!
! Printing summary of the results
! 
             if ( debug ) then
               write(*,'(X,A,20I3)') 'Total number of entities',nagg
               write(*,'(X,A,20I3)') 'Aggregates of each type ',nmol
               write(*,'(X,A,20I3)') 'New basis of molecules  ',imol
               write(*,*)
             end if
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
           call build_adj(nnode,adj,thr,trr%pos,trr%NATOMS,sys%nat,    &
                          (/trr%box(1,1),trr%box(2,2),trr%box(3,3)/),  &
                                                                  debug)
!
           if ( debug ) then 
             write(*,*)
             write(*,*) 'Molecule-based Adjacency matrix'
             write(*,*) '-------------------------------'
             write(*,*)
             do j = 1, trr%NATOMS/sys%nat
               write(*,'(3X,13I2)') (adj(i,j),i=1,trr%NATOMS/sys%nat)
             end do
             write(*,*)
           end if
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
! Printing summary of the results
!
       pop(:) = pop(:)/nsteps
       write(*,'(1X,A)') 'Output information'
       write(*,'(1X,18("-"))')
       write(*,*)
       write(*,'(1X,A,20F6.2)') 'Average populations :',pop(:)
!~        write(*,'(1X,20(X,F6.2))') pop(:)
       write(*,*)
!
! Deallocate memory
!

       deallocate(adj,imol,iagg,itag,nmol,pop)
!
       close(uniout)
!
! Printing timing
!
       call system_clock(t2)	      
       tcpu = dble(t2-t1)/dble(count_rate)
!
       write(*,'(1X,A,3(I4,A))') 'Total CPU time     ',                &
                                  int(tcpu/(60*60)),' hours',          &
                                  mod(int(tcpu/60),60),' minutes',     &
                                  mod(int(tcpu),60),' seconds'  
       write(*,'(1X,53("-"))')
!
       write(*,'(1X,A,3(I4,A))') 'Total reading time ',                &
                                  int(tread/(60*60)),' hours',         &
                                  mod(int(tread/60),60),' minutes',    &
                                  mod(int(tread),60),' seconds'  
       write(*,'(1X,A,3(I4,A))') 'Total building time',                &
                                  int(tadj/(60*60)),' hours',          &
                                  mod(int(tadj/60),60),' minutes',     &
                                  mod(int(tadj),60),' seconds'  
       write(*,'(1X,A,3(I4,A))') 'Total BFS time     ',                &
                                  int(tbfs/(60*60)),' hours',          &
                                  mod(int(tbfs/60),60),' minutes',     &
                                  mod(int(tbfs),60),' seconds'  
       write(*,'(1X,A,3(I4,A))') 'Total sorting time ',                &
                                  int(tsort/(60*60)),' hours',         &
                                  mod(int(tsort/60),60),' minutes',    &
                                  mod(int(tsort),60),' seconds'  
       write(*,*)
! Printing finishing date     
       call print_end()
!
       end program aggregate
!
!======================================================================!
!
       subroutine command_line(inp,conf,outp,nprint,maxagg,thr,maxdis, &
                               debug)
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
       character(len=leninp),intent(out)         ::  inp     !  Input file name
       character(len=lenout),intent(out)         ::  outp    !  Input file name
       character(len=leninp),intent(out)         ::  conf    !  Structure file name
       real(kind=8),intent(out)                  ::  thr     !  Threshold distance
       real(kind=8),intent(out)                  ::  maxdis  !  Screening distance
       integer,intent(out)                       ::  maxagg  !  Maximum aggregate size
       integer,intent(out)                       ::  nprint  !  Populations printing steps interval
       logical,intent(out)                       ::  debug   !  Debug mode
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
       inp    = 'md.xtc'
       conf   = 'conf.gro'
       outp   = 'md.dat'
       nprint = 1
       maxagg = 10
       maxdis = 1.5d0
       thr    = 0.25d0
       debug  = .FALSE.
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
           case ('-f','-file','-files','--file','--files') ! FLAG: change to -t,--trajectory
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
           case ('-thr','--thr','--threshold')
             call get_command_argument(i,next,status=io)
             read(next,*) thr
             i = i + 1
           case ('-a','-maxagg','--maxagg','--maximum-size')
             call get_command_argument(i,next,status=io)
             read(next,*) maxagg
             i = i + 1
           case ('-d','-maxdis','--maxdis','--maximum-distance')
             call get_command_argument(i,next,status=io)
             read(next,*) maxdis
             i = i + 1
           case ('-n','-nprint','--nprint','--print-steps')
             call get_command_argument(i,next,status=io)
             read(next,*) nprint
             i = i + 1
           case ('-v','--debug','--verbose')
             debug = .TRUE.
             i = i + 1
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
       write(*,'(2X,A)') '-f,--file             Trajectory file name'
       write(*,'(2X,A)') '-c,--configuration    Configuration file name'
       write(*,'(2X,A)') '-p,--populations      Populations file name'
       write(*,'(2X,A)') '-n,--nprint           Populations printi'//  &
                                                     'ng steps interval'
       write(*,'(2X,A)') '-a,--maxagg           Maximum aggregate size'
       write(*,'(2X,A)') '-d,--maxdis           Screening distance'
       write(*,'(2X,A)') '-thr,--threshold      Distance threshold'
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
       character(len=lenarg)                    ::  straux  !  Auxiliary string
       character(len=5)                         ::  aux     !
       integer                                  ::  io      !  Input/Output status
       integer                                  ::  i,j,k   !  Indexes
!
       open(unit=uniinp,file=trim(conf),action='read',       &
            status='old',iostat=io)
      if ( io .ne. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Missing input file'
         write(*,*)
         write(*,'(3X,3(A))')   'Input file ',trim(conf),    &
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
       return
       end subroutine read_gro
!
!======================================================================!
!
       subroutine build_adj(debug,nnode,adj,thr,maxdis,coord,natconf,  &
                            natmol,box)
!
       use geometry
!
       implicit none
!
! Input/output variables
!
       real(kind=4),dimension(3),intent(in)          ::  box      !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=4),dimension(3,natconf),intent(in)  ::  coord    !  Simulation box !FLAG: kind=8 to kind=4
       real(kind=8),intent(in)                       ::  thr      !  Threshold distance
       real(kind=8),intent(in)                       ::  maxdis   !  Screening distance
       logical,dimension(nnode,nnode),intent(out)    ::  adj      !  Adjacency matrix
       integer,intent(in)                            ::  nnode    !  Number of residues
       integer,intent(in)                            ::  natconf  !  Total number of atoms
       integer,intent(in)                            ::  natmol   !  Atoms per residue
       logical,intent(in)                            ::  debug    !  Debug mode
!
! Local variables
!
       logical,dimension(nnode)                      ::  check    !
       real(kind=4),dimension(3)                     ::  r        !  Minimum image vector !FLAG: kind=8 to kind=4
       real(kind=8)                                  ::  dist     !  Minimum image distance
       integer                                       ::  iat,jat  !  Atom indexes
       integer                                       ::  irenum   !  Residue index
       integer                                       ::  jrenum   !  Residue index
       integer                                       ::  i,j      !  Indexes
       integer                                       ::  ii,jj    !  Indexes
!
! Building adjacency matrix
!
       adj(:,:) = .FALSE.
       dist     = thr + 0.1
       i = 1
       do irenum = 1, nnode-1
         check(:) = .TRUE.
         do iat = 1, natmol
           if ( debug ) then
             write(*,*)
             write(*,'(3X,3(A,X,I4,X))') 'RESIDUE',irenum,'atom',iat,  &
                                                                   'i',i
             write(*,'(3X,29("-"))')
           end if
           jj = irenum*natmol + 1
           do jrenum = irenum + 1, nnode
             if ( (.not. adj(jrenum,irenum)) .and. check(jrenum) ) then
               j    = jj
               jat  = 1
               dist = thr + 0.1
               do while ( (jat.le.natmol) .and. (dist.le.maxdis) .and. &
                          (dist.ge.thr) ) 
                 r    = minimgvec(coord(:,i),coord(:,j),box)
                 dist = sqrt(dot_product(r,r))
!
                 if ( debug ) then
                   write(*,'(1X,4(A,X,I4,X),A,F6.3)') 'Comparing w'//  &
                         'ith residue',jrenum,'atom',jat,'j',j,'jj',   &
                                                          jj,'dist',dist
                 end if
!
                 if ( dist .lt. thr ) then
                   adj(jrenum,irenum) = .TRUE.
                   adj(irenum,jrenum) = .TRUE.
                 end if
                 j   = j + 1
                 jat = jat + 1
               end do
               if ( dist .gt. maxdis ) check(jrenum) = .FALSE.
             end if
             jj = jj + natmol
           end do
           i = i + 1
         end do
       end do
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
       integer,dimension(nnode),intent(out)          ::  imol    !  Molecule identifier
       integer,dimension(nnode),intent(out)          ::  iagg    !  Aggregate identifier
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
           if ( debug ) then
             write(*,*)
             write(*,'(X,A,I4)') 'Starting outer loop cycle in node',  &
                                                                   inode
             write(*,*)          '====================================='
             write(*,*)
             write(*,'(X,A,20L3)') 'Visited               :',notvis
             write(*,'(X,A,20I3)') 'Queue                 :',queue
             write(*,*) 'New aggregate found'
             write(*,'(X,A,I3)')   'Number of aggregates  :',nagg
             write(*,'(X,A,20I3)') 'Identifier            :',iagg
             write(*,*)
           end if
!
! Inner loop over the queue elements
!
           do while ( iqueue .lt. nqueue )
! Saving actual element in the queue
             knode = queue(iqueue)
! Checking if node k was already visited
!
               if ( debug ) then
                 write(*,*)
                 write(*,'(X,A,I4)') 'Taking from the queue node',knode
                 write(*,*)          '-------------------------------'
                 write(*,*)
                 write(*,'(3X,A,I4)') 'Starting inner loop'
                 write(*,'(3X,A)')    '-------------------'
                 write(*,*)
               end if
! Check the connection between actual queue element and the rest of nodes
               do jnode = inode + 1, nnode
                 if ( debug ) write(*,'(4X,A,I4,A)') 'Checking if '//  &
                                             'node',jnode,' was visited'
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
                   if ( debug ) then
                     write(*,*)
                     write(*,'(4X,A)')      'New connection found'
                     write(*,'(4X,A,I4)')   'Adding to the queue e'//  &
                                                          'lement',jnode
                     write(*,*)
                     write(*,'(4X,A,20I3)') 'Queue                '//  &
                                                              ' :',queue
                     write(*,'(4X,A,I3)')   'Number of aggregates '//  & 
                                                               ' :',nagg
                     write(*,'(4X,A,20I3)') 'Identifier           '//  &
                                                               ' :',iagg
                     write(*,*)
                   end if
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
       if ( debug ) then
         write(*,*)
         write(*,*) 'BFS completed'
         write(*,*) '-------------'
         write(*,*)
         write(*,'(X,A,20I3)') 'imol-raw  : ',imol
         write(*,'(X,A,20I3)') 'iagg-raw  : ',iagg
         write(*,'(X,A,20I3)') 'itag-raw  : ',itag
         write(*,*)
       end if
!
       deallocate(notvis,queue)
!
       return
       end subroutine blockdiag
!
!======================================================================!
