!======================================================================!
!
       module utils
       implicit none
!
       contains
!
!======================================================================!
!
       character(len=32) function print_host()
!
       implicit none
!
       call hostnm(print_host)

       return
       end function print_host
!
!======================================================================!
!
       subroutine print_start()
!
       implicit none
!
       include 'info.h'
!
       character(len=lencmd)  ::  cmd         !  Command executed
!
       write(*,'(1X,A)') 'Starting program at '//fdate()//             &
                         ' on '//print_host() 
       write(*,*)  
!
       write(*,'(2X,A)') 'Executing version '//trim(version)
       call get_command(cmd) 
       write(*,'(2X,A)') 'Command executed:'
       write(*,'(4X,A)') trim(cmd)
       write(*,*)   
!
       return
       end subroutine print_start
!
!======================================================================!
!
       subroutine print_end()
!
       implicit none
!
       write(*,'(1X,4(A))') 'Finishing program at ', fdate(),          &
                            ' on ', print_host()
       write(*,*)
!     
       call exit(0)
!
       return
       end subroutine print_end
!
!======================================================================!
!
       subroutine check_arg(opt,io,arg,cmd)
!
       implicit none
!
       character(len=*),intent(in)  ::  opt
       character(len=*),intent(in)  ::  arg
       character(len=*),intent(in)  ::  cmd
       integer,intent(in)           ::  io
!
       if ( (io .ne. 0) .or. (opt(1:1) .eq. '-') ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  No argument introduced'//      &
                              ' for command-line'//   &
                              ' option'
         write(*,*)
         write(*,'(4X,A)')    trim(cmd)
         write(*,*)
         write(*,'(3X,2(A))') 'Argument missing for command'//        &
                              '-line option  :  ', arg
         write(*,'(2X,68("="))')
         write(*,*)
         call exit(0)
       end if
!
       return
       end subroutine check_arg
!
!======================================================================!
!
       subroutine read_string(i,lenstr,inp,nfile,arg,cmd)
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
       integer,intent(inout)              ::  i       !  Argument index
       integer,intent(in)                 ::  lenstr  !  Input length
       character(len=lenstr),intent(out)  ::  inp     !  Input file name
       integer,intent(out)                ::  nfile   !  Number of input files
       character(len=lenarg),intent(in)   ::  arg     !  Argument read
       character(len=lencmd),intent(in)   ::  cmd     !  Command executed
! Local variables
       character(len=lenarg)              ::  next    !  Next argument to be read
       integer                            ::  io      !  Status
!
       call get_command_argument(i,next,status=io)
       call check_arg(next,io,arg,cmd)
       nfile = 1
       inp   = next
       do
         call get_command_argument(i+nfile,next,status=io)
         if ( (io .ne. 0) .or. (next(1:1) .eq. '-') ) exit
         nfile = nfile + 1
         inp   = trim(inp)//' '//trim(next)
       end do
       i = i + nfile
!
       return
       end subroutine read_string
!
!======================================================================!
!
       subroutine read_realvec(i,n,box)
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
       integer,intent(inout)                  ::  i       !  Argument index
       integer,intent(in)                     ::  n       !  Vector dimension
       real(kind=8),dimension(n),intent(out)  ::  box     !  Double precision vector
! Local variables
       character(len=lenarg)                  ::  next    !  Next argument to be read
       integer                                ::  io      !  Status
       integer                                ::  k       !  Index
!
       do k = 1, n
         call get_command_argument(i,next,status=io)
         read(next,*) box(k)
         i = i + 1
       end do
!
       return
       end subroutine read_realvec
!
!======================================================================!
!
       subroutine read_intvec(i,n,nbox)
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
       integer,intent(inout)             ::  i       !  Argument index
       integer,intent(in)                ::  n       !  Vector dimension
       integer,dimension(n),intent(out)  ::  nbox    !  Integer vector
! Local variables
       character(len=lenarg)             ::  next    !  Next argument to be read
       integer                           ::  io      !  Status
       integer                           ::  k       !  Index
!
       do k = 1, n
         call get_command_argument(i,next,status=io)
         read(next,*) nbox(k)
         i = i + 1
       end do
!
       return
       end subroutine read_intvec
!
!======================================================================!
!
       end module utils
!
!======================================================================!
