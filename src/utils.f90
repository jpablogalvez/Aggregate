!======================================================================!
!
       module utils
!
       use units,      only:  uniinp
!
       use printings,  only:  print_end
!
       implicit none
!
       contains
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
         write(*,'(3X,A)')    'ERROR:  No argument introduced'//       &
                              ' for command-line option'
         write(*,*)
         write(*,'(4X,A)')    trim(cmd)
         write(*,*)
         write(*,'(3X,A)') 'Argument missing for option  :  '//trim(arg)
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       return
       end subroutine check_arg
!
!======================================================================!
!
       subroutine read_string(i,lenstr,inp,nfile,arg,cmd)
!
       use lengths,  only:  lenarg
!
       implicit none
!
! Input/output variables
!
       integer,intent(inout)              ::  i       !  Argument index
       integer,intent(in)                 ::  lenstr  !  Input length
       character(len=lenstr),intent(out)  ::  inp     !  Input file name
       integer,intent(out)                ::  nfile   !  Number of input files
       character(len=*),intent(in)        ::  arg     !  Argument read
       character(len=*),intent(in)        ::  cmd     !  Command executed
!
! Local variables
!
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
       subroutine read_srealvec(i,n,box)
!
       use lengths,  only:  lenarg
!
       implicit none
!
! Input/output variables
       integer,intent(inout)                  ::  i       !  Argument index
       integer,intent(in)                     ::  n       !  Vector dimension
       real(kind=4),dimension(n),intent(out)  ::  box     !  Double precision vector
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
       end subroutine read_srealvec
!
!======================================================================!
!
       subroutine read_drealvec(i,n,box)
!
       use lengths,  only:  lenarg
!
       implicit none
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
       end subroutine read_drealvec
!
!======================================================================!
!
       subroutine read_intvec(i,n,ivec)
!
       use lengths,  only:  lenarg
!
       implicit none
!
! Input/output variables
!
       integer,intent(inout)             ::  i       !  Argument index
       integer,intent(in)                ::  n       !  Vector dimension
       integer,dimension(n),intent(out)  ::  ivec    !  Integer vector
!
! Local variables
!
       character(len=lenarg)             ::  next    !  Next argument to be read
       integer                           ::  io      !  Status
       integer                           ::  k       !  Index
!
       do k = 1, n
         call get_command_argument(i,next,status=io)
         read(next,*) ivec(k)
         i = i + 1
       end do
!
       return
       end subroutine read_intvec
!
!======================================================================!
!
       subroutine chkcomment(line,key)
!
       implicit none
!
       character(len=*),intent(in)   ::  line
       character(len=*),intent(out)  ::  key
!
       integer                       ::  posi
! Removing white spaces at the beggining
       key = adjustl(line) 
! Removing comments at the end
       posi = index(key,'!')           
       if ( posi .gt. 0 ) key = key(:posi-1)
!
       return
       end subroutine chkcomment
!
!======================================================================!
!
       subroutine chkkeyarg(key,line,arg)
!
       implicit none
!
       character(len=*),intent(in)     ::  key
       character(len=*),intent(inout)  ::  line
       character(len=*),intent(out)    ::  arg
!
       integer                         ::  posi
!
       posi = scan(line,' ') 
       if ( posi .ne. 0 ) then 
         arg  = line(:posi-1)
!
         line = line(posi+1:)  
         line = adjustl(line)
       else
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Not argument introduced for in'//  &
                                                           'put keyword'
         write(*,*) 
         write(*,'(3X,A)') 'Argument missing for keyword  :  '//       &
                                                               trim(key)
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       end if
!
       return
       end subroutine chkkeyarg
!
!======================================================================!
!
       subroutine chklineopt(line,key,arg)
!
       use lengths,  only:  lenarg,lenline
!
       implicit none
!
       character(len=lenline),intent(in)     ::  line
       character(len=lenline),intent(inout)  ::  key
       character(len=lenarg),intent(out)     ::  arg
!
       integer                               ::  posi
! Removing white spaces at the beggining
       key = adjustl(line)
! Removing comments at the end
       posi = scan(key,'!') 
       if ( posi .gt. 0 ) key = key(:posi-1)
! Saving the keyword and the values separately
       posi = scan(key,'=') 
       if ( posi .ne. 0 ) then 
         arg = key(posi+1:)  
         key = key(:posi-1)
       end if
! Changing uppercase letters by lowercase letters 
       key = lowercase(key)
!
       return
       end subroutine chklineopt
!
!======================================================================!
!
       subroutine unksect(key,blck)
!
       implicit none
!
       character(len=*),intent(in)  ::  key
       character(len=*),intent(in)  ::  blck
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Unknown section in block '//trim(blck)        
       write(*,*) 
       write(*,'(3X,A)') 'Section '//trim(key)//' not known'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine unksect
!
!======================================================================!
!
       subroutine unkopt(opt,sect,blck)
!
       implicit none
!
       character(len=*),intent(in)  ::  opt
       character(len=*),intent(in)  ::  sect
       character(len=*),intent(in)  ::  blck
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Unknown option in section '//        &
                                    trim(sect)//' of block '//trim(blck)
       write(*,*) 
       write(*,'(3X,A)') 'Option '//trim(opt)//' not known'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine unkopt
!
!======================================================================!
!
       subroutine unkkeysect(key,sect)
!
       implicit none
!
       character(len=*),intent(in)  ::  key
       character(len=*),intent(in)  ::  sect
!
       write(*,*)
       write(*,'(2X,68("="))')
!~        write(*,'(3X,A)') 'ERROR:  Unknown keyword'//                   &
!~                                      ' in option 'trim(opt)//          &
!~                                      ' of section '//trim(sect)//      &
!~                                      ' of block '//trim(blck)
       write(*,'(3X,A)') 'ERROR:  Unknown keyword in section '//       &
                                                              trim(sect)
       write(*,*) 
       write(*,'(3X,A)') 'Keyword '//trim(key)//' not known'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine unkkeysect
!
!======================================================================!
!
       subroutine errkey(aux,sect)
!
       implicit none
!
       character(len=*),intent(in)  ::  sect
       character(len=*),intent(in)  ::  aux
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Keyword badly introduced'
       write(*,*) 
       write(*,'(3X,A)') 'Please, check '//aux//' '//trim(sect)//      &
                                                    ' in the input file'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine errkey
!
!======================================================================!
!
       subroutine errkeyoptsectblck(key,opt,sect,blck)
!
       implicit none
!
       character(len=*),intent(in)  ::  key
       character(len=*),intent(in)  ::  opt
       character(len=*),intent(in)  ::  sect
       character(len=*),intent(in)  ::  blck
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Unknown keyword from input file'
       write(*,*) 
       write(*,'(3X,A)') 'Keyword '//trim(key)//' for option '//       &
                          trim(opt)//'of section '//trim(sect)//       &
                                    'in block'//trim(blck)//' not known'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()  
!
       return
       end subroutine errkeyoptsectblck
!
!======================================================================!
!
       subroutine errkeychar(aux,ch,sect)
!
       implicit none
!
       character(len=*),intent(in)  ::  sect
       character(len=*),intent(in)  ::  aux
       character(len=*),intent(in)  ::  ch
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Keyword badly introduced'
       write(*,*) 
       write(*,'(3X,A)') 'Group name specified not defined in **MO'//  &
                                                'LREP block: '//trim(ch)
       write(*,'(3X,A)') 'Please, check '//aux//' '//trim(sect)//      &
                                                    ' in the input file'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine errkeychar
!
!======================================================================!
!
       subroutine endblck(blck)
!
       implicit none
!
       character(len=*),intent(in)  ::  blck
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Input file finished before block'//  &
                                                         ' was complete'
       write(*,*) 
       write(*,'(3X,A)') 'Block '//trim(blck)//' is not complete'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine endblck
!
!======================================================================!
!
       subroutine endsect(sect)
!
       implicit none
!
       character(len=*),intent(in)  ::  sect
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Input file finished before secti'//  &
                                                       'on was complete'
       write(*,*) 
       write(*,'(3X,A)') 'Section '//trim(sect)//' is not complete'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine endsect
!
!======================================================================!
!
       subroutine endopt(opt)
!
       implicit none
!
       character(len=*),intent(in)  ::  opt
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Input file finished before optio'//  &
                                                        'n was complete'
       write(*,*) 
       write(*,'(3X,A)') 'Option '//trim(opt)//' is not complete'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine endopt
!
!======================================================================!
!
       subroutine endkey(key)
!
       implicit none
!
       character(len=*),intent(in)  ::  key
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Input file finished before keywo'//  &
                                                       'rd was complete'
       write(*,*) 
       write(*,'(3X,A)') 'Keyword '//trim(key)//' is not complete'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine endkey
!
!======================================================================!
!
       subroutine endinp(sel,flag)
!
       implicit none
!
       character(len=*),intent(in)  ::  sel
       character(len=*),intent(in)  ::  flag
!
       select case ( sel )
         case ('key')
           call endkey(flag)
         case ('opt')
           call endopt(flag)
         case ('sect')
           call endsect(flag)
         case ('blck')
           call endblck(flag)
         case default
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Subroutine ENDINP called inc'//  &
                                                              'orrectly'
           write(*,*) 
           write(*,'(3X,A)') 'Input file is not complete'
           write(*,'(2X,68("="))')
           write(*,*) 
           call print_end()           
       end select
!
       return
       end subroutine endinp
!
!======================================================================!
!
       function uppercase(aux) result(str)
!
       implicit none
!
       character(len=*),intent(in)   ::  aux
       character(len=len_trim(aux))  ::  str
!
       integer                       ::  i
!
! Changing lowercase letters by uppercase letters 
!
       str = aux 
!
       do i = 1, len_trim(aux)
         select case(str(i:i))
           case('a':'z')
             str(i:i) = achar(iachar(str(i:i))-32)
         end select
       end do 
!    
       return
       end function uppercase
!
!======================================================================!
!
       function lowercase(aux) result(str)
!
       implicit none
!
       character(len=*),intent(in)   ::  aux
       character(len=len_trim(aux))  ::  str
!
       integer                       ::  i
!
! Changing uppercase letters by lowercase letters 
!
       str = aux
!
       do i = 1, len_trim(aux)
         select case(str(i:i))
           case('A':'Z')
             str(i:i) = achar(iachar(str(i:i))+32)
         end select
       end do 
!     
       return
       end function lowercase
!
!======================================================================!
!
       subroutine findline(key,sel,flag)
!
       use lengths,  only:  lenline
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(out)  ::  key   !
       character(len=*),intent(in)   ::  sel   !
       character(len=*),intent(in)   ::  flag  !

!
! Local variables
!
       character(len=lenline)        ::  line  !
       integer                       ::  io    !  Input/Output status
!
! Reading MOLREP block sections 
!
       do
         read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
         if ( io /= 0 ) call endinp(sel,flag)
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
         return
       end do
!
       return
       end subroutine findline
!
!======================================================================!
!
! FINDCV - Find Character Vector 
!
! This subroutine returns the location of the last element in the array
!  CVEC of size N with the value given in the CVAL argument
!
       function findcv(n,cvec,cval) result(posi)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)                        ::  n
       character(len=*),dimension(n),intent(in)  ::  cvec
       character(len=*),intent(in)               ::  cval
       integer                                   ::  posi
!
! Local variables
!
       integer                                   ::  i
! 
! Finding string CVAL in the array CVEC
!
       posi = 0
       do i = 1, n
         if ( trim(adjustl(cvec(i))) .eq. trim(adjustl(cval)) ) then
           posi = i
         end if
       end do 
!     
       return
       end function findcv
!
!======================================================================!
!
       subroutine errfindcvopt(opt,sect)
!
       implicit none
!
       character(len=*),intent(in)  ::  opt
       character(len=*),intent(in)  ::  sect
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Keyword badly introduced'
       write(*,*) 
       write(*,'(3X,A)') 'Please, check option '//opt//' of section '  &
                                            //sect//' in the input file'
       write(*,'(2X,68("="))')
       write(*,*) 
       call print_end()
!
       return
       end subroutine errfindcvopt
!
!======================================================================!
!
       subroutine rndmseed(seed)
!
       implicit none
!
! Declaration of the in/out variables
!
       logical,intent(in)                ::  seed    !
!
! Declaration of the local variables
!
       integer,dimension(:),allocatable  ::  rndm    !
       integer,dimension(8)              ::  values  !
       integer                           ::  n       ! 
       integer                           ::  i,j     ! 
!
! Initializing the random seed
!       
       if ( seed ) then
         call random_seed(size=n)
         allocate(rndm(n))
!
         call date_and_time(VALUES=values)
!
         do i = 1, n
           rndm(i) = 0
           do j = 1, 8
             rndm(i) = rndm(i) + values(j)*i
           end do
         end do
!
         call random_seed(put=rndm)
       else
         call random_seed(size=n)
!
         allocate(rndm(n))
!
         do i = 1, n
           rndm(i) = i*n
         end do
!
         call random_seed(put=rndm)
       end if
!
       deallocate(rndm)
!
       return
       end subroutine rndmseed
!
!======================================================================!
!
       subroutine znum2atname(znum,atname)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)            ::  znum    !  Atomic number
       character(len=5),intent(out)  ::  atname  !  Atom names
! 
!
!
       select case ( znum )
         case ( 1 )
           atname = 'H'
         case ( 2 )
           atname = 'He'
         case ( 3 )
           atname = 'Li'
         case ( 4 )
           atname = 'Be'
         case ( 5 )
           atname = 'B'
         case ( 6 )
           atname = 'C'
         case ( 7 )
           atname = 'N'
         case ( 8 )
           atname = 'O'
         case ( 9 )
           atname = 'F'
         case ( 10 )
           atname = 'Ne'
         case ( 17 )
           atname = 'Cl'
         case ( 18 )
           atname = 'Ar'
         case ( 36 )
           atname = 'Kr'
         case default
           write(*,*) znum, 'Not yet!'
           write(*,*) 
           call print_end()
       end select
!
       return
       end subroutine znum2atname
!
!======================================================================!
!
       subroutine znum2atmass(znum,atmass)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)        ::  znum    !  Atomic number
       real(kind=8),intent(out)  ::  atmass  !  Atom mass
! 
!
!
       select case ( znum )
         case ( 1 )
           atmass = 1.00783   ! Gaussian value
!~          case ( 2 )
!~            atmass = 'He'
!~          case ( 3 )
!~            atmass = 'Li'
!~          case ( 4 )
!~            atmass = 'Be'
!~          case ( 5 )
!~            atmass = 'B'
         case ( 6 )
           atmass = 12.0d0    ! Gaussian value
!~          case ( 7 )
!~            atmass = 'N'
         case ( 8 )
           atmass = 15.99491  ! Gaussian value
!~          case ( 9 )
!~            atmass = 'F'
!~          case ( 10 )
!~            atmass = 'Ne'
!~          case ( 17 )
!~            atmass = 'Cl'
!~          case ( 18 )
!~            atmass = 'Ar'
!~          case ( 36 )
!~            atmass = 'Kr'
         case default
           write(*,*) znum, 'Not yet!'
           write(*,*) 
           call print_end()
       end select
!
       return
       end subroutine znum2atmass
!
!======================================================================!
!
       subroutine countlines(inp,uni,n)
!
       use printings, only: print_missinp
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp
       integer,intent(in)           ::  uni
       integer,intent(out)          ::  n
!
! Local variables
!
       integer                      ::  io
!
! Counting the number of lines up to the end of the file
!
       open(unit=uni,file=inp,action='read',status='old',iostat=io)
!
       if ( io .ne. 0 ) call print_missinp(inp)
!
       n = 0
       do
         read(uni,*,iostat=io)
         if ( io .ne. 0 ) exit
         n = n + 1
       end do
!
       close(uni)
!
       return
       end subroutine countlines
!
!======================================================================!
!
       subroutine find_key(uni,key,line,iost)
!
       use lengths, only: lenline
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  
       character(len=lenline),intent(out)  ::  line    !  Line read
       integer,intent(out)                 ::  iost    !  Reading status
       integer,intent(in)                  ::  uni     !  Reading status
!
! Local variables
!
       integer                             ::  keylen  !  Keyword length
!
! Finding the first line starting with the specified keyword 
!
       keylen = len_trim(key)
       iost   = 0
!
       do
         read(uni,'(A)',iostat=iost) line
!~ write(*,*)'reading:',trim(line),iost
         if ( iost /= 0 ) exit
         if ( len_trim(line) .eq. 0 )  cycle
         if ( len_trim(line) .lt. keylen ) cycle
!
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) return
       end do
!~ write(*,*) 'out'
!
       iost = 1
!
       return
       end subroutine find_key
!
!======================================================================!
!
       subroutine below_key(uni,key,iost)
!
       use lengths, only: lenline
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  
       integer,intent(out)                 ::  iost    !  Reading status
       integer,intent(in)                  ::  uni     !  Reading status
!
! Local variables
!
character(len=100) ::  ioerrmsg
       character(len=lenline)              ::  line    !  Line read
       integer                             ::  keylen  !  Keyword length
!
! Finding the first line starting with the specified keyword 
!
       keylen = len_trim(key)
       iost   = 0
!
       do
         read(uni,'(A)',iostat=iost,iomsg=ioerrmsg) line
!~ write(*,*)'reading:',trim(line),iost!,trim(ioerrmsg)
         if ( iost /= 0 ) exit
         if ( len_trim(line) .eq. 0 )  cycle
         if ( len_trim(line) .lt. keylen ) cycle
!
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) return
       end do
!~ write(*,*) 'out'
!
       iost = 1
!
       return
       end subroutine below_key
!
!======================================================================!
!
       subroutine find_last(uni,key,line,iost)
!
       use lengths, only: lenline
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  
       character(len=lenline),intent(out)  ::  line    !  Line read
       integer,intent(in)                  ::  uni     !
       integer,intent(out)                 ::  iost    !  Reading status
!
! Local variables
!
       character(len=lenline)              ::  straux  !  Auxliary string
       integer                             ::  keylen  !  Keyword length
       integer                             ::  io      !  Input/Output status
!
! Finding the last line starting with the specified keyword 
!
       keylen = len_trim(key)
       iost   = 1
!
       do
         read(uni,'(A)',iostat=io) straux
         if ( io /= 0 ) exit
         straux = adjustl(straux)
         if ( straux(:keylen) .eq. trim(key) ) then
           iost = 0 
           line = straux
         end if
       end do
!
       return
       end subroutine find_last
!
!======================================================================!
!
       subroutine print_badread(inp,key)
!
       use lengths, only: leninp
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)       ::  key     !  Input file name
       character(len=leninp),intent(in)  ::  inp     !  Input file name
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)')      'ERROR:  Missing information in the '//  &
                              'input file'
       write(*,*)
       write(*,'(3X,2(A))')   'Keyword    : ',trim(adjustl(key))
       write(*,'(3X,2(A))')   'Input file : ',trim(inp)
       write(*,'(2X,68("="))')
       write(*,*)
       call exit(0)
!
       return
       end subroutine print_badread
!
!======================================================================!
!
       end module utils
!
!======================================================================!
