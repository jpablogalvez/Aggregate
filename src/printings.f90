!======================================================================!
!
       module printings
!
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
       use lengths,  only:  lencmd
       use version
!
       implicit none
!
       character(len=lencmd)  ::  cmd  !  Command executed
!
       write(*,'(1X,A)') 'Starting program at '//fdate()//             &
                         ' on '//print_host() 
       write(*,*)  
!
       call print_version()
       write(*,*)  
!
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
       write(*,'(1X,A)') 'Finishing program at '//fdate()//' on '//    &
                                                            print_host()
       write(*,*)
!     
       call exit(0)
!
       return
       end subroutine print_end
!
!======================================================================!
!
       subroutine print_time(uni,blnk,key,lenin,time)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       real(kind=8),intent(in)      ::  time     !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
       integer,intent(in)           ::  lenin    !   
!
! Local variables
!
       character(len=64)            ::  straux   !
       character(len=256)           ::  fmt1     !
       integer                      ::  lenmid   !
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk - 1   ! FLAG: check if iaux is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,5(1X,I2,1X,A))'
!
       write(uni,fmt1) key,int(time/(60*60*24)),'days',                &
                              int(time/(60*60)),'h',                   &
                           mod(int(time/60),60),'min',                 &
                              mod(int(time),60),'sec',                 &  
                      int(100*(time-int(time))),'msec'
!
       return
       end subroutine print_time
!
!======================================================================!
!
       subroutine print_speed(uni,blnk,key,lenin,time,cputime)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       real(kind=8),intent(in)      ::  time     !
       real(kind=8),intent(in)      ::  cputime  !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
       integer,intent(in)           ::  lenin    !   
!
! Local variables
!
       character(len=64)            ::  straux   !
       character(len=256)           ::  fmt1     !
       real(kind=8)                 ::  speed    !
       integer                      ::  lenmid   !
!
       speed = cputime/time
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk - 1   ! FLAG: check if iaux is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,5(1X,I2,1X,A),F4.1,A)'
!
       write(uni,fmt1) key,int(time/(60*60*24)),'days',                &
                              int(time/(60*60)),'h',                   &
                           mod(int(time/60),60),'min',                 &
                              mod(int(time),60),'sec',                 &  
                      int(100*(time-int(time))),'msec (',speed,')'
!
       return
       end subroutine print_speed
!
!======================================================================!
!
       subroutine print_missinp(inp)
!
       implicit none
!
       character(len=*),intent(in)  ::  inp
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Missing input file'
       write(*,*)
       write(*,'(3X,A)') 'Input file '//trim(inp)//' not found in '//  &
                                                 'the current directory'
       write(*,'(2X,68("="))')
       write(*,*)
       call print_end()
!
       return
       end subroutine print_missinp
!
!======================================================================!
!
       subroutine print_inp(uni,inp,lfin)
!
       use lengths, only: lenline
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp
       integer,intent(in)           ::  uni
       integer,intent(in)           ::  lfin
!
! Local variables
!
       character(len=16)            ::  aux
       character(len=lenline)       ::  line
       integer                      ::  io
!
! Printing input data file
!
       write(aux,*) lfin
       aux = adjustl(aux)
       aux = '('//trim(aux)//'("."))'
!
       open(unit=uni,file=trim(inp),action='read',status='old',        &
            iostat=io)
       if ( io .ne. 0 ) call print_missinp(inp)
!
       write(6,'(2X,A)') 'Dump of input data file'
       write(*,*)
       write(*,aux)
!
       do
         read(uni,'(A)',iostat=io) line 
         if ( io /= 0 ) exit
         write(*,'(A)') trim(line)
       end do
!
       write(*,aux)
       write(*,*)
!
       close(uni)
!
       return
       end subroutine print_inp
!
!======================================================================!
!
       subroutine line_sp(uni,blnk,key,lenin,sep,sizedp,sval,lenfin)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       character(len=*),intent(in)  ::  sizedp   !  
       character(len=*),intent(in)  ::  sep      !    
       real(kind=4),intent(in)      ::  sval     !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
       integer,intent(in)           ::  lenin    !   
       integer,intent(in)           ::  lenfin   !   
!
! Local variables
!
       character(len=64)            ::  straux   !
       character(len=256)           ::  fmt1     !
       integer                      ::  iaux     !
       integer                      ::  lenmid   !
       integer                      ::  io       !
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk    ! FLAG: check if iaux is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,("'//sep//'"),'
!
       io = scan(sizedp,'.')
       if ( io .eq. 0 ) then 
          write(*,*)
          write(*,'(2X,68("="))')
          write(*,'(3X,A)') 'ERROR:  Subroutine LINE_DP called inc'//  &
                                                              'orrectly'
          write(*,*) 
          write(*,'(3X,A)') 'Error while printing information'
          write(*,'(2X,68("="))')
          write(*,*) 
          call print_end() 
       end if
!
       iaux = io - 2
       write(straux,*) iaux
       straux = adjustl(straux)
       straux = '(I'//trim(straux)//')'
!
       read(sizedp(2:io-1),straux) iaux
       iaux = lenfin - lenin - iaux - len(sep)  ! FLAG: check if iaux is negative
       write(straux,*) iaux
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,'//sizedp//')'
!
       write(uni,fmt1) key,sval
!
       return
       end subroutine line_sp
!
!======================================================================!
!
       subroutine line_dp(uni,blnk,key,lenin,sep,sizedp,dval,lenfin)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       character(len=*),intent(in)  ::  sizedp   !  
       character(len=*),intent(in)  ::  sep      !    
       real(kind=8),intent(in)      ::  dval     !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
       integer,intent(in)           ::  lenin    !   
       integer,intent(in)           ::  lenfin   !   
!
! Local variables
!
       character(len=64)            ::  straux   !
       character(len=256)           ::  fmt1     !
       integer                      ::  iaux     !
       integer                      ::  lenmid   !
       integer                      ::  io       !
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk    ! FLAG: check if iaux is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,("'//sep//'"),'
!
       io = scan(sizedp,'.')
       if ( io .eq. 0 ) then 
          write(*,*)
          write(*,'(2X,68("="))')
          write(*,'(3X,A)') 'ERROR:  Subroutine LINE_DP called inc'//  &
                                                              'orrectly'
          write(*,*) 
          write(*,'(3X,A)') 'Error while printing information'
          write(*,'(2X,68("="))')
          write(*,*) 
          call print_end() 
       end if
!
       iaux = io - 2
       write(straux,*) iaux
       straux = adjustl(straux)
       straux = '(I'//trim(straux)//')'
!
       read(sizedp(2:io-1),straux) iaux
       iaux = lenfin - lenin - iaux - len(sep)  ! FLAG: check if iaux is negative
       write(straux,*) iaux
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,'//sizedp//')'
!
       write(uni,fmt1) key,dval
!
       return
       end subroutine line_dp
!
!======================================================================!
!
       subroutine line_dvec(uni,blnk,key,lenin,sep,ndim,vecblnk,       &
                            sizedp,dvec,lenfin)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)              ::  key      !   
       character(len=*),intent(in)              ::  sizedp   !  
       character(len=*),intent(in)              ::  sep      !    
       real(kind=8),dimension(ndim),intent(in)  ::  dvec     !
       integer,intent(in)                       ::  ndim     !   
       integer,intent(in)                       ::  blnk     !   
       integer,intent(in)                       ::  vecblnk  !   
       integer,intent(in)                       ::  uni      !   
       integer,intent(in)                       ::  lenin    !   
       integer,intent(in)                       ::  lenfin   !   
!
! Local variables
!
       character(len=256)                       ::  fmt1     !
       character(len=64)                        ::  straux   !
       integer                                  ::  iaux     !
       integer                                  ::  lenmid   !
       integer                                  ::  io       !
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk    ! FLAG: check if iaux is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,("'//sep//'"),'
!
       io = scan(sizedp,'.')
       if ( io .eq. 0 ) then 
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Subroutine LINE_DVEC called in'//  &
                                                             'correctly'
         write(*,*) 
         write(*,'(3X,A)') 'Error while printing information'
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end() 
       end if
!
       iaux = io - 2
       write(straux,*) iaux
       straux = adjustl(straux)
       straux = '(I'//trim(straux)//')'
!
       read(sizedp(2:io-1),straux) iaux
       iaux = lenfin - lenin - ndim*iaux - len(sep) - ndim*vecblnk  ! FLAG: check if iaux is negative
       write(straux,*) iaux
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,'
!
       write(straux,*) ndim
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'('
!
       write(straux,*) vecblnk
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,'//sizedp//'))'
!
       write(uni,fmt1) key,dvec
!
       return
       end subroutine line_dvec
!
!======================================================================!
!
       subroutine line_ivec(uni,blnk,key,lenin,sep,ndim,vecblnk,       &
                            sizeint,ivec,lenfin)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key      !   
       character(len=*),intent(in)         ::  sizeint  !  
       character(len=*),intent(in)         ::  sep      !    
       integer,dimension(ndim),intent(in)  ::  ivec     !
       integer,intent(in)                  ::  ndim     !   
       integer,intent(in)                  ::  blnk     !   
       integer,intent(in)                  ::  vecblnk  !   
       integer,intent(in)                  ::  uni      !   
       integer,intent(in)                  ::  lenin    !   
       integer,intent(in)                  ::  lenfin   !   
!
! Local variables
!
       character(len=256)                  ::  fmt1     !
       character(len=64)                   ::  straux   !
       integer                             ::  iaux     !
       integer                             ::  lenmid   !
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk    ! FLAG: check if iaux is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,("'//sep//'"),'
!
!
       read(sizeint(2:),*) iaux
       iaux = lenfin - lenin - ndim*iaux - len(sep) - ndim*vecblnk  ! FLAG: check if iaux is negative
       write(straux,*) iaux
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,'
!
       write(straux,*) ndim
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'('
!
       write(straux,*) vecblnk
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,'//sizeint//'))'
!
       write(uni,fmt1) key,ivec
!
       return
       end subroutine line_ivec
!
!======================================================================!
!
       subroutine line_int(uni,blnk,key,lenin,sep,sizeint,ival,lenfin)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       character(len=*),intent(in)  ::  sizeint  !   
       character(len=*),intent(in)  ::  sep      !   
       integer,intent(in)           ::  ival     !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
       integer,intent(in)           ::  lenin    !   
       integer,intent(in)           ::  lenfin   !   
!
! Local variables
!
       character(len=256)           ::  fmt1     !
       character(len=64)            ::  straux   !
       integer                      ::  iaux     !
       integer                      ::  lenmid   !
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk             ! FLAG: check if lenmid is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,("'//sep//'"),'
!
       read(sizeint(2:),*) iaux
       iaux = lenfin - lenin - iaux - len(sep)      ! FLAG: check if iaux is negative
       write(straux,*) iaux
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,'//sizeint//')'
!
       write(uni,fmt1) key,ival
!
       return
       end subroutine line_int
!
!======================================================================!
!
       subroutine line_str(uni,blnk,key,lenin,sep,strval,lenfin)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       character(len=*),intent(in)  ::  strval   !
       character(len=*),intent(in)  ::  sep      !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
       integer,intent(in)           ::  lenin    !   
       integer,intent(in)           ::  lenfin   !   
!
! Local variables
!
       character(len=256)           ::  fmt1     !
       character(len=64)            ::  straux   !
       integer                      ::  iaux     !
       integer                      ::  lenmid   !
!
       write(straux,'(I4)') blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A,'
!
       lenmid = lenin - len(key) - blnk                ! FLAG: check if iaux is negative
       write(straux,*) lenmid
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,("'//sep//'"),'
!
       iaux = lenfin - lenin - len(strval) - len(sep)      ! FLAG: check if iaux is negative
       write(straux,*) iaux
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'X,A)'
!
       write(uni,fmt1) key,strval
!
       return
       end subroutine line_str
!
!======================================================================!
!
       subroutine line_log(uni,blnk,key,lenin,sep,logval,lenfin)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       character(len=*),intent(in)  ::  sep      !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
       integer,intent(in)           ::  lenin    !   
       integer,intent(in)           ::  lenfin   !  
       logical,intent(in)           ::  logval   ! 
!
       if ( logval ) then
         call line_str(uni,blnk,key,lenin,sep,'YES',lenfin)
       else
         call line_str(uni,blnk,key,lenin,sep,'NO',lenfin)
       end if
!
       return
       end subroutine line_log
!
!======================================================================!
!
       subroutine print_title(uni,blnk,key,sub)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       character(len=*),intent(in)  ::  sub      !
       integer,intent(in)           ::  blnk     !   
       integer,intent(in)           ::  uni      !   
!
! Local variables
!
       character(len=256)           ::  fmt1     !
       character(len=64)            ::  straux   !
       integer                      ::  iaux     !
!
! Printing title line
!
       write(straux,*) blnk
       straux = adjustl(straux)
!
       fmt1 = '('//trim(straux)//'X,A)'
!
       write(uni,fmt1) key
!
! Printing highlighting line
!
       fmt1 = '('//trim(straux)//'X,'
!
       iaux = floor(real(len(key))/len(sub))   ! FLAG: check when len(sub) gt 1
       write(straux,*) iaux
       straux = adjustl(straux)
!
       fmt1 = trim(fmt1)//trim(straux)//'("'//sub//'"))'
!
       write(uni,fmt1)
!
       return
       end subroutine print_title
!
!======================================================================!
!
       subroutine print_titleint(uni,blnk1,key,blnk2,sizeint,ival,sub)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  key      !   
       character(len=*),intent(in)  ::  sizeint  !   
       character(len=*),intent(in)  ::  sub      !
       integer,intent(in)           ::  blnk1    !   
       integer,intent(in)           ::  blnk2    !   
       integer,intent(in)           ::  ival     !   
       integer,intent(in)           ::  uni      !   
!
! Local variables
!
       character(len=256)           ::  fmt1     !
       character(len=64)            ::  straux1  !
       character(len=64)            ::  straux2  !
       integer                      ::  iaux     !
!
! Printing information line
!
       write(straux1,*) blnk1
       straux1 = adjustl(straux1)
!
       fmt1 = '('//trim(straux1)//'X,A,'
!
       write(straux2,*) blnk2
       straux2 = adjustl(straux2)
!
       fmt1 = trim(fmt1)//trim(straux2)//'X,'//sizeint//')'
!
       write(uni,fmt1) key,ival
!
! Printing highlighting line
!
       fmt1 = '('//trim(straux1)//'X,'
!
       read(sizeint(2:),*) iaux
       iaux = iaux + blnk2 + floor(real(len(key))/len(sub))   ! FLAG: check if iaux is negative
       write(straux2,*) iaux
       straux2 = adjustl(straux2)
!
       fmt1 = trim(fmt1)//trim(straux2)//'("'//sub//'"))'
!
       write(uni,fmt1)
!
       return
       end subroutine print_titleint

!
!======================================================================!
!
       subroutine print_ivec(i,n,A)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)               ::  i       !
       integer,intent(in)               ::  n       !
       integer,intent(in),dimension(n)  ::  A       !
!
! Local variables
!
       integer,parameter                ::  num=10  !
       integer                          ::  ilower  !
       integer                          ::  iupper  !
       integer                          ::  j       !
!
       do ilower = 1, n, num
         iupper = min(ilower + num - 1,n)
         write(*,'(1X,10(X,I6))') (i+j,j=ilower,iupper)
         write(*,'(1X,10(X,I6))') (A(j),j=ilower,iupper)
         write(*,*)
       end do
!
       return
       end subroutine print_ivec
!
!======================================================================!
!
       subroutine print_dvec(i,n,A)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)                    ::  i       !
       integer,intent(in)                    ::  n       !
       real(kind=8),intent(in),dimension(n)  ::  A       !
!
! Local variables
!
       integer,parameter                     ::  num=10  !
       integer                               ::  ilower  !
       integer                               ::  iupper  !
       integer                               ::  j       !
!
       do ilower = 1, n, num
         iupper = min(ilower + num - 1,n)
         write(*,'(1X,10(X,I6))') (i+j,j=ilower,iupper)
         write(*,'(1X,10(X,I6))') (A(j),j=ilower,iupper)
         write(*,*)
       end do
!
       return
       end subroutine print_dvec
!
!======================================================================!
!
       subroutine print_info(i,n,A,B,C,str1,str2,str3)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)               ::  i       !
       integer,intent(in)               ::  n       !
       integer,intent(in),dimension(n)  ::  A,B,C   !
       character(len=*),intent(in)      ::  str1    !   
       character(len=*),intent(in)      ::  str2    !   
       character(len=*),intent(in)      ::  str3    !   
!
! Local variables
!
       integer,parameter                ::  num=10  !
       integer                          ::  ilower  !
       integer                          ::  iupper  !
       integer                          ::  j       !
!
       do ilower = 1, n, num
         iupper = min(ilower + num - 1,n)
         write(*,'(11X,10(X,I6))')    (i+j,j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str1,(A(j),j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str2,(B(j),j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str3,(C(j),j=ilower,iupper)
         write(*,*)
       end do
!
       return
       end subroutine print_info
!
!======================================================================!
!
       subroutine print_dictionary(i,n,A,m,B,str1,str2)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)                 ::  i       !
       integer,intent(in)                 ::  n       !
       integer,intent(in)                 ::  m       !
       integer,intent(in),dimension(n)    ::  A       !
       integer,intent(in),dimension(m,n)  ::  B       !
       character(len=*),intent(in)        ::  str1    !   
       character(len=*),intent(in)        ::  str2    !   
!
! Local variables
!
       integer,parameter                ::  num=10  !
       integer                          ::  ilower  !
       integer                          ::  iupper  !
       integer                          ::  j,k     !
!
       do ilower = 1, n, num
         iupper = min(ilower + num - 1,n)
         write(*,'(11X,10(X,I6))')    (i+j,j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str1,(A(j),j=ilower,iupper)
         do k = 1, m
           write(*,'(1X,A10,10(X,I6))') str2,(B(k,j),j=ilower,iupper)
         end do
         write(*,*)
       end do
!
       return
       end subroutine print_dictionary
!
!======================================================================!
!
       subroutine print_test(n,i,m,A,B,C,D,str1,str2,str3,str4)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)               ::  i       !
       integer,intent(in)               ::  n       !
       integer,intent(in)               ::  m       !
       integer,intent(in),dimension(n)  ::  A       !
       integer,intent(in),dimension(m)  ::  B,C,D   !
       character(len=*),intent(in)      ::  str1    !   
       character(len=*),intent(in)      ::  str2    !   
       character(len=*),intent(in)      ::  str3    !   
       character(len=*),intent(in)      ::  str4    !   
!
! Local variables
!
       integer,parameter                ::  num=10  !
       integer                          ::  ilower  !
       integer                          ::  iupper  !
       integer                          ::  j       !
!
       do ilower = 1, n, num
         iupper = min(ilower + num - 1,n)
         write(*,'(11X,10(X,I6))')    (j,j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str1,(A(j),j=ilower,iupper)
         write(*,*)
       end do
!
       do ilower = 1, m, num
         iupper = min(ilower + num - 1,m)
         write(*,'(11X,10(X,I6))')    (i+j,j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str2,(B(j),j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str3,(C(j),j=ilower,iupper)
         write(*,'(1X,A10,10(X,I6))') str4,(D(j),j=ilower,iupper)
         write(*,*)
       end do
!
       return
       end subroutine print_test
!
!======================================================================!
!
       end module printings
!
!======================================================================!
