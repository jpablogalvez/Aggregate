!======================================================================!
!
       module input_section
!
       implicit none
!
       private
       public  ::  read_inp,                                           &
                   read_gro
!
       contains
!
!======================================================================!
!
       subroutine read_gro(conf,sys)
!
       use datatypes
       use lengths, only: leninp
       use printings
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  conf     !  Structure file name
       type(groinp),intent(inout)        ::  sys      !  System information
!
! Local variables
!
       character(len=20)                 ::  straux
       character(len=20)                 ::  aux
       integer                           ::  io      !  Input/Output status
       integer                           ::  i,j       !  Index
!
! Reading Gromacs configuration file
! ----------------------------------
!
       open(unit=uniinp,file=trim(conf),action='read',                 &
            status='old',iostat=io)
!
       if ( io .ne. 0 ) call print_missinp(conf)
!
       read(uniinp,'(A)') sys%title
       read(uniinp,*)     sys%nat
!
       allocate(sys%renum(sys%nat)   ,  &
                sys%rename(sys%nat)  ,  &
                sys%atname(sys%nat)  ,  &
                sys%atnum(sys%nat)   ,  &
                sys%mass(sys%nat))
!
       do i = 1, sys%nat       ! TODO: check if lines are correctly read
         read(uniinp,'(I5,2A5,I5,3F8.3)') sys%renum(i),   &
                                          sys%rename(i),  &
                                          sys%atname(i),  &
                                          sys%atnum(i)
       end do
       read(uniinp,*) sys%latvec
!
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
!~ !
         sys%atname(j) = straux(:5)
         straux = uppercase(straux)
!
         select case ( straux )
           case ( 'H' )
             sys%mass(j) = 1.007825d0
           case ( 'HE' )
             sys%mass(j)   = 4.002602d0  ! Not exact
             sys%atname(j) = 'He'
           case ( 'LI' )
             sys%mass(j) = 6.941d0     ! Not exact
             sys%atname(j) = 'Li'
           case ( 'BE' )
             sys%mass(j) = 9.012182d0  ! Not exact
             sys%atname(j) = 'Be'
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
             sys%atname(j) = 'Ne'
           case ( 'CL' )
             sys%mass(j) = 35.453d0    ! Not exact
             sys%atname(j) = 'Cl'
           case ( 'AR' )
             sys%mass(j) = 39.948d0    ! Not exact
             sys%atname(j) = 'Ar'
           case ( 'KR' )
             sys%mass(j) = 83.798d0    ! Not exact
             sys%atname(j) = 'Kr'
           case default
             write(*,*) straux, 'Not yet!'
             call exit(0)
         end select
       end do
!
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
       subroutine read_inp(inp,nat,tgrp,grptag,body,grps,subg,atms,    &
                           mbody,mgrps,msubg,matms,thr,thr2,thrang,    &
                           neiang)
!
       use lengths, only: leninp,lentag,lenline
       use printings
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)                  ::  inp     !  General input file name
       character(len=leninp),intent(out)                 ::  tgrp    !  Groups file title
       character(len=lentag),dimension(nat),intent(out)  ::  grptag  !  Names of the groups
       real(kind=8),dimension(nat,nat),intent(out)       ::  thr     !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(out)       ::  thr2    !  Distance threshold
       real(kind=8),dimension(nat,nat),intent(out)       ::  thrang  !  Angle threshold
       integer,dimension(nat),intent(out)                ::  neiang  !  First neighbour index
       integer,dimension(nat),intent(out)                ::  body    !  Number of groups in each body
       integer,dimension(nat),intent(out)                ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(out)                ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(out)                ::  atms    !  Atoms identifier
       integer,intent(out)                               ::  mbody   !  Number of bodies
       integer,intent(out)                               ::  mgrps   !  Number of groups
       integer,intent(out)                               ::  msubg   !  Number of subgroups
       integer,intent(out)                               ::  matms   !  
       integer,intent(in)                                ::  nat     !  Number of atoms in the molecule
!
! Local variables
!
       character(len=lenline)                            ::  line    !
       character(len=lenline)                            ::  key     !
       integer                                           ::  io      !  Input/Output status
       integer                                           ::  i       !  Indexes
!
! Setting defaults
!
       thr(:,:)  = 0.0d0
       thr2(:,:) = 0.0d0
!
       thrang(:,:) = 0.0d0
!
       neiang(:) = 0
!
       tgrp = 'MOLREP Title'
!
       mbody    = nat
       mgrps    = nat
       msubg    = nat
       matms    = nat
!
       do i = 1, nat
         write(grptag(i),'(I8)') i
         grptag(i) = 'Atom-'//trim(adjustl(grptag(i)))
         body(i)   = i
         grps(i)   = i
         subg(i)   = i
         atms(i)   = i
       end do
!
! Reading general input file
! --------------------------
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(inp)
!
! Reading input file block 
!   
       do
! Reading input file line
         read(uniinp,'(A)',iostat=io) line 
! If the end of the input file is reached exit
         if ( io /= 0 ) exit
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
           case ('**MOLREP','**MOLREPRE')
!~              write(*,*) 
!~              write(*,*) 'Reading **MOLREP block'
!~              write(*,*) 
!
             call findline(line,'blck','**MOLREP')
!
             call read_molrep(line,'**MOLREP',nat,tgrp,grptag,body,    &
                              grps,subg,atms,mbody,mgrps,msubg,matms)      
!      
           case ('**THRESH','**THRESHOLD')
!~              write(*,*) 
!~              write(*,*) 'Reading **THRESHOLD block'
!~              write(*,*) 
!
             call findline(line,'blck','**THRESHOLD')
!
             call read_threshold(line,'**THRESHOLD',nat,thr,mgrps,     &
                                 grptag)
!
           case ('**INTERTHRESH','**INTERTHRESHOLD')
!~              write(*,*) 
!~              write(*,*) 'Reading **THRESHOLD block'
!~              write(*,*) 
!
             call findline(line,'blck','**INTERTHRESHOLD')
!
             call read_interthreshold(line,'**INTERTHRESHOLD',nat,thr2,&
                                      mgrps,grptag)
!
           case ('**ANGLES','**THRANG','**THRANGLE','**THREANGLES')
!~              write(*,*) 
!~              write(*,*) 'Reading **ANGLES block'
!~              write(*,*) 
!
             call findline(line,'blck','**ANGLES')
!
             call read_angles(line,'**ANGLES',nat,thrang,mgrps,msubg,  &
                              grptag,neiang)
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
       subroutine read_molrep(key,blck,nat,tgrp,grptag,body,grps,subg, &
                              atms,mbody,mgrps,msubg,matms)
!
       use lengths, only: leninp,lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                         ::  blck    !  Block name
       character(len=lenline),intent(inout)                ::  key     !
       character(len=leninp),intent(inout)                 ::  tgrp    !  Groups file title
       integer,intent(in)                                  ::  nat     !  Number of atoms in the molecule
       character(len=lentag),dimension(nat),intent(inout)  ::  grptag  !  Names of the groups
       integer,dimension(nat),intent(inout)                ::  body    !  Number of groups in each body
       integer,dimension(nat),intent(inout)                ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(inout)                ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(inout)                ::  atms    !  Atoms identifier
       integer,intent(inout)                               ::  mbody   !  Number of bodies
       integer,intent(inout)                               ::  mgrps   !  Number of groups
       integer,intent(inout)                               ::  msubg   !  Number of subgroups
       integer,intent(inout)                               ::  matms   !  
!
! Local variables
!
       integer                                             ::  posi    !
!
! Reading MOLREP block sections 
! -----------------------------
!
       mbody    = 0
       mgrps    = 0
       msubg    = 0
       matms    = 0
!
       body(:)  = 0
       grps(:)  = 0
       subg(:)  = 0
       atms(:)  = 0
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.TITLE')
!~              write(*,*) 
!~              write(*,*) '  Reading .TITLE option'
!~              write(*,*)
!
             read(uniinp,*) tgrp
             tgrp = adjustl(tgrp)
!
             call findline(key,'blck','**MOLREP')
!
           case ('*BODY')
             mbody = mbody + 1
!
             call findline(key,'sect','*BODY')
!
!~              write(*,*) 
!~              write(*,*) '  Reading *BODY section'
!~              write(*,*)
!
             call read_body(key,'*BODY',blck,nat,grptag,body,grps,     &
                            subg,atms,mbody,mgrps,msubg,matms)
!
           case ('**END')
!~              write(*,*) 
!~              write(*,*) 'Exiting from **MOLREP block'
!~              write(*,*)
             return
           case default
             call unksect(key,blck)
         end select  
       end do
!
       return
       end subroutine read_molrep
!
!======================================================================!
!
       subroutine read_body(key,sect,blck,nat,grptag,body,grps,subg,   &
                            atms,mbody,mgrps,msubg,matms)
!
       use lengths, only: lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                         ::  sect    !  Section name
       character(len=*),intent(in)                         ::  blck    !  Block name
       character(len=lenline),intent(inout)                ::  key     !  
       integer,intent(in)                                  ::  nat     !  Number of atoms in the molecule
       character(len=lentag),dimension(nat),intent(inout)  ::  grptag  !  Names of the groups
       integer,dimension(nat),intent(inout)                ::  body    !  Number of groups in each body
       integer,dimension(nat),intent(inout)                ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(inout)                ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(inout)                ::  atms    !  Atoms identifier
       integer,intent(inout)                               ::  mbody   !  Number of bodies
       integer,intent(inout)                               ::  mgrps   !  Number of groups
       integer,intent(inout)                               ::  msubg   !  Number of subgroups
       integer,intent(inout)                               ::  matms   !
!
! Local variables
!
       integer                                             ::  posi    !
!
! Reading BODY section options 
! ----------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different block sections       
         select case (key)
           case ('.GRP','.GROUP')
             body(mbody) = body(mbody) + 1
!
!~              write(*,*) 
!~              write(*,*) '    Reading .GRP option'
!~              write(*,*)
!
             call read_grps('.GRP',nat,grptag,grps,subg,atms,mgrps,    &
                            msubg,matms)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.SGRP','.SUBGROUP')
             body(mbody) = body(mbody) + 1
!
!~              write(*,*) 
!~              write(*,*) '    Reading .SGRP option'
!~              write(*,*)
             call read_grps('.SGRP',nat,grptag,grps,subg,atms,mgrps,   &
                            msubg,matms)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect,blck)
         end select  
       end do
!
       return
       end subroutine read_body
!
!======================================================================!
!
       subroutine read_grps(opt,nat,grptag,grps,subg,atms,mgrps,       &
                            msubg,matms)
!
       use lengths, only: leninp,lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                         ::  opt     !  Option name
       integer,intent(in)                                  ::  nat     !  Number of atoms in the molecule
       character(len=lentag),dimension(nat),intent(inout)  ::  grptag  !  Names of the groups
       integer,dimension(nat),intent(inout)                ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(inout)                ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(inout)                ::  atms    !  Atoms identifier
       integer,intent(inout)                               ::  mgrps   !  Number of groups
       integer,intent(inout)                               ::  msubg   !  NUmber of subgroups
       integer,intent(inout)                               ::  matms   !
!
! Local variables
!
       character(len=lenline)                              ::  line    ! 
       character(len=lenline)                              ::  key     ! 
       character(len=leninp)                               ::  arg     !  
       integer                                             ::  natgrp  !
       integer                                             ::  posi    ! 
       integer                                             ::  io      !  Input/Output status
       integer                                             ::  old     !
!
! Reading GRPS option keywords 
! ----------------------------
!
       do
         read(uniinp,'(A)',iostat=io) line  
!
         if ( io /= 0 ) call endopt(opt)
! If reads a white line then reads the next line
         if ( len_trim(line) == 0 ) cycle
! Processing the line read 
         call chkcomment(line,key)
! If the line just contains a comment reads the next line
         if ( len_trim(key) == 0 ) cycle
!
! Procesing the keywords line
!
         mgrps = mgrps + 1  
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
             call errkey('option',opt)
           end if
! Saving the arguments
           select case (key)
             case ('atoms')
               call chkkeyarg(key,line,arg)
!
               read(arg,*) natgrp
!
               if ( trim(opt) .eq. '.GRP' ) then
                 msubg       = msubg + 1
                 grps(mgrps) = 1
                 subg(msubg) = natgrp
               else if ( trim(opt) .eq. '.SGRP' ) then
                 old   = msubg
                 msubg = msubg + natgrp
                 grps(mgrps) = natgrp
                 subg(old+1:msubg) = 1
               end if    
!
               read(uniinp,*) atms(matms+1:matms+natgrp)
!
               matms = matms + natgrp
!        
             case ('name')
               call chkkeyarg(key,line,arg)
               grptag(mgrps) = arg(:lentag)
!
             case default
               call unkkeysect(key,opt)
           end select  
!
           key = line
!
         end do
       end do
!
       return
       end subroutine read_grps
!
!======================================================================!
!
       subroutine read_threshold(key,blck,nat,thr,mgrps,grptag)
!
       use lengths, only: lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                      ::  blck    !
       character(len=lenline),intent(inout)             ::  key     !
       real(kind=8),dimension(nat,nat),intent(inout)    ::  thr     !
       character(len=lentag),dimension(nat),intent(in)  ::  grptag  !
       integer,intent(in)                               ::  nat     !
       integer,intent(in)                               ::  mgrps   !  Number of groups

!
! Local variables
!
       character(len=lentag)                            ::  caux1   !
       character(len=lentag)                            ::  caux2   !
       real(kind=8)                                     ::  daux    !
       integer                                          ::  iaux1   !
       integer                                          ::  iaux2   !
       integer                                          ::  posi    !
!
! Reading THRESHOLD block options 
! -------------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different section options       
         select case (key)
           case ('.THR')
!~              write(*,*) 
!~              write(*,*) 'Reading .THR option'
!~              write(*,*)
!
             read(uniinp,*) daux    ! FLAG: check if a value is introduced
!
             thr(:,:) = daux
!
             call findline(key,'blck','**THRESHOLD')
!
           case ('.VALUES')
!~              write(*,*) 
!~              write(*,*) 'Reading .VALUES option'
!~              write(*,*)
!
             call findline(key,'opt','.VALUES')
!
             do
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               caux1 = key(:posi-1)
               iaux1 = findcv(mgrps,grptag,caux1)  
               if ( iaux1 .eq. 0 ) call errkeychar('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
!
               caux2 = key(:posi-1)
               iaux2 = findcv(mgrps,grptag,caux2)  
               if ( iaux1 .eq. 0 ) call errkeychar('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               key = key(:posi-1)
               read(key,*) daux
!
               thr(iaux1,iaux2) = daux
               thr(iaux2,iaux1) = daux
!
               call findline(key,'opt','.VALUES')
!
               if ( (uppercase(key(1:1)).eq.'*') .or.                  &
                                     (uppercase(key(1:1)).eq.'.') ) exit
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
       subroutine read_interthreshold(key,blck,nat,thr,mgrps,grptag)
!
       use lengths, only: lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                         ::  blck    !
       character(len=lenline),intent(inout)                ::  key     !
       real(kind=8),dimension(nat,nat),intent(inout)       ::  thr     !
       character(len=lentag),dimension(nat),intent(in)     ::  grptag  !
       integer,intent(in)                                  ::  nat     !
       integer,intent(in)                                  ::  mgrps   !  Number of groups

!
! Local variables
!
       character(len=lentag)                               ::  caux1   !
       character(len=lentag)                               ::  caux2   !
       real(kind=8)                                        ::  daux    !
       integer                                             ::  iaux1   !
       integer                                             ::  iaux2   !
       integer                                             ::  posi    !
!
! Reading THRESHOLD block options 
! -------------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different section options       
         select case (key)
           case ('.THR')
!~              write(*,*) 
!~              write(*,*) 'Reading .THR option'
!~              write(*,*)
!
             read(uniinp,*) daux    ! FLAG: check if a value is introduced
!
             thr(:,:) = daux
!
             call findline(key,'blck','**INTERTHRESHOLD')
!
           case ('.VALUES')
!~              write(*,*) 
!~              write(*,*) 'Reading .VALUES option'
!~              write(*,*)
!
             call findline(key,'opt','.VALUES')
!
             do
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               caux1 = key(:posi-1)
               iaux1 = findcv(mgrps,grptag,caux1)  
               if ( iaux1 .eq. 0 ) call errkeychar('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
!
               caux2 = key(:posi-1)
               iaux2 = findcv(mgrps,grptag,caux2)  
               if ( iaux1 .eq. 0 ) call errkeychar('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               key = key(:posi-1)
               read(key,*) daux
!
               thr(iaux1,iaux2) = daux
               thr(iaux2,iaux1) = daux
!
               call findline(key,'opt','.VALUES')
!
               if ( (uppercase(key(1:1)).eq.'*') .or.                  &
                                     (uppercase(key(1:1)).eq.'.') ) exit
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
       end subroutine read_interthreshold
!
!======================================================================!
!
       subroutine read_angles(key,blck,nat,thr,mgrps,msubg,grptag,     &
                              neiang)
!
       use parameters
       use lengths, only: lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                         ::  blck    !
       character(len=lenline),intent(inout)                ::  key     !
       real(kind=8),dimension(nat,nat),intent(inout)       ::  thr     !
       character(len=lentag),dimension(nat),intent(in)     ::  grptag  !
       integer,dimension(nat),intent(inout)                ::  neiang  !
       integer,intent(in)                                  ::  nat     !
       integer,intent(in)                                  ::  mgrps   !  Number of groups
       integer,intent(in)                                  ::  msubg   !  Number of subgroups

!
! Local variables
!
       character(len=lentag)                               ::  caux1   !
       character(len=lentag)                               ::  caux2   !
       real(kind=8)                                        ::  daux    !
       integer                                             ::  iaux1   !
       integer                                             ::  iaux2   !
       integer                                             ::  posi    !
!
! Reading THRESHOLD block options 
! -------------------------------
!
       do
! Changing lowercase letters by uppercase letters 
         key = uppercase(key)
! Keeping just the first string
         posi = index(key,' ')           
         if ( posi .gt. 0 ) key = key(:posi-1)
! Reading the different section options       
         select case (key)
           case ('.NEI','.NEIGHBOUR','.NEIGHBOURS')
!~              write(*,*) 
!~              write(*,*) 'Reading .NEIGHBOURS option'
!~              write(*,*)
!
             read(uniinp,*) neiang(:msubg)    ! FLAG: check if a value is introduced
!
             call findline(key,'blck','**ANGLES')
!
           case ('.THR')
!~              write(*,*) 
!~              write(*,*) 'Reading .THR option'
!~              write(*,*)
!
             read(uniinp,*) daux    ! FLAG: check if a value is introduced
!
             if ( daux .gt. zero ) thr(:,:) = pi - daux*pi/180.0d0
!
             call findline(key,'blck','**ANGLES')
!
           case ('.VALUES')
!~              write(*,*) 
!~              write(*,*) 'Reading .VALUES option'
!~              write(*,*)
!
             call findline(key,'opt','.VALUES')
!
             do
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               caux1 = key(:posi-1)
               iaux1 = findcv(mgrps,grptag,caux1)  
               if ( iaux1 .eq. 0 ) call errkeychar('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
!
               caux2 = key(:posi-1)
               iaux2 = findcv(mgrps,grptag,caux2)  
               if ( iaux1 .eq. 0 ) call errkeychar('option','.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               key = key(:posi-1)
               read(key,*) daux
!~ write(*,*) 'INPUT VALUE',daux
!
               thr(iaux1,iaux2) = pi - daux*pi/180.0d0
               thr(iaux2,iaux1) = pi - daux*pi/180.0d0
!~ write(*,*) 'THRE VALUE',thr(iaux2,iaux1)*180.0d0/pi
!
               call findline(key,'opt','.VALUES')
!
               if ( (uppercase(key(1:1)).eq.'*') .or.                  &
                                     (uppercase(key(1:1)).eq.'.') ) exit
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
       end subroutine read_angles
!
!======================================================================!
!
       end module input_section
!
!======================================================================!
