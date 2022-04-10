!======================================================================!
!
       module input_section
       implicit none
!
       include 'info.h'
!
       private
       public  ::  read_inp,                                           &
                   read_gro
!
       contains
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
! ----------------------------------
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
       subroutine read_inp(inp,nat,tgrp,grptag,body,grps,subg,igrps,   &
                           nbody,ngrps,nsubg,thr)
!
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)             ::  inp     !  General input file name
       character(len=leninp),intent(out)            ::  tgrp    !  Groups file title
       character(len=8),dimension(nat),intent(out)  ::  grptag  !  Names of the groups
       real(kind=8),dimension(nat,nat),intent(out)  ::  thr     !  Distance threshold
       integer,dimension(nat),intent(out)           ::  body    !  Number of groups in each body
       integer,dimension(nat),intent(out)           ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(out)           ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(out)           ::  igrps   !  Atoms identifier
       integer,intent(out)                          ::  nbody   !  Number of bodies
       integer,intent(out)                          ::  ngrps   !  Number of groups
       integer,intent(out)                          ::  nsubg   !  Number of subgroups
       integer,intent(in)                           ::  nat     !  Number of atoms in the molecule
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
       nbody    = nat
       ngrps    = nat
       nsubg    = nat
!
       do i = 1, nat
         write(grptag(i),'(I8)') i
         grptag(i) = 'Atom-'//trim(adjustl(grptag(i)))
         body(i)   = i
         grps(i)   = i
         subg(i)   = i
         igrps(i)  = i
       end do
!
! Reading general input file
! --------------------------
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
       if ( io .ne. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  Missing input file'
         write(*,*)
         write(*,'(3X,A)')    'Input file '//trim(inp)//' not foun'//  &
                                            'd in the current directory'
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
           case ('**MOLREP','**MOLREPRE')
!~              write(*,*) 
!~              write(*,*) 'Reading **MOLREP block'
!~              write(*,*) 
!
             call findline(line,'blck','**MOLREP')
!
             call read_molrep(line,'**MOLREP',nat,tgrp,grptag,body,    &
                              grps,subg,igrps,nbody,ngrps,nsubg)      
!      
           case ('**THRESH','**THRESHOLD')
!~              write(*,*) 
!~              write(*,*) 'Reading **THRESHOLD block'
!~              write(*,*) 
!
             call findline(line,'blck','**THRESHOLD')
!
             call read_threshold(line,'**THRESHOLD',nat,thr,ngrps,     &
                                 grptag)
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
                              igrps,nbody,ngrps,nsubg)
!
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                  ::  blck    !  Block name
       character(len=leninp),intent(inout)          ::  key     !
       character(len=leninp),intent(out)            ::  tgrp    !  Groups file title
       integer,intent(in)                           ::  nat     !  Number of atoms in the molecule
       character(len=8),dimension(nat),intent(out)  ::  grptag  !  Names of the groups
       integer,dimension(nat),intent(out)           ::  body    !  Number of groups in each body
       integer,dimension(nat),intent(out)           ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(out)           ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(out)           ::  igrps   !  Atoms identifier
       integer,intent(out)                          ::  nbody   !  Number of bodies
       integer,intent(out)                          ::  ngrps   !  Number of groups
       integer,intent(out)                          ::  nsubg   !  Number of subgroups
!
! Local variables
!
       character(len=leninp)                        ::  line    !
       integer                                      ::  posi    !
       integer                                      ::  aux     !
       integer                                      ::  io      !  Input/Output status
       integer                                      ::  i,j,k   !  Indexes
!
! Reading MOLREP block sections 
! -----------------------------
!  
       tgrp = 'MOLREP Title'
!
       body(:)  = 0
       grps(:)  = 0
       subg(:)  = 0
       igrps(:) = 0
!
       nbody = 0
       ngrps = 0
       nsubg = 0
!
       aux = 0
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
             nbody = nbody + 1
!
             call findline(key,'sect','*BODY')
!
!~              write(*,*) 
!~              write(*,*) '  Reading *BODY section'
!~              write(*,*)
!
             call read_body(key,'*BODY',nat,grptag,body,grps,subg,     &
                            igrps,nbody,ngrps,nsubg,aux)
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
       subroutine read_body(key,sect,nat,grptag,body,grps,subg,igrps,  &
                            nbody,ngrps,nsubg,aux)
!
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                    ::  sect    !  Section name
       character(len=leninp),intent(inout)            ::  key     !  
       integer,intent(in)                             ::  nat     !  Number of atoms in the molecule
       character(len=8),dimension(nat),intent(inout)  ::  grptag  !  Names of the groups
       integer,dimension(nat),intent(inout)           ::  body    !  Number of groups in each body
       integer,dimension(nat),intent(inout)           ::  grps    !  Number of subgroups in each group
       integer,dimension(nat),intent(inout)           ::  subg    !  Number of atoms in each subgroup
       integer,dimension(nat),intent(inout)           ::  igrps   !  Atoms identifier
       integer,intent(inout)                          ::  nbody   !  Number of bodies
       integer,intent(inout)                          ::  ngrps   !  Number of groups
       integer,intent(inout)                          ::  nsubg   !  Number of subgroups
       integer,intent(inout)                          ::  aux     !
!
! Local variables
!
       character(len=leninp)                          ::  line    !
       integer                                        ::  posi    !
       integer                                        ::  io      !  Input/Output status
       integer                                        ::  i,j,k   !  Indexes
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
             body(nbody) = body(nbody) + 1
!
!~              write(*,*) 
!~              write(*,*) '    Reading .GRP option'
!~              write(*,*)
!
             call read_grps('.GRP',nat,grptag,grps,subg,igrps,ngrps,   &
                            nsubg,aux)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case ('.SGRP','.SUBGROUP')
             body(nbody) = body(nbody) + 1
!
!~              write(*,*) 
!~              write(*,*) '    Reading .SGRP option'
!~              write(*,*)
             call read_grps('.SGRP',nat,grptag,grps,subg,igrps,ngrps,  &
                            nsubg,aux)
!
             call findline(key,'sect',sect)
             if ( key(1:1) .eq. '*' ) return
!
           case default
             call unkopt(key,sect)
         end select  
       end do
!
       return
       end subroutine read_body
!
!======================================================================!
!
       subroutine read_grps(opt,nat,grptag,grps,subg,igrps,ngrps,      &
                            nsubg,aux)
!
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                    ::  opt     !  Option name
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
! Reading GRPS option keywords 
! ----------------------------
!
       do
         read(uniinp,'(A)',iostat=io) line  
!
         if ( io /= 0 ) call endopt(opt)
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
                 nsubg       = nsubg + 1
                 grps(ngrps) = 1
                 subg(nsubg) = natgrp
               else if ( trim(opt) .eq. '.SGRP' ) then
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
               call unkkey(key,opt)
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
       subroutine read_threshold(key,blck,nat,thr,ngrps,grptag)
!
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                    ::  blck    !
       character(len=leninp),intent(inout)            ::  key     !
       real(kind=8),dimension(nat,nat),intent(inout)  ::  thr     !
       character(len=8),dimension(nat),intent(in)     ::  grptag  !
       integer,intent(in)                             ::  nat     !
       integer,intent(in)                             ::  ngrps   !  Number of groups

!
! Local variables
!
       character(len=leninp)                          ::  line    !
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
               iaux1 = findcv(ngrps,grptag,caux1)  
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
               iaux2 = findcv(ngrps,grptag,caux2)  
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
       end module input_section
!
!======================================================================!
