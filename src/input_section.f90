!======================================================================!
!
       module input_section
!
       implicit none
!
       private
       public  ::  read_inp
!
       contains
!
!======================================================================!
!
       subroutine read_inp(inp,ntype,rep,mat,thr,thr2,thrang)
!
       use datatypes,  only:  repre
       use lengths,    only:  leninp,lentag,lenline
       use printings
       use utils
!
       implicit none
!
! Input/output variables
!
       type(repre),dimension(ntype),intent(inout)        ::  rep     !  Topological representations
       character(len=leninp),intent(in)                  ::  inp     !  General input file name
       real(kind=8),dimension(mat,mat),intent(out)       ::  thr     !  Distance threshold
       real(kind=8),dimension(mat,mat),intent(out)       ::  thr2    !  Distance threshold
       real(kind=8),dimension(mat,mat),intent(out)       ::  thrang  !  Angle threshold
       integer,intent(in)                                ::  ntype   !
       integer,intent(in)                                ::  mat     !
!
! Local variables
!
       character(len=lentag),dimension(mat)              ::  grptag   !  Names of the groups
       character(len=lenline)                            ::  line    !
       character(len=lenline)                            ::  key     !
       integer                                           ::  imol    !
       integer                                           ::  itag    !
       integer                                           ::  io      !  Input/Output status
       integer                                           ::  i,j     !  Indexes
!
! Setting defaults
!
       thr(:,:)  = 0.0d0
       thr2(:,:) = 0.0d0
!
       thrang(:,:) = 0.0d0
!
       do j = 1, ntype
!
         rep(j)%tgrp = 'MOLREP Title'
!
         rep(j)%mbody    = rep(j)%nat
         rep(j)%mgrps    = rep(j)%nat
         rep(j)%msubg    = rep(j)%nat
         rep(j)%matms    = rep(j)%nat
! 
         do i = 1, rep(j)%nat
           write(rep(j)%grptag(i),'(I8)') i
           rep(j)%grptag(i) = 'Atom-'//trim(adjustl(rep(j)%grptag(i)))
           rep(j)%body(i)   = i
           rep(j)%grps(i)   = i
           rep(j)%subg(i)   = i
           rep(j)%atms(i)   = i
         end do
!
       end do
!
! Reading general input file
! --------------------------
!
       imol = 1
       itag = 0
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
             call read_molrep(line,'**MOLREP',rep(imol)%nat,           & 
                              rep(imol)%tgrp,rep(imol)%grptag,         &
                              rep(imol)%nbody,rep(imol)%ngrps,         &
                              rep(imol)%nsubg,rep(imol)%atms,          &
                              rep(imol)%mbody,rep(imol)%mgrps,         &
                              rep(imol)%msubg,rep(imol)%matms)      
!            
             grptag(itag+1:itag+rep(imol)%mgrps) = rep(imol)%grptag(:rep(imol)%mgrps)
!
             itag = itag + rep(imol)%mgrps
             imol = imol + 1 
!      
           case ('**THRESH','**THRESHOLD')
!~              write(*,*) 
!~              write(*,*) 'Reading **THRESHOLD block'
!~              write(*,*) 
!
             call findline(line,'blck','**THRESHOLD')
!
             call read_threshold(line,'**THRESHOLD',mat,thr,itag,grptag(:itag))
!
           case ('**INTERTHRESH','**INTERTHRESHOLD')
!~              write(*,*) 
!~              write(*,*) 'Reading **THRESHOLD block'
!~              write(*,*) 
!
             call findline(line,'blck','**INTERTHRESHOLD')
!
             call read_interthreshold(line,'**INTERTHRESHOLD',         &
                                      mat,thr2,itag,grptag(:itag))
!
           case ('**ANGLES','**THRANG','**THRANGLE','**THREANGLES')
!~              write(*,*) 
!~              write(*,*) 'Reading **ANGLES block'
!~              write(*,*) 
!
             call findline(line,'blck','**ANGLES')
!
             call read_angles(line,'**ANGLES',mat,thrang,itag,grptag(:itag))
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
       subroutine read_threshold(key,blck,nat,thr,ntag,grptag)
!
       use lengths, only: lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lentag),dimension(ntag),intent(in)    ::  grptag  !
       character(len=lenline),intent(inout)                ::  key     !
       character(len=*),intent(in)                         ::  blck    !
       real(kind=8),dimension(nat,nat),intent(inout)       ::  thr     !
       integer,intent(in)                                  ::  nat     !
       integer,intent(in)                                  ::  ntag    !

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
               iaux1 = findcv(ntag,grptag,caux1)
               if ( iaux1 .eq. 0 ) call errkeychar('option',caux1,'.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
!
               caux2 = key(:posi-1)
               iaux2 = findcv(ntag,grptag,caux2)
               if ( iaux2 .eq. 0 ) call errkeychar('option',caux2,'.VALUES')
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
!~ write(*,*) 'adding',iaux1,iaux2,daux
               thr(iaux1,iaux2) = daux
               thr(iaux2,iaux1) = daux
!
               call findline(key,'opt','.VALUES')
!
               if ( (uppercase(key(1:1)).eq.'*') .or.                  &
                                     (uppercase(key(1:1)).eq.'.') ) exit

             end do
!~ stop
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
       subroutine read_interthreshold(key,blck,nat,thr,ntag,grptag)
!
       use lengths, only: lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lentag),dimension(ntag),intent(in)    ::  grptag  !
       character(len=lenline),intent(inout)                ::  key     !
       character(len=*),intent(in)                         ::  blck    !
       real(kind=8),dimension(nat,nat),intent(inout)       ::  thr     !
       integer,intent(in)                                  ::  nat     !
       integer,intent(in)                                  ::  ntag    !

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
               iaux1 = findcv(ntag,grptag,caux1)  
               if ( iaux1 .eq. 0 ) call errkeychar('option',caux1,'.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
               posi  = scan(key,' ') 
!
               caux2 = key(:posi-1)
               iaux2 = findcv(ntag,grptag,caux2)  
               if ( iaux2 .eq. 0 ) call errkeychar('option',caux2,'.VALUES')
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
       subroutine read_angles(key,blck,nat,thr,ntag,grptag)
!
       use parameters
       use lengths,    only:  lentag,lenline
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=lentag),dimension(ntag),intent(in)    ::  grptag  !
       character(len=lenline),intent(inout)                ::  key     !
       character(len=*),intent(in)                         ::  blck    !
       real(kind=8),dimension(nat,nat),intent(inout)       ::  thr     !
       integer,intent(in)                                  ::  nat     !
       integer,intent(in)                                  ::  ntag    !

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
!~            case ('.NEI','.NEIGHBOUR','.NEIGHBOURS')
!             write(*,*) 
!             write(*,*) 'Reading .NEIGHBOURS option'
!             write(*,*)
!~ !
!~              read(uniinp,*) neiang(:msubg)    ! FLAG: check if a value is introduced
!~ !
!~              call findline(key,'blck','**ANGLES')
!~ !
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
!
! Extract position of first label
!
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               caux1 = key(:posi-1) 
               iaux1 = findcv(ntag,grptag,caux1)
               if ( iaux1 .eq. 0 ) call errkeychar('option',caux1,'.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
! Extract position of second label
!
               posi  = scan(key,' ') 
!
               caux2 = key(:posi-1)
               iaux2 = findcv(ntag,grptag,caux2)
               if ( iaux2 .eq. 0 ) call errkeychar('option',caux2,'.VALUES')
!
               key   = key(posi+1:)   
               key   = adjustl(key)
               if ( len_trim(key) .eq. 0 ) call errkey('option',       &
                                                       '.VALUES')
!
! Extract threshold value
!
               posi  = scan(key,' ') 
               if ( posi .eq. 0 ) call errkey('option','.VALUES')
!
               key = key(:posi-1)
               read(key,*) daux
!
               thr(iaux1,iaux2) = pi - daux*pi/180.0d0
               thr(iaux2,iaux1) = pi - daux*pi/180.0d0
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
