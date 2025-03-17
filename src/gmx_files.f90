!======================================================================!
!
       module gmx_files
!
       use utils,    only:  find_key,below_key,find_last,print_badread
       use lengths,  only:  leninp,lenarg,lenline,lenlab
!
       implicit none
!
       private
       public  ::  read_gro,                                           &
                   top_parser,                                         &
!~                    read_top,                                           &
                   count_moltype
!
       contains
!
!======================================================================!
!
       subroutine read_gro(uni,conf,sys)
!
       use datatypes,  only:  groinp
       use printings,  only:  print_missinp
       use utils,      only:  uppercase
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  conf     !  Structure file name
       type(groinp),intent(inout)        ::  sys      !  System information
       integer,intent(in)                ::  uni      !  
!
! Local variables
!
       integer                           ::  io       !  Input/Output status
       integer                           ::  i        !  Index
!
! Reading Gromacs configuration file
! ----------------------------------
!
       open(unit=uni,file=trim(conf),action='read',                 &
            status='old',iostat=io)
!
       if ( io .ne. 0 ) call print_missinp(conf)
!
       read(uni,'(A)') sys%title
       read(uni,*)     sys%nat
!
       allocate(sys%renum(sys%nat)   ,  &
                sys%rename(sys%nat)  ,  &
                sys%atname(sys%nat)  ,  &
                sys%atnum(sys%nat)   ,  &
                sys%mass(sys%nat))
!
       do i = 1, sys%nat       ! TODO: check if lines are correctly read (fixed format)
         read(uni,'(I5,2A5,I5,3F8.3)') sys%renum(i),   &
                                       sys%rename(i),  &
                                       sys%atname(i),  &
                                       sys%atnum(i)
       end do
       read(uni,*) sys%latvec
!
       close(uni)
!
       return
       end subroutine read_gro
!
!======================================================================!
!
       subroutine top_parser(uni,top,ntype,nnode)
!
       use systeminf,  only:  sys
       use printings,  only:  print_missinp,print_end
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)             ::  top      !
       integer,dimension(ntype),intent(out)         ::  nnode    !  
       integer,intent(in)                           ::  uni      !  
       integer,intent(in)                           ::  ntype    !  
!
! Local variables
!
       character(len=lenline)                       ::  line     !  Line read
       character(len=lenline)                       ::  key      !  Line read
       character(len=10)                            ::  str1     !  
       character(len=10)                            ::  str2     !  
       character(len=10)                            ::  str3     !  
       integer                                      ::  keylen   !
       integer                                      ::  io       !
       integer                                      ::  i,j      !
!
! Reading Gromacs topology information
! ------------------------------------
!
       open(unit=uni,file=trim(top),action='read',                     &    
            status='old',form='formatted',iostat=io)
       if ( io .ne. 0 ) call print_missinp(top)
!
! First round to allocate information
! ...................................
!
       do i = 1, ntype    
!
! Counting atoms information
!
         call below_key(uni,'[ atoms ]',io)
         if ( io .ne. 0 ) call print_badread(top,'[ atoms ]')
!
         sys(i)%nat = 0
         do
           read(uni,'(A)',iostat=io) line
           if ( io /= 0 ) exit
           if ( len_trim(line) .eq. 0 ) cycle
           line = adjustl(line)
           if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
           if ( line(1:1) .eq. '[' ) exit  ! Finished when reads [ bonds ] or [ moleculetype ]
           sys(i)%nat = sys(i)%nat + 1
         end do
!
         allocate(sys(i)%renum(sys(i)%nat)   ,  &
                  sys(i)%rename(sys(i)%nat)  ,  &
                  sys(i)%atname(sys(i)%nat)  ,  &
                  sys(i)%atnum(sys(i)%nat)   ,  &
                  sys(i)%mass(sys(i)%nat))
!
! Counting bonds information
!
         if ( sys(i)%nat .gt. 1 ) then ! Read bonds only when dealing with molecules
!
           sys(i)%nbond = 0
           do
             read(uni,'(A)',iostat=io) line
             if ( io /= 0 ) exit
             if ( len_trim(line) .eq. 0 ) cycle
             line = adjustl(line)
             if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
             if ( line(1:1) .eq. '[' ) exit
             sys(i)%nbond = sys(i)%nbond + 1
           end do
!
           allocate(sys(i)%ibond(2,sys(i)%nbond))
!
         else
!
           sys(i)%nbond = 0
!
         end if
!
       end do
!
! Second round to read information
! ................................
!
! Reading moleculetype information
!
       rewind(uni)
!
       i = 0
       do while ( i .le. ntype )
!
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit                                                 ! 
!~ write(*,*) trim(line)
         if ( len_trim(line) .eq. 0 )  cycle                                 !  Skip white lines
         line = adjustl(line)
         if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle !  Skip comments
!
! Reading moleculetype information
!
         key = '[ moleculetype ]'
         keylen = len_trim(key)
!
         if ( len_trim(line) .ge. keylen ) then
           if ( line(:keylen) .eq. trim(key) ) then
!
!~ write(*,*) 'Reading moleculetype information'
!
             i = i + 1
!
             do
               read(uni,'(A)',iostat=io) line
               if ( io /= 0 ) call print_badread(top,'[ moleculetype ]')
!~ write(*,*) trim(line)
               if ( len_trim(line) .eq. 0 ) cycle
               line = adjustl(line)
               if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
               exit
             end do
!
             read(line,*) sys(i)%title 
!
           end if
         end if
!
! Reading atoms information
!
         key = '[ atoms ]'  ! TODO: check if atoms section is missing
         keylen = len_trim(key)
!
         if ( len_trim(line) .ge. keylen ) then
           if ( line(:keylen) .eq. trim(key) ) then
!
!~ write(*,*) 'Reading atoms information'
!
             j = 0
             do while ( j .lt. sys(i)%nat )
!
               read(uni,'(A)',iostat=io) line
               if ( io /= 0 ) call print_badread(top,'[ atoms ]')
!~ write(*,*) trim(line)
               if ( len_trim(line) .eq. 0 ) cycle
               line = adjustl(line)
               if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
!
               j = j + 1
!
               read(line,*) sys(i)%atnum(j),str1,sys(i)%renum(j),      &
                            sys(i)%rename(j),sys(i)%atname(j),         &
                            str2,str3,sys(i)%mass(j) 
!
             end do
!
! Reading bonds information
!
             if ( sys(i)%nat .gt. 1 ) then ! Read bonds only when dealing with molecules
!
!~ write(*,*) 'Reading bonds information'
! 
               key = '[ bonds ]'
               keylen = len_trim(key)
!
               do
!
                 read(uni,'(A)',iostat=io) line
                 if ( io /= 0 ) call print_badread(top,'[ bonds ]')
!~ write(*,*) trim(line)
                 if ( len_trim(line) .eq. 0 )  cycle                                 !  Skip white lines
                 line = adjustl(line)
                 if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle !  Skip comments
!
                 if ( len_trim(line) .ge. keylen ) then
                   if ( line(:keylen) .eq. trim(key) ) then
                     exit
                   end if
                 end if
!
               end do
!
               j = 0
               do while ( j .lt. sys(i)%nbond )
!
                 read(uni,'(A)',iostat=io) line
                 if ( io /= 0 ) call print_badread(top,'[ bonds ]')
!~ write(*,*) trim(line)
                 if ( len_trim(line) .eq. 0 ) cycle
                 line = adjustl(line)
                 if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
!
                 j = j + 1
!
                 read(line,*) sys(i)%ibond(:,j)
!
               end do
!
             end if
!
           end if
         end if
!
       end do
!
!~        do i = 1, ntype
!~ write(*,*)
!~ write(*,*) 'Reading moleculetype information'
!~          key = '[ moleculetype ]'
!~          keylen = len_trim(key)
!~          io     = 1
!~ !
!~          do
!~            read(uni,'(A)',iostat=io) line
!~ write(*,*)'reading:',trim(line),io!,trim(ioerrmsg)
!~            if ( io /= 0 ) exit
!~            if ( len_trim(line) .eq. 0 )  cycle
!~            if ( len_trim(line) .lt. keylen ) cycle
!~ !
!~            line = adjustl(line)
!~            if ( line(:keylen) .eq. trim(key) ) then
!~              io = 0
!~              exit
!~            end if
!~          end do
!~ !
!~          if ( io .ne. 0 ) call print_badread(top,'[ moleculetype ]')
!~ !
!~          do
!~            read(uni,'(A)',iostat=io) line
!~            if ( io /= 0 ) exit
!~            if ( len_trim(line) .eq. 0 ) cycle
!~            line = adjustl(line)
!~            if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
!~            exit
!~          end do
!~ !
!~          read(line,*) sys(i)%title   
!~ !
!~ ! Reading atoms information
!~ !
!~ write(*,*) 'Reading atoms information'
!~          key = '[ atoms ]'
!~          keylen = len_trim(key)
!~          io     = 1
!~ !
!~          do
!~            read(uni,'(A)',iostat=io) line
!~ write(*,*)'reading:',trim(line),io!,trim(ioerrmsg)
!~            if ( io /= 0 ) exit
!~            if ( len_trim(line) .eq. 0 )  cycle
!~            if ( len_trim(line) .lt. keylen ) cycle
!~ !
!~            line = adjustl(line)
!~            if ( line(:keylen) .eq. trim(key) ) then
!~              io = 0
!~              exit
!~            end if
!~          end do
!~ !
!~          if ( io .ne. 0 ) call print_badread(top,'[ atoms ]')
!~ !
!~          j = 0
!~          do while ( j .lt. sys(i)%nat )
!~ !
!~            read(uni,'(A)',iostat=io) line
!~            if ( io /= 0 ) exit
!~            if ( len_trim(line) .eq. 0 ) cycle
!~            line = adjustl(line)
!~            if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
!~ !
!~            j = i + 1
!~ !
!~            read(line,*) sys(i)%atnum(j),str1,sys(i)%renum(j),          &
!~                         sys(i)%rename(j),sys(i)%atname(j),             &
!~                         str2,str3,sys(i)%mass(j) 
!~ !
!~          end do
!~ !
!~ ! Reading bonds information
!~ !
!~          if ( sys(i)%nat .gt. 1 ) then ! Read bonds only when dealing with molecules
!~ !
!~            j = 0
!~            do while ( j .lt. sys(i)%nbond )
!~ !
!~              read(uni,'(A)',iostat=io) line
!~              if ( io /= 0 ) exit
!~              if ( len_trim(line) .eq. 0 ) cycle
!~              line = adjustl(line)
!~              if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
!~ !
!~              j = j + 1
!~ !
!~              read(line,*) sys(i)%ibond(:,j)
!~ !
!~            end do
!~ !
!~          end if
!~ !
!~        end do
!
! Reading molecules information
!
       rewind(uni)
!
       call below_key(uni,'[ molecules ]',io)
       if ( io .ne. 0 ) call print_badread(top,'[ molecules ]')
!
       i = 0
       do while ( i .lt. ntype )
!
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).eq.';') .or. (trim(line(1:1)).eq.'#') ) cycle
!
         i = i + 1
!
         read(line,*) str1,nnode(i)
!
         if ( trim(str1) .ne. trim(sys(i)%title) ) then
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Molecule type not correct'
           write(*,*)
           write(*,'(3X,A)') 'Molecule type defined: '//trim(sys(i)%title)
           write(*,'(3X,A)') 'Please, check the following definiti'//  &
                                     'on in [ molecules ]: '//trim(str1)
           write(*,*)
           write(*,'(3X,A)') 'Possibly, the molecules are defined '//  &
                             'in a different order, or there has b'//  &
                                               'een a misspelling error'
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,*)  
           call print_end()
         end if
!
       end do
!
       close(uni)
!
       return
       end subroutine top_parser
!
!======================================================================!
!
!~        subroutine read_top(uni,nat,top,ntype,itype,ctype,attype,ptype, &
!~                            ljeps,ljsig,charge,nbfunc,comrule,genpair,  &
!~                            fudgelj,fudgeqq,resname,nrexcl)
!~ !
!~        use printings,  only:  print_missinp
!~ !
!~        implicit none
!~ !
!~ ! Input/output variables
!~ !
!~        character(len=leninp),intent(in)                  ::  top      !
!~        integer,intent(in)                                ::  nat      !
!~        integer,intent(in)                                ::  uni      !  
!~ !
!~        character(len=lenlab),intent(out)                 ::  genpair  !
!~        real(kind=8),intent(out)                          ::  fudgelj  !
!~        real(kind=8),intent(out)                          ::  fudgeqq  !
!~        integer,intent(out)                               ::  nbfunc   !
!~        integer,intent(out)                               ::  comrule  !
!~        integer,intent(out)                               ::  nrexcl   !
!~ !
!~        character(len=lenlab),intent(inout)               ::  resname  !
!~        character(len=1),dimension(nat),intent(out)       ::  ptype    !
!~        character(len=lenlab),dimension(nat),intent(out)  ::  ctype    !
!~        character(len=lenlab),dimension(nat),intent(out)  ::  attype   !  
!~        real(kind=8),dimension(nat),intent(out)           ::  ljeps    !
!~        real(kind=8),dimension(nat),intent(out)           ::  ljsig    !
!~        real(kind=8),dimension(nat),intent(out)           ::  charge   !
!~        integer,dimension(nat),intent(out)                ::  itype    !
!~        integer,intent(out)                               ::  ntype    !
!~ !
!~ ! Local variables
!~ !
!~        character(len=lenline)                            ::  line     !  Line read
!~        integer                                           ::  io       !
!~ !
!~ ! Reading Gromacs topology information
!~ ! ------------------------------------
!~ !
!~        open(unit=uni,file=trim(top),action='read',                     &    
!~             status='old',iostat=io)
!~ !
!~        if ( io .ne. 0 ) call print_missinp(top)
!~ !
!~ ! Reading defaults section
!~ !
!~        call read_defaults(uni,top,nbfunc,comrule,genpair,fudgelj,fudgeqq)
!~ !
!~        rewind(uni)
!~ !
!~ ! Reading atomtypes section
!~ ! 
!~        call read_attype(uni,top,nat,ntype,ptype,attype,ljeps,ljsig)
!~ !
!~        rewind(uni)
!~ !
!~ ! Reading moleculetype section
!~ ! 
!~        call read_moltype(uni,top,resname,nrexcl)
!~ !
!~        rewind(uni)
!~ !
!~ ! Reading atoms section
!~ ! 
!~        call read_atoms(uni,top,nat,ntype,ctype,itype,charge)
!~ !
!~        close(uni)
!~ !
!~        return
!~        end subroutine read_top
!
!======================================================================!
!
       subroutine read_defaults(uni,top,nbfunc,comrule,genpair,        &
                                fudgelj,fudgeqq)
!
       implicit none
!
!
! Input/output variables
!
       integer,intent(in)                                ::  uni      !  
       character(len=leninp),intent(in)                  ::  top      !
!
       character(len=lenlab),intent(out)                 ::  genpair  !
       real(kind=8),intent(out)                          ::  fudgelj  !
       real(kind=8),intent(out)                          ::  fudgeqq  !
       integer,intent(out)                               ::  nbfunc   !
       integer,intent(out)                               ::  comrule  !
!
! Local variables
!
       character(len=lenline)                            ::  line     !  Line read
       integer                                           ::  io       !
!
! Reading defaults section
!
       line = '[none]'
       call find_key(uni,'[ defaults ]',line,io)
       if ( io .ne. 0 ) call print_badread(top,'[ defaults ]')
!
       do
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
         exit
       end do
!
       read(line,*) nbfunc,comrule,genpair,fudgelj,fudgeqq
!
       return
       end subroutine read_defaults
!
!======================================================================!
!
       subroutine read_attype(uni,top,nat,ntype,ptype,attype,bond,     &
                              mass,charge,ljeps,ljsig)
!
       implicit none
!
!
! Input/output variables
!
       character(len=leninp),intent(in)                  ::  top      !
       integer,intent(in)                                ::  uni      !  
       integer,intent(in)                                ::  nat      !  
!
       character(len=1),dimension(nat),intent(out)       ::  ptype    !
       character(len=lenlab),dimension(nat),intent(out)  ::  attype   !  
       real(kind=8),dimension(nat),intent(out)           ::  bond     !
       real(kind=8),dimension(nat),intent(out)           ::  mass     !
       real(kind=8),dimension(nat),intent(out)           ::  charge   !
       real(kind=8),dimension(nat),intent(out)           ::  ljeps    !
       real(kind=8),dimension(nat),intent(out)           ::  ljsig    !
       integer,intent(out)                               ::  ntype    !  
!
! Local variables
!
       character(len=lenline)                            ::  line     !  Line read
       integer                                           ::  io       !
!
! Reading atomtypes section
!
       line = '[none]'
       call find_key(uni,'[ atomtypes ]',line,io)
       if ( io .ne. 0 ) call print_badread(top,'[ atomtypes ]')
!
       ntype = 0
       do
!
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
         if ( line(1:1) .eq. '[' ) exit
!
         ntype = ntype + 1
!
         read(line,*) attype(ntype),bond(ntype),mass(ntype),           &  ! TODO: breaks with different format
                      charge(ntype),ptype(ntype),                      &
                      ljsig(ntype),ljeps(ntype) 
!
       end do
!
       return
       end subroutine read_attype
!
!======================================================================!
!
       subroutine count_moltype(uni,top,ntype)
!
       use printings, only:  print_missinp
!
       implicit none
!
!
! Input/output variables
!
       character(len=leninp),intent(in)                  ::  top     !
       integer,intent(in)                                ::  uni     !  
       integer,intent(out)                               ::  ntype   !
!
! Local variables
!
       character(len=lenline)                            ::  line    !  Line read
       integer                                           ::  io      !
!
! Counting moleculetype section
!
       open(unit=uni,file=trim(top),action='read',                   &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(top)
!
       line = '[none]'
       call find_key(uni,'[ moleculetype ]',line,io)
       if ( io .ne. 0 ) call print_badread(top,'[ moleculetype ]')
!
       ntype = 1
       do
         line = '[none]'
         call find_key(uni,'[ moleculetype ]',line,io)
         if ( io .ne. 0 ) exit
!
         ntype = ntype + 1
       end do
!
       close(uni)
!
       return
       end subroutine count_moltype
!
!======================================================================!
!
       subroutine read_moltype(uni,top,resname,nrexcl)
!
       implicit none
!
!
! Input/output variables
!
       integer,intent(in)                                ::  uni      !  
       character(len=leninp),intent(in)                  ::  top      !
!
       character(len=lenlab),intent(inout)               ::  resname  !
       integer,intent(out)                               ::  nrexcl   !
!
! Local variables
!
       character(len=lenline)                            ::  line     !  Line read
       character(len=lenline)                            ::  str1     !  
       integer                                           ::  io       !
!
! Reading moleculetype section
!
       line = '[none]'
       call find_key(uni,'[ moleculetype ]',line,io)
       if ( io .ne. 0 ) call print_badread(top,'[ moleculetype ]')
!
       do
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
         exit
       end do
!
       read(line,*) str1,nrexcl
!
       if ( len_trim(resname) .eq. 0 ) resname = trim(adjustl(str1))
!
       return
       end subroutine read_moltype
!
!======================================================================!
!
       subroutine read_atoms(uni,top,nat,attype,ntype,ctype,itype,     &
                             atnr,resnr,residue,atom,cgnr,charge,mass)
!
       use utils,     only: findcv
       use printings, only: print_end
!
       implicit none
!
!
! Input/output variables
!
       character(len=leninp),intent(in)                  ::  top      !
       character(len=lenlab),dimension(nat),intent(in)   ::  attype   !
       integer,intent(in)                                ::  uni      !  
       integer,intent(in)                                ::  nat      !  
       integer,intent(in)                                ::  ntype    !
!
       character(len=lenlab),dimension(nat),intent(out)  ::  residue  !
       character(len=lenlab),dimension(nat),intent(out)  ::  atom     !
       character(len=lenlab),dimension(nat),intent(out)  ::  ctype    !
       real(kind=8),dimension(nat),intent(out)           ::  charge   !
       real(kind=8),dimension(nat),intent(out)           ::  mass     !
       integer,dimension(nat),intent(out)                ::  cgnr     !
       integer,dimension(nat),intent(out)                ::  atnr     !
       integer,dimension(nat),intent(out)                ::  resnr    !
       integer,dimension(nat),intent(out)                ::  itype    !
!
! Local variables
!
       character(len=lenline)                            ::  line     !  Line read
       integer                                           ::  io       !
       integer                                           ::  i        !
!
! Reading atoms section
!
       line = '[none]'
       call find_key(uni,'[ atoms ]',line,io)
       if ( io .ne. 0 ) call print_badread(top,'[ atoms ]')
!
       itype(:) = 0
!
       i = 0
       do while ( i .lt. nat )
!
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
!
         i = i + 1
!
         read(line,*) atnr(i),ctype(i),resnr(i),residue(i),atom(i),    &
                                               cgnr(i),charge(i),mass(i) ! TODO: read top mass or use qc output mass
!
         io = findcv(ntype,attype(:ntype),ctype(i))
!
         if ( io .le. 0 ) then
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Atomtype not specified'
           write(*,*)
           write(*,'(3X,A)') 'Please, check the following atom: '//    &
                                                                ctype(i)
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,*)  
           call print_end()
         end if
!
         itype(i) = io
!
       end do
!
       return
       end subroutine read_atoms
!
!======================================================================!
!
       subroutine count_bond(uni,intop,nbond)
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  intop   !
       integer,intent(out)               ::  nbond   !
       integer,intent(in)                ::  uni     !
!
! Local variables
!
       character(len=lenline)            ::  line     !  Line read
       integer                           ::  io       !
!
! Reading Gromacs bonds information
! ---------------------------------
!
! Counting the number of bonds in the molecule
!
       line = '[none]'
       call find_key(uni,'[ bonds ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ bonds ]')
!
       nbond = 0
       do
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
         if ( line(1:1) .eq. '[' ) exit
         nbond = nbond + 1
       end do
!
       return
       end subroutine count_bond
!
!======================================================================!
!
       subroutine read_bond(uni,intop,nbond,ibond,fbond,rbond,kbond)
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)           ::  intop   !
       real(kind=8),dimension(nbond),intent(out)  ::  rbond   !
       real(kind=8),dimension(nbond),intent(out)  ::  kbond   !
       integer,dimension(2,nbond),intent(out)     ::  ibond   !
       integer,dimension(nbond),intent(out)       ::  fbond   !
       integer,intent(in)                         ::  nbond   !
       integer,intent(in)                         ::  uni     !
!
! Local variables
!
       character(len=lenline)                     ::  line     !  Line read
       integer                                    ::  io       !
       integer                                    ::  i        !
!
! Reading Gromacs bonds information
! ---------------------------------
!
! Reading bonds section  ! TODO: only valid for harmonic bond and G96 bond interactions
! 
       line = '[none]'
       call find_key(uni,'[ bonds ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ bonds ]')
!
       i = 0
       do while ( i .lt. nbond )
!
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
!
         i = i + 1
!
         read(line,*) ibond(:,i),fbond(i),rbond(i),kbond(i)
!
       end do
!
       return
       end subroutine read_bond
!
!======================================================================!
!
       subroutine read_system(uni,top,sysname)
!
       implicit none
!
!
! Input/output variables
!
       integer,intent(in)                                ::  uni      !  
       character(len=leninp),intent(in)                  ::  top      !
       character(len=lenline),intent(inout)              ::  sysname  !
!
! Local variables
!
       character(len=lenline)                            ::  line     !  Line read
       integer                                           ::  io       !
!
! Reading moleculetype section
!
       line = '[none]'
       call find_key(uni,'[ system ]',line,io)
       if ( io .ne. 0 ) call print_badread(top,'[ system ]')
!
       do
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
         exit
       end do
!
       sysname = line
!
       return
       end subroutine read_system
!
!======================================================================!
!
       subroutine read_molecules(uni,top,ntype,nnode,resname)
!
       use printings,  only:  print_end
!
       implicit none
!
!
! Input/output variables
!
       integer,intent(in)                                  ::  uni      !    
       character(len=leninp),intent(in)                    ::  top      !
       integer,intent(in)                                  ::  ntype    !
       integer,dimension(ntype),intent(out)                ::  nnode    !
       character(len=lenlab),dimension(ntype),intent(in)   ::  resname  !
!
! Local variables
!
       character(len=lenline)                              ::  line     !  Line read
       character(len=lenlab)                               ::  str      !  
       integer                                             ::  io       !
       integer                                             ::  i        !
!
! Reading moleculetype section
!
       line = '[none]'
       call find_key(uni,'[ molecules ]',line,io)
       if ( io .ne. 0 ) call print_badread(top,'[ molecules ]')
!
       i = 0
       do while ( i .lt. ntype )
!
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (trim(line(1:1)).ne.';') .or. (trim(line(1:1)).ne.'#') ) cycle
!
         i = i + 1
!
         read(line,*) str,nnode(i)
!
         if ( trim(str) .ne. trim(resname(i)) ) then
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Molecule type not correct'
           write(*,*)
           write(*,'(3X,A)') 'Please, check the following definiti'//  &
                                      'on in [ molecules ]: '//trim(str)
           write(*,*)
           write(*,'(3X,A)') 'Possibly, the molecules are defined '//  &
                             'in a different order, or there has b'//  &
                                               'een a misspelling error'
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,*)  
           call print_end()
         end if
!
       end do
!
       return
       end subroutine read_molecules
!
!======================================================================!
!
       end module gmx_files
!
!======================================================================!
