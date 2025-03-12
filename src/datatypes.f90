!======================================================================!
!
       module datatypes
!
       use lengths,  only:  leninp,lentag
!
       implicit none
!
       type groinp
         character(len=32)                          ::  fname   !  Input file
         character(len=72)                          ::  title   !  Title
         integer                                    ::  nat     !  Number of atoms
         real(kind=8)                               ::  totm    !  Total mass 
         integer,dimension(:),allocatable           ::  renum   !  Residue number
         character(len=5),dimension(:),allocatable  ::  rename  !  Residue name
         character(len=5),dimension(:),allocatable  ::  atname  !  Atom name
         integer,dimension(:),allocatable           ::  atnum   !  Atom number
         real(kind=8),dimension(:),allocatable      ::  mass    !  Atom mass
         real(kind=8),dimension(:,:),allocatable    ::  coord   !  Coordinates
         real(kind=8),dimension(:,:),allocatable    ::  vel     !  Velocities
         real(kind=8),dimension(3)                  ::  latvec  !  Box vectors
         integer,dimension(:,:),allocatable         ::  ibond   !  Bond information
         integer                                    ::  nbond   !  Number of bonds
         logical,dimension(:,:),allocatable         ::  adj     !  Adjacency matrix
       end type groinp  
!
       type repre
         character(len=leninp)                           ::  tgrp     !  Groups file title
         character(len=lentag),dimension(:),allocatable  ::  grptag   !  Names of the groups
         integer,dimension(:),allocatable                ::  body     !  Number of groups in each body
         integer,dimension(:),allocatable                ::  nbody    !  Number of groups in each body
         integer,dimension(:),allocatable                ::  ibody    !  Number of groups in each body
         integer,dimension(:),allocatable                ::  grps     !  Number of subgroups in each group
         integer,dimension(:),allocatable                ::  ngrps    !  Number of subgroups in each group
         integer,dimension(:),allocatable                ::  igrps    !  Number of subgroups in each group
         integer,dimension(:),allocatable                ::  subg     !  Number of atoms in each subgroup
         integer,dimension(:),allocatable                ::  nsubg    !  Number of atoms in each subgroup
         integer,dimension(:),allocatable                ::  isubg    !  Number of atoms in each subgroup
         integer,dimension(:),allocatable                ::  atms     !  Atoms identifier
         integer,dimension(:),allocatable                ::  neiang   !
         integer                                         ::  mbody    !  Number of bodies
         integer                                         ::  mgrps    !  Number of groups
         integer                                         ::  msubg    !  Number of subgroups
         integer                                         ::  matms    !  Number of interacting atoms in the monomer
         integer                                         ::  nat      !  Number of atoms in the monomer
         integer                                         ::  iat      !  Starting atom
       end type repre
!
       end module datatypes
!
!======================================================================!
