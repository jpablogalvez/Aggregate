!======================================================================!
!
       module systeminf
!
       use datatypes
       use xdr, only: xtcfile
!
       implicit none
!
! System information
!
       type(groinp),dimension(:),allocatable   ::  sys      !  Monomer information
       type(repre),dimension(:),allocatable    ::  rep      !  Topological representations
       type(adjinf),dimension(:),allocatable   ::  adjgrps  !  Templates in the grps representations
       type(adjinf),dimension(:),allocatable   ::  adjbody  !  Templates in the bodies representations
       type(xtcfile)                           ::  xtcf     !  Trajectory information
!
! System size information
!
       integer,dimension(:),allocatable       ::  nnode     !  Total number of molecules of each size
       integer,dimension(:),allocatable       ::  inode     !  Initial molecule of each size
       integer,dimension(:),allocatable       ::  nat       !  Number of atoms in each molecule type
       integer,dimension(:),allocatable       ::  iat       !  Initial atoms of each molecule type
       integer,dimension(:),allocatable       ::  natms     !
       integer,dimension(:),allocatable       ::  iatms     !
       integer,dimension(:),allocatable       ::  ngrps     !  Number of active atoms in each representation type
       integer,dimension(:),allocatable       ::  igrps     !  Initial active atom of each representation type
       integer,dimension(:),allocatable       ::  mmon      !  Number of monomers in the aggregate
       integer,dimension(:),allocatable       ::  mgrpsmon  !  Total number of grps per aggregate type
       integer,dimension(:),allocatable       ::  mbodymon  !  Total number of bodies per aggregate type
       integer,dimension(:,:),allocatable     ::  nmon      !  Number of monomers of each type
       integer,dimension(:,:),allocatable     ::  imon      !  Initial monomers of each type
       integer,dimension(:,:),allocatable     ::  ngrpsmon  !  Number of grps per monomer type
       integer,dimension(:,:),allocatable     ::  nbodymon  !  Number of bodies per monomer type
       integer,dimension(:,:),allocatable     ::  igrpsmon  !  Accumulation of grps per monomer type
       integer,dimension(:,:),allocatable     ::  ibodymon  !  Accumulation of bodies per monomer type
       integer                                ::  mtype     !  Number of molecule types
       integer                                ::  mnode     !  Total number of molecules
       integer                                ::  maxat     !  Total number of atoms in the system
       integer                                ::  mat       !  Total number of atoms in the monomers
       integer                                ::  matms     !  Total number of groups in the system
       integer                                ::  mgrps     !  Total number of groups in the monomers
!
       end module systeminf
!
!======================================================================!
