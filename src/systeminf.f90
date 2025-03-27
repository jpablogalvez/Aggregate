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
       type(groinp),dimension(:),allocatable  ::  sys   !  Monomer information
       type(repre),dimension(:),allocatable   ::  rep   !  Topological representations
       type(xtcfile)                          ::  xtcf  !  Trajectory information
!
! System size information
!
       integer,dimension(:),allocatable       ::  nnode    !  Total number of molecules of each size
       integer,dimension(:),allocatable       ::  inode    !  Initial molecule of each size
       integer,dimension(:),allocatable       ::  nat      !  Number of atoms in each molecule type
       integer,dimension(:),allocatable       ::  iat      !  Initial atoms of each molecule type
       integer,dimension(:),allocatable       ::  natms    !
       integer,dimension(:),allocatable       ::  iatms    !
       integer,dimension(:),allocatable       ::  ngrps    !  Number of active atoms in each representation type
       integer,dimension(:),allocatable       ::  igrps    !  Initial active atom of each representation type
       integer                                ::  mtype    !  Number of molecule types
       integer                                ::  mnode    !  Total number of molecules
       integer                                ::  maxat    !  Total number of atoms in the system
       integer                                ::  mat      !  Total number of atoms in the molecules
       integer                                ::  matms    !  Total number of subgroups in the molecules
       integer                                ::  mgrps    !  
!
       end module systeminf
!
!======================================================================!
