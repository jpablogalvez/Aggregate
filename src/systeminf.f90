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
       end module systeminf
!
!======================================================================!
