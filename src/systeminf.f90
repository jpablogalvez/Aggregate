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
       type(groinp)   ::  sys   !  Monomer information
       type(xtcfile)  ::  xtcf  !  Trajectory information
!
       end module systeminf
!
!======================================================================!
