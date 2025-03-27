!======================================================================!
!
       module thresholds
!
       implicit none
!
! Interaction criteria information
! 
       real(kind=4),dimension(:,:),allocatable  ::  thr     !  Distance threshold
       real(kind=4),dimension(:,:),allocatable  ::  thrang  !  Angles threshold
       real(kind=4)                             ::  neidis  !  Screening distance
       integer,dimension(:),allocatable         ::  neiang  !
!
       end module thresholds
!
!======================================================================!
