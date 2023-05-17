!======================================================================!
!
       module thresholds
!
       implicit none
!
! Interaction criteria information
! 
       real(kind=8),dimension(:,:),allocatable  ::  thr     !  Distance threshold
       real(kind=8),dimension(:,:),allocatable  ::  thr2    !  Distance threshold
       real(kind=8),dimension(:,:),allocatable  ::  thrang  !  Angles threshold
       real(kind=8)                             ::  neidis  !  Screening distance
       integer,dimension(:),allocatable         ::  neiang  !  First neighbour index
!
       end module thresholds
!
!======================================================================!
