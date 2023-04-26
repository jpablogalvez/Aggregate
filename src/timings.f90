!======================================================================!
!
       module timings
!
       implicit none
!
       save
!
! Declaration of time control variables
!
       real(kind=8)  ::  tcpu     !  Total CPU time
       real(kind=8)  ::  tread    !  Total reading time
       real(kind=8)  ::  tadj     !  Total adjacency matrix building time
       real(kind=8)  ::  tbfs     !  Total BFS time
       real(kind=8)  ::  tsort    !  Total sorting time
       real(kind=8)  ::  tscrn    !  Total screening time
       real(kind=8)  ::  tlife    !  Total lifetimes time
       real(kind=8)  ::  tpim     !  Total PIM analysis time
       real(kind=8)  ::  tconf    !  Total conformational analysis time
!
! Declaration of system_clock variables 
!
       integer       ::  count_rate
       integer       ::  count_max 
!
       end module timings
!
!======================================================================!
