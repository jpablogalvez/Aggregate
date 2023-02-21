!
! Declaration of time control variables
!
       integer       ::  t1,t2    !  CPU times
       integer       ::  t1read   !  Initial reading time
       integer       ::  t2read   !  Final reading time
       integer       ::  t1adj    !  Initial building time
       integer       ::  t2adj    !  Final building time
       integer       ::  t1bfs    !  Initial BFS time
       integer       ::  t2bfs    !  Final BFS time
       integer       ::  t1sort   !  Initial sorting time
       integer       ::  t2sort   !  Final sorting time
       integer       ::  t1pim    !  Initial PIM analysis time
       integer       ::  t2pim    !  Final PIM analysis time
       real(kind=8)  ::  tcpu     !  Total CPU time
       real(kind=8)  ::  tread    !  Total reading time
       real(kind=8)  ::  tadj     !  Total adjacency matrix building time
       real(kind=8)  ::  tbfs     !  Total BFS time
       real(kind=8)  ::  tsort    !  Total sorting time
       real(kind=8)  ::  tpim     !  Total PIM analysis time
!
! Declaration of system_clock variables 
!
       integer       ::  count_rate
       integer       ::  count_max 
