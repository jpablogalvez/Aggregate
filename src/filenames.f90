!======================================================================!
!
       module filenames
!
       use lengths,  only:  leninp,lenout    
!
       implicit none
!
! Input/output files
!
       character(len=leninp)  ::  traj    ! Trajectory file name
       character(len=leninp)  ::  weight  ! Weights file name
       character(len=leninp)  ::  conf    ! Configuration file name
       character(len=leninp)  ::  inp     ! General input file name
       character(len=lenout)  ::  outp    ! Output file name
!
       end module filenames
!
!======================================================================!
