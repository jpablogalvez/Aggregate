!======================================================================!
!
       module properties
!
       implicit none
!
! Average properties
!
       real(kind=8),dimension(:,:,:),allocatable       ::  pim      !  Pairwise interaction matrix
       real(kind=8),dimension(:),allocatable           ::  pop      !  Populations
       real(kind=8),dimension(:),allocatable           ::  conc     !  Concentrations
       real(kind=8),dimension(:),allocatable           ::  frac     !  Molar fractions
       real(kind=8)                                    ::  cin      !  Stechiometric concentration
       real(kind=8)                                    ::  volu     !  Simulation box volume
!
       end module properties
!
!======================================================================!
